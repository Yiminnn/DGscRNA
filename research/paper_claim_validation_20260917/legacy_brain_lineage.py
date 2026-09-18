"""Map panel names and curated labels onto ONE brain lineage vocabulary.

The curated labels in brain_GBM are 7 CL terms. Panel names come from CellMarker 2.0
brain rows and are finer. Both sides must land in the same vocabulary or macro-F1 is
meaningless. Anything unmappable becomes 'Other' and is scored as a miss, never dropped
silently.
"""
import re
LINEAGES = ["myeloid","neoplastic","OPC","oligodendrocyte","astrocyte","neuron","vascular","Other"]

_PANEL = {
    # myeloid / immune
    "microglial cell":"myeloid","m1 microglial cell":"myeloid","m2 microglial cell":"myeloid",
    "macrophage":"myeloid","m1 macrophage":"myeloid","m2 macrophage":"myeloid",
    "myeloid cell":"myeloid","dendritic cell":"myeloid","monocyte":"myeloid",
    "conventional dendritic cell 2a(cdc2a)":"myeloid","conventional dendritic cell 2b(cdc2b)":"myeloid",
    "pre-dendritic cell(pre-dc)":"myeloid","plasmacytoid dendritic cell":"myeloid",
    "t cell":"myeloid","cd4+ t cell":"myeloid","cd8+ t cell":"myeloid","b cell":"myeloid",
    "regulatory t(treg) cell":"myeloid","natural killer cell":"myeloid","mast cell":"myeloid",
    # malignant / progenitor-like
    "cancer cell":"neoplastic","cancer stem cell":"neoplastic","stem cell":"neoplastic",
    "mesenchymal cell":"neoplastic","radial glial cell":"neoplastic",
    "neuronal progenitor cell":"neoplastic","neural stem cell":"neoplastic",
    "progenitor cell":"neoplastic","glioma stem cell":"neoplastic",
    # glia
    "oligodendrocyte precursor cell":"OPC","oligodendrocyte progenitor cell":"OPC",
    "oligodendrocyte":"oligodendrocyte","mature oligodendrocyte":"oligodendrocyte",
    "newly formed oligodendrocyte":"oligodendrocyte",
    "astrocyte":"astrocyte","fibrous astrocyte":"astrocyte","protoplasmic astrocyte":"astrocyte",
    # neurons
    "neuron":"neuron","excitatory neuron":"neuron","inhibitory neuron":"neuron",
    "interneuron":"neuron","gabaergic neuron":"neuron","glutamatergic neuron":"neuron",
    "dopaminergic neuron":"neuron","purkinje cell":"neuron",
    # vasculature / stroma
    "endothelial cell":"vascular","pericyte":"vascular","fibroblast":"vascular",
    "smooth muscle cell":"vascular","vascular lymphangioblast":"vascular",
    "mural cell":"vascular","ependymal cell":"Other","choroid plexus cell":"Other",
}
_CURATED = {
    "myeloid cell":"myeloid","neoplastic cell":"neoplastic",
    "oligodendrocyte precursor cell":"OPC","oligodendrocyte":"oligodendrocyte",
    "astrocyte":"astrocyte","neuron":"neuron","vascular lymphangioblast":"vascular",
}

def _clean(s):
    s = str(s)
    m = re.search(r"\+([^+]+)$", s)          # strip '<db>_<src>+<tissue>+' prefix
    if m: s = m.group(1)
    return re.sub(r"\s+", " ", s.strip().lower())

# Pattern rules for the long tail. Ordered: first match wins, so put specific before generic.
# 'Lake et al.Science.Ex*/In*' are excitatory/inhibitory cortical neuron subtypes from that paper.
_PATTERNS = [
    (r"lake et al.*\.(ex|in)\d", "neuron"),
    (r"\b(glioblastoma|glioma)\b", "neoplastic"),
    (r"\bmalignant\b", "neoplastic"),
    (r"radi[ac]l glial", "neoplastic"),
    (r"intermediate progenitor", "neoplastic"),
    (r"\b(nkt|natural killer|cytotoxic|lymphocyte|memory t|exhausted|helper t|gamma delta)\b", "myeloid"),
    (r"\bt cell\b|\bb cell\b|\bplasma cell\b", "myeloid"),
    (r"microglia", "myeloid"),
    (r"macrophage|monocyte|dendritic", "myeloid"),
    # granulocytes ARE myeloid (CL:0000763 'myeloid cell' subsumes the granulocyte
    # lineage). Without this rule a neutrophil panel scored as an error even when the
    # gold label was 'myeloid cell': 641 cells in CM_bonemarrow. Megakaryocyte and
    # erythroid stay 'Other' -- CL nests them under myeloid too, but the curators of
    # this brain dataset plainly did not mean them by 'myeloid cell'.
    (r"neutrophil|granulocyte|\bMDSC\b|myeloid-derived suppressor", "myeloid"),
    (r"neuron|neural cell", "neuron"),
    (r"oligodendrocyte precursor|oligodendrocyte progenitor|\bopc\b", "OPC"),
    (r"oligodendrocyte", "oligodendrocyte"),
    (r"astrocyt", "astrocyte"),
    (r"endothelial|pericyte|mural|smooth muscle|fibroblast|stromal", "vascular"),
]

def panel_lineage(label):
    if label is None or str(label) in ("nan","Undecided","Unknown","noise",""):
        return "Other"
    c = _clean(label)
    if c in _PANEL: return _PANEL[c]
    for pat, lin in _PATTERNS:
        if re.search(pat, c): return lin
    return "Other"

def curated_lineage(label):
    return _CURATED.get(_clean(label), "Other")

def panel_source(label):
    m = re.search(r"_(normal|cancer)\+", str(label))
    return m.group(1) if m else None
