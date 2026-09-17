"""Harmonize CL-ontology GT labels and DG-scRNA/scType marker-set predictions
to a common major-lineage vocabulary for cross-method benchmarking."""
import re
LINEAGES = ["T cell","B cell","NK cell","Plasma cell","Myeloid","Dendritic cell",
            "Mast cell","Endothelial","Fibroblast/Stromal","Epithelial","Glial/Neuronal"]
# tissue context: which parenchymal lineage malignant cells belong to per dataset
MALIGNANT_TO = {"breast_TNBC":"Epithelial","colorectal":"Epithelial","kidney_ccRCC":"Epithelial",
                "brain_GBM":"Glial/Neuronal","blood_DLBCL":"B cell"}
def _norm(s): return re.sub(r'[^a-z0-9 ]',' ',str(s).lower())
def to_lineage(label, dataset=None):
    """Map a single label string to a major lineage. dataset optional (for malignant context)."""
    s=_norm(label)
    # malignant / tumor cells -> tissue parenchyma
    if any(k in s for k in ["malignant","neoplastic","abnormal","tumor","cancer stem","tumour"]):
        return MALIGNANT_TO.get(dataset,"Epithelial")
    # immune
    if "mast" in s: return "Mast cell"
    if "plasmacytoid dendritic" in s or "pdc" in s: return "Dendritic cell"
    if "dendritic" in s or re.search(r'\bdc\b',s): return "Dendritic cell"
    if "plasma cell" in s or "plasmablast" in s: return "Plasma cell"
    if "nk t" in s or "nkt" in s: return "T cell"           # NKT -> T
    if "natural killer" in s or re.search(r'\bnk\b',s): return "NK cell"
    if re.search(r'\bt cell',s) or "t follicular" in s or "regulatory t" in s or \
       "cd4" in s or "cd8" in s or "thymocyte" in s or re.search(r'\bth\b',s) or \
       "t helper" in s or "cytotoxic t" in s or re.search(r'\btreg\b',s): return "T cell"
    if re.search(r'\bb cell',s) or "b lymph" in s or "germinal center" in s: return "B cell"
    if any(k in s for k in ["macrophage","monocyte","myeloid","kupffer","microglia",
                             "neutrophil","granulocyte","dendritic","langerhans","osteoclast"]):
        return "Myeloid"
    # endothelial
    if "endothel" in s or "lymphangio" in s or "vascular lymph" in s: return "Endothelial"
    # stromal / mesenchymal
    if any(k in s for k in ["fibroblast","pericyte","smooth muscle","perivascular","stromal",
                            "myofibroblast","stellate","mesenchym","myoepithelial"]):
        return "Fibroblast/Stromal"
    # glia / neuron
    if any(k in s for k in ["astrocyte","oligodendrocyte","neuron","opc","glia","glial",
                            "microglial","radial glia","neural","schwann","ependymal"]):
        return "Glial/Neuronal"
    # epithelial (broad: colonocyte, goblet, luminal, basal, secretory, parietal, mucous, acinar, ductal, etc.)
    if any(k in s for k in ["epitheli","colonocyte","goblet","enterocyte","enteroendocrine",
                            "crypt","tuft","secretory","luminal","basal","paneth","club","ciliated",
                            "parietal","mucous","acinar","ductal","follicular cell","hepatocyte",
                            "keratinocyte","alveolar","pneumocyte","best4","glandular","chief cell"]):
        return "Epithelial"
    return "Other"
def predname_to_lineage(pred, dataset=None):
    """DG-scRNA/scType predictions look like 'CellMarker_AllTissues_normal+Blood+CD8 T cell'
    or scType lineage names. Take the trailing lineage token then map."""
    p=str(pred)
    if p in ("Unknown","Undecided","nan","Unassigned","",None): return p if p in ("Unknown","Undecided") else "Unknown"
    tok = p.split("+")[-1] if "+" in p else p
    return to_lineage(tok, dataset)
if __name__=="__main__":
    # coverage test against all GT labels
    gt=['B cell','BEST4+ colonocyte','CD4-positive, alpha-beta T cell','CD4-positive, alpha-beta memory T cell',
    'CD8-positive, alpha-beta cytotoxic T cell','T cell','T follicular helper cell','abnormal cell',
    'activated CD8-positive, alpha-beta T cell','astrocyte','basal cell of epithelium of lobular bronchiole',
    'basal-myoepithelial cell of mammary gland','capillary endothelial cell','colon goblet cell','colonocyte',
    'conventional dendritic cell','early colonocyte','endothelial cell','enteroendocrine cell of colon',
    'fibroblast of breast','intestinal crypt stem cell of colon','luminal epithelial cell of mammary gland',
    'macrophage','malignant cell','mast cell','mature NK T cell','mucous cell of stomach','muscle fibroblast',
    'myeloid cell','myoepithelial cell','natural killer cell','neoplastic cell','neuron','oligodendrocyte',
    'oligodendrocyte precursor cell','parietal cell','pericyte','perivascular cell','plasma cell',
    'plasmacytoid dendritic cell','plasmacytoid dendritic cell, human','regulatory T cell',
    'salivary gland glandular cell','secretory cell','tuft cell of colon','vascular associated smooth muscle cell',
    'vascular lymphangioblast']
    unmapped=[]
    for g in gt:
        L=to_lineage(g,"colorectal" if "colon" in g.lower() else None)
        if L=="Other": unmapped.append(g)
    print("GT labels ->", {g:to_lineage(g) for g in gt})
    print("\nUNMAPPED (Other):", unmapped if unmapped else "NONE — full coverage")
