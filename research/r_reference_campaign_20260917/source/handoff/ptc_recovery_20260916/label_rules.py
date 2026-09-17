"""Name-only PTC ontology rules, frozen before evaluation; no TCR outcomes are used."""
import re

UNKNOWN={'','unknown','undecided','no annotation','unassigned','unlabeled','nan','na','none'}
TISSUE_PREFIXES={'Thymus','Lymph','Lymph node','Lymphoid tissue','Thyroid','Blood','Epithelium'}

def simplify(name):
    name=str(name)
    fields=name.split('+')
    if name.startswith('CellMarker_') and len(fields)>=4:
        name='+'.join(fields[3:])
    elif fields[0]=='NCOMMREFF' and len(fields)>1:
        name='+'.join(fields[1:])
    elif fields[0]=='cancer' and len(fields)>3:
        name='+'.join(fields[3:])
    elif fields[0] in TISSUE_PREFIXES and len(fields)>2:
        name='+'.join(fields[2:])
    return re.sub(r'\s+',' ',name.replace('.',' ').replace('_',' ')).strip()

def broad_lineage(name):
    s=simplify(name).lower()
    if s in UNKNOWN:return 'Unknown'
    if re.search(r'natural killer t|\bnkt\b',s):return 'NKT'
    if re.search(r'natural killer|\bnk\b',s):return 'NK'
    if re.search(r'\bt(?:\s|\(|-|/|[0-9]+\b|$)|\btreg\b|\bmait\b|\bctl\b',s):return 'T'
    if re.search(r'\bplasma(?:\s|$)|plasmablast|antibody secreting',s):return 'Plasma'
    if re.search(r'\bb(?:\s|\(|-|$)|\bbreg\b',s):return 'B'
    if re.search(r'macrophag|monocyt|myeloid|dendritic|neutrophil|granulocyte|basophil|mast cell|langerhans|microglia',s):return 'Myeloid'
    if re.search(r'endothelial|lymphatic',s):return 'Endothelial'
    if re.search(r'fibroblast|smooth muscle|pericyte|stromal|mesenchymal|myoid',s):return 'Stromal'
    if re.search(r'epithelial|follicular|basal cell|club cell|ciliated|ionocyte|lonocyte|goblet|secretory|mucous|hillock|deuterosomal|tuft|brush cell|neuroendocrine|cholangiocyte|corneocyte',s):
        # CD8 intraepithelial has no explicit T identity in its name: keep ambiguous.
        if 'cd8' in s:return 'Lymphoid_ambiguous'
        return 'Epithelial'
    if re.search(r'lymphocyte|lymphoid|immune cell|leukocyte|peripheral blood mononuclear|precursor memory',s):return 'Lymphoid_ambiguous'
    if re.search(r'cancer|malignant',s):return 'Tumor_unspecified'
    return 'Other'

def strict_T(name):
    return broad_lineage(name)=='T'

def legacy_broad_T_sensitivity(name):
    # Explicit reconstruction of the old permissive sensitivity; never strict T identity.
    return broad_lineage(name) in {'T','NKT','NK','Lymphoid_ambiguous'}

def T_subtype(name):
    if broad_lineage(name)!='T':return broad_lineage(name)
    s=simplify(name).lower()
    if re.search('treg|regulatory',s):return 'Treg'
    if re.search('gamma|γδ|delta',s):return 'gamma_delta_T'
    if 'mait' in s or 'mucosa-associated invariant' in s:return 'MAIT'
    if 'cd8' in s:return 'CD8_T'
    if 'cd4' in s:return 'CD4_T'
    if 'cytotoxic' in s:return 'Cytotoxic_T_unspecified'
    return 'T_unspecified'

MODULES={
 'T_core':['CD3D','CD3E','TRAC','TRBC1','TRBC2','CD247'],
 'NK_cytotoxic':['NKG7','GNLY','KLRD1','FCGR3A'],
 'B':['MS4A1','CD79A','CD79B','CD37'],
 'Plasma':['MZB1','JCHAIN','SDC1'],
 'Myeloid':['LYZ','LST1','TYROBP','FCER1G'],
 'Endothelial':['PECAM1','VWF','CLDN5','KDR'],
 'Epithelial':['EPCAM','KRT8','KRT18','KRT19'],
 'Thyroid':['TG','TPO','IYD','SLC5A5'],
 'Stromal':['COL1A1','COL1A2','DCN','LUM','ACTA2','TAGLN'],
 'Proliferation':['MKI67','TOP2A','STMN1','TYMS'],
 'Interferon':['ISG15','IFIT1','IFIT3','MX1','OAS1','STAT1'],
 'Stress':['FOS','JUN','HSPA1A','DDIT3'],
 'CD4':['CD4','IL7R','CCR7'],
 'CD8':['CD8A','CD8B','CCL5'],
 'Treg':['FOXP3','IL2RA','CTLA4']}

MAPPING_NOTES={
 'primary':'Explicit T names only; NKT,NK and ambiguous lymphocytes are separate.',
 'CD8 intraepithelial cell':'Ambiguous without explicit T identity; do not optimize mapping against TCR.',
 'NCOMMREFF+NK Cells':'Preserve NK label despite CD8A/NKG7 marker ambiguity.',
 'generic cancer/medullary/cortical names':'Do not infer epithelial identity solely from generic labels.',
 'legacy_broad':'T+NKT+NK+ambiguous lymphoid is a permissive sensitivity, not biological ground truth.',
 'RNA_support':'T_core at least2 genes detected; descriptive expression support, not independent truth or used in fitting.',
 'coarse_concordance':'S2/S3 source-label concordance only; not independent multi-class accuracy.'}
