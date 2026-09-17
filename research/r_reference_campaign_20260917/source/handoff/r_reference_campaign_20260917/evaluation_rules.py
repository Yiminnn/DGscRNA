"""Prespecified semantic equivalences. No metric-based relabelling or cluster matching."""
import re
from functools import lru_cache
UNKNOWN={'unknown','undecided','no_annotation','unassigned','nan',''}
SYNONYMS={
 'alpha':'alpha cell','pancreatic a cell':'alpha cell',
 'beta':'beta cell','type b pancreatic cell':'beta cell',
 'delta':'delta cell','pancreatic d cell':'delta cell',
 'gamma':'gamma cell','pp':'gamma cell','pp cell':'gamma cell','pancreatic pp cell':'gamma cell',
 'pancreatic gamma cell':'gamma cell','epsilon':'epsilon cell','pancreatic epsilon cell':'epsilon cell',
 'pancreatic alpha cell':'alpha cell','pancreatic beta cell':'beta cell','pancreatic delta cell':'delta cell',
 'pancreatic acinar cell':'acinar cell','acinar':'acinar cell','ductal':'ductal cell','pancreatic ductal cell':'ductal cell',
 'endothelial':'endothelial cell','mast':'mast cell','macrophage':'macrophage','schwann':'schwann cell',
 't cell':'t cell','plasmocyte':'plasma cell','b cell (plasmocyte)':'plasma cell',
 'nk cell':'natural killer cell','nkt cell':'natural killer t cell','natural killer t (nkt) cell':'natural killer t cell',
 'erythrocyte':'erythroid cell','red blood cell (erythrocyte)':'erythroid cell',
 'hspcs':'hematopoietic stem and progenitor cell',
 'hematopoietic stem/progenitor cell':'hematopoietic stem and progenitor cell',
 'cb cd34+':'cd34+ progenitor cell',
 'monocyte-derived dendritic cell':'monocyte-derived dendritic cell',
 'plasmacytoid dendritic cell(pdc)':'plasmacytoid dendritic cell',
 'plasmacytoid dendritic cell (pdc)':'plasmacytoid dendritic cell',
}

@lru_cache(maxsize=65536)
def canonical(label):
    s=str(label).strip().lower().replace('_',' ')
    for a,b in [('α','alpha'),('β','beta'),('γ','gamma'),('δ','delta'),('ε','epsilon')]:s=s.replace(a,b)
    s=re.sub(r'\s+',' ',s)
    s=re.sub(r'\s*\((?:alpha|beta|gamma|delta|epsilon) cell\)','',s)
    s=re.sub(r'\s*\((?:nk|nkt)\)','',s)
    s=re.sub(r'\bcells\b','cell',s)
    s=re.sub(r'\bmonocytes\b','monocyte',s)
    s=re.sub(r'\bprogenitors\b','progenitor',s)
    s=re.sub(r'\bery(throcytes)\b','erythrocyte',s)
    if s in UNKNOWN:return 'Unknown'
    return SYNONYMS.get(s,s)

@lru_cache(maxsize=65536)
def broad(label,dataset=None):
    """Secondary common-lineage endpoint; retained separately from curated labels."""
    s=canonical(label)
    if s=='Unknown':return s
    # Preserve the key multiclass endocrine and immune distinctions.
    for name in ['alpha','beta','gamma','delta','epsilon']:
        if s==name+' cell' or ('pancrea' in s and name in s):return name.capitalize()+' cell'
    if 'acinar' in s:return 'Acinar cell'
    if 'ductal' in s:return 'Ductal cell'
    if 'stellate' in s or s=='psc cell':return 'Stellate cell'
    if 'schwann' in s:return 'Schwann cell'
    if 'megakaryo' in s or 'platelet' in s:return 'Megakaryocyte/platelet'
    if 'erythro' in s or 'red blood' in s:return 'Erythroid'
    if 'hematopoietic' in s or s=='hspcs' or 'cd34' in s:return 'Hematopoietic progenitor'
    if 'malignant' in s or 'neoplastic' in s or 'cancer cell' in s or 'cancer stem' in s or s=='abnormal cell':return 'Malignant'
    if 'nkt' in s or 'natural killer t' in s or 'nk t' in s:return 'NKT cell'
    if 'cd4' in s and ('t ' in s or 't-' in s):return 'CD4 T cell'
    if 'cd8' in s and ('t ' in s or 't-' in s):return 'CD8 T cell'
    if 'plasmacytoid' in s and 'dendritic' in s:return 'Plasmacytoid DC'
    if 'monocyte-derived dendritic' in s:return 'Monocyte-derived DC'
    if 'cd14' in s and 'monocyte' in s:return 'CD14 monocyte'
    if 'cd16' in s and 'monocyte' in s:return 'CD16 monocyte'
    if 'astrocyte' in s:return 'Astrocyte'
    if 'oligodendrocyte precursor' in s or s=='opc':return 'OPC'
    if 'oligodendrocyte' in s:return 'Oligodendrocyte'
    if 'neuron' in s:return 'Neuron'
    # Never collapse every unmatched label into a mutually correct "Other" class.
    import sys
    from pathlib import Path
    p=str(Path(__file__).resolve().parents[1])
    if p not in sys.path:sys.path.insert(0,p)
    from harmonize import to_lineage
    value=to_lineage(s,dataset)
    return 'unmapped:'+s if value=='Other' else value

@lru_cache(maxsize=65536)
def ptc_general(label):
    if label in {'Unknown','Undecided','No_Annotation'}:return label
    fields=str(label).split('+')
    if fields[0]=='NCOMMREFF':return fields[1]
    if fields[0]=='cancer':return '+'.join(fields[3:])
    if fields[0] in ['CellMarker_normal','CellMarker_cancer']:return '+'.join(fields[3:])
    return '+'.join(fields[2:]) if len(fields)>1 else str(label)
