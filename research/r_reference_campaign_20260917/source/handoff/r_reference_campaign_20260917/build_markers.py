"""Freeze tissue-informed CellMarker contexts before scoring any benchmark output."""
import hashlib,json,os,re
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
CONTEXTS={
 'brain_GBM':(['Brain'],['Blood','Blood vessel','Lymph node'],r'glioblast|glioma'),
 'breast_TNBC':(['Breast','Mammary gland'],['Blood','Lymph node','Blood vessel','Epithelium','Adipose tissue'],r'breast|mammary'),
 'colorectal':(['Intestine','Colon','Gastrointestinal tract'],['Liver','Lung','Peritoneum','Blood','Lymph node','Blood vessel','Epithelium'],r'colorect|colon|rectal'),
 'kidney_ccRCC':(['Kidney'],['Blood','Lymph node','Blood vessel','Epithelium'],r'renal|kidney'),
 'blood_DLBCL':(['Lymph node','Lymph','Lymphoid tissue','Blood'],['Skin','Intestine','Stomach','Airway','Bronchi','Bronchus','Blood vessel','Bone marrow','Spleen'],r'lymphoma|DLBCL'),
 'pancreas':(['Pancreas'],['Blood','Lymph node','Blood vessel','Epithelium','Intestine'],None),
 'immune_ALL_human':(['Blood','Bone marrow'],['Lymph node','Lymph','Lymphoid tissue','Thymus','Spleen'],None),
}

def hcl_primary(scope):
    rules=[('Adrenal',['Suprarenal gland','Adrenal gland']),('CordBlood',['Blood','Bone marrow','Umbilical cord']),
     ('BoneMarrow',['Bone marrow']),('SpinalCord',['Spinal cord','Brain']),('TemporalLobe|Cerebellum|Brain',['Brain']),
     ('Lung',['Lung']),('Intestine|Colon|Rectum|JeJunum|Duodenum|Epityphlon|Ileum',['Intestine','Colon','Gastrointestinal tract']),
     ('Kidney|Ureter',['Kidney']),('Pleura',['Pleura','Lung']),('Pancreas',['Pancreas']),('Muscle',['Muscle','Skeletal muscle']),
     ('Liver',['Liver']),('PeripheralBlood',['Blood']),('Spleen',['Spleen']),('Stomach',['Stomach']),
     ('MaleGonad',['Testis','Gonad']),('FemaleGonad',['Ovary','Gonad']),('Omentum',['Adipose tissue','Peritoneum']),
     ('Thyroid',['Thyroid']),('Esophagus',['Esophagus']),('Trachea',['Trachea','Airway']),('Chorionic|Placenta',['Placenta']),
     ('Gallbladder',['Gall bladder','Bile duct']),('Artery',['Artery','Blood vessel']),('Bladder',['Bladder']),
     ('Cervix',['Cervix','Uterine cervix']),('Heart',['Heart']),('Uterus',['Uterus','Endometrium']),
     ('Skin',['Skin','Epidermis']),('Fallopiantube',['Oviduct']),('Rib|Calvaria',['Bone']),('Thymus',['Thymus']),
     ('Prostate',['Prostate']),('Eyes',['Eye']),('HESC',['Embryo']),('Adipose',['Adipose tissue'])]
    for pattern,tissues in rules:
        if re.search(pattern,scope):return tissues
    raise ValueError(f'Unmapped HCL tissue {scope}')

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import pandas as pd
    src=ROOT/'handoff/refdb/Cell_marker_Human.xlsx'
    df=pd.read_excel(src).fillna('')
    for c in df:df[c]=df[c].astype(str).str.strip()
    df=df[(df.species.str.lower()=='human') & df.Symbol.ne('') & df.cell_name.ne('')].copy()
    df['native']=df.apply(lambda r:'CM2+'+'+'.join([r.cell_type,r.tissue_class,r.cancer_type or 'unspecified',r.cell_name]),axis=1)
    panels={};metadata=[]
    for name,g in df.groupby('native',sort=True):
        genes=sorted({s.strip() for value in g.Symbol for s in re.split('[,;]',value) if s.strip() and s.strip().lower()!='nan'})
        panels[name]=genes
        row=g.iloc[0]
        metadata.append(dict(native=name,cell_name=row.cell_name,tissue=row.tissue_class,context=row.cell_type,
          cancer_type=row.cancer_type,ontology_ids=';'.join(sorted(set(g.cellontology_id)-{''})),
          PMIDs=';'.join(sorted(set(g.PMID)-{''})),denominator=len(genes)))
    dest=OUT/'markers';dest.mkdir(exist_ok=True)
    pd.DataFrame(metadata).to_csv(dest/'native_panel_metadata.csv',index=False)
    units=[]
    inventory=json.loads((OUT/'dataset_inventory.json').read_text())
    for key in inventory:
        scopes=sorted(inventory[key]['metadata']['tissue_sample']['values']) if key=='HCL' else ['whole']
        for scope in scopes:
            unit=key if scope=='whole' else 'HCL__'+re.sub('[^A-Za-z0-9]+','_',scope)
            if key=='HCL':
                primary=hcl_primary(scope);support=['Blood','Blood vessel','Lymph node']
                if scope.startswith('Fetal') or scope in ['HESC','ChorionicVillus','Placenta']:support+=['Embryo','Placenta']
                disease=None
            else:primary,support,disease=CONTEXTS['pancreas' if key in ['baron_human','muraro','segerstolpe','xin'] else key]
            primary_mask=df.tissue_class.isin(primary)
            normal=df.cell_type.eq('Normal cell')
            selections={'CM2_primary_normal':primary_mask & normal}
            if disease:
                selections['CM2_primary_disease']=primary_mask & ~normal & df.cancer_type.str.contains(disease,case=False,regex=True)
                selections['CM2_primary_all_context']=primary_mask
            for tissue in support:
                if tissue in primary:continue
                mask=df.tissue_class.eq(tissue)
                if disease is None:mask &= normal
                selections['CM2_'+re.sub('[^A-Za-z0-9]+','_',tissue)]=mask
            union=df.tissue_class.isin(primary+support)
            if disease is None:union &= normal
            selections['CM2_related_union']=union
            selections['CM2_AllHuman']=pd.Series(True,index=df.index)
            libraries={};empty=[]
            for lib,mask in selections.items():
                names=sorted(df.loc[mask,'native'].unique())
                if not names:empty.append(lib);continue
                libraries[lib]={name:panels[name] for name in names}
            (dest/(unit+'.json')).write_text(json.dumps(libraries,ensure_ascii=False,indent=1)+'\n')
            units.append(dict(unit=unit,dataset=key,scope=scope,primary=primary,support=support,disease_pattern=disease,
                libraries={k:len(v) for k,v in libraries.items()},empty_contexts_not_fitted=empty,
                selection_basis='prespecified sampled anatomy and normal/disease state; no prediction metrics or expression filtering',
                full_gene_denominators_preserved=True))
    manifest=dict(source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
      script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),job=os.environ['SLURM_JOB_ID'],
      source_url='https://bio-bigdata.hrbmu.edu.cn/CellMarker2.0/index.html',units=units,
      marker_truth_independence='Marker evidence can cite source studies; source PMIDs retained. Not an independent-reference generalization claim.')
    (dest/'marker_roster.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print('Frozen marker rosters',len(units),'units; native panels',len(panels),flush=True)

if __name__=='__main__':run()
