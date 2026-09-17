"""Verify and evaluate completed terminal conditions; fitting never imports this module."""
import json
import os
from pathlib import Path
import sys
from ptc_common import BASE,GROUPS,require_slurm,sha,utc,write_json,task_list,geometry_dir
from label_rules import broad_lineage,legacy_broad_T_sensitivity
from verify_terminal import verify_terminal

def evaluation_task(index):
    if index<496:
      t=task_list()[index];gd=geometry_dir(t)
      parts=[(gd/m/'score_RNA',gd/m,{**t,'clusterer':m,'correction':'per_sample_python'}) for m in t['clusterers']]
      return gd,parts,'single'
    if index<502:
      i=index-496;group=['MTN','TUT'][i//3];correction=['NONE','CCA','HARMONY'][i%3]
      gd=BASE/'pooled'/group/correction;parts=[]
      for space in ['PCA30','UMAP2']:
        for method in ['SNN','HDBSCAN_R']:
          for assay in (['RNA','integrated'] if correction=='CCA' else ['RNA']):
            part=gd/f'{space}_{method}'
            context=dict(group=group,correction=correction,space=space,clusterer=method,
              feature='group_selected2000',input_space='pca30' if space=='UMAP2' else 'direct',
              dr='UMAP' if space=='UMAP2' else 'PCA',dim=2 if space=='UMAP2' else 30,seed=42,
              family='legacy_integrated_scoring_and_DL' if assay=='integrated' else 'controlled_batch')
            parts.append((part/f'score_{assay}',part,context))
      return gd,parts,'pooled'
    i=index-502;sample=sum(GROUPS.values(),[])[i];group='MTN' if i<4 else 'TUT'
    gd=BASE/'matched_samples'/sample;parts=[]
    for space in ['PCA30','UMAP2']:
      for method in ['SNN','HDBSCAN_R']:
        part=gd/f'{space}_{method}'
        context=dict(sample=sample,group=group,correction='per_sample_NONE',space=space,clusterer=method,
           feature='group_selected2000',input_space='pca30' if space=='UMAP2' else 'direct',
           dr='UMAP' if space=='UMAP2' else 'PCA',dim=2 if space=='UMAP2' else 30,seed=42,family='matched_R_single')
        parts.append((part/'score_RNA',part,context))
    return gd,parts,'matched'

def metrics(pred,ref):
    import numpy as np
    from sklearn.metrics import f1_score,adjusted_rand_score
    names=np.unique(pred)
    mapping={v:broad_lineage(v) for v in names}
    broad=np.array([mapping[v] for v in pred],dtype=str)
    strict=broad=='T'
    permissive=np.isin(broad,['T','NKT','NK','Lymphoid_ambiguous'])
    unknown=broad=='Unknown'
    y2=ref.S2_terminal_broad.to_numpy(dtype=str);y3=ref.S3_terminal_broad.to_numpy(dtype=str)
    out=dict(n_cells=len(pred),n_called=int((~unknown).sum()),coverage=float((~unknown).mean()),
       n_strict_T=int(strict.sum()),n_NK=int((broad=='NK').sum()),n_NKT=int((broad=='NKT').sum()),
       n_ambiguous_lymphoid=int((broad=='Lymphoid_ambiguous').sum()),
       S2_broad_concordance=float((broad==y2).mean()),S3_broad_concordance=float((broad==y3).mean()),
       S2_macro_F1_reference_present=float(f1_score(y2,broad,labels=np.unique(y2),average='macro',zero_division=0)),
       S3_macro_F1_reference_present=float(f1_score(y3,broad,labels=np.unique(y3),average='macro',zero_division=0)),
       S2_annotation_ARI=float(adjusted_rand_score(y2,broad)))
    for tag,col in [('productive','TCR_cell_high_confidence_productive_TCR'),
                    ('any_contig','TCR_any_filtered_contig'),('S3_supplied','TCR_S3_supplied'),
                    ('paired_TRA_TRB','TCR_paired_productive_TRA_TRB')]:
      detected=ref[col].to_numpy(dtype=bool)
      for mapping,predicted in [('strict',strict),('permissive',permissive)]:
        tp=int((predicted&detected).sum());fn=int((~predicted&detected).sum());fp=int((predicted&~detected).sum())
        prefix=tag+'_'+mapping
        out.update({prefix+'_TP':tp,prefix+'_FN':fn,prefix+'_predicted_no_detection':fp,
           prefix+'_recall':tp/(tp+fn) if tp+fn else np.nan,
           prefix+'_detection_yield':tp/(tp+fp) if tp+fp else np.nan,
           prefix+'_apparent_F1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.0})
    support=ref.T_RNA_support_ge2.to_numpy(dtype=bool)
    undetected=~ref.TCR_cell_high_confidence_productive_TCR.to_numpy(dtype=bool)
    for name,mask in [('predicted_T',strict),('undetected_predicted_T',undetected&strict),('Unknown',unknown)]:
      out[name+'_n']=int(mask.sum())
      out[name+'_RNA_T_support']=float(support[mask].mean()) if mask.any() else np.nan
      out[name+'_median_nCount']=float(np.median(ref.nCount_RNA.to_numpy()[mask])) if mask.any() else np.nan
      out[name+'_median_nFeature']=float(np.median(ref.nFeature_RNA.to_numpy()[mask])) if mask.any() else np.nan
      out[name+'_median_percent_mt']=float(np.median(ref['percent.mt'].to_numpy()[mask])) if mask.any() else np.nan
    return out,broad

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    import torch
    from sklearn.metrics import adjusted_rand_score,homogeneity_completeness_v_measure
    torch.set_num_threads(4);torch.set_num_interop_threads(1)
    index=int(sys.argv[1] if len(sys.argv)>1 else os.environ['SLURM_ARRAY_TASK_ID'])
    gd,parts,mode=evaluation_task(index)
    assert (gd/'ANNOTATION_COMPLETE').exists()
    assert (BASE/'evaluation_reference/COMPLETE').exists()
    ref=pd.read_csv(BASE/'evaluation_reference/reference_cells.csv.gz',index_col=0,keep_default_na=False)
    dest=BASE/'evaluations'/f'task_{index:03d}';dest.mkdir(parents=True,exist_ok=True)
    if (dest/'COMPLETE').exists():
      m=json.loads((dest/'manifest.json').read_text());assert m['source_sha256']==sha(Path(__file__))
      assert (dest/'COMPLETE').read_text().strip()==sha(dest/'manifest.json');return
    rows=[];conditions=[];cluster_rows=[];compositions=[]
    for score,part,context in parts:
      assert (part/'CLUSTER_COMPLETE').read_text().strip()==sha(part/'cluster_manifest.json')
      clmeta=json.loads((part/'cluster_manifest.json').read_text())
      assert clmeta['clusters_sha256']==sha(part/'clusters.csv')
      assert (score/'REFINEMENT_COMPLETE').exists()
      sm=json.loads((score/'score_manifest.json').read_text())
      assert (score/'SCORE_COMPLETE').read_text().strip()==sha(score/'score_manifest.json')
      seeds=pd.read_csv(score/'initial_calls.csv.gz',keep_default_na=False,dtype=str).set_index('cell_id')
      cluster=pd.read_csv(part/'clusters.csv',dtype=str).set_index('cell_id').loc[seeds.index,'cluster']
      assert seeds.index.is_unique and set(cluster.index)==set(seeds.index)
      cr=ref.loc[seeds.index]
      base={k:context.get(k) for k in ['group','feature','input_space','dr','dim','seed','family','correction','clusterer','space','geometry_id']}
      base.update(mode=mode,scoring_assay=sm['assay'],precision_recovery=bool(clmeta.get('precision_recovery',False)),
                  task_index=index,partition_path=str(part.relative_to(BASE)))
      sample_positions=[(s,np.flatnonzero(cr['sample'].to_numpy()==s)) for s in cr['sample'].unique()]
      for sample,idx in sample_positions:
        rr=cr.iloc[idx];labels=cluster.iloc[idx].to_numpy()
        y=rr.S2_terminal_native.to_numpy()
        h,c,v=homogeneity_completeness_v_measure(y,labels)
        noise='0' if context['clusterer']=='HDBSCAN_R' else '-1'
        cluster_rows.append(dict(**base,sample=sample,patient=str(rr.patient.iloc[0]),n_cells=len(rr),
          cluster_S2_native_ARI=float(adjusted_rand_score(y,labels)),
          cluster_S2_native_V=float(v),n_clusters_in_sample=len(np.unique(labels)),
          noise_fraction=float((labels==noise).mean()) if 'HDBSCAN' in context['clusterer'] else 0.0))
      for aid,arm in sm['arms'].items():
        terminal=score/'refinement'/aid
        tm,z=verify_terminal(terminal,seeds[aid].to_numpy(dtype=str),seeds.index)
        assert tm['score_manifest_sha256']==sha(score/'score_manifest.json') and tm['arm']==arm
        cid=str(terminal.relative_to(BASE))
        info=dict(**base,condition_id=cid,library=arm['library'],cutoff=arm['cutoff'],
          DL_input=arm['DL_input'],cache_key=tm['cache_key'],dl_status=tm['dl_status'],
          actual_model_trained=tm['model_training_executed_for_cache'],n_known=tm['n_known'],
          n_pool=tm['n_pool'],n_training_classes=tm['n_training_classes'],terminal_sha256=tm['terminal_sha256'])
        conditions.append(info)
        for stage,key in [('marker_only_ablation','initial'),('terminal_DL090','final090'),('terminal_DL070_sensitivity','final070')]:
          pred=z[key]
          for sample,idx in sample_positions:
            rr=cr.iloc[idx]
            result,broad=metrics(pred[idx],rr)
            rows.append(dict(**info,sample=sample,patient=str(rr.patient.iloc[0]),stage=stage,**result))
            values,counts=np.unique(broad,return_counts=True)
            compositions.extend(dict(condition_id=cid,sample=sample,patient=str(rr.patient.iloc[0]),stage=stage,
                                     broad_lineage=name,n=int(n),fraction=float(n/len(rr))) for name,n in zip(values,counts))
        z.close()
    pd.DataFrame(rows).to_csv(dest/'sample_metrics.csv.gz',index=False)
    pd.DataFrame(conditions).to_csv(dest/'condition_verification.csv.gz',index=False)
    pd.DataFrame(cluster_rows).to_csv(dest/'cluster_metrics.csv.gz',index=False)
    pd.DataFrame(compositions).to_csv(dest/'lineage_composition.csv.gz',index=False)
    outputs={p.name:sha(p) for p in dest.glob('*.csv.gz')}
    m=dict(status='completed',task_index=index,mode=mode,n_conditions=len(conditions),n_metric_rows=len(rows),
      source_sha256=sha(Path(__file__)),verifier_sha256=sha(Path(__file__).with_name('verify_terminal.py')),
      label_rule_sha256=sha(Path(__file__).with_name('label_rules.py')),
      reference_manifest_sha256=sha(BASE/'evaluation_reference/manifest.json'),outputs=outputs,
      job=os.environ['SLURM_JOB_ID'],completed_at=utc())
    write_json(dest/'manifest.json',m);(dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')
    print('Verified and evaluated',len(conditions),'terminal conditions in task',index,flush=True)

if __name__=='__main__':run()
