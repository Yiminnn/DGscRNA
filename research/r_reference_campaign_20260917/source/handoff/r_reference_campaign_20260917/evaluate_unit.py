"""Complete-grid evaluation with frozen semantic rules and explicit abstentions."""
import hashlib,json,os,sys
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
from evaluation_rules import canonical,broad,ptc_general

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    from sklearn.metrics import f1_score,accuracy_score,balanced_accuracy_score,confusion_matrix,precision_recall_fscore_support,roc_auc_score,adjusted_rand_score,normalized_mutual_info_score,fowlkes_mallows_score
    prep=Path(sys.argv[1]);dataset=sys.argv[2]
    dest=prep/'evaluation';dest.mkdir(exist_ok=True)
    units=list(prep.glob('*/score_manifest.json'))
    assert len(units)==4,(str(prep),len(units))
    metrics=[];perclass=[];donors=[];clusters=[];statuses=[];sources={}
    is_ptc=dataset=='PTC'
    if is_ptc:
        reference=pd.read_csv(ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline/paper_baseline_reference.csv.gz',keep_default_na=False).set_index('cell_id')
        tcr=pd.read_csv(ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline/original_DG_binary_pairs_for_R.csv.gz',keep_default_na=False).set_index('cell_id')
        t_names=set(json.loads((ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline/historical_T_names_from_vignette.json').read_text()))
        sys.path.insert(0,str(ROOT/'handoff/ptc_recovery_20260916'))
        from label_rules import strict_T,broad_lineage
        from functools import lru_cache
        strict_T=lru_cache(maxsize=65536)(strict_T)
        broad_lineage=lru_cache(maxsize=65536)(broad_lineage)
    else:
        unit=prep.parent.name
        reference=pd.read_csv(OUT/'inputs'/dataset/unit/'evaluation_only.csv.gz',keep_default_na=False).set_index('cell_id')
        panelmeta=pd.read_csv(OUT/'markers/native_panel_metadata.csv',keep_default_na=False)
        native_names=dict(zip(panelmeta.native,panelmeta.cell_name))
    for manifest in sorted(units):
        source=manifest.parent;m=json.loads(manifest.read_text())
        sources[str(manifest)]=sha(manifest)
        cells=pd.read_csv(source/'cells.csv',keep_default_na=False,dtype=str)
        ids=cells.cell_id.to_numpy()
        assert len(set(ids))==len(ids) and set(ids)<=set(reference.index)
        ref=reference.loc[ids]
        if is_ptc:
            yy=tcr.loc[ids,'truth'].astype(int).to_numpy()
            y_native=ref.paper_final_native.to_numpy(dtype=str)
            y_coarse=np.array([broad_lineage(s) for s in y_native])
            scopes=[('all',np.ones(len(ids),dtype=bool))]
            for col in ['group','sample']:
                for value in ref[col].unique():scopes.append((str(value),(ref[col]==value).to_numpy()))
        else:
            y=ref.truth.map(canonical).to_numpy(dtype=str)
            y_coarse=np.array([broad(s,dataset) for s in ref.truth])
            donor_values=ref.donor.astype(str).to_numpy()
        if (source/'clusters.csv').exists():cl=pd.read_csv(source/'clusters.csv',dtype=str).set_index('cell_id').loc[ids,'cluster'].to_numpy()
        else:cl=cells.cluster.to_numpy()
        cluster_truth=y_coarse if is_ptc else y
        clusters.append(dict(route=source.name,n_cells=len(ids),n_clusters=len(set(cl)),
            evaluation='archived_label_concordance' if is_ptc else 'curated_label_clustering',
            ARI=adjusted_rand_score(cluster_truth,cl),NMI=normalized_mutual_info_score(cluster_truth,cl),
            FMI=fowlkes_mallows_score(cluster_truth,cl)))
        for aid,arm in m['arms'].items():
            td=source/'terminal'/aid
            assert (td/'TERMINAL_COMPLETE').exists(),str(td)
            tm=json.loads((td/'training_manifest.json').read_text())
            assert (td/'TERMINAL_COMPLETE').read_text().strip()==sha(td/'training_manifest.json')
            assert tm['provenance']['score_manifest_sha256']==sha(manifest)
            z=np.load(td/'terminal.npz',allow_pickle=False)
            assert len(z['initial'])==len(ids)
            known=z['initial']!='Undecided'
            assert np.array_equal(z['initial'][known],z['final090'][known])
            statuses.append(dict(route=source.name,arm_id=aid,library=arm['library'],cutoff=arm['cutoff'],
                dl_status=tm['dl_status'],training_executed=tm['training_executed'],n_known=tm['n_known'],
                n_pool=tm['n_pool'],elapsed_seconds=tm['elapsed_seconds']))
            for stage in ['initial','final090','final070']:
                native=z[stage].astype(str)
                unresolved=np.isin(native,['Unknown','Undecided','No_Annotation'])
                common=dict(dataset=dataset,unit=prep.parent.name if not is_ptc else prep.name,route=source.name,
                   arm_id=aid,library=arm['library'],cutoff=arm['cutoff'],stage=stage,dl_status=tm['dl_status'])
                if is_ptc:
                    historical=np.array([ptc_general(v) in t_names for v in native],dtype=int)
                    strict=np.array([strict_T(v) for v in native],dtype=int)
                    concordance=np.array([broad_lineage(v) for v in native])
                    for scope,mask in scopes:
                        for rule,pred in [('paper_broad_T_compatibility',historical),('strict_T_name_rule',strict)]:
                            tn,fp,fn,tp=confusion_matrix(yy[mask],pred[mask],labels=[0,1]).ravel()
                            metrics.append(dict(**common,scope=scope,endpoint=rule,n_cells=int(mask.sum()),
                              F1_nonT=f1_score(yy[mask],pred[mask],pos_label=0,zero_division=0),
                              F1_T=f1_score(yy[mask],pred[mask],pos_label=1,zero_division=0),
                              AUC_binary=roc_auc_score(yy[mask],pred[mask]) if len(set(yy[mask]))==2 else None,
                              accuracy=accuracy_score(yy[mask],pred[mask]),TP=int(tp),FP=int(fp),FN=int(fn),TN=int(tn),
                              unknown_fraction=float(unresolved[mask].mean()),
                              saved_native_concordance=float((native[mask]==y_native[mask]).mean()),
                              saved_broad_lineage_concordance=float((concordance[mask]==y_coarse[mask]).mean())))
                else:
                    semantic=np.array([canonical(native_names.get(v,v)) for v in native])
                    coarse=np.array([broad(native_names.get(v,v),dataset) for v in native])
                    for level,truth,pred in [('curated_semantic',y,semantic),('common_lineage',y_coarse,coarse)]:
                        labels=sorted(set(truth))
                        metrics.append(dict(**common,endpoint=level,n_cells=len(ids),n_truth_classes=len(labels),
                          macro_F1=f1_score(truth,pred,labels=labels,average='macro',zero_division=0),
                          weighted_F1=f1_score(truth,pred,labels=labels,average='weighted',zero_division=0),
                          accuracy=accuracy_score(truth,pred),unknown_fraction=float(unresolved.mean()),
                          predictions_outside_truth_vocabulary=float((~np.isin(pred,labels)).mean())))
                        if level=='curated_semantic' and stage in ['initial','final090']:
                            precision,recall,f1,support=precision_recall_fscore_support(truth,pred,labels=labels,zero_division=0)
                            for lab,pp,rr,ff,nn in zip(labels,precision,recall,f1,support):
                                perclass.append(dict(**common,cell_type=lab,precision=pp,recall=rr,F1=ff,support=int(nn)))
                            all_labels=sorted(set(truth)|set(pred))
                            conf=pd.DataFrame(confusion_matrix(truth,pred,labels=all_labels),index=all_labels,columns=all_labels)
                            cdir=dest/'confusion'/source.name;cdir.mkdir(parents=True,exist_ok=True)
                            conf.to_csv(cdir/(aid+'_'+stage+'.csv.gz'))
                            for donor in sorted(set(donor_values)):
                                mask=donor_values==donor
                                donor_labels=sorted(set(truth[mask]))
                                donors.append(dict(**common,donor=donor,n_cells=int(mask.sum()),
                                  macro_F1=f1_score(truth[mask],pred[mask],labels=donor_labels,average='macro',zero_division=0),
                                  weighted_F1=f1_score(truth[mask],pred[mask],labels=donor_labels,average='weighted',zero_division=0),
                                  accuracy=accuracy_score(truth[mask],pred[mask]),unknown_fraction=float(unresolved[mask].mean())))
    pd.DataFrame(metrics).to_csv(dest/'metrics.csv.gz',index=False)
    pd.DataFrame(statuses).to_csv(dest/'terminal_statuses.csv',index=False)
    pd.DataFrame(clusters).to_csv(dest/'clustering_metrics.csv',index=False)
    if perclass:pd.DataFrame(perclass).to_csv(dest/'per_class.csv.gz',index=False)
    if donors:pd.DataFrame(donors).to_csv(dest/'per_donor.csv.gz',index=False)
    summary=dict(status='complete',dataset=dataset,unit_directory=str(prep),annotation_conditions=len(statuses),
      clustering_conditions=len(clusters),no_refit_or_label_selection_in_evaluation=True,
      evaluation_rule_sha256=sha(Path(__file__).with_name('evaluation_rules.py')),sources=sources,
      job=os.environ['SLURM_JOB_ID'],script_sha256=sha(__file__),
      notes=['Initial = no-DL ablation; final090 = primary terminal endpoint; final070 = confidence sensitivity.',
       'Curated semantic labels and common lineages are distinct resolutions; all Unknowns retained.',
       'PTC TCR-undetected is not established non-T ground truth; saved labels measure concordance.'])
    (dest/'evaluation_manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
    (dest/'COMPLETE').write_text(sha(dest/'evaluation_manifest.json')+'\n')
    print('EVALUATED',dataset,prep,len(statuses),'conditions',flush=True)

if __name__=='__main__':run()
