"""Additional independent whole-unit checks; no model creation or prediction."""
from collections import Counter
from pathlib import Path
import json,os
import common as c

def manual_classes(y,p,labels):
    import numpy as np
    out=[]
    for label in labels:
        tp=int(((y==label)&(p==label)).sum());fp=int(((y!=label)&(p==label)).sum());fn=int(((y==label)&(p!=label)).sum())
        out.append(dict(label=label,precision=tp/(tp+fp) if tp+fp else 0.,recall=tp/(tp+fn) if tp+fn else 0.,
            F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.,support=tp+fn))
    return out

def partition_scores(a,b):
    """Contingency-count formulas, independent of sklearn metric calls."""
    import numpy as np
    n=len(a);pairs=Counter(zip(a,b));ra=Counter(a);rb=Counter(b)
    if n<2:return dict(ARI=1.,NMI=1.,FMI=0.)
    choose=lambda x:x*(x-1)/2
    joint=sum(choose(v) for v in pairs.values());aa=sum(choose(v) for v in ra.values());bb=sum(choose(v) for v in rb.values())
    expected=aa*bb/choose(n);denom=(aa+bb)/2-expected
    ari=(joint-expected)/denom if denom else 1.
    mi=max(0.,sum(v/n*np.log(v*n/(ra[x]*rb[y])) for (x,y),v in pairs.items()))
    ha=-sum(v/n*np.log(v/n) for v in ra.values());hb=-sum(v/n*np.log(v/n) for v in rb.values())
    nmi=mi/((ha+hb)/2) if ha+hb else 1.
    fmi=(joint/aa*joint/bb)**.5 if aa and bb else 0.
    return dict(ARI=ari,NMI=nmi,FMI=fmi)

def close(a,b,field):
    import numpy as np
    if a is None:assert b is None or np.isnan(float(b)),field
    else:assert np.isclose(float(a),float(b),rtol=0,atol=1e-12,equal_nan=True),(field,a,b)

def check_figures(directory,conditions):
    from PIL import Image
    m=c.js(directory/'manifest.json');assert c.checked(directory)
    assert c.sha(m['display'])==m['display_sha256']
    expected={name+'.'+ext for name in c.ROUTES for ext in ['png','pdf']}
    assert set(m['files'])==expected and len(m['conditions'])==13 and m['every_candidate_plotted'] is True
    assert {Path(r['condition']['dest']).name for r in m['conditions']}==set(c.ROUTES)
    expected_conditions={item['route']:item for item in conditions}
    assert {item['condition']['route']:item['condition'] for item in m['conditions']}==expected_conditions
    for record in m['conditions']:
        route=Path(record['condition']['dest']);tm=route/'terminal/L00_mean/terminal_manifest.json'
        assert set(record['files'])=={route.name+'.png',route.name+'.pdf'}
        assert c.sha(tm)==record['terminal_manifest_sha256']
        t=c.js(tm);assert record['training_executed']==t['training_executed'] and record['dl_status']==t['dl_status']
    for name,digest in m['files'].items():
        p=directory/name;assert c.sha(p)==digest and p.stat().st_size>100
        if p.suffix=='.png':
            with Image.open(p) as image:assert min(image.size)>200;image.verify()
        else:
            with p.open('rb') as f:assert f.read(5)==b'%PDF-'
    return m

def solver_audit(unit,cfg,rep,contract):
    allowed=contract['allowed_run_producers'];assert rep['source_sha256'] in allowed
    assert cfg['source_sha256'] in allowed and c.js(unit/'fit_manifest.json')['source_sha256'] in allowed
    assert rep['reference_labels_used_for_fit'] is False
    if cfg['space']!='ICA2':
        expected={'noDR':{},'PCA2':dict(n_components=2,svd_solver='randomized',iterated_power=7,random_state=42),
            'FA2':dict(n_components=2,max_iter=5000,tol=.01,random_state=42),
            'Isomap2':dict(n_components=2,n_neighbors=15,n_jobs=1),
            'TSNE2':dict(n_components=2,perplexity=30,init='pca',learning_rate='auto',max_iter=1000,random_state=42,n_jobs=1)}
        if cfg['space']=='UMAP2':
            for k,v in dict(seed=42,n_neighbors=30,n_components=2,metric='cosine',min_dist=.3,threads=1).items():assert rep[k]==v,(k,rep[k])
        else:
            for k,v in expected[cfg['space']].items():assert rep['params'][k]==v,(k,rep['params'][k])
        assert not any(w.startswith('ConvergenceWarning:') for w in rep.get('warnings',[]))
        return dict(space=cfg['space'],producer=rep['source_sha256'])
    import numpy as np
    policy=c.js(c.POLICY);params=rep['params']
    for k,v in dict(n_components=2,tol=1e-4,random_state=42,whiten='unit-variance',whiten_solver='svd',fun='logcosh',fun_args=None,w_init=None).items():assert params[k]==v,(k,params[k])
    assert rep['sklearn_version']=='1.9.0'
    assert not any(w.startswith('ConvergenceWarning:') for w in rep['warnings'])
    oldfail=unit/'CONVERGENCE_FAILURE.json'
    if rep['source_sha256']==contract['v8_run_sha256']:
        attempts=rep['attempts'];assert rep['policy_sha256']==c.sha(c.POLICY)
        assert rep['comparator_policy_id']==policy['policy_id'] and rep['adaptive_helper_sha256']==contract['adaptive_helper_sha256']
        assert len(attempts)>=1 and len(attempts)<=4
        assert [(a['algorithm'],a['cap']) for a in attempts]==[(p['algorithm'],p['max_iter']) for p in policy['attempts'][:len(attempts)]]
        assert [a['accepted'] for a in attempts]==[False]*(len(attempts)-1)+[True]
        for a in attempts:
            assert c.sha(a['path'])==a['sha256'];m=c.js(a['path'])
            assert m['reference_labels_used'] is False and m['policy_sha256']==c.sha(c.POLICY)
            assert m['producer_source_sha256']==contract['v8_run_sha256'] and m['input_geometry']==cfg['geometry']
            assert m['adaptive_helper_sha256']==contract['adaptive_helper_sha256'] and m['sklearn_version']=='1.9.0'
            assert m['algorithm']==a['algorithm'] and m['cap']==a['cap'] and m['accepted']==a['accepted'] and m['finite'] is True
            assert m['n_iter_']==a['n_iter'] and m['warnings']==a['warnings'] and m['fixed_point_residuals']==a['fixed_point_residuals']
            for key,value in params.items():
                if key not in ['algorithm','max_iter']:assert m['params'][key]==value,key
            assert m['params']['algorithm']==a['algorithm'] and m['params']['max_iter']==a['cap']
            if not a['accepted']:assert a['algorithm']=='parallel' and any(w.startswith('ConvergenceWarning:') for w in a['warnings'])
        canonical='converged' if len(attempts)==1 else 'not_converged'
        assert rep['canonical_parallel5000_status']==canonical
        assert rep['actual_solver']==attempts[-1]['algorithm']==params['algorithm']
        assert rep['actual_iteration_cap']==attempts[-1]['cap']==params['max_iter']
        assert rep['n_iter_']==attempts[-1]['n_iter'] and rep['warnings']==attempts[-1]['warnings']
        assert rep['fixed_point_residuals']==attempts[-1]['fixed_point_residuals']
        assert rep['fallback_used']==(rep['actual_solver']=='deflation')
        if rep['fallback_used']:
            assert rep['n_iter_']<rep['actual_iteration_cap']
            residuals=rep['fixed_point_residuals'];assert len(residuals)==2 and np.isfinite(residuals).all() and max(residuals)<1e-4
    else:
        assert params['algorithm']=='parallel' and params['max_iter']==5000 and 0<rep['n_iter_']<=5000
        attempts=[];canonical='converged'
    return dict(sample=cfg['sample'],budget=cfg['budget'],space='ICA2',space_display='ICA2 adaptive',
        actual_solver=rep.get('actual_solver',params['algorithm']),actual_iteration_cap=rep.get('actual_iteration_cap',params['max_iter']),
        canonical_parallel5000_status=canonical,fallback_used=rep.get('fallback_used',False),n_iter=rep['n_iter_'],
        producer_source_sha256=rep['source_sha256'],policy_sha256=c.sha(c.POLICY),representation_manifest_sha256=c.sha(unit/'representation.json'),
        adaptive_figure_manifest_sha256=c.sha(unit/'figures_adaptive_v1/manifest.json'),attempts_json=json.dumps(attempts,sort_keys=True),
        old_failure_path=str(oldfail) if oldfail.exists() else None,old_failure_sha256=c.sha(oldfail) if oldfail.exists() else None)

def terminal_state(td,initial,pred,tm,cells,contract):
    import numpy as np
    import pandas as pd
    train=c.js(td/'training_manifest.json');assert c.sha(td/'training_manifest.json')==tm['training_manifest_sha256']
    assert train['terminal_valid'] is True and tm['terminal_valid'] is True
    assert train['dl_status']==tm['dl_status'] and train['training_executed']==tm['training_executed']
    for name,digest in train['outputs'].items():assert c.sha(td/name)==digest,name
    sig=train['provenance']['input_signature'];assert sig['training_source_sha256']==contract['refinement_source_sha256']
    assert train['params']==sig['params']==contract['DL_params']
    assert train['provenance']['driver_source_sha256']==contract['terminal_driver_sha256']
    assert train['provenance']['reference_labels_used_for_fit'] is False
    assert sig['DL_sha256']==tm['DL_sha256'] and sig['n_cells']==len(cells)
    assert sig['n_features']==tm['DL_features']==train['input_width']
    assert sig['torch_threads']==4 and sig['torch_version']=='2.5.1'
    assert sig['cell_order_sha256']==__import__('hashlib').sha256('\n'.join(cells).encode()).hexdigest()
    with np.load(td/'terminal.npz',allow_pickle=False) as z:
        init=np.asarray(initial,dtype=str);known=np.flatnonzero(init!='Undecided');pool=np.flatnonzero(init=='Undecided')
        assert np.array_equal(z['initial'],init) and np.array_equal(z['pool_indices'],pool)
        assert np.array_equal(z['classes'],np.asarray(sorted(set(init))))
        assert train['classes']==z['classes'].tolist()
        assert len(known)==tm['n_known']==train['n_known'] and len(pool)==tm['n_pool']==train['n_pool']
        assert len(set(init[known]))==tm['n_training_classes']==train['n_training_classes']
        signature_bytes=json.dumps(sig,sort_keys=True).encode()
        cache_key=__import__('hashlib').sha256(signature_bytes+b'\0'+json.dumps(init.tolist(),ensure_ascii=False).encode()).hexdigest()
        assert cache_key==tm['cache_key']==train['provenance']['cache_key']
        for key in ['initial','final090','final070','lineage']:assert np.array_equal(z[key],pred[key])
        for key in ['final090','final070']:assert np.array_equal(z[key][known],init[known])
        for key in ['final090','final070']:assert len(z[key])==len(cells)
        confidence=pd.to_numeric(pred.confidence_rounded,errors='coerce').to_numpy()
        np.testing.assert_allclose(confidence,z['confidence_rounded'],rtol=0,atol=1e-7,equal_nan=True)
        assert np.isnan(z['confidence_rounded'][known]).all()
        assert np.array_equal(z['predicted'][known],np.full(len(known),''))
        assert np.array_equal(z['lineage'][known],np.full(len(known),'retained_initial'))
        if len(pool)==0:expected='no_op_all_initially_known'
        elif len(known)==0:expected='no_known_labels_archived_Undecided_terminal'
        elif len(known)<2:expected='structural_insufficient_known_split'
        else:expected='trained_single_known_class' if len(set(init[known]))==1 else 'trained'
        assert tm['dl_status']==expected
        if tm['training_executed']:
            assert expected in ['trained','trained_single_known_class'] and train['prediction_executed'] is True
            probs=z['probabilities'];assert probs.shape==(len(pool),len(z['classes'])) and np.isfinite(probs).all()
            assert (probs>=0).all() and (probs<=1).all();np.testing.assert_allclose(probs.sum(1),1,rtol=1e-5,atol=1e-6)
            rounded=np.asarray([round(v,4) for v in probs.max(1)],dtype=np.float32)
            calls=z['classes'][probs.argmax(1)]
            assert np.array_equal(z['confidence_rounded'][pool],rounded) and np.array_equal(z['predicted'][pool],calls)
            for key,threshold in [('final090',.9),('final070',.7)]:assert np.array_equal(z[key][pool],np.where(rounded>=threshold,calls,'Unknown'))
            assert np.array_equal(z['lineage'][pool],np.where(rounded>=.9,'DL_accepted','DL_low_confidence'))
            import torch
            assert torch.__version__==sig['torch_version']
            order=torch.randperm(len(known),generator=torch.Generator().manual_seed(42)).numpy();ntrain=round(len(known)*.9)
            assert np.array_equal(z['train_indices'],known[order[:ntrain]]) and np.array_equal(z['validation_indices'],known[order[ntrain:]])
            assert train['n_train']==ntrain and train['n_validation']==len(known)-ntrain
            history=c.js(td/'training_history.json');assert len(history)==10 and [h['epoch'] for h in history]==list(range(1,11))
            assert all(h['n']==ntrain and np.isfinite(h['mean_loss']) and h['mean_loss']>=0 and 0<=h['accuracy']<=1 for h in history)
        else:
            assert train['prediction_executed'] is False and z['probabilities'].shape==(0,len(z['classes']))
            assert len(z['train_indices'])==len(z['validation_indices'])==0
            assert np.isnan(z['confidence_rounded']).all() and np.all(z['predicted']=='')
            assert np.array_equal(z['lineage'][pool],np.full(len(pool),'unresolved'))
            for key in ['final090','final070']:
                expected_final=init.copy()
                if expected=='structural_insufficient_known_split':expected_final[pool]='Unknown'
                assert np.array_equal(z[key],expected_final)
        assert tm['n_final_called090']==int((~np.isin(z['final090'],['Unknown','Undecided'])).sum())
        assert tm['n_final_called070']==int((~np.isin(z['final070'],['Unknown','Undecided'])).sum())
    return dict(dl_status=expected,training_executed=tm['training_executed'],n_known=len(known),n_pool=len(pool))
