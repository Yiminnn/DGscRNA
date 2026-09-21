"""Truth-blind adaptive ICA comparator; canonical failures remain explicit.

Native sklearn calls are used without changing the default successful 5000 path.
This is a documented adaptive comparator, not a claim of canonical ICA parity
for the units which require deflation.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import os

POLICY_NAME='ICA2 adaptive (parallel; deflation fallback)'
POLICY_ID='ICA2_parallel_caps_then_deflation_v2'


def policy(core):
    path=core.CAMPAIGN/'protocol/embedding_convergence_repair_20260920_v2.json'
    data=json.loads(path.read_text())
    assert data['policy_id']==POLICY_ID
    assert data['original_protocol_sha256']==core.sha(core.CAMPAIGN/'protocol/embedding.json')
    assert data['attempts']==[
        {'algorithm':'parallel','max_iter':5000},
        {'algorithm':'parallel','max_iter':50000},
        {'algorithm':'parallel','max_iter':100000},
        {'algorithm':'deflation','max_iter':100000},
    ]
    return data,path


def deflation_residuals(model,x):
    """One independent native fixed-point step for each deflated component.

Unit-variance scaling of _unmixing is removed before testing the original
per-component stopping residual. No target labels or fitted annotations enter.
"""
    import numpy as np
    from sklearn.decomposition._fastica import _logcosh, _gs_decorrelation
    white=(model.whitening_ @ (x-model.mean_).T)*np.sqrt(len(x))
    actual=model._unmixing.copy()
    actual/=np.linalg.norm(actual,axis=1)[:,None]
    values=[]
    for j in range(2):
        gx,dg=_logcosh(actual[j]@white,{'alpha':1.0})
        nxt=(white*gx).mean(axis=1)-dg.mean()*actual[j]
        _gs_decorrelation(nxt,actual,j)
        nxt/=np.linalg.norm(nxt)
        values.append(float(abs(abs(nxt@actual[j])-1)))
    return values


def fit(core,cfg,x):
    import numpy as np
    import sklearn
    from sklearn.decomposition import FastICA
    rules,path=policy(core)
    assert sklearn.__version__==rules['fixed_parameters']['sklearn_version']=='1.9.0'
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')
    attempts_dir=Path(cfg['dest'])/'ica_attempts'/(
        str(os.environ.get('SLURM_JOB_ID','unknown'))+'_'+stamp)
    attempts_dir.mkdir(parents=True,exist_ok=False)
    attempts=[]
    for number,item in enumerate(rules['attempts']):
        algorithm,cap=item['algorithm'],item['max_iter']
        # First construction matches the v6 expression exactly, including its
        # pinned defaults. Later retries change only declared algorithm/cap.
        if number==0:
            model=FastICA(n_components=2,max_iter=5000,tol=1e-4,whiten='unit-variance',random_state=42)
        else:
            model=FastICA(n_components=2,max_iter=cap,tol=1e-4,whiten='unit-variance',
                          random_state=42,algorithm=algorithm)
        z,metadata=core.fit_model(model,x)
        metadata['sklearn_version']=sklearn.__version__
        warned=any(w.startswith('ConvergenceWarning:') for w in metadata['warnings'])
        finite=bool(np.isfinite(z).all())
        residuals=None
        if algorithm=='deflation' and finite:
            residuals=deflation_residuals(model,x)
        converged=(not warned and finite and (algorithm!='deflation' or (
            model.n_iter_<cap and np.isfinite(residuals).all() and max(residuals)<1e-4)))
        record=dict(**metadata,attempt=number+1,finite=finite,accepted=bool(converged),
                    fixed_point_residuals=residuals,algorithm=algorithm,cap=cap,
                    input_geometry=str(cfg['geometry']),reference_labels_used=False,
                    producer_source_sha256=core.sha(core.CODE/'run.py'),
                    adaptive_helper_sha256=core.sha(Path(__file__)),
                    policy_sha256=core.sha(path))
        attempt_path=attempts_dir/f'{number+1}_{algorithm}_{cap}.json'
        core.write_json(attempt_path,record)
        attempts.append(dict(path=str(attempt_path),sha256=core.sha(attempt_path),
                             algorithm=algorithm,cap=cap,n_iter=model.n_iter_,
                             accepted=bool(converged),finite=finite,
                             warnings=metadata['warnings'],fixed_point_residuals=residuals))
        if converged:
            metadata.update(comparator_policy=POLICY_NAME,comparator_policy_id=POLICY_ID,
                policy_sha256=core.sha(path),actual_solver=algorithm,actual_iteration_cap=cap,
                canonical_parallel5000_status='converged' if number==0 else 'not_converged',
                fallback_used=algorithm=='deflation',attempts=attempts,
                fixed_point_residuals=residuals,adaptive_helper_sha256=core.sha(Path(__file__)))
            return z,metadata
        if not finite or (algorithm=='parallel' and not warned):
            raise RuntimeError('Invalid ICA output without the declared convergence-warning retry trigger')
    core.write_json(attempts_dir/'ALL_ATTEMPTS_FAILED.json',dict(status='failed',
        attempts=attempts,comparator_policy=POLICY_NAME,policy_sha256=core.sha(path),
        last_iterate_accepted=False,reference_labels_used=False))
    raise RuntimeError('Adaptive ICA exhausted all declared solvers; no invalid iterate accepted')
