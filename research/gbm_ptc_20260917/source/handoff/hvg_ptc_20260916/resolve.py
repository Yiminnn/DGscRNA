"""Verify structurally unavailable conditions without changing frozen fits or annotations."""
import json
from common import OUT, require_slurm, sha, utc, write_json


def structural_resolution(sample, geometry, arm):
    require_slurm()
    import numpy as np
    prep=OUT/'prepared'/sample
    pm=json.loads((prep/'manifest.json').read_text())
    if not (prep/'PREPARED').exists():
        return None
    assert (prep/'PREPARED').read_text().strip()==sha(prep/'manifest.json')
    g,aid=geometry['geometry_id'],arm['arm_id']
    dest=OUT/'fits'/sample/g/aid
    resolution=dict(sample=sample,geometry_id=g,arm_id=aid,timestamp=utc(),
                    final_valid=False,terminal_valid=False,scientific_parameters_changed=False,
                    original_error_preserved=True,prepared_manifest_sha256=sha(prep/'manifest.json'))
    unavailable=pm.get('unavailable_features',{}).get(geometry['feature'])
    if unavailable:
        assert unavailable['stage']=='VST_loess' and unavailable['status']=='structural_failure'
        resolution.update(status='structural_preprocessing_failure',reason='VST_LOESS_near_singularity',
                          evidence=unavailable,partition_available=False)
    elif (dest/'manifest.json').exists() and not (dest/'COMPLETE').exists():
        am=json.loads((dest/'manifest.json').read_text())
        if am.get('status')!='failed':
            return None
        error=am.get('error','')
        if ('Could not calculate statistics for groups' in error and 'only contain one sample' in error
                and (dest/'clusters.npy').exists()):
            cl=np.load(dest/'clusters.npy',allow_pickle=False)
            groups,counts=np.unique(cl[cl!=-1],return_counts=True)
            singleton=groups[counts==1]
            assert len(singleton)>0, 'Singleton failure without an actual singleton cluster'
            assert len(cl)==pm['n_cells']
            resolution.update(status='structural_annotation_failure',reason='legacy_DEG_cannot_score_singleton_clusters',
                singleton_cluster_ids=singleton.tolist(),partition_available=True,
                cluster_sha256=sha(dest/'clusters.npy'),original_fit_manifest_sha256=sha(dest/'manifest.json'),error=error)
        elif 'ill-defined empirical covariance' in error:
            evidence=OUT/'numerical_failure_checks'/sample/g/aid/'verification.json'
            if not evidence.exists():return None
            check=json.loads(evidence.read_text())
            assert check['original_fit_manifest_sha256']==sha(dest/'manifest.json')
            assert check['embedding_sha256']==sha(dest.parent/'embedding.npy')
            assert check['parameters_changed'] is False and len(check['attempts'])>=3
            if not check['reproducible_covariance_failure']:return None
            resolution.update(status='numerical_clustering_failure',reason='GMM_covariance_failure_reproduced_without_parameter_changes',
                partition_available=False,error=error,evidence_path=str(evidence),evidence_sha256=sha(evidence),
                original_fit_manifest_sha256=sha(dest/'manifest.json'))
        else:
            return None
    else:
        return None
    out=OUT/'structural_resolutions'/sample/g/f'{aid}.json'
    write_json(out,resolution)
    return resolution
