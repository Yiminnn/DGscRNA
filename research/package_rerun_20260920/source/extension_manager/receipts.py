"""Read-only receipt and checksum validation; never imports scientific libraries."""
from pathlib import Path
from common import need, sha, load, receipt, PARENT

EVALUATOR = '424fb152cee5d95990b0eae724761601199230cb55512e22bd8de5cfcdfd25a9'
NO_TRAIN = {'no_op_all_initially_known', 'no_known_labels_archived_Undecided_terminal',
            'structural_insufficient_known_split'}
PRIMARY, BACKUP = 'CM2_glioma_other', 'CM2_primary_all_context'


def validate(task, release_hash, *, deep=False, pilot=False):
    """Return a complete immutable proof bundle, or raise. No missing-as-success."""
    output = Path(task['output']).resolve()
    validation_root = Path(task.get('verification_output', output)).resolve()
    default_representation = pilot and task['family'] == 'representation' and task.get('pilot_mode') == 'default'
    artifacts = {}

    def remember(path):
        path = Path(path).resolve()
        need(path.is_relative_to(output) or path.is_relative_to(validation_root), f'Acceptance escapes its unit: {path}')
        artifacts[str(path)] = sha(path)
        return artifacts[str(path)]

    def read_receipt(directory, name, flag):
        data = receipt(directory, name, flag)
        remember(Path(directory) / name); remember(Path(directory) / flag)
        return data

    def declared_payload(path, digest):
        path = Path(path).resolve()
        need(path.is_relative_to(output) or path.is_relative_to(validation_root), 'Scientific payload escapes unit')
        need(path.is_file(), 'Scientific payload disappeared')
        if deep:
            need(sha(path) == digest, f'Scientific payload checksum changed: {path}')
        artifacts[str(path)] = digest

    def identity(data):
        need(data.get('status') == ('passed' if default_representation else 'passed_exact') and data.get('sample') == task['sample']
             and data.get('budget') == task['budget'], 'Wrong/incomplete extension receipt')

    def evaluate(parent, parity_name, n):
        parent = Path(parent)
        proof_path = parent / parity_name
        proof = load(proof_path); remember(proof_path); identity(proof)
        need(len(proof['conditions']) == n, 'Parity condition count differs')
        keys = lambda rows: [(r['route'], r['library'], str(r['cutoff'])) for r in rows]
        need(len(set(keys(proof['conditions']))) == n and
             all(r['status'] == 'passed_exact' for r in proof['conditions']), 'Missing/duplicate exact terminal proof')
        evaluation = read_receipt(parent / 'evaluation', 'manifest.json', 'COMPLETE')
        need(evaluation['status'] == 'completed' and evaluation['sample'] == task['sample']
             and evaluation['budget'] == task['budget'] and evaluation['n_conditions'] == n
             and evaluation['n_valid'] == evaluation['n_threshold_rows'] == n * 2,
             'Lfine acceptance incomplete or identifies another unit')
        need(evaluation['script_sha256'] == EVALUATOR
             and evaluation['core_provider_sha256'] == sha(PARENT / 'evaluate_core_lfine.py')
             and evaluation['frozen_semantics'] == load(PARENT / 'evaluate_core_lfine_sources.json'),
             'Lfine semantics changed')
        spec_path = parent / 'evaluation_spec.json'
        spec = load(spec_path)
        need(evaluation['specification_sha256'] == remember(spec_path), 'Evaluation specification changed')
        expected_proof = dict(path=str(proof_path), sha256=sha(proof_path))
        need(evaluation['parity_receipt'] == spec['parity_receipt'] == expected_proof,
             'Lfine does not bind independent parity')
        need(spec['sample'] == task['sample'] and spec['budget'] == task['budget']
             and spec['expected_conditions'] == n and set(keys(spec['conditions'])) == set(keys(proof['conditions'])),
             'Terminal semantic roster differs')
        need(evaluation['outputs']['metrics.csv.gz'] == remember(parent / 'evaluation/metrics.csv.gz'),
             'Lfine metric checksum changed')
        proofs = dict(zip(keys(proof['conditions']), proof['conditions']))
        for condition in spec['conditions']:
            terminal = Path(condition['terminal_directory']).resolve()
            need(terminal.is_relative_to(parent) or (task['family'] == 'A1' and terminal.is_relative_to(output)),
                 'Terminal directory escapes evaluated condition')
            tm = read_receipt(terminal, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
            digest = sha(terminal / 'terminal_manifest.json')
            key = (condition['route'], condition['library'], str(condition['cutoff']))
            need(proofs[key]['actual_manifest_sha256'] == digest and tm['terminal_valid'] is True,
                 'Terminal differs from exact proof')
            need(tm['arm']['library'] == key[1] and str(tm['arm']['cutoff']) == key[2], 'Terminal arm differs')
            need(tm['dl_status'] in NO_TRAIN or tm['dl_status'].startswith('trained'), 'Invalid terminal state')
            need(bool(tm['training_executed']) == tm['dl_status'].startswith('trained'), 'Training state inconsistent')
            declared_payload(terminal / 'terminal.npz', tm['terminal_sha256'])
            if 'training_manifest_sha256' in tm:
                need(remember(terminal / 'training_manifest.json') == tm['training_manifest_sha256'], 'Training manifest changed')
                for filename, digest in load(terminal / 'training_manifest.json')['outputs'].items():
                    declared_payload(terminal / filename, digest)
        for path, wanted in evaluation['terminal_artifacts'].items():
            need(remember(path) == wanted, 'Evaluated terminal metadata/predictions changed')
        return dict(parity_sha256=sha(proof_path), evaluation_manifest_sha256=sha(parent / 'evaluation/manifest.json'),
                    metrics_sha256=evaluation['outputs']['metrics.csv.gz'], terminal_conditions=n, lfine_rows=n * 2)

    family = task['family']
    if family == 'A1':
        space = task['configuration']['space']
        record = read_receipt(validation_root, 'manifest.json', 'COMPLETE'); identity(record)
        spaces = record['spaces']
        need((spaces == [space] or (pilot and space in spaces and len(spaces) == 7))
             and record['n_representations'] == len(spaces)
             and record['n_partitions'] == record['n_terminal_conditions'] == len(spaces) * 13
             and record['n_Lfine_threshold_rows'] == len(spaces) * 26,
             'A1 verification roster differs')
        fit = read_receipt(output, 'fit_manifest.json', 'FIT_COMPLETE')
        need(fit['status'] == 'fit_complete_evaluation_pending' and fit['config']['space'] == space
             and all(fit['config'][k] == task[k] for k in ['sample', 'budget'])
             and len(fit['config']['conditions']) == 13 and fit['config']['reference_labels_used_for_fit'] is False,
             'A1 fit contract differs')
        read_receipt(output, 'representation.json', 'REPRESENTATION_COMPLETE')
        comparison_path = validation_root / 'R_comparison.json'
        need(record['R_comparison_sha256'] == remember(comparison_path), 'A1 R proof changed')
        comparison = load(comparison_path)
        need(len(comparison) == len(spaces) * 13 and all(c['DEG_exact'] and c['density_exact'] for c in comparison),
             'A1 independent R comparison failed')
        proof = record['evaluation_proofs'][space]
        evaluated = evaluate(validation_root / space, 'parity.json', 13)
        for name, actual in [('parity', validation_root / space / 'parity.json'),
                             ('specification', validation_root / space / 'evaluation_spec.json'),
                             ('manifest', validation_root / space / 'evaluation/manifest.json'),
                             ('metrics', validation_root / space / 'evaluation/metrics.csv.gz')]:
            need(proof[name + '_path'] == str(actual) and proof[name + '_sha256'] == sha(actual), 'A1 linked proof differs')
        parity = load(validation_root / space / 'parity.json')
        need(parity['actual_fit_manifest_sha256'] == sha(output / 'fit_manifest.json'), 'A1 parity not bound to current fit')
        for path, digest in record['source_hashes'].items():
            need(sha(path) == digest, 'A1 verification source changed')
        accepted = [evaluated]
    elif family == 'learning':
        accepted = []
        epochs = task.get('pilot_epochs', 30)
        primary = read_receipt(output / PRIMARY, 'acceptance.json', 'LEARNING_VERIFIED_COMPLETE')
        libraries = [PRIMARY] if primary['training_executed'] else [PRIMARY, BACKUP]
        if primary['training_executed']:
            need(not (output / BACKUP).exists(), 'Unpermitted backup marker after successful primary training')
        for library in libraries:
            directory = output / library
            record = primary if library == PRIMARY else read_receipt(directory, 'acceptance.json', 'LEARNING_VERIFIED_COMPLETE')
            identity(record)
            need(record['gate_sha256'] == release_hash and record['epochs'] == epochs
                 and record['model_seed'] == task['configuration']['model_seed'] and record['library'] == library,
                 'Learning configuration differs')
            training = read_receipt(directory, 'training_manifest.json', 'COMPLETE')
            need(record['training_manifest_sha256'] == sha(directory / 'training_manifest.json'), 'Training manifest changed')
            for filename, digest in training['outputs'].items():
                declared_payload(directory / filename, digest)
            for checkpoint in training['checkpoints']:
                cp = directory / f'epoch{checkpoint["epoch"]:02d}'
                declared_payload(cp / 'terminal.npz', checkpoint['terminal_sha256'])
                declared_payload(cp / 'model_state.pt', checkpoint['model_sha256'])
            expected_epochs = [e for e in [5, 10, 20, 30] if e <= epochs] if record['training_executed'] else [0]
            need(record['training_executed'] or record['dl_status'] in NO_TRAIN, 'Invalid no-training endpoint')
            need([c['epoch'] for c in record['conditions']] == expected_epochs
                 and [e['epoch'] for e in record['Lfine']] == expected_epochs, 'Incomplete learning checkpoints')
            for epoch, claimed in zip(expected_epochs, record['Lfine']):
                checkpoint = directory / f'epoch{epoch:02d}' if epoch else directory
                result = evaluate(checkpoint, 'checkpoint_parity.json', 1)
                need(claimed['manifest'] == str(checkpoint / 'evaluation/manifest.json')
                     and claimed['manifest_sha256'] == result['evaluation_manifest_sha256']
                     and claimed['metrics_sha256'] == result['metrics_sha256'] and claimed['valid_threshold_rows'] == 2,
                     'Learning evaluation summary differs')
                accepted.append(result)
            for path, digest in record['source_hashes'].items():
                need(sha(path) == digest, 'Learning source differs from accepted execution')
    else:
        names = dict(geometry=('geometry_acceptance.json', 'GEOMETRY_VERIFIED_COMPLETE', 'geometry_parity.json', 8),
                     mlp=('manifest.json', 'COMPLETE', 'parity.json', 8),
                     neighbors=('acceptance.json', 'NEIGHBOR_VERIFIED_COMPLETE', 'parity.json', 2),
                     no_cluster=('acceptance.json', 'COMPLETE', 'parity.json', 5),
                     representation=('manifest.json', 'COMPLETE', 'parity.json', 2))
        need(family in names, 'Optional extension family is not activated/reviewed')
        name, flag, parity_name, n = names[family]
        if default_representation:
            n = 8
        record = read_receipt(output, name, flag); identity(record)
        need(record['gate_sha256'] == release_hash and record['terminal_conditions'] == n, 'Wrong release/condition count')
        accepted = [evaluate(output, parity_name, n)]
        evaluation = accepted[0]
        parity_hash = record['parity_receipt']['sha256'] if family == 'geometry' else record['parity_sha256']
        need(parity_hash == evaluation['parity_sha256'], 'Acceptance parity checksum differs')
        if family == 'neighbors':
            need(record['configuration'] == task['configuration']['name'] and
                 record['Lfine']['manifest_sha256'] == evaluation['evaluation_manifest_sha256']
                 and record['Lfine']['metrics_sha256'] == evaluation['metrics_sha256'], 'Neighbor acceptance differs')
        else:
            key = 'evaluation_manifest_sha256' if family == 'no_cluster' else 'lfine_manifest_sha256'
            need(record[key] == evaluation['evaluation_manifest_sha256'], 'Lfine acceptance checksum differs')
        if family == 'geometry':
            need(record['pilot_only'] is pilot and record['adapter_sha256'] == sha(PARENT / 'geometry_control.py'),
                 'Geometry pilot/source mismatch')
        if family == 'mlp':
            need(record['control'] == task['configuration']['control'] and record['params'] == task['configuration']['params']
                 and record['adapter_sources_sha256'] == sha(PARENT / 'mlp_controls_sources.json'), 'MLP control differs')
        if family == 'no_cluster':
            need(record['source_lock_sha256'] == sha(PARENT / 'no_cluster_adapter/SOURCE_LOCK.json'), 'A2 sources differ')
        if family == 'representation':
            mode = 'default' if default_representation else ('changed' if pilot else 'task')
            need(record['mode'] == mode and record['source_lock_sha256'] == sha(PARENT / 'representation_sources.json'),
                 'Representation source or mode differs')
            if default_representation:
                configs = record['configurations']
                need(len(configs) == 4 and {(x['space'],x['method']) for x in configs} ==
                     {('PCA30','SNN'),('PCA30','HDBSCAN_R'),('UMAP2','SNN'),('UMAP2','HDBSCAN_R')}
                     and all(x['sample'] == task['sample'] and x['budget'] == task['budget'] and x['minPts'] == 50
                             and x['resolution'] == .5 and x['embedding_seed'] == 42 for x in configs),
                     'Default representation configuration differs')
                need(record['archived_representation_all_artifacts_exact'] is True
                     and record['annotation_and_marker_density_exact_to_core'] is True
                     and record['core_membership_diagnostic_classification'] == 'independently_reproduced_CSV_roundtrip_difference',
                     'Default diagnostic exception is not the independently verified roundtrip contract')
                differences = [x['diagnostic_roundtrip'] for x in record['stage_checks'] if x['diagnostic_roundtrip']]
                audit = load(PARENT / 'representation_sources.json')['files']['roundtrip_audit']
                need(len(differences) == 2 and sha(audit['path']) == audit['sha256'] and
                     all(d['RDS_and_CSV_memberships_each_exact_to_independent_recomputation'] is True
                         and d['causal_audit_sha256'] == audit['sha256'] for d in differences), 'Default causal audit differs')
            else:
                need(record['configurations'] == [task['configuration']], 'Representation configuration differs')
    return dict(status='all_required_extension_receipts_verified', index=task['index'], family=family,
                sample=task['sample'], budget=task['budget'], output=str(output), artifacts=artifacts,
                terminal_conditions=sum(e['terminal_conditions'] for e in accepted),
                lfine_rows=sum(e['lfine_rows'] for e in accepted), evaluations=accepted)
