"""Bounded SLURM-only tests; every scheduler mutation is mocked."""
from pathlib import Path
from unittest import mock
import copy
import json
import os
import sys
import subprocess
import tempfile
import unittest
import common as c
import manage as m
import receipts


class Contracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        c.safe_identity()
        cls.temporary = tempfile.TemporaryDirectory(prefix='extension_manager_test_')
        cls.root = Path(cls.temporary.name)
        base = c.load(c.HERE / 'CAMPAIGN_CANDIDATE_v2.json')
        base['launchers'] = {name:c.record(c.HERE / name) for name in c.LAUNCHERS}
        cls.gate_path = cls.root / 'candidate.json'
        c.write(cls.gate_path, base)
        cls.gate, cls.tasks, cls.release, cls.campaign_root = c.gate_metadata(cls.gate_path, candidate=True)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def state(self):
        return dict(schema=1, campaign_gate_sha256='a' * 64, controller_sha256=c.sha(m.__file__),
                    limits=c.LIMITS, batch_size=64, task_count=743, submissions=[], accepted={})

    def rows(self, running=0, pending=0, job='900'):
        return [dict(job=f'{job}_{i}', state='RUNNING' if i < running else 'PENDING', name='other')
                for i in range(running + pending)]

    def test_queue_reserves_and_running_reserve(self):
        self.assertTrue(m.may_submit(self.rows(190, 500), 64))
        self.assertFalse(m.may_submit(self.rows(224), 64))
        self.assertFalse(m.may_submit(self.rows(190, 582), 64))
        self.assertTrue(m.may_submit(self.rows(223), 64))

    def test_full_roster_profiles_and_noDR_dependencies(self):
        core_tasks = c.load(c.core.checked_file(self.release['task_manifest']))
        default = self.campaign_root / 'extensions/representation_pilots/TKU4163_hvg2000_default_v2/manifest.json'
        tasks = c.task_roster(self.release, core_tasks, full=True, representation_default=default)
        self.assertEqual(len(tasks), 2617)
        self.assertEqual(tasks[923]['family'], 'A1')
        self.assertEqual(tasks[923]['extension_dependencies'], [])
        self.assertEqual(tasks[924]['extension_dependencies'], [923])
        self.assertEqual(tasks[930]['extension_dependencies'], [])
        self.assertEqual(tasks[931]['extension_dependencies'], [930])
        self.assertTrue(all(t['profile'] == 'A1' for t in tasks[923:]))
        self.assertTrue(all(t['profile'] == 'standard' for t in tasks[:923]))
        self.assertTrue(all(Path(t['verification_output']).is_relative_to(self.campaign_root / 'extensions/A1_verification')
                            for t in tasks[923:]))

    def test_worker_and_launcher_smoke(self):
        for name in ['worker.py', 'manage.py', 'build_candidate.py']:
            result = subprocess.run([sys.executable, '-s', str(c.HERE / name), '--help'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('usage:', result.stdout)
        result = subprocess.run(['/bin/bash', '-n', str(c.HERE / 'array.sbatch')], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_adaptive_throttle_bounds(self):
        rows = self.rows(190) + self.rows(16, 48, job='123')
        decision = m.throttle_target(rows, '123')
        self.assertEqual(decision['target'], 49)
        self.assertEqual(decision['unrelated_running'], 190)
        self.assertEqual(m.throttle_target([], '123')['target'], 64)
        blocked = m.throttle_target(self.rows(240), '123')
        self.assertEqual(blocked['target'], 8)
        self.assertTrue(blocked['reserve_temporarily_unavailable'])

    def test_scontrol_only_owned_array(self):
        state = self.state()
        state['submissions'] = [dict(status='submitted', job_id='123', concurrency=16)]
        rows = self.rows(16, 48, job='123') + self.rows(190)
        control = self.root / 'throttle'; control.mkdir()
        with mock.patch.object(m, 'queue', return_value=rows), mock.patch.object(m.subprocess, 'run') as run:
            run.return_value = subprocess.CompletedProcess([], 0, '', '')
            m.adapt_throttle(state, control, rows)
        self.assertEqual(run.call_args.args[0], ['scontrol', 'update', 'JobId=123', 'ArrayTaskThrottle=49'])
        self.assertEqual(state['submissions'][0]['current_throttle'], 49)
        self.assertTrue(state['throttle_decisions'][0]['scientific_parameters_unchanged'])

    def test_crash_after_scontrol_forces_idempotent_recovery(self):
        state = self.state()
        state['submissions'] = [dict(status='submitted', job_id='124', concurrency=16, current_throttle=16)]
        state['throttle_intent'] = dict(array_job='124', target=64, previous=16)
        # Scheduler may actually be at64, while persisted current_throttle is16.
        # New capacity target is also16, so comparing only cached state is wrong.
        rows = self.rows(16, 48, job='124') + self.rows(223)
        control = self.root / 'throttle_recovery'; control.mkdir()
        with mock.patch.object(m, 'queue', return_value=rows), mock.patch.object(m.subprocess, 'run') as run:
            run.return_value = subprocess.CompletedProcess([], 0, '', '')
            m.adapt_throttle(state, control, rows)
        self.assertEqual(run.call_args.args[0][-1], 'ArrayTaskThrottle=16')
        self.assertNotIn('throttle_intent', state)
        self.assertEqual(state['throttle_decisions'][0]['reconciles_unresolved_intent']['target'], 64)

    def test_duplicate_task_registration_rejected(self):
        state = self.state()
        state['submissions'] = [dict(status='submitted', job_id=j, indices=[0], concurrency=16) for j in ['1', '2']]
        with self.assertRaisesRegex(RuntimeError, 'Duplicate scientific'):
            m.registry(state, 'a' * 64, self.tasks, 64)

    def test_state_source_tamper_rejected(self):
        state = self.state(); state['controller_sha256'] = '0' * 64
        with self.assertRaisesRegex(RuntimeError, 'state/options/source'):
            m.registry(state, 'a' * 64, self.tasks, 64)

    def test_finished_job_recovers_lost_sbatch_reply(self):
        state = self.state()
        state['submission_intent'] = dict(name='unique', created_at='2026-09-21T03:00:00+00:00', indices=[0], concurrency=16)
        control = self.root / 'recover'; control.mkdir()
        with mock.patch.object(m, 'accounting', return_value={'555_0': dict(job='555_0', name='unique', state='COMPLETED')}):
            self.assertTrue(m.recover_intent(state, control, []))
        self.assertEqual(state['submissions'][0]['job_id'], '555')
        self.assertNotIn('submission_intent', state)

    def test_ambiguous_sbatch_reply_never_retries(self):
        state = self.state()
        state['submission_intent'] = dict(name='uncertain', created_at=c.utc(), indices=[0], concurrency=16,
                                          explicit_capacity_rejection=False)
        with mock.patch.object(m, 'accounting', return_value={}):
            self.assertFalse(m.recover_intent(state, self.root, []))
        self.assertEqual(state['submissions'], [])
        self.assertIn('submission_intent', state)

    def test_submit_intent_persisted_before_mutation(self):
        state = self.state(); control = self.root / 'submit'; control.mkdir()
        task = dict(self.tasks[0], output=str(self.root / 'fresh_output'))
        gate = dict(self.gate)
        calls = []
        def fake_run(argv, **kwargs):
            saved = c.load(control / 'state.json')
            self.assertEqual(saved['submission_intent']['argv'], argv)
            self.assertTrue(argv[0] == 'sbatch' and any(s == '--array=0%16' for s in argv))
            calls.append(argv)
            return subprocess.CompletedProcess(argv, 0, '998\n', '')
        with mock.patch.object(c, 'gate_metadata'), mock.patch.object(c, 'dependency_proofs', return_value=[{}]), \
             mock.patch.object(m, 'queue', return_value=[]), mock.patch.object(m.subprocess, 'run', side_effect=fake_run):
            self.assertTrue(m.submit(state, control, self.gate_path, gate, [task], self.release, self.root, [0]))
        self.assertEqual(len(calls), 1)
        self.assertEqual(state['submissions'][0]['job_id'], '998')

    def test_generic_qos_failure_is_not_capacity_retry(self):
        state = self.state(); control = self.root / 'qos'; control.mkdir()
        task = dict(self.tasks[0], output=str(self.root / 'fresh_qos'))
        response = subprocess.CompletedProcess([], 1, '', 'Batch job submission failed: Job violates accounting/QOS policy')
        with mock.patch.object(c, 'gate_metadata'), mock.patch.object(c, 'dependency_proofs', return_value=[{}]), \
             mock.patch.object(m, 'queue', return_value=[]), mock.patch.object(m.subprocess, 'run', return_value=response):
            m.submit(state, control, self.gate_path, self.gate, [task], self.release, self.root, [0])
        self.assertFalse(state['submission_intent']['explicit_capacity_rejection'])

    def test_oom_preserved_without_retry(self):
        state = self.state(); state['submissions'] = [dict(status='submitted', job_id='77', indices=[0])]
        row = dict(job='77_0', state='OUT_OF_MEMORY', exit_code='0:125', max_rss='64000M')
        busy, failed = m.failures_and_busy(state, self.tasks, self.root, {'77_0': row}, [])
        self.assertFalse(busy); self.assertEqual(failed[0]['scheduler']['max_rss'], '64000M')
        self.assertFalse(failed[0]['automatic_retry'])

    def test_finished_without_receipt_not_complete(self):
        state = self.state(); state['submissions'] = [dict(status='submitted', job_id='78', indices=[0])]
        busy, failed = m.failures_and_busy(state, self.tasks, self.root,
                                         {'78_0': dict(state='COMPLETED', exit_code='0:0')}, [])
        self.assertFalse(busy); self.assertEqual(len(failed), 1)

    def test_gate_changed_tasks_rejected(self):
        gate = copy.deepcopy(self.gate); gate['families']['geometry'] = 362
        path = self.root / 'tampered_gate.json'; c.write(path, gate)
        with self.assertRaisesRegex(RuntimeError, 'scope or scheduler limits'):
            c.gate_metadata(path, candidate=True)

    def test_candidate_cannot_launch(self):
        with self.assertRaisesRegex(RuntimeError, 'not reviewed'):
            c.gate_metadata(self.gate_path)

    def geometry_task(self):
        task = next(t for t in self.tasks if t['family'] == 'geometry' and t['sample'] == 'TKU4163' and t['budget'] == 'hvg2000')
        return dict(task, output=str(self.campaign_root / 'pilots/geometry_fixed_DL2000_TKU4163_hvg2000'))

    def test_lfine_metric_tamper_rejected(self):
        task = self.geometry_task(); original = receipts.sha
        def tampered(path):
            return '0' * 64 if str(path).endswith('/evaluation/metrics.csv.gz') else original(path)
        with mock.patch.object(receipts, 'sha', side_effect=tampered):
            with self.assertRaisesRegex(RuntimeError, 'metric checksum'):
                receipts.validate(task, self.gate['release_gate']['sha256'], pilot=True)

    def test_terminal_payload_tamper_rejected_deep(self):
        task = self.geometry_task(); original = receipts.sha
        with mock.patch.object(receipts, 'sha', side_effect=lambda p:'0' * 64 if str(p).endswith('/terminal.npz') else original(p)):
            with self.assertRaisesRegex(RuntimeError, 'payload checksum'):
                receipts.validate(task, self.gate['release_gate']['sha256'], deep=True, pilot=True)

    def test_real_A1_single_space_receipts_and_old_reference(self):
        old = c.load(c.PARENT / 'embedding_adapter/tasks.json')[1]
        base = self.campaign_root / 'embedding_adapter_v1'
        verify = base / 'verification_PCA2_early'
        spec = c.load(verify / 'PCA2/evaluation_spec.json')
        actual = Path(spec['conditions'][0]['terminal_directory']).parents[2]
        task = dict(index=924, family='A1', sample=old['sample'], budget=old['budget'], configuration=old,
                    output=str(actual), verification_output=str(verify))
        proof = receipts.validate(task, self.gate['release_gate']['sha256'], deep=True, pilot=True)
        self.assertEqual(proof['terminal_conditions'], 13); self.assertEqual(proof['lfine_rows'], 26)
        self.assertTrue(c.historical_A1_ready(task, deep=True))

    def representation_default(self):
        return dict(index=-1, family='representation', sample='TKU4163', budget='hvg2000', configuration={},
                    output=str(self.campaign_root / 'extensions/representation_pilots/TKU4163_hvg2000_default_v2'),
                    pilot_mode='default')

    def test_representation_default_deep_contract(self):
        proof = receipts.validate(self.representation_default(), self.gate['release_gate']['sha256'], deep=True, pilot=True)
        self.assertEqual(proof['terminal_conditions'], 8)
        self.assertEqual(proof['lfine_rows'], 16)

    def test_representation_default_payload_tamper_rejected(self):
        original = receipts.sha
        with mock.patch.object(receipts, 'sha', side_effect=lambda p:'0' * 64 if str(p).endswith('/terminal.npz') else original(p)):
            with self.assertRaisesRegex(RuntimeError, 'payload checksum'):
                receipts.validate(self.representation_default(), self.gate['release_gate']['sha256'], deep=True, pilot=True)

    def test_cached_receipt_detects_stat_change(self):
        root = self.root / 'stat_test'; unit = c.worker_receipt_root(root, 0); unit.mkdir(parents=True)
        artifact = root / 'data.txt'; artifact.write_text('before')
        value = dict(status='extension_independently_verified', campaign_gate_sha256='a' * 64,
                     release_gate_sha256='b' * 64, worker_sha256=c.sha(c.HERE / 'worker.py'), task_index=0,
                     family='geometry', sample='X', budget='hvg2000', result=dict(artifacts={str(artifact):c.sha(artifact)}))
        c.write(unit / 'acceptance.json', value)
        (unit / 'EXTENSION_VERIFIED_COMPLETE').write_text(c.sha(unit / 'acceptance.json'))
        previous = dict(sha256=c.sha(unit / 'acceptance.json'), artifact_stats=m.artifact_stats([str(artifact)]))
        artifact.write_text('after')
        with self.assertRaisesRegex(RuntimeError, 'artifact stat changed'):
            m.acceptance(root, dict(index=0, family='geometry', sample='X', budget='hvg2000'), 'a' * 64, 'b' * 64, previous)


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Contracts)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    report = dict(status='passed' if result.wasSuccessful() else 'failed', tests=result.testsRun,
                  failures=len(result.failures), errors=len(result.errors), scheduler_mutations=0,
                  source_hashes={str(c.HERE / name):c.sha(c.HERE / name) for name in c.LAUNCHERS},
                  selftest_sha256=c.sha(__file__), job=os.environ['SLURM_JOB_ID'], step=os.environ.get('SLURM_STEP_ID'))
    c.write(c.HERE / 'selftest_result.json', report)
    raise SystemExit(0 if result.wasSuccessful() else 1)
