"""Regression tests for the observed SLURM partial throttle-update response."""
from pathlib import Path
from unittest import mock
import copy
import json
import os
import subprocess
import tempfile
import unittest
import common as c
import manage_v2 as m
import throttle_v2 as t


class Throttle(unittest.TestCase):
    def setUp(self):
        c.safe_identity()
        self.temp = tempfile.TemporaryDirectory(prefix='pkg_throttle_v2_')
        self.control = Path(self.temp.name)
        self.state = dict(submissions=[dict(status='submitted', job_id='123', name='owned',
            concurrency=16, current_throttle=16, indices=[0,1])], accepted={})

    def tearDown(self):
        self.temp.cleanup()

    def rows(self, pending=True, job='123'):
        return [dict(job=job+'_0', state='RUNNING', name='owned'),
                dict(job=job+'_1', state='PENDING' if pending else 'RUNNING', name='owned')]

    def response(self, cap=64, state='RUNNING', job='123', owner='yimin(1)'):
        return subprocess.CompletedProcess([],0,
            f'JobId=456 ArrayJobId={job} ArrayTaskId=0 ArrayTaskThrottle={cap} JobState={state} JobName=owned UserId={owner}\n','')

    def update(self, code=1, error='123_0-1,4: Job has already finished\n'):
        return subprocess.CompletedProcess([],code,'',error)

    def run_case(self, responses, queues=None):
        before = self.rows()
        with mock.patch.object(m,'queue',side_effect=queues or [before,before]), \
             mock.patch.object(t.subprocess,'run',side_effect=responses) as calls:
            m.adapt_throttle(self.state,self.control,before)
        return calls

    def test_skip_without_own_pending(self):
        rows=self.rows(False)+self.rows(True,'999')
        with mock.patch.object(m,'queue',return_value=rows),mock.patch.object(t.subprocess,'run') as call:
            m.adapt_throttle(self.state,self.control,rows)
        call.assert_not_called()
        self.assertEqual(self.state['last_capacity_decision']['action'],'skipped_no_pending_members')

    def test_partial_success_exact_readback(self):
        self.run_case([self.update(),self.response()])
        self.assertEqual(self.state['submissions'][0]['current_throttle'],64)
        self.assertNotIn('throttle_intent',self.state)
        self.assertEqual(self.state['last_capacity_decision']['resolution'],'readback_verified_target')

    def test_mismatched_readback_retains_intent(self):
        with self.assertRaisesRegex(RuntimeError,'differs from requested'):
            self.run_case([self.update(),self.response(16)])
        self.assertIn('throttle_intent',self.state)

    def test_mixed_error_is_not_swallowed(self):
        with self.assertRaisesRegex(RuntimeError,'unrecognized failure'):
            self.run_case([self.update(error='123_0: Job has already finished\nAccess/permission denied\n')])
        self.assertIn('throttle_intent',self.state)

    def test_missing_throttle_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError,'missing its numeric'):
            self.run_case([self.update(),self.response(cap='')])

    def test_wrong_array_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError,'unowned array'):
            self.run_case([self.update(),self.response(job='999')])

    def test_unrecognized_owner_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError,'unowned array'):
            self.run_case([self.update(),self.response(owner='other(2)')])

    def test_all_finish_between_snapshots(self):
        self.run_case([self.update(),self.response(61,'COMPLETED')],queues=[self.rows(),[]])
        self.assertEqual(self.state['last_capacity_decision']['resolution'],'superseded_no_pending_members')
        self.assertFalse(self.state['last_capacity_decision']['successful_throttle_update_claimed'])
        self.assertEqual(self.state['submissions'][0]['current_throttle'],61)

    def test_crash_intent_reissues_even_when_cached_target_matches(self):
        self.state['throttle_intent']=dict(array_job='123',target=64,previous=16)
        rows=self.rows()+[dict(job=f'999_{i}',state='RUNNING',name='other') for i in range(223)]
        calls=self.run_case([self.update(0,''),self.response(16)],queues=[rows,rows])
        self.assertEqual(calls.call_args_list[0].args[0][-1],'ArrayTaskThrottle=16')
        self.assertNotIn('throttle_intent',self.state)

    def test_completed_intent_retired_before_next_wave(self):
        self.state['throttle_intent']=dict(array_job='123',target=61,returncode=1,stderr='123_0: Job has already finished')
        self.state['accepted']={'0':dict(sha='a'),'1':dict(sha='b')}
        history={f'123_{i}':dict(job=f'123_{i}',state='COMPLETED',exit_code='0:0') for i in range(2)}
        t.retire_completed_intent(self.state,self.control,history,[],c)
        self.assertNotIn('throttle_intent',self.state)
        self.state['submissions'].append(dict(status='submitted',job_id='456',name='owned',indices=[2],concurrency=16))
        with mock.patch.object(m,'queue',return_value=self.rows(False,'456')),mock.patch.object(t.subprocess,'run') as call:
            m.adapt_throttle(self.state,self.control,self.rows(False,'456'))
        call.assert_not_called()

    def test_incomplete_accounting_cannot_retire_intent(self):
        self.state['throttle_intent']=dict(array_job='123',target=61)
        with self.assertRaisesRegex(RuntimeError,'successful scientific acceptance'):
            t.retire_completed_intent(self.state,self.control,{},[],c)
        self.assertIn('throttle_intent',self.state)

    def test_migration_preserves_full_state_and_is_idempotent(self):
        import manage as old
        self.state.update(schema=1,campaign_gate_sha256='a'*64,controller_sha256=c.sha(old.__file__),
            limits=c.LIMITS,batch_size=64,task_count=2,submission_intent=dict(id='unchanged-uuid'))
        self.state['accepted']={'0':dict(sha='accepted-original')}
        before=copy.deepcopy(self.state)
        contract=dict(path='reviewed',sha256='b'*64)
        # Simulate a crash after writing the original state snapshot.
        c.write(self.control/'controller_v2_migration_prior.json',self.state)
        with mock.patch.object(m,'controller_contract',return_value=contract):
            m.migrate_scheduler_state(self.state,self.control,'a'*64,[{},{}],64)
            once=copy.deepcopy(self.state)
            m.migrate_scheduler_state(self.state,self.control,'a'*64,[{},{}],64)
        self.assertEqual(self.state,once)
        self.assertEqual(c.load(self.control/'controller_v2_migration_prior.json'),before)
        for key in ('submissions','accepted','submission_intent','campaign_gate_sha256','limits'):
            self.assertEqual(self.state[key],before[key])

    def test_exact_source_derivation(self):
        deriv=c.load(c.HERE/'CONTROLLER_V2_DERIVATION.json')
        text=Path(deriv['derivative']).read_text()
        self.assertEqual(c.sha(deriv['derivative']),deriv['derivative_sha256'])
        for edit in reversed(deriv['literal_replacements']):
            self.assertEqual(text.count(edit['new']),1)
            text=text.replace(edit['new'],edit['old'])
        self.assertEqual(text,Path(deriv['original']).read_text())


if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Throttle))
    record=dict(status='passed' if result.wasSuccessful() else 'failed',tests=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),scheduler_mutations=0,
        source_hashes={str(c.HERE/n):c.sha(c.HERE/n) for n in ['manage_v2.py','throttle_v2.py','selftest_v2.py','CONTROLLER_V2_DERIVATION.json']},
        job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'))
    c.write(c.HERE/'selftest_v2_result.json',record)
    raise SystemExit(0 if result.wasSuccessful() else 1)
