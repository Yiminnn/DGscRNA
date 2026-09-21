"""Publish the strict four-unit release gate from existing independent proofs only."""
import common as c
c.contract()
proofs={}
for task in c.pilot_tasks():
    assert c.valid_report(task),task
    path=c.reportpath(task)/'manifest.json';m=c.js(path)
    assert (m['sample'],m['budget'],m['space'])==c.taskkey(task)
    assert m['n_conditions']==13 and m['n_metric_rows']==39 and m['no_model_fits'] is True
    assert m['job'] and m['validation_contract_sha256']==c.sha(c.OUT/'contract.json')
    proofs[str(path)]=c.sha(path)
result=dict(status='passed',validation_contract_sha256=c.sha(c.OUT/'contract.json'),proofs=proofs,
    published_at=c.utc(),scope='Three PCA2 and one ICA2 independent validator pilots only; full A1 remains incomplete')
path=c.OUT/'PILOT_ACCEPTANCE.json'
if path.exists():c.pilot_accepted()
else:c.write(path,result)
print(path)
