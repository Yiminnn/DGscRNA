"""Bind the four exact completed Lfine pilot proofs; metadata only."""
import common as c
contract=c.contract();proofs={}
for task in contract['pilot_units']:
    ep=c.OUT/'evaluation'/task['sample']/task['budget']/task['space'];vp=c.OUT/'verification'/task['sample']/task['budget']/task['space']
    em=c.checked_outputs(ep);vm=c.checked_outputs(vp)
    assert em['status']=='completed' and em['n_valid']==39 and vm['status']=='passed' and vm['n_unavailable']==0
    proofs[str(vp/'manifest.json')]=c.sha(vp/'manifest.json')
path=c.OUT/'PILOT_ACCEPTANCE.json'
record=dict(status='passed',contract_sha256=c.sha(c.OUT/'contract.json'),proofs=proofs,published_at=c.utc(),scope='Four Lfine replay/count-check pilots only')
if path.exists():assert c.js(path)['proofs']==proofs and c.js(path)['contract_sha256']==record['contract_sha256']
else:c.write(path,record)
print(path)
