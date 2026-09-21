"""Close the root-reviewed full extension gate after all actual pilot proofs."""
from pathlib import Path
import argparse
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'extension_manager'))
import common as c


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    c.safe_identity()
    candidate, tasks, _, _ = c.gate_metadata(args.candidate.resolve(), candidate=True)
    c.need(candidate['status'] == 'candidate_pending_independent_review'
           and candidate['scope'] == 'full_2617' and len(tasks) == 2617,
           'Only the complete independently validated candidate can be reviewed')
    review_path = HERE / 'extension_manager/independent_review.json'
    c.need(c.sha(review_path) == '0a2ba03861daaf312bc13790ec675ffe0f9c5e569d4257d0d886339056124452',
           'Independent review receipt differs from the reviewed source')
    review = c.load(review_path)
    c.need(review['status'] == 'passed_independent_source_review'
           and review['remaining_material_findings'] == [], 'Unresolved independent review findings')
    for path, digest in review['source_hashes'].items():
        c.need(c.sha(path) == digest, 'Reviewed scheduler/worker changed: ' + path)
    test = review['selftest_receipt']
    c.need(c.sha(test['path']) == test['sha256']
           == 'd1722f6c14ec7f139b2151db7173bea4cc5894cffd2c8b4c4d20192fcee5b2de'
           and test['tests'] == 22, 'Required scheduler/receipt regression proof changed')
    c.need({p['role'] for p in candidate['pilots']} == {
        'geometry', 'mlp_original', 'mlp_seed0', 'B_neighbors', 'B_learning', 'A2',
        'representation_default', 'representation_changed', 'A1_full91'}, 'Incomplete pilot coverage')
    c.need(len(candidate['pilot_validation']) == 17, 'Expected complete deep pilot audits')
    target = args.out.resolve()
    c.need(not target.exists(), 'Preserve existing reviewed gate')
    candidate.update(status='reviewed_extension_campaign', reviewed_at=c.utc(),
        root_review=dict(reviewer='/root', candidate=c.record(args.candidate),
            independent_review=c.record(review_path), root_review_source=c.record(__file__),
            actual_full_pilot_receipts_required=True, all_2617_families_included=True,
            public_package_gate_or_installed_wheel_changed=False,
            scientific_campaign_complete=False))
    c.write(target, candidate)
    c.gate_metadata(target)
    print('REVIEWED_EXTENSION_GATE', target, c.sha(target), len(tasks))


if __name__ == '__main__':
    main()
