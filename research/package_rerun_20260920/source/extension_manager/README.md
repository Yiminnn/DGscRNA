# Frozen extension campaign

The unified roster has 2,617 units: geometry 363, MLP 54, Reviewer B neighbors
66 and learning 18, no-cluster 242, representation 180, and A1 1,694. A full
campaign gate is built only after all source locks and complete pilot proofs
pass. The earlier 743-unit candidates are nonlaunchable review artifacts.

The existing release gate, installed package, core launchers and accepted core
directories remain unchanged. Every worker requires its own current
`PACKAGE_VERIFIED_COMPLETE` core dependencies. Geometry additionally requires
the same sample's HVG2000 core. A1 additionally waits for all 13 historical
reference terminals and for that sample/budget's new noDR unit before any of
its other six spaces can run; noDR therefore owns the shared geometry export.

## Scheduling and recovery

- One active array contains at most 64 tasks, all using the same resource
  profile. Standard controls use 4 CPUs / 64 GB / 8 hours; A1 uses 4 CPUs /
  32 GB / 24 hours. No scheduler requeue is requested.
- Initial throttle is 16. Every 60 seconds, the controller may update only its
  own array's throttle to 8–64 using the expanded all-user queue. Admission
  leaves 64 queue slots below 900 and targets 16 running slots for core below
  the global 256 ceiling. Unrelated controllers can consume that reserve;
  status explicitly reports when it is unavailable. No unrelated job is changed.
- Each `sbatch` has a durable UUID/name intent before submission. A lost reply
  is reconciled through both queue and accounting, including finished jobs.
  Only explicit submit-count rejection can be retried. Timeout or generic QOS
  failures remain unresolved until reviewed. An unresolved throttle intent is
  recovered by an idempotent update to the newly calculated target.
- State and status are atomic JSON. A directory lock prevents concurrent
  controllers. Inspect its owner/job before manually resolving a stale lock.
  A `STOP` file stops new scheduling without cancelling any jobs.
- Failed/OOM/scientifically invalid units preserve their outputs, worker phase,
  logs and `sacct` state/exit/RSS records. No automatic scientific retry or
  parameter change occurs. A worker's durable unit directory prevents a
  duplicate/requeued process from fitting again.

## Acceptance

Workers run the frozen adapter, independently validate all required parity and
Lfine receipts, then record `EXTENSION_VERIFIED_COMPLETE`. Learning acceptance
includes every required checkpoint and the backup marker only after a legal
primary no-training state. Every accepted result binds its actual array job,
task, core dependencies, source gate and scientific payload checksums.

A1 fit directories are `extensions/A1/<sample>/<budget>/<space>`. Their frozen
task `verification_output` points to
`extensions/A1_verification/<sample>/<budget>/<space>`; only these two explicit
roots may contain A1 acceptance artifacts. The verifier's receipt, not
`FIT_COMPLETE` alone, establishes A1 completion.

Workers record file size/mtime/ctime/inode snapshots after their full checksum
audit. The controller checks this unchanged snapshot and the full receipt chain
at first acceptance. Later polls compare the small acceptance hash and artifact
snapshots; a changed snapshot stops for review. At campaign completion it
performs one complete payload checksum audit. `GBM_EXTENSIONS_COMPLETE` is
written only when all 2,617 receipts pass and every registered array task ended
successfully. There is no automatic PTC launch.

## Build and launch

`build_candidate.py --full` requires the unchanged release gate, the final
representation default and changed manifests, and the full A1 91-terminal
manifest. It creates a candidate with exact launcher/source/pilot/task hashes.
The root reviews the candidate and tests before creating the actual gate with
status `reviewed_extension_campaign`. Candidate status is rejected by both
controller and worker. Never edit frozen files to extend a running campaign;
use separately reviewed versioned paths and control namespaces if a later
campaign is needed.

Only the root launches `manage.py --campaign <actual gate> --watch` in SLURM.
`--once` performs one scheduling iteration; it is not a dry-run. Use
`selftest.py` for mocked scheduler and tamper tests. All scientific operations,
including deep checksum audits and tests, run through SLURM.
