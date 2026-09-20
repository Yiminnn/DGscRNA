# Resource ledger schema and boundaries

`registry_inventory.csv` lists the finite top-level registries and hashes examined; no recursive cache/results scan. `job_provenance.csv` links explicit job IDs to recorded commands or source metadata. `accounting_queries.json` records own-user-only `sacct --user=yimin --duplicates --array` queries. Array task and retry accounting are retained, with DBIndex/SLUID disambiguation. `accounting_raw.csv` preserves raw allocation and step records.

`job_attempt_ledger.csv` uses allocation records as cost units, or individually recorded steps in a shared interactive allocation. It does not add batch/extern step CPU-hours on top of the parent allocation. `allocated_cpu_hours = ElapsedRaw × AllocCPUS /3600`; measuredCPU derives from recorded TotalCPU, with UserCPU/SystemCPU retained. These quantities are different. Actual scheduler AllocCPUS may exceed the program's requested threads for memory scheduling. Failed/requeued attempts remain visible. Bare shared7204040 costs are excluded; only explicitly identified steps can be attributed.

`batch_maxrss_kib` and `max_step_maxrss_kib` are reported maxima, not aggregate fork/process peak. Never multiply or sum them to manufacture a true simultaneous peak. `phase_cost_summary.csv` keeps combined pipeline jobs unsplit unless source timers prove separation. Sum of job wall times is separate from calendar span. Running/pending attempts are provisional.

`terminal_execution_evidence.json` keeps requested condition counts, inherited training flags, fresh training, cache reuse, and unique models distinct. Missing cache uniqueness stays null. The old public/PTC grid is not counted as a fresh followup fit. `separate_scaling_evidence.json` records the45 independent CPU scaling runs, maximum120k, separately from complete grid-search cost. GPU memory isN/A for the CPU implementation.

The first ledger is deliberately partial while the new campaign runs and historical cache/combined-stage coverage is unresolved. A green computation job is not closure of reviewerF cost coverage.

Phase labels in this first inventory are provisional keyword classifications from explicit submission commands or recorded source paths. They are not instrumented stage timers, and do not close the prepare/integration/embedding/DEG/scoring/DL cost split. Family-level recorded accounting totals are independent of this provisional phase classification. Scheduler timestamp strings retain the cluster's local timezone; collection timestamps are explicitlyUTC.
