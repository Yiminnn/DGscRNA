# A1 adaptive ICA implementation v8

This version defines a truth-blind **adaptive ICA comparator**: native sklearn
parallel FastICA at 5,000, 50,000 and 100,000 iterations, followed only on
nonconvergence by native deflation FastICA at 100,000 iterations. Fixed settings
remain 2 components, seed 42, tolerance 1e-4, unit-variance/SVD whitening and the
default logcosh contrast. The first finite converged result is accepted.

The first parallel call is unchanged. Deflation is a different solver variant,
not an assertion of canonical parallel equivalence. Its convergence also requires
n_iter < cap and an independent per-component fixed-point residual below 1e-4.
Failed attempts remain in immutable, timestamped attempt files. No author labels,
annotation metrics, clusters or refinement outputs are read for solver choice.

Only successful v5/v6 caches with explicitly pinned producers may be consumed.
Their producing source is retained. Separate consumer receipts and immutable
producer-manifest snapshots describe reuse. Native R export/HDBSCAN/scoring and
DL/refinement files remain byte-identical. `embedding.json` is not modified.

This snapshot does not authorize production activation. Root reviews the v2
protocol, successful ICA regression, native failing-case recovery and source
compatibility before replacing a dispatcher. Existing v5 pilot jobs are untouched.

Presentation changes are reviewed separately. Existing figures referenced by
notebooks/receipts must not be overwritten: new adaptive-labelled figures need a
separate version/path. Main summaries must call every ICA row an adaptive policy,
with actual solver/cap/attempt fields; canonical failures remain visible as
technical failures. The copied selector is the frozen v6 baseline; the separately
hardened final selector must be integrated by root before final cohort release.

V8 adds only cache-integrity validation: actual embedding CSV hash and cell order,
all seven score-input fingerprint fields against current prepared inputs, and
terminal prediction CSV hashes. Completed 13-candidate rosters are explicit.
Scoring, DL and adaptive ICA arithmetic are unchanged from v7; its successful
regression is preserved and hash-linked. No v7 production outputs were produced.
