# Adversarial review — per-species QUERY / READ path (`taffy lift --bigwig` vector branch)

Scope: `taf_lift.c` `bigwig_lift_window` vector path (`bw->vecN > 0`, ~L1723-1874) +
its use of the libBigWig vector read API. Distrusting the apes-only (N=8, T<2^32,
sparse) validation; hunting N=577 / T=72.5e9 / dense / whole-chromosome failure.

## Verdict
**No blocking or high-severity correctness bug found in the read path for the
intended use** (lift to a LEAF target, dense N-vector, single `uni0` 64-bit axis).
The apes validation is weak, but the logic generalizes: 64-bit types are clean
throughout, the denominator correctly excludes insertion columns, the
coalesce/two-pointer is column-sorted and correct, and `--species` indexing is
verified untransposed. The real at-scale risks are **resource (memory)** and one
**latent inconsistency that bites only the planned sparse format** — plus a
misleading comment. Details below, severity-ranked, confirmed first.

## Confirmed

### C1 (MED, resource) — unbounded `bin_acc` + single-cluster `ov` at fine bins
- `bin_acc = st_calloc(nbins * OUT_N, sizeof(double))` (L1766) is allocated upfront.
  A whole 250 Mb chromosome at `--bin 1000`, all-N (OUT_N=577) =
  250000 × 577 × 8 ≈ **1.15 GB**, before any data is read. `st_calloc` aborts the
  process on OOM (sonLib), so this is a hard failure, not graceful.
- `bwGetOverlappingIntervalsVec(lo,hi)` (L1789) returns one cluster's records ×N
  floats in a single allocation. `COALESCE_GAP=16k` means a dense window collapses
  to few clusters; worst case one cluster spans the window → ~317k records × 577 ×
  4 ≈ **731 MB** in `ov->value` alone.
- Not a correctness bug, and the browser uses COARSE bins (nbins~2000 → ~9 MB), so
  intended use is fine. But the CLI imposes no bound: a naive whole-chrom
  `--bin 1000` all-N query can need ~2 GB and abort on a constrained node.
- Fix sketch: cap/stream bins, or warn when `nbins*OUT_N*8` exceeds a threshold;
  document that whole-chrom all-N wants coarse `-B`.

### C2 (LOW now / HIGH if sparse format ships) — NaN denominator: vector ≠ scalar
- Vector (L1811-1816): the per-component sum is guarded by `if (!isnan(val))`, but
  `cnt += w` (L1816) is **outside** that guard → a NaN-valued record still counts
  its columns in the denominator (diluting coverage toward 0).
- Scalar (L1920): `if (w > 0 && !isnan(value[q])) { sum += ...; cnt += w; }` →
  NaN records are **excluded** from `cnt`.
- So with any NaN in the data, the two paths use different denominators and the
  `Σ_c(vector) == scalar` identity (the whole basis of the parent's correctness
  check) breaks. Masked today because the DENSE writer stores `0.0`, never NaN.
- This is explicitly relevant: the design memo plans "sparse later if ~150 GB
  bites" — a sparse/fill format that surfaces NaN for absent components would make
  this a live bug. Fix: move `cnt += w` inside the same condition the scalar uses
  (or define NaN-handling once, shared).

### C3 (LOW, doc) — stale comment claims the read uses the zoom path
- L1716-1722: "Read via the ZOOM path (bwStatsVec) -- a raw vector fetch would
  pull N floats × every column". The code does exactly the "raw vector fetch" it
  warns against: `bwGetOverlappingIntervalsVec` over RAW base records (L1789).
  `bwStatsVec` is never called on this path. Misleads the at-scale perf story (the
  whole-chrom ~1 s cost is raw base-record I/O, not zoom-summarized). Fix the
  comment (or, if a zoom-coarse read is the intended future, that's the open
  optimization, not what's implemented).

## Suspected / lower

### S1 (LOW-MED, robustness) — `ov == NULL` conflates "no overlap" with "read error"
- L1799 `if (ov != NULL)`: `bwGetOverlappingIntervalsVec` returns NULL for both an
  empty range AND an internal/IO error. On error the cluster silently yields 0
  coverage (`cnt==0`, counted in `n_na`) and the track shows zeros, no error
  surfaced. Same in the scalar path. At 577-way / whole-chrom scale a transient
  read failure becomes a silent partial-zero track. Hard to distinguish without a
  libBigWig error channel; at minimum the large `n_na` could be promoted to a
  warning when it's a high fraction of runs.

### S2 (LOW) — denominator is COVERED columns (`cnt`), not run length (`L`)
- `cnt = Σ w` over overlapping records (L1816), used as the per-run denominator
  (L1821). Correct for a LEAF target: every column of `[cs,ce)` has the target →
  a record exists → `cnt == L`. For an ANCESTOR target, or any target whose
  columns can lack a per-species record (zero-suppressed bins), `cnt < L` and the
  fraction is inflated. The lift doesn't restrict `-g` to leaf genomes. Intended
  use is leaf targets, so low; worth a guard or a doc note.

### S3 (LOW) — uncovered run (`cnt==0`) emits 0, not "no data"
- A run that maps into a bin but finds no record (`cnt==0`) still adds `o` to
  `bin_lv` (L1835), so the bin is emitted as 0 rather than omitted. Observed on
  the apes 500 Mb-SLICE test bw (target regions mapping beyond the built column
  range show as 0). For a FULL 577-way build every leaf column is covered, so this
  won't occur in production — but partial/truncated builds silently read as
  zero-coverage instead of absent.

### S4 (LOW) — minor hardening gaps
- No `ov->N == Nsp` check before `vq = ov->value + q*Nsp` (L1810); always equal in
  practice (same file), but undefended.
- `.names` with MORE than `Nsp` lines is not flagged — only the first `Nsp` are
  read (L1737), so a stray leading/junk line would mislabel every column. Too-few
  names IS caught (L1744). Names > 1023 chars desync (`linebuf[1024]`).
- `--species` on a SCALAR (`vecN==0`) bigWig is silently ignored (falls through to
  the scalar total-depth path); no error.

## Verified clean (coverage notes)
- **64-bit / overflow:** `iv.start/end/t_start` are `int64_t`; `lo/hi` int64 → cast
  `(uint64_t)` for the API; `ov->start/end/l` are `uint64_t` (T=72.5e9 < 2^63, casts
  to int64 safe); value offset `(size_t)q * Nsp`; `bin_acc[bi*OUT_N+oc]` all int64.
  No 32-bit truncation on the read path.
- **Sort / clustering:** `tui_query` qsorts by `start` (tui.c L1840) → the
  coalesce (L1782) and two-pointer (L1800-1817) column-ordering assumption holds
  (the "t-sorted" note in tui.c L680 is a different function). `p` resets per
  cluster, advances monotonically; each interval summed once; records spanning two
  intervals weighted into each; 1kb records can't span a >16k cluster gap.
- **Insertion-column leak:** the denominator/numerator for run `[cs,ce)` come only
  from records intersected with `[cs,ce)` via the two-pointer — the coalesced fetch
  range `[lo,hi)` reads extra records but they are NOT summed into a run unless they
  overlap it, so insertion columns never enter `cnt` or the sum. The memory's
  warning ("merging insertion gaps breaks the denominator") is respected: this path
  uses the BASE tui intervals, not the chained .tui.
- **`--species` indexing:** verified empirically — `--species hg38` = 0.9913 =
  all-N col(hg38); `--species GCA_028885655.2` = 0.9858 = its column; bogus name →
  clean error. `.names` order == component order == all-N column order.
- **Edges:** empty region `a==b` → header only, no crash; region past seq end →
  empty output, no crash; missing `.names` → clean error (the lift self-guards; the
  bench driver's pre-check is belt-and-suspenders), no segfault.
- **`-B` guard:** `bin_size<=0` rejected before the call (L2194) and `-B` parse
  requires `[1,2^60]` → `N` divisor never 0.
- **Reverse strand:** `g_lo/g_hi` from `t_start` (= b-1 for rev) and `L`
  (L1824-1825) matches the scalar path; binned aggregation is orientation-agnostic.

## Bottom line
The read path is correct and 64-bit-clean for the production case (leaf target,
dense 577-vector, full-axis build). Before trusting it at 577/whole-chrom scale:
(1) bound or document the `bin_acc`/`ov` memory for fine-binned all-N queries (C1);
(2) fix the NaN denominator asymmetry before any sparse/NaN-bearing format ships
(C2); (3) correct the stale `bwStatsVec` comment (C3). None of these block the
current dense build + bench.
