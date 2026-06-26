# Adversarial review — per-species DEPTH WRITE path (at N=577 / T>2³² / sharded / high -T)

Scope: `taf_depth.c` (`--perSpecies --bigwig`: UniBW, gerp_flush_bin[_vec], score_one_block,
phase-C binner, teardown), `taf_merge_bigwig.c`, and taffy's use of the libBigWig vector
WRITE API (libBigWig internals reviewed separately in REVIEW_vector_format.md; spot-checked
here only where taffy's usage could trip a latent 64-bit bug).

Method: read-only code trace. The existing apes validation (N=8, T<2³², small, sparse) was
treated as proving nothing about the real target; every finding below is reasoned at
N=577 / T=72.5e9 / sharded build / -T 32.

**Bottom line: no BLOCKING or HIGH bug found after genuinely trying to break it.** Two MED
items (one correctness-of-an-auxiliary-output at non-default `--minLeaves`, one memory
quantification) and two LOW/informational. The riskiest things at scale — cross-shard
component alignment and 64-bit coordinate handling — are sound.

---

## CONFIRMED (code trace)

### M1 [MED] `--depth` text alongside `--perSpecies` diverges from the standalone scalar `--depth` when `--minLeaves > 1`
`gerp_flush_bin_vec` (taf_depth.c:287-288) writes the `--depth` bedGraph value as
`vsum/cnt` where `vsum = Σ_c sumvec[c]` — the **ungated** total leaf-presence count
(every leaf present at every column, accumulated at L431-433, explicitly NOT min_leaves-gated).
The standalone scalar track (`gerp_flush_bin`, L257-260) emits `bin_sum/cnt` where `bin_sum`
sums `depth` only for **scored** columns (`scored = depth >= min_leaves`, L425).
- At `--minLeaves 1` (the build/bench default) the two are identical (a depth-0 column adds
  nothing to either; a depth-d≥1 column adds d to both) — so the Σ_c==scalar gate held, and
  there is **no impact on the actual 577 build**.
- At `--minLeaves ≥ 2` they differ: a depth-1 column is unscored (adds 0 to `bin_sum`) but its
  one present leaf adds 1 to `vsum`. So the `--depth`-beside-`--perSpecies` mean reads higher.
- The inline comment L274-275 ("= sum of the per-leaf counts / cnt, matching the scalar depth
  track") over-claims; it only matches at minLeaves=1.
Fix: either gate `vsum`/emit `bin_sum` for the `--depth` sibling, or soften the comment to
"matches the scalar depth track only at --minLeaves 1". The per-species **bigWig** itself is
unaffected (counts are ungated by design; the lift divides by full span).

### M2 [MED] `bin_sum_vec` peak memory is N×=577 the scalar and persists at the high-water mark
`results[i].bin_sum_vec` is `bin_cap × n_leaves × int64` per batch slot (L414-415), grown
geometrically and **never shrunk** for the run. `bin_cap` reaches the largest single block's
bin span any slot ever saw. With `batch_cap = 4·n_threads`:
- peak ≈ `4·n_threads × (max_block_columns / bin_size) × N × 8` bytes.
- At N=577, -T 32 (batch_cap=128), bin=1000: a 1 Mb single block → ~590 MB; a 10 Mb block → ~5.9 GB.
Plus a fixed `UniBW.values` buffer = `UNIBW_BATCH(65536) × N × 4` = **151 MB at N=577** (L200),
and the merge loads a whole shard via `bwGetOverlappingIntervalsVec` (~232 MB for the default
`-s 1e8` at N=577; scales with `-s`).
Not a leak and not a correctness bug — but the SLURM `--mem` must cover it. The 24 G/shard
default is comfortable for typical TAF block sizes; the risk is a pathological very-long single
block at high -T. The code already flags `bin_sum_vec` as uncapped (L28-ish comment). Recommend:
document the `--mem ≈ f(max_block, N, -T)` relationship, or cap/stream `bin_sum_vec` for huge blocks.

---

## LOW / informational

### L1 Zoom levels capped at uint32_max bin size (makeZoomLevels, bwWrite.c:956)
`zoom` (the zoom-level bin width, stored as uint32) is capped at `((uint32_t)-1)/multiplier`.
For T=72.5e9 the coarsest zoom bin is ≤ ~4.29 Gb → ~17 bins for the whole axis. **Not a
correctness bug** (zoom-record positions are u64 @L970 comment; `nBases` u32 fits since a
zoom bin spans ≤ its u32-capped level). It just bounds whole-axis coarse-read depth — exactly
the gap the separately-planned "zoom-coarse-read lever" addresses. No action for the write path.

### L2 Row-0-gap `continue` is defensive-only (score_one_block, L362)
`if (rb == '-') continue;` skips a row-0 gap column from binning while `uni_col`/`block_start_col`
still count it (uni_col += column_number, L961). For a valid universal MAF row-0 is gap-free so
this never fires. If it ever did, `bin_cnt` would undercount those columns (narrower emitted
interval) but the per-species **aggregate** stays correct (the lift sums over the full query bin
and divides by full span). Defensive; no impact for valid input.

---

## Tried hard to break, could NOT (negative findings — confidence builders)

- **Cross-shard component alignment (the #1 sharding risk): SOUND.** `leaf_id` is the gerp-tree
  post-order leaf index, assigned in `gt_walk` (gerp.c:65-69) and resolved per row by
  `stHash_search(gt->leaf_by_label, genome)` — a pure function of the tree, NOT of block row
  order. Every shard parses the same `# hal` header → identical tree → identical leaf_id and
  identical `gerp_tree_leaf_name(gt, c)`. So all shards emit the same component for a species,
  `merge-bigwig` combines matching components, and `.names` (written from the same order, L872-875)
  labels correctly. A row-order-derived index here would have silently scrambled all 577 tracks;
  it isn't.
- **64-bit clean end to end (T>2³²).** bin index = `(block_start_col+col)/bin_size` (int64),
  `start = bin*bin_size` (int64, max ~7.25e10 fits), `unibw_add` casts to uint64 (L231-232).
  Merge read path is u64 throughout: `bwGetOverlappingIntervalsVec(...,uint64_t start,end)`
  (bigWig.h:394), R-tree `baseStart/baseEnd` u64 with u64 comparisons (bwValues.c:147-241),
  Vec record decode `start=*(uint64_t*)p; end=*(uint64_t*)(p+8)` (bwValues.c Core L53-54),
  `pushIntervalsVec(uint64_t,uint64_t,...)`. No 32-bit truncation of any column position, query
  bound, or byte offset in the write/merge path. (The apes T<2³² could not have exercised this.)
- **Parallelism (-T).** Phase B writes `results[i]` (per batch slot) using `ts[omp_get_thread_num()]`
  (per thread); OpenMP `for` gives distinct `i` per iteration, and Phase C runs after the
  parallel-for barrier. No two threads touch the same `results[i]` or the same `ts[t]`; the
  running binner (`cur_*`) is serial Phase-C only. No race, no lost update, no double-count.
- **Shard / columnRange edges.** `LO` must be bin-aligned (L776) and `HI` bin-aligned unless
  `HI==T` (L822); since global bins are `col/bin_size`-aligned, no bin straddles a shard seam —
  each bin's columns live wholly in one shard. The iterator clips a boundary block into disjoint
  sub-blocks (verified split-vs-single on apes, per L944-951). `merge-bigwig` rejects overlap/
  out-of-order (`iv->start[0] < prev_end`, L133), N/T mismatch (L125), handles empty shards
  (`iv->l==0`), and uses `bwAddIntervalsVec` only for the first non-empty batch then
  `bwAppendIntervalsVec` (L141-146) with `(size_t)off*N` chunk math (no 32-bit overflow).
- **Sentinel + final flush.** First flush is `cur_bin=-1, cur_cnt=0` → both flushers early-return
  on `cnt<=0` (L257/L279), so no bogus negative-start interval. Final post-loop flush (L1070)
  emits the last open bin; empty input → no-op.
- **float precision.** Emitted counts ≤ bin_size (=1000) ≪ 2²⁴, exact as float (`(float)sumvec[c]`,
  L292). bin_cnt ≤ bin_size and columns are disjoint across blocks, so no over-count inflates it.

---

## Recommendations (none blocking)
1. M1: gate the `--depth`-beside-`--perSpecies` value or fix the comment (cosmetic at the default).
2. M2: document `--mem` vs (max block, N, -T); consider bounding `bin_sum_vec` if very long blocks appear.
3. Optional: add ONE at-scale regression that the apes test cannot give — a tiny synthetic with a
   universal column > 2³² (e.g. craft a 2-shard build straddling 2³² + merge, assert merged==single).
   This is the one property with zero current coverage.
