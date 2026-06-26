# Focused review: `--bigwig` 64-bit writer + gating

Scope: the new 64-bit bigWig writer and its rejection gating, as committed in
taffy `7c198ae` and libBigWig `c8d7481`. Reviewed: `UniBW` batch writer in
`taf_depth.c`, the `BIGWIG64_MAGIC` read/gate path in `libBigWig/bwRead.c`, and
the zoom-pyramid construction in `libBigWig/bwWrite.c`.

Three distinct confirmed issues after dedup (the two reported instances of the
error-propagation bug are the same defect). No uncertain items remain.

## 1. Confirmed bugs (prioritized)

### B1 (medium) — mid-stream full-batch write failure is swallowed; can exit 0 with a truncated/corrupt `.bw`
- **Where:** `taf_depth.c:246` (`gerp_flush_bin` discards `unibw_add`'s `int`
  return); call sites `taf_depth.c:921` and `:962`; supporting
  `unibw_flush` `:197-206`, `unibw_add` `:208-215`; the only fatal-setting path
  is `unibw_close` at `:1003`.
- **Mechanism:** `unibw_add` returns `1` only when a *full* 64k batch
  (`n == UNIBW_BATCH`) fails to flush. `gerp_flush_bin` is `void` and ignores
  that return. `unibw_flush` resets `u->n = 0` (dropping the failed 65536-entry
  batch) and merely logs to stderr. `unibw_close` only flushes/propagates the
  *final partial* batch (`n < 65536`). So any earlier full-batch write failure
  (e.g. disk-full during `bwAppendIntervals`) drops 65536 intervals, keeps
  `wrote=true`, and continues; if the trailing flush then succeeds the process
  returns 0 with a silently truncated/corrupt bigWig.
- **Severity rationale:** I/O-error path only — never affects happy-path runs;
  round-trip and bedGraph parity are unaffected. Downgraded from "always exits
  0": the dominant real failure (persistent disk-full) also fails the close-time
  tail flush, which *does* set `fatal=1` and gives a nonzero exit. The truly
  silent exit-0 window is a *transient* mid-stream write error that recovers
  before close — real but narrow. Still a genuine error-propagation defect, and
  relevant here because the per-species/577 bigWigs are multi-GB and disk-full
  is a flagged recurring hazard.
- **Fix:** make `gerp_flush_bin` return `int`, propagate `unibw_add`'s return at
  `:246` up through both call sites (`:921`, `:962`) into `fatal`; on failure
  stop and ideally mark/unlink the partial output so no consumer mistakes it for
  complete.

### B2 (low) — `bwIsBigWig()` magic check is inverted vs. the fork's gate (latent footgun, no live caller)
- **Where:** `libBigWig/bwRead.c:315` (mirror:
  `taffy/submodules/libBigWig/bwRead.c:315`); gate at `bwRead.c:131-135`;
  exported via `bigWig.h:284`.
- **Mechanism:** `bwIsBigWig()` still tests 32-bit `BIGWIG_MAGIC` (returns 1 for
  it, 0 for `BIGWIG64_MAGIC`), while `bwHdrRead` now *rejects* `BIGWIG_MAGIC` and
  only accepts `BIGWIG64_MAGIC`. So `bwIsBigWig` recognizes exactly the files
  `bwOpen("r")` will refuse, and returns 0 for the valid 64-bit files this fork
  produces and reads.
- **Severity rationale:** grep shows **zero** callers of `bwIsBigWig`/`bbIsBigBed`
  in taffy outside the submodule — no runtime path is affected. Latent footgun
  in an exported API only.
- **Fix:** change `bwRead.c:315` to compare against `BIGWIG64_MAGIC`.

### B3 (low) — zoom pyramid: the single global-final zoom bin per level reads sum/mean = 0
- **Where:** `libBigWig/bwWrite.c:968` (conditional `ZSUM` write) and `:1077-1078`
  (`free` of sum/sumsq with no final flush); within
  `updateInterval`/`constructZoomLevels`.
- **Mechanism:** `ZSUM` is written only retroactively on a bin transition
  (`updateInterval` else-branch at `:968`). The add-new-entry path (`:982-993`)
  sets `ZN/ZMIN/ZMAX` but never `ZSUM/ZSQ`. The final bin has no successor
  interval to trigger the retroactive write, and `constructZoomLevels` frees the
  buffers (`:1077-1078`) without a terminal flush. Buffers are `calloc`'d
  (`:1014/1016`), so the unwritten `ZSUM` reads back as 0 → wrong sum/mean for
  that one bin at zoom resolution. `ZN/ZMIN/ZMAX` and positions stay correct.
- **Scope:** narrower than "last bin of each level" — `sum[k]` is flushed at each
  chromosome boundary (tid mismatch in `overlapsInterval` `:918` → else-branch
  `:968`), so only the *single global-final* bin per zoom level is affected.
  Empirically: data-level deficit -127560 over the trailing ~16kb, -0.035% over
  the full range (-32% in the isolated bin).
- **Not introduced by this PR:** structurally identical to upstream libBigWig;
  the 64-bit fork only changed record layout/widths (`:931-946`). Pre-existing
  upstream, hence low for this change (medium is defensible on absolute-
  correctness grounds).
- **Fix (optional, upstreamable):** flush the trailing `sum/sumsq` into the
  buffer's final entry before the `free` at `:1077-1078`.

## 2. Uncertain items

None. All four reported instances resolved to confirmed; instances #1 and #2
are the same defect (B1).

## 3. Verdict

The 64-bit writer and rejection gating are **functionally sound on the happy
path** — round-trip, bedGraph parity, and 64-bit magic gating all hold, and the
two structural defects (B2, B3) are low-impact latent issues, one with no live
caller and one pre-existing upstream and confined to a single zoom bin.

The one issue worth landing before phase 2 builds the vector writer on top is
**B1**: the vector writer will inherit the same batch-flush / error-propagation
pattern, so fixing the dropped-return-value plumbing now (return `int` from
`gerp_flush_bin`, propagate into `fatal`, mark/unlink partial output) prevents
the same silent-corruption-on-write-failure mode from being duplicated into the
multi-GB vector builds where disk-full is a real hazard. B2 and B3 are safe to
defer (B2 a one-line magic fix; B3 an upstream zoom flush), and zoom is usable
as-is.
