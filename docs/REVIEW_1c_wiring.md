# Adversarial review: phase-1c 64-bit wiring + chunking removal

Scope: the depth write path was already proven bit-identical by parity. This
review targets the two open questions of phase-1c: (a) the lift **read** path
over the 64-bit libBigWig fork, and (b) whether dropping the 2e9 chunking is
complete and safe. Two distinct bugs were confirmed; both are correctness
defects on *standard / pre-1c* bigWig inputs. Neither touches the parity-proven
write path.

## 1. Confirmed bugs (prioritized)

### BLOCKING — fork misparses standard UCSC 32-bit bigWigs (R-tree + interval coords)
- **Location:** `~/dev/libBigWig/bwValues.c` `readRTreeIdx` (L43, L45) and
  `bwGetRTreeNode` (L84-103), consumed by `~/dev/taffy/taf_lift.c`.
- **Defect:** the fork reads `baseStart`/`baseEnd` as `uint64_t` with **no magic
  or version gating** in the index parse. UCSC `wigToBigWig` writes them as
  32-bit. The standard R-tree header is 48 bytes; the fork reads 56
  (magic4 + blockSize4 + nItems8 + chrIdxStart4 + **baseStart8** + chrIdxEnd4 +
  **baseEnd8** + idxSize8 + nItemsPerSlot4 + pad4). The extra 8 bytes shift
  `rootOffset` and every subsequent node entry is over-read.
- **Reachability:** the shipped depth workflow shells out to UCSC `wigToBigWig`;
  two sample bigwigs carry the standard 32-bit `BIGWIG_MAGIC`. The fork's
  header check (`bwRead.c:131`) accepts `BIGWIG_MAGIC` too, so the bad file
  *opens*. Running the compiled fork reader on a valid UCSC bigwig segfaults
  while UCSC `bigWigInfo` reads the same file correctly — the file is valid; the
  fork misparses it.
- **Impact:** real-world standard input to `taffy lift --bigwig` yields a
  crashing or all-zero track.
- **Fix:** gate the 64-bit interpretation on `magic == BIGWIG64_MAGIC`. When the
  file carries the standard magic, parse `baseStart`/`baseEnd` (and the node
  entries in `bwGetRTreeNode`) as 32-bit, or reject the file with a clear error.
  Do not read 64-bit fields off a 32-bit file unconditionally.

### MEDIUM — pre-1c multi-chrom universal bigWig silently returns a partial-zero track
- **Location:** `~/dev/taffy/taf_lift.c:1681-1704` (axis guard), `:1720`
  (`uchrom = chrom[0]`), `:1731` (`base = 0`); enabled by
  `~/dev/libBigWig/bwRead.c:131` (accepts both magics).
- **Defect:** the old universal-depth bigWig was multi-chrom — `uni0,uni1,...,uniK`
  with `chunk = col/2e9`, `localpos = col%2e9`, `uni0` capped at `[0,2e9)`. The
  new single-axis reader's guard (L1682) inspects only `chrom[0]` (`"uni0"`), so
  `integer_axis = true` and the old file passes. The reader then sets `base = 0`
  (L1731), `uchrom = chrom[0] = "uni0"` (L1720), and queries
  `bwGetOverlappingIntervals(bw, "uni0", lo, hi)` with **absolute** universal
  columns from `tui_query`. For columns `>= 2e9` (which live on `uni1..uniK`),
  `uni0` has no data → `bwGetOverlappingIntervals` returns NULL → the
  `if (ov != NULL)` guard yields `sum=0, cnt=0 → d=0`.
- **Impact:** columns `< 2e9` read correctly; the **majority** of a large genome
  (T=72.5e9 for the 577-way) silently returns 0. No crash. This is a
  partial-zero track (not the uniformly-zero one the original claim's wording
  emphasized).
- **Severity rationale:** fires only on an obsolete pre-1c input; new
  single-chrom `uni0` files read correctly and the writer changed in lockstep.
- **Fix:** add a magic-or-layout check (e.g. `magic == BIGWIG64_MAGIC` and/or
  `nKeys == 1` with a single `uni0` axis) and **reject** a multi-chrom
  `uni0,uni1,...` file, matching the project's "old .tui REJECTED" policy. This
  converts a silent miscompute into a clean error.

## 2. Uncertain items

None. Both confirmed claims were re-verified against current source
(`bwValues.c:21-59`, `bwRead.c:124-144`, `taf_lift.c:1675-1744`).

## 3. Verdict

**Not safe to commit as-is.** The write path is parity-proven, and the
chunking-removal logic on a *correctly-formatted single-axis 64-bit* input is
sound (`base=0`, absolute columns on `uni0`, no straddle pre-split needed). But
the **lift read path is NOT robust**: it silently mishandles two input classes
that are reachable today because `bwRead.c:131` opens both magics —

1. (blocking) standard UCSC 32-bit bigWigs, which the shipped `wigToBigWig`
   workflow actually produces, crash/zero out; and
2. (medium) pre-1c multi-chrom universal bigWigs return a majority-zero track.

Both share one root cause: the reader accepts files it cannot correctly
interpret instead of gating on `BIGWIG64_MAGIC` / single-`uni0` layout. **Add
magic/layout gating (reject-or-correctly-parse) before commit.** After that gate
is in place, the wiring and chunking removal are correct.
