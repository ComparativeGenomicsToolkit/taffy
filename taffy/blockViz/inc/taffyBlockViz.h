/*
 * taffyBlockViz: C API mirroring HAL's halBlockViz for use by genome
 * browsers (or other consumers) that need block-level alignment data
 * out of a universal .tui index.  All structs / functions follow the
 * hal_X / halX shapes so a thin compat shim (typedefs + #defines)
 * gets existing hal blockViz callers running against taffy.
 *
 * SCOPE -- THIS IS THE INITIAL BASIC-INTEGRATION SUBSET:
 *   Supported:
 *     - taffyOpen / taffyClose / taffyCloseGenome
 *     - taffyGetSpecies                       (list genomes)
 *     - taffyGetChroms                        (list chroms in a genome)
 *     - taffyGetBlocksInTargetRange[_filterByChrom]
 *     - taffyGetDna                           (STUB: returns "NNN..."
 *                                              of requested length)
 *     - taffyGetGenomeMetadata                (returns NULL today)
 *     - taffyFree* for every returned list
 *
 *   NOT supported (will set errStr and return NULL / -1):
 *     - tReversed != 0
 *     - mapBackAdjacencies != 0
 *     - coalescenceLimitName != NULL          (taffy has no analog)
 *     - dupMode != TAFFY_QUERY_AND_TARGET_DUPS (only the default
 *                                              "include all paralogs"
 *                                              shape is wired today)
 *     - LOD-related calls                     (no LOD concept in .tui)
 *     - MAF emit                              (not in initial cut)
 *
 * `errStr`:  if non-NULL, *errStr is set to a malloc'd error message
 * (caller frees) on failure.  If NULL, the call still fails the same
 * way but the message is dropped on the floor (matching halBlockViz).
 *
 * COORDINATES: BED-like, 0-based, half-open, forward-strand-relative
 * for both qStart and tStart.  `strand` in taffy_block_t encodes the
 * alignment strand of qChrom relative to the target reference.
 */

#ifndef _TAFFY_BLOCK_VIZ_H
#define _TAFFY_BLOCK_VIZ_H

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* keep integer type definition in one place; matches HAL exactly. */
typedef long taffy_int_t;

/** range of coordinates in the target */
struct taffy_target_range_t {
    struct taffy_target_range_t *next;
    taffy_int_t tStart;
    taffy_int_t size;
};

/** paralogous ranges in the *target* genome (e.g. self-alignment of
 * the target that landed in the same column).  Browsers render these
 * as a separate track ("blue lines") so they're split out from the
 * main mappedBlocks list. */
struct taffy_target_dupe_list_t {
    struct taffy_target_dupe_list_t *next;
    taffy_int_t id;
    struct taffy_target_range_t *tRange;
    char *qChrom;
};

/** Per-chain summary returned alongside the block lists.  One entry per
 * chainId that appears in mappedBlocks, in score-descending order.  Use
 * the id field to look up which mappedBlocks (via taffy_block_t.chainId)
 * belong to this chain.  totalBp is the chain's covered bp; totalScore
 * is its score; nAlns is its aln count.  Coverage (wide) queries report
 * one entry per query chromosome (matching the per-qChrom chainIds). */
struct taffy_chain_summary_t {
    struct taffy_chain_summary_t *next;
    taffy_int_t id;
    taffy_int_t totalScore;
    taffy_int_t totalBp;
    taffy_int_t nAlns;
};

/** Mapped + target-dupe block results from a single
 * taffyGetBlocksInTargetRange() call. */
struct taffy_block_results_t {
    struct taffy_block_t *mappedBlocks;
    struct taffy_target_dupe_list_t *targetDupeBlocks;
    struct taffy_chain_summary_t *chainSummaries;   /* score-desc; one per emitted chainId */
};

/** One contiguous aligned block (one gap-free run of qChrom that
 * lifts to a contiguous run on the target reference).
 * Coordinates are forward-strand-relative.
 */
struct taffy_block_t {
    struct taffy_block_t *next;
    char *qChrom;
    taffy_int_t tStart;
    taffy_int_t qStart;
    taffy_int_t size;
    char strand;                /* '+' or '-' */
    /* Opaque, always-positive (1-based) group id for snake-trace
     * rendering: blocks of the same chainId form one snake.  Stable
     * WITHIN one taffyGetBlocksInTargetRange result; NOT stable across
     * queries.  Use it instead of inferring adjacency from
     * qChrom+strand+tStart (which collides on paralog chains sharing
     * qChrom+strand).
     *   - Detail (narrow) queries: one id per real alignment chain (an
     *     orthologous arm is one id; off-diagonal paralogs get their own).
     *   - Coverage (wide) queries that hit the positional-coarsen path:
     *     one id PER DIAGONAL -- the per-bin coverage blocks are grouped by
     *     the same collinear rule as detail, so an orthologous arm is one
     *     row and each paralog copy (a different diagonal that dominates
     *     some bins) gets its own, instead of every copy of a qChrom
     *     winding through a single per-chromosome row.
     * chainSummaries carries the matching id->score/bp entries either way.
     */
    taffy_int_t chainId;
    char *qSequence;            /* query DNA, NULL unless seqMode != NO_SEQUENCES */
    char *tSequence;            /* target DNA, NULL unless seqMode != NO_SEQUENCES */
};

/** Genome (species) metadata. */
struct taffy_species_t {
    struct taffy_species_t *next;
    char *name;
    taffy_int_t length;         /* sum of chrom lengths */
    taffy_int_t numChroms;
    char *parentName;           /* "" (no tree-walk needed for browser blocks) */
    double parentBranchLength;  /* 0.0 */
};

/** Chromosome / sequence within a genome. */
struct taffy_chromosome_t {
    struct taffy_chromosome_t *next;
    char *name;
    taffy_int_t length;
};

/** Optional key/value metadata; today taffyGetGenomeMetadata returns NULL. */
struct taffy_metadata_t {
    struct taffy_metadata_t *next;
    char *key;
    char *value;
};

/** seqMode -- whether to fetch sequence into qSequence / tSequence on
 * each block.  Same values as HAL's hal_seqmode_type_t.
 * Initial cut: only HAL_NO_SEQUENCES is honoured; the others set errStr. */
typedef enum {
    TAFFY_NO_SEQUENCES      = 0,
    TAFFY_HARD_MASKED_SEQ   = 1,
    TAFFY_SOFT_MASKED_SEQ   = 2,
    TAFFY_LOD0_SEQ          = 3
} taffy_seqmode_type_t;

/** dupMode -- which paralog edges to include.  Initial cut: only
 * TAFFY_QUERY_AND_TARGET_DUPS is honoured. */
typedef enum {
    TAFFY_NO_DUPS               = 0,
    TAFFY_QUERY_DUPS            = 1,
    TAFFY_QUERY_AND_TARGET_DUPS = 2
} taffy_dup_type_t;

/* ------------------------------------------------------------------ */
/* Lifecycle                                                           */
/* ------------------------------------------------------------------ */

/** Open a universal .tui file; returns an opaque non-negative handle
 * on success or -1 on failure (with errStr set if non-NULL).
 *
 * @param tuiPath path to the .tui sidecar OR to its companion .maf.gz
 *  (we accept either and resolve the .tui via the same convention as
 *  `taffy lift -i`).
 */
int taffyOpen(const char *tuiPath, char **errStr);

/** Close a handle returned by taffyOpen().  Returns 0 on success, -1
 * on failure. */
int taffyClose(int taffyHandle, char **errStr);

/** Release any cached per-genome lift table held by `genomeName` on
 * the open handle.  Useful to drop memory after a browser pan moves
 * to a new species.  Returns 0 on success / -1 on failure. */
int taffyCloseGenome(int taffyHandle, const char *genomeName, char **errStr);

/* ------------------------------------------------------------------ */
/* Chain tuning (per handle)                                           */
/* ------------------------------------------------------------------ */

/** Configure the chaining cost / gap-break parameters used by
 * taffyGetBlocksInTargetRange when partitioning visited runs into
 * chains (the per-handle defaults are 0 / 1 / 10 Mb, matching
 * TAFFY_CHAIN_DEFAULT_{OPEN,EXTEND,MAX_GAP}).
 *
 * Semantics of each param:
 *   chain_open    -- fixed cost per chain join.  0 = always pays to
 *                    join collinear runs.  Higher values prevent
 *                    short alignments from chaining across gaps.
 *   chain_extend  -- per-bp cost added to chain_open scaled by
 *                    (q_gap + t_gap).  1 = a 100-bp gap costs 100
 *                    (on top of chain_open).
 *   max_gap_length -- any candidate join whose q_gap or t_gap exceeds
 *                    this is rejected outright -- the chain breaks.
 *                    Pass INT64_MAX to disable the gap-break.
 *
 * Empirical sensitivity on apes universal MAF (sweep across 4 hg38 chr5
 * regions, max_gap from 1 KB to 100 MB): primary-chain selection is
 * INSENSITIVE to max_gap_length over 5 orders of magnitude -- the
 * gap-cost (chain_extend * gap_bp) already breaks chains at appropriate
 * points before any hard cutoff fires.  chain_extend IS the
 * sensitive knob -- ranging it from 0 to 1e4 changes the secondary-
 * chain partition count by ~4x while leaving primary unchanged.  At
 * snake-track scales, default 10 Mb max_gap is functionally equivalent
 * to INT64_MAX; the hard cutoff is a safety net for pathological
 * (e.g. cross-chromosome) cases.
 *
 * Pass -1 in any field to leave that field unchanged.  All other
 * negative values are rejected (returns -1 with errStr).  Returns 0
 * on success.  Thread-safe (takes the same lock as the query path).
 */
int taffySetChainParams(int taffyHandle,
                        int64_t chain_open,
                        int64_t chain_extend,
                        int64_t max_gap_length,
                        char **errStr);

/** Read back the currently-configured chain parameters for the
 * handle.  Any non-NULL output pointer receives the current value;
 * NULL pointers are ignored (call with NULLs for ones you don't
 * care about).  Returns 0 on success, -1 on invalid handle. */
int taffyGetChainParams(int taffyHandle,
                        int64_t *chain_open,
                        int64_t *chain_extend,
                        int64_t *max_gap_length,
                        char **errStr);

/** Coverage budget: the number of bins the wide-zoom (coarsened) path
 * tiles the query into, so coverage output is bounded at ~this many blocks
 * (one per occupied bin).  Default 1000 (= HAL's halSnakeTrack.c
 * NUM_LEVELS).  The detail (narrow) path is deliberately NOT bounded by
 * this: it emits EVERY block in the span, so a fixed locus renders
 * identically regardless of where the window edges fall (panning-stable) --
 * the span itself bounds the count.  The detail<->coverage switch is
 * span-based: coverage kicks in when (span / this) exceeds a coverage-bin
 * width (~10 kb), i.e. at chromosome-scale zoom, NOT when a window happens
 * to contain many blocks.  Lower it for a coarser overview, raise it for a
 * finer coverage tiling.  Must be >= 1 (rejected with errStr set otherwise).
 *
 * Returns 0 on success, -1 on invalid handle or invalid n. */
int taffySetMaxOutputBlocks(int taffyHandle, int64_t n, char **errStr);

/** Read back the currently-configured cap.  *n receives the value;
 * NULL is accepted as a no-op.  Returns 0 on success, -1 on invalid
 * handle. */
int taffyGetMaxOutputBlocks(int taffyHandle, int64_t *n, char **errStr);

/** Noise filter: drop output blocks too small to render meaningfully.
 * ON by default with conservative fractions (spanFrac=0.001 ~ 1/1000 px,
 * relFrac=0.1), so it declutters out of the box; pass BOTH 0 to disable.  A
 * block is dropped iff its size is BOTH < spanFrac of the query window AND <
 * relFrac of the largest output block -- equivalently, < min(window*spanFrac,
 * maxBlock*relFrac).
 *
 * Rationale: the query window span stands in for the track's pixel width,
 * so passing spanFrac ~ 1/pixels makes this a sub-pixel test with no pixel
 * count in the API.  The relative term is self-protecting: a uniformly-small
 * region survives (nothing is small relative to its own max) while slivers
 * beside a real feature drop.  Dropped blocks are simply not emitted, so the
 * snake renderer bridges the surviving neighbours -- removing a sliver also
 * removes its hairline connectors.  The largest block always survives, so the
 * result is never emptied.
 *
 * Applies to both detail and coverage output (at wide zoom the coverage bins
 * are themselves ~1/max_output_blocks of the window, so the sparse ones are
 * the clutter this removes).  Both fractions must be > 0 to engage;
 * spanFrac=1 makes it relative-only, relFrac=1 makes it window-only.
 *
 * Returns 0 on success, -1 on invalid handle or out-of-[0,1] fraction
 * (set errStr).  Thread-safe. */
int taffySetMinBlockFilter(int taffyHandle, double spanFrac, double relFrac,
                           char **errStr);

/** Read back the min-block filter fractions.  Either out-pointer may be
 * NULL.  Returns 0 on success, -1 on invalid handle. */
int taffyGetMinBlockFilter(int taffyHandle, double *spanFrac, double *relFrac,
                           char **errStr);

/** Per-handle run floor: drop input alignment runs shorter than this (bp)
 * BEFORE chaining.  On wide overview queries this collapses the chain DP
 * over the sub-pixel micro-alignment tail (a whole-chromosome view is
 * ~90% sub-100bp runs / <1% of bp), making segdup-dense chromosomes
 * ~6-10x faster.  It slightly relaxes paralogy dedup (a few extra paralog
 * blocks survive) -- sub-pixel, acceptable for overviews.
 *
 * minRun = -1 (default) = AUTO per query span (span/500000: ~100-500 bp at
 * whole-chrom, 0 below ~500 kb so detail queries stay byte-exact); 0 = off;
 * > 0 = explicit floor.  Returns 0 on success, -1 on invalid handle or
 * minRun < -1 (set errStr).  Thread-safe. */
int taffySetMinRunSize(int taffyHandle, int64_t minRun, char **errStr);

/** Read back the run floor (-1 = auto).  Out-pointer may be NULL.
 * Returns 0 on success, -1 on invalid handle. */
int taffyGetMinRunSize(int taffyHandle, int64_t *minRun, char **errStr);

/** Per-handle overlap-fraction threshold for the chain paralogy filter
 * applied after taffy_chain has assigned chain ids.  Chains are walked
 * in score-desc order and accepted iff their union-of-aln q-coverage
 * overlaps the already-kept chains' q-coverage by at most this
 * fraction of the candidate chain's own q-bp.
 *
 * Semantics:
 *   frac = 0.5   (DEFAULT)   balanced: drop a chain only when >50% of
 *                            its q-coverage is redundant with a kept
 *                            chain.  Filters true paralogs while keeping
 *                            real ortholog coverage; inversions and
 *                            disjoint-q chains stay regardless of strand.
 *   frac = 0.0               strict: drop on ANY q-overlap.  WRONG for
 *                            the universal .tui -- the alignment is split
 *                            into many chains with tiny boundary
 *                            overlaps, so this culls ~17% of real
 *                            ortholog coverage (e.g. hs1->hg38).  Use
 *                            only on a pre-merged / unfragmented source.
 *   frac in (0,1]            accept chains overlapping up to this
 *                            fraction of their own q-bp.
 *   frac = -1.0              disable the filter entirely (browser sees
 *                            every chain, modulo the max_output_blocks
 *                            budget).
 *
 * Returns 0 on success, -1 on invalid handle or out-of-range frac
 * (set errStr).  Thread-safe. */
int taffySetChainOverlapFrac(int taffyHandle, double frac, char **errStr);

/** Read back the currently-configured overlap-frac threshold.  *frac
 * receives the value; NULL is accepted as a no-op.  Returns 0 on
 * success, -1 on invalid handle. */
int taffyGetChainOverlapFrac(int taffyHandle, double *frac, char **errStr);

/* ------------------------------------------------------------------ */
/* Free functions for the returned linked lists.                       */
/* ------------------------------------------------------------------ */

void taffyFreeBlockResults(struct taffy_block_results_t *results);
void taffyFreeBlocks(struct taffy_block_t *block);
void taffyFreeTargetDupeLists(struct taffy_target_dupe_list_t *dupes);
void taffyFreeSpeciesList(struct taffy_species_t *species);
void taffyFreeChromList(struct taffy_chromosome_t *chromosome);
void taffyFreeMetadataList(struct taffy_metadata_t *metadata);

/* ------------------------------------------------------------------ */
/* Block query (the hot path).                                         */
/* ------------------------------------------------------------------ */

/** Return blocks of `qSpecies` aligned within `tChrom`:[tStart, tEnd)
 * on `tSpecies`.  Output is a linked list of taffy_block_t.
 *
 * Two regimes, chosen by ZOOM (query span), NOT by block count -- so a
 * fixed locus renders identically as you pan a fixed zoom (no blinking):
 *   - DETAIL (narrow span): every aligned block is returned at full
 *     resolution; bounded by local alignment density, not max_output_blocks.
 *   - COVERAGE (wide span, once span/max_output_blocks exceeds a ~10 kb
 *     bin): the span is tiled into <= max_output_blocks absolute-aligned
 *     bins, one coverage block per occupied bin; per-dupe detail is
 *     dropped.
 *
 * Pipeline:
 *   1. The tSpecies range is mapped to universal-column intervals; for
 *      each, the qSpecies runs are visited and chained + collinear-merged
 *      in windows (cheap, bounded memory) into merged blocks.
 *   2. Blocks are grouped into chains by collinearity (one diagonal per
 *      chain); the overlap-frac paralogy filter (taffySetChainOverlapFrac,
 *      default 0.5) drops chains mostly redundant with a higher-scoring
 *      kept chain.
 *   3. Within each chain, redundant nested fragments are merged/dropped and
 *      partial overlaps trimmed so blocks tile the reference WITHOUT
 *      overlapping -- a snake renderer then has no fake loop-backs.
 *   4. Detail emits every surviving block; coverage bins them (grouping
 *      bins per diagonal).  Each taffy_block_t carries a real
 *      (tStart, qStart, size, strand) and a chainId; routing follows dupMode.
 *
 * Implication for snake-track callers: at paralog-rich detail queries
 * (segdup hotspots) you'll see all the paralog dupes; pass
 * dupMode = TAFFY_NO_DUPS for primary-only.
 *
 * @param taffyHandle      handle from taffyOpen
 * @param qSpecies         genome to fetch blocks FROM
 * @param tSpecies         genome whose coordinates [tStart, tEnd) frame
 *                         the query range (any genome the .tui knows --
 *                         leaf or ancestor)
 * @param tChrom           chromosome name in tSpecies (bare, no genome prefix)
 * @param tStart           half-open range start
 * @param tEnd             half-open range end (0 = chrom length)
 * @param tReversed        UNSUPPORTED -- must be 0
 * @param seqMode          UNSUPPORTED for !=TAFFY_NO_SEQUENCES (stub today)
 * @param dupMode          UNSUPPORTED for !=TAFFY_QUERY_AND_TARGET_DUPS
 * @param mapBackAdjacencies UNSUPPORTED -- must be 0
 * @param coalescenceLimitName UNSUPPORTED -- must be NULL
 * @param errStr           error message destination (caller-frees) or NULL
 *
 * @return malloc'd taffy_block_results_t (free with taffyFreeBlockResults),
 *         or NULL on error.
 */
struct taffy_block_results_t *taffyGetBlocksInTargetRange(
    int taffyHandle,
    const char *qSpecies, const char *tSpecies,
    const char *tChrom, taffy_int_t tStart, taffy_int_t tEnd,
    taffy_int_t tReversed,
    taffy_seqmode_type_t seqMode,
    taffy_dup_type_t dupMode,
    int mapBackAdjacencies,
    const char *coalescenceLimitName,
    char **errStr);

/** Same as taffyGetBlocksInTargetRange but only emits blocks whose
 * qChrom matches `qChrom` (still on `qSpecies`).  qChrom == NULL is
 * equivalent to the unfiltered call. */
struct taffy_block_results_t *taffyGetBlocksInTargetRange_filterByChrom(
    int taffyHandle,
    const char *qSpecies, const char *tSpecies,
    const char *tChrom, taffy_int_t tStart, taffy_int_t tEnd,
    taffy_int_t tReversed,
    taffy_seqmode_type_t seqMode,
    taffy_dup_type_t dupMode,
    int mapBackAdjacencies,
    const char *qChrom,
    const char *coalescenceLimitName,
    char **errStr);

/* ------------------------------------------------------------------ */
/* Introspection                                                       */
/* ------------------------------------------------------------------ */

/** List of all genomes in the .tui (sorted alphabetically).
 * `parentName` and `parentBranchLength` are stubbed to "" / 0.0 in
 * this initial cut.  Genome names -- including any dotted version
 * suffix (NCBI accessions like "GCA_028858775.2") -- come straight
 * from the .tui's "O g" genome-record roster, which every .tui
 * carries. */
struct taffy_species_t *taffyGetSpecies(int taffyHandle, char **errStr);

/** List of chromosomes in the named genome (sorted alphabetically). */
struct taffy_chromosome_t *taffyGetChroms(int taffyHandle, const char *speciesName, char **errStr);

/** Fetch DNA in [start, end) of `chromName` on `speciesName`.
 *
 * STUB IN THIS INITIAL CUT: always returns a malloc'd string of
 * (end - start) 'N' characters and sets no error.  Browser code that
 * pulls sequence to render won't crash; output is just opaque.  Will
 * be replaced with a real fetch once we wire up DNA access. */
char *taffyGetDna(int taffyHandle, const char *speciesName, const char *chromName,
                  taffy_int_t start, taffy_int_t end, char **errStr);

/** Per-genome metadata as a linked list.  Returns NULL today (.tui
 * does not yet carry per-genome metadata). */
struct taffy_metadata_t *taffyGetGenomeMetadata(int taffyHandle, const char *genomeName, char **errStr);

#ifdef __cplusplus
}
#endif

#endif /* _TAFFY_BLOCK_VIZ_H */
