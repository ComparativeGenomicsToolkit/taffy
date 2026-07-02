/*
 * taffy summary: per-reference bigMafSummary (bed3+4) computed from a universal
 * MAF/TAF, for ONE reference genome at a time.
 *
 * Replicates UCSC mafToBigMafSummary's pairwise scoring + region merging (kent
 * src/hg/utils/mafToBigMafSummary/mafToBigMafSummary.c, scorer from
 * src/lib/mafScore.c), but driven from taffy's universal alignment instead of a
 * reference-anchored MAF.  For each reference (master) row in a block we emit,
 * for every other LEAF species present, a summary record scored master-vs-that-
 * species in the reference's coordinates.  Ancestor rows are dropped.
 *
 * SINGLE-COVER: a reference-anchored MAF is single-coverage.  The universal MAF
 * is single-cover *by base* already -- the cactus-hal2maf --universal
 * --noRefDupes invariant (tui.h) guarantees every base of every genome maps to
 * at most one column.  What it does NOT prevent is several rows of the SAME
 * genome in one block (paralogs aligned to the same ancestral column -- both
 * the reference and the species exhibit this).  Two consequences, both handled
 * here:
 *   - a species on N rows would otherwise emit N records at the SAME reference
 *     position; we single-cover the species in-block by keeping only its
 *     best-scoring row (--chainOverlapFrac 0, the default).
 *   - the reference on N rows is N DISJOINT reference positions (disjoint by the
 *     --noRefDupes invariant); we treat EACH as a master (the old first-row-only
 *     pick silently dropped the rest).
 * The result is one record per (reference base, species) = the reference-
 * anchored single-coverage the real bigMafSummary has.  No cross-block chaining
 * is needed (the reference is already single-cover by base).
 *
 * MEMORY / SCALE: universal block order is the ROW-0 ANCESTOR's order, in which a
 * LEAF reference is badly shuffled (measured ~40% chromosome-jumps between
 * consecutive reference rows).  The raw scored records are ~200-300 GB per DEEP
 * reference (hg38-class), ~90 TB across a 577-way -- far too much to spill to disk
 * (per-reference OR per-column both hit this).  So we NEVER spill: each reference
 * accumulates its records IN RAM (RefAcc -- one interleaved array, chromId per
 * record) and COMPACTS incrementally -- once it grows past a --mergeMem-derived
 * budget we sort by (chrom,src,start) and window-bin IN PLACE, a ~50-70x collapse
 * to the merged (bed) size.  Peak RAM ~= the resident references' merged size, not
 * the raw; scratch is ZERO.  Tune with -M / --refSubset (fewer refs per pass).
 *
 * The region-merge had to change: the kent adaptive mergeGap+minSize merge is
 * ORDER-DEPENDENT (a run flushed at minSize can't re-absorb a record the shuffle
 * delivers inside it later -> overlaps + score corruption), so it CANNOT run
 * incrementally on the shuffled scan.  Instead each (chrom, src, floor(start/
 * minSize)) fixed window becomes ONE coverage-weighted run -- order-independent and
 * exactly associative, so incremental compaction == one big merge.  Result is
 * BROWSER-EQUIVALENT at ~minSize resolution, not byte-identical to kent (mergeGap
 * is bridged only WITHIN a window).  Output is (chrom,start)-sorted, bedToBigBed-
 * ready.  Score is a DOUBLE (the coverage-weighted average needs the precision).
 *
 * Output is bed3+4 (chrom start end src score leftStatus rightStatus), the same
 * schema as kent's mafSummary.as.  Feed straight to
 * `bedToBigBed -type=bed3+4 -as=mafSummary.as -tab`.
 *
 *  Released under the MIT license, see LICENSE.txt
 */

#define _GNU_SOURCE   /* mkstemp/fdopen */
#include "taf.h"
#include "tui.h"
#include "gerp.h"
#include "line_iterator.h"
#include "sonLib.h"
#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <omp.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>

/* ----- defaults: match kent mafToBigMafSummary ----- */
static int mergeGap   = 500;
static int minSize    = 10000;
static int maxSize    = 50000;
static int64_t minSeqSize = 1;
/* per-species in-block dedup: 0.0 = single-cover (keep each species' best
 * paralog row, matching the reference-anchored summary); 1.0 = keep all rows
 * (multi-cover; the region merge then length-weighted-averages them). */
static double  chainOverlapFrac = 0.0;
static int     nThreads = 0;   /* 0 => ncores */
/* per-chrom sort+merge holds the whole chrom in RAM (recs + merged); with many
 * threads on a scaffold-dense reference that sums past the box.  Cap concurrent
 * in-RAM chrom bytes to this budget (workers wait their turn; a single chrom
 * larger than the budget still runs, alone). */
static int64_t mergeMemBudget = 8LL << 30;   /* 8 GiB */

/* ===================================================================== */
/* HOXD70 pairwise scorer -- ported verbatim from kent src/lib/mafScore.c */
/* (mafScoreUseTraditional + mafScoreRangeMultiz, pairwise specialization) */
/* ===================================================================== */
#define NACHARS 128
#define DASH '-'
typedef int ss_t[NACHARS][NACHARS];
static ss_t ss;
static int gop[256], gtype[128];
static const unsigned char nchars[] = "ACGT";
static const int HOXD70_sym[4][4] = {
    {  91, -114,  -31, -123 },
    {-114,  100, -125,  -31 },
    { -31, -125,  100, -114 },
    {-123,  -31, -114,   91 },
};
#define SS(c,d) ss[(unsigned char)(c)][(unsigned char)(d)]
#define GTY(c)  gtype[(unsigned char)(c)]
#define GAP(w,x,y,z) gop[(GTY(w)<<6)+(GTY(x)<<4)+(GTY(y)<<2)+GTY(z)]

static void DNA_scores(ss_t s, const int matrix[4][4]) {
    int i, j, bad, a, b, A, B;
    for (i = 0; i < NACHARS; ++i)
        for (j = 0; j < NACHARS; ++j)
            s[i][j] = -100;
    for (i = 0; i < 4; ++i) {
        A = nchars[i]; a = tolower(A);
        for (j = 0; j < 4; ++j) {
            B = nchars[j]; b = tolower(B);
            s[A][B] = s[a][B] = s[A][b] = s[a][b] = matrix[i][j];
        }
    }
    bad = -1000;
    for (i = 0; i < NACHARS; ++i)
        s['X'][i] = s[i]['X'] = s['x'][i] = s[i]['x'] = bad;
}

static void gap_costs(int gap_open) {
    int i, X, D;
    for (i = 0; i < 128; ++i) gtype[i] = 0;
    D = DASH; gtype[D] = 1;
    for (i = 0; i < 256; ++i) gop[i] = 0;
    X = 'A';
    GAP(X,X,X,D) = gap_open;
    GAP(X,X,D,X) = gap_open;
    GAP(X,D,D,X) = gap_open;
    GAP(D,X,X,D) = gap_open;
    GAP(D,D,X,D) = gap_open;
    GAP(D,D,D,X) = gap_open;
}

static void scorer_init(void) {   /* mafScoreUseTraditional */
    int i;
    DNA_scores(ss, HOXD70_sym);
    int E = 30, O = 400;
    for (i = 0; i < 128; ++i) ss[i][DASH] = ss[DASH][i] = -E;
    ss[DASH][DASH] = 0;
    gap_costs(O);
}

/* Pairwise [0,1] score of master vs other over a block's text columns.
 * `master`/`other` are the two rows' alignment strings (length textSize);
 * masterBases is the master row's non-gap base count.  Mirrors kent
 * scorePairwise(): score the columns covering masterBases reference bases,
 * normalize by masterBases, rescale (score+100)/200, clamp [0,1]. */
static double score_pairwise(const char *master, const char *other,
                             int64_t textSize, int64_t masterBases) {
    if (masterBases <= 0) return 0.0;
    int64_t deltaB = masterBases, endT;
    for (endT = 0; endT < textSize; ++endT) {
        if (deltaB <= 0) break;
        if (master[endT] != DASH) deltaB -= 1;
    }
    double score = 0.0;
    for (int64_t i = 0; i < endT; ++i) {
        char br = master[i], bi = other[i];
        score += SS(br, bi);
        if (i > 0) score -= GAP(master[i-1], other[i-1], br, bi);
    }
    score = score / (double) masterBases;
    score = (score - (-100.0)) * (1.0 / 200.0);
    if (score < 0.0) score = 0.0;
    if (score > 1.0) score = 1.0;
    return score;
}

/* ===================================================================== */
/* name interning: string <-> uint32 id (per-id reverse name array)       */
/* ===================================================================== */
typedef struct { stHash *map; char **names; int64_t n, cap; } Interner;

static Interner *interner_new(void) {
    Interner *in = st_calloc(1, sizeof(Interner));
    in->map = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    in->cap = 256; in->n = 0; in->names = st_malloc(in->cap * sizeof(char *));
    return in;
}
/* returns id, adding s on first sight. NOT thread-safe (scan is single-thread). */
static uint32_t intern_id(Interner *in, const char *s) {
    void *v = stHash_search(in->map, (void *) s);
    if (v != NULL) return (uint32_t)(intptr_t) v - 1;
    char *copy = stString_copy(s);
    uint32_t id = (uint32_t) in->n;
    stHash_insert(in->map, copy, (void *)(intptr_t)(id + 1));
    if (in->n == in->cap) { in->cap *= 2; in->names = st_realloc(in->names, in->cap * sizeof(char *)); }
    in->names[in->n++] = copy;
    return id;
}
static const char *intern_name(Interner *in, uint32_t id) { return in->names[id]; }  /* read-only after scan */
static void interner_free(Interner *in) {
    /* keys (the copies) are freed by the map's key destructor; names[] alias them */
    stHash_destruct(in->map); free(in->names); free(in);
}

/* ===================================================================== */
/* packed records + per-reference-chrom buckets                           */
/* ===================================================================== */
typedef struct { uint32_t chromId; int64_t start, end; uint32_t srcId; double score; } PackedRec;  /* start/end int64: VGP sequences exceed 4.29 Gb (uint32 wraps) */

static int cmp_compact(const void *a, const void *b) {   /* (chromId, srcId, start) -- compaction */
    const PackedRec *x = a, *y = b;
    if (x->chromId != y->chromId) return x->chromId < y->chromId ? -1 : 1;
    if (x->srcId != y->srcId) return x->srcId < y->srcId ? -1 : 1;
    if (x->start != y->start) return x->start < y->start ? -1 : 1;
    return 0;
}
static const int64_t *g_lexrank;                          /* chromId -> lexicographic rank (finalize) */
static int cmp_output(const void *a, const void *b) {     /* (chrom lex, start, srcId) -- bed order */
    const PackedRec *x = a, *y = b;
    int64_t rx = g_lexrank[x->chromId], ry = g_lexrank[y->chromId];
    if (rx != ry) return rx < ry ? -1 : 1;
    if (x->start != y->start) return x->start < y->start ? -1 : 1;
    if (x->srcId != y->srcId) return x->srcId < y->srcId ? -1 : 1;
    return 0;
}

/* Bin records into fixed minSize-aligned windows per (chrom, src) -- an
 * order-INDEPENDENT summary that survives the shuffled scan + incremental
 * compaction.  Each (chrom, src, floor(start/minSize)) window becomes ONE run
 * [min_start,max_end] with a coverage-weighted score = Sum(score*len)/(max_end-
 * min_start).  Reads a (chromId,srcId,start)-sorted `r` of n records/runs; a prior
 * run re-bins by its OWN window (its score*span re-supplying its total), so this is
 * exactly associative + idempotent.
 *
 * WHY not the kent adaptive mergeGap+minSize merge: that is order-DEPENDENT.  A leaf
 * reference is badly shuffled in universal (row-0 ancestor) order, so a run flushed
 * at minSize during one compaction cannot absorb a record the shuffle delivers
 * INSIDE it later -> overlapping runs + score corruption (observed: scores > 1, 6x
 * rows).  Fixed-grid binning never flushes-before-complete, so it is safe; result is
 * browser-equivalent at ~minSize resolution.  Merges IN PLACE (each emit writes
 * r[out] with out < i, no transient buffer); returns the run count. */
static int64_t region_merge_inplace(PackedRec *r, int64_t n) {
    int64_t out = 0, win = (minSize > 0 ? minSize : 1);
    int have = 0; uint32_t pch = 0, psrc = 0; int64_t pstart = 0, pend = 0, pwin = 0; double total = 0;
#define EMIT() do { r[out].chromId = pch; r[out].start = pstart; r[out].end = pend; \
        r[out].srcId = psrc; r[out].score = total / (double) (pend - pstart); out++; } while (0)
    for (int64_t i = 0; i < n; i++) {
        uint32_t ch = r[i].chromId, src = r[i].srcId; int64_t start = r[i].start, end = r[i].end; double score = r[i].score;
        int64_t w = start / win;
        if (have && ch == pch && src == psrc && w == pwin) {
            total += score * (double) (end - start);   /* additive -> associative */
            if (end > pend) pend = end;                 /* start monotone (sorted); end may not be */
        } else {
            if (have) EMIT();                            /* out < i here -- r[out] already consumed */
            pch = ch; psrc = src; pwin = w; pstart = start; pend = end; total = score * (double) (end - start); have = 1;
        }
    }
    if (have) EMIT();
#undef EMIT
    return out;
}

/* A reference's records, ALL chroms interleaved in ONE array (chromId carried in
 * each record).  Self-compacting: sort by (chromId,srcId,start) + window-bin in
 * place, realloc'd DOWN to the merged size.  ONE allocation per reference (not per
 * chrom): releasing the raw is a single mmap shrink, with none of the per-chrom
 * grow/shrink churn that otherwise balloons the allocator to tens of GB. */
typedef struct {
    PackedRec *recs; int64_t n, cap;
    int64_t    baseline;                     /* n right after the last compaction */
} RefAcc;

static void refacc_append(RefAcc *a, uint32_t chromId, int64_t start, int64_t end,
                          uint32_t srcId, double score) {
    if (a->n == a->cap) {
        a->cap = a->cap ? a->cap * 2 : 1024;
        a->recs = st_realloc(a->recs, (size_t) a->cap * sizeof(PackedRec));
    }
    a->recs[a->n].chromId = chromId; a->recs[a->n].start = start; a->recs[a->n].end = end;
    a->recs[a->n].srcId = srcId; a->recs[a->n].score = score;
    a->n++;
}

/* Sort by (chromId,srcId,start) + window-bin IN PLACE, keeping the capacity.  n
 * drops to the merged count; later appends refill the same buffer.  Do NOT realloc
 * down + regrow: that churns the allocator (every freed chunk is kept in the arena,
 * ballooning RSS to ~25x the working set).  So the one array grows to ~floorRecs
 * ONCE and stays flat -- RSS tracks the compaction threshold, not the raw. */
static void refacc_compact(RefAcc *a) {
    if (a->n <= 1) { a->baseline = a->n; return; }
    qsort(a->recs, (size_t) a->n, sizeof(PackedRec), cmp_compact);
    a->n = region_merge_inplace(a->recs, a->n);
    a->baseline = a->n;
}

/* Amortized-doubling trigger: compact once the reference has grown to 2x its last
 * merged size, so peak per reference stays ~2x its merged (bed) footprint -- never
 * the raw.  The floor (a slice of --mergeMem, the total in-RAM budget) amortizes the
 * sort.  Total peak RSS ~= 2x the resident references' merged size; tune via -M /
 * --refSubset (fewer refs/pass), not this. */
static void refacc_maybe_compact(RefAcc *a) {
    int64_t floorRecs = mergeMemBudget / (64 * (int64_t) sizeof(PackedRec));
    if (floorRecs < (1 << 20)) floorRecs = 1 << 20;       /* >= ~1M records (~40 MB) */
    if (a->n >= floorRecs && a->n >= 2 * a->baseline)
        refacc_compact(a);
}

static void refacc_free(RefAcc *a) {
    free(a->recs);
    a->recs = NULL; a->n = a->cap = a->baseline = 0;
}

/* Close a stream we WROTE to, aborting on any accumulated write error (ferror)
 * or a flush/close failure -- a silent short write (ENOSPC, quota) must never
 * pass as a valid-but-truncated output. */
static void wclose(FILE *f, const char *what) {
    int err = ferror(f);
    if (fclose(f) != 0) err = 1;
    if (err) { fprintf(stderr, "taffy summary: write error on %s\n", what); exit(1); }
}

/* ===================================================================== */
/* per-block scoring (streams packed records into per-chrom buckets)      */
/* ===================================================================== */
typedef struct {
    GerpTree *gt; const char *ref; size_t refLen;
    Interner *chromInt, *srcInt;
    RefAcc *acc;
    int64_t raw_n, n_with_ref, n_split_skipped;
} ScoreCtx;

static void score_block(ScoreCtx *cx, Alignment *aln) {
    int single_cover = (chainOverlapFrac < 1.0);
    for (Alignment_Row *ref = aln->row; ref != NULL; ref = ref->n_row) {
        const char *rg = gerp_tree_resolve_genome(cx->gt, ref->sequence_name);
        if (rg == NULL || strcmp(rg, cx->ref) != 0) continue;       /* only reference rows are masters */
        if (ref->length < minSeqSize || ref->length <= 0 || ref->sequence_length < minSeqSize) continue;
        cx->n_with_ref++;
        int64_t rstart = ref->strand ? ref->start
            : ref->sequence_length - (ref->start + ref->length);
        int64_t rend = rstart + ref->length;
        const char *chrom = (strlen(ref->sequence_name) > cx->refLen + 1)
            ? ref->sequence_name + cx->refLen + 1 : ref->sequence_name;
        uint32_t chromId = intern_id(cx->chromInt, chrom);
        if (ref->length > maxSize) cx->n_split_skipped++;

        /* best score per leaf species (single-cover the species in-block) */
        stHash *best = single_cover
            ? stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, free) : NULL;
        for (Alignment_Row *row = aln->row; row != NULL; row = row->n_row) {
            if (row == ref || row->length == 0) continue;
            const char *g = gerp_tree_resolve_genome(cx->gt, row->sequence_name);
            if (g == NULL || strcmp(g, cx->ref) == 0) continue;     /* unknown genome or reference paralog */
            if (gerp_tree_is_ancestor(cx->gt, g)) continue;         /* ancestor */
            double sc = score_pairwise(ref->bases, row->bases, aln->column_number, ref->length);
            uint32_t srcId = intern_id(cx->srcInt, g);
            if (single_cover) {
                const char *gname = intern_name(cx->srcInt, srcId);  /* stable key */
                double *cur = stHash_search(best, (void *) gname);
                if (cur == NULL) { cur = st_malloc(sizeof(double)); *cur = sc; stHash_insert(best, (void *) gname, cur); }
                else if (sc > *cur) *cur = sc;
            } else {
                refacc_append(cx->acc, chromId, rstart, rend, srcId, sc); cx->raw_n++;
            }
        }
        if (single_cover) {
            stHashIterator *hit = stHash_getIterator(best);
            char *gname;
            while ((gname = stHash_getNext(hit)) != NULL) {
                double *cur = stHash_search(best, gname);
                uint32_t srcId = intern_id(cx->srcInt, gname);
                refacc_append(cx->acc, chromId, rstart, rend, srcId, *cur); cx->raw_n++;
            }
            stHash_destructIterator(hit);
            stHash_destruct(best);
        }
        refacc_maybe_compact(cx->acc);
    }
}

/* ===================================================================== */
/* fan-out (--allRefs): every leaf genome is a master in ONE scan          */
/* ===================================================================== */
/* Each reference accumulates its scored records in RAM (per-chrom), self-compacting
 * (sort+region-merge) so its footprint tracks the merged bed size, never the raw --
 * no scan-temp, no disk spill.  In --allRefs the OpenMP scan partitions references
 * across threads (srcId%nt), so each thread owns its refs' RefAcc lock-free, and
 * compaction is per-reference so it stays lock-free. */
typedef struct {
    char     *name;        /* genome name (owned) */
    size_t    nameLen;
    Interner *chromInt;    /* this reference's chroms (chromId carried in the records) */
    RefAcc    acc;         /* in-RAM per-chrom records, compacted on the fly */
    int64_t   raw_n, n_with_ref, n_split_skipped;
} RefState;

typedef struct {
    GerpTree *gt;
    Interner *srcInt;      /* shared across references */
    stHash   *refMap;      /* genome name -> RefState* (key aliases ->name) */
    stList   *refList;     /* RefStates in first-seen order */
    const char *tmpDir;
    /* --refSubset j/M (with --allRefs): master only leaves whose name-hash %M==j,
     * so the file is read M times and each job masters ~1/M of the references
     * (node-local raw stays ~1/M, no shared scratch).  mod<=1 = master all. */
    int          refSubsetMod, refSubsetIdx;
} AllRefsCtx;

static RefState *get_refstate(AllRefsCtx *cx, const char *g) {
    RefState *rs = stHash_search(cx->refMap, (void *) g);
    if (rs == NULL) {
        rs = st_calloc(1, sizeof(RefState));   /* acc zeroed: chrom=NULL, live=0 */
        rs->name = stString_copy(g);
        rs->nameLen = strlen(g);
        rs->chromInt = interner_new();
        stHash_insert(cx->refMap, rs->name, rs);   /* key aliases rs->name */
        stList_append(cx->refList, rs);
    }
    return rs;
}


/* djb2 string hash -- deterministic across processes so every --refSubset job
 * partitions the reference roster identically (no shared roster needed). */
static uint64_t name_hash(const char *s) {
    uint64_t h = 5381;
    for (; *s; s++) h = h * 33u + (unsigned char) *s;
    return h;
}

/* Write <dir>/<ref>.chrom.sizes from the .tui sequence lengths (genome prefix
 * stripped to match the bed's chrom names) for bedToBigBed.  No-op if tui==NULL. */
static void write_chrom_sizes(Tui *tui, const char *ref, const char *dir) {
    if (tui == NULL) return;
    stHash *sl = tui_genome_seq_lengths(tui, ref);
    if (sl == NULL) { fprintf(stderr, "  WARNING: %s not in .tui roster -- no chrom.sizes\n", ref); return; }
    char path[PATH_MAX];
    snprintf(path, sizeof path, "%s/%s.chrom.sizes", dir, ref);
    FILE *cs = fopen(path, "w");
    if (cs == NULL) { fprintf(stderr, "taffy summary: cannot write %s\n", path); exit(1); }
    size_t glen = strlen(ref);
    stHashIterator *sit = stHash_getIterator(sl);
    char *fn;
    while ((fn = stHash_getNext(sit)) != NULL) {
        int64_t len = (int64_t) (intptr_t) stHash_search(sl, fn);
        const char *chrom = (strlen(fn) > glen + 1 && strncmp(fn, ref, glen) == 0 && fn[glen] == '.')
                            ? fn + glen + 1 : fn;
        fprintf(cs, "%s\t%" PRId64 "\n", chrom, len);
    }
    stHash_destructIterator(sit);
    wclose(cs, path);
    stHash_destruct(sl);
}

/* Mirror of score_block, but EVERY leaf row is a master routed to its own
 * RefState.  For a given reference the emitted record SET equals the single-ref
 * score_block(ref=that genome) set, so each per-reference bed matches a
 * single-ref run (compared as a sorted set). */
static void score_block_allrefs(AllRefsCtx *cx, Alignment *aln) {
    int single_cover = (chainOverlapFrac < 1.0);
    /* PRE-PASS (sequential): populate the SHARED interners.  intern_id is NOT
     * thread-safe, so the parallel scoring below must only READ them: srcInt
     * gets every leaf species, refMap/refList get every master genome. */
    int64_t n_leaf = 0, n_master = 0;
    for (Alignment_Row *r = aln->row; r != NULL; r = r->n_row) {
        const char *g = gerp_tree_resolve_genome(cx->gt, r->sequence_name);
        if (g == NULL || gerp_tree_is_ancestor(cx->gt, g)) continue;   /* leaf only */
        intern_id(cx->srcInt, g); n_leaf++;
        if (cx->refSubsetMod > 1 &&
            (int) (name_hash(g) % (uint64_t) cx->refSubsetMod) != cx->refSubsetIdx) continue;  /* not my subset */
        if (r->length < minSeqSize || r->length <= 0 || r->sequence_length < minSeqSize) continue;
        get_refstate(cx, g); n_master++;   /* one RefState per master genome */
    }
    /* PARALLEL scoring -- only when the block is heavy enough to outweigh the
     * OpenMP fork/join (the dense root/backbone blocks dominate; the many small
     * blocks stay sequential).  Thread t owns the reference genomes with interned
     * srcId %% nt == t: it scores ALL master rows of its genomes and appends ONLY
     * into those genomes' RefAcc -> zero contention, and the shared interners/refMap
     * are read-only here.  (Each RefAcc merges its own records -> -T8 == -T1 as a
     * set; the compaction boundary is browser-equivalent, not byte-identical.) */
    int nt = nThreads > 0 ? nThreads : omp_get_max_threads();
    if (nt < 1) nt = 1;
    /* per-block work below which the OpenMP fork/join isn't worth it; tunable
     * (env override, parsed once -- this fn runs on the single scan thread). */
    static int64_t scan_par_min = -1;
    if (scan_par_min < 0) {
        const char *e = getenv("TAFFY_SUMMARY_SCAN_PAR_MIN");
        scan_par_min = (e != NULL) ? atoll(e) : 100000;
    }
    int do_par = (nt > 1 && n_master >= 2 &&
                  (int64_t) n_master * n_leaf * aln->column_number >= scan_par_min);
    #pragma omp parallel num_threads(nt) if(do_par)
    {
        int tid = omp_get_thread_num();
        int cur_nt = omp_get_num_threads();   /* 1 when !do_par */
    for (Alignment_Row *ref = aln->row; ref != NULL; ref = ref->n_row) {
        const char *rg = gerp_tree_resolve_genome(cx->gt, ref->sequence_name);
        if (rg == NULL || gerp_tree_is_ancestor(cx->gt, rg)) continue;   /* leaf masters only */
        if (cx->refSubsetMod > 1 &&
            (int) (name_hash(rg) % (uint64_t) cx->refSubsetMod) != cx->refSubsetIdx) continue;  /* not my subset */
        if (ref->length < minSeqSize || ref->length <= 0 || ref->sequence_length < minSeqSize) continue;
        if (cur_nt > 1 &&
            (int)(intern_id(cx->srcInt, rg) % (uint32_t) cur_nt) != tid) continue;  /* not my genome */
        RefState *rs = get_refstate(cx, rg);
        rs->n_with_ref++;
        int64_t rstart = ref->strand ? ref->start
            : ref->sequence_length - (ref->start + ref->length);
        int64_t rend = rstart + ref->length;
        const char *chrom = (strlen(ref->sequence_name) > rs->nameLen + 1)
            ? ref->sequence_name + rs->nameLen + 1 : ref->sequence_name;
        uint32_t chromId = intern_id(rs->chromInt, chrom);   /* per-genome interner: this genome's sole thread */
        if (ref->length > maxSize) rs->n_split_skipped++;

        stHash *best = single_cover
            ? stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, free) : NULL;
        for (Alignment_Row *row = aln->row; row != NULL; row = row->n_row) {
            if (row == ref || row->length == 0) continue;
            const char *g = gerp_tree_resolve_genome(cx->gt, row->sequence_name);
            if (g == NULL || strcmp(g, rg) == 0) continue;     /* unknown genome or reference paralog */
            if (gerp_tree_is_ancestor(cx->gt, g)) continue;    /* ancestor */
            double sc = score_pairwise(ref->bases, row->bases, aln->column_number, ref->length);
            uint32_t srcId = intern_id(cx->srcInt, g);          /* read-only (pre-interned) */
            if (single_cover) {
                const char *gname = intern_name(cx->srcInt, srcId);  /* stable key */
                double *cur = stHash_search(best, (void *) gname);
                if (cur == NULL) { cur = st_malloc(sizeof(double)); *cur = sc; stHash_insert(best, (void *) gname, cur); }
                else if (sc > *cur) *cur = sc;
            } else {
                refacc_append(&rs->acc, chromId, rstart, rend, srcId, sc); rs->raw_n++;
            }
        }
        if (single_cover) {
            stHashIterator *hit = stHash_getIterator(best);
            char *gname;
            while ((gname = stHash_getNext(hit)) != NULL) {
                double *cur = stHash_search(best, gname);
                uint32_t srcId = intern_id(cx->srcInt, gname);   /* read-only (pre-interned) */
                refacc_append(&rs->acc, chromId, rstart, rend, srcId, *cur); rs->raw_n++;
            }
            stHash_destructIterator(hit);
            stHash_destruct(best);
        }
        refacc_maybe_compact(&rs->acc);
    }
    }
}

/* ===================================================================== */
/* finalize: final compaction + write the reference's bed3+4              */
/* ===================================================================== */
static Interner *g_chrom_for_sort;
static int cmp_chrom_name(const void *a, const void *b) {   /* chroms lexicographically (bedToBigBed wants -k1,1) */
    return strcmp(intern_name(g_chrom_for_sort, *(const uint32_t *) a),
                  intern_name(g_chrom_for_sort, *(const uint32_t *) b));
}

/* Final window-bin, then write bed3+4 in (chrom,start)-sorted order (chroms
 * lexicographic, then start, srcId).  Returns the row count.  Shared by the -r path
 * and each --allRefs reference (called sequentially per reference, so the
 * g_chrom_for_sort / g_lexrank globals are set race-free). */
static int64_t refacc_finalize(RefAcc *a, Interner *chromInt, Interner *srcInt, FILE *out) {
    refacc_compact(a);
    if (a->n == 0) return 0;
    /* chrom output order = lexicographic by name -> chromId -> lex rank */
    int nChrom = (int) chromInt->n;
    if (nChrom < 1) nChrom = 1;
    uint32_t *order = st_malloc((size_t) nChrom * sizeof(uint32_t));
    for (int i = 0; i < nChrom; i++) order[i] = (uint32_t) i;
    g_chrom_for_sort = chromInt;
    qsort(order, (size_t) nChrom, sizeof(uint32_t), cmp_chrom_name);
    int64_t *lex = st_malloc((size_t) nChrom * sizeof(int64_t));
    for (int i = 0; i < nChrom; i++) lex[order[i]] = i;
    g_lexrank = lex;
    qsort(a->recs, (size_t) a->n, sizeof(PackedRec), cmp_output);
    for (int64_t k = 0; k < a->n; k++)
        fprintf(out, "%s\t%" PRId64 "\t%" PRId64 "\t%s\t%g\t\t\n",
                intern_name(chromInt, a->recs[k].chromId), a->recs[k].start, a->recs[k].end,
                intern_name(srcInt, a->recs[k].srcId), a->recs[k].score);
    free(order); free(lex);
    return a->n;
}

int taf_summary_main(int argc, char *argv[]) {
    char *inputFile = NULL, *outputFile = NULL, *refGenome = NULL, *tmpDir = NULL;
    int allRefs = 0;
    int refSubsetIdx = 0, refSubsetMod = 1;              /* --refSubset j/M: per-reference fan-out */

    static struct option lopts[] = {
        {"mergeGap",        required_argument, 0, 'g'},
        {"minSize",         required_argument, 0, 's'},
        {"maxSize",         required_argument, 0, 'S'},
        {"minSeqSize",      required_argument, 0, 'q'},
        {"chainOverlapFrac",required_argument, 0, 'F'},
        {"tmp",             required_argument, 0, 'T'},
        {"threads",         required_argument, 0, 'p'},
        {"mergeMem",        required_argument, 0, 1005},
        {"allRefs",         no_argument,       0, 'A'},
        {"refSubset",       required_argument, 0, 1006},
        {"help",            no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "i:r:o:g:s:S:q:F:T:p:Ah", lopts, NULL)) != -1) {
        switch (c) {
            case 'i': inputFile  = optarg; break;
            case 'r': refGenome  = optarg; break;
            case 'o': outputFile = optarg; break;
            case 'A': allRefs    = 1; break;
            case 1006:   /* --refSubset j/M: master only the refs with name-hash %M==j */
                if (sscanf(optarg, "%d/%d", &refSubsetIdx, &refSubsetMod) != 2 || refSubsetMod <= 0 || refSubsetIdx < 0 || refSubsetIdx >= refSubsetMod) {
                    fprintf(stderr, "taffy summary: --refSubset must be j/M with 0<=j<M\n"); return 1;
                }
                break;
            case 'g': mergeGap   = atoi(optarg); break;
            case 's': minSize    = atoi(optarg); break;
            case 'S': maxSize    = atoi(optarg); break;
            case 'q': minSeqSize = atoll(optarg); break;
            case 'F': chainOverlapFrac = atof(optarg); break;
            case 'T': tmpDir     = optarg; break;
            case 'p': nThreads   = atoi(optarg); break;
            case 1005: {   /* --mergeMem G : per-chrom merge RAM budget in GiB */
                double g = atof(optarg);
                if (g <= 0) { fprintf(stderr, "taffy summary: --mergeMem must be > 0 (GiB)\n"); return 1; }
                mergeMemBudget = (int64_t)(g * (double)(1LL << 30));
                break;
            }
            case 'h':
            default:
                fprintf(stderr,
                  "usage: taffy summary -i <universal .taf/.maf[.gz]> -r <refGenome> [-o out.bed]\n"
                  "  Emit a per-reference bigMafSummary (bed3+4) from a universal alignment.\n"
                  "  Output is (chrom,start)-sorted; feed straight to bedToBigBed -type=bed3+4.\n"
                  "  -i FILE   universal MAF/TAF (its `# hal` header tree names leaves vs ancestors)\n"
                  "  -r NAME   reference genome (the 'master'); its leaf neighbours are scored\n"
                  "  -o FILE   output bed (default stdout); with --allRefs, an output DIRECTORY\n"
                  "  --allRefs  ONE scan emits a bed per LEAF genome to <-o DIR>/<genome>.bed\n"
                  "            (every leaf is a master; -r is ignored)\n"
                  "  --refSubset j/M  with --allRefs: master only 1/M of the references (by\n"
                  "            name-hash); the file is read M times -> per-reference fan-out with no\n"
                  "            shared scratch (raw stays node-local). Also writes <ref>.chrom.sizes\n"
                  "  --tmp DIR   accepted for compatibility; the merge is now fully IN-RAM\n"
                  "            (self-compacting per reference) -- no scratch files are written\n"
                  "  --threads N (ncores)  parallel scan scoring (--allRefs) + final per-chrom merge\n"
                  "  --mergeMem G (8)  in-RAM compaction budget (GiB): a reference compacts\n"
                  "            (merges to runs) as it grows, bounding peak to ~2x its merged size\n"
                  "  --chainOverlapFrac F (%.2f)  per-species in-block dedup: 0=single-cover\n"
                  "            (best paralog row, matches reference-anchored summary), 1=keep all\n"
                  "  --mergeGap N (%d)  --minSize N (%d)  --maxSize N (%d)  --minSeqSize N (%" PRId64 ")\n",
                  chainOverlapFrac, mergeGap, minSize, maxSize, minSeqSize);
                return (c == 'h') ? 0 : 1;
        }
    }
    if (inputFile == NULL) {
        fprintf(stderr, "taffy summary: -i is required\n");
        return 1;
    }
    if (allRefs) {
        if (outputFile == NULL) {
            fprintf(stderr, "taffy summary: --allRefs requires -o <output directory>\n");
            return 1;
        }
    } else if (refGenome == NULL) {
        fprintf(stderr, "taffy summary: -r is required (or use --allRefs)\n");
        return 1;
    }
    if (refSubsetMod > 1 && !allRefs) {
        fprintf(stderr, "taffy summary: --refSubset requires --allRefs\n"); return 1;
    }
    if (chainOverlapFrac < 0.0 || chainOverlapFrac > 1.0) {
        fprintf(stderr, "taffy summary: --chainOverlapFrac must be in [0,1]\n");
        return 1;
    }
    if (tmpDir == NULL) { tmpDir = getenv("TMPDIR"); if (tmpDir == NULL) tmpDir = "/tmp"; }
    if (nThreads <= 0) { long nc = sysconf(_SC_NPROCESSORS_ONLN); nThreads = (nc > 0) ? (int) nc : 1; }

    scorer_init();

    /* Single-threaded bgzf decode ON PURPOSE: htslib's multithreaded bgzf_mt
     * (what LI_set_bgzf_threads(>1) enables) leaks ~unbounded on a large file
     * (measured 6.4 GB on apes at -p 8 vs 246 MB at -p 1 -- catastrophic at
     * 577-way).  The scan's real parallelism is the OpenMP scoring (--threads),
     * not decode.  If decode becomes the bottleneck, fix bgzf_mt teardown then. */
    LI_set_bgzf_threads(1);

    FILE *in_fh = fopen(inputFile, "r");
    if (in_fh == NULL) { fprintf(stderr, "taffy summary: cannot open %s\n", inputFile); return 1; }
    LI *li = LI_construct(in_fh);
    int input_format = check_input_format(LI_peek_at_next_line(li));
    if (input_format != 0 && input_format != 1) {
        fprintf(stderr, "taffy summary: input must be MAF or TAF\n"); return 1;
    }
    bool rle = false;
    Tag *header_tag = (input_format == 0) ? taf_read_header_2(li, &rle) : maf_read_header(li);

    Tag *hal = tag_find(header_tag, (char *) TAF_HAL_TREE_KEY);
    if (hal == NULL) {
        fprintf(stderr, "taffy summary: input header has no `# hal` tree (need it to "
                        "tell leaves from ancestors)\n");
        return 1;
    }
    GerpTree *gt = gerp_tree_construct(hal->value);
    if (gt == NULL) { fprintf(stderr, "taffy summary: failed to parse `# hal` tree\n"); return 1; }
    if (!allRefs && gerp_tree_is_ancestor(gt, refGenome)) {
        fprintf(stderr, "taffy summary: -r '%s' is an internal (ancestor) node, not a genome\n", refGenome);
        return 1;
    }

    Interner *srcInt = interner_new();   /* species ids -- shared across references */

    if (allRefs) {
        /* -------- fan-out: ONE scan, every leaf genome is a master accumulating
         * its records in RAM (self-compacting); then write each to <outDir>/
         * <genome>.bed.  Each per-reference bed equals the single-ref `-r` run as a
         * set (compaction boundaries are browser-equivalent, not byte-identical). */
        if (mkdir(outputFile, 0777) != 0 && errno != EEXIST) {
            fprintf(stderr, "taffy summary: cannot create output dir %s\n", outputFile); return 1;
        }
        stHash *refMap = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
        stList *refList = stList_construct();
        AllRefsCtx acx = { gt, srcInt, refMap, refList, tmpDir };
        acx.refSubsetMod = refSubsetMod; acx.refSubsetIdx = refSubsetIdx;
        Alignment *prev = NULL;
        int64_t n_blocks = 0;
        while (1) {
            Alignment *aln = (input_format == 0) ? taf_read_block(prev, rle, li) : maf_read_block(li);
            if (prev != NULL) { alignment_destruct(prev, 1); prev = NULL; }
            if (aln == NULL) break;
            n_blocks++;
            score_block_allrefs(&acx, aln);
            prev = aln;
        }
        stHash_destruct(refMap);   /* lookups done; RefStates live on in refList */

        int nRefs = stList_length(refList);
        int64_t total_out = 0;
        char path[PATH_MAX];
        Tui *atui = NULL;   /* <inputFile>.tui -> per-ref chrom.sizes for bedToBigBed */
        { char tp[PATH_MAX]; snprintf(tp, sizeof tp, "%s.tui", inputFile); atui = tui_load(tp);
          if (atui == NULL) fprintf(stderr, "taffy summary --allRefs: no %s -- no chrom.sizes written\n", tp); }
        for (int i = 0; i < nRefs; i++) {
            RefState *rs = stList_get(refList, i);
            snprintf(path, sizeof path, "%s/%s.bed", outputFile, rs->name);
            FILE *rout = fopen(path, "w");
            if (rout == NULL) { fprintf(stderr, "taffy summary: cannot write %s\n", path); return 1; }
            int64_t no = refacc_finalize(&rs->acc, rs->chromInt, srcInt, rout);
            wclose(rout, path);
            write_chrom_sizes(atui, rs->name, outputFile);
            total_out += no;
            fprintf(stderr, "  %s: %" PRId64 " master rows, %" PRId64 " raw recs -> %" PRId64 " merged rows\n",
                    rs->name, rs->n_with_ref, rs->raw_n, no);
            refacc_free(&rs->acc);
            interner_free(rs->chromInt); free(rs->name); free(rs);
        }
        if (atui != NULL) tui_destruct(atui);
        fprintf(stderr, "taffy summary --allRefs: %d references, %" PRId64 " blocks, "
                "%" PRId64 " total merged rows (single-cover F=%.2f, %d threads)\n",
                nRefs, n_blocks, total_out, chainOverlapFrac, nThreads);
        stList_destruct(refList);
    } else {
        FILE *out = (outputFile == NULL) ? stdout : fopen(outputFile, "w");
        if (out == NULL) { fprintf(stderr, "taffy summary: cannot write %s\n", outputFile); return 1; }

        /* -------- single reference: stream blocks, score, single-cover species
         * in-block, accumulate + self-compact in RAM.  Keep only the previous block
         * alive (TAF deltas reference it); free it on the next read. */
        Interner *chromInt = interner_new();
        RefAcc acc = { NULL, 0, 0, 0 };
        ScoreCtx cx = { gt, refGenome, strlen(refGenome), chromInt, srcInt, &acc, 0, 0, 0 };
        Alignment *prev = NULL;
        int64_t n_blocks = 0;
        while (1) {
            Alignment *aln = (input_format == 0) ? taf_read_block(prev, rle, li) : maf_read_block(li);
            if (prev != NULL) { alignment_destruct(prev, 1); prev = NULL; }
            if (aln == NULL) break;
            n_blocks++;
            score_block(&cx, aln);
            prev = aln;
        }
        int64_t raw_n = cx.raw_n, n_with_ref = cx.n_with_ref, n_split_skipped = cx.n_split_skipped;

        int64_t n_out = refacc_finalize(&acc, chromInt, srcInt, out);
        if (out != stdout) wclose(out, outputFile);
        else if (fflush(out) != 0 || ferror(out)) { fprintf(stderr, "taffy summary: write error on stdout\n"); exit(1); }

        int nt = nThreads; if (nt < 1) nt = 1;
        fprintf(stderr, "taffy summary: %" PRId64 " blocks (%" PRId64 " reference master rows, "
                "single-cover F=%.2f), %" PRId64 " raw records -> %" PRId64 " merged summary rows "
                "(%d chroms, %d threads)\n",
                n_blocks, n_with_ref, chainOverlapFrac, raw_n, n_out, (int) chromInt->n, nt);
        if (n_split_skipped > 0)
            fprintf(stderr, "taffy summary: NOTE: %" PRId64 " reference master row(s) had span > "
                    "maxSize=%d and were NOT split (cosmetic; universal-MAF blocks are normally small)\n",
                    n_split_skipped, maxSize);

        refacc_free(&acc);
        interner_free(chromInt);
    }

    interner_free(srcInt);
    gerp_tree_destruct(gt);
    LI_destruct(li);
    if (header_tag) tag_destruct(header_tag);
    return 0;
}
