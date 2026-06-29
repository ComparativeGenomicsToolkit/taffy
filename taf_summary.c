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
 * consecutive reference rows), so a global reorder is unavoidable and merging in
 * scan order does not work.  But the reorder need not be one giant external sort
 * (~2.9 billion records / ~200 GB temp at 577-way).  Instead we BUCKET BY
 * REFERENCE CHROM: each scored record is written, as a fixed-size packed binary
 * struct, to a per-reference-chrom temp file.  After the scan each chrom's bucket
 * fits in RAM, so we sort+merge it in-memory (no external sort), in PARALLEL
 * across chroms, then concatenate the per-chrom outputs in chrom order.  Peak RAM
 * is bounded by (largest chrom's records) x threads, not the total; disk is the
 * packed record volume (interned int ids + packed structs ~= 1/3 the old text).
 *
 * The score is stored as a DOUBLE in the packed record (not float): the region
 * merge length-weighted-averages scores, and float would diverge from the
 * external-sort reference output under %g formatting.  Output is byte-identical
 * (as a record SET) to the external-sort version; row ORDER differs (we emit
 * (chrom,start)-sorted, which is bedToBigBed-ready -- no post-sort needed).
 *
 * Output is bed3+4 (chrom start end src score leftStatus rightStatus), the same
 * schema as kent's mafSummary.as.  Feed straight to
 * `bedToBigBed -type=bed3+4 -as=mafSummary.as -tab`.
 *
 *  Released under the MIT license, see LICENSE.txt
 */

#define _GNU_SOURCE   /* mkstemp/fdopen */
#include "taf.h"
#include "gerp.h"
#include "line_iterator.h"
#include "sonLib.h"
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
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
typedef struct { uint32_t chromId, start, end, srcId; double score; } PackedRec;  /* 24B, score 8-aligned */

typedef struct {
    FILE   *wfp;                 /* write handle during the scan (lazy) */
    char    in_path[PATH_MAX];   /* packed bucket temp */
    char    out_path[PATH_MAX];  /* merged bed temp (set by worker; empty if none) */
    int64_t n_recs;              /* raw records written */
    int64_t merged_n;            /* merged records (set by worker) */
} Bucket;

typedef struct { Bucket *arr; int n, cap; } BucketSet;

static Bucket *get_bucket(BucketSet *bs, uint32_t id, const char *tmpDir) {
    if ((int) id >= bs->cap) {
        int nc = bs->cap ? bs->cap : 64;
        while (nc <= (int) id) nc *= 2;
        bs->arr = st_realloc(bs->arr, (size_t) nc * sizeof(Bucket));
        for (int i = bs->cap; i < nc; i++) {
            bs->arr[i].wfp = NULL; bs->arr[i].n_recs = 0; bs->arr[i].merged_n = 0;
            bs->arr[i].in_path[0] = 0; bs->arr[i].out_path[0] = 0;
        }
        bs->cap = nc;
    }
    if ((int) id >= bs->n) bs->n = (int) id + 1;
    Bucket *b = &bs->arr[id];
    if (b->wfp == NULL) {
        snprintf(b->in_path, PATH_MAX, "%s/taffy_sum_bk_XXXXXX", tmpDir);
        int fd = mkstemp(b->in_path);
        if (fd < 0) { fprintf(stderr, "taffy summary: cannot create bucket temp in %s (raise ulimit -n?)\n", tmpDir); exit(1); }
        b->wfp = fdopen(fd, "wb");
        if (b->wfp == NULL) { fprintf(stderr, "taffy summary: fdopen bucket failed\n"); exit(1); }
    }
    return b;
}

static inline void put_rec(Bucket *b, uint32_t chromId, uint32_t start, uint32_t end, uint32_t srcId, double score) {
    PackedRec r = { chromId, start, end, srcId, score };
    fwrite(&r, sizeof r, 1, b->wfp);
    b->n_recs++;
}

/* ===================================================================== */
/* per-block scoring (streams packed records into per-chrom buckets)      */
/* ===================================================================== */
typedef struct {
    GerpTree *gt; const char *ref; size_t refLen;
    Interner *chromInt, *srcInt;
    BucketSet *bs; const char *tmpDir;
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
        Bucket *bk = get_bucket(cx->bs, chromId, cx->tmpDir);
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
                put_rec(bk, chromId, (uint32_t) rstart, (uint32_t) rend, srcId, sc); cx->raw_n++;
            }
        }
        if (single_cover) {
            stHashIterator *hit = stHash_getIterator(best);
            char *gname;
            while ((gname = stHash_getNext(hit)) != NULL) {
                double *cur = stHash_search(best, gname);
                uint32_t srcId = intern_id(cx->srcInt, gname);
                put_rec(bk, chromId, (uint32_t) rstart, (uint32_t) rend, srcId, *cur); cx->raw_n++;
            }
            stHash_destructIterator(hit);
            stHash_destruct(best);
        }
    }
}

/* ===================================================================== */
/* fan-out (--allRefs): every leaf genome is a master in ONE scan          */
/* ===================================================================== */
/* During the scan each reference keeps ONE open write stream (all its chroms
 * mixed; chromId is carried in the packed record) -- N_refs streams total, NOT
 * per-(ref,chrom) (that opened tens of thousands of buffered streams across all
 * refs at once -> OOM).  Post-processing re-buckets ONE reference at a time, so
 * the per-chrom fan-out (and its memory) is bounded to a single reference. */
typedef struct {
    char     *name;        /* genome name (owned) */
    size_t    nameLen;
    Interner *chromInt;    /* this reference's chroms (chromId carried in packed recs) */
    FILE     *wfp;         /* ONE scan-temp write stream (lazy) */
    char      tmp_path[PATH_MAX];
    int64_t   n_recs;      /* packed records written to wfp */
    int64_t   raw_n, n_with_ref, n_split_skipped;
} RefState;

typedef struct {
    GerpTree *gt;
    Interner *srcInt;      /* shared across references */
    stHash   *refMap;      /* genome name -> RefState* (key aliases ->name) */
    stList   *refList;     /* RefStates in first-seen order */
    const char *tmpDir;
} AllRefsCtx;

static RefState *get_refstate(AllRefsCtx *cx, const char *g) {
    RefState *rs = stHash_search(cx->refMap, (void *) g);
    if (rs == NULL) {
        rs = st_calloc(1, sizeof(RefState));   /* wfp=NULL, n_recs=0, tmp_path[0]=0 */
        rs->name = stString_copy(g);
        rs->nameLen = strlen(g);
        rs->chromInt = interner_new();
        stHash_insert(cx->refMap, rs->name, rs);   /* key aliases rs->name */
        stList_append(cx->refList, rs);
    }
    return rs;
}

/* append one packed record to this reference's single scan-temp (lazy open) */
static inline void put_ref_rec(RefState *rs, const char *tmpDir, uint32_t chromId,
                               uint32_t start, uint32_t end, uint32_t srcId, double score) {
    if (rs->wfp == NULL) {
        snprintf(rs->tmp_path, PATH_MAX, "%s/taffy_sum_ref_XXXXXX", tmpDir);
        int fd = mkstemp(rs->tmp_path);
        if (fd < 0) { fprintf(stderr, "taffy summary: cannot create per-reference temp in %s\n", tmpDir); exit(1); }
        rs->wfp = fdopen(fd, "wb");
        if (rs->wfp == NULL) { fprintf(stderr, "taffy summary: fdopen per-reference temp failed\n"); exit(1); }
    }
    PackedRec r = { chromId, start, end, srcId, score };
    fwrite(&r, sizeof r, 1, rs->wfp);
    rs->n_recs++;
}

/* Mirror of score_block, but EVERY leaf row is a master routed to its own
 * RefState.  For a given reference the emitted record SET equals the single-ref
 * score_block(ref=that genome) set, so each per-reference bed matches a
 * single-ref run (compared as a sorted set). */
static void score_block_allrefs(AllRefsCtx *cx, Alignment *aln) {
    int single_cover = (chainOverlapFrac < 1.0);
    for (Alignment_Row *ref = aln->row; ref != NULL; ref = ref->n_row) {
        const char *rg = gerp_tree_resolve_genome(cx->gt, ref->sequence_name);
        if (rg == NULL || gerp_tree_is_ancestor(cx->gt, rg)) continue;   /* leaf masters only */
        if (ref->length < minSeqSize || ref->length <= 0 || ref->sequence_length < minSeqSize) continue;
        RefState *rs = get_refstate(cx, rg);
        rs->n_with_ref++;
        int64_t rstart = ref->strand ? ref->start
            : ref->sequence_length - (ref->start + ref->length);
        int64_t rend = rstart + ref->length;
        const char *chrom = (strlen(ref->sequence_name) > rs->nameLen + 1)
            ? ref->sequence_name + rs->nameLen + 1 : ref->sequence_name;
        uint32_t chromId = intern_id(rs->chromInt, chrom);
        if (ref->length > maxSize) rs->n_split_skipped++;

        stHash *best = single_cover
            ? stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, free) : NULL;
        for (Alignment_Row *row = aln->row; row != NULL; row = row->n_row) {
            if (row == ref || row->length == 0) continue;
            const char *g = gerp_tree_resolve_genome(cx->gt, row->sequence_name);
            if (g == NULL || strcmp(g, rg) == 0) continue;     /* unknown genome or reference paralog */
            if (gerp_tree_is_ancestor(cx->gt, g)) continue;    /* ancestor */
            double sc = score_pairwise(ref->bases, row->bases, aln->column_number, ref->length);
            uint32_t srcId = intern_id(cx->srcInt, g);
            if (single_cover) {
                const char *gname = intern_name(cx->srcInt, srcId);  /* stable key */
                double *cur = stHash_search(best, (void *) gname);
                if (cur == NULL) { cur = st_malloc(sizeof(double)); *cur = sc; stHash_insert(best, (void *) gname, cur); }
                else if (sc > *cur) *cur = sc;
            } else {
                put_ref_rec(rs, cx->tmpDir, chromId, (uint32_t) rstart, (uint32_t) rend, srcId, sc); rs->raw_n++;
            }
        }
        if (single_cover) {
            stHashIterator *hit = stHash_getIterator(best);
            char *gname;
            while ((gname = stHash_getNext(hit)) != NULL) {
                double *cur = stHash_search(best, gname);
                uint32_t srcId = intern_id(cx->srcInt, gname);
                put_ref_rec(rs, cx->tmpDir, chromId, (uint32_t) rstart, (uint32_t) rend, srcId, *cur); rs->raw_n++;
            }
            stHash_destructIterator(hit);
            stHash_destruct(best);
        }
    }
}

/* ===================================================================== */
/* per-chrom in-memory sort + region merge (run in parallel)              */
/* ===================================================================== */
typedef struct { uint32_t start, end, srcId; double score; } MergedRec;

static int cmp_packed(const void *a, const void *b) {   /* by (srcId, start) */
    const PackedRec *x = a, *y = b;
    if (x->srcId != y->srcId) return x->srcId < y->srcId ? -1 : 1;
    if (x->start != y->start) return x->start < y->start ? -1 : 1;
    return 0;
}
static int cmp_merged(const void *a, const void *b) {   /* by (start, srcId) */
    const MergedRec *x = a, *y = b;
    if (x->start != y->start) return x->start < y->start ? -1 : 1;
    if (x->srcId != y->srcId) return x->srcId < y->srcId ? -1 : 1;
    return 0;
}

typedef struct {
    BucketSet *bs; Interner *chromInt, *srcInt; const char *tmpDir;
    int next; pthread_mutex_t lock;
} WorkPool;

/* read bucket `cid` into RAM, sort by (srcId,start), kent region-merge per src,
 * sort merged by (start,srcId), write (chrom,start)-ordered bed3+4 to a temp. */
static void process_chrom(WorkPool *wp, int cid) {
    Bucket *b = &wp->bs->arr[cid];
    if (b->n_recs == 0) { b->merged_n = 0; return; }

    FILE *f = fopen(b->in_path, "rb");
    if (f == NULL) { fprintf(stderr, "taffy summary: cannot reopen bucket %s\n", b->in_path); exit(1); }
    int64_t n = b->n_recs;
    PackedRec *recs = st_malloc((size_t) n * sizeof(PackedRec));
    if (fread(recs, sizeof(PackedRec), (size_t) n, f) != (size_t) n) {
        fprintf(stderr, "taffy summary: short read on bucket %s\n", b->in_path); exit(1);
    }
    fclose(f);
    remove(b->in_path);
    qsort(recs, (size_t) n, sizeof(PackedRec), cmp_packed);

    MergedRec *m = NULL; int64_t mn = 0, mcap = 0;
#define PUSH(s,e,sr,sc) do { \
        if (mn == mcap) { mcap = mcap ? mcap * 2 : 1024; m = st_realloc(m, (size_t) mcap * sizeof(MergedRec)); } \
        m[mn].start = (s); m[mn].end = (e); m[mn].srcId = (sr); m[mn].score = (sc); mn++; \
    } while (0)
    int have = 0, first = 1; uint32_t psrc = 0, pstart = 0, pend = 0; double pscore = 0;
    for (int64_t i = 0; i < n; i++) {
        uint32_t src = recs[i].srcId, start = recs[i].start, end = recs[i].end; double score = recs[i].score;
        if (first || src != psrc) {
            if (have) { PUSH(pstart, pend, psrc, pscore); have = 0; }
            psrc = src; first = 0;
        }
        if (have && ((int64_t) start + 1 - (int64_t) pend < mergeGap)) {
            double total = pscore * (double)((int64_t) pend - (int64_t) pstart)
                         + score  * (double)((int64_t) end  - (int64_t) start);
            pend = end;
            pscore = total / (double)((int64_t) pend - (int64_t) pstart);
        } else {
            if (have) PUSH(pstart, pend, psrc, pscore);
            pstart = start; pend = end; pscore = score; have = 1;
        }
        if ((int64_t) pend - (int64_t) pstart > minSize) { PUSH(pstart, pend, psrc, pscore); have = 0; }
    }
    if (have) PUSH(pstart, pend, psrc, pscore);
#undef PUSH
    free(recs);

    qsort(m, (size_t) mn, sizeof(MergedRec), cmp_merged);

    snprintf(b->out_path, PATH_MAX, "%s/taffy_sum_out_XXXXXX", wp->tmpDir);
    int ofd = mkstemp(b->out_path);
    if (ofd < 0) { fprintf(stderr, "taffy summary: cannot create out temp in %s\n", wp->tmpDir); exit(1); }
    FILE *of = fdopen(ofd, "w");
    const char *chrom = intern_name(wp->chromInt, (uint32_t) cid);
    for (int64_t i = 0; i < mn; i++)
        fprintf(of, "%s\t%" PRIu32 "\t%" PRIu32 "\t%s\t%g\t\t\n",
                chrom, m[i].start, m[i].end, intern_name(wp->srcInt, m[i].srcId), m[i].score);
    fclose(of);
    free(m);
    b->merged_n = mn;
}

static void *worker(void *arg) {
    WorkPool *wp = arg;
    for (;;) {
        pthread_mutex_lock(&wp->lock);
        int cid = (wp->next < wp->bs->n) ? wp->next++ : -1;
        pthread_mutex_unlock(&wp->lock);
        if (cid < 0) break;
        process_chrom(wp, cid);
    }
    return NULL;
}

/* concat order: reference chroms lexicographically (bedToBigBed wants -k1,1 -k2,2n) */
static Interner *g_chrom_for_sort;
static int cmp_chrom_name(const void *a, const void *b) {
    return strcmp(intern_name(g_chrom_for_sort, *(const uint32_t *) a),
                  intern_name(g_chrom_for_sort, *(const uint32_t *) b));
}

/* flush+close buckets, parallel per-chrom sort+merge, then concat the per-chrom
 * outputs in chrom-name order to `out`.  Returns the merged row count.  Shared
 * by the single-ref path and each --allRefs reference (so a reference's bed is
 * produced identically either way). */
static int64_t finalize_reference(BucketSet *bs, Interner *chromInt, Interner *srcInt,
                                  const char *tmpDir, int nThreads, FILE *out) {
    for (int i = 0; i < bs->n; i++)
        if (bs->arr[i].wfp != NULL) {
            if (fclose(bs->arr[i].wfp) != 0) {
                fprintf(stderr, "taffy summary: error writing bucket %s\n", bs->arr[i].in_path); exit(1);
            }
            bs->arr[i].wfp = NULL;
        }

    WorkPool wp = { bs, chromInt, srcInt, tmpDir, 0, PTHREAD_MUTEX_INITIALIZER };
    int nt = nThreads; if (nt > bs->n) nt = bs->n; if (nt < 1) nt = 1;
    pthread_t *tids = st_malloc((size_t) nt * sizeof(pthread_t));
    for (int t = 0; t < nt; t++) pthread_create(&tids[t], NULL, worker, &wp);
    for (int t = 0; t < nt; t++) pthread_join(tids[t], NULL);
    free(tids);

    uint32_t *order = st_malloc((size_t) (bs->n > 0 ? bs->n : 1) * sizeof(uint32_t));
    for (int i = 0; i < bs->n; i++) order[i] = (uint32_t) i;
    g_chrom_for_sort = chromInt;
    qsort(order, (size_t) bs->n, sizeof(uint32_t), cmp_chrom_name);
    int64_t n_out = 0;
    char buf[1 << 16];
    for (int i = 0; i < bs->n; i++) {
        Bucket *b = &bs->arr[order[i]];
        if (b->merged_n == 0 || b->out_path[0] == 0) continue;
        FILE *cf = fopen(b->out_path, "r");
        if (cf == NULL) { fprintf(stderr, "taffy summary: cannot read out temp %s\n", b->out_path); exit(1); }
        size_t r;
        while ((r = fread(buf, 1, sizeof buf, cf)) > 0) fwrite(buf, 1, r, out);
        fclose(cf);
        remove(b->out_path);
        n_out += b->merged_n;
    }
    free(order);
    return n_out;
}

/* re-read a reference's single scan-temp, route its records into per-chrom
 * buckets (bounded to THIS reference's chroms), then finalize as usual.  Lets
 * --allRefs keep only N_refs streams open during the scan and one reference's
 * chrom fan-out during post -- exactly the single-ref path's footprint. */
static int64_t rebucket_and_finalize(RefState *rs, Interner *srcInt,
                                     const char *tmpDir, int nThreads, FILE *out) {
    if (rs->wfp != NULL) {
        if (fclose(rs->wfp) != 0) { fprintf(stderr, "taffy summary: error writing %s\n", rs->tmp_path); exit(1); }
        rs->wfp = NULL;
    }
    BucketSet bs = { NULL, 0, 0 };
    if (rs->n_recs > 0) {
        FILE *f = fopen(rs->tmp_path, "rb");
        if (f == NULL) { fprintf(stderr, "taffy summary: cannot reopen %s\n", rs->tmp_path); exit(1); }
        enum { RBUF = 1 << 16 };
        PackedRec *rbuf = st_malloc((size_t) RBUF * sizeof(PackedRec));
        size_t got;
        while ((got = fread(rbuf, sizeof(PackedRec), RBUF, f)) > 0)
            for (size_t i = 0; i < got; i++) {
                Bucket *bk = get_bucket(&bs, rbuf[i].chromId, tmpDir);
                put_rec(bk, rbuf[i].chromId, rbuf[i].start, rbuf[i].end, rbuf[i].srcId, rbuf[i].score);
            }
        free(rbuf);
        fclose(f);
        remove(rs->tmp_path);
    }
    int64_t n_out = finalize_reference(&bs, rs->chromInt, srcInt, tmpDir, nThreads, out);
    free(bs.arr);
    return n_out;
}

int taf_summary_main(int argc, char *argv[]) {
    char *inputFile = NULL, *outputFile = NULL, *refGenome = NULL, *tmpDir = NULL;
    int allRefs = 0;

    static struct option lopts[] = {
        {"mergeGap",        required_argument, 0, 'g'},
        {"minSize",         required_argument, 0, 's'},
        {"maxSize",         required_argument, 0, 'S'},
        {"minSeqSize",      required_argument, 0, 'q'},
        {"chainOverlapFrac",required_argument, 0, 'F'},
        {"tmp",             required_argument, 0, 'T'},
        {"threads",         required_argument, 0, 'p'},
        {"allRefs",         no_argument,       0, 'A'},
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
            case 'g': mergeGap   = atoi(optarg); break;
            case 's': minSize    = atoi(optarg); break;
            case 'S': maxSize    = atoi(optarg); break;
            case 'q': minSeqSize = atoll(optarg); break;
            case 'F': chainOverlapFrac = atof(optarg); break;
            case 'T': tmpDir     = optarg; break;
            case 'p': nThreads   = atoi(optarg); break;
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
                  "  --tmp DIR (TMPDIR|/tmp)  scratch dir for the packed per-chrom buckets (needs\n"
                  "            room for ~the packed record volume); bounds peak RAM\n"
                  "  --threads N (ncores)  parallel per-chrom sort+merge workers\n"
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
    if (chainOverlapFrac < 0.0 || chainOverlapFrac > 1.0) {
        fprintf(stderr, "taffy summary: --chainOverlapFrac must be in [0,1]\n");
        return 1;
    }
    if (tmpDir == NULL) { tmpDir = getenv("TMPDIR"); if (tmpDir == NULL) tmpDir = "/tmp"; }
    if (nThreads <= 0) { long nc = sysconf(_SC_NPROCESSORS_ONLN); nThreads = (nc > 0) ? (int) nc : 1; }

    /* one open fd per reference chrom during the scan -- raise the soft limit */
    struct rlimit rl;
    if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur < rl.rlim_max) {
        rl.rlim_cur = rl.rlim_max; setrlimit(RLIMIT_NOFILE, &rl);
    }

    scorer_init();

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
        /* -------- fan-out: ONE scan, every leaf genome is a master routed to its
         * own per-chrom buckets; then finalize each to <outDir>/<genome>.bed.
         * Each per-reference bed equals the single-ref `-r <genome>` run. */
        if (mkdir(outputFile, 0777) != 0 && errno != EEXIST) {
            fprintf(stderr, "taffy summary: cannot create output dir %s\n", outputFile); return 1;
        }
        stHash *refMap = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
        stList *refList = stList_construct();
        AllRefsCtx acx = { gt, srcInt, refMap, refList, tmpDir };
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
        for (int i = 0; i < nRefs; i++) {
            RefState *rs = stList_get(refList, i);
            snprintf(path, sizeof path, "%s/%s.bed", outputFile, rs->name);
            FILE *rout = fopen(path, "w");
            if (rout == NULL) { fprintf(stderr, "taffy summary: cannot write %s\n", path); return 1; }
            int64_t no = rebucket_and_finalize(rs, srcInt, tmpDir, nThreads, rout);
            fclose(rout);
            total_out += no;
            fprintf(stderr, "  %s: %" PRId64 " master rows, %" PRId64 " raw recs -> %" PRId64 " merged rows\n",
                    rs->name, rs->n_with_ref, rs->raw_n, no);
            interner_free(rs->chromInt); free(rs->name); free(rs);
        }
        fprintf(stderr, "taffy summary --allRefs: %d references, %" PRId64 " blocks, "
                "%" PRId64 " total merged rows (single-cover F=%.2f, %d threads)\n",
                nRefs, n_blocks, total_out, chainOverlapFrac, nThreads);
        stList_destruct(refList);
    } else {
        FILE *out = (outputFile == NULL) ? stdout : fopen(outputFile, "w");
        if (out == NULL) { fprintf(stderr, "taffy summary: cannot write %s\n", outputFile); return 1; }

        /* -------- single reference: stream blocks, score, single-cover species
         * in-block, route to per-chrom buckets.  Keep only the previous block
         * alive (TAF deltas reference it); free it on the next read. */
        Interner *chromInt = interner_new();
        BucketSet bs = { NULL, 0, 0 };
        ScoreCtx cx = { gt, refGenome, strlen(refGenome), chromInt, srcInt, &bs, tmpDir, 0, 0, 0 };
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

        int64_t n_out = finalize_reference(&bs, chromInt, srcInt, tmpDir, nThreads, out);
        if (out != stdout) fclose(out);

        int nt = nThreads; if (nt > bs.n) nt = bs.n; if (nt < 1) nt = 1;
        fprintf(stderr, "taffy summary: %" PRId64 " blocks (%" PRId64 " reference master rows, "
                "single-cover F=%.2f), %" PRId64 " raw records -> %" PRId64 " merged summary rows "
                "(%d chroms, %d threads)\n",
                n_blocks, n_with_ref, chainOverlapFrac, raw_n, n_out, bs.n, nt);
        if (n_split_skipped > 0)
            fprintf(stderr, "taffy summary: NOTE: %" PRId64 " reference master row(s) had span > "
                    "maxSize=%d and were NOT split (cosmetic; universal-MAF blocks are normally small)\n",
                    n_split_skipped, maxSize);

        free(bs.arr);
        interner_free(chromInt);
    }

    interner_free(srcInt);
    gerp_tree_destruct(gt);
    LI_destruct(li);
    if (header_tag) tag_destruct(header_tag);
    return 0;
}
