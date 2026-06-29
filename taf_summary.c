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
 * MEMORY: universal block order != reference order, so the per-record entries
 * must be sorted by (src, chrom, start) before the region merge.  Buffering them
 * all in RAM scales with the SUMMARY size (the 28-species fish generates enough
 * raw records to OOM a 31 GB box).  Instead we STREAM raw records to a temp file,
 * EXTERNAL-`sort` it (spills to disk under --tmp), and stream-merge the sorted
 * result -- peak RAM is then a small constant (the merge's one pending block +
 * IO buffers), independent of record count.  Disk use is ~the raw record volume;
 * point --tmp at a roomy filesystem.
 *
 * Output is bed3+4 (chrom start end src score leftStatus rightStatus), the same
 * schema as kent's mafSummary.as, grouped by (src, chrom).  bedToBigBed wants
 * (chrom, start) order, so pipe through `sort -k1,1 -k2,2n` then `bedToBigBed
 * -type=bed3+4 -as=mafSummary.as -tab`.
 *
 *  Released under the MIT license, see LICENSE.txt
 */

#define _GNU_SOURCE   /* mkstemp/popen/fdopen */
#include "taf.h"
#include "gerp.h"
#include "line_iterator.h"
#include "sonLib.h"
#include <ctype.h>
#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
/* record streaming (to temp) + streaming per-(src,chrom) region merge    */
/* ===================================================================== */

static char *intern(stHash *pool, const char *s) {
    char *v = stHash_search(pool, (void *) s);
    if (v == NULL) { v = stString_copy(s); stHash_insert(pool, v, v); }
    return v;
}

/* bed3+4; empty left/right status (kent leaves them blank for these) */
static void emit(FILE *f, const char *chrom, const char *src, int64_t start, int64_t end, double score) {
    fprintf(f, "%s\t%" PRId64 "\t%" PRId64 "\t%s\t%g\t\t\n", chrom, start, end, src, score);
}

/* one raw record to the temp: src \t chrom \t start \t end \t score.
 * src/chrom FIRST so external `sort -k1 -k2 -k3n` groups by (src,chrom) and
 * orders by start, matching the merge's expectations.  score at full %.17g so
 * the merge sees the exact doubles (no rounding before length-weighted avg). */
static void rec_to_tmp(FILE *tmp, const char *src, const char *chrom,
                       int64_t start, int64_t end, double score) {
    fprintf(tmp, "%s\t%s\t%" PRId64 "\t%" PRId64 "\t%.17g\n", src, chrom, start, end, score);
}

/* Stream-merge the EXTERNALLY-sorted temp (records grouped by (src,chrom),
 * ordered by start) using kent mafToBigMafSummary's plain region merge: merge
 * adjacent within mergeGap (length-weighted average), output a block once it
 * grows past minSize, flush at each group boundary.  Constant memory. */
static int64_t stream_merge(FILE *sorted, FILE *out) {
    int64_t n_out = 0;
    char line[4096];
    char src[1024], chrom[1024], psrc[1024] = "\1", pchrom[1024] = "\1";
    int have = 0; int64_t pstart = 0, pend = 0; double pscore = 0;
    while (fgets(line, sizeof line, sorted) != NULL) {
        int64_t start, end; double score;
        if (sscanf(line, "%1023[^\t]\t%1023[^\t]\t%" SCNd64 "\t%" SCNd64 "\t%lf",
                   src, chrom, &start, &end, &score) != 5)
            continue;
        if (strcmp(src, psrc) != 0 || strcmp(chrom, pchrom) != 0) {
            if (have) { emit(out, pchrom, psrc, pstart, pend, pscore); n_out++; have = 0; }
            strcpy(psrc, src); strcpy(pchrom, chrom);
        }
        if (have && (start + 1 - pend < mergeGap)) {
            double total = pscore * (double)(pend - pstart) + score * (double)(end - start);
            pend = end;
            pscore = total / (double)(pend - pstart);
        } else {
            if (have) { emit(out, pchrom, psrc, pstart, pend, pscore); n_out++; }
            pstart = start; pend = end; pscore = score; have = 1;
        }
        if (pend - pstart > minSize) { emit(out, pchrom, psrc, pstart, pend, pscore); n_out++; have = 0; }
    }
    if (have) { emit(out, pchrom, psrc, pstart, pend, pscore); n_out++; }
    return n_out;
}

/* per-block scoring context (streams raw records to cx->tmp) */
typedef struct {
    GerpTree *gt; const char *ref; size_t refLen;
    stHash *pool;
    FILE *tmp; int64_t raw_n;
    int64_t n_with_ref, n_split_skipped;
} ScoreCtx;

/* For EACH reference (master) row in the block, emit one record per leaf
 * species, scored master-vs-species in the reference's coords.  Single-cover
 * each species in-block (best-scoring paralog row) unless --chainOverlapFrac>=1. */
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
        char *chromI = intern(cx->pool, chrom);
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
            char *gi = intern(cx->pool, g);
            if (single_cover) {
                double *cur = stHash_search(best, gi);
                if (cur == NULL) { cur = st_malloc(sizeof(double)); *cur = sc; stHash_insert(best, gi, cur); }
                else if (sc > *cur) *cur = sc;
            } else {
                rec_to_tmp(cx->tmp, gi, chromI, rstart, rend, sc); cx->raw_n++;
            }
        }
        if (single_cover) {
            stHashIterator *hit = stHash_getIterator(best);
            char *gi;
            while ((gi = stHash_getNext(hit)) != NULL) {
                double *cur = stHash_search(best, gi);
                rec_to_tmp(cx->tmp, gi, chromI, rstart, rend, *cur); cx->raw_n++;
            }
            stHash_destructIterator(hit);
            stHash_destruct(best);
        }
    }
}

int taf_summary_main(int argc, char *argv[]) {
    char *inputFile = NULL, *outputFile = NULL, *refGenome = NULL, *tmpDir = NULL;

    static struct option lopts[] = {
        {"mergeGap",        required_argument, 0, 'g'},
        {"minSize",         required_argument, 0, 's'},
        {"maxSize",         required_argument, 0, 'S'},
        {"minSeqSize",      required_argument, 0, 'q'},
        {"chainOverlapFrac",required_argument, 0, 'F'},
        {"tmp",             required_argument, 0, 'T'},
        {"help",            no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "i:r:o:g:s:S:q:F:T:h", lopts, NULL)) != -1) {
        switch (c) {
            case 'i': inputFile  = optarg; break;
            case 'r': refGenome  = optarg; break;
            case 'o': outputFile = optarg; break;
            case 'g': mergeGap   = atoi(optarg); break;
            case 's': minSize    = atoi(optarg); break;
            case 'S': maxSize    = atoi(optarg); break;
            case 'q': minSeqSize = atoll(optarg); break;
            case 'F': chainOverlapFrac = atof(optarg); break;
            case 'T': tmpDir     = optarg; break;
            case 'h':
            default:
                fprintf(stderr,
                  "usage: taffy summary -i <universal .taf/.maf[.gz]> -r <refGenome> [-o out.bed]\n"
                  "  Emit a per-reference bigMafSummary (bed3+4) from a universal alignment.\n"
                  "  Output is grouped by (src,chrom); pipe through `sort -k1,1 -k2,2n` then\n"
                  "  bedToBigBed -type=bed3+4 for the browser.\n"
                  "  -i FILE   universal MAF/TAF (its `# hal` header tree names leaves vs ancestors)\n"
                  "  -r NAME   reference genome (the 'master'); its leaf neighbours are scored\n"
                  "  -o FILE   output bed (default stdout)\n"
                  "  --tmp DIR (TMPDIR|/tmp)  scratch dir for the external sort spill (needs room\n"
                  "            for ~the raw record volume); bounds peak RAM regardless of size\n"
                  "  --chainOverlapFrac F (%.2f)  per-species in-block dedup: 0=single-cover\n"
                  "            (best paralog row, matches reference-anchored summary), 1=keep all\n"
                  "  --mergeGap N (%d)  --minSize N (%d)  --maxSize N (%d)  --minSeqSize N (%" PRId64 ")\n",
                  chainOverlapFrac, mergeGap, minSize, maxSize, minSeqSize);
                return (c == 'h') ? 0 : 1;
        }
    }
    if (inputFile == NULL || refGenome == NULL) {
        fprintf(stderr, "taffy summary: -i and -r are required\n");
        return 1;
    }
    if (chainOverlapFrac < 0.0 || chainOverlapFrac > 1.0) {
        fprintf(stderr, "taffy summary: --chainOverlapFrac must be in [0,1]\n");
        return 1;
    }
    if (tmpDir == NULL) { tmpDir = getenv("TMPDIR"); if (tmpDir == NULL) tmpDir = "/tmp"; }

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
    if (gerp_tree_is_ancestor(gt, refGenome)) {
        fprintf(stderr, "taffy summary: -r '%s' is an internal (ancestor) node, not a genome\n", refGenome);
        return 1;
    }

    FILE *out = (outputFile == NULL) ? stdout : fopen(outputFile, "w");
    if (out == NULL) { fprintf(stderr, "taffy summary: cannot write %s\n", outputFile); return 1; }

    /* temp file for raw records (streamed during the scan, external-sorted after) */
    char tmpl[PATH_MAX];
    snprintf(tmpl, sizeof tmpl, "%s/taffy_summary_XXXXXX", tmpDir);
    int tfd = mkstemp(tmpl);
    if (tfd < 0) { fprintf(stderr, "taffy summary: cannot create temp in %s\n", tmpDir); return 1; }
    FILE *recf = fdopen(tfd, "w");
    if (recf == NULL) { fprintf(stderr, "taffy summary: fdopen temp failed\n"); close(tfd); remove(tmpl); return 1; }

    /* -------- stream blocks: score each, single-covering species in-block --------
     * The universal --noRefDupes invariant (tui.h) makes the reference single-
     * cover by base, so no cross-block chaining is needed.  We keep only the
     * previous block alive (TAF coordinate deltas reference it) and free it once
     * the next block is read.  Raw records go to recf, not RAM. */
    stHash *pool = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    ScoreCtx cx = { gt, refGenome, strlen(refGenome), pool, recf, 0, 0, 0 };
    Alignment *prev = NULL;
    int64_t n_blocks = 0;
    while (1) {
        Alignment *aln = (input_format == 0) ? taf_read_block(prev, rle, li) : maf_read_block(li);
        if (prev != NULL) { alignment_destruct(prev, 1); prev = NULL; }   /* consumed by the read above */
        if (aln == NULL) break;
        n_blocks++;
        score_block(&cx, aln);
        prev = aln;
    }
    int64_t raw_n = cx.raw_n, n_with_ref = cx.n_with_ref, n_split_skipped = cx.n_split_skipped;
    if (fclose(recf) != 0) { fprintf(stderr, "taffy summary: error writing temp %s\n", tmpl); remove(tmpl); return 1; }

    /* external sort by (src, chrom, start), then stream-merge -- constant RAM.
     * -t<TAB> so the fields split on tab; -T puts the sort's own spill on tmpDir. */
    char cmd[2 * PATH_MAX + 128];
    snprintf(cmd, sizeof cmd, "LC_ALL=C sort -S 1G -t'\t' -k1,1 -k2,2 -k3,3n -T '%s' '%s'", tmpDir, tmpl);
    FILE *sorted = popen(cmd, "r");
    if (sorted == NULL) { fprintf(stderr, "taffy summary: cannot spawn sort\n"); remove(tmpl); return 1; }
    int64_t n_out = stream_merge(sorted, out);
    int prc = pclose(sorted);
    remove(tmpl);
    if (prc != 0) {
        fprintf(stderr, "taffy summary: external sort failed (exit %d) -- check --tmp '%s' has space\n", prc, tmpDir);
        if (out != stdout) fclose(out);
        return 1;
    }

    if (out != stdout) fclose(out);
    fprintf(stderr, "taffy summary: %" PRId64 " blocks (%" PRId64 " reference master rows, "
            "single-cover F=%.2f), %" PRId64 " raw records -> %" PRId64 " merged summary rows\n",
            n_blocks, n_with_ref, chainOverlapFrac, raw_n, n_out);
    if (n_split_skipped > 0)
        fprintf(stderr, "taffy summary: NOTE: %" PRId64 " reference master row(s) had span > "
                "maxSize=%d and were NOT split (cosmetic; universal-MAF blocks are normally small)\n",
                n_split_skipped, maxSize);

    stHash_destruct(pool);
    gerp_tree_destruct(gt);
    LI_destruct(li);
    if (header_tag) tag_destruct(header_tag);
    return 0;
}
