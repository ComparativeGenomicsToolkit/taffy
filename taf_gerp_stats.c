/*
 * taffy gerp-rank: post-process a gerp RS wig + depth wig into a
 * depth-corrected + percentile-ranked wig.
 *
 * Three streaming passes (each is one walk of the joint wig pair):
 *   Pass 1: per-depth (count, sum, sum_sq) -> mean(d) + stddev(d).
 *   Pass 2: z-score per column -> histogram of z values.
 *   Pass 3: emit (chrom, pos, percentile) wig (and/or zscore wig).
 *
 * Each pass batches N tuples and dispatches the per-line work to an
 * OpenMP worker pool.  bgzf_mt parallelises the decompress on input.
 *
 * Input wig contract: variableStep, one chrom per variableStep header,
 * data lines of "<pos> <value>".  Depth wig is a positional superset of
 * the RS wig (every non-gap row-0 base in the original universal MAF);
 * RS wig has only the scored positions.  Tool asserts the RS positions
 * are all present in the depth wig in the same chrom order.
 */

#include "taf.h"
#include "gerp_rank.h"
#include "sonLib.h"
#include <ctype.h>
#include <getopt.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////////////////////////////////////////////////////////
// Streaming wig reader.  Reads one wig line at a time.  variableStep
// headers update the current chrom; data lines yield (pos, value).
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    LI    *li;
    char  *cur_chrom;        // owned strdup; replaced on each variableStep
    bool   eof;
    // Pending peeked tuple (set by `_advance` when it's already buffered).
    bool    have_pending;
    int64_t pending_pos;
    double  pending_val;
    char   *pending_chrom;   // borrowed: points into cur_chrom (do not free)
} WigStream;

static WigStream *wig_stream_open(const char *path) {
    FILE *fh = fopen(path, "r");
    if (fh == NULL) {
        fprintf(stderr, "taffy gerp-rank: cannot open wig file: %s\n", path);
        return NULL;
    }
    WigStream *ws = st_calloc(1, sizeof(WigStream));
    ws->li = LI_construct(fh);
    return ws;
}

static void wig_stream_close(WigStream *ws) {
    if (ws == NULL) return;
    if (ws->li != NULL) LI_destruct(ws->li);
    free(ws->cur_chrom);
    free(ws);
}

// Parse one wig line into ws's pending slot, or set ws->eof if we ran out.
// Returns true if a data tuple is now available in `pending_*`; false on
// EOF.  Skips comment / blank / fixedStep lines and updates cur_chrom on
// variableStep.
static bool wig_stream_next(WigStream *ws) {
    if (ws->eof) return false;
    while (1) {
        char *line = LI_get_next_line(ws->li);
        if (line == NULL) { ws->eof = true; ws->have_pending = false; return false; }
        // Skip leading whitespace.
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0' || *p == '\n' || *p == '#') { free(line); continue; }
        // variableStep header.
        if (strncmp(p, "variableStep", 12) == 0) {
            // Look for "chrom=<name>".
            char *c = strstr(p, "chrom=");
            if (c == NULL) {
                fprintf(stderr, "taffy gerp-rank: variableStep line missing chrom: %s\n", p);
                free(line);
                ws->eof = true;
                return false;
            }
            c += 6;
            char *e = c;
            while (*e && *e != ' ' && *e != '\t' && *e != '\n' && *e != '\r') e++;
            *e = '\0';
            free(ws->cur_chrom);
            ws->cur_chrom = stString_copy(c);
            free(line);
            continue;
        }
        if (strncmp(p, "fixedStep", 9) == 0) {
            fprintf(stderr, "taffy gerp-rank: fixedStep wig is not supported -- "
                            "convert with wigToBedGraph or similar first.\n");
            free(line);
            ws->eof = true;
            return false;
        }
        if (strncmp(p, "track", 5) == 0 || strncmp(p, "browser", 7) == 0) {
            // UCSC header line.  Tolerant skip.
            free(line);
            continue;
        }
        // Data line: "<pos> <value>".
        char *endp = NULL;
        long long pos = strtoll(p, &endp, 10);
        if (endp == p) {
            // Garbage line; skip but warn.
            st_logDebug("taffy gerp-rank: skipping unparseable wig line: %.40s\n", p);
            free(line);
            continue;
        }
        while (*endp == ' ' || *endp == '\t') endp++;
        char *valp = endp;
        double val = strtod(valp, &endp);
        if (endp == valp) {
            st_logDebug("taffy gerp-rank: missing value on wig line: %.40s\n", p);
            free(line);
            continue;
        }
        if (ws->cur_chrom == NULL) {
            fprintf(stderr, "taffy gerp-rank: data line before any variableStep header\n");
            free(line);
            ws->eof = true;
            return false;
        }
        ws->have_pending  = true;
        ws->pending_pos   = (int64_t)pos;
        ws->pending_val   = val;
        ws->pending_chrom = ws->cur_chrom;
        free(line);
        return true;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Joint streamer: drives the RS wig; advances the depth wig in lockstep
// to find the same (chrom, pos).  Returns false on EOF of RS wig.
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    WigStream *rs;
    WigStream *depth;       // may be NULL if no depth file was supplied
    // Track when chrom changes so the caller can emit variableStep
    // headers in pass 3.  Set on every successful tuple read; toggled
    // to false after the caller acknowledges via _ack_chrom().
    bool       chrom_changed;
    const char *last_emitted_chrom;
} JointWig;

static JointWig *joint_open(const char *rs_path, const char *depth_path) {
    WigStream *rs = wig_stream_open(rs_path);
    if (rs == NULL) return NULL;
    WigStream *de = NULL;
    if (depth_path != NULL) {
        de = wig_stream_open(depth_path);
        if (de == NULL) { wig_stream_close(rs); return NULL; }
    }
    JointWig *j = st_calloc(1, sizeof(JointWig));
    j->rs    = rs;
    j->depth = de;
    return j;
}

static void joint_close(JointWig *j) {
    if (j == NULL) return;
    wig_stream_close(j->rs);
    wig_stream_close(j->depth);
    free(j);
}

// Advance the depth stream until its current tuple matches (chrom, pos)
// from the RS stream.  Returns true on match.  Returns false (and prints
// an error) if we run off the end of depth OR pass `pos` without hitting
// it (depth was expected to be a superset).
static bool advance_depth_to(JointWig *j, const char *chrom, int64_t pos,
                             int64_t *out_depth) {
    if (j->depth == NULL) { *out_depth = 0; return true; }
    while (1) {
        if (!j->depth->have_pending && !wig_stream_next(j->depth)) {
            fprintf(stderr, "taffy gerp-rank: depth wig EOF before RS position "
                            "%s:%" PRIi64 "\n", chrom, pos);
            return false;
        }
        if (strcmp(j->depth->pending_chrom, chrom) != 0) {
            // Different chrom: discard this depth tuple, keep advancing.
            j->depth->have_pending = false;
            continue;
        }
        if (j->depth->pending_pos < pos) {
            j->depth->have_pending = false;
            continue;
        }
        if (j->depth->pending_pos == pos) {
            *out_depth = (int64_t)j->depth->pending_val;
            j->depth->have_pending = false;
            return true;
        }
        // depth.pos > pos -> the RS position isn't in the depth wig.
        fprintf(stderr, "taffy gerp-rank: depth wig has no entry for RS position "
                        "%s:%" PRIi64 " (next depth: %s:%" PRIi64 "). "
                        "The two wigs must come from the same gerp run.\n",
                chrom, pos,
                j->depth->pending_chrom, j->depth->pending_pos);
        return false;
    }
}

// Pull the next (chrom, pos, rs, depth) tuple.  Returns false on EOF of RS.
// On error inside depth alignment, returns false AND sets *err = 1.
static bool joint_next(JointWig *j, const char **chrom_out, int64_t *pos_out,
                       double *rs_out, int64_t *depth_out, int *err) {
    *err = 0;
    if (!j->rs->have_pending && !wig_stream_next(j->rs)) return false;
    const char *chrom = j->rs->pending_chrom;
    int64_t     pos   = j->rs->pending_pos;
    double      rs    = j->rs->pending_val;
    j->rs->have_pending = false;

    int64_t depth = 0;
    if (!advance_depth_to(j, chrom, pos, &depth)) { *err = 1; return false; }

    *chrom_out = chrom;
    *pos_out   = pos;
    *rs_out    = rs;
    *depth_out = depth;
    if (j->last_emitted_chrom == NULL ||
        strcmp(j->last_emitted_chrom, chrom) != 0) {
        j->chrom_changed       = true;
        j->last_emitted_chrom  = chrom;
    } else {
        j->chrom_changed = false;
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////
// Pass 1: per-depth stats.  Stream the joint wig; tally (count, sum,
// sum_sq) per integer depth.  No depth-correction yet.
/////////////////////////////////////////////////////////////////////////////

static int pass1_collect_stats(JointWig *j, DepthStats *ds_global,
                               int n_threads, int batch_size, int64_t *n_total) {
    // Batch of (depth, rs) tuples.  Workers populate per-thread DepthStats.
    int64_t *b_depth = st_malloc((size_t)batch_size * sizeof(int64_t));
    double  *b_rs    = st_malloc((size_t)batch_size * sizeof(double));
    DepthStats **per_thread = st_malloc((size_t)n_threads * sizeof(DepthStats *));
    for (int t = 0; t < n_threads; t++) per_thread[t] = depth_stats_construct(ds_global->max_depth);

    int err = 0;
    int64_t total = 0;
    while (!err) {
        int n = 0;
        while (n < batch_size) {
            const char *chrom = NULL;
            int64_t pos = 0;
            double  rs  = 0.0;
            int64_t d   = 0;
            int     e   = 0;
            if (!joint_next(j, &chrom, &pos, &rs, &d, &e)) { err = e; break; }
            b_depth[n] = d;
            b_rs[n]    = rs;
            n++;
        }
        if (n == 0) break;
        total += n;
        #pragma omp parallel for schedule(static) num_threads(n_threads)
        for (int i = 0; i < n; i++) {
#ifdef _OPENMP
            int t = omp_get_thread_num();
#else
            int t = 0;
#endif
            depth_stats_observe(per_thread[t], b_depth[i], b_rs[i]);
        }
    }
    for (int t = 0; t < n_threads; t++) {
        depth_stats_merge(ds_global, per_thread[t]);
        depth_stats_destruct(per_thread[t]);
    }
    free(per_thread);
    free(b_depth);
    free(b_rs);
    *n_total = total;
    return err;
}

/////////////////////////////////////////////////////////////////////////////
// Pass 2: histogram of (depth-corrected) z-scores.
/////////////////////////////////////////////////////////////////////////////

static int pass2_collect_histogram(JointWig *j, const DepthStats *ds,
                                   Histogram *h_global, int n_threads,
                                   int batch_size) {
    int64_t *b_depth = st_malloc((size_t)batch_size * sizeof(int64_t));
    double  *b_rs    = st_malloc((size_t)batch_size * sizeof(double));
    Histogram **per_thread = st_malloc((size_t)n_threads * sizeof(Histogram *));
    for (int t = 0; t < n_threads; t++)
        per_thread[t] = histogram_construct(h_global->n_bins, h_global->min_val, h_global->max_val);

    int err = 0;
    while (!err) {
        int n = 0;
        while (n < batch_size) {
            const char *chrom = NULL;
            int64_t pos = 0;
            double  rs  = 0.0;
            int64_t d   = 0;
            int     e   = 0;
            if (!joint_next(j, &chrom, &pos, &rs, &d, &e)) { err = e; break; }
            b_depth[n] = d;
            b_rs[n]    = rs;
            n++;
        }
        if (n == 0) break;
        #pragma omp parallel for schedule(static) num_threads(n_threads)
        for (int i = 0; i < n; i++) {
#ifdef _OPENMP
            int t = omp_get_thread_num();
#else
            int t = 0;
#endif
            double z = depth_stats_zscore(ds, b_depth[i], b_rs[i]);
            histogram_observe(per_thread[t], z);
        }
    }
    for (int t = 0; t < n_threads; t++) {
        histogram_merge(h_global, per_thread[t]);
        histogram_destruct(per_thread[t]);
    }
    free(per_thread);
    free(b_depth);
    free(b_rs);
    return err;
}

/////////////////////////////////////////////////////////////////////////////
// Pass 3: emit output wig with the chosen (z, percentile, or both).
/////////////////////////////////////////////////////////////////////////////

typedef enum { EMIT_PERCENTILE, EMIT_ZSCORE, EMIT_BOTH } EmitMode;

// Reuse the per-block GerpBuf-style buffers from taf_gerp.c via a copy.
typedef struct {
    char  *buf;
    size_t len;
    size_t cap;
} OutBuf;

static void outbuf_init(OutBuf *b, size_t cap)  { b->cap = cap; b->buf = st_malloc(cap); b->len = 0; }
static void outbuf_destroy(OutBuf *b)           { free(b->buf); b->buf = NULL; b->len = b->cap = 0; }
static void outbuf_reset(OutBuf *b)             { b->len = 0; }
static inline void outbuf_reserve(OutBuf *b, size_t n) {
    if (b->len + n <= b->cap) return;
    while (b->len + n > b->cap) b->cap *= 2;
    b->buf = st_realloc(b->buf, b->cap);
}
static inline void outbuf_putn(OutBuf *b, const char *s, size_t n) {
    outbuf_reserve(b, n);
    memcpy(b->buf + b->len, s, n);
    b->len += n;
}

static int pass3_emit(JointWig *j, const DepthStats *ds, const Histogram *h,
                      LW *out, EmitMode emit_mode, int n_threads, int batch_size) {
    int64_t *b_depth = st_malloc((size_t)batch_size * sizeof(int64_t));
    double  *b_rs    = st_malloc((size_t)batch_size * sizeof(double));
    int64_t *b_pos   = st_malloc((size_t)batch_size * sizeof(int64_t));
    char   **b_chrom = st_malloc((size_t)batch_size * sizeof(char *));
    bool    *b_first_of_chrom = st_calloc((size_t)batch_size, sizeof(bool));
    OutBuf  *b_out   = st_calloc((size_t)batch_size, sizeof(OutBuf));
    for (int i = 0; i < batch_size; i++) outbuf_init(&b_out[i], 256);

    int err = 0;
    while (!err) {
        int n = 0;
        while (n < batch_size) {
            const char *chrom = NULL;
            int64_t pos = 0;
            double  rs  = 0.0;
            int64_t d   = 0;
            int     e   = 0;
            bool was_changed_prev = j->chrom_changed;  // unused; we re-check below
            (void) was_changed_prev;
            if (!joint_next(j, &chrom, &pos, &rs, &d, &e)) { err = e; break; }
            b_chrom[n] = (char *) chrom;  // borrowed; cur_chrom in WigStream
            b_pos[n]   = pos;
            b_rs[n]    = rs;
            b_depth[n] = d;
            b_first_of_chrom[n] = j->chrom_changed;
            n++;
        }
        if (n == 0) break;

        #pragma omp parallel for schedule(static) num_threads(n_threads)
        for (int i = 0; i < n; i++) {
            outbuf_reset(&b_out[i]);
            if (b_first_of_chrom[i]) {
                outbuf_putn(&b_out[i], "variableStep chrom=", 19);
                outbuf_putn(&b_out[i], b_chrom[i], strlen(b_chrom[i]));
                outbuf_putn(&b_out[i], "\n", 1);
            }
            double z = depth_stats_zscore(ds, b_depth[i], b_rs[i]);
            double p = histogram_percentile(h, z);
            char line[96];
            int ln = 0;
            switch (emit_mode) {
                case EMIT_PERCENTILE:
                    ln = snprintf(line, sizeof(line), "%" PRIi64 " %.2f\n", b_pos[i], p);
                    break;
                case EMIT_ZSCORE:
                    ln = snprintf(line, sizeof(line), "%" PRIi64 " %.4f\n", b_pos[i], z);
                    break;
                case EMIT_BOTH:
                    ln = snprintf(line, sizeof(line), "%" PRIi64 "\t%.4f\t%.2f\n", b_pos[i], z, p);
                    break;
            }
            if (ln > 0) outbuf_putn(&b_out[i], line, (size_t)ln);
        }

        for (int i = 0; i < n; i++) {
            if (b_out[i].len > 0) LW_putn(out, b_out[i].buf, b_out[i].len);
        }
    }
    for (int i = 0; i < batch_size; i++) outbuf_destroy(&b_out[i]);
    free(b_out);
    free(b_first_of_chrom);
    free(b_chrom);
    free(b_pos);
    free(b_rs);
    free(b_depth);
    return err;
}

/////////////////////////////////////////////////////////////////////////////
// Driver
/////////////////////////////////////////////////////////////////////////////

static void usage(void) {
    fprintf(stderr, "taffy gerp-rank [options]\n");
    fprintf(stderr, "Depth-correct and percentile-rank the RS wig from `taffy gerp`.\n");
    fprintf(stderr, "-i --inputFile : RS wig from `taffy gerp -o` (required)\n");
    fprintf(stderr, "-D --depthFile : Depth wig from `taffy gerp -D` (required for --depthCorrect=zscore)\n");
    fprintf(stderr, "-o --outputFile : Output wig (writes stdout if omitted)\n");
    fprintf(stderr, "-c --useCompression : Bgzip the output\n");
    fprintf(stderr, "--depthCorrect=MODE : none | sqrt | zscore (default zscore)\n");
    fprintf(stderr, "--emit=MODE : percentile | zscore | both (default percentile)\n");
    fprintf(stderr, "--histBins N : Histogram bin count for percentile lookup (default 1000)\n");
    fprintf(stderr, "--histRange MIN MAX : Score range covered by histogram (default -10 10)\n");
    fprintf(stderr, "--minBucketN N : Min observations per depth bucket before falling back (default 30)\n");
    fprintf(stderr, "--maxDepth N : Largest expected depth (default 4096; auto-detected if too small)\n");
    fprintf(stderr, "-T --threads N : Worker threads + bgzf I/O (default 1)\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help\n");
}

enum {
    OPT_DEPTH_CORRECT = 256,
    OPT_EMIT,
    OPT_HIST_BINS,
    OPT_HIST_RANGE,
    OPT_MIN_BUCKET_N,
    OPT_MAX_DEPTH,
};

typedef enum { DC_NONE, DC_SQRT, DC_ZSCORE } DepthCorrect;

int taf_gerp_rank_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    char *logLevelString = NULL;
    char *inputFile  = NULL;
    char *depthFile  = NULL;
    char *outputFile = NULL;
    bool use_compression = false;
    DepthCorrect depth_correct = DC_ZSCORE;
    EmitMode emit_mode = EMIT_PERCENTILE;
    int64_t hist_bins = 1000;
    double  hist_min = -10.0, hist_max = 10.0;
    int64_t min_bucket_n = 30;
    int64_t max_depth = 4096;
    int n_threads = 1;
    int batch_size = 4096;

    while (1) {
        static struct option long_options[] = {
            { "logLevel",       required_argument, 0, 'l' },
            { "inputFile",      required_argument, 0, 'i' },
            { "depthFile",      required_argument, 0, 'D' },
            { "outputFile",     required_argument, 0, 'o' },
            { "useCompression", no_argument,       0, 'c' },
            { "depthCorrect",   required_argument, 0, OPT_DEPTH_CORRECT },
            { "emit",           required_argument, 0, OPT_EMIT },
            { "histBins",       required_argument, 0, OPT_HIST_BINS },
            { "histRange",      required_argument, 0, OPT_HIST_RANGE },
            { "minBucketN",     required_argument, 0, OPT_MIN_BUCKET_N },
            { "maxDepth",       required_argument, 0, OPT_MAX_DEPTH },
            { "threads",        required_argument, 0, 'T' },
            { "help",           no_argument,       0, 'h' },
            { 0, 0, 0, 0 }
        };
        int idx = 0;
        int64_t key = getopt_long(argc, argv, "l:i:D:o:cT:h", long_options, &idx);
        if (key == -1) break;
        switch (key) {
            case 'l': logLevelString = optarg; break;
            case 'i': inputFile  = optarg; break;
            case 'D': depthFile  = optarg; break;
            case 'o': outputFile = optarg; break;
            case 'c': use_compression = true; break;
            case 'T': n_threads  = atoi(optarg); break;
            case OPT_DEPTH_CORRECT:
                if      (strcmp(optarg, "none")    == 0) depth_correct = DC_NONE;
                else if (strcmp(optarg, "sqrt")    == 0) depth_correct = DC_SQRT;
                else if (strcmp(optarg, "zscore")  == 0) depth_correct = DC_ZSCORE;
                else { fprintf(stderr, "taffy gerp-rank: unknown --depthCorrect %s\n", optarg); return 1; }
                break;
            case OPT_EMIT:
                if      (strcmp(optarg, "percentile") == 0) emit_mode = EMIT_PERCENTILE;
                else if (strcmp(optarg, "zscore")     == 0) emit_mode = EMIT_ZSCORE;
                else if (strcmp(optarg, "both")       == 0) emit_mode = EMIT_BOTH;
                else { fprintf(stderr, "taffy gerp-rank: unknown --emit %s\n", optarg); return 1; }
                break;
            case OPT_HIST_BINS:    hist_bins    = atol(optarg); break;
            case OPT_HIST_RANGE:
                // Format "<min>:<max>" -- two getopts would be cleaner but
                // long_options doesn't easily support nargs=2.
                {
                    char *colon = strchr(optarg, ':');
                    if (colon == NULL) {
                        fprintf(stderr, "taffy gerp-rank: --histRange requires MIN:MAX\n");
                        return 1;
                    }
                    *colon = '\0';
                    hist_min = atof(optarg);
                    hist_max = atof(colon + 1);
                }
                break;
            case OPT_MIN_BUCKET_N: min_bucket_n = atol(optarg); break;
            case OPT_MAX_DEPTH:    max_depth    = atol(optarg); break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }
    if (n_threads < 1) n_threads = 1;
    if (hist_bins < 2 || hist_max <= hist_min) {
        fprintf(stderr, "taffy gerp-rank: invalid histogram config\n");
        return 1;
    }
    if (inputFile == NULL) {
        fprintf(stderr, "taffy gerp-rank: -i required\n");
        usage();
        return 1;
    }
    if (depth_correct == DC_ZSCORE && depthFile == NULL) {
        fprintf(stderr, "taffy gerp-rank: --depthCorrect=zscore needs -D <depth.wig>\n");
        return 1;
    }
    if (depth_correct == DC_SQRT && depthFile == NULL) {
        fprintf(stderr, "taffy gerp-rank: --depthCorrect=sqrt needs -D <depth.wig>\n");
        return 1;
    }

    st_setLogLevelFromString(logLevelString);
    LI_set_bgzf_threads(n_threads);

    // --- Pass 1: collect per-depth stats. -----------------------------
    JointWig *j1 = joint_open(inputFile, (depth_correct != DC_NONE) ? depthFile : NULL);
    if (j1 == NULL) return 1;
    DepthStats *ds = depth_stats_construct(max_depth);
    int64_t n_total = 0;
    st_logInfo("taffy gerp-rank: pass 1 / 3 (per-depth stats)\n");
    int rc = pass1_collect_stats(j1, ds, n_threads, batch_size, &n_total);
    joint_close(j1);
    if (rc != 0) return 1;
    depth_stats_finalize(ds, min_bucket_n);
    st_logInfo("taffy gerp-rank: pass 1 done -- %" PRIi64 " columns in %" PRIi64 "s\n",
               n_total, (int64_t)(time(NULL) - startTime));

    if (depth_correct == DC_NONE) {
        // With no depth correction, the "z-score" is the raw RS.  Wire
        // the DepthStats so that stddev=1, mean=0 across all buckets:
        for (int64_t d = 0; d <= ds->max_depth; d++) {
            ds->mean[d]        = 0.0;
            ds->stddev[d]      = 1.0;
            ds->fallback_to[d] = d;
        }
    } else if (depth_correct == DC_SQRT) {
        // sqrt: scale RS by sqrt(depth).  Reuse the zscore plumbing by
        // setting mean=0 and stddev = 1/sqrt(d) (then z = RS*sqrt(d)).
        for (int64_t d = 0; d <= ds->max_depth; d++) {
            ds->mean[d]        = 0.0;
            ds->stddev[d]      = (d > 0) ? 1.0 / sqrt((double) d) : 1.0;
            ds->fallback_to[d] = d;
        }
    }
    // DC_ZSCORE: leave finalize()'s output as-is.

    // --- Pass 2: histogram of corrected values. ------------------------
    JointWig *j2 = joint_open(inputFile, (depth_correct != DC_NONE) ? depthFile : NULL);
    if (j2 == NULL) { depth_stats_destruct(ds); return 1; }
    Histogram *h = histogram_construct(hist_bins, hist_min, hist_max);
    st_logInfo("taffy gerp-rank: pass 2 / 3 (z-score histogram, %" PRIi64
               " bins in [%g, %g))\n", hist_bins, hist_min, hist_max);
    rc = pass2_collect_histogram(j2, ds, h, n_threads, batch_size);
    joint_close(j2);
    if (rc != 0) { depth_stats_destruct(ds); histogram_destruct(h); return 1; }
    histogram_finalize(h);
    if (h->n_clipped_lo > 0 || h->n_clipped_hi > 0) {
        st_logInfo("taffy gerp-rank: %" PRIi64 " values clipped below %g, "
                   "%" PRIi64 " clipped at/above %g (of %" PRIi64 " total)\n",
                   h->n_clipped_lo, h->min_val, h->n_clipped_hi, h->max_val, h->n_total);
    }
    st_logInfo("taffy gerp-rank: pass 2 done in %" PRIi64 "s total\n",
               (int64_t)(time(NULL) - startTime));

    // --- Pass 3: emit. -------------------------------------------------
    FILE *out_fh = (outputFile == NULL) ? stdout : fopen(outputFile, "w");
    if (out_fh == NULL) {
        fprintf(stderr, "taffy gerp-rank: cannot open output: %s\n", outputFile);
        depth_stats_destruct(ds); histogram_destruct(h);
        return 1;
    }
    LW *out = LW_construct(out_fh, use_compression);
    JointWig *j3 = joint_open(inputFile, (depth_correct != DC_NONE) ? depthFile : NULL);
    if (j3 == NULL) { LW_destruct(out, false); depth_stats_destruct(ds); histogram_destruct(h); return 1; }
    st_logInfo("taffy gerp-rank: pass 3 / 3 (emit)\n");
    rc = pass3_emit(j3, ds, h, out, emit_mode, n_threads, batch_size);
    joint_close(j3);
    LW_destruct(out, false);
    if (outputFile != NULL) fclose(out_fh);

    st_logInfo("taffy gerp-rank: done in %" PRIi64 "s\n",
               (int64_t)(time(NULL) - startTime));

    depth_stats_destruct(ds);
    histogram_destruct(h);
    return rc;
}
