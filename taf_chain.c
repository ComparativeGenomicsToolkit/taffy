/*
 * taffy tui-chain: build a "chained" .tui from an existing .tui by
 * merging consecutive gap-free (or near-gap-free) runs into longer
 * chains within each chunk.  Output preserves the .tui format (same
 * schema, same reader, same lookup path) but with N x fewer runs per
 * chunk -- the lift hot path (tui_genome_lift_visit_runs) iterates
 * runs, so reducing the count gives a proportional wall-clock win.
 *
 * Sidecar-style: lives alongside the source .tui as <input>.chained.tui,
 * does NOT replace it.  The base .tui remains the precise, run-accurate
 * truth for tools that need per-run fidelity (taffy view, future per-
 * base queries).  taffy lift / blockViz can prefer the chained sidecar
 * for speed when present.
 *
 * Chain semantics:
 *   Two runs A and B (within the SAME genome, SAME seq, SAME strand)
 *   are chainable iff
 *     g_gap = B.g_start - (A.g_start + A.length)   in  [0, maxGap]
 *     t_gap depends on strand:
 *       + strand: B.t_start - (A.t_start + A.length)         in [0, maxGap]
 *       - strand: A.t_start - (B.t_start + B.length)         in [0, maxGap]
 *   The chained run's (g_start, length) covers the SPAN of both inputs.
 *   For + strand, t_start = A.t_start.  For - strand, t_start = B.t_start
 *   (the smaller value, since the target axis decreases with g).
 *
 * The TRADE-OFF: internal gaps inside a chain are not represented in
 * the encoded length; queries at positions within an internal gap will
 * report a target position that doesn't reflect a real source alignment.
 * --maxGap is the user's lever -- 0 = lossless (only chain truly
 * contiguous runs), 1000 = axtChain-equivalent (close-species default).
 *
 *  Released under the MIT license, see LICENSE.txt
 */

#define _POSIX_C_SOURCE 200809L

#include "taf.h"
#include "tui.h"
#include "sonLib.h"

#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <malloc.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

static void usage(void) {
    fprintf(stderr,
        "taffy tui-chain [options] -i <input.maf|input.tui>\n"
        "\n"
        "Build a chained .tui from an existing .tui by merging consecutive\n"
        "gap-free runs into longer chains within each chunk.  Output keeps\n"
        "the same .tui schema, so existing readers (taffy view, blockViz,\n"
        "taffy lift) work unchanged.  The chained file is a sidecar -- the\n"
        "original .tui is not modified.\n"
        "\n"
        "Options:\n"
        "  -i --inputFile FILE     Base MAF/TAF whose <input>.tui will be\n"
        "                          chained, OR with --fromTui, an existing\n"
        "                          .tui directly.  REQUIRED.\n"
        "  -o --outputFile FILE    Output chained .tui path.\n"
        "                          DEFAULT: <input>.chained.tui\n"
        "  -G --maxGap INT         Max gap (in bp) allowed inside a chain on\n"
        "                          BOTH the universal-column axis and the\n"
        "                          target genome axis.  0 = lossless (only\n"
        "                          merge contiguous runs); higher values\n"
        "                          give more aggressive compression at the\n"
        "                          cost of internal precision.  Typical\n"
        "                          values: 100 (conservative), 1000 (axt-\n"
        "                          Chain-equivalent), 10000 (very aggressive).\n"
        "                          DEFAULT 1000.\n"
        "     --fromTui            Treat -i as the source .tui itself.\n"
        "  -l --logLevel STR       critical|info|debug.  DEFAULT critical.\n"
        "  -h --help               Print this help.\n");
}

/* ------------------------------------------------------------------ */
/* In-memory chained run.  Identical shape to TuiRun but seq is an    */
/* owned strdup (the TuiRun.seq pointer borrowed from gl is reset when */
/* gl is destructed at end of genome).                                 */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *seq;     /* borrowed from gl during streaming; later owned */
    int64_t     g_start;
    int64_t     length;
    int64_t     t_start;
    int         strand;
} ChainRun;

typedef struct {
    int64_t     max_gap;
    stHash     *active;        /* const char* seq -> ChainRun* (active chain) */
    stList     *out;           /* list of ChainRun*, in arrival (g) order */
    int64_t     n_in;          /* input runs streamed */
    int64_t     n_out;         /* chained runs produced (after flush) */
} ChainCtx;

/* Visitor: extend the active chain for this seq, or emit-and-restart. */
static void chain_visit_cb(const TuiRun *r, void *user) {
    ChainCtx *cx = (ChainCtx *)user;
    cx->n_in++;
    ChainRun *cur = (ChainRun *)stHash_search(cx->active, (void *)r->seq);
    if (cur == NULL) {
        cur = (ChainRun *)st_malloc(sizeof(ChainRun));
        cur->seq     = r->seq;
        cur->g_start = r->g_start;
        cur->length  = r->length;
        cur->t_start = r->t_start;
        cur->strand  = r->strand;
        stHash_insert(cx->active, (void *)r->seq, cur);
        return;
    }
    /* Can we extend cur with r? */
    if (cur->strand == r->strand) {
        int64_t g_gap = r->g_start - (cur->g_start + cur->length);
        if (g_gap >= 0 && g_gap <= cx->max_gap) {
            int64_t t_gap;
            int64_t new_t_start = cur->t_start;
            if (cur->strand) {
                /* + strand: target advances with g; cur.t_start is the
                 * smallest t (left edge); new chain extends to the right. */
                t_gap = r->t_start - (cur->t_start + cur->length);
            } else {
                /* - strand: target decreases as g advances.  cur was emitted
                 * at the larger g (= smaller t); r is at even larger g (=
                 * even smaller t).  So new chain's t_start = r->t_start. */
                t_gap = cur->t_start - (r->t_start + r->length);
                new_t_start = r->t_start;
            }
            if (t_gap >= 0 && t_gap <= cx->max_gap) {
                cur->t_start = new_t_start;
                cur->length  = (r->g_start + r->length) - cur->g_start;
                return;
            }
        }
    }
    /* Cannot chain: emit cur, start fresh.  Note we append cur (not a
     * copy) and replace the hash slot with a new allocation. */
    stList_append(cx->out, cur);
    cur = (ChainRun *)st_malloc(sizeof(ChainRun));
    cur->seq     = r->seq;
    cur->g_start = r->g_start;
    cur->length  = r->length;
    cur->t_start = r->t_start;
    cur->strand  = r->strand;
    stHash_insert(cx->active, (void *)r->seq, cur);
}

/* After streaming a genome, move any active chains from `active` into
 * `out` (in seq-iteration order; we resort below). */
static void chain_flush_active(ChainCtx *cx) {
    stHashIterator *it = stHash_getIterator(cx->active);
    void *k;
    while ((k = stHash_getNext(it)) != NULL) {
        ChainRun *cur = (ChainRun *)stHash_search(cx->active, k);
        if (cur) stList_append(cx->out, cur);
    }
    stHash_destructIterator(it);
    stHash_destruct(cx->active);
    cx->active = stHash_construct();
    cx->n_out = stList_length(cx->out);
}

/* Sort comparator: (seq string, g_start).  Use strcmp -- pointer-only
 * compare was a footgun (any code that duplicates seq strings, even
 * once, splits a single logical seq into many "different" groups). */
static int chain_cmp(const void *a, const void *b) {
    const ChainRun *p = *(const ChainRun *const *)a;
    const ChainRun *q = *(const ChainRun *const *)b;
    int c = strcmp(p->seq, q->seq);
    if (c != 0) return c;
    if (p->g_start < q->g_start) return -1;
    if (p->g_start > q->g_start) return 1;
    return 0;
}

typedef struct {
    char       *full_name;
    int64_t     seq_idx_in_dir;
    ChainRun  **runs;             /* owned ChainRun pointers (move ownership) */
    int64_t     n_runs;
} SeqGroup;

/* SeqGroup compare: by seq_idx_in_dir asc (matches d-record order). */
static int seq_group_cmp(const void *a, const void *b) {
    const SeqGroup *p = (const SeqGroup *)a;
    const SeqGroup *q = (const SeqGroup *)b;
    if (p->seq_idx_in_dir < q->seq_idx_in_dir) return -1;
    if (p->seq_idx_in_dir > q->seq_idx_in_dir) return 1;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Per-genome emit: write S/C/R chunks honoring TUI caps.              */
/* Via the shared tui_write_sequence.  Each S-record is one seq's      */
/* chained runs, packed into 1..N chunks (TUI_CHUNK_RUNS / G_MAX caps). */
/* ------------------------------------------------------------------ */

/* Bucket chains by seq (with genome prefix re-prepended), produce
 * one SeqGroup per active seq, sorted by seq_idx_in_dir for stable
 * output ordering.  Caller frees groups[] + their inner arrays. */
static int64_t bucket_by_seq(stList *chains, const char *genome,
                             stHash *seq_to_dir_idx,
                             SeqGroup **groups_out) {
    /* Sort chains by (seq, g_start) so seq groups are contiguous. */
    int64_t n = stList_length(chains);
    if (n == 0) { *groups_out = NULL; return 0; }
    ChainRun **arr = (ChainRun **)st_malloc((size_t)n * sizeof(ChainRun *));
    for (int64_t i = 0; i < n; i++) arr[i] = (ChainRun *)stList_get(chains, i);
    qsort(arr, (size_t)n, sizeof(ChainRun *), chain_cmp);

    /* Pass 1: count distinct seqs. */
    int64_t n_g = 0;
    for (int64_t i = 0; i < n; ) {
        int64_t j = i + 1;
        while (j < n && arr[j]->seq == arr[i]->seq) j++;
        n_g++; i = j;
    }
    SeqGroup *groups = (SeqGroup *)st_malloc((size_t)(n_g ? n_g : 1) * sizeof(SeqGroup));
    int64_t gi = 0;
    size_t glen = strlen(genome);
    for (int64_t i = 0; i < n; ) {
        int64_t j = i + 1;
        while (j < n && arr[j]->seq == arr[i]->seq) j++;
        size_t slen = strlen(arr[i]->seq);
        char *full = (char *)st_malloc(glen + 1 + slen + 1);
        memcpy(full, genome, glen);
        full[glen] = '.';
        memcpy(full + glen + 1, arr[i]->seq, slen + 1);
        int64_t *idx_p = (int64_t *)stHash_search(seq_to_dir_idx, full);
        if (idx_p == NULL) {
            st_errAbort("tui-chain: seq %s not in directory map (this is a bug)",
                        full);
        }
        groups[gi].full_name      = full;
        groups[gi].seq_idx_in_dir = *idx_p;
        groups[gi].n_runs         = j - i;
        groups[gi].runs           = (ChainRun **)st_malloc((size_t)(j - i) * sizeof(ChainRun *));
        for (int64_t k = i; k < j; k++) groups[gi].runs[k - i] = arr[k];
        gi++;
        i = j;
    }
    free(arr);
    qsort(groups, (size_t)n_g, sizeof(SeqGroup), seq_group_cmp);
    *groups_out = groups;
    return n_g;
}

/* Encode and write one genome's chained runs via the shared tui_write_sequence
 * (the single source of truth for the on-disk S/C/R layout). */
static void write_genome(OneFile *of, SeqGroup *groups, int64_t n_groups,
                         int64_t *c_ord_emit, int64_t *s_ord_counter,
                         stHash *capture) {
    for (int64_t gi = 0; gi < n_groups; gi++) {
        SeqGroup *g = &groups[gi];
        /* Approximate seq_len = max(t_start + length) across this seq's chains. */
        int64_t max_t_end = 0;
        for (int64_t i = 0; i < g->n_runs; i++) {
            ChainRun *cr = g->runs[i];
            int64_t te = cr->t_start + cr->length;
            if (te > max_t_end) max_t_end = te;
        }
        /* s_ord bookkeeping for the directory's S-ordinal (implicit in write
         * order; the C records carry no per-chunk parent-S-ord). */
        ++(*s_ord_counter);
        int64_t s_ord = *s_ord_counter;
        if (capture != NULL) {
            int64_t *p = (int64_t *)st_malloc(sizeof(int64_t));
            *p = s_ord;
            stHash_insert(capture, stString_copy(g->full_name), p);
        }

        /* Chunk the chained runs (honoring TUI_CHUNK_RUNS + TUI_CHUNK_G_MAX)
         * into the shared write descriptor, then emit via tui_write_sequence. */
        TuiWriteChunk *wc = NULL;
        int64_t nwc = 0, wc_cap = 0;
        int64_t i = 0;
        while (i < g->n_runs) {
            int64_t start = i;
            int64_t cm = 0;
            int64_t g_min = g->runs[start]->g_start;
            int64_t g_max = g_min + g->runs[start]->length;
            int64_t t_min = g->runs[start]->t_start;
            int64_t t_max = t_min + g->runs[start]->length;
            while (i < g->n_runs && cm < TUI_CHUNK_RUNS) {
                ChainRun *cr = g->runs[i];
                int64_t new_g_max = (cr->g_start + cr->length > g_max)
                                      ? cr->g_start + cr->length : g_max;
                /* Honour the g-span cap, BUT always accept at least one
                 * run per chunk -- a single chain may itself exceed
                 * TUI_CHUNK_G_MAX (gap-merged across long regions), and
                 * rejecting it would infinite-loop in the outer while. */
                if (cm > 0 && new_g_max - g_min > TUI_CHUNK_G_MAX) break;
                int64_t te = cr->t_start + cr->length;
                if (cr->t_start < t_min) t_min = cr->t_start;
                if (te > t_max) t_max = te;
                g_max = new_g_max;
                cm++; i++;
            }
            /* Encode triples: (t_start, g_start, lenc) for tui_encode_runs. */
            int64_t *buf = (int64_t *)st_malloc((size_t)cm * 3 * sizeof(int64_t));
            for (int64_t k = 0; k < cm; k++) {
                ChainRun *cr = g->runs[start + k];
                buf[3*k + 0] = cr->t_start;
                buf[3*k + 1] = cr->g_start;
                buf[3*k + 2] = (cr->length << 1) | (cr->strand ? 0 : 1);
            }
            int64_t raw_len = 0, def_len = 0;
            uint8_t *def = tui_encode_runs(buf, cm, &raw_len, &def_len);
            free(buf);

            if (nwc == wc_cap) {
                wc_cap = wc_cap ? wc_cap * 2 : 8;
                wc = (TuiWriteChunk *)st_realloc(wc, wc_cap * sizeof(TuiWriteChunk));
            }
            wc[nwc].g_min = g_min; wc[nwc].g_max = g_max;
            wc[nwc].t_min = t_min; wc[nwc].t_max = t_max;
            wc[nwc].raw_len = raw_len; wc[nwc].def = def; wc[nwc].def_len = def_len;
            nwc++;
        }
        tui_write_sequence(of, g->full_name, max_t_end, wc, nwc, c_ord_emit);
        for (int64_t k = 0; k < nwc; k++) free(wc[k].def);
        free(wc);
        /* Each ChainRun was malloc'd in chain_visit_cb and owned via
         * SeqGroup.runs; free them now that they've been encoded. */
        for (int64_t i = 0; i < g->n_runs; i++) free(g->runs[i]);
        free(g->runs);
        free(g->full_name);
    }
    free(groups);
}

/* ------------------------------------------------------------------ */
/* Directory.  Name-sorted; written after all S/C/R.                     */
/* ------------------------------------------------------------------ */

typedef struct {
    char    *name;
    int64_t  length;
    int64_t  dir_idx;
} DirEntry;

static int dir_entry_cmp(const void *a, const void *b) {
    return strcmp(((const DirEntry *)a)->name, ((const DirEntry *)b)->name);
}

/* ------------------------------------------------------------------ */
/* Main                                                                 */
/* ------------------------------------------------------------------ */

int taf_chain_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    char    *input_file     = NULL;
    char    *output_file    = NULL;
    int64_t  max_gap        = 1000;
    int      from_tui       = 0;
    char    *log_level      = NULL;

    enum { OPT_FROM_TUI = 256 };
    static struct option long_options[] = {
        { "inputFile",   required_argument, 0, 'i' },
        { "outputFile",  required_argument, 0, 'o' },
        { "maxGap",      required_argument, 0, 'G' },
        { "fromTui",     no_argument,       0, OPT_FROM_TUI },
        { "logLevel",    required_argument, 0, 'l' },
        { "help",        no_argument,       0, 'h' },
        { 0, 0, 0, 0 }
    };
    while (1) {
        int oi = 0;
        int k = getopt_long(argc, argv, "i:o:G:l:h", long_options, &oi);
        if (k == -1) break;
        switch (k) {
            case 'i': input_file  = optarg; break;
            case 'o': output_file = optarg; break;
            case 'G': max_gap     = atoll(optarg); break;
            case OPT_FROM_TUI:    from_tui = 1; break;
            case 'l': log_level   = optarg; break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }
    if (input_file == NULL) {
        fprintf(stderr, "ERROR: -i required\n");
        usage(); return 1;
    }
    if (max_gap < 0) {
        fprintf(stderr, "ERROR: --maxGap must be >= 0\n");
        return 1;
    }
    st_setLogLevelFromString(log_level);

    /* Resolve paths. */
    char *src_tui_path;
    if (from_tui) {
        src_tui_path = stString_copy(input_file);
    } else {
        src_tui_path = tui_path(input_file);
    }
    char *out_path_owned = NULL;
    const char *out_path = output_file;
    if (out_path == NULL) {
        /* DEFAULT: <input>.chained.tui (or strip trailing .tui first). */
        size_t slen = strlen(src_tui_path);
        const char *suffix = ".chained.tui";
        if (slen >= 4 && strcmp(src_tui_path + slen - 4, ".tui") == 0) {
            out_path_owned = (char *)st_malloc(slen + strlen(suffix) + 1);
            memcpy(out_path_owned, src_tui_path, slen - 4);
            strcpy(out_path_owned + slen - 4, suffix);
        } else {
            out_path_owned = (char *)st_malloc(slen + strlen(suffix) + 1);
            strcpy(out_path_owned, src_tui_path);
            strcat(out_path_owned, suffix);
        }
        out_path = out_path_owned;
    }

    /* Load source .tui. */
    Tui *src = tui_load(src_tui_path);
    if (src == NULL) {
        fprintf(stderr, "tui-chain: tui_load failed for %s\n", src_tui_path);
        free(src_tui_path); free(out_path_owned);
        return 1;
    }
    int64_t T_src = tui_total_columns(src);
    st_logInfo("tui-chain: source .tui has T = %" PRIi64 " universal columns; "
               "maxGap = %" PRIi64 "\n", T_src, max_gap);

    /* Load the genome roster from the source's g-records. */
    int64_t n_genomes = 0;
    TuiGenomeInfo *roster = tui_genome_names(src_tui_path, &n_genomes);
    if (roster == NULL || n_genomes == 0) {
        fprintf(stderr,
                "tui-chain: source .tui has no g-records "
                "(corrupt, or not built by `taffy index -u`?)\n");
        tui_destruct(src); free(src_tui_path); free(out_path_owned);
        return 1;
    }
    st_logInfo("tui-chain: %" PRIi64 " genomes from .tui g-records\n", n_genomes);

    /* Build the directory (name-sorted, same as source). */
    stHash *seqlens = tui_sequence_lengths(src_tui_path);
    if (seqlens == NULL) {
        fprintf(stderr, "tui-chain: tui_sequence_lengths failed\n");
        tui_genome_info_free(roster, n_genomes);
        tui_destruct(src); free(src_tui_path); free(out_path_owned);
        return 1;
    }
    int64_t n_seqs = stHash_size(seqlens);
    DirEntry *dir = (DirEntry *)st_malloc((size_t)n_seqs * sizeof(DirEntry));
    int64_t di = 0;
    stHashIterator *sit = stHash_getIterator(seqlens);
    char *seqname;
    while ((seqname = (char *)stHash_getNext(sit)) != NULL) {
        dir[di].name = seqname;
        dir[di].length = (int64_t)(intptr_t)stHash_search(seqlens, seqname);
        di++;
    }
    stHash_destructIterator(sit);
    qsort(dir, (size_t)n_seqs, sizeof(DirEntry), dir_entry_cmp);
    for (int64_t i = 0; i < n_seqs; i++) dir[i].dir_idx = i;

    stHash *seq_to_dir_idx = stHash_construct3(
        stHash_stringKey, stHash_stringEqualKey, NULL, free);
    for (int64_t i = 0; i < n_seqs; i++) {
        int64_t *p = (int64_t *)st_malloc(sizeof(int64_t));
        *p = i;
        stHash_insert(seq_to_dir_idx, dir[i].name, p);
    }

    /* Open output .tui. */
    char blurb[256];
    snprintf(blurb, sizeof(blurb),
             "chained from %s with maxGap=%" PRIi64, src_tui_path, max_gap);
    OneSchema *schema = NULL;
    OneFile *of = tui_open_write(out_path, "taffy", "tui-chain", blurb, &schema);
    if (of == NULL) {
        stHash_destruct(seq_to_dir_idx); free(dir);
        stHash_destruct(seqlens);
        tui_genome_info_free(roster, n_genomes);
        tui_destruct(src); free(src_tui_path); free(out_path_owned);
        tui_atexit_disarm();
        return 1;
    }

    /* t: total columns + format version (chaining preserves universal coords). */
    tui_write_header(of, T_src);

    /* X: copy the source's X-record (universal-column -> .maf.gz file
     * position).  Chain merging doesn't change universal columns, so the
     * source positions remain valid.  Lets `taffy view -r` work against
     * the chained .tui as long as the original .maf.gz is reachable. */
    {
        int64_t xn = tui_idx_n(src);
        if (xn > 0) {
            int64_t x_raw = 0, x_def = 0;
            uint8_t *xdef = tui_encode_idx(tui_idx_cols(src), tui_idx_fpos(src),
                                           xn, &x_raw, &x_def);
            oneInt(of, 0) = x_raw;
            oneInt(of, 1) = xn;
            oneWriteLine(of, 'X', x_def, xdef);
            free(xdef);
        } else {
            /* Source had no X-track; emit a stub so the schema-required
             * X line is present.  view -r will fail, lift will work. */
            int64_t c0 = 0, f0 = 0, x_raw = 0, x_def = 0;
            uint8_t *xdef = tui_encode_idx(&c0, &f0, 1, &x_raw, &x_def);
            oneInt(of, 0) = x_raw;
            oneInt(of, 1) = 1;
            oneWriteLine(of, 'X', x_def, xdef);
            free(xdef);
        }
    }

    /* capture: full_name -> s_ord, used by deferred d-record writer.
     * See the deferred d: directory writer below. */
    stHash *capture = stHash_construct3(
        stHash_stringKey, stHash_stringEqualKey, free, free);

    int64_t c_ord_emit = 0;
    int64_t s_ord_emit = 0;
    int64_t total_in   = 0, total_out = 0;
    time_t per_genome_t0 = time(NULL);
    for (int64_t gi = 0; gi < n_genomes; gi++) {
        const char *gname = roster[gi].name;
        time_t g_t0 = time(NULL);

        TuiGenomeLift *gl = tui_genome_lift_load(src, gname);
        if (gl == NULL) {
            st_logInfo("tui-chain: skip %s (no chunks)\n", gname);
            continue;
        }

        ChainCtx cx = {
            .max_gap = max_gap,
            .active  = stHash_construct(),
            .out     = stList_construct(),
            .n_in    = 0,
            .n_out   = 0,
        };
        tui_genome_lift_stream_runs(gl, chain_visit_cb, &cx);
        chain_flush_active(&cx);

        int64_t n_chains = stList_length(cx.out);
        if (n_chains == 0) {
            stList_destruct(cx.out);
            stHash_destruct(cx.active);
            tui_genome_lift_destruct(gl);
            st_logInfo("tui-chain: %s -> 0 chains\n", gname);
            continue;
        }

        /* bucket_by_seq reads cr->seq (gl-owned interned strings) via
         * strcmp; gl must stay alive until bucket completes. */
        SeqGroup *groups = NULL;
        int64_t n_groups = bucket_by_seq(cx.out, gname, seq_to_dir_idx, &groups);
        stList_destruct(cx.out);
        stHash_destruct(cx.active);
        tui_genome_lift_destruct(gl);

        write_genome(of, groups, n_groups, &c_ord_emit, &s_ord_emit, capture);
        /* write_genome frees each ChainRun struct + the groups[].runs
         * array + groups[].full_name + groups[].  No leftovers to clean. */

        total_in  += cx.n_in;
        total_out += cx.n_out;

        malloc_trim(0);
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        st_logInfo("tui-chain: %s: %" PRIi64 " runs in -> %" PRIi64
                   " chains out (%.1fx) in %" PRIi64 " s [rss %ld MB]\n",
                   gname, cx.n_in, cx.n_out,
                   cx.n_in > 0 ? (double)cx.n_in / (double)cx.n_out : 0.0,
                   (int64_t)(time(NULL) - g_t0),
                   ru.ru_maxrss / 1024);
    }
    st_logInfo("tui-chain: per-genome work done in %" PRIi64 " s; "
               "%" PRIi64 " runs in -> %" PRIi64 " chains out (%.1fx)\n",
               (int64_t)(time(NULL) - per_genome_t0),
               total_in, total_out,
               total_in > 0 ? (double)total_in / (double)total_out : 0.0);

    /* d: directory.  Written AFTER S/C/R, same deferred-pointer pattern as
     * the capture map built above.  active seqs get
     * their captured s_ord; inactive get sentinel -1. */
    int64_t active_seqs = 0;
    for (int64_t i = 0; i < n_seqs; i++) {
        int64_t *ord_p = (int64_t *)stHash_search(capture, dir[i].name);
        int64_t stored_ord = ord_p ? (*ord_p - 1) : (int64_t)(-1);
        oneInt(of, 1) = stored_ord;
        oneInt(of, 2) = dir[i].length;
        oneWriteLine(of, 'd', strlen(dir[i].name), (void *)dir[i].name);
        if (ord_p) active_seqs++;
    }
    st_logInfo("tui-chain: wrote %" PRIi64 " d-records (%" PRIi64
               " active, %" PRIi64 " inactive)\n",
               (int64_t)n_seqs, active_seqs, (int64_t)n_seqs - active_seqs);

    /* g: roster (preserve from source/list). */
    for (int64_t gi = 0; gi < n_genomes; gi++) {
        oneInt(of, 1) = roster[gi].total_bp;
        oneInt(of, 2) = roster[gi].n_chroms;
        oneWriteLine(of, 'g', strlen(roster[gi].name), (void *)roster[gi].name);
    }

    stHash_destruct(capture);
    oneFileClose(of);
    oneSchemaDestroy(schema);
    tui_atexit_disarm();

    stHash_destruct(seq_to_dir_idx);
    free(dir);
    stHash_destruct(seqlens);
    tui_genome_info_free(roster, n_genomes);
    tui_destruct(src);
    free(src_tui_path);
    free(out_path_owned);

    st_logInfo("tui-chain: total wall %" PRIi64 " s; wrote %s\n",
               (int64_t)(time(NULL) - startTime), out_path);
    return 0;
}
