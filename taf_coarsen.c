/*
 * taffy tui-coarsen: build a level-of-detail (LOD) .tui from an existing
 * .tui by aggregating runs into per-(g_seq, g_bin, uc_bin) bp cells.
 *
 * The output is a regular .tui (same schema, same reader): each output
 * "run" encodes one cell as
 *     t_start  = g_bin_idx * B        (position on g_seq in g_bp)
 *     g_start  = uc_bin_idx * B       (position in universal columns)
 *     length   = bp_in_cell           (coverage; usually <= B)
 *     strand   = '+'                  (LOD collapses strand)
 *
 * Pair the output with a HAL-style .lod.txt mapping (span-threshold ->
 * file) and let the browser pick the right LOD per query span; blockViz
 * needs no changes to read these.
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
        "taffy tui-coarsen [options] -i <input.maf|input.tui> -B <binSize>\n"
        "\n"
        "Build a level-of-detail (LOD) .tui from an existing .tui by\n"
        "aggregating runs into per-(g_seq, g_bin, uc_bin) bp cells.\n"
        "The output is a regular .tui readable by every taffy reader\n"
        "(view, lift, blockViz); same schema, same code path.  Pair with\n"
        "a HAL-style .lod.txt for browser LOD selection.\n"
        "\n"
        "Options:\n"
        "  -i --inputFile FILE     Base MAF/TAF whose <input>.tui will be\n"
        "                          coarsened, OR with --fromTui, an existing\n"
        "                          .tui directly.  REQUIRED.\n"
        "  -o --outputFile FILE    Output LOD .tui path.\n"
        "                          DEFAULT: <input>.lod<B>.tui\n"
        "  -B --binSize INT        LOD bin width in bp.  Applied symmetrically\n"
        "                          to BOTH the g_seq axis and the universal-\n"
        "                          column axis.  REQUIRED.  Typical values:\n"
        "                          1000, 10000, 100000, 1000000.\n"
        "     --minCellBp INT      Noise floor: drop cells with bp < this.\n"
        "                          DEFAULT 1 (keep everything).\n"
        "     --fromTui            Treat -i as the source .tui itself (for\n"
        "                          chaining a coarser LOD from a finer one).\n"
        "     --lodTxt PATH        Append/update a HAL-style \"<B>\\t<path>\\n\"\n"
        "                          line at PATH (replaces any existing line\n"
        "                          for the same B).\n"
        "     --genomeList FILE    Read genome names from FILE (one per line;\n"
        "                          '#' comments OK).  Overrides the g-record\n"
        "                          roster requirement -- use this when the\n"
        "                          source .tui predates the g-record schema\n"
        "                          (b8f6811) and rebuilding isn't practical.\n"
        "                          Each name MUST be a real d-line prefix on\n"
        "                          the source .tui or that genome is skipped.\n"
        "                          TRANSITIONAL flag; remove once all .tui\n"
        "                          files have g records.\n"
        "     --fillBins           Emit length=B for every output cell instead\n"
        "                          of length=bp_coverage.  Each aligned bin\n"
        "                          renders as a full bin-width block in a\n"
        "                          snake viewer; the per-bin coverage signal\n"
        "                          in length is lost but the visual continuity\n"
        "                          improves at zoom-out.  Overrides\n"
        "                          --minBinFrac.\n"
        "     --minBinFrac F       For each cell emit length=max(bp, B*F).\n"
        "                          F in [0.0, 1.0]; 0.0 = pure coverage\n"
        "                          (default), 1.0 = same as --fillBins.\n"
        "                          Use a small fraction (e.g. 0.25) to keep\n"
        "                          the coverage signal but give thin bins a\n"
        "                          visible minimum width.  DEFAULT 0.0.\n"
        "     --maxParalogsPerBin K\n"
        "                          For each (seq, uc_bin), keep at most K\n"
        "                          cells (the K highest-bp).  At LOD, source\n"
        "                          paralogs scattered across many g_bins\n"
        "                          share a uc_bin and stack vertically in\n"
        "                          the snake renderer; this cap bounds the\n"
        "                          worst-case stack depth.  Pair with\n"
        "                          --fillBins for clean bin-aligned rows.\n"
        "                          DEFAULT 0 (no cap).\n"
        "  -l --logLevel STR       critical|info|debug.  DEFAULT critical.\n"
        "  -h --help               Print this help.\n");
}

/* ------------------------------------------------------------------ */
/* Cell hash: (seq, g_bin, uc_bin) -> bp                              */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *seq;       /* gl-owned, pointer-stable */
    int64_t     g_bin;
    int64_t     uc_bin;
} CellKey;

static uint64_t cell_key_hash(const void *k) {
    const CellKey *c = (const CellKey *)k;
    /* Fibonacci-style mix; same pattern as taf_lift.c BinKey. */
    return ((uint64_t)(uintptr_t)c->seq) * 0x9E3779B97F4A7C15ULL
         ^ ((uint64_t)c->g_bin * 0xC2B2AE3D27D4EB4FULL)
         ^ ((uint64_t)c->uc_bin);
}

static int cell_key_equal(const void *a, const void *b) {
    const CellKey *p = (const CellKey *)a;
    const CellKey *q = (const CellKey *)b;
    return p->seq == q->seq && p->g_bin == q->g_bin && p->uc_bin == q->uc_bin;
}

/* ------------------------------------------------------------------ */
/* Run splitter: emit per-(seq, g_bin, uc_bin) cell with bp           */
/* ------------------------------------------------------------------ */

typedef struct {
    int64_t  B;
    stHash  *cells;          /* CellKey* -> int64_t* bp */
} CoarsenCtx;

static inline void cell_add(CoarsenCtx *cx, const char *seq,
                            int64_t g_bin, int64_t uc_bin, int64_t bp) {
    CellKey probe = { seq, g_bin, uc_bin };
    int64_t *bp_p = (int64_t *)stHash_search(cx->cells, &probe);
    if (bp_p != NULL) { *bp_p += bp; return; }
    CellKey *k = (CellKey *)st_malloc(sizeof(CellKey));
    *k = probe;
    int64_t *v = (int64_t *)st_malloc(sizeof(int64_t));
    *v = bp;
    stHash_insert(cx->cells, k, v);
}

/* Split a single TuiRun across both uc_bin AND g_bin boundaries; emit
 * one cell update per sub-segment. */
static void coarsen_visit_cb(const TuiRun *r, void *user) {
    CoarsenCtx *cx = (CoarsenCtx *)user;
    int64_t B = cx->B;
    int64_t uc_lo = r->g_start;
    int64_t uc_hi = r->g_start + r->length;
    int64_t c = uc_lo;
    while (c < uc_hi) {
        int64_t uc_bin = c / B;
        int64_t uc_bin_end = (uc_bin + 1) * B;
        int64_t c_uc_end = uc_hi < uc_bin_end ? uc_hi : uc_bin_end;
        /* Within this uc_bin slice, walk g_bin sub-segments. */
        int64_t cur = c;
        while (cur < c_uc_end) {
            int64_t t_pos = r->strand
                ? r->t_start + (cur - r->g_start)
                : r->t_start + r->length - 1 - (cur - r->g_start);
            int64_t g_bin = t_pos / B;
            int64_t bp_left_in_g_bin;
            if (r->strand) {
                /* + strand: t_pos increases with cur. */
                int64_t g_bin_end_pos = (g_bin + 1) * B;
                bp_left_in_g_bin = g_bin_end_pos - t_pos;
            } else {
                /* - strand: t_pos decreases as cur advances. */
                int64_t g_bin_start_pos = g_bin * B;
                bp_left_in_g_bin = t_pos - g_bin_start_pos + 1;
            }
            int64_t cur_end = cur + bp_left_in_g_bin;
            if (cur_end > c_uc_end) cur_end = c_uc_end;
            int64_t bp = cur_end - cur;
            if (bp > 0) cell_add(cx, r->seq, g_bin, uc_bin, bp);
            cur = cur_end;
        }
        c = c_uc_end;
    }
}

/* ------------------------------------------------------------------ */
/* Paralog dedup (optional, --maxParalogsPerBin)                       */
/* ------------------------------------------------------------------ */

/* Walk the cells hash; within each (seq, uc_bin) group keep only the
 * top max_per_bin cells by bp; free the rest.
 *
 * Motivation: at LOD, source paralogs that get fragmented across many
 * g_bins all share the same uc_bin and emit as separate cells.  In a
 * snake viewer these stack vertically -- one row per cell at the
 * crowded uc_bin -- which the user can't render gracefully.  Capping
 * at K limits worst-case stack height to ~K + cross-seq paralog count.
 *
 * Runs in O(N log N) on cell count.  Walks the hash once into an array,
 * sorts by (seq, uc_bin, -bp), then sweeps groups and frees beyond K.
 */
typedef struct {
    CellKey *k;
    int64_t  bp;
} DedupEntry;

static int dedup_entry_cmp(const void *a, const void *b) {
    const DedupEntry *p = (const DedupEntry *)a;
    const DedupEntry *q = (const DedupEntry *)b;
    /* (seq ptr, uc_bin asc, bp desc) */
    if (p->k->seq < q->k->seq) return -1;
    if (p->k->seq > q->k->seq) return 1;
    if (p->k->uc_bin < q->k->uc_bin) return -1;
    if (p->k->uc_bin > q->k->uc_bin) return 1;
    if (p->bp > q->bp) return -1;
    if (p->bp < q->bp) return 1;
    return 0;
}

static int64_t cell_dedup_paralogs(stHash *cells, int64_t max_per_bin) {
    int64_t n = stHash_size(cells);
    if (n == 0 || max_per_bin <= 0) return 0;
    DedupEntry *arr = (DedupEntry *)st_malloc((size_t)n * sizeof(DedupEntry));
    stHashIterator *it = stHash_getIterator(cells);
    int64_t i = 0;
    CellKey *k;
    while ((k = (CellKey *)stHash_getNext(it)) != NULL) {
        int64_t *bp_p = (int64_t *)stHash_search(cells, k);
        arr[i].k  = k;
        arr[i].bp = *bp_p;
        i++;
    }
    stHash_destructIterator(it);
    qsort(arr, (size_t)n, sizeof(DedupEntry), dedup_entry_cmp);
    int64_t removed = 0;
    int64_t j = 0;
    while (j < n) {
        int64_t end = j + 1;
        while (end < n && arr[end].k->seq    == arr[j].k->seq
                       && arr[end].k->uc_bin == arr[j].k->uc_bin) end++;
        /* Group is [j, end).  Keep first max_per_bin (already sorted by
         * bp desc within group); drop the rest. */
        for (int64_t m = j + max_per_bin; m < end; m++) {
            /* removeAndFreeKey calls the registered key destructor (free),
             * and returns the value pointer for us to free. */
            int64_t *v = (int64_t *)stHash_removeAndFreeKey(cells, arr[m].k);
            free(v);
            removed++;
        }
        j = end;
    }
    free(arr);
    return removed;
}

/* ------------------------------------------------------------------ */
/* Per-genome emit: sort cells -> chunk -> encode -> write S/C/R     */
/* ------------------------------------------------------------------ */

typedef struct {
    int64_t g_bin;
    int64_t uc_bin;
    int64_t bp;
} EmitCell;

typedef struct {
    char       *full_name;        /* "<genome>.<seq>" -- owns; freed in write_genome */
    int64_t     seq_idx_in_dir;   /* parent_S_ord = seq_idx_in_dir + 1 */
    EmitCell   *cells;
    int64_t     n_cells;
} SeqGroup;

static int seq_group_cmp(const void *a, const void *b) {
    const SeqGroup *p = (const SeqGroup *)a;
    const SeqGroup *q = (const SeqGroup *)b;
    if (p->seq_idx_in_dir < q->seq_idx_in_dir) return -1;
    if (p->seq_idx_in_dir > q->seq_idx_in_dir) return  1;
    return 0;
}

static int emit_cell_cmp(const void *a, const void *b) {
    const EmitCell *p = (const EmitCell *)a;
    const EmitCell *q = (const EmitCell *)b;
    if (p->uc_bin < q->uc_bin) return -1;
    if (p->uc_bin > q->uc_bin) return  1;
    if (p->g_bin  < q->g_bin)  return -1;
    if (p->g_bin  > q->g_bin)  return  1;
    return 0;
}

/* Drain the cell hash for one genome into per-seq sorted groups.  Caller
 * frees groups + cells.  Returns n_groups; *groups_out is malloc'd.
 *
 * Note: cell keys hold the BARE seq name (gl strips the "genome." prefix
 * at load -- tui.c:2020), but the global directory map is keyed by the
 * FULL "genome.seq" d-line name.  Reconstruct the full key here. */
static int64_t drain_cells(stHash *cells, int64_t min_cell_bp,
                           const char *genome_name,
                           stHash *seq_to_dir_idx,
                           SeqGroup **groups_out) {
    /* Bucket cells by seq pointer. */
    /* Value destructor frees the per-seq stList, which in turn frees its
     * EmitCells (constructed with stList_construct3(0, free) below).
     * Plain stHash_construct() registers NULL destructors and leaks every
     * EmitCell per genome (~40 B/cell × millions of cells across a run). */
    stHash *by_seq = stHash_construct2(NULL, (void (*)(void *))stList_destruct);
    stHashIterator *it = stHash_getIterator(cells);
    CellKey *k;
    while ((k = (CellKey *)stHash_getNext(it)) != NULL) {
        int64_t *bp_p = (int64_t *)stHash_search(cells, k);
        if (*bp_p < min_cell_bp) continue;
        stList *l = (stList *)stHash_search(by_seq, (void *)k->seq);
        if (l == NULL) {
            l = stList_construct3(0, free);
            stHash_insert(by_seq, (void *)k->seq, l);
        }
        EmitCell *e = (EmitCell *)st_malloc(sizeof(EmitCell));
        e->g_bin = k->g_bin; e->uc_bin = k->uc_bin; e->bp = *bp_p;
        stList_append(l, e);
    }
    stHash_destructIterator(it);

    /* Materialize SeqGroup[] */
    int64_t n_g = stHash_size(by_seq);
    SeqGroup *groups = (SeqGroup *)st_malloc((size_t)(n_g ? n_g : 1) * sizeof(SeqGroup));
    int64_t gi = 0;
    it = stHash_getIterator(by_seq);
    void *seq_key;
    while ((seq_key = stHash_getNext(it)) != NULL) {
        stList *l = (stList *)stHash_search(by_seq, seq_key);
        int64_t n_cells = stList_length(l);
        EmitCell *arr = (EmitCell *)st_malloc((size_t)n_cells * sizeof(EmitCell));
        for (int64_t i = 0; i < n_cells; i++) arr[i] = *(EmitCell *)stList_get(l, i);
        qsort(arr, (size_t)n_cells, sizeof(EmitCell), emit_cell_cmp);
        /* Build the full "genome.seq" name (gl strips the prefix in
         * r->seq; reader expects full names in S records to match
         * d-line entries). */
        char *full = (char *)st_malloc(strlen(genome_name) + 1 +
                                       strlen((const char *)seq_key) + 1);
        sprintf(full, "%s.%s", genome_name, (const char *)seq_key);
        int64_t *idx_p = (int64_t *)stHash_search(seq_to_dir_idx, full);
        if (idx_p == NULL) {
            st_errAbort("tui-coarsen: seq %s not in directory map (this is a bug)",
                        full);
        }
        groups[gi].full_name = full;
        groups[gi].seq_idx_in_dir = *idx_p;
        groups[gi].cells = arr;
        groups[gi].n_cells = n_cells;
        gi++;
    }
    stHash_destructIterator(it);
    stHash_destruct(by_seq);
    qsort(groups, (size_t)n_g, sizeof(SeqGroup), seq_group_cmp);
    *groups_out = groups;
    return n_g;
}

/* Write one genome's S/C/R records.  Mirrors phase2_genome_write at
 * tui.c:1131-1158.  Builds chunks honoring TUI_CHUNK_RUNS and
 * TUI_CHUNK_G_MAX caps so the reader's per-chunk t-range skip works.
 *
 * Output-length policy: by default length=bp (the raw coverage in
 * the cell).  fill_bins=1 emits length=B (full bin-width).  Otherwise
 * length=max(bp, B*min_bin_frac).  See --fillBins / --minBinFrac help.
 *
 * s_ord_counter is the running tally of S-records written across the
 * whole file (1-indexed; incremented before each S-record).  capture
 * (full_name -> int64_t* s_ord) records the actual s_ord for each
 * written S, used later to populate the d-record's parent-S pointer
 * with the real ordinal (not the dir_index, which doesn't match the
 * genome-roster emission order). */
static void write_genome(OneFile *of, SeqGroup *groups, int64_t n_groups,
                         int64_t B, int64_t *c_ord_emit,
                         int64_t *s_ord_counter, stHash *capture,
                         int fill_bins, double min_bin_frac) {
    int64_t min_len = fill_bins ? B : (int64_t)(B * min_bin_frac);
    if (min_len > B) min_len = B;
    if (min_len < 0) min_len = 0;
    for (int64_t gi = 0; gi < n_groups; gi++) {
        SeqGroup *g = &groups[gi];
        /* S line: seq name + length-in-coarse-units = number of g_bins.
         * Reader treats this as the "seq length" field.  Coarse units. */
        int64_t max_g_bin = 0;
        for (int64_t i = 0; i < g->n_cells; i++)
            if (g->cells[i].g_bin > max_g_bin) max_g_bin = g->cells[i].g_bin;
        oneInt(of, 1) = (max_g_bin + 1) * B;   /* approximate seq_len upper bound */
        ++(*s_ord_counter);                    /* this S becomes ordinal *s_ord_counter */
        oneWriteLine(of, 'S', strlen(g->full_name), (void *)g->full_name);
        int64_t s_ord = *s_ord_counter;
        /* Capture (full_name -> s_ord) so the d-record writer (post-loop)
         * can point this seq's d-line at the right S-record.  The dir
         * index does NOT match s_ord here: S-records come out in
         * genome-roster order, while d-records are name-sorted globally. */
        if (capture != NULL) {
            int64_t *p = (int64_t *)st_malloc(sizeof(int64_t));
            *p = s_ord;
            stHash_insert(capture, stString_copy(g->full_name), p);
        }
        int64_t i = 0;
        while (i < g->n_cells) {
            /* Build one chunk's triples.  cells[] is sorted by uc_bin
             * then g_bin -- the "g_start" of a coarse run is uc_bin * B,
             * monotonically non-decreasing.  Chunk closes on RUNS cap or
             * when (uc_bin_max - uc_bin_min) * B would exceed
             * TUI_CHUNK_G_MAX. */
            int64_t start = i;
            int64_t cm = 0;
            int64_t uc_bin_min = g->cells[start].uc_bin;
            int64_t uc_bin_max = uc_bin_min;
            int64_t g_bin_min = g->cells[start].g_bin;
            int64_t g_bin_max = g_bin_min;
            while (i < g->n_cells && cm < TUI_CHUNK_RUNS) {
                int64_t new_uc_max = g->cells[i].uc_bin > uc_bin_max
                                     ? g->cells[i].uc_bin : uc_bin_max;
                if ((new_uc_max - uc_bin_min) * B > TUI_CHUNK_G_MAX) break;
                if (g->cells[i].g_bin > g_bin_max) g_bin_max = g->cells[i].g_bin;
                if (g->cells[i].g_bin < g_bin_min) g_bin_min = g->cells[i].g_bin;
                uc_bin_max = new_uc_max;
                cm++; i++;
            }
            /* Encode triples. */
            int64_t *buf = (int64_t *)st_malloc((size_t)cm * 3 * sizeof(int64_t));
            for (int64_t k = 0; k < cm; k++) {
                EmitCell *e = &g->cells[start + k];
                int64_t len = e->bp;
                if (len < min_len) len = min_len;
                buf[3*k + 0] = e->g_bin * B;     /* t_start */
                buf[3*k + 1] = e->uc_bin * B;    /* g_start */
                buf[3*k + 2] = (len << 1) | 0;   /* lenc: + strand */
            }
            int64_t raw_len = 0, def_len = 0;
            uint8_t *def = tui_encode_runs(buf, cm, &raw_len, &def_len);
            free(buf);

            ++(*c_ord_emit);
            oneInt(of, 0) = uc_bin_min * B;
            oneInt(of, 1) = uc_bin_max * B + B;     /* g_max exclusive-ish */
            oneInt(of, 2) = s_ord;
            oneInt(of, 3) = *c_ord_emit;
            oneInt(of, 4) = g_bin_min * B;          /* t_min */
            oneInt(of, 5) = g_bin_max * B + B;      /* t_max */
            oneWriteLine(of, 'C', 0, NULL);
            oneInt(of, 0) = raw_len;
            oneWriteLine(of, 'R', def_len, def);
            free(def);
        }
        free(g->cells);
        free(g->full_name);
    }
    free(groups);
}

/* ------------------------------------------------------------------ */
/* .lod.txt update                                                    */
/* ------------------------------------------------------------------ */

/* Append or replace the line for `bin_size` in PATH.  Format:
 *   <bin_size>\t<lod_path>\n
 * Atomic via temp+rename.  Existing lines for other bin_sizes preserved. */
static void update_lod_txt(const char *path, int64_t bin_size, const char *lod_path) {
    char tmp[4096];
    snprintf(tmp, sizeof(tmp), "%s.tmp.%d", path, (int)getpid());
    FILE *out = fopen(tmp, "w");
    if (out == NULL) {
        fprintf(stderr, "tui-coarsen: cannot write %s: %s\n", tmp, strerror(errno));
        return;
    }
    int wrote_this = 0;
    FILE *in = fopen(path, "r");
    if (in != NULL) {
        char *line = NULL; size_t cap = 0; ssize_t n;
        while ((n = getline(&line, &cap, in)) > 0) {
            int64_t b = 0;
            if (sscanf(line, "%" SCNi64, &b) == 1 && b == bin_size) {
                fprintf(out, "%" PRIi64 "\t%s\n", bin_size, lod_path);
                wrote_this = 1;
            } else {
                fputs(line, out);
            }
        }
        free(line);
        fclose(in);
    }
    if (!wrote_this) fprintf(out, "%" PRIi64 "\t%s\n", bin_size, lod_path);
    fclose(out);
    if (rename(tmp, path) != 0) {
        fprintf(stderr, "tui-coarsen: rename %s -> %s failed: %s\n",
                tmp, path, strerror(errno));
    }
}

/* ------------------------------------------------------------------ */
/* Build per-genome seq directory: name -> dir-index                  */
/* ------------------------------------------------------------------ */

typedef struct {
    char    *name;        /* full "genome.seq" key */
    int64_t  length;
    int64_t  dir_idx;     /* position in name-sorted directory */
} DirEntry;

static int dir_entry_cmp(const void *a, const void *b) {
    return strcmp(((const DirEntry *)a)->name, ((const DirEntry *)b)->name);
}

/* ------------------------------------------------------------------ */
/* Main                                                               */
/* ------------------------------------------------------------------ */

int taf_coarsen_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    char    *input_file     = NULL;
    char    *output_file    = NULL;
    int64_t  bin_size       = 0;
    int64_t  min_cell_bp    = 1;
    int      from_tui       = 0;
    char    *lod_txt        = NULL;
    char    *genome_list    = NULL;
    char    *log_level      = NULL;
    /* Output-length policy.  Each cell carries the raw bp of source
     * alignment that fell in (g_bin, uc_bin).  Default emits length=bp,
     * which makes a partially-aligned bin render as a small block at
     * the bin start in a snake viewer -- visually sparse.  --fillBins
     * emits length=B so every cell with any alignment renders as a
     * full bin-width block (lose coverage signal, gain visual
     * continuity).  --minBinFrac F emits length=max(bp, B*F): keeps
     * the variable-length signal but enforces a visible minimum.
     * Mutually exclusive in effect: fillBins overrides minBinFrac. */
    int      fill_bins      = 0;
    double   min_bin_frac   = 0.0;
    /* Per-uc_bin paralog cap.  At LOD, source paralogs that get scattered
     * across many g_bins all share the same uc_bin and stack vertically
     * in a snake viewer (each pile is one row of the chain).  This cap
     * keeps at most the top-K cells (by bp) per (seq, uc_bin), removing
     * the LOD-introduced paralog explosion.  0 = no cap (default). */
    int64_t  max_paralogs_per_bin = 0;

    enum { OPT_MIN_CELL_BP = 256, OPT_FROM_TUI, OPT_LOD_TXT, OPT_GENOME_LIST,
           OPT_FILL_BINS, OPT_MIN_BIN_FRAC, OPT_MAX_PARALOGS_PER_BIN };
    static struct option long_options[] = {
        { "inputFile",          required_argument, 0, 'i' },
        { "outputFile",         required_argument, 0, 'o' },
        { "binSize",            required_argument, 0, 'B' },
        { "minCellBp",          required_argument, 0, OPT_MIN_CELL_BP },
        { "fromTui",            no_argument,       0, OPT_FROM_TUI },
        { "lodTxt",             required_argument, 0, OPT_LOD_TXT },
        { "genomeList",         required_argument, 0, OPT_GENOME_LIST },
        { "fillBins",           no_argument,       0, OPT_FILL_BINS },
        { "minBinFrac",         required_argument, 0, OPT_MIN_BIN_FRAC },
        { "maxParalogsPerBin",  required_argument, 0, OPT_MAX_PARALOGS_PER_BIN },
        { "logLevel",           required_argument, 0, 'l' },
        { "help",               no_argument,       0, 'h' },
        { 0, 0, 0, 0 }
    };

    while (1) {
        int oi = 0;
        int k = getopt_long(argc, argv, "i:o:B:l:h", long_options, &oi);
        if (k == -1) break;
        switch (k) {
            case 'i': input_file  = optarg; break;
            case 'o': output_file = optarg; break;
            case 'B': bin_size    = atoll(optarg); break;
            case OPT_MIN_CELL_BP:  min_cell_bp = atoll(optarg); break;
            case OPT_FROM_TUI:     from_tui = 1; break;
            case OPT_LOD_TXT:      lod_txt = optarg; break;
            case OPT_GENOME_LIST:  genome_list = optarg; break;
            case OPT_FILL_BINS:    fill_bins = 1; break;
            case OPT_MIN_BIN_FRAC: min_bin_frac = atof(optarg); break;
            case OPT_MAX_PARALOGS_PER_BIN: max_paralogs_per_bin = atoll(optarg); break;
            case 'l': log_level   = optarg; break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }
    if (input_file == NULL || bin_size <= 0) {
        if (input_file == NULL) fprintf(stderr, "ERROR: -i required\n");
        if (bin_size   <= 0)    fprintf(stderr, "ERROR: -B (positive) required\n");
        usage();
        return 1;
    }
    if (min_cell_bp < 1) {
        fprintf(stderr, "ERROR: --minCellBp must be >= 1\n");
        return 1;
    }
    if (min_bin_frac < 0.0 || min_bin_frac > 1.0) {
        fprintf(stderr, "ERROR: --minBinFrac must be in [0.0, 1.0]\n");
        return 1;
    }
    if (max_paralogs_per_bin < 0) {
        fprintf(stderr, "ERROR: --maxParalogsPerBin must be >= 0\n");
        return 1;
    }
    st_setLogLevelFromString(log_level);

    /* Resolve source .tui path. */
    char *src_tui_path;
    if (from_tui) {
        src_tui_path = stString_copy(input_file);
    } else {
        src_tui_path = tui_path(input_file);
    }
    /* Default output: <input>.lod<B>.tui */
    char *out_path_owned = NULL;
    const char *out_path;
    if (output_file != NULL) {
        out_path = output_file;
    } else {
        out_path_owned = (char *)st_malloc(strlen(input_file) + 64);
        sprintf(out_path_owned, "%s.lod%" PRIi64 ".tui", input_file, bin_size);
        out_path = out_path_owned;
    }

    Tui *src = tui_load(src_tui_path);
    if (src == NULL) {
        fprintf(stderr, "tui-coarsen: cannot open %s\n", src_tui_path);
        free(src_tui_path); free(out_path_owned);
        return 1;
    }
    int64_t T_src = tui_total_columns(src);
    st_logInfo("tui-coarsen: source .tui has T = %" PRIi64 " universal columns; "
               "bin size B = %" PRIi64 "\n", T_src, bin_size);

    /* Load genome roster.  Two sources, in order of preference:
     *   1. --genomeList FILE  (transitional override; one name per line,
     *      '#' comments OK).  Use when the .tui predates the g-record
     *      schema and rebuilding isn't practical.  Each name must be a
     *      real d-line prefix; mismatches are silently skipped.
     *   2. The .tui's g records (the proper long-term path).
     * If neither is available we refuse to coarsen -- guessing genome
     * boundaries from d-line first-dot heuristics breaks on dotted
     * accessions (GCA_028885655.2 vs GCA_028885655). */
    int64_t n_genomes = 0;
    TuiGenomeInfo *roster = NULL;
    if (genome_list != NULL) {
        /* Parse the genome list file: one name per line; skip blanks/
         * '#' comments.  total_bp / n_chroms are unknown here (we only
         * need names for coarsening), so populate as 0 -- they're only
         * used when copying the g-record roster forward into the LOD
         * file, and a missing g records source means the LOD won't have
         * accurate totals either way (acceptable for the transitional
         * path; documented in --genomeList help text). */
        FILE *gf = fopen(genome_list, "r");
        if (gf == NULL) {
            fprintf(stderr, "tui-coarsen: cannot read --genomeList %s: %s\n",
                    genome_list, strerror(errno));
            tui_destruct(src);
            free(src_tui_path); free(out_path_owned);
            return 1;
        }
        int64_t cap = 32;
        roster = (TuiGenomeInfo *)st_malloc((size_t)cap * sizeof(TuiGenomeInfo));
        char *line = NULL; size_t lcap = 0; ssize_t got;
        while ((got = getline(&line, &lcap, gf)) > 0) {
            while (got > 0 && (line[got-1] == '\n' || line[got-1] == '\r'
                              || line[got-1] == ' '  || line[got-1] == '\t')) {
                line[--got] = 0;
            }
            char *s = line;
            while (*s == ' ' || *s == '\t') s++;
            if (*s == 0 || *s == '#') continue;
            if (n_genomes == cap) {
                cap *= 2;
                roster = (TuiGenomeInfo *)st_realloc(roster,
                                                    (size_t)cap * sizeof(TuiGenomeInfo));
            }
            roster[n_genomes].name     = stString_copy(s);
            roster[n_genomes].total_bp = 0;
            roster[n_genomes].n_chroms = 0;
            n_genomes++;
        }
        free(line);
        fclose(gf);
        st_logInfo("tui-coarsen: %" PRIi64 " genomes from --genomeList %s\n",
                   n_genomes, genome_list);
    } else {
        roster = tui_genome_names(src_tui_path, &n_genomes);
        if (roster == NULL || n_genomes == 0) {
            fprintf(stderr,
                    "tui-coarsen: source .tui has no genome roster (g records).\n"
                    "  Either rebuild the .tui with a newer taffy (commit b8f6811\n"
                    "  or later) for the proper genome roster, OR pass\n"
                    "  --genomeList <file> with one genome name per line as a\n"
                    "  transitional workaround.  Refusing to coarsen without one.\n");
            tui_destruct(src);
            free(src_tui_path); free(out_path_owned);
            return 1;
        }
        st_logInfo("tui-coarsen: %" PRIi64 " genomes from .tui g-records\n", n_genomes);
    }

    /* Pull seq lengths for the directory.  Build name-sorted directory
     * (mirrors tui_create's tui.c:1441-1456 d-block). */
    stHash *seqlens = tui_sequence_lengths(src_tui_path);
    if (seqlens == NULL) {
        fprintf(stderr, "tui-coarsen: tui_sequence_lengths failed\n");
        tui_genome_info_free(roster, n_genomes);
        tui_destruct(src);
        free(src_tui_path); free(out_path_owned);
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

    /* seq_name -> dir_idx, for parent_S_ord stamping in C records. */
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
             "LOD B=%" PRIi64 " from %s", bin_size, src_tui_path);
    OneSchema *schema = NULL;
    OneFile *of = tui_open_write(out_path, "taffy", "tui-coarsen", blurb, &schema);
    if (of == NULL) {
        stHash_destruct(seq_to_dir_idx);
        free(dir);
        stHash_destruct(seqlens);
        tui_genome_info_free(roster, n_genomes);
        tui_destruct(src);
        free(src_tui_path); free(out_path_owned);
        tui_atexit_disarm();
        return 1;
    }

    /* t: total columns (coarse units = ceil(T_src / B)). */
    int64_t T_coarse = (T_src + bin_size - 1) / bin_size;
    oneInt(of, 0) = T_coarse;
    oneWriteLine(of, 't', 0, NULL);

    /* X: stub single (0,0) anchor.  LOD has no underlying MAF to anchor;
     * structural compliance only.  taffy view -r on a LOD .tui would
     * fail loudly, which is what we want. */
    {
        int64_t c0 = 0, f0 = 0;
        int64_t x_raw = 0, x_def = 0;
        uint8_t *xdef = tui_encode_idx(&c0, &f0, 1, &x_raw, &x_def);
        oneInt(of, 0) = x_raw;
        oneInt(of, 1) = 1;
        oneWriteLine(of, 'X', x_def, xdef);
        free(xdef);
    }

    /* d-records are written LAST -- see the rationale block below the
     * per-genome loop.  Briefly: each d-record stores the s_ord of its
     * matching S-record, but the s_ord assignment depends on the order
     * S-records are emitted (which follows the genome-roster order, not
     * the name-sorted dir order).  So we collect the (name -> s_ord)
     * mapping while writing S-records, then sort + emit d-records after. */

    /* full_name -> int64_t* s_ord, populated by write_genome's S writes. */
    stHash *capture = stHash_construct3(
        stHash_stringKey, stHash_stringEqualKey, free, free);

    /* Per-genome coarsen + write.  Single-threaded for v1; OpenMP is a
     * follow-up (the parallel-for-ordered template at tui.c:1494-1512). */
    int64_t c_ord_emit = 0;
    int64_t s_ord_emit = 0;
    time_t per_genome_t0 = time(NULL);
    for (int64_t gi = 0; gi < n_genomes; gi++) {
        const char *gname = roster[gi].name;
        time_t g_t0 = time(NULL);

        TuiGenomeLift *gl = tui_genome_lift_load(src, gname);
        if (gl == NULL) {
            st_logInfo("tui-coarsen: skip %s (no chunks)\n", gname);
            continue;
        }

        CoarsenCtx cx = {
            .B = bin_size,
            .cells = stHash_construct3(cell_key_hash, cell_key_equal, free, free),
        };
        tui_genome_lift_stream_runs(gl, coarsen_visit_cb, &cx);

        int64_t n_cells_raw = stHash_size(cx.cells);
        if (max_paralogs_per_bin > 0) {
            int64_t dropped = cell_dedup_paralogs(cx.cells, max_paralogs_per_bin);
            if (dropped > 0) {
                st_logDebug("tui-coarsen: %s: dedup dropped %" PRIi64
                            " of %" PRIi64 " cells (max=%" PRIi64 ")\n",
                            gname, dropped, n_cells_raw, max_paralogs_per_bin);
            }
        }

        SeqGroup *groups = NULL;
        int64_t n_groups = drain_cells(cx.cells, min_cell_bp,
                                       gname, seq_to_dir_idx, &groups);
        stHash_destruct(cx.cells);
        tui_genome_lift_destruct(gl);

        if (n_groups == 0) {
            free(groups);
            st_logInfo("tui-coarsen: %s -> 0 surviving cells (all below floor)\n", gname);
            continue;
        }

        write_genome(of, groups, n_groups, bin_size, &c_ord_emit,
                     &s_ord_emit, capture,
                     fill_bins, min_bin_frac);

        /* Diagnostic: per-genome RSS + glibc heap trim.  On the 577-way
         * we saw RSS grow ~30 MB/genome despite every per-genome alloc
         * being freed; suspecting glibc ptmalloc arena retention from
         * the millions-of-tiny-allocs workload (CellKey + int64_t +
         * EmitCell per cell).  malloc_trim(0) hands free chunks at the
         * top of the arena back to the OS; if RSS plateaus after this
         * line, fragmentation was the cause.  Cheap (microseconds for
         * arenas of this shape), so keep it on. */
        malloc_trim(0);
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        st_logInfo("tui-coarsen: %s: %" PRIi64 " cells (raw) over %" PRIi64
                   " seqs in %" PRIi64 " s [rss %ld MB]\n",
                   gname, n_cells_raw, n_groups,
                   (int64_t)(time(NULL) - g_t0),
                   ru.ru_maxrss / 1024);
    }
    st_logInfo("tui-coarsen: per-genome work done in %" PRIi64 " s\n",
               (int64_t)(time(NULL) - per_genome_t0));

    /* d: directory.  Written AFTER S/C/R because each d-record stores
     * the s_ord of its matching S-record (read back via oneGoto('S', ord)),
     * and that ordinal is only known once write_genome has run.
     *
     * Dir order is name-sorted globally across all genomes (binary-search
     * lookup requirement; see tui.h:97 + tui_find_d_lower_bound).  S-records
     * are emitted in genome-roster order, so we cannot match s_ord = dir_idx.
     * Instead we look each dir entry up in `capture`, populated by the per-
     * genome S writes, and use its actual s_ord.  Inactive seqs (no
     * alignment, so no S-record was written) get stored ord = -1 so the
     * reader's oneGoto('S', 0) silently fails and the genome lift skips
     * that seq cleanly.
     *
     * ONElib is order-agnostic on object types: the '&' footer index records
     * file offsets per ordinal regardless of where in the file each object
     * was written, so d-records after S/C/R index correctly.  Verified
     * against the read-side binary search in tui_find_d_lower_bound. */
    int64_t active_seqs = 0;
    for (int64_t i = 0; i < n_seqs; i++) {
        int64_t *ord_p = (int64_t *)stHash_search(capture, dir[i].name);
        int64_t stored_ord = ord_p ? (*ord_p - 1) : (int64_t)(-1);
        oneInt(of, 1) = stored_ord;
        oneInt(of, 2) = dir[i].length;
        oneWriteLine(of, 'd', strlen(dir[i].name), (void *)dir[i].name);
        if (ord_p) active_seqs++;
    }
    st_logInfo("tui-coarsen: wrote %" PRIi64 " d-records (%" PRIi64
               " active, %" PRIi64 " inactive sentinel ord=-1)\n",
               (int64_t)n_seqs, active_seqs, (int64_t)n_seqs - active_seqs);

    /* g: genome roster -- copy over from source. */
    for (int64_t gi = 0; gi < n_genomes; gi++) {
        oneInt(of, 1) = roster[gi].total_bp;
        oneInt(of, 2) = roster[gi].n_chroms;
        oneWriteLine(of, 'g', strlen(roster[gi].name), (void *)roster[gi].name);
    }

    stHash_destruct(capture);

    oneFileClose(of);
    oneSchemaDestroy(schema);
    tui_atexit_disarm();

    if (lod_txt != NULL) {
        update_lod_txt(lod_txt, bin_size, out_path);
        st_logInfo("tui-coarsen: updated %s\n", lod_txt);
    }

    stHash_destruct(seq_to_dir_idx);
    free(dir);
    stHash_destruct(seqlens);
    tui_genome_info_free(roster, n_genomes);
    tui_destruct(src);
    free(src_tui_path);
    free(out_path_owned);

    st_logInfo("tui-coarsen: total wall %" PRIi64 " s; wrote %s\n",
               (int64_t)(time(NULL) - startTime), out_path);
    return 0;
}
