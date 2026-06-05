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
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
    stHash *by_seq = stHash_construct();      /* seq -> stList of EmitCell* */
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
 * TUI_CHUNK_G_MAX caps so the reader's per-chunk t-range skip works. */
static void write_genome(OneFile *of, SeqGroup *groups, int64_t n_groups,
                         int64_t B, int64_t *c_ord_emit) {
    for (int64_t gi = 0; gi < n_groups; gi++) {
        SeqGroup *g = &groups[gi];
        /* S line: seq name + length-in-coarse-units = number of g_bins.
         * Reader treats this as the "seq length" field.  Coarse units. */
        int64_t max_g_bin = 0;
        for (int64_t i = 0; i < g->n_cells; i++)
            if (g->cells[i].g_bin > max_g_bin) max_g_bin = g->cells[i].g_bin;
        oneInt(of, 1) = (max_g_bin + 1) * B;   /* approximate seq_len upper bound */
        oneWriteLine(of, 'S', strlen(g->full_name), (void *)g->full_name);
        int64_t s_ord = g->seq_idx_in_dir + 1;
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
                buf[3*k + 0] = e->g_bin * B;     /* t_start */
                buf[3*k + 1] = e->uc_bin * B;    /* g_start */
                buf[3*k + 2] = (e->bp << 1) | 0; /* lenc: + strand */
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
    char    *log_level      = NULL;

    enum { OPT_MIN_CELL_BP = 256, OPT_FROM_TUI, OPT_LOD_TXT };
    static struct option long_options[] = {
        { "inputFile",  required_argument, 0, 'i' },
        { "outputFile", required_argument, 0, 'o' },
        { "binSize",    required_argument, 0, 'B' },
        { "minCellBp",  required_argument, 0, OPT_MIN_CELL_BP },
        { "fromTui",    no_argument,       0, OPT_FROM_TUI },
        { "lodTxt",     required_argument, 0, OPT_LOD_TXT },
        { "logLevel",   required_argument, 0, 'l' },
        { "help",       no_argument,       0, 'h' },
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
            case OPT_MIN_CELL_BP: min_cell_bp = atoll(optarg); break;
            case OPT_FROM_TUI:    from_tui = 1; break;
            case OPT_LOD_TXT:     lod_txt = optarg; break;
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

    /* Load genome roster.  Hard-error on legacy .tui without g records;
     * coarsening needs the correct genome list. */
    int64_t n_genomes = 0;
    TuiGenomeInfo *roster = tui_genome_names(src_tui_path, &n_genomes);
    if (roster == NULL || n_genomes == 0) {
        fprintf(stderr,
                "tui-coarsen: source .tui has no genome roster (g records).\n"
                "  Rebuild the .tui with a newer taffy (commit b8f6811 or later)\n"
                "  so genome resolution is exact.  Refusing to coarsen.\n");
        tui_destruct(src);
        free(src_tui_path); free(out_path_owned);
        return 1;
    }
    st_logInfo("tui-coarsen: %" PRIi64 " genomes\n", n_genomes);

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

    /* d: directory, name-sorted, seqLen field = original base-pair length
     * (we keep this faithful so taffy stats etc. still report sensible
     * sizes; the COARSE seq-length lives implicitly in chunk g_max). */
    for (int64_t i = 0; i < n_seqs; i++) {
        oneInt(of, 1) = i;                   /* S-ord */
        oneInt(of, 2) = dir[i].length;       /* seq length in base pairs */
        oneWriteLine(of, 'd', strlen(dir[i].name), (void *)dir[i].name);
    }

    /* Per-genome coarsen + write.  Single-threaded for v1; OpenMP is a
     * follow-up (the parallel-for-ordered template at tui.c:1494-1512). */
    int64_t c_ord_emit = 0;
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

        write_genome(of, groups, n_groups, bin_size, &c_ord_emit);

        st_logInfo("tui-coarsen: %s: %" PRIi64 " cells (raw) over %" PRIi64
                   " seqs in %" PRIi64 " s\n",
                   gname, n_cells_raw, n_groups,
                   (int64_t)(time(NULL) - g_t0));
    }
    st_logInfo("tui-coarsen: per-genome work done in %" PRIi64 " s\n",
               (int64_t)(time(NULL) - per_genome_t0));

    /* g: genome roster -- copy over from source. */
    for (int64_t gi = 0; gi < n_genomes; gi++) {
        oneInt(of, 1) = roster[gi].total_bp;
        oneInt(of, 2) = roster[gi].n_chroms;
        oneWriteLine(of, 'g', strlen(roster[gi].name), (void *)roster[gi].name);
    }

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
