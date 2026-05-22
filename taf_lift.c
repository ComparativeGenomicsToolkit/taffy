/*
 * taffy lift: lift a .wig annotation from the universal MAF's ancestor
 * (row-0) coordinates into a target leaf genome's coordinates, using the
 * .tui universal-column index.  Per-record sparse lift (one wig record at
 * a time): each input record is translated via tui_query (anc -> universal
 * column) + tui_genome_lift_column (column -> target leaf coord).
 *
 *  Released under the MIT license, see LICENSE.txt
 */

#include "taf.h"
#include "tui.h"
#include "sonLib.h"
#include <getopt.h>
#include <time.h>

static void usage(void) {
    fprintf(stderr, "taffy lift [options]\n");
    fprintf(stderr, "Lift a .wig annotation from universal MAF row-0 (ancestor) coords to a target leaf genome via the .tui index\n");
    fprintf(stderr, "-i --inputFile [FILE_NAME] : REQUIRED Path to the universal MAF/TAF (its .tui sidecar is expected at <input>.tui)\n");
    fprintf(stderr, "-w --wig       [FILE_NAME] : REQUIRED Input .wig in ancestor row-0 coords (chrom = full genome.sequence)\n");
    fprintf(stderr, "-g --genome    [STRING]    : REQUIRED Target genome name (e.g. hg38)\n");
    fprintf(stderr, "-o --outputFile [FILE_NAME] : Output wig (default stdout)\n");
    fprintf(stderr, "-l --logLevel  : Set log level\n");
    fprintf(stderr, "-h --help      : Print this help\n");
}

typedef struct {
    int64_t seq_idx;
    int64_t pos;
    double  value;
} LiftRec;

static int liftrec_cmp(const void *a, const void *b) {
    const LiftRec *x = a, *y = b;
    if (x->seq_idx != y->seq_idx) return x->seq_idx < y->seq_idx ? -1 : 1;
    return x->pos < y->pos ? -1 : x->pos > y->pos ? 1 : 0;
}

int taf_lift_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);
    char *logLevelString = NULL;
    char *maf_file       = NULL;
    char *wig_file       = NULL;
    char *target_genome  = NULL;
    char *output_file    = NULL;

    while (1) {
        static struct option long_options[] = {
            { "logLevel",   required_argument, 0, 'l' },
            { "inputFile",  required_argument, 0, 'i' },
            { "wig",        required_argument, 0, 'w' },
            { "genome",     required_argument, 0, 'g' },
            { "outputFile", required_argument, 0, 'o' },
            { "help",       no_argument,       0, 'h' },
            { 0, 0, 0, 0 }
        };
        int idx = 0;
        int key = getopt_long(argc, argv, "l:i:w:g:o:h", long_options, &idx);
        if (key == -1) break;
        switch (key) {
            case 'l': logLevelString = optarg; break;
            case 'i': maf_file       = optarg; break;
            case 'w': wig_file       = optarg; break;
            case 'g': target_genome  = optarg; break;
            case 'o': output_file    = optarg; break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }
    st_setLogLevelFromString(logLevelString);

    if (maf_file == NULL || wig_file == NULL || target_genome == NULL) {
        fprintf(stderr, "ERROR: -i, -w, and -g are required\n");
        usage();
        return 1;
    }

    // The .tui lives at <maf>.tui (same convention as tai_path).
    char *tui_p = tui_path(maf_file);
    Tui *tui = tui_load(tui_p);
    if (tui == NULL) {
        fprintf(stderr, "ERROR: failed to open .tui at %s\n", tui_p);
        free(tui_p);
        return 1;
    }
    st_logInfo("Loaded .tui (%" PRIi64 " columns) in %" PRIi64 " s\n",
               tui_total_columns(tui), (int64_t)(time(NULL) - startTime));

    time_t t0 = time(NULL);
    TuiGenomeLift *gl = tui_genome_lift_load(tui, target_genome);
    if (gl == NULL) {
        fprintf(stderr, "ERROR: genome '%s' not found in .tui\n", target_genome);
        tui_destruct(tui);
        free(tui_p);
        return 1;
    }
    st_logInfo("Loaded %" PRIi64 " runs for target genome '%s' in %" PRIi64 " s\n",
               tui_genome_lift_n_runs(gl), target_genome, (int64_t)(time(NULL) - t0));

    // String-table for the output target-sequence names so the output
    // records can be sorted by (seq_idx, pos) cheaply.
    stHash *seqtab = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    stList *seqs   = stList_construct3(0, NULL);
    LiftRec *out = NULL;
    int64_t out_n = 0, out_cap = 0;
    int64_t n_in = 0, n_lifted = 0, n_no_column = 0, n_no_genome = 0;

    // Stream-parse the wig directly instead of going through wig_parse.
    // wig_parse builds an stHash<seq, stHash<coord -> double*>> which (a) costs
    // ~150 MB per million records and (b) uses (void*)coord as the inner-hash
    // key -- when coord == 0 the key is NULL, indistinguishable from the
    // stHash_getNext end-of-iteration sentinel, so any wig containing
    // ancestor position 0 silently iterates to zero records.  Streaming
    // avoids both issues and is closer to constant memory.
    FILE *wf = fopen(wig_file, "r");
    if (wf == NULL) {
        fprintf(stderr, "ERROR: failed to open wig %s\n", wig_file);
        tui_genome_lift_destruct(gl);
        tui_destruct(tui);
        free(tui_p);
        return 1;
    }
    LI *wli = LI_construct(wf);

    // Wig parser state.
    char    *wig_chrom = NULL;     // borrowed-then-owned chrom from latest header
    int      wig_fixed = 0;        // 0 = variableStep, 1 = fixedStep
    int64_t  wig_fs_pos = 0;       // next 0-based coord for fixedStep
    int64_t  wig_fs_step = 1;

    // Single-slot cache of the most-recently-loaded source-seq runs.  Refreshed
    // on each wig header line.  Wigs almost always group records by chrom so
    // this hits 100% after the first record of each chrom.
    char    *cur_seq  = NULL;
    int64_t *cur_runs = NULL;      // 3 * cur_n int64s: (t_start, g_start, lenc)
    int64_t  cur_n    = 0;
    int      cur_seq_in_tui = 0;

    t0 = time(NULL);
    char *line;
    while ((line = LI_get_next_line(wli)) != NULL) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0' || *p == '#') { free(line); continue; }

        int is_fixed = (strncmp(p, "fixedStep",   9) == 0);
        int is_var   = (strncmp(p, "variableStep", 12) == 0);
        if (is_fixed || is_var) {
            // Header line: tokenize the rest as key=value pairs.
            char *chrom = NULL;
            int64_t fs_start = 1, fs_step = 1;
            while (*p && *p != ' ' && *p != '\t') p++;  // skip keyword
            while (1) {
                while (*p == ' ' || *p == '\t') p++;
                if (!*p) break;
                char *tok = p;
                while (*p && *p != ' ' && *p != '\t') p++;
                char saved = *p; *p = '\0';
                char *eq = strchr(tok, '=');
                if (eq) {
                    *eq = '\0';
                    char *key = tok, *val = eq + 1;
                    if      (strcmp(key, "chrom") == 0) chrom = val;
                    else if (strcmp(key, "start") == 0) fs_start = atoll(val);
                    else if (strcmp(key, "step")  == 0) fs_step  = atoll(val);
                    // span is intentionally ignored for now; the writer that
                    // we care about (GERP & friends) emits step=1 + span=1.
                }
                if (saved == 0) break;
                *p++ = saved;
            }
            wig_fixed   = is_fixed;
            wig_fs_pos  = fs_start - 1;   // 1-based -> 0-based
            wig_fs_step = fs_step;
            free(wig_chrom);
            wig_chrom = chrom ? stString_copy(chrom) : NULL;
            // Refresh the cached source-seq runs on chrom change.
            if (wig_chrom != NULL &&
                (cur_seq == NULL || strcmp(wig_chrom, cur_seq) != 0)) {
                free(cur_seq); free(cur_runs);
                cur_seq = stString_copy(wig_chrom);
                cur_runs = tui_load_seq_runs(tui, wig_chrom, &cur_n);
                cur_seq_in_tui = (cur_runs != NULL);
            }
            free(line);
            continue;
        }

        // Data line.
        if (wig_chrom == NULL) { free(line); continue; }
        int64_t coord;
        double  value;
        if (wig_fixed) {
            value = atof(p);
            coord = wig_fs_pos;
            wig_fs_pos += wig_fs_step;
        } else {
            char *q = p;
            while (*q && *q != ' ' && *q != '\t') q++;
            if (!*q) { free(line); continue; }
            *q = '\0';
            coord = atoll(p) - 1;          // wig is 1-based
            value = atof(q + 1);
        }
        n_in++;
        free(line);

        if (!cur_seq_in_tui) { n_no_column++; continue; }

        // Binary-search cur_runs (sorted by t_start) for the run covering
        // coord, then map to a universal column.
        int64_t lo = 0, hi = cur_n;
        while (lo < hi) {
            int64_t mid = (lo + hi) / 2;
            if (cur_runs[3*mid + 0] <= coord) lo = mid + 1; else hi = mid;
        }
        int64_t i = lo - 1;
        if (i < 0) { n_no_column++; continue; }
        int64_t t    = cur_runs[3*i + 0];
        int64_t g    = cur_runs[3*i + 1];
        int64_t lenc = cur_runs[3*i + 2];
        int64_t len  = lenc >> 1;
        int     rev  = (int)(lenc & 1);
        if (coord >= t + len) { n_no_column++; continue; }
        int64_t column = rev ? g + (t + len - 1 - coord)
                             : g + (coord - t);

        const char *tseq_name = NULL;
        int64_t tpos = 0;
        bool tstrand = true;
        if (!tui_genome_lift_column(gl, column, &tseq_name, &tpos, &tstrand)) {
            n_no_genome++;
            continue;
        }
        void *v = stHash_search(seqtab, (void *)tseq_name);
        int64_t seq_idx;
        if (v == NULL) {
            seq_idx = stList_length(seqs);
            stList_append(seqs, (void *)tseq_name);
            stHash_insert(seqtab, (void *)tseq_name, (void *)(intptr_t)(seq_idx + 1));
        } else {
            seq_idx = (intptr_t)v - 1;
        }
        if (out_n == out_cap) {
            out_cap = out_cap ? out_cap * 2 : 1024;
            out = st_realloc(out, out_cap * sizeof(LiftRec));
        }
        out[out_n].seq_idx = seq_idx;
        out[out_n].pos     = tpos;
        out[out_n].value   = value;
        out_n++;
        n_lifted++;
    }
    LI_destruct(wli);
    fclose(wf);
    free(wig_chrom);

    st_logInfo("Lifted %" PRIi64 "/%" PRIi64 " records in %" PRIi64 " s "
               "(no-column=%" PRIi64 ", no-genome=%" PRIi64 ")\n",
               n_lifted, n_in, (int64_t)(time(NULL) - t0),
               n_no_column, n_no_genome);

    // Sort by (seq_idx, pos) so we can emit variableStep blocks.
    qsort(out, out_n, sizeof(LiftRec), liftrec_cmp);

    FILE *fh = output_file ? fopen(output_file, "w") : stdout;
    if (fh == NULL) {
        fprintf(stderr, "ERROR: failed to open output %s\n", output_file);
        free(out);
        free(cur_seq); free(cur_runs);
        stList_destruct(seqs);
        stHash_destruct(seqtab);
        tui_genome_lift_destruct(gl);
        tui_destruct(tui);
        free(tui_p);
        return 1;
    }
    int64_t cur_out_seq = -1;
    for (int64_t i = 0; i < out_n; i++) {
        if (out[i].seq_idx != cur_out_seq) {
            cur_out_seq = out[i].seq_idx;
            fprintf(fh, "variableStep chrom=%s\n", (const char*)stList_get(seqs, cur_out_seq));
        }
        // wig is 1-based by spec; our internal coord is 0-based.
        fprintf(fh, "%" PRIi64 " %g\n", out[i].pos + 1, out[i].value);
    }
    if (output_file) fclose(fh);

    free(out);
    free(cur_seq);
    free(cur_runs);
    stList_destruct(seqs);
    stHash_destruct(seqtab);
    tui_genome_lift_destruct(gl);
    tui_destruct(tui);
    free(tui_p);
    st_logInfo("Total wall: %" PRIi64 " s\n", (int64_t)(time(NULL) - startTime));
    return 0;
}
