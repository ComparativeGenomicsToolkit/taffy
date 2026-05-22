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

    t0 = time(NULL);
    stHash *wig = wig_parse(wig_file, "", 1);    // zero-based internally
    if (wig == NULL) {
        fprintf(stderr, "ERROR: failed to parse wig %s\n", wig_file);
        tui_genome_lift_destruct(gl);
        tui_destruct(tui);
        free(tui_p);
        return 1;
    }
    st_logInfo("Parsed wig in %" PRIi64 " s\n", (int64_t)(time(NULL) - t0));

    // Build a string-table for the target sequence names produced by the
    // lift so the output records can be sorted by (seq_idx, pos) cheaply.
    stHash *seqtab = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    stList *seqs   = stList_construct3(0, NULL);
    LiftRec *out = NULL;
    int64_t out_n = 0, out_cap = 0;
    int64_t n_in = 0, n_lifted = 0, n_no_column = 0, n_no_genome = 0;

    t0 = time(NULL);
    stHashIterator *it = stHash_getIterator(wig);
    char *anc_seq;
    while ((anc_seq = stHash_getNext(it)) != NULL) {
        stHash *inner = stHash_search(wig, anc_seq);
        if (inner == NULL) continue;
        stHashIterator *jt = stHash_getIterator(inner);
        void *kv;
        while ((kv = stHash_getNext(jt)) != NULL) {
            int64_t coord = (int64_t)kv;
            double  value = *(double*)stHash_search(inner, kv);
            n_in++;

            int64_t n_iv = 0;
            TuiInterval *iv = tui_query(tui, anc_seq, coord, coord + 1, &n_iv);
            if (n_iv == 0 || iv == NULL) {
                free(iv);
                n_no_column++;
                continue;
            }
            int64_t column = iv[0].start;
            free(iv);

            const char *tseq_name = NULL;
            int64_t tpos = 0;
            bool tstrand = true;
            if (!tui_genome_lift_column(gl, column, &tseq_name, &tpos, &tstrand)) {
                n_no_genome++;
                continue;
            }

            // String-table the target sequence name (gl owns the string;
            // we just need a stable integer index for sort).
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
        stHash_destructIterator(jt);
    }
    stHash_destructIterator(it);
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
        stList_destruct(seqs);
        stHash_destruct(seqtab);
        stHash_destruct(wig);
        tui_genome_lift_destruct(gl);
        tui_destruct(tui);
        free(tui_p);
        return 1;
    }
    int64_t cur_seq = -1;
    for (int64_t i = 0; i < out_n; i++) {
        if (out[i].seq_idx != cur_seq) {
            cur_seq = out[i].seq_idx;
            fprintf(fh, "variableStep chrom=%s\n", (const char*)stList_get(seqs, cur_seq));
        }
        // wig is 1-based by spec; our internal coord is 0-based.
        fprintf(fh, "%" PRIi64 " %g\n", out[i].pos + 1, out[i].value);
    }
    if (output_file) fclose(fh);

    free(out);
    stList_destruct(seqs);
    stHash_destruct(seqtab);
    stHash_destruct(wig);
    tui_genome_lift_destruct(gl);
    tui_destruct(tui);
    free(tui_p);
    st_logInfo("Total wall: %" PRIi64 " s\n", (int64_t)(time(NULL) - startTime));
    return 0;
}
