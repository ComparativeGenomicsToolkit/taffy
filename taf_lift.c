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
#include <errno.h>
#include <signal.h>
#include <unistd.h>

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

/////////////////////////////////////////////////////////////////////////////
// External sort: when the in-memory buffer of LiftRec hits SPILL_BUDGET, we
// qsort it and write the sorted records to a fresh spill file, then reset
// the buffer.  At end-of-input, the residual is flushed the same way; the
// final wig output is produced by a k-way heap-merge of all spills.  If
// nothing ever spilled (small lift) we fall back to the in-memory emit.
//
// LiftRec.seq_idx is a process-global index into the `seqs` stList, so
// records remain consistently ordered across spills.
//
// SPILL_BUDGET picked so 60+ spills are unusual at vertebrate scale: 256 MB
// per spill (~11M records) gives ~6 spills for a 70M-record fish-of-tetra
// lift, ~250 spills for a 2.5B-record full-gene-set hg38 lift -- both well
// within "open file" limits and merge-heap costs.
// Override at runtime with TAFFY_LIFT_BUDGET_MB.
/////////////////////////////////////////////////////////////////////////////
#define LIFT_SPILL_BUDGET_DEFAULT_MB 256

static int64_t lift_spill_budget_bytes(void) {
    const char *env = getenv("TAFFY_LIFT_BUDGET_MB");
    int64_t mb = LIFT_SPILL_BUDGET_DEFAULT_MB;
    if (env != NULL && *env != 0) {
        char *end = NULL;
        long v = strtol(env, &end, 10);
        if (end != env && v > 0) mb = (int64_t)v;
    }
    return mb * 1024LL * 1024LL;
}

// Tracked spill paths so we can rm() on normal exit AND on
// atexit/signal aborts (mirrors the tui.c cleanup pattern).  All spills
// live inside a single per-process mkdtemp'd directory (lift_spill_dir)
// to avoid TOCTOU symlink attacks on /tmp and adjacent-to-output name
// collisions across concurrent invocations.
static char **lift_spill_paths   = NULL;
static int    lift_spill_n       = 0;
static int    lift_spill_cap     = 0;
static int    lift_atexit_armed  = 0;
static char  *lift_spill_dir     = NULL;  // owned; rmdir'd on cleanup

// Normal-path cleanup: unlink each spill file, rmdir the spill directory,
// free the path strings, clear the tracker.  Calls free() / remove() /
// rmdir() which are NOT async-signal-safe -- only invoke on normal exit.
static void lift_remove_spills(void) {
    for (int i = 0; i < lift_spill_n; i++) {
        if (lift_spill_paths[i] != NULL) {
            (void) remove(lift_spill_paths[i]);
            free(lift_spill_paths[i]);
            lift_spill_paths[i] = NULL;
        }
    }
    free(lift_spill_paths);
    lift_spill_paths = NULL;
    lift_spill_n = lift_spill_cap = 0;
    if (lift_spill_dir != NULL) {
        (void) rmdir(lift_spill_dir);
        free(lift_spill_dir);
        lift_spill_dir = NULL;
    }
}

static void lift_atexit_handler(void) {
    if (lift_atexit_armed) lift_remove_spills();
}

// Signal-path cleanup: ONLY async-signal-safe primitives.  unlink() and
// rmdir() are on POSIX's safe list, free()/remove() are NOT (free can
// deadlock if the signal arrived mid-malloc; remove() isn't async-safe).
// We don't free the path strings -- the _exit() below takes the process
// down before any leak matters.  _exit(128+sig) gives the standard
// "killed by signal" exit code and skips both atexit handlers and stdio
// flushing (also unsafe in a signal context).
static void lift_signal_handler(int sig) {
    for (int i = 0; i < lift_spill_n; i++) {
        if (lift_spill_paths[i] != NULL) {
            (void) unlink(lift_spill_paths[i]);
        }
    }
    if (lift_spill_dir != NULL) {
        (void) rmdir(lift_spill_dir);
    }
    _exit(128 + sig);
}

static void lift_arm_cleanup(void) {
    if (lift_atexit_armed) return;
    lift_atexit_armed = 1;
    atexit(lift_atexit_handler);
    signal(SIGINT,  lift_signal_handler);
    signal(SIGTERM, lift_signal_handler);
    signal(SIGHUP,  lift_signal_handler);
}

// Create the per-process spill directory on first call.  Uses mkdtemp so
// the directory name is unique (avoids races with other taffy instances
// AND can't be pre-created by a hostile actor on shared systems --
// mkdtemp creates with mode 0700 and follows no symlinks).  Returns 0 on
// success, nonzero on mkdtemp failure.
static int lift_ensure_spill_dir(const char *output_file) {
    if (lift_spill_dir != NULL) return 0;
    const char *tmpdir = getenv("TMPDIR");
    if (tmpdir == NULL) tmpdir = "/tmp";
    char *tmpl = st_malloc(strlen(output_file ? output_file : tmpdir) + 32);
    if (output_file) {
        sprintf(tmpl, "%s.spills.XXXXXX", output_file);
    } else {
        sprintf(tmpl, "%s/taffy_lift.XXXXXX", tmpdir);
    }
    if (mkdtemp(tmpl) == NULL) {
        fprintf(stderr, "ERROR: failed to create spill directory %s: %s\n",
                tmpl, strerror(errno));
        free(tmpl);
        return 1;
    }
    lift_spill_dir = tmpl;
    return 0;
}

// Compose <lift_spill_dir>/spill.NNNNN.  Caller owns the returned string.
static char *lift_spill_path(int idx) {
    char *p = st_malloc(strlen(lift_spill_dir) + 32);
    sprintf(p, "%s/spill.%05d", lift_spill_dir, idx);
    return p;
}

// qsort the in-memory buffer, write all records to a new spill file (binary
// LiftRec back to back), record its path in the cleanup list.  Returns 0
// on success, nonzero on write/open failure.
static int lift_flush_spill(LiftRec *buf, int64_t n, const char *output_file) {
    if (n <= 0) return 0;
    qsort(buf, n, sizeof(LiftRec), liftrec_cmp);
    lift_arm_cleanup();
    if (lift_ensure_spill_dir(output_file) != 0) return 1;

    if (lift_spill_n == lift_spill_cap) {
        lift_spill_cap = lift_spill_cap ? lift_spill_cap * 2 : 16;
        lift_spill_paths = st_realloc(lift_spill_paths, lift_spill_cap * sizeof(char *));
    }
    char *path = lift_spill_path(lift_spill_n);
    FILE *fh = fopen(path, "wb");
    if (fh == NULL) {
        fprintf(stderr, "ERROR: failed to open spill file %s: %s\n", path, strerror(errno));
        free(path);
        return 1;
    }
    size_t w = fwrite(buf, sizeof(LiftRec), (size_t)n, fh);
    int rc = fclose(fh);
    if (w != (size_t)n || rc != 0) {
        fprintf(stderr, "ERROR: short write or close failure on spill %s\n", path);
        (void) remove(path);
        free(path);
        return 1;
    }
    lift_spill_paths[lift_spill_n++] = path;
    return 0;
}

// Per-spill read cursor for the k-way merge.  Always exposes `head` as the
// next-to-emit record while not eof; spill_pop advances (refilling the
// read buffer as needed) and returns 0 on exhaustion.  The merge heap is
// keyed on `head`, never on buf indices, so it stays valid through
// advances without the heap key going stale.
#define SPILL_READ_BATCH 4096
typedef struct {
    FILE   *fh;
    LiftRec buf[SPILL_READ_BATCH];
    int     buf_n;
    int     buf_i;
    int     eof;
    LiftRec head;
} Spill;

static int spill_pop(Spill *s) {
    if (s->buf_i >= s->buf_n) {
        if (s->eof) return 0;
        size_t got = fread(s->buf, sizeof(LiftRec), SPILL_READ_BATCH, s->fh);
        s->buf_n = (int)got;
        s->buf_i = 0;
        if (got == 0) { s->eof = 1; return 0; }
    }
    s->head = s->buf[s->buf_i++];
    return 1;
}

static int heap_cmp(const Spill *a, const Spill *b) {
    return liftrec_cmp(&a->head, &b->head);
}

static void heap_sift_down(Spill **h, int n, int i) {
    while (1) {
        int l = 2*i + 1, r = 2*i + 2, best = i;
        if (l < n && heap_cmp(h[l], h[best]) < 0) best = l;
        if (r < n && heap_cmp(h[r], h[best]) < 0) best = r;
        if (best == i) break;
        Spill *tmp = h[i]; h[i] = h[best]; h[best] = tmp;
        i = best;
    }
}

static void heap_build(Spill **h, int n) {
    for (int i = n/2 - 1; i >= 0; i--) heap_sift_down(h, n, i);
}

// k-way merge: open each spill, pop the first record (priming `head`),
// build a min-heap, repeatedly emit head[0], pop next into the top entry
// and sift -- or evict if that spill is now exhausted.
static int lift_merge_spills(char **paths, int n_paths, FILE *out_fh, stList *seqs) {
    Spill *all = st_calloc(n_paths, sizeof(Spill));
    Spill **heap = st_malloc(n_paths * sizeof(Spill *));
    int heap_n = 0;
    for (int i = 0; i < n_paths; i++) {
        all[i].fh = fopen(paths[i], "rb");
        if (all[i].fh == NULL) {
            fprintf(stderr, "ERROR: failed to open spill %s for merge: %s\n",
                    paths[i], strerror(errno));
            for (int j = 0; j < i; j++) if (all[j].fh) fclose(all[j].fh);
            free(all); free(heap);
            return 1;
        }
        if (spill_pop(&all[i])) heap[heap_n++] = &all[i];
    }
    heap_build(heap, heap_n);

    int64_t cur_out_seq = -1;
    while (heap_n > 0) {
        Spill *top = heap[0];
        if (top->head.seq_idx != cur_out_seq) {
            cur_out_seq = top->head.seq_idx;
            fprintf(out_fh, "variableStep chrom=%s\n",
                    (const char *)stList_get(seqs, cur_out_seq));
        }
        fprintf(out_fh, "%" PRIi64 " %g\n", top->head.pos + 1, top->head.value);
        if (spill_pop(top)) {
            heap_sift_down(heap, heap_n, 0);
        } else {
            // Exhausted -- swap last entry in, shrink heap, re-sift.
            heap[0] = heap[--heap_n];
            if (heap_n > 0) heap_sift_down(heap, heap_n, 0);
        }
    }
    for (int i = 0; i < n_paths; i++) if (all[i].fh) fclose(all[i].fh);
    free(all); free(heap);
    return 0;
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
    st_logInfo("Indexed %" PRIi64 " chunks for target genome '%s' in %" PRIi64 " s\n",
               tui_genome_lift_n_chunks(gl), target_genome, (int64_t)(time(NULL) - t0));

    // String-table for the output target-sequence names so the output
    // records can be sorted by (seq_idx, pos) cheaply.
    stHash *seqtab = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    stList *seqs   = stList_construct3(0, NULL);
    LiftRec *out = NULL;
    int64_t out_n = 0, out_cap = 0;
    int64_t n_in = 0, n_lifted = 0, n_no_column = 0, n_no_genome = 0;
    const int64_t spill_budget = lift_spill_budget_bytes();
    const int64_t spill_budget_n = spill_budget / (int64_t)sizeof(LiftRec);

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
                    if      (strcmp(key, "chrom") == 0) {
                        // strdup immediately: val points into the line buffer,
                        // and the *p++ = saved restoration below would put the
                        // original whitespace back over val's NUL terminator
                        // -- after which val would read past its intended end
                        // (e.g. on "chrom=X start=1 step=1", val would become
                        // "X start=1 step=1" once we restore the separator).
                        free(chrom);
                        chrom = stString_copy(val);
                    }
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
            // `chrom` is already a heap copy (strdup'd from the line buffer
            // at parse time); move ownership directly into wig_chrom.
            free(wig_chrom);
            wig_chrom = chrom;
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

        // Universal column -> target genome bases.  One source position can
        // lift to MULTIPLE target positions (paralogs of the same target
        // genome aligning to the same ancestral column).  Stack buffer for
        // the common case (0 or 1 paralog); fall back to heap on overflow.
        TuiGenomeMatch stack_m[32];
        TuiGenomeMatch *m = stack_m;
        TuiGenomeMatch *heap_m = NULL;
        int nm = tui_genome_lift_column(gl, column, stack_m, 32);
        if (nm == 0) { n_no_genome++; continue; }
        if (nm > 32) {
            heap_m = st_malloc((size_t)nm * sizeof(TuiGenomeMatch));
            tui_genome_lift_column(gl, column, heap_m, nm);
            m = heap_m;
        }
        for (int k = 0; k < nm; k++) {
            const char *tseq_name = m[k].seq;
            int64_t     tpos      = m[k].pos;
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
                // Double until we hit the spill budget, then stop growing --
                // the post-append spill-flush below will reset out_n and we
                // never need cap > spill_budget_n.  (The previous version had
                // a "safety re-add" branch that could push cap above the
                // budget for one realloc on tiny budgets; not needed since
                // initial cap is 1024 <= any sensible budget.)
                int64_t new_cap = out_cap ? out_cap * 2 : 1024;
                if (new_cap > spill_budget_n) new_cap = spill_budget_n;
                out_cap = new_cap;
                out = st_realloc(out, out_cap * sizeof(LiftRec));
            }
            out[out_n].seq_idx = seq_idx;
            out[out_n].pos     = tpos;
            out[out_n].value   = value;
            out_n++;
            n_lifted++;
            // External-sort spill: cap memory at spill_budget_n records.
            if (out_n >= spill_budget_n) {
                if (lift_flush_spill(out, out_n, output_file) != 0) {
                    free(heap_m);
                    free(out); free(cur_seq); free(cur_runs); free(wig_chrom);
                    LI_destruct(wli); fclose(wf);
                    stList_destruct(seqs); stHash_destruct(seqtab);
                    tui_genome_lift_destruct(gl); tui_destruct(tui); free(tui_p);
                    return 1;
                }
                out_n = 0;
            }
        }
        free(heap_m);
    }
    LI_destruct(wli);
    fclose(wf);
    free(wig_chrom);

    st_logInfo("Lifted %" PRIi64 "/%" PRIi64 " records in %" PRIi64 " s "
               "(no-column=%" PRIi64 ", no-genome=%" PRIi64 ", spills=%d)\n",
               n_lifted, n_in, (int64_t)(time(NULL) - t0),
               n_no_column, n_no_genome, lift_spill_n);

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

    if (lift_spill_n == 0) {
        // Nothing ever spilled -- emit directly from the in-memory buffer.
        qsort(out, out_n, sizeof(LiftRec), liftrec_cmp);
        int64_t cur_out_seq = -1;
        for (int64_t i = 0; i < out_n; i++) {
            if (out[i].seq_idx != cur_out_seq) {
                cur_out_seq = out[i].seq_idx;
                fprintf(fh, "variableStep chrom=%s\n",
                        (const char*)stList_get(seqs, cur_out_seq));
            }
            fprintf(fh, "%" PRIi64 " %g\n", out[i].pos + 1, out[i].value);
        }
    } else {
        // Flush the residual to a final spill, then k-way merge all spills.
        if (out_n > 0 && lift_flush_spill(out, out_n, output_file) != 0) {
            if (output_file) fclose(fh);
            free(out); free(cur_seq); free(cur_runs);
            stList_destruct(seqs); stHash_destruct(seqtab);
            tui_genome_lift_destruct(gl); tui_destruct(tui); free(tui_p);
            return 1;
        }
        out_n = 0;
        time_t t_merge = time(NULL);
        int rc = lift_merge_spills(lift_spill_paths, lift_spill_n, fh, seqs);
        st_logInfo("Merged %d spill(s) in %" PRIi64 " s\n",
                   lift_spill_n, (int64_t)(time(NULL) - t_merge));
        if (rc != 0) {
            if (output_file) fclose(fh);
            free(out); free(cur_seq); free(cur_runs);
            stList_destruct(seqs); stHash_destruct(seqtab);
            tui_genome_lift_destruct(gl); tui_destruct(tui); free(tui_p);
            lift_remove_spills();
            return 1;
        }
        lift_remove_spills();
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
