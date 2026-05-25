/*
 * taffy lift: lift a .wig annotation from the universal MAF's ancestor
 * (row-0) coordinates into a target leaf genome's coordinates, using the
 * .tui universal-column index.  Per-record sparse lift (one wig record at
 * a time): each input record is translated via tui_query (anc -> universal
 * column) + tui_genome_lift_column (column -> target leaf coord).
 *
 *  Released under the MIT license, see LICENSE.txt
 */

/* mkdtemp() is in stdlib.h under _POSIX_C_SOURCE>=200809L; set the macro
 * before any system headers are pulled in by the project headers below. */
#define _POSIX_C_SOURCE 200809L

#include "taf.h"
#include "tui.h"
#include "sonLib.h"
#include <getopt.h>
#include <time.h>
#include <errno.h>
#include <signal.h>
#include <unistd.h>
#include <stdlib.h>

static void usage(void) {
    fprintf(stderr, "taffy lift [options]\n");
    fprintf(stderr, "Lift a .wig annotation from universal MAF row-0 (ancestor) coords to a target leaf genome via the .tui index\n");
    fprintf(stderr, "-i --inputFile [FILE_NAME] : REQUIRED Path to the universal MAF/TAF (its .tui sidecar is expected at <input>.tui)\n");
    fprintf(stderr, "-w --wig       [FILE_NAME] : REQUIRED Input .wig in ancestor row-0 coords (chrom = full genome.sequence)\n");
    fprintf(stderr, "-g --genome    [STRING]    : REQUIRED Target genome name (e.g. hg38)\n");
    fprintf(stderr, "-o --outputFile [FILE_NAME] : Output wig (default stdout)\n");
    fprintf(stderr, "-m --memCap    [SIZE]      : Max in-memory output buffer before spilling to disk; suffix K/M/G allowed (default 2G).  Env TAFFY_LIFT_BUDGET_MB also honored (CLI wins).\n");
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

// LSD radix sort LiftRec[] by composite key (seq_idx << 32) | pos.
// 8 bits per pass, 4-5 passes suffice at vertebrate scale (seq_idx in
// the low hundreds = 8 bits, pos in the low billions = ~28-32 bits).
// O(N) total -- on a 60M-record sort, this is ~3 s vs qsort's ~10 s
// (no function-pointer comparator overhead, sequential memory access).
// Cost: one auxiliary LiftRec[] (= memory doubles during the sort);
// this is paid once at end-of-input before the emit phase, both arrays
// freed before the wig is written.  Assumes both seq_idx and pos fit
// in uint32 -- safe for any realistic universal MAF (genome with >4B
// sequences in a single target, or a chromosome >4 Gb, would need a
// wider key).
//
// Wrapper liftrec_sort() dispatches by N: small buffers (under
// RADIX_THRESHOLD) use qsort directly because the radix sort's
// per-call malloc of a >MB tmp buffer hits glibc's mmap path on
// every call and dominates at small N -- the in-memory fast path
// sorts ~60M records once and radix wins big; the spill flush path
// sorts ~budget/24 records per call (default 2G/24 = ~83M, still
// radix territory; but at -m 1M the per-call buffer is ~44K records
// and qsort is faster).
#define RADIX_THRESHOLD 1000000

// Output emit format break-even.  variableStep (~14 B/row: "pos val\n")
// vs fixedStep (~50 B header + ~2 B/row: "val\n").  Equal at ~6 records.
// Used by both the in-memory emit path and the spill-merge lookahead.
#define FS_THRESHOLD 6

static void liftrec_radix_sort(LiftRec *arr, int64_t n) {
    if (n < 2) return;
    uint64_t max_key = 0;
    for (int64_t i = 0; i < n; i++) {
        // Guard the composite-key truncation assumption (see comment on
        // liftrec_radix_sort): if either field exceeds 32 bits the
        // composite key overflows and the sort silently misorders.
        // Catch it loudly here rather than corrupt output downstream.
        assert(arr[i].seq_idx >= 0 && arr[i].seq_idx < (1LL << 32));
        assert(arr[i].pos     >= 0 && arr[i].pos     < (1LL << 32));
        uint64_t k = ((uint64_t)arr[i].seq_idx << 32) | (uint32_t)arr[i].pos;
        if (k > max_key) max_key = k;
    }
    int n_passes = 0;
    while (n_passes < 8 && (max_key >> (n_passes * 8)) != 0) n_passes++;
    if (n_passes == 0) return;

    LiftRec *tmp = st_malloc((size_t)n * sizeof(LiftRec));
    LiftRec *src = arr, *dst = tmp;
    int64_t count[257];
    for (int pass = 0; pass < n_passes; pass++) {
        memset(count, 0, sizeof(count));
        int shift = pass * 8;
        for (int64_t i = 0; i < n; i++) {
            uint64_t k = ((uint64_t)src[i].seq_idx << 32) | (uint32_t)src[i].pos;
            count[((k >> shift) & 0xff) + 1]++;
        }
        for (int b = 0; b < 256; b++) count[b + 1] += count[b];
        for (int64_t i = 0; i < n; i++) {
            uint64_t k = ((uint64_t)src[i].seq_idx << 32) | (uint32_t)src[i].pos;
            uint8_t b = (k >> shift) & 0xff;
            dst[count[b]++] = src[i];
        }
        LiftRec *swap = src; src = dst; dst = swap;
    }
    if (src != arr) memcpy(arr, src, (size_t)n * sizeof(LiftRec));
    free(tmp);
}

// Pick the right algorithm for the input size.  Radix has fixed per-call
// overhead (a >MB malloc on the in-memory path's hot loop, glibc mmap'd
// at sizes above 128 KB) that's amortized at large N but loses to qsort
// at small N -- the spill-flush path with -m 1M lives there.
static void liftrec_sort(LiftRec *arr, int64_t n) {
    if (n < RADIX_THRESHOLD) qsort(arr, n, sizeof(LiftRec), liftrec_cmp);
    else                     liftrec_radix_sort(arr, n);
}

/////////////////////////////////////////////////////////////////////////////
// External sort: when the in-memory buffer of LiftRec hits the budget, we
// qsort it and write the sorted records to a fresh spill file, then reset
// the buffer.  At end-of-input, the residual is flushed the same way; the
// final wig output is produced by a k-way heap-merge of all spills.  If
// nothing ever spilled (small lift) we fall back to the in-memory emit.
//
// LiftRec.seq_idx is a process-global index into the `seqs` stList, so
// records remain consistently ordered across spills.
//
// Budget precedence: --memCap CLI > TAFFY_LIFT_BUDGET_MB env > default.
// Default of 4 GB handles most single-chrom and full-gene-set lifts in
// memory; the spill path only kicks in for whole-genome dense lifts (e.g.
// chr1 step=1 hg38 -> ape ~5 GB, full canonical gene set ~30 GB).
/////////////////////////////////////////////////////////////////////////////
#define LIFT_SPILL_BUDGET_DEFAULT_BYTES (2LL * 1024 * 1024 * 1024)

// Parse "12345", "256M", "4G", "1024K" -> bytes; -1 on parse error.  Empty
// or whitespace-only input is also an error (callers treat -1 as "absent
// CLI override" and fall back to env/default; if --memCap was passed
// explicitly with garbage we still error out before we get here).
static int64_t lift_parse_size(const char *s) {
    if (s == NULL || *s == 0) return -1;
    char *end = NULL;
    long long v = strtoll(s, &end, 10);
    if (end == s || v < 0) return -1;
    int64_t mult = 1;
    if (*end != 0) {
        if (*(end + 1) != 0) return -1;
        switch (*end) {
            case 'k': case 'K': mult = 1024LL; break;
            case 'm': case 'M': mult = 1024LL * 1024; break;
            case 'g': case 'G': mult = 1024LL * 1024 * 1024; break;
            default: return -1;
        }
    }
    return (int64_t)v * mult;
}

// Resolve the budget: CLI override (parsed at argv-time, passed in here as
// >0 or -1 for "not set") wins over the env var, wins over the default.
static int64_t lift_spill_budget_bytes(int64_t cli_override_bytes) {
    if (cli_override_bytes > 0) return cli_override_bytes;
    const char *env = getenv("TAFFY_LIFT_BUDGET_MB");
    if (env != NULL && *env != 0) {
        char *end = NULL;
        long v = strtol(env, &end, 10);
        if (end != env && v > 0) return (int64_t)v * 1024LL * 1024LL;
    }
    return LIFT_SPILL_BUDGET_DEFAULT_BYTES;
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
    liftrec_sort(buf, n);
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
// Pop the heap's smallest record into *out, advance the source spill,
// rebalance the heap.  Returns 1 if a record was popped, 0 if heap empty.
static int heap_pop_into(Spill **heap, int *heap_n, LiftRec *out) {
    if (*heap_n == 0) return 0;
    *out = heap[0]->head;
    if (spill_pop(heap[0])) {
        heap_sift_down(heap, *heap_n, 0);
    } else {
        heap[0] = heap[--(*heap_n)];
        if (*heap_n > 0) heap_sift_down(heap, *heap_n, 0);
    }
    return 1;
}

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

    // Streaming-collapse emit: same logic as the in-memory fast path,
    // but with a small ring buffer for the lookahead needed to decide
    // fixedStep vs variableStep.  Refill from the heap as the buffer
    // drains.
    LiftRec buf[FS_THRESHOLD];
    int buf_n = 0;
    int64_t cur_var_seq = -1;

    // Initial fill.
    while (buf_n < FS_THRESHOLD && heap_n > 0) {
        if (heap_pop_into(heap, &heap_n, &buf[buf_n])) buf_n++;
    }

    while (buf_n > 0) {
        int rl = 1;
        while (rl < buf_n &&
               buf[rl].seq_idx == buf[0].seq_idx &&
               buf[rl].pos == buf[rl - 1].pos + 1) {
            rl++;
        }
        if (rl >= FS_THRESHOLD) {
            int64_t cur_seq = buf[0].seq_idx;
            int64_t next_pos = buf[0].pos + 1;
            fprintf(out_fh, "fixedStep chrom=%s start=%" PRIi64 " step=1\n",
                    (const char *)stList_get(seqs, cur_seq), buf[0].pos + 1);
            for (int k = 0; k < rl; k++) {
                fprintf(out_fh, "%g\n", buf[k].value);
                next_pos = buf[k].pos + 1;
            }
            // Slide remaining records down.
            for (int k = rl; k < buf_n; k++) buf[k - rl] = buf[k];
            buf_n -= rl;
            // Keep extending the run by pulling from the heap as long
            // as the next-smallest stays consecutive with this fixedStep.
            // We must check buf[0] (refilled if empty) first since the
            // heap top alone isn't enough -- buf may still hold a record
            // that's earlier in the merge order than the new heap top.
            while (1) {
                if (buf_n == 0 && heap_n > 0) {
                    if (heap_pop_into(heap, &heap_n, &buf[0])) buf_n = 1;
                }
                if (buf_n == 0) break;
                if (buf[0].seq_idx != cur_seq || buf[0].pos != next_pos) break;
                fprintf(out_fh, "%g\n", buf[0].value);
                next_pos++;
                // Slide one off the front.
                for (int k = 1; k < buf_n; k++) buf[k - 1] = buf[k];
                buf_n--;
            }
            cur_var_seq = -1;
        } else {
            // Emit buf[0] as variableStep.
            if (buf[0].seq_idx != cur_var_seq) {
                fprintf(out_fh, "variableStep chrom=%s\n",
                        (const char *)stList_get(seqs, buf[0].seq_idx));
                cur_var_seq = buf[0].seq_idx;
            }
            fprintf(out_fh, "%" PRIi64 " %g\n", buf[0].pos + 1, buf[0].value);
            for (int k = 1; k < buf_n; k++) buf[k - 1] = buf[k];
            buf_n--;
        }
        // Refill buf to FS_THRESHOLD if more records remain.
        while (buf_n < FS_THRESHOLD && heap_n > 0) {
            if (heap_pop_into(heap, &heap_n, &buf[buf_n])) buf_n++;
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
    int64_t mem_cap_cli  = -1;     // -1 sentinel = "not set on CLI"

    while (1) {
        static struct option long_options[] = {
            { "logLevel",   required_argument, 0, 'l' },
            { "inputFile",  required_argument, 0, 'i' },
            { "wig",        required_argument, 0, 'w' },
            { "genome",     required_argument, 0, 'g' },
            { "outputFile", required_argument, 0, 'o' },
            { "memCap",     required_argument, 0, 'm' },
            { "help",       no_argument,       0, 'h' },
            { 0, 0, 0, 0 }
        };
        int idx = 0;
        int key = getopt_long(argc, argv, "l:i:w:g:o:m:h", long_options, &idx);
        if (key == -1) break;
        switch (key) {
            case 'l': logLevelString = optarg; break;
            case 'i': maf_file       = optarg; break;
            case 'w': wig_file       = optarg; break;
            case 'g': target_genome  = optarg; break;
            case 'o': output_file    = optarg; break;
            case 'm':
                mem_cap_cli = lift_parse_size(optarg);
                if (mem_cap_cli <= 0) {
                    fprintf(stderr, "ERROR: bad --memCap '%s' (expected e.g. 4G, 256M, 1024K, or a byte count)\n", optarg);
                    return 1;
                }
                break;
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
    const int64_t spill_budget = lift_spill_budget_bytes(mem_cap_cli);
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
        liftrec_sort(out, out_n);
        // Collapse consecutive same-seq positions into fixedStep blocks
        // (one header, just values).  For non-runs or short runs, keep
        // variableStep ("pos val" per row, one header per seq).  See
        // FS_THRESHOLD comment near liftrec_radix_sort.
        int64_t cur_var_seq = -1;  // -1 = no variableStep header open
        int64_t i = 0;
        while (i < out_n) {
            int64_t j = i;
            while (j + 1 < out_n &&
                   out[j+1].seq_idx == out[i].seq_idx &&
                   out[j+1].pos == out[j].pos + 1) {
                j++;
            }
            int64_t rl = j - i + 1;
            if (rl >= FS_THRESHOLD) {
                fprintf(fh, "fixedStep chrom=%s start=%" PRIi64 " step=1\n",
                        (const char*)stList_get(seqs, out[i].seq_idx),
                        out[i].pos + 1);
                for (int64_t k = i; k <= j; k++) fprintf(fh, "%g\n", out[k].value);
                cur_var_seq = -1;  // fixedStep closed; force new header for any var follow-up
            } else {
                if (out[i].seq_idx != cur_var_seq) {
                    fprintf(fh, "variableStep chrom=%s\n",
                            (const char*)stList_get(seqs, out[i].seq_idx));
                    cur_var_seq = out[i].seq_idx;
                }
                for (int64_t k = i; k <= j; k++) {
                    fprintf(fh, "%" PRIi64 " %g\n", out[k].pos + 1, out[k].value);
                }
            }
            i = j + 1;
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
