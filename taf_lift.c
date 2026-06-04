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
#include "chain.h"
#include <time.h>
#include <errno.h>
#include <signal.h>
#include <unistd.h>
#include <stdlib.h>

static void usage(void) {
    fprintf(stderr, "taffy lift [options]\n");
    fprintf(stderr, "Lift a .wig OR a BED annotation from universal MAF row-0 (ancestor) coords to a target leaf genome via the .tui index\n");
    fprintf(stderr, "-i --inputFile [FILE_NAME] : REQUIRED Path to the universal MAF/TAF (its .tui sidecar is expected at <input>.tui)\n");
    fprintf(stderr, "-w --wig       [FILE_NAME] : Input .wig in ancestor row-0 coords (chrom = full genome.sequence) -- OR --\n");
    fprintf(stderr, "-b --bed       [FILE_NAME] : Input BED3+ in ancestor row-0 coords (chrom = full genome.sequence).  Exactly one of -w / -b is required.\n");
    fprintf(stderr, "-g --genome    [STRING]    : REQUIRED Target genome name (e.g. hg38)\n");
    fprintf(stderr, "-o --outputFile [FILE_NAME] : Output (wig if input was wig, BED if input was BED; default stdout)\n");
    fprintf(stderr, "-m --memCap    [SIZE]      : (wig mode) Max in-memory output buffer before spilling to disk; suffix K/M/G allowed (default 2G).  Env TAFFY_LIFT_BUDGET_MB also honored (CLI wins).\n");
    fprintf(stderr, "-G --maxGap    [INT]       : (bed mode) Merge two adjacent target rows when the gap between them is <= INT bp.  Default 0 (touch / overlap only).  Max 2^60.  The merged interval spans the un-lifted gap region (no UCSC liftOver / halLiftover analogue).  Downstream tools that count target coverage will overcount by up to INT * number-of-merges.\n");
    fprintf(stderr, "-S --minSize   [INT]       : (bed mode) Drop output rows whose target length is < INT bp.  Default 1 (no filter).  Applied AFTER --maxGap so small fragments that would coalesce don't get filtered prematurely; reported in the summary log as 'dropped < --minSize'.\n");
    fprintf(stderr, "-F --fast                  : (bed mode) Use chunk-iteration lift instead of per-column open/close.  O(runs_in_range) vs O(columns_in_range); 10-1000x faster on chromosome-scale queries.  Output is equivalent modulo merge order (use with --maxGap N for browser-style block collapse).\n");
    fprintf(stderr, "-B --bin       [INT]       : (bed mode, requires --fast) Emit coarse-grained bedGraph: for every N-bp window on the TARGET genome, the value is the total source-bp lifted into that window.  Output is `seq<TAB>bin_start<TAB>bin_end<TAB>bp_covered`, sorted by (seq, bin).  Skips the per-row merge + filter machinery -- mutually exclusive with --maxGap / --minSize.  Intended for browser chromosome-scale tracks where per-base accuracy isn't needed.\n");
    fprintf(stderr, "-C --chainFilter [INT]     : (bed mode, requires --fast; mutex with --bin) Keep only the top N chains per BED record.  Buffers the visited runs, runs the taffy chainer (chain_open=0, chain_extend=1, max_gap=10Mb), and replays survivors through --maxGap / --minSize.  N=1 = primary chain only; useful for browser-snake-style filtering of paralog noise.  Drops are counted in the summary log.  Memory is O(n_runs_per_BED_record); a whole-chromosome record on a large alignment can buffer millions of runs (~50 bytes each).  Narrow the input BED if the working set must be bounded.\n");
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

/////////////////////////////////////////////////////////////////////////////
// BED lift path: input BED3+ in source-genome coords -> output BED with the
// same shape (BED3 in -> BED3 out; BED4+ in -> BED4+ out with target coords
// substituted for cols 0-2 and the target strand recomputed for col 5).
//
// For each input interval:
//   1. tui_query() maps the source seq range to a set of universal-column
//      intervals (sorted, disjoint, half-open).
//   2. Within each column interval, walk column by column.  For each column
//      tui_genome_lift_column() returns up to N target matches (paralogs).
//   3. Maintain a small fixed array of "open" target intervals.  Each match
//      either extends an existing open (same seq+strand, next consecutive
//      target pos) or starts a new one.  Opens not extended on a given
//      column are closed + emitted; opens still active at the end of a
//      column-interval are also closed (a gap in the universal-column
//      walk corresponds to non-monotonic / non-orthologous target coords
//      where we can't continue the run).
//
// Output is emitted in source-walk order, NOT globally sorted by target.
// Users wanting target-sorted output can pipe through `sort -k1,1 -k2,2n`.
/////////////////////////////////////////////////////////////////////////////

#define BED_MAX_OPEN 32   // max simultaneous open target intervals (= paralog cap)

typedef struct {
    const char *seq;      // borrowed from gl->seq_names; pointer-stable per gl
    int64_t     start;    // 0-based inclusive
    int64_t     end;      // 0-based exclusive
    int         strand;   // 1 = '+', 0 = '-'
    int         active;
} OpenBedInterval;

// Pending output line that may still merge with the next close_opens()
// flush.  Bucket key: (seq, output-strand) when emit_strand is true; (seq)
// alone when false.  In the latter case different alignment strands abut
// freely on the target, which is how halLiftover collapses BED3/4/5 output.
// `touched` is a monotonically-increasing serial assigned on each insert /
// extend; the LRU eviction picks the slot with the smallest `touched`.
typedef struct {
    const char *seq;     // borrowed from gl->seq_names; pointer-stable
    int64_t     start;
    int64_t     end;
    int64_t     touched;      // serial for LRU eviction
    char        out_strand;   // '+' / '-' if emit_strand, else 0
    int         active;
} PendingBed;

// Cap on the pending-merge buffer.  Sized to make LRU eviction rare on
// realistic paralog-dense workloads -- the prior 256-slot cap allowed
// enough LRU evictions on dense regions (apes hg38 chr5 paralog hotspot
// is the canonical case) that default and --fast paths lost DIFFERENT
// merge opportunities and emitted the same bp coverage with a slightly
// different row packing (off-by-a-few rows at chunk boundaries).
//
// Sizing:  the buffer holds open target intervals that might still merge
// with the next incoming run.  Empirical concurrent-slot counts:
//   apes / rodent / fish typical regions :   ~30-50
//   apes hg38 chr5 paralog hotspot       :  ~200-500
//   pathological whole-chrom (mexican-tetra
//     49 Mb single BED record)           :  >1M  (LRU eviction fires
//                                               here at any feasible cap;
//                                               documented as accepted,
//                                               bp exact)
// 1024 covers the realistic hotspot cases.  Going higher (8192/16384) does
// close the residual 1-3 split-vs-merged rows that LRU still produces, but
// costs 4-7x wall clock because pending_push / pending_cascade do linear
// O(p_cap) scans on every push -- the buffer is small enough at 1024 to
// stay in L1/L2 cache.  When the pool fills the LRU slot is flushed to
// disk and reused; the merge window shrinks but bp coverage is preserved.
#define PENDING_MAX 1024

// Flush one pending row to disk.  Caller clears `pe->active` after.
// `min_size` (>=1) drops rows whose length is below threshold without
// emitting; the drop is counted in *n_filtered_p (if non-NULL) so the
// final summary log can report how many rows the user's --minSize cap
// suppressed (otherwise a low output count looks like 'nothing mapped'
// when in fact it 'mapped but was filtered').
//
// Caller must decide whether min_size should fire here: the LRU-evict
// path in pending_push passes min_size=1 (i.e. emit unconditionally),
// because mid-line eviction is forced by capacity pressure, not by an
// actual end-of-merge-window event.  Only flush_pending at end-of-bed-
// line passes the user's min_size, because by then no further merges
// are possible.
static void emit_pending(FILE *fh, PendingBed *pe, int64_t min_size,
                         const char *name, const char *score,
                         int64_t *n_out_p, int64_t *n_filtered_p) {
    if (!pe->active) return;
    if ((pe->end - pe->start) < min_size) {
        if (n_filtered_p) (*n_filtered_p)++;
        return;
    }
    fprintf(fh, "%s\t%" PRIi64 "\t%" PRIi64, pe->seq, pe->start, pe->end);
    if (name) {
        fprintf(fh, "\t%s", name);
        if (score) {
            fprintf(fh, "\t%s", score);
            if (pe->out_strand) fprintf(fh, "\t%c", pe->out_strand);
        }
    }
    fputc('\n', fh);
    if (n_out_p) (*n_out_p)++;
}

// Try to merge a new tile into pending[]; if no abutting slot exists,
// allocate a fresh slot (growing pending[] in-place up to PENDING_MAX);
// if pending[] is at the cap, evict the LRU slot to disk and reuse it.
// Returns the (possibly reallocated) pending pointer; caller updates
// *p_cap.
//
// Three passes:
//   1. find a slot whose interval overlaps or is within max_gap bp of
//      the new tile on the same (seq, strand), via pending_overlap_or_
//      gap.  On match, extend the slot to the span union and cascade
//      to absorb other slots that now connect.  Bumps `touched` so
//      the slot survives the next eviction.
//   2. if no slot connected and pending[] has an inactive slot or
//      capacity to grow (cap < PENDING_MAX), allocate a fresh one.
//   3. if pending[] is at PENDING_MAX and all slots are active, find the
//      slot with the smallest `touched` (the LRU), flush it to disk, and
//      reuse the slot.  Linear scan; PENDING_MAX is bounded.
//
// Slots otherwise flush at end-of-bed-line via `flush_pending`.  Memory
// is bounded by PENDING_MAX × sizeof(PendingBed) ≈ 16 KB; correctness is
// preserved (an LRU slot that goes idle was about to be flushed at
// end-of-line anyway -- its content writes out the same).
// Cascade: after extending slot `hit`, keep absorbing any other slot
// whose edge now abuts hit's edge.  Required because the universal MAF
// can emit the SAME source region from two anchor blocks (one "main" +
// one "patch" covering 1-4 bp gaps between main's runs).  Without
// cascading, the main rows get pushed first (forming non-abutting slots
// in the pending buffer), then the patch rows extend each slot's end to
// abut the NEXT slot's start -- but no further merge fires, so the
// output keeps the per-block fragmentation the user shouldn't see.
// Unified merge predicate: a's interval and [b_start, b_end) overlap,
// touch, or have a between-edges gap <= max_gap (both directions).
// max_gap=0 reduces to "touch or overlap" -- one bp wider than strict
// equality, but in practice the lift never pushes two intervals that
// overlap on the same (seq, strand) at max_gap=0 (open[] tracks
// non-overlapping runs per paralog), so default behavior is byte-
// identical to the pre-patch strict-abut version.  With max_gap > 0
// the predicate also catches the "contained bridge" case: a bridge
// that arrives AFTER its containing slot has grown across a gap will
// be absorbed instead of emitted as a phantom sub-row.
//
// The merged interval is [min(a.start, b_start), max(a.end, b_end));
// at max_gap > 0 it spans the un-lifted gap region as if it were
// lifted -- this is the design choice of "collapse to syntenic block."
// (No UCSC liftOver / halLiftover analogue: liftOver gates at the
// chain-block level via -minBlocks; halLiftover has no --maxGap.)
// Downstream tools that count target-side coverage will overcount by
// up to max_gap × (number-of-merges); --minSize applies to the merged
// row, not the originals.
static int pending_overlap_or_gap(const PendingBed *a, int64_t b_start,
                                  int64_t b_end, int64_t max_gap) {
    // max_gap is CLI-bounded at parse time so b_end + max_gap and
    // a->end + max_gap cannot overflow int64 (chromosome coords are
    // sub-Tb; max_gap is capped at 1<<60).
    return b_end + max_gap >= a->start && a->end + max_gap >= b_start;
}

// Cascade: after a slot extends, absorb any other slot it now
// overlaps-or-gaps with.  Order-independent: the predicate is
// symmetric and the merge op (min start, max end) is monotone, so
// pending[]'s final state is independent of slot-scan order.
// Terminates: each merged=1 iteration deactivates one slot (j) and
// hit is never deactivated, so cascade runs <= cap iterations.
static void pending_cascade(PendingBed *pending, int cap, int hit,
                            int emit_strand, int64_t max_gap) {
    int merged;
    do {
        merged = 0;
        for (int j = 0; j < cap; j++) {
            if (j == hit || !pending[j].active) continue;
            if (pending[j].seq != pending[hit].seq) continue;
            if (emit_strand && pending[j].out_strand != pending[hit].out_strand) continue;
            if (pending_overlap_or_gap(&pending[hit], pending[j].start,
                                       pending[j].end, max_gap)) {
                if (pending[j].start < pending[hit].start) pending[hit].start = pending[j].start;
                if (pending[j].end   > pending[hit].end)   pending[hit].end   = pending[j].end;
                pending[j].active = 0;
                merged = 1;
                break;
            }
        }
    } while (merged);
}

static PendingBed *pending_push(FILE *fh, PendingBed *pending, int *p_cap,
                                int64_t *p_touch,
                                const char *seq, int64_t start, int64_t end,
                                char out_strand, int emit_strand,
                                int64_t max_gap, int64_t min_size,
                                const char *name, const char *score,
                                int64_t *n_out_p) {
    int64_t now = ++(*p_touch);
    // 1. extend-on-abut (or overlap, or within max_gap, both directions);
    //    cascade-merge after.  At max_gap=0 this matches strict-abut for
    //    non-overlapping intervals; with max_gap > 0 a previously-extended
    //    slot can also absorb a CONTAINED bridge that arrives later
    //    (otherwise the bridge would emit as a phantom sub-row inside the
    //    merged region).  See pending_overlap_or_gap above.
    for (int i = 0; i < *p_cap; i++) {
        if (!pending[i].active) continue;
        if (pending[i].seq != seq) continue;
        if (emit_strand && pending[i].out_strand != out_strand) continue;
        if (pending_overlap_or_gap(&pending[i], start, end, max_gap)) {
            if (start < pending[i].start) pending[i].start = start;
            if (end   > pending[i].end)   pending[i].end   = end;
            pending[i].touched = now;
            pending_cascade(pending, *p_cap, i, emit_strand, max_gap);
            return pending;
        }
    }
    // 2a. inactive slot in existing buffer
    int slot = -1;
    for (int i = 0; i < *p_cap; i++) {
        if (!pending[i].active) { slot = i; break; }
    }
    // 2b. grow up to the cap (doubling)
    if (slot < 0 && *p_cap < PENDING_MAX) {
        int new_cap = (*p_cap) * 2;
        if (new_cap > PENDING_MAX) new_cap = PENDING_MAX;
        pending = st_realloc(pending, (size_t)new_cap * sizeof(PendingBed));
        for (int i = *p_cap; i < new_cap; i++) pending[i].active = 0;
        slot = *p_cap;
        *p_cap = new_cap;
    }
    // 3. at the cap and all slots full -> evict LRU
    if (slot < 0) {
        int victim = 0;
        int64_t min_touch = pending[0].touched;
        for (int i = 1; i < *p_cap; i++) {
            if (pending[i].touched < min_touch) {
                min_touch = pending[i].touched;
                victim = i;
            }
        }
        // LRU eviction is forced by PENDING_MAX capacity, not by an
        // actual end-of-merge-window.  Emit unconditionally (min_size=1)
        // so a small fragment that hasn't yet had the chance to grow
        // isn't silently dropped by the size filter.  The user-visible
        // --minSize filter applies only at flush_pending below.
        emit_pending(fh, &pending[victim], 1, name, score, n_out_p, NULL);
        slot = victim;
    }
    pending[slot] = (PendingBed){ .seq = seq, .start = start, .end = end,
                                  .touched = now,
                                  .out_strand = out_strand, .active = 1 };
    return pending;
}

// Close open BED intervals matching the predicate; rows feed `pending[]`
// for adjacent-merge before they hit disk.  `used_this_col` is a per-call
// mask: if NULL, every active open is closed; if non-NULL, only active
// opens NOT in the mask are closed (i.e. those that didn't get extended
// this column).  Output strand for each closed interval is the XOR of
// `input_strand_sign` (+1 / -1) and the open's alignment strand;
// `emit_strand` controls whether the strand column appears at all (set
// false for BED3/4/5 in -> no strand out, or input strand '.').
static void close_opens(FILE *fh, OpenBedInterval *open, int n_open,
                        const int *used_this_col, int input_strand_sign,
                        int emit_strand,
                        PendingBed **p_pending, int *p_cap, int64_t *p_touch,
                        int64_t max_gap, int64_t min_size,
                        const char *name, const char *score,
                        int64_t *n_out_p) {
    for (int s = 0; s < n_open; s++) {
        if (!open[s].active) continue;
        if (used_this_col && used_this_col[s]) continue;
        char out_strand = 0;
        if (emit_strand) {
            int sign = (open[s].strand ? +1 : -1) * input_strand_sign;
            out_strand = (sign >= 0) ? '+' : '-';
        }
        *p_pending = pending_push(fh, *p_pending, p_cap, p_touch,
                                  open[s].seq, open[s].start, open[s].end,
                                  out_strand, emit_strand,
                                  max_gap, min_size,
                                  name, score, n_out_p);
        open[s].active = 0;
    }
}

// Drain all pending rows to disk (end-of-bed-line, or any point where
// the caller has no more tiles to add).  Leaves `pending[]` empty.
// `min_size` (>=1) filters at emit; rows whose length < min_size are
// counted into *n_filtered_p instead of emitted, so the BED summary
// log can distinguish "filtered" from "unmapped".  This is the only
// site where --minSize fires -- the LRU-evict path in pending_push
// bypasses it (see comment there).
static void flush_pending(FILE *fh, PendingBed *pending, int n_pending,
                          int64_t min_size,
                          const char *name, const char *score,
                          int64_t *n_out_p, int64_t *n_filtered_p) {
    for (int i = 0; i < n_pending; i++) {
        emit_pending(fh, &pending[i], min_size, name, score,
                     n_out_p, n_filtered_p);
        pending[i].active = 0;
    }
}

// --- Chunk-iteration BED lift (--fast) ----------------------------------
//
// The default bed_lift_main_impl below walks every source column in the
// query range, asking "which run(s) cover this column?" via the per-column
// open/close machinery.  For a whole-chromosome query (~10^8 columns)
// that's the dominant cost.
//
// This variant inverts the loop: for each col_iv returned by tui_query,
// walk the chunks whose [g_min, g_max) intersects the col_iv, decode each
// (once, cached), and emit each contained run as a target-coord row
// directly.  Cost is O(n_runs_in_range) instead of O(n_columns_in_range),
// typically 10^3-10^5x fewer iterations on chromosome-scale queries.
//
// The pending_push/cascade/flush_pending machinery is reused so --maxGap
// / --minSize work the same.  Output rows should be equivalent to the
// column-walk modulo merge order: every chunk-walk emit is exactly one
// gap-free target run, which is what an open/close cycle produces in
// the column-walk -- the two paths agree on the set of pre-merge rows,
// then pending_push collapses them identically.
//
// Strand math (clip a run to col_iv [c_lo, c_hi)):
//   clip_start = max(g_start, c_lo);  clip_end = min(g_start+length, c_hi);
//   '+' strand: t_out = [t_start + (clip_start - g_start),
//                        t_start + (clip_end   - g_start))
//   '-' strand: t_out = [t_start + length - (clip_end   - g_start),
//                        t_start + length - (clip_start - g_start))
// Output BED strand: (run.strand=='+' XOR input.strand=='-').

// Context passed through the tui visitor callback for the chunk lift.
typedef struct {
    FILE       *fh;
    PendingBed *pending;
    int         pending_cap;
    int64_t     pending_touch;
    int64_t     c_lo, c_hi;
    int         input_strand_sign;
    int         emit_strand;
    int64_t     max_gap, min_size;
    const char *name;
    const char *score;
    int64_t    *n_out_p;
    // --bin coarse-grained bedGraph mode.  When bin_size > 0 the visitor
    // does NOT touch pending[]/pending_push at all -- it accumulates a
    // flat array of (target_seq, bin_idx, bp) tuples here, which is
    // sorted, merged, and emitted as bedGraph at end-of-lift.
    int64_t     bin_size;
    struct BinEntry {
        const char *seq;   // borrowed from gl->seq_names; pointer-stable per gl
        int64_t     bin;   // bin_idx = target_pos / bin_size
        int64_t     bp;    // source bp contributed to this bin by one clipped run
    } *bins;
    size_t      n_bins, cap_bins;
    // --chainFilter: keep only the top-N highest-scoring chains.  When
    // chain_filter > 0 the visitor BUFFERS each clipped run instead of
    // calling pending_push; bed_lift_chunk_impl runs the chain at the
    // end of each BED record, picks the top N chains by total_score,
    // and replays surviving entries through pending_push so --maxGap /
    // --minSize / strand merge still fire on the kept subset.
    int64_t     chain_filter;
    struct LiftRun {
        const char *seq;          // target seq (gl-owned, pointer stable)
        int64_t     c_start, c_end; // column range (chain q-axis)
        int64_t     t_start, t_end; // target genome pos (chain t-axis, forward)
        int         strand;         // run strand: +1 = forward, -1 = reverse
        char        out_strand;     // already-XOR'd output strand char ('+'/'-'/0)
    } *run_buf;
    size_t      n_run_buf, cap_run_buf;
} ChunkLiftCtx;

// Visitor callback: clip the run to [c_lo, c_hi), map to target coords
// honouring strand, and push the (target_seq, t_start, t_end) interval
// through pending_push (so --maxGap merging + --minSize filter still
// fire at flush time).  Strand math:
//   '+': target advances with source.  Output [t + (cs - gs), t + (ce - gs))
//   '-': target advances inverse.       Output [t + len - (ce - gs),
//                                                t + len - (cs - gs))
// Output strand bit = run.strand XOR input.strand sign (matches the
// open/close path's per-column derivation exactly).
static void chunk_lift_visit_cb(const TuiRun *r, void *user) {
    ChunkLiftCtx *cx = (ChunkLiftCtx *)user;
    int64_t r_end = r->g_start + r->length;
    int64_t cs = r->g_start > cx->c_lo ? r->g_start : cx->c_lo;
    int64_t ce = r_end < cx->c_hi     ? r_end     : cx->c_hi;
    if (cs >= ce) return;                          // no overlap (defensive)
    int64_t t_out_start, t_out_end;
    if (r->strand) {
        t_out_start = r->t_start + (cs - r->g_start);
        t_out_end   = r->t_start + (ce - r->g_start);
    } else {
        t_out_start = r->t_start + r->length - (ce - r->g_start);
        t_out_end   = r->t_start + r->length - (cs - r->g_start);
    }
    // --bin coarse-grained bedGraph path: accumulate per-(target_seq, bin)
    // bp directly, skip the pending_push / open-merge machinery entirely.
    // Strand is irrelevant for bp coverage so we ignore emit_strand here.
    if (cx->bin_size > 0) {
        int64_t bs = cx->bin_size;
        int64_t bin_lo = t_out_start / bs;
        int64_t bin_hi = (t_out_end - 1) / bs;
        for (int64_t b = bin_lo; b <= bin_hi; b++) {
            int64_t b_start = b * bs;
            int64_t b_end   = b_start + bs;
            int64_t cs2 = (t_out_start > b_start) ? t_out_start : b_start;
            int64_t ce2 = (t_out_end   < b_end)   ? t_out_end   : b_end;
            int64_t bp = ce2 - cs2;
            if (bp <= 0) continue;
            if (cx->n_bins == cx->cap_bins) {
                cx->cap_bins = cx->cap_bins ? cx->cap_bins * 2 : 4096;
                cx->bins = st_realloc(cx->bins, cx->cap_bins * sizeof(*cx->bins));
            }
            cx->bins[cx->n_bins].seq = r->seq;
            cx->bins[cx->n_bins].bin = b;
            cx->bins[cx->n_bins].bp  = bp;
            cx->n_bins++;
        }
        return;
    }
    char out_strand = 0;
    if (cx->emit_strand) {
        int sign = (r->strand ? +1 : -1) * cx->input_strand_sign;
        out_strand = (sign >= 0) ? '+' : '-';
    }
    // --chainFilter path: buffer for the post-record chain pass.
    // Don't touch pending_push; replay happens after chaining.
    if (cx->chain_filter > 0) {
        if (cx->n_run_buf == cx->cap_run_buf) {
            cx->cap_run_buf = cx->cap_run_buf ? cx->cap_run_buf * 2 : 4096;
            cx->run_buf = st_realloc(cx->run_buf,
                                     cx->cap_run_buf * sizeof(*cx->run_buf));
        }
        struct LiftRun *L = &cx->run_buf[cx->n_run_buf++];
        L->seq        = r->seq;
        L->c_start    = cs;
        L->c_end      = ce;
        L->t_start    = t_out_start;
        L->t_end      = t_out_end;
        L->strand     = r->strand ? +1 : -1;
        L->out_strand = out_strand;
        return;
    }
    cx->pending = pending_push(cx->fh, cx->pending, &cx->pending_cap,
                               &cx->pending_touch,
                               r->seq, t_out_start, t_out_end,
                               out_strand, cx->emit_strand,
                               cx->max_gap, cx->min_size,
                               cx->name, cx->score, cx->n_out_p);
}

// Comparator for the bin-accumulator array: sort by (seq name lexically,
// bin index ascending).  seq names are interned in gl->seq_names so
// pointer-equal seq strings compare equal under strcmp too -- the
// post-sort merge can use either pointer or strcmp equality.
static int bin_entry_cmp(const void *a, const void *b) {
    const struct BinEntry *x = a, *y = b;
    int c = strcmp(x->seq, y->seq);
    if (c) return c;
    return (x->bin < y->bin) ? -1 : (x->bin > y->bin) ? 1 : 0;
}

static int bed_lift_chunk_impl(Tui *tui, TuiGenomeLift *gl,
                               const char *bed_file, const char *output_file,
                               int64_t max_gap, int64_t min_size,
                               int64_t bin_size,
                               int64_t chain_filter) {
    FILE *bf = fopen(bed_file, "r");
    if (bf == NULL) {
        fprintf(stderr, "ERROR: failed to open bed %s: %s\n", bed_file, strerror(errno));
        return 1;
    }
    FILE *fh = output_file ? fopen(output_file, "w") : stdout;
    if (fh == NULL) {
        fprintf(stderr, "ERROR: failed to open output %s: %s\n", output_file, strerror(errno));
        fclose(bf);
        return 1;
    }

    char *line = NULL;
    size_t cap = 0;
    int64_t lineno = 0, n_in = 0, n_out = 0, n_unmapped = 0, n_filtered = 0;
    int64_t n_chain_filtered = 0;  // runs dropped by chainFilter (only meaningful when chain_filter > 0)
    time_t t0 = time(NULL);
    ChunkLiftCtx cx;
    memset(&cx, 0, sizeof(cx));
    cx.fh = fh;
    cx.max_gap = max_gap;
    cx.min_size = min_size;
    cx.bin_size = bin_size;
    cx.chain_filter = chain_filter;
    cx.n_out_p = &n_out;

    ssize_t got;
    while ((got = getline(&line, &cap, bf)) > 0) {
        lineno++;
        if (line[0] == '#' || line[0] == '\n' || line[0] == 0) continue;
        while (got > 0 && (line[got-1] == '\n' || line[got-1] == '\r')) line[--got] = 0;
        char *fields[6] = {0};
        int n_fields = 0;
        char *p = line;
        while (n_fields < 6 && p) {
            fields[n_fields++] = p;
            char *t = strchr(p, '\t');
            if (t) { *t = 0; p = t + 1; } else p = NULL;
        }
        if (n_fields < 3) {
            fprintf(stderr, "WARN: skipping bed line %" PRIi64 " (need >=3 cols, got %d)\n",
                    lineno, n_fields);
            continue;
        }
        const char *chrom  = fields[0];
        int64_t start = atoll(fields[1]);
        int64_t end   = atoll(fields[2]);
        const char *name   = (n_fields >= 4) ? fields[3] : NULL;
        const char *score  = (n_fields >= 5) ? fields[4] : NULL;
        const char *strand = (n_fields >= 6) ? fields[5] : NULL;
        cx.input_strand_sign = (strand && strand[0] == '-') ? -1 : +1;
        cx.emit_strand = (strand != NULL && strand[0] != '.');
        cx.name = name;
        cx.score = score;

        n_in++;
        int64_t n_iv = 0;
        TuiInterval *iv = tui_query(tui, chrom, start, end, &n_iv);
        if (iv == NULL || n_iv == 0) {
            free(iv);
            n_unmapped++;
            continue;
        }
        if (cx.bin_size == 0 && cx.pending_cap == 0) {
            cx.pending_cap = BED_MAX_OPEN;
            cx.pending = st_calloc((size_t)cx.pending_cap, sizeof(PendingBed));
        }
        if (cx.bin_size == 0)
            for (int s = 0; s < cx.pending_cap; s++) cx.pending[s].active = 0;

        for (int64_t k = 0; k < n_iv; k++) {
            cx.c_lo = iv[k].start;
            cx.c_hi = iv[k].end;
            tui_genome_lift_visit_runs(gl, cx.c_lo, cx.c_hi,
                                       chunk_lift_visit_cb, &cx);
        }

        // --chainFilter post-visit pass: chain the buffered runs of THIS
        // record, pick the top-N chains by total_score, and replay the
        // survivors through pending_push so --maxGap / --minSize / strand
        // merge still fire on the kept subset.
        if (cx.chain_filter > 0 && cx.n_run_buf > 0) {
            int64_t n_buf = (int64_t) cx.n_run_buf;
            TaffyAln *alns = st_malloc((size_t) n_buf * sizeof(TaffyAln));
            for (int64_t i = 0; i < n_buf; i++) {
                struct LiftRun *L = &cx.run_buf[i];
                alns[i].q_name  = chrom;
                alns[i].q_start = L->c_start;
                alns[i].q_end   = L->c_end;
                alns[i].t_name  = L->seq;
                alns[i].t_start = L->t_start;
                alns[i].t_end   = L->t_end;
                alns[i].strand  = L->strand;
                alns[i].score   = L->c_end - L->c_start;
                alns[i].user    = (void *) (intptr_t) i;   // buf index
            }
            TaffyChainCostParams cost = {
                TAFFY_CHAIN_DEFAULT_OPEN, TAFFY_CHAIN_DEFAULT_EXTEND
            };
            int64_t *chain_id = st_calloc((size_t) n_buf, sizeof(int64_t));
            TaffyChainInfo *chains_out = NULL;
            int64_t n_chains_out = 0;
            taffy_chain(alns, n_buf,
                        taffy_chain_default_gap_cost, &cost,
                        TAFFY_CHAIN_DEFAULT_MAX_GAP,
                        chain_id, &chains_out, &n_chains_out);

            // Mark the top-N chain ids by score; taffy_chain returns
            // chains_out sorted desc, so chains_out[0..min(N,nc)) are
            // the winners.  Use a flat bool[] keyed on (max id + 1) for
            // O(1) membership.
            int64_t max_id = 0;
            for (int64_t k = 0; k < n_chains_out; k++)
                if (chains_out[k].id > max_id) max_id = chains_out[k].id;
            char *keep_chain = st_calloc((size_t)(max_id + 1), sizeof(char));
            int64_t n_keep_chains = cx.chain_filter < n_chains_out
                                       ? cx.chain_filter : n_chains_out;
            for (int64_t k = 0; k < n_keep_chains; k++)
                keep_chain[chains_out[k].id] = 1;

            // taffy_chain SORTS alns in place by (q_name, strand, q_start),
            // destroying the input order.  Walk alns post-sort and use
            // aln.user (the original buf index) to recover the LiftRun.
            // Set up pending state if not already.
            if (cx.pending_cap == 0) {
                cx.pending_cap = BED_MAX_OPEN;
                cx.pending = st_calloc((size_t)cx.pending_cap, sizeof(PendingBed));
            }
            for (int s = 0; s < cx.pending_cap; s++) cx.pending[s].active = 0;

            for (int64_t i = 0; i < n_buf; i++) {
                if (!keep_chain[chain_id[i]]) { n_chain_filtered++; continue; }
                int64_t bi = (int64_t)(intptr_t) alns[i].user;
                struct LiftRun *L = &cx.run_buf[bi];
                cx.pending = pending_push(cx.fh, cx.pending, &cx.pending_cap,
                                          &cx.pending_touch,
                                          L->seq, L->t_start, L->t_end,
                                          L->out_strand, cx.emit_strand,
                                          cx.max_gap, cx.min_size,
                                          cx.name, cx.score, cx.n_out_p);
            }

            free(keep_chain);
            free(chain_id);
            free(chains_out);
            free(alns);

            // Reset buffer for the next BED record (keep allocation).
            cx.n_run_buf = 0;
        }

        if (cx.bin_size == 0) {
            // BED-row mode: flush this record's open intervals.  Bin mode
            // accumulates ACROSS bed records and flushes once at the end.
            flush_pending(cx.fh, cx.pending, cx.pending_cap, min_size,
                          name, score, &n_out, &n_filtered);
            cx.pending_touch = 0;
        }
        free(iv);
    }

    if (cx.bin_size > 0) {
        // bedGraph emit: sort all bin entries by (seq, bin), then merge-
        // and-emit consecutive same-key entries.  All BED input records
        // contribute to one combined bedGraph (lift "this whole region of
        // source -> coverage on target"), so output is *across* records,
        // not per-record.
        qsort(cx.bins, cx.n_bins, sizeof(*cx.bins), bin_entry_cmp);
        size_t i = 0;
        while (i < cx.n_bins) {
            const char *seq = cx.bins[i].seq;
            int64_t bin = cx.bins[i].bin, bp_total = 0;
            size_t j = i;
            // Pointer equality is safe here: seq strings are interned in
            // gl->seq_names, so equal strings (post-sort neighbours) point
            // to the same buffer.  Falls back to strcmp implicitly via the
            // sort -- if pointer-equal fails, the bin index also differs.
            while (j < cx.n_bins && cx.bins[j].seq == seq && cx.bins[j].bin == bin) {
                bp_total += cx.bins[j].bp;
                j++;
            }
            fprintf(fh, "%s\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\n",
                    seq, bin * cx.bin_size, (bin + 1) * cx.bin_size, bp_total);
            n_out++;
            i = j;
        }
        free(cx.bins);
    }

    free(cx.pending);
    free(cx.run_buf);
    free(line);
    fclose(bf);
    if (output_file) fclose(fh);

    if (cx.bin_size > 0) {
        st_logInfo("BED lift (--fast --bin %" PRIi64 "): %" PRIi64 " input -> %" PRIi64 " bedGraph rows "
                   "(%" PRIi64 " unmapped) in %" PRIi64 " s\n",
                   cx.bin_size, n_in, n_out, n_unmapped,
                   (int64_t)(time(NULL) - t0));
    } else if (cx.chain_filter > 0) {
        st_logInfo("BED lift (--fast --chainFilter %" PRIi64 "): %" PRIi64 " input -> "
                   "%" PRIi64 " output intervals (%" PRIi64 " unmapped, "
                   "%" PRIi64 " dropped < --minSize, %" PRIi64 " dropped by chain) "
                   "in %" PRIi64 " s\n",
                   cx.chain_filter, n_in, n_out, n_unmapped, n_filtered,
                   n_chain_filtered, (int64_t)(time(NULL) - t0));
    } else {
        st_logInfo("BED lift (--fast): %" PRIi64 " input -> %" PRIi64 " output intervals "
                   "(%" PRIi64 " unmapped, %" PRIi64 " dropped < --minSize) "
                   "in %" PRIi64 " s\n",
                   n_in, n_out, n_unmapped, n_filtered,
                   (int64_t)(time(NULL) - t0));
    }
    return 0;
}

static int bed_lift_main_impl(Tui *tui, TuiGenomeLift *gl,
                              const char *bed_file, const char *output_file,
                              int64_t max_gap, int64_t min_size) {
    FILE *bf = fopen(bed_file, "r");
    if (bf == NULL) {
        fprintf(stderr, "ERROR: failed to open bed %s: %s\n", bed_file, strerror(errno));
        return 1;
    }
    FILE *fh = output_file ? fopen(output_file, "w") : stdout;
    if (fh == NULL) {
        fprintf(stderr, "ERROR: failed to open output %s: %s\n", output_file, strerror(errno));
        fclose(bf);
        return 1;
    }

    // Open-interval tracker.  BED_MAX_OPEN is the starting stack size --
    // enough for typical paralog density (a handful) but real apes/primate
    // chr22 paralog hotspots produce 40+ matches at a single column, so we
    // grow heap-backed as needed.  Cached across input BED lines.
    OpenBedInterval *open = NULL;
    int *used_this_col = NULL;
    PendingBed *pending = NULL;
    int open_cap = 0;
    int pending_cap = 0;
    int64_t pending_touch = 0;          // monotone serial for LRU eviction
    char *line = NULL;
    size_t cap = 0;
    int64_t lineno = 0, n_in = 0, n_out = 0, n_unmapped = 0, n_filtered = 0;
    time_t t0 = time(NULL);

    ssize_t got;
    while ((got = getline(&line, &cap, bf)) > 0) {
        lineno++;
        if (line[0] == '#' || line[0] == '\n' || line[0] == 0) continue;
        // Strip trailing newline / CR.
        while (got > 0 && (line[got-1] == '\n' || line[got-1] == '\r')) line[--got] = 0;
        // Tokenize by tab.  Cols 0-2 required; 3 (name), 4 (score), 5 (strand)
        // optional.  Cols >= 6 (thickStart etc.) ignored.
        char *fields[6] = {0};
        int n_fields = 0;
        char *p = line;
        while (n_fields < 6 && p) {
            fields[n_fields++] = p;
            char *t = strchr(p, '\t');
            if (t) { *t = 0; p = t + 1; } else p = NULL;
        }
        if (n_fields < 3) {
            fprintf(stderr, "WARN: skipping bed line %" PRIi64 " (need >=3 cols, got %d)\n",
                    lineno, n_fields);
            continue;
        }
        const char *chrom  = fields[0];
        int64_t start = atoll(fields[1]);
        int64_t end   = atoll(fields[2]);
        const char *name   = (n_fields >= 4) ? fields[3] : NULL;
        const char *score  = (n_fields >= 5) ? fields[4] : NULL;
        const char *strand = (n_fields >= 6) ? fields[5] : NULL;
        int input_strand_sign = (strand && strand[0] == '-') ? -1 : +1;
        // BED convention: strand of '.' means "unknown / unstranded".  Treat
        // it as no-strand-column (skip the strand on output) rather than
        // forcing a definite output direction derived purely from the
        // alignment.  All other strand values produce a definite output.
        int emit_strand = (strand != NULL && strand[0] != '.');

        n_in++;

        // tui_query expects the FULL "<genome>.<chrom>" key.
        int64_t n_iv = 0;
        TuiInterval *iv = tui_query(tui, chrom, start, end, &n_iv);
        if (iv == NULL || n_iv == 0) {
            free(iv);
            n_unmapped++;
            continue;
        }

        // Lazily (re)allocate open / used_this_col / pending to the
        // initial cap on the first BED line; later columns grow open[]
        // (and pending[] in lock-step for new buckets) on demand.
        // pending[] holds rows that have been closed but might still merge
        // with an abutting close later in the same BED line; its capacity
        // tracks the worst-case unmerged row count per line (not just the
        // simultaneously-open count, which is much smaller), so it grows
        // independently of open_cap.
        if (open_cap == 0) {
            open_cap = BED_MAX_OPEN;
            open = st_malloc((size_t)open_cap * sizeof(OpenBedInterval));
            used_this_col = st_malloc((size_t)open_cap * sizeof(int));
        }
        if (pending_cap == 0) {
            pending_cap = BED_MAX_OPEN;
            pending = st_calloc((size_t)pending_cap, sizeof(PendingBed));
        }
        for (int s = 0; s < open_cap; s++) open[s].active = 0;
        for (int s = 0; s < pending_cap; s++) pending[s].active = 0;

        // Match buffer is reused ACROSS iv[k]s -- holding it inside the iv
        // loop forced a free()+malloc() per column-interval AND threw away
        // any growth, hurting nothing here but cleaner this way.
        TuiGenomeMatch *m = NULL;
        int m_cap = 0;
        for (int64_t k = 0; k < n_iv; k++) {
            int64_t c_lo = iv[k].start;
            int64_t c_hi = iv[k].end;
            // Per-column lift; group consecutive same-(seq,strand) hits.
            // After each column we close opens that weren't extended.
            for (int64_t c = c_lo; c < c_hi; c++) {
                if (m_cap == 0) {
                    m_cap = BED_MAX_OPEN;
                    m = st_malloc((size_t)m_cap * sizeof(TuiGenomeMatch));
                }
                int nm = tui_genome_lift_column(gl, c, m, m_cap);
                if (nm == 0) {
                    // No target hit at this column -- close all opens.
                    close_opens(fh, open, open_cap, NULL,
                                input_strand_sign, emit_strand,
                                &pending, &pending_cap, &pending_touch,
                                max_gap, min_size,
                                name, score, &n_out);
                    continue;
                }
                if (nm > m_cap) {
                    // More paralogs than buf holds (apes-density chr22
                    // paralog hotspots can have 40+ at one column).
                    // Grow + re-query so we see them all.
                    m_cap = nm;
                    m = st_realloc(m, (size_t)m_cap * sizeof(TuiGenomeMatch));
                    tui_genome_lift_column(gl, c, m, m_cap);
                }
                // open[]/used[] slot pool.  Worst case at this column is
                // (active opens NOT extended) + (new matches needing slots),
                // both bounded by `prior_active + nm`.  We can't cheaply count
                // active opens here, so just ensure the pool grows on demand
                // in the inactive-slot search below when needed.
                // For each match, find an extending open or open a new one.
                // Extension semantics differ by strand:
                //   '+' (m.strand=1): target pos advances; next expected =
                //                     open.end, and on hit we open.end++.
                //   '-' (m.strand=0): target pos decreases as source advances
                //                     (reverse-complement alignment); next
                //                     expected = open.start - 1, and on hit
                //                     we open.start-- to grow the interval
                //                     downward.
                //
                // tui_genome_lift_column returns matches in an order
                // determined by chunk/run layout, not column content (see
                // chunk_collect in tui.c).  A given paralog pair (A,B)
                // shows up in the same relative order at every column they
                // co-occur, so an open's slot remains stable across the
                // paralog's whole run.
                memset(used_this_col, 0, (size_t)open_cap * sizeof(int));
                for (int i = 0; i < nm; i++) {
                    int m_strand = (m[i].strand ? 1 : 0);
                    int slot = -1;
                    for (int s = 0; s < open_cap; s++) {
                        if (!open[s].active || used_this_col[s]) continue;
                        if (open[s].seq != m[i].seq) continue;
                        if (open[s].strand != m_strand) continue;
                        int extends = m_strand
                            ? (open[s].end == m[i].pos)
                            : (open[s].start == m[i].pos + 1);
                        if (extends) { slot = s; break; }
                    }
                    if (slot >= 0) {
                        if (m_strand) open[slot].end++;
                        else          open[slot].start--;
                        used_this_col[slot] = 1;
                    } else {
                        // Find an inactive slot; grow the pool if none exists
                        // (prior column's opens may all still be active and
                        // not yet extended-and-marked this column, e.g. when
                        // every paralog starts a fresh run at this column).
                        for (int s = 0; s < open_cap; s++) {
                            if (!open[s].active) { slot = s; break; }
                        }
                        if (slot < 0) {
                            int new_cap = open_cap * 2;
                            open = st_realloc(open, (size_t)new_cap * sizeof(OpenBedInterval));
                            used_this_col = st_realloc(used_this_col, (size_t)new_cap * sizeof(int));
                            for (int s = open_cap; s < new_cap; s++) {
                                open[s].active = 0;
                                used_this_col[s] = 0;
                            }
                            slot = open_cap;
                            open_cap = new_cap;
                        }
                        open[slot] = (OpenBedInterval){
                            .seq = m[i].seq,
                            .start = m[i].pos,
                            .end = m[i].pos + 1,
                            .strand = (m[i].strand ? 1 : 0),
                            .active = 1,
                        };
                        used_this_col[slot] = 1;
                    }
                }
                // Close opens that weren't extended this column.
                close_opens(fh, open, open_cap, used_this_col,
                            input_strand_sign, emit_strand,
                            &pending, &pending_cap, &pending_touch,
                            max_gap, min_size,
                            name, score, &n_out);
            }
            // No close_opens at iv[k]/iv[k+1] boundary: opens carry over so
            // that target intervals which happen to remain contiguous across
            // a source-side column gap (the source has a run break here but
            // the target's own run spans the gap with a matching coordinate
            // step) merge into one output BED line instead of two.  The
            // next column's extend-vs-close check handles correctness:
            // extension fires only when end==m.pos for '+' (resp. start==
            // m.pos+1 for '-'), so if the target's coord doesn't actually
            // line up after the gap, the open closes at the first column of
            // iv[k+1] and a fresh one opens in its place -- same output as
            // before, just deferred by one column.  In practice (rodents +
            // fish vs apes) this cuts output BED line count ~100-200x by
            // collapsing the per-source-block fragmentation that the
            // universal-column run topology introduced.
        }
        // End of this BED input line -- close any still-open intervals,
        // then drain the pending-merge buffer to disk.
        close_opens(fh, open, open_cap, NULL,
                    input_strand_sign, emit_strand,
                    &pending, &pending_cap, &pending_touch,
                    max_gap, min_size,
                    name, score, &n_out);
        flush_pending(fh, pending, pending_cap, min_size,
                      name, score, &n_out, &n_filtered);
        // Reset the touch counter between BED input lines so we don't drift
        // toward int64 overflow on million-line inputs (mod-arithmetic LRU
        // would also work but reset is simpler -- per-line scope is fine
        // since flush_pending also clears active=0 on every slot).
        pending_touch = 0;
        free(m);
        free(iv);
    }
    free(open); free(used_this_col); free(pending);
    free(line);
    fclose(bf);
    if (output_file) fclose(fh);

    st_logInfo("BED lift: %" PRIi64 " input -> %" PRIi64 " output intervals "
               "(%" PRIi64 " unmapped, %" PRIi64 " dropped < --minSize) "
               "in %" PRIi64 " s\n",
               n_in, n_out, n_unmapped, n_filtered,
               (int64_t)(time(NULL) - t0));
    return 0;
}

int taf_lift_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);
    char *logLevelString = NULL;
    char *maf_file       = NULL;
    char *wig_file       = NULL;
    char *bed_file       = NULL;
    char *target_genome  = NULL;
    char *output_file    = NULL;
    int64_t mem_cap_cli  = -1;     // -1 sentinel = "not set on CLI"
    int64_t max_gap      = 0;      // --maxGap (bed mode), 0 = strict abut
    int64_t min_size     = 1;      // --minSize (bed mode), 1 = no filter
    int     fast_mode    = 0;      // --fast (bed mode), 0 = legacy column walk
    int64_t bin_size     = 0;      // --bin (bed mode, requires --fast), 0 = off
    int64_t chain_filter = 0;      // --chainFilter N (bed mode, requires --fast), 0 = off

    while (1) {
        static struct option long_options[] = {
            { "logLevel",   required_argument, 0, 'l' },
            { "inputFile",  required_argument, 0, 'i' },
            { "wig",        required_argument, 0, 'w' },
            { "bed",        required_argument, 0, 'b' },
            { "genome",     required_argument, 0, 'g' },
            { "outputFile", required_argument, 0, 'o' },
            { "memCap",     required_argument, 0, 'm' },
            { "maxGap",     required_argument, 0, 'G' },
            { "minSize",    required_argument, 0, 'S' },
            { "fast",       no_argument,       0, 'F' },
            { "bin",        required_argument, 0, 'B' },
            { "chainFilter", required_argument, 0, 'C' },
            { "help",       no_argument,       0, 'h' },
            { 0, 0, 0, 0 }
        };
        int idx = 0;
        int key = getopt_long(argc, argv, "l:i:w:b:g:o:m:G:S:FB:C:h", long_options, &idx);
        if (key == -1) break;
        switch (key) {
            case 'l': logLevelString = optarg; break;
            case 'i': maf_file       = optarg; break;
            case 'w': wig_file       = optarg; break;
            case 'b': bed_file       = optarg; break;
            case 'g': target_genome  = optarg; break;
            case 'o': output_file    = optarg; break;
            case 'm':
                mem_cap_cli = lift_parse_size(optarg);
                if (mem_cap_cli <= 0) {
                    fprintf(stderr, "ERROR: bad --memCap '%s' (expected e.g. 4G, 256M, 1024K, or a byte count)\n", optarg);
                    return 1;
                }
                break;
            case 'G': {
                // strtoll instead of atoll: detect garbage / overflow.
                // Upper-bound cap so `end + max_gap` cannot wrap int64
                // (chrom coords are sub-Tb, 1<<60 is well clear of any
                // realistic sum).
                errno = 0;
                char *end_p = NULL;
                long long v = strtoll(optarg, &end_p, 10);
                if (errno == ERANGE || end_p == optarg || *end_p != 0 ||
                    v < 0 || v > ((long long)1 << 60)) {
                    fprintf(stderr, "ERROR: --maxGap must be in [0, 2^60] "
                                    "(got '%s')\n", optarg);
                    return 1;
                }
                max_gap = (int64_t)v;
                break;
            }
            case 'S': {
                errno = 0;
                char *end_p = NULL;
                long long v = strtoll(optarg, &end_p, 10);
                if (errno == ERANGE || end_p == optarg || *end_p != 0 ||
                    v < 1) {
                    fprintf(stderr, "ERROR: --minSize must be >= 1 "
                                    "(got '%s')\n", optarg);
                    return 1;
                }
                min_size = (int64_t)v;
                break;
            }
            case 'F': fast_mode = 1; break;
            case 'B': {
                errno = 0;
                char *end_p = NULL;
                long long v = strtoll(optarg, &end_p, 10);
                if (errno == ERANGE || end_p == optarg || *end_p != 0 ||
                    v < 1 || v > ((long long)1 << 60)) {
                    fprintf(stderr, "ERROR: --bin must be in [1, 2^60] "
                                    "(got '%s')\n", optarg);
                    return 1;
                }
                bin_size = (int64_t)v;
                break;
            }
            case 'C': {
                errno = 0;
                char *end_p = NULL;
                long long v = strtoll(optarg, &end_p, 10);
                if (errno == ERANGE || end_p == optarg || *end_p != 0 ||
                    v < 1 || v > ((long long)1 << 31)) {
                    fprintf(stderr, "ERROR: --chainFilter must be in [1, 2^31] "
                                    "(got '%s')\n", optarg);
                    return 1;
                }
                chain_filter = (int64_t)v;
                break;
            }
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }
    st_setLogLevelFromString(logLevelString);

    if (maf_file == NULL || target_genome == NULL) {
        fprintf(stderr, "ERROR: -i and -g are required\n");
        usage();
        return 1;
    }
    if ((wig_file == NULL) == (bed_file == NULL)) {
        fprintf(stderr, "ERROR: exactly one of -w/--wig or -b/--bed must be specified\n");
        usage();
        return 1;
    }
    // --maxGap / --minSize / --fast / --bin only apply to BED output (rows can
    // merge / drop / chunk-iterate / bin); wig is per-position, no "row" concept.
    if (wig_file != NULL && (max_gap != 0 || min_size != 1 || fast_mode || bin_size > 0)) {
        fprintf(stderr, "ERROR: --maxGap / --minSize / --fast / --bin are BED-mode only "
                        "(use with -b/--bed, not -w/--wig)\n");
        return 1;
    }
    if (bin_size > 0 && !fast_mode) {
        fprintf(stderr, "ERROR: --bin requires --fast\n");
        return 1;
    }
    if (bin_size > 0 && (max_gap != 0 || min_size != 1)) {
        fprintf(stderr, "ERROR: --bin is mutually exclusive with --maxGap / --minSize "
                        "(bin output is fixed-width bedGraph; merge / drop don't apply)\n");
        return 1;
    }
    if (chain_filter > 0 && !fast_mode) {
        fprintf(stderr, "ERROR: --chainFilter requires --fast\n");
        return 1;
    }
    if (chain_filter > 0 && bin_size > 0) {
        fprintf(stderr, "ERROR: --chainFilter is mutually exclusive with --bin "
                        "(bin output isn't per-block, so chain filtering doesn't apply)\n");
        return 1;
    }
    if (chain_filter > 0 && wig_file != NULL) {
        fprintf(stderr, "ERROR: --chainFilter is BED-mode only\n");
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

    // BED mode: separate code path that emits target-coord intervals in
    // source-walk order (no global sort, no spill machinery -- typical
    // BED inputs are intervals, not per-base, so output volume is
    // manageable in O(open paralog count) memory).
    if (bed_file != NULL) {
        int rc = fast_mode
            ? bed_lift_chunk_impl(tui, gl, bed_file, output_file, max_gap, min_size, bin_size, chain_filter)
            : bed_lift_main_impl (tui, gl, bed_file, output_file, max_gap, min_size);
        tui_genome_lift_destruct(gl);
        tui_destruct(tui);
        free(tui_p);
        st_logInfo("Total wall: %" PRIi64 " s\n", (int64_t)(time(NULL) - startTime));
        return rc;
    }

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
