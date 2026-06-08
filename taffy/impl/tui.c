/*
 * tui: build a Universal Column Index (.tui ONEcode container) for a
 * `cactus-hal2maf --universal` MAF/TAF.  See taffy/inc/tui.h.
 *
 *  Released under the MIT license, see LICENSE.txt
 */

#include "taf.h"
#include "tai.h"
#include "tui.h"
#include "sonLib.h"
#include "ONElib.h"
#include <errno.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <signal.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <zlib.h>

// ONEcode schema. Exactly one list field per line (STRING / INT_LIST); other
// fields are scalar INT.  The directory keys on the full genome.sequence name
// (the writer resolves the "genome" prefix via genome_of() only to group
// per-genome spills; the reader does no genome resolution at all).
//
// R holds one sequence's merged runs as a per-sequence structure-of-arrays
// delta+varint blob, zlib-deflated (see encode_runs).  The scalar INT is the
// inflated blob length (sizes the inflate buffer); the STRING list is the
// deflated bytes.  Absolute (t,g,len) triples defeat ONElib's Huffman; the
// SoA gap|gsk|lenc form is ~5x smaller than absolute end to end.
//
// X is the universal-column -> file-offset anchor index (the .tai-equivalent);
// strictly increasing, sampled every TUI_IDX_BLOCK universal columns at
// coordinate-bearing block starts.  Stored like R as a SoA delta+varint blob.
// The query-time extractor binary-searches it to seek the underlying stream
// to the nearest anchor <= the queried universal column.
static const char *TUI_SCHEMA =
    "P 3 tui\n"
    "D t 1 3 INT\n"                          // total columns T (global)
    "D X 3 3 INT 3 INT 6 STRING\n"           // univ-col index: inflatedLen, nRec, deflate(SoA)
    "O d 3 6 STRING 3 INT 3 INT\n"           // dir: seqName, S-ordinal, seqLen -- O so oneGoto works
    "O S 2 6 STRING 3 INT\n"                 // sequence object: seqName, seqLen
    "O C 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT\n"  // chunk header:
                                                   //   field 0 g_min, 1 g_max,
                                                   //   2 parent S-ord (1-based),
                                                   //   3 self c_ord  (1-based),
                                                   //   4 t_min,  5 t_max
                                                   // t_min/t_max are the per-source-seq
                                                   // coordinate bounds covered by this
                                                   // chunk's runs -- the reader uses
                                                   // them to skip chunks that don't
                                                   // overlap a source-side query
                                                   // range, avoiding the zlib decode
                                                   // of irrelevant R blobs.  Old .tui
                                                   // files written before this field
                                                   // existed only have 4 INTs on C;
                                                   // the reader checks the schema's
                                                   // field count at open and falls
                                                   // back to "decode every chunk" in
                                                   // that case.
    "D R 2 3 INT 6 STRING\n"                 // runs (one per chunk): inflatedLen, deflate(blob)
    "O g 3 6 STRING 3 INT 3 INT\n";          // genome roster (per resolved genome):
                                             //   name, total_bp (sum of seqLens),
                                             //   n_chroms.  Written at the end of
                                             //   the .tui after all per-seq d/S/C/R
                                             //   are emitted.  Old .tui (pre-roster)
                                             //   files lack this line type entirely;
                                             //   tui_genome_names() detects the
                                             //   absence via of->info['g']==NULL and
                                             //   returns NULL so callers can fall
                                             //   back to a heuristic.

// Writer-side chunk caps.  A chunk closes when EITHER trigger fires:
//
//   * TUI_CHUNK_RUNS  -- max runs per chunk.  Picked to keep each R blob a
//                        few MB compressed (a typical run encodes to ~25 B
//                        after delta+varint+deflate at 64k runs/chunk ≈
//                        ~1.6 MB).  Dense targets close on this trigger,
//                        giving ~250 chunks per gorilla chr1 -- the lazy-
//                        load granularity for narrow queries.
//
//   * TUI_CHUNK_G_MAX -- max universal-column span per chunk.  Sparse /
//                        divergent targets (e.g. tarpon vs mexican tetra)
//                        produce 65k runs scattered across millions of
//                        universal columns; without this cap each chunk's
//                        g-range is huge, the reader's outer walk-back in
//                        tui_genome_lift_column never early-exits, and a
//                        whole-chrom lift walks ~138 of 318 chunks per
//                        column.  Closing at a g-span cap keeps chunks
//                        TIGHT regardless of run density, so the early-
//                        exit fires after 1-3 steps as designed.  Dense
//                        targets hit TUI_CHUNK_RUNS first and never see
//                        this cap.  Old .tui files (without this cap)
//                        keep working at their current speed -- reader
//                        code is unchanged.
//
// At vertebrate scale on a 49 Mb tetra chrom -> tarpon lift, this dropped
// wall from 196 s to <expected 5-10 s; on tetra-close targets unchanged.
#define TUI_CHUNK_RUNS  65536
#define TUI_CHUNK_G_MAX 1000000      // 1 M universal columns

char *tui_path(const char *maf_path) {
    char *p = st_malloc(strlen(maf_path) + 5);
    sprintf(p, "%s.tui", maf_path);
    return p;
}

/////////////////////////////////////////////////////////////////////////////
// Crash/abort cleanup of phase-1 spill files.
//
// tui_cleanup() removes spills on the normal return path, but mid-phase-1
// st_errAbort() (unresolvable genome name, disk-full write, '-' strand
// row-0, etc.), assert() failures, and user SIGINT all bypass it -- and
// these spills are multi-GB at vertebrate scale.  We track every spill
// path in a file-scope list as it's created (`spill_for`, `idx_anchor`),
// and we register an atexit + signal handler that walks the list and
// remove()s each.  Also tracks the .tui output once oneFileOpenWrite has
// touched it, so a half-written .tui doesn't masquerade as a real index.
//
// On the normal exit path tui_atexit_disarm() is called before tui_cleanup()
// frees the path strings -- then the handler is a no-op.
//
// Single-threaded only.  taffy index runs tui_create once per process; this
// is not safe to call from multiple threads concurrently.
/////////////////////////////////////////////////////////////////////////////
static stList *atexit_spill_paths = NULL;  // owned strdup'd path strings
static char *atexit_tui_path = NULL;        // owned strdup'd .tui output path

// Normal-exit atexit handler: walks the spill list, unlinks each file,
// frees memory.  Runs in process context so stList_destruct (which calls
// free) is safe.  Called explicitly on st_errAbort and other graceful
// failure paths, and registered via atexit() as a backstop.
static void tui_atexit_cleanup(void) {
    if (atexit_spill_paths != NULL) {
        int64_t n = stList_length(atexit_spill_paths);
        for (int64_t i = 0; i < n; i++) {
            (void) unlink((char *) stList_get(atexit_spill_paths, i));
        }
        stList_destruct(atexit_spill_paths);
        atexit_spill_paths = NULL;
    }
    if (atexit_tui_path != NULL) {
        (void) unlink(atexit_tui_path);
        free(atexit_tui_path);
        atexit_tui_path = NULL;
    }
}

// Signal-path cleanup: ONLY async-signal-safe primitives.  unlink() and
// _exit() are on POSIX's safe list; free() and remove() are NOT (free()
// can deadlock if the signal arrived mid-malloc; remove() isn't on the
// async-safe list on every platform).  We READ but never MUTATE the
// stList -- stList_get is array-indexed access with no allocation in
// the hot path, so it's safe enough in practice.  We don't free path
// strings here; the _exit() below takes the process down before any
// leak matters.  _exit(128+sig) gives the standard "killed by signal"
// exit code and skips atexit handlers and stdio flushing (also unsafe
// in a signal context).  Mirrors the taf_lift.c lift_signal_handler
// pattern.
static void tui_signal_handler(int sig) {
    if (atexit_spill_paths != NULL) {
        int64_t n = stList_length(atexit_spill_paths);
        for (int64_t i = 0; i < n; i++) {
            (void) unlink((char *) stList_get(atexit_spill_paths, i));
        }
    }
    if (atexit_tui_path != NULL) {
        (void) unlink(atexit_tui_path);
    }
    _exit(128 + sig);
}

// Add `path` to the cleanup-on-abort list.  On first call, registers the
// atexit handler and a few signal handlers.  Caller retains ownership of the
// passed string; we strdup.
static void tui_atexit_track_spill(const char *path) {
    if (atexit_spill_paths == NULL) {
        atexit_spill_paths = stList_construct3(0, free);
        atexit(tui_atexit_cleanup);
        signal(SIGINT,  tui_signal_handler);
        signal(SIGTERM, tui_signal_handler);
        signal(SIGHUP,  tui_signal_handler);
    }
    stList_append(atexit_spill_paths, stString_copy(path));
}

// Track the .tui output path so a half-written file is removed on crash.
// Called after oneFileOpenWriteNew() succeeds; cleared by tui_atexit_disarm().
// Exported (extern in tui.h) so secondary writers (taf_coarsen) can use the
// same crash-cleanup safety net.
void tui_atexit_track_tui(const char *path) {
    if (atexit_tui_path != NULL) { free(atexit_tui_path); }
    atexit_tui_path = stString_copy(path);
}

// Disarm the crash cleanup on the normal success path, BEFORE tui_cleanup()
// frees its path strings.  Idempotent.  Exported for taf_coarsen.
void tui_atexit_disarm(void) {
    if (atexit_spill_paths != NULL) {
        stList_destruct(atexit_spill_paths);
        atexit_spill_paths = NULL;
    }
    if (atexit_tui_path != NULL) {
        free(atexit_tui_path);
        atexit_tui_path = NULL;
    }
}

// directory component of a path (malloc'd); "." if none. portable (no libgen)
static char *dir_of(const char *path) {
    const char *slash = strrchr(path, '/');
    if (slash == NULL) return stString_copy(".");
    if (slash == path) return stString_copy("/");
    char *d = st_malloc(slash - path + 1);
    memcpy(d, path, slash - path);
    d[slash - path] = '\0';
    return d;
}

// basename component of a path (pointer into path; no alloc)
static const char *base_of(const char *path) {
    const char *slash = strrchr(path, '/');
    return slash ? slash + 1 : path;
}

// Resolve the genome name of a full "genome.sequence" name (caller frees).
//
// <=1 '.'  : unambiguous -> split at the (only) '.', or the whole name if none.
//            (covers ancestors like `MuridaeAnc3.MuridaeAnc3refChr0` and any
//            plain `genome.chr`); no name map needed even if one is supplied.
// >1  '.'  : ambiguous (e.g. accessions `GCF_947179515.1.NC_067472.1`).  A
//            genome_name_map is REQUIRED; resolution reuses taf.c's
//            extract_genome_name (walks each '.' offset, returns the first
//            prefix that is a key of the map).  Missing map or no match is a
//            fatal error -- silently mis-splitting a coordinate index is worse.
static char *genome_of(const char *seq_name, stHash *gmap) {
    int ndots = 0;
    for (const char *s = seq_name; (s = strchr(s, '.')) != NULL; s++) ndots++;
    if (ndots <= 1) {
        const char *dot = strchr(seq_name, '.');
        int64_t n = dot ? (dot - seq_name) : (int64_t)strlen(seq_name);
        char *g = st_malloc(n + 1);
        memcpy(g, seq_name, n);
        g[n] = '\0';
        return g;
    }
    if (gmap == NULL) {
        st_errAbort("tui: sequence name '%s' has more than one '.' so its "
                    "genome/sequence split is ambiguous; supply a genome name "
                    "list with `taffy index -n <file>`", seq_name);
    }
    char *g = extract_genome_name(seq_name, NULL, gmap);  // reuse taf.c parser
    if (g == NULL) {
        st_errAbort("tui: could not resolve a genome for '%s' using the "
                    "provided name map (no listed genome is a prefix of it)",
                    seq_name);
    }
    return g;
}

// Build a genome-name membership set (keys, dummy values) from a Newick string
// -- every node label, leaf and internal (e.g. hal2maf's "# hal" tree gives the
// full genome set including ancestors).  A label is a maximal run of non-
// "(),:;"/non-space chars that is NOT a branch length (i.e. not the run right
// after a ':').  Same stHash shape extract_genome_name / genome_of expect.
// Grammar is the plain hal2maf form: NO quoted labels and NO NHX comments
// ([&&NHX...] / '...' would be mis-tokenized); labels must be < sizeof(tok)
// (over-long is a hard error, not a silent truncate).
static stHash *genome_set_from_newick(const char *nwk) {
    stHash *s = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    char tok[8192];
    int n = 0, after_colon = 0;
    for (const char *p = nwk; ; p++) {
        char c = *p;
        int delim = (c == '\0' || c == '(' || c == ')' || c == ',' ||
                     c == ':' || c == ';' || c == ' ' || c == '\t' ||
                     c == '\n' || c == '\r');
        if (!delim) {
            if (n >= (int)sizeof(tok) - 1) {
                st_errAbort("tui: # hal tree label exceeds %d chars; cannot "
                            "resolve genome names from it", (int)sizeof(tok) - 1);
            }
            tok[n++] = c;
        } else {
            if (n > 0) {
                tok[n] = '\0';
                if (!after_colon && stHash_search(s, tok) == NULL) {
                    stHash_insert(s, stString_copy(tok), (void *)1);
                }
                n = 0;
            }
            if (c == ':') after_colon = 1;          // next run is a branch length
            else if (c != '\0') after_colon = 0;    // '(' ',' ')' ';' end it
            if (c == '\0') break;
        }
    }
    return s;
}

// One distinct "genome.sequence" name plus its gmap-resolved genome.  The
// genome is computed once (genome_of) and cached here so Phase 2 never
// re-resolves per sequence, and -- crucially -- so the sort key is the TRUE
// genome, not a first-dot guess.  That makes "each genome's sequences form one
// contiguous block in seqs[]" true BY CONSTRUCTION (was only an assumption
// when ordering split at the first '.').
typedef struct { char *seq; char *genome; } SeqKey;

// THE global order: genome-major (true resolved genome), then full sequence
// name.  seqs[] uses this; runs[] (run_cmp) and the Phase-2 forward cursor use
// strcmp on the full name, which within one genome's spill is exactly this
// order's within-genome part -- so the cursor can't skip/misattribute.
static int seqkey_cmp(const void *a, const void *b) {
    const SeqKey *x = a, *y = b;
    int c = strcmp(x->genome, y->genome);
    return c ? c : strcmp(x->seq, y->seq);
}

/////////////////////////////////////////////////////////////////////////////
// Varint helpers (LEB128) -- used by both the spill format (below) and the
// per-sequence run codec further down.
/////////////////////////////////////////////////////////////////////////////

static inline uint64_t zigzag(int64_t v) { return ((uint64_t)v << 1) ^ (uint64_t)(v >> 63); }
static inline int64_t unzigzag(uint64_t z) { return (int64_t)(z >> 1) ^ -(int64_t)(z & 1); }

static inline size_t put_uvarint(uint8_t *p, uint64_t v) {
    size_t n = 0;
    while (v >= 0x80) { p[n++] = (uint8_t)(v | 0x80); v >>= 7; }
    p[n++] = (uint8_t)v;
    return n;
}
static inline uint64_t get_uvarint(const uint8_t **pp) {
    const uint8_t *p = *pp;
    uint64_t v = 0; int s = 0; uint8_t b;
    do { b = *p++; v |= (uint64_t)(b & 0x7f) << s; s += 7; } while (b & 0x80);
    *pp = p;
    return v;
}

// Spill file binary record format (one stream of records, no header):
//
//   tag byte:
//     'N' (0x4E) -- new-name: tag, uvarint seq_idx, uvarint name_len, name_bytes
//                   Emitted ONCE the first time a (seq_name) is seen in this
//                   spill.  Subsequent data records reference the seq_idx.
//     'D' (0x44) -- data: tag, uvarint seq_idx, uvarint t_start, uvarint g_start,
//                   uvarint length, byte strand ('+' or '-')
//
// Replaces the prior text format ("%s\t%lld\t%lld\t%lld\t%c\n") which averaged
// ~70-100 bytes/record on vertebrate data (seq_name 15-30 chars + 3 decimal
// int64s + separators).  Binary is ~10-15 bytes/record; observed amplification
// vs decompressed MAF drops from ~80x to ~10x.  No correctness change; on-disk
// only.  The Run struct returned by load_genome_runs is unchanged so callers
// don't need updates.

// Write one spilled run.  `seen_seqs` is a per-spill stHash<seq_name, int64_t*>
// owned by the caller (Phase1::spill_seqs); we lazily insert seq_idx on first
// sight and write the 'N' dictionary record before the 'D' data record.
static void spill_run(FILE *fh, stHash *seen_seqs, const char *seq_name,
                      int64_t t_start, int64_t g_start, int64_t length, bool strand) {
    int64_t *p_idx = stHash_search(seen_seqs, (void *)seq_name);
    int64_t seq_idx;
    if (p_idx == NULL) {
        seq_idx = (int64_t)stHash_size(seen_seqs);
        int64_t *new_idx = st_malloc(sizeof(int64_t));
        *new_idx = seq_idx;
        stHash_insert(seen_seqs, stString_copy(seq_name), new_idx);

        // Emit 'N' record: tag + varint seq_idx + varint name_len + name bytes.
        size_t nl = strlen(seq_name);
        uint8_t hdr[1 + 10 + 10];  // tag + 2 varints (max 10 bytes each)
        size_t hn = 0;
        hdr[hn++] = 'N';
        hn += put_uvarint(hdr + hn, (uint64_t)seq_idx);
        hn += put_uvarint(hdr + hn, (uint64_t)nl);
        if (fwrite(hdr, 1, hn, fh) != hn ||
            (nl > 0 && fwrite(seq_name, 1, nl, fh) != nl)) {
            st_errAbort("tui: failed writing N record to spill (disk full / write error)");
        }
    } else {
        seq_idx = *p_idx;
    }

    // Emit 'D' record: tag + 4 varints + strand byte.  Caps at 1 + 4*10 + 1.
    uint8_t buf[1 + 4 * 10 + 1];
    size_t n = 0;
    buf[n++] = 'D';
    n += put_uvarint(buf + n, (uint64_t)seq_idx);
    n += put_uvarint(buf + n, (uint64_t)t_start);
    n += put_uvarint(buf + n, (uint64_t)g_start);
    n += put_uvarint(buf + n, (uint64_t)length);
    buf[n++] = strand ? '+' : '-';
    if (fwrite(buf, 1, n, fh) != n) {
        st_errAbort("tui: failed writing D record to spill (disk full / write error)");
    }
}

// Emit the maximal gap-free stretches of one row of one block.  `g` is the
// global column of this block's first column.
static void emit_row_runs(FILE *spill, stHash *seen_seqs, Alignment_Row *row,
                          int64_t g, int64_t column_number) {
    const char *b = row->bases;
    int64_t pre = 0;            // non-gap bases of this row before column c
    int64_t c = 0;
    while (c < column_number) {
        if (b[c] == '-') { c++; continue; }
        int64_t c0 = c;
        while (c < column_number && b[c] != '-') c++;
        int64_t slen = c - c0;                       // gap-free stretch [c0, c)
        int64_t t_start, g_start = g + c0;
        if (row->strand) {                           // '+'
            t_start = row->start + pre;
        } else {                                     // '-' : MAF-native start
            int64_t fwd_end = row->sequence_length - row->start;
            t_start = fwd_end - pre - slen;
        }
        spill_run(spill, seen_seqs, row->sequence_name, t_start, g_start, slen, row->strand);
        pre += slen;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Phase 1: one streaming scan, spill runs per-genome (column order, O(1) RAM)
/////////////////////////////////////////////////////////////////////////////

// Per-genome spill state.  At vertebrate scale (1k+ genomes in the # hal
// tree) we can't keep one FILE* per genome open for all of phase 1 -- Linux
// soft RLIMIT_NOFILE defaults to 1024 and we'd hit EMFILE.  Instead phase 1
// holds a bounded LRU pool of currently-open FILE*s (max_open below) and
// re-opens evicted spills in append mode on next write.  Each genome's
// SpillEnt is created on first sighting (next_spill_id++ assigns its file
// name) and lives until phase-1 cleanup.
typedef struct SpillEnt {
    char  *path;          // we own; created once on first sighting of this genome
    FILE  *fh;            // NULL when evicted from the open-FH pool; reopened on next write
    struct SpillEnt *prev_lru, *next_lru;  // LRU list (MRU=head, LRU=tail)
} SpillEnt;

typedef struct {
    stHash *spill_ents;   // genomeName -> SpillEnt*  (lifetime of phase 1)
    SpillEnt *lru_head;   // MRU end (most recently opened/touched)
    SpillEnt *lru_tail;   // LRU end (next eviction candidate)
    int     n_open;       // count of SpillEnt with fh != NULL
    int     max_open;     // cap; eviction starts when n_open == max_open
    stHash *spill_seqs;   // genomeName -> stHash<seq_name, int64_t*>  (per-spill
                          //                seq-name dictionary; written inline
                          //                as 'N' records, used by reader to
                          //                rehydrate seq names)
    stHash *seq_len;      // "genome.seq" -> int64_t* seqLen
    stList *seq_keys;     // distinct "genome.seq", first-seen order
    int64_t T;            // total columns
    const char *tmp_dir;  // where spill files go
    const char *out_base; // basename of out path (for unique spill names)
    int next_spill_id;
    stHash *gmap;         // optional name map for >1-dot genome resolution
    // Universal-column -> file_pos index (the .tai-equivalent, monotone key).
    // Anchors sampled every TUI_IDX_BLOCK columns at coordinate-bearing block
    // starts (MAF: every block).  Spilled (col<TAB>file_pos), column-ordered.
    FILE   *idx_fh;
    char   *idx_path;
    int64_t idx_last_col;
    int64_t idx_n;        // anchors written so far (0 => none yet)
} Phase1;

#define TUI_IDX_BLOCK 10000   // universal columns between index anchors (== tai default)

// Record a (universal column -> file_pos) anchor for the block that the NEXT
// maf/taf_read_block will return.  `col` = p1->T at the block start (universal
// column of its first column); `file_pos` = LI_tell(li) BEFORE the read
// (exactly the tai_create_maf/taf convention).  TAF: only anchor on a
// coordinate-bearing line (so tai_resync_taf_line can resync there); MAF:
// every block start is a valid anchor.  Always anchor the first block.
static void idx_anchor(Phase1 *p1, LI *li, int is_maf, bool rle) {
    int64_t col = p1->T;
    if (p1->idx_n != 0 && col - p1->idx_last_col < TUI_IDX_BLOCK) return;
    char *peek = LI_peek_at_next_line(li);
    if (peek == NULL) return;
    // anchor = start of the block's first line (the peeked line): exactly the
    // tai_create_taf convention so LI_seek + tai_resync_taf_line work.
    int64_t file_pos = LI_get_position(li);
    if (!is_maf) {              // TAF: only a full-coordinate (resync-able) line
        stList *toks = stString_split(peek);
        int ok = tai_taf_line_is_anchor(toks, rle);
        stList_destruct(toks);
        if (!ok) return;
    }
    if (p1->idx_fh == NULL) {
        p1->idx_path = stString_print("%s/%s.tuiIdx.%ld", p1->tmp_dir,
                                      p1->out_base, (long)getpid());
        tui_atexit_track_spill(p1->idx_path);
        p1->idx_fh = fopen(p1->idx_path, "w");
        if (p1->idx_fh == NULL) {
            fprintf(stderr, "tui: cannot open idx spill %s\n", p1->idx_path);
            exit(1);
        }
    }
    if (fprintf(p1->idx_fh, "%" PRIi64 "\t%" PRIi64 "\n", col, file_pos) < 0) {
        st_errAbort("tui: failed writing idx spill (disk full / write error)");
    }
    p1->idx_last_col = col;
    p1->idx_n++;
}

// Assert the row-0 invariants the universal lift relies on: '+' strand and
// gap-free (row-0 base count == column count).  Both are guaranteed by
// `cactus-hal2maf --universal` (maxRefGap=0); a violation here is bad input.
static void assert_row0_universal(const Alignment_Row *row0, int64_t cn) {
    if (!row0->strand) {
        st_errAbort("tui: row-0 sequence '%s' is on '-' strand; a universal "
                    "MAF row-0 must be '+'", row0->sequence_name);
    }
    if (row0->length != cn) {
        st_errAbort("tui: row-0 '%s' is not gap-free (length %" PRIi64
                    " != column_number %" PRIi64 "); --universal / "
                    "maxRefGap==0 expected", row0->sequence_name,
                    row0->length, cn);
    }
}

// LRU helpers: unlink/insert at head.  Caller invariant: `e` is not on
// the list when lru_push_head is called.
static void lru_unlink(Phase1 *p1, SpillEnt *e) {
    if (e->prev_lru) e->prev_lru->next_lru = e->next_lru;
    else             p1->lru_head           = e->next_lru;
    if (e->next_lru) e->next_lru->prev_lru = e->prev_lru;
    else             p1->lru_tail           = e->prev_lru;
    e->prev_lru = e->next_lru = NULL;
}
static void lru_push_head(Phase1 *p1, SpillEnt *e) {
    e->prev_lru = NULL;
    e->next_lru = p1->lru_head;
    if (p1->lru_head) p1->lru_head->prev_lru = e;
    p1->lru_head = e;
    if (p1->lru_tail == NULL) p1->lru_tail = e;
}

// Evict the LRU spill (fclose its FH; ent stays in the hash so its path is
// remembered for reopen).  Never evicts `pinned` (the ent we're about to
// open / write to).  Returns 1 on success.
//
// Invariant: every ent in the LRU list has fh != NULL; spill_for only calls
// this in the cold branch where pinned->fh == NULL, so pinned is by
// construction not on the list and the `v == pinned` defense is dead code.
// Keep it -- correctness invariants drift, and the cost is one pointer
// compare per eviction.
static int spill_evict_one(Phase1 *p1, SpillEnt *pinned) {
    SpillEnt *v = p1->lru_tail;
    if (v == pinned) v = v->prev_lru;
    if (v == NULL) return 0;
    if (fclose(v->fh) != 0) {
        st_errAbort("tui: spill close (LRU eviction) failed: %s -- %s",
                    v->path, strerror(errno));
    }
    v->fh = NULL;
    lru_unlink(p1, v);
    p1->n_open--;
    return 1;
}

static FILE *spill_for(Phase1 *p1, const char *genome) {
    SpillEnt *e = stHash_search(p1->spill_ents, (void *)genome);
    if (e == NULL) {
        // first sighting of this genome: allocate the ent and pick a path.
        // File creation itself is deferred to the fopen below (still on this
        // call, but conceptually separate so the EMFILE-eviction loop applies
        // uniformly to first-open and reopen cases).
        e = st_calloc(1, sizeof(SpillEnt));
        e->path = stString_print("%s/%s.tuiSpill.%ld.%d", p1->tmp_dir,
                                  p1->out_base, (long)getpid(),
                                  p1->next_spill_id++);
        tui_atexit_track_spill(e->path);
        stHash_insert(p1->spill_ents, stString_copy(genome), e);
        // Per-spill seq-name dictionary; values are heap int64_t* seq indices.
        stHash *seen = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
                                         free, free);
        stHash_insert(p1->spill_seqs, stString_copy(genome), seen);
    }
    if (e->fh == NULL) {
        // Cold ent: need to (re)open.  Make room in the FH pool first.
        while (p1->n_open >= p1->max_open) {
            if (!spill_evict_one(p1, e)) break;   // pool empty besides `e`
        }
        // "ab" creates on first open and appends on reopen (no truncation).
        e->fh = fopen(e->path, "ab");
        if (e->fh == NULL) {
            st_errAbort("tui: cannot open spill %s -- %s (n_open=%d max_open=%d)",
                        e->path, strerror(errno), p1->n_open, p1->max_open);
        }
        p1->n_open++;
        lru_push_head(p1, e);
    } else if (p1->lru_head != e) {
        // Warm ent but not MRU: bump to head so it survives the next eviction.
        lru_unlink(p1, e);
        lru_push_head(p1, e);
    }
    return e->fh;
}

static void note_seq(Phase1 *p1, const char *seq_name, int64_t seq_len) {
    if (stHash_search(p1->seq_len, (void *)seq_name) == NULL) {
        int64_t *v = st_malloc(sizeof(int64_t));
        *v = seq_len;
        char *k = stString_copy(seq_name);
        stHash_insert(p1->seq_len, k, v);
        stList_append(p1->seq_keys, k);   // shares the hash's key allocation
    }
}

static void tui_phase1_block(Phase1 *p1, Alignment *aln) {
    int64_t cn = aln->column_number;
    if (aln->row == NULL) { p1->T += cn; return; }  // degenerate a-line, no rows
    assert_row0_universal(aln->row, cn);
    for (Alignment_Row *row = aln->row; row != NULL; row = row->n_row) {
        char *gname = genome_of(row->sequence_name, p1->gmap);
        note_seq(p1, row->sequence_name, row->sequence_length);
        FILE *spill = spill_for(p1, gname);
        stHash *seen = stHash_search(p1->spill_seqs, (void *)gname);
        emit_row_runs(spill, seen, row, p1->T, cn);
        free(gname);
    }
    p1->T += cn;
}

/////////////////////////////////////////////////////////////////////////////
// Phase 2: per-genome in-RAM sort by (seq, t_start), write ONEcode container
//
// Memory: bounded PER GENOME by that genome's spill size -- not absolutely.
// At vertebrate scale the biggest leaf (a Mus chr1 / sleep mr-style spill)
// can carry ~30M runs at ~40 B/Run = ~1.2 GB peak.  Acceptable on a build
// host with the kind of RAM budget --universal already demands (the writer
// holds Phase 1 spill FHs and Phase 2 sort arenas together, never both at
// the same scale).
/////////////////////////////////////////////////////////////////////////////

typedef struct { char *seq; int64_t t, g, len; char strand; } Run;

// Directory-line write key: links a name-sorted output position back to the
// genome-major S-ordinal of the sequence (so the reader can binary-search the
// d-lines by name and `oneGoto S` straight to the matching object).
typedef struct { int64_t idx; const char *name; } DKey;
static int dkey_cmp(const void *a, const void *b) {
    return strcmp(((const DKey *)a)->name, ((const DKey *)b)->name);
}

// runs[] hold one genome's runs (one spill), so full-name strcmp is exactly
// the within-genome order seqs[] uses (genome prefix is constant here).
static int run_cmp(const void *a, const void *b) {
    const Run *x = a, *y = b;
    int s = strcmp(x->seq, y->seq);
    if (s != 0) return s;
    return (x->t < y->t) ? -1 : (x->t > y->t) ? 1 : 0;
}

// Compare two (t, g, lenc) triples by the g field (the universal column).
// Used by the writer to reorder a sequence's runs from t-sorted (the
// merge-loop's natural order) to g-sorted (so per-chunk g ranges are
// tight and the reader can filter chunks effectively).
static int triple_cmp_by_g(const void *a, const void *b) {
    const int64_t *x = a, *y = b;
    return (x[1] < y[1]) ? -1 : (x[1] > y[1]) ? 1 : 0;
}

// Compare two (t, g, lenc) triples by the t field.  Used by
// tui_load_seq_runs / tui_query to restore the t-sorted view callers
// (source-coord lookups in taffy lift) expect, after concatenating
// per-chunk g-sorted blobs.
static int triple_cmp_by_t(const void *a, const void *b) {
    const int64_t *x = a, *y = b;
    return (x[0] < y[0]) ? -1 : (x[0] > y[0]) ? 1 : 0;
}

// Load every run record of one genome spill into a Run array (caller frees
// each runs[k].seq + the runs[] array itself).  See spill_run for the on-disk
// format: alternating 'N' (seq dictionary) and 'D' (data) records.
//
// We mmap the spill; the per-spill seq dictionary is built locally as we walk
// 'N' records.  Each 'D' record's seq pointer borrows directly from the dict
// (NOT stString_copy'd per run) -- the dict is returned to the caller so it
// outlives runs[] and gets freed via *dict_out instead of per-run.
//
// Why interning matters: pre-fix this function stString_copy'd dict[idx]
// PER RUN.  For a giant-ancestor genome (~1B runs at 577-way root scale)
// that's ~30 GB of redundant string copies per worker.  Combined with
// n_threads concurrent workers in phase 2, this is what was OOMing the
// 1.5 TB index host on the 577-way at high thread counts.
//
// Caller contract: free dict_out[0..dict_cap_out) then free(dict_out)
// AFTER you're done with runs[].  Order matters -- runs[i].seq is a
// borrowed pointer into dict_out[i], so freeing dict first would leave
// runs holding dangling pointers.
static Run *load_genome_runs(const char *path, int64_t *n_out,
                             char ***dict_out, int64_t *dict_cap_out) {
    int fd = open(path, O_RDONLY);
    if (fd < 0) { *n_out = 0; return NULL; }
    struct stat sb;
    if (fstat(fd, &sb) < 0) { close(fd); *n_out = 0; return NULL; }
    if (sb.st_size == 0) { close(fd); *n_out = 0; return NULL; }
    uint8_t *base = mmap(NULL, (size_t)sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (base == MAP_FAILED) {
        st_errAbort("tui: mmap of spill %s failed", path);
    }

    // Per-spill local dictionary: seq_idx -> seq_name (caller-freed).  Sized
    // to seq_idx as 'N' records grow it.
    int64_t dict_cap = 1024;
    char **dict = st_malloc(dict_cap * sizeof(char *));
    for (int64_t i = 0; i < dict_cap; i++) dict[i] = NULL;

    int64_t cap = 1024, n = 0;
    Run *runs = st_malloc(cap * sizeof(Run));

    const uint8_t *p = base, *end = base + sb.st_size;
    while (p < end) {
        uint8_t tag = *p++;
        if (tag == 'N') {
            uint64_t idx = get_uvarint(&p);
            uint64_t nl  = get_uvarint(&p);
            if ((int64_t)idx >= dict_cap) {
                int64_t old = dict_cap;
                while (dict_cap <= (int64_t)idx) dict_cap *= 2;
                dict = st_realloc(dict, dict_cap * sizeof(char *));
                for (int64_t i = old; i < dict_cap; i++) dict[i] = NULL;
            }
            char *nm = st_malloc((size_t)nl + 1);
            memcpy(nm, p, (size_t)nl);
            nm[nl] = '\0';
            dict[(int64_t)idx] = nm;
            p += nl;
        } else if (tag == 'D') {
            if (n == cap) { cap *= 2; runs = st_realloc(runs, cap * sizeof(Run)); }
            uint64_t idx = get_uvarint(&p);
            uint64_t t   = get_uvarint(&p);
            uint64_t g   = get_uvarint(&p);
            uint64_t len = get_uvarint(&p);
            char st = (char)*p++;
            if ((int64_t)idx >= dict_cap || dict[(int64_t)idx] == NULL) {
                st_errAbort("tui: corrupt spill %s: D-record references "
                            "undefined seq_idx %" PRIi64, path, (int64_t)idx);
            }
            Run *r = &runs[n++];
            r->seq    = dict[(int64_t)idx];   // borrowed; dict outlives runs
            r->t      = (int64_t)t;
            r->g      = (int64_t)g;
            r->len    = (int64_t)len;
            r->strand = st;
        } else {
            st_errAbort("tui: corrupt spill %s: unknown record tag 0x%02x "
                        "at offset %" PRIi64, path, (int)tag,
                        (int64_t)(p - 1 - base));
        }
    }

    // Hand the dict off to the caller -- runs[i].seq is a borrowed pointer
    // into it, so it has to outlive the runs[] sort + per-seq processing.
    munmap(base, (size_t)sb.st_size);

    *n_out        = n;
    *dict_out     = dict;
    *dict_cap_out = dict_cap;
    return runs;
}

/////////////////////////////////////////////////////////////////////////////
// Per-sequence run codec: zigzag-delta + LEB128 varint, then zlib deflate.
/////////////////////////////////////////////////////////////////////////////

// (defined earlier in the file -- see "Varint helpers" section near the top)

// Encode m runs (buf = m triples [t, g, lenc=(len<<1|rev)], absolute) and
// zlib-deflate them.  Caller frees the returned buffer; *raw_len = inflated
// length, *def_len = deflated length.
//
// Structure-of-arrays: three concatenated varint streams, not interleaved.
// Per run (prev_end_t = pt+pl, prev_end_g = pg+pl; pt,pg,pl carry the prior
// run's absolute t,g and length; all seeded 0 so run 0 stores absolutes):
//   gap  = t - prev_end_t   (== 0 for ~99% of splits: the sequence stays
//                             forward-contiguous, only OTHER lineages' columns
//                             intervened -> this stream deflates to ~nothing)
//   gsk  = g - prev_end_g   (intervening universal columns; the real signal)
//   lenc = (len<<1)|strand
// The delta math never uses strand -> the codec is a strand-agnostic
// bijection; strand only rides in lenc for the query formula.  Measured ~22%
// smaller than the interleaved form (SoA lets zlib exploit each stream's very
// different statistics, esp. the all-zero gap stream).
//
// Raw layout: header [uvarint m, uvarint |gap bytes|, uvarint |gsk bytes|]
// then gap||gsk||lenc.  raw_len bounds the inflate buffer; m bounds decode.
/* Exported via tui.h (internal-but-shared with taf_coarsen). */
uint8_t *tui_encode_runs(const int64_t *buf, int64_t m,
                         int64_t *raw_len, int64_t *def_len) {
    uint8_t *G = st_malloc((size_t)(m * 10 + 1));   // gap stream
    uint8_t *K = st_malloc((size_t)(m * 10 + 1));   // gsk stream
    uint8_t *L = st_malloc((size_t)(m * 10 + 1));   // lenc stream
    size_t gn = 0, kn = 0, ln = 0;
    int64_t pt = 0, pg = 0, pl = 0;
    for (int64_t k = 0; k < m; k++) {
        int64_t t = buf[3*k+0], g = buf[3*k+1], lenc = buf[3*k+2];
        gn += put_uvarint(G + gn, zigzag(t - (pt + pl)));
        kn += put_uvarint(K + kn, zigzag(g - (pg + pl)));
        ln += put_uvarint(L + ln, (uint64_t)lenc);
        pt = t; pg = g; pl = lenc >> 1;
    }
    uint8_t hdr[30];
    size_t hn = 0;
    hn += put_uvarint(hdr + hn, (uint64_t)m);
    hn += put_uvarint(hdr + hn, (uint64_t)gn);
    hn += put_uvarint(hdr + hn, (uint64_t)kn);
    size_t rn = hn + gn + kn + ln;
    uint8_t *raw = st_malloc(rn ? rn : 1);
    memcpy(raw, hdr, hn);
    memcpy(raw + hn, G, gn);
    memcpy(raw + hn + gn, K, kn);
    memcpy(raw + hn + gn + kn, L, ln);
    free(G); free(K); free(L);

    uLongf cap = compressBound((uLong)rn);
    uint8_t *def = st_malloc(cap ? cap : 1);
    uLongf dl = cap;
    if (compress2(def, &dl, raw, (uLong)rn, 9) != Z_OK) {
        st_errAbort("tui: zlib compress2 failed");
    }
    free(raw);
    *raw_len = (int64_t)rn;
    *def_len = (int64_t)dl;
    return def;
}

// Inverse of encode_runs: inflate, split the three streams by the header, and
// rebuild the m absolute triples into out[3*m].  Returns m.
static int64_t decode_runs(const uint8_t *def, int64_t def_len,
                           int64_t raw_len, int64_t *out, int64_t out_cap) {
    uint8_t *raw = st_malloc((size_t)(raw_len ? raw_len : 1));
    uLongf rl = (uLongf)raw_len;
    if (uncompress(raw, &rl, def, (uLong)def_len) != Z_OK || (int64_t)rl != raw_len) {
        st_errAbort("tui: zlib uncompress failed");
    }
    const uint8_t *h = raw;
    int64_t m = (int64_t)get_uvarint(&h);
    int64_t gn = (int64_t)get_uvarint(&h);
    int64_t kn = (int64_t)get_uvarint(&h);
    const uint8_t *gp = h, *kp = h + gn, *lp = h + gn + kn;
    int64_t pt = 0, pg = 0, pl = 0;
    for (int64_t k = 0; k < m; k++) {
        int64_t t = pt + pl + unzigzag(get_uvarint(&gp));
        int64_t g = pg + pl + unzigzag(get_uvarint(&kp));
        int64_t lenc = (int64_t)get_uvarint(&lp);
        assert(3 * k + 2 < out_cap);
        out[3*k+0] = t; out[3*k+1] = g; out[3*k+2] = lenc;
        pt = t; pg = g; pl = lenc >> 1;
    }
    free(raw);
    return m;
}

// Universal-column index codec.  n anchors (col, file_pos), both STRICTLY
// increasing in anchor order, so plain non-negative uvarint deltas (no
// zigzag).  Two SoA streams C|F; header = uvarint(n), uvarint(|C bytes|);
// one zlib deflate.  Same shape as encode_runs but tuned for the index pairs.
/* Exported via tui.h (internal-but-shared with taf_coarsen). */
uint8_t *tui_encode_idx(const int64_t *col, const int64_t *fpos, int64_t n,
                        int64_t *raw_len, int64_t *def_len) {
    uint8_t *C = st_malloc((size_t)(n * 10 + 1));
    uint8_t *F = st_malloc((size_t)(n * 10 + 1));
    size_t cn = 0, fn = 0;
    int64_t pc = 0, pf = 0;
    for (int64_t k = 0; k < n; k++) {
        cn += put_uvarint(C + cn, (uint64_t)(col[k]  - pc));   // >= 0 (sorted)
        fn += put_uvarint(F + fn, (uint64_t)(fpos[k] - pf));
        pc = col[k]; pf = fpos[k];
    }
    uint8_t hdr[20];
    size_t hn = 0;
    hn += put_uvarint(hdr + hn, (uint64_t)n);
    hn += put_uvarint(hdr + hn, (uint64_t)cn);
    size_t rn = hn + cn + fn;
    uint8_t *raw = st_malloc(rn ? rn : 1);
    memcpy(raw, hdr, hn);
    memcpy(raw + hn, C, cn);
    memcpy(raw + hn + cn, F, fn);
    free(C); free(F);
    uLongf cap = compressBound((uLong)rn);
    uint8_t *def = st_malloc(cap ? cap : 1);
    uLongf dl = cap;
    if (compress2(def, &dl, raw, (uLong)rn, 9) != Z_OK) {
        st_errAbort("tui: zlib compress2 failed (X)");
    }
    free(raw);
    *raw_len = (int64_t)rn;
    *def_len = (int64_t)dl;
    return def;
}

static int64_t decode_idx(const uint8_t *def, int64_t def_len, int64_t raw_len,
                          int64_t *col, int64_t *fpos) {
    uint8_t *raw = st_malloc((size_t)(raw_len ? raw_len : 1));
    uLongf rl = (uLongf)raw_len;
    if (uncompress(raw, &rl, def, (uLong)def_len) != Z_OK || (int64_t)rl != raw_len) {
        st_errAbort("tui: zlib uncompress failed (X)");
    }
    const uint8_t *h = raw;
    int64_t n  = (int64_t)get_uvarint(&h);
    int64_t cb = (int64_t)get_uvarint(&h);
    const uint8_t *cp = h, *fp = h + cb;
    int64_t pc = 0, pf = 0;
    for (int64_t k = 0; k < n; k++) {
        pc += (int64_t)get_uvarint(&cp); col[k]  = pc;
        pf += (int64_t)get_uvarint(&fp); fpos[k] = pf;
    }
    free(raw);
    return n;
}

// Remove spill files and free phase-1 state.  Shared by the success and the
// error-return paths so a failed run never leaks spill files on disk.
// `tree_map` (the # hal-derived genome set, or NULL) is freed here too so its
// lifetime is tied to one place regardless of which exit path is taken.
static void tui_cleanup(Phase1 *p1, SeqKey *seqks, int64_t n_seqs,
                        char *eff_tmp, stHash *tree_map) {
    if (p1->spill_ents != NULL) {
        stHashIterator *hit = stHash_getIterator(p1->spill_ents);
        char *gk;
        while ((gk = stHash_getNext(hit)) != NULL) {
            SpillEnt *e = stHash_search(p1->spill_ents, gk);
            if (e->fh != NULL) {
                fclose(e->fh);
                e->fh = NULL;
            }
            remove(e->path);
            free(e->path);
            free(e);
        }
        stHash_destructIterator(hit);
        stHash_destruct(p1->spill_ents);
    }
    if (p1->idx_fh != NULL) fclose(p1->idx_fh);
    if (p1->idx_path != NULL) { remove(p1->idx_path); free(p1->idx_path); }
    for (int64_t i = 0; i < n_seqs; i++) free(seqks[i].genome);  // .seq owned by seq_len
    free(seqks);
    if (p1->spill_seqs != NULL) stHash_destruct(p1->spill_seqs);
    stHash_destruct(p1->seq_len);
    stList_destruct(p1->seq_keys);
    if (tree_map != NULL) stHash_destruct(tree_map);
    free(eff_tmp);
}

/////////////////////////////////////////////////////////////////////////////
// Phase 2 worker / writer split.  The per-genome body (load spill -> sort
// -> per-seq colinear merge -> chunk -> encode) is pure CPU once we have
// the genome's spill path and per-seqks slen lookup; it has no shared
// state.  The writer side (oneWriteLine S/C/R into the OneFile, plus the
// file-position c_ord_emit counter that gets stamped into each C record)
// must stay single-threaded.  We split the loop into:
//
//   phase2_genome_work(gi, ...) -> Phase2Genome*   [worker, parallel]
//   phase2_genome_write(of, g, &c_ord_emit)        [writer, serial]
//
// so an OpenMP parallel-for-ordered can run N workers concurrently while
// the ordered region serialises writes in seqks order (matching the
// pre-committed d-record S-ordinals).
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    int64_t g_min, g_max;
    int64_t t_min, t_max;
    int64_t raw_len;
    uint8_t *def;
    int64_t def_len;
    int64_t cm;              // # runs in chunk (for self-check + bookkeeping)
} Phase2Chunk;

typedef struct {
    char  *seq_name;         // borrowed from seqks (do NOT free here)
    int64_t slen;
    int64_t seqks_idx;       // 0-based position in seqks (s_ord = seqks_idx + 1)
    Phase2Chunk *chunks;     // length n_chunks
    int64_t n_chunks;
} Phase2Seq;

typedef struct {
    Phase2Seq *seqs;         // length n_seqs
    int64_t n_seqs;
    int64_t nr;              // total runs loaded from this genome's spill (for log)
} Phase2Genome;

// Per-genome scheduling metadata, computed once before the parallel loop.
typedef struct {
    int64_t seqks_start, seqks_end;  // half-open range in seqks[]
    const char *gname;               // borrowed (seqks[seqks_start].genome)
    const char *spill_path;          // NULL if this genome had no rows
} P2GenomeRange;

// Worker: build all encoded output for `gi`'s genome.  Pure CPU + the spill
// file read.  No shared mutable state touched: spill_path + slen_by_idx +
// the seqks[] window are all read-only.  Result is heap-allocated; the
// writer consumes and frees it.
static Phase2Genome *phase2_genome_work(const P2GenomeRange *gr,
                                        const SeqKey *seqks,
                                        const int64_t *slen_by_idx) {
    int64_t nr = 0;
    Run *runs = NULL;
    char  **seq_dict = NULL;
    int64_t  seq_dict_cap = 0;
    if (gr->spill_path != NULL) {
        runs = load_genome_runs(gr->spill_path, &nr, &seq_dict, &seq_dict_cap);
        if (runs != NULL) qsort(runs, nr, sizeof(Run), run_cmp);
    }

    Phase2Genome *g = st_calloc(1, sizeof(Phase2Genome));
    g->n_seqs   = gr->seqks_end - gr->seqks_start;
    g->seqs     = st_calloc(g->n_seqs ? g->n_seqs : 1, sizeof(Phase2Seq));
    g->nr       = nr;

    int64_t rc = 0;  // forward cursor into sorted runs[]
    for (int64_t i = gr->seqks_start; i < gr->seqks_end; i++) {
        char *sk = seqks[i].seq;
        Phase2Seq *ps = &g->seqs[i - gr->seqks_start];
        ps->seq_name = sk;
        ps->slen      = slen_by_idx[i];
        ps->seqks_idx = i;

        // Find this seq's runs slice in sorted runs[].
        while (rc < nr && strcmp(runs[rc].seq, sk) < 0) rc++;
        int64_t a = rc;
        int64_t bnd = a;
        while (bnd < nr && strcmp(runs[bnd].seq, sk) == 0) bnd++;
        int64_t cnt = bnd - a;
        int64_t *buf = st_malloc((cnt ? 3 * cnt : 1) * sizeof(int64_t));

        // Colinear merge.  Identical logic to the serial path; kept inline
        // rather than factored so this file's only behaviour change is the
        // worker/writer split.
        int64_t m = 0;
        for (int64_t k = a; k < bnd; k++) {
            Run *r = &runs[k];
            int64_t rt = r->t, rg = r->g, rl = r->len;
            char rs = r->strand;
            if (m > 0) {
                int64_t *pp = &buf[3 * (m - 1)];
                int64_t pt = pp[0], pg = pp[1], pl = pp[2] >> 1;
                char ps2 = (pp[2] & 1) ? '-' : '+';
                if (ps2 == rs && rt == pt + pl) {
                    if (rs == '+' && rg == pg + pl) {
                        pp[2] = (pl + rl) << 1;
                        continue;
                    }
                    if (rs == '-' && rg == pg - rl) {
                        pp[1] = rg;
                        pp[2] = ((pl + rl) << 1) | 1;
                        continue;
                    }
                }
            }
            buf[3*m+0] = rt;
            buf[3*m+1] = rg;
            buf[3*m+2] = (rl << 1) | (rs == '-' ? 1 : 0);
            m++;
        }
        qsort(buf, m, 3 * sizeof(int64_t), triple_cmp_by_g);

        // Chunk + encode.  Same chunk-boundary computation as the serial
        // path; we just store the encoded blob + metadata instead of
        // emitting C/R immediately.
        Phase2Chunk *chunks = NULL;
        int64_t n_chunks = 0, chunk_cap = 0;
        int64_t chunk_pos = 0;
        while (chunk_pos < m) {
            int64_t cm_max = m - chunk_pos;
            if (cm_max > TUI_CHUNK_RUNS) cm_max = TUI_CHUNK_RUNS;
            int64_t *cb = &buf[3 * chunk_pos];
            int64_t g_min = cb[1];
            int64_t g_max = cb[1] + (cb[2] >> 1);
            int64_t t_min = cb[0];
            int64_t t_max = cb[0] + (cb[2] >> 1);
            int64_t cm = 1;
            for (; cm < cm_max; cm++) {
                int64_t t  = cb[3*cm + 0];
                int64_t gv = cb[3*cm + 1];
                int64_t le = (cb[3*cm + 2] >> 1);
                int64_t cand_g_max = (gv + le > g_max) ? gv + le : g_max;
                int64_t cand_g_min = (gv < g_min) ? gv : g_min;
                if (cand_g_max - cand_g_min > TUI_CHUNK_G_MAX) break;
                g_max = cand_g_max;
                g_min = cand_g_min;
                if (t < t_min) t_min = t;
                if (t + le > t_max) t_max = t + le;
            }
            if (n_chunks == chunk_cap) {
                chunk_cap = chunk_cap ? chunk_cap * 2 : 4;
                chunks = st_realloc(chunks, chunk_cap * sizeof(Phase2Chunk));
            }
            Phase2Chunk *c = &chunks[n_chunks++];
            c->g_min = g_min; c->g_max = g_max;
            c->t_min = t_min; c->t_max = t_max;
            c->cm = cm;
            c->def = tui_encode_runs(cb, cm, &c->raw_len, &c->def_len);
            // Per-chunk codec self-check (same as serial path).
            int64_t *chk = st_malloc((cm ? 3 * cm : 1) * sizeof(int64_t));
            int64_t dm = decode_runs(c->def, c->def_len, c->raw_len, chk, 3 * cm);
            assert(dm == cm && memcmp(chk, cb, (size_t)(3 * cm) * sizeof(int64_t)) == 0);
            free(chk);
            chunk_pos += cm;
        }
        ps->chunks    = chunks;
        ps->n_chunks  = n_chunks;
        free(buf);
        rc = bnd;
    }

    // runs[].seq are borrowed pointers into seq_dict; free runs first then
    // the dict.  (Order reversed and you free pointers that runs still
    // holds, then read them in the next free() iteration's debug builds.)
    free(runs);
    if (seq_dict != NULL) {
        for (int64_t i = 0; i < seq_dict_cap; i++) free(seq_dict[i]);
        free(seq_dict);
    }
    return g;
}

// Writer: drain one genome's worker result into the OneFile.  Single-
// threaded by construction (the OpenMP `ordered` region in the loop).
// c_ord_emit is the file-position C ordinal counter; we stamp it INTO each
// C record (4th field) so the reader doesn't depend on ONElib's
// oneObject(C) accumulator (off-by-one across cross-type oneGoto in this
// lib version).
static void phase2_genome_write(OneFile *of, Phase2Genome *g, int64_t *c_ord_emit) {
    for (int64_t i = 0; i < g->n_seqs; i++) {
        Phase2Seq *ps = &g->seqs[i];
        oneInt(of, 1) = ps->slen;
        oneWriteLine(of, 'S', strlen(ps->seq_name), (void *)ps->seq_name);
        int64_t s_ord = ps->seqks_idx + 1;
        for (int64_t k = 0; k < ps->n_chunks; k++) {
            Phase2Chunk *c = &ps->chunks[k];
            ++(*c_ord_emit);
            oneInt(of, 0) = c->g_min;
            oneInt(of, 1) = c->g_max;
            oneInt(of, 2) = s_ord;
            oneInt(of, 3) = *c_ord_emit;
            oneInt(of, 4) = c->t_min;
            oneInt(of, 5) = c->t_max;
            oneWriteLine(of, 'C', 0, NULL);
            oneInt(of, 0) = c->raw_len;
            oneWriteLine(of, 'R', c->def_len, c->def);
            free(c->def);
        }
        free(ps->chunks);
    }
    free(g->seqs);
    free(g);
}

/* Open a new .tui in WRITE mode, install the crash-cleanup atexit hook,
 * stamp the provenance line.  Returns the OneFile* on success (caller
 * owns); writes *schema_out (caller-owned via oneSchemaDestroy after
 * oneFileClose).  Returns NULL on open failure -- caller is responsible
 * for any pre-allocated state cleanup and calling tui_atexit_disarm().
 *
 * Exported via tui.h (internal-but-shared with taf_coarsen) so any tool
 * writing a .tui-format file uses the same boilerplate and the same
 * crash safety net.  prog/what/blurb populate the OneCode provenance
 * record (oneAddProvenance) -- same call shape tui_create uses. */
OneFile *tui_open_write(const char *out_path, const char *prog,
                        const char *what, const char *blurb,
                        OneSchema **schema_out) {
    OneSchema *schema = oneSchemaCreateFromText(TUI_SCHEMA);
    OneFile *of = oneFileOpenWriteNew(out_path, schema, "tui", true, 1);
    if (of == NULL) {
        fprintf(stderr, "tui: cannot write %s\n", out_path);
        oneSchemaDestroy(schema);
        if (schema_out) *schema_out = NULL;
        return NULL;
    }
    tui_atexit_track_tui(out_path);
    oneAddProvenance(of, prog, what, blurb, 0);
    if (schema_out) *schema_out = schema;
    return of;
}

int tui_create(LI *li, const char *out_path, const char *tmp_dir,
                stHash *genome_name_map, int n_threads) {
    int input_format = check_input_format(LI_peek_at_next_line(li));  // 0=TAF 1=MAF
    if (input_format != 0 && input_format != 1) {
        fprintf(stderr, "tui: input must be MAF or TAF\n");
        return 1;
    }
    bool rle = 0;
    Tag *tag = (input_format == 0) ? taf_read_header_2(li, &rle) : maf_read_header(li);
    // Genome resolution map precedence: an explicit -n list overrides; else
    // fall back to the genome set in the "# hal" tree comment if present.
    stHash *eff_map = genome_name_map, *tree_map = NULL;
    if (eff_map == NULL) {
        Tag *h = tag_find(tag, (char *) TAF_HAL_TREE_KEY);
        if (h != NULL) {
            tree_map = genome_set_from_newick(h->value);
            eff_map = tree_map;
            st_logInfo("tui: using genome set from the # hal tree (%" PRIi64
                       " genomes)\n", (int64_t) stHash_size(tree_map));
        }
    }
    tag_destruct(tag);

    // spill dir: explicit --tmpDir, else the output file's own directory
    char *eff_tmp = (tmp_dir != NULL && tmp_dir[0] != '\0')
                        ? stString_copy(tmp_dir) : dir_of(out_path);

    // Raise RLIMIT_NOFILE soft to hard.  Phase 1 cycles one FD per open spill
    // (one per genome in the cap) plus a handful for stdin/out/err, the MAF
    // stream, and the idx spill.  At vertebrate scale (1k+ genomes) the
    // default Linux soft limit (1024) leaves no headroom; the user's 1153-
    // genome 577-way run failed at spill #1019 with this exact pattern.
    // Raising soft to hard is a no-op on systems already configured for
    // hi-FD work and free headroom elsewhere.
    struct rlimit rl = { .rlim_cur = 1024, .rlim_max = 1024 };  // safe defaults
    int max_open_cap = 1024;            // safe-fallback cap if getrlimit fails
    if (getrlimit(RLIMIT_NOFILE, &rl) == 0) {
        rlim_t want = rl.rlim_max;
        if (rl.rlim_cur < want) {
            rl.rlim_cur = want;
            (void) setrlimit(RLIMIT_NOFILE, &rl);   // best-effort; ignore EPERM
            (void) getrlimit(RLIMIT_NOFILE, &rl);   // re-read what we actually got
        }
        // Reserve 32 FDs for stdin/out/err + idx + maf + libc misc; cap at a
        // sane ceiling so we don't burn the whole process FD table on spills.
        rlim_t avail = (rl.rlim_cur > 32) ? rl.rlim_cur - 32 : 0;
        if (avail > 8192) avail = 8192;
        if (avail < 8)    avail = 8;          // LRU below cap=8 thrashes badly
        max_open_cap = (int) avail;
    }
    // Env-var override -- useful for tests (force the LRU into play on small
    // fixtures) and for diagnosing thrash vs FD-pressure issues in the field.
    const char *env_cap = getenv("TAFFY_TUI_MAX_OPEN");
    if (env_cap != NULL && env_cap[0] != '\0') {
        char *end = NULL;
        long v = strtol(env_cap, &end, 10);
        if (end != env_cap && v >= 4 && v <= 65536) max_open_cap = (int) v;
    }

    Phase1 p1;
    memset(&p1, 0, sizeof(p1));
    p1.tmp_dir = eff_tmp;
    p1.out_base = base_of(out_path);
    p1.next_spill_id = 0;
    p1.gmap = eff_map;
    // spill_ents value freed inline by tui_cleanup (the SpillEnt struct +
    // its path string are not standalone-freeable via a single callback).
    p1.spill_ents = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    p1.max_open = max_open_cap;
    p1.n_open   = 0;
    p1.lru_head = p1.lru_tail = NULL;
    // per-spill seq-name dictionaries: values are inner stHashes (destroyed
    // here on cleanup so their per-record int64_t* + seq_name copies go too).
    p1.spill_seqs = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
                                      free, (void (*)(void *))stHash_destruct);
    p1.seq_len = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    p1.seq_keys = stList_construct();
    p1.T = 0;
    p1.idx_fh = NULL; p1.idx_path = NULL; p1.idx_last_col = 0; p1.idx_n = 0;
    st_logInfo("tui: phase-1 spill FH pool: max_open=%d (RLIMIT_NOFILE soft=%lu)\n",
               p1.max_open, (unsigned long) rl.rlim_cur);

    // idx_anchor captures the block's anchor offset itself (LI_get_position =
    // start of the peeked block-first line; the exact tai_create_taf
    // convention so LI_seek + the shared readers behave identically at query).
    time_t t_total_start = time(NULL);
    time_t t_phase1_start = t_total_start;
    st_logInfo("tui: starting phase 1 (streaming %s scan, spilling per-genome runs)\n",
               input_format == 0 ? "TAF" : "MAF");
    int64_t block_count = 0;
    Alignment *aln, *p_aln = NULL;
    // Periodic in-phase log -- the only signal a user has during a multi-hour
    // run.  Time-gated rather than block-gated so the cadence doesn't depend
    // on block size (vertebrate-scale blocks vary by ~3 orders of magnitude).
    // The interval is short for quick smoke tests + still useful for the
    // 27-hour 577-way build (~30 messages/hour).
    //
    // Per-iteration cost is just `++tick_ctr & MASK == 0` (one int increment
    // + AND + branch, all branch-predicted in the common case).  The actual
    // time(NULL) check happens only once per TUI_TICK_INTERVAL blocks.  This
    // matters at vertebrate scale: a 1B-block input previously made 1B vDSO
    // time() calls per block-emit (~5-10 s total) -- with the throttle that
    // drops to ~1M calls (~10 ms), well under the noise floor of the I/O
    // dominated phase-1 loop.  At 10k blocks/s the time() check fires
    // ~10/s, still way more frequently than the 600 s log threshold.
    time_t t_last_log = t_phase1_start;
    const int64_t LOG_EVERY_SEC = 600;
    const int64_t TUI_TICK_INTERVAL = 1024;       // power-of-2 for bitmask
    const int64_t TUI_TICK_MASK     = TUI_TICK_INTERVAL - 1;
    int64_t tick_ctr = 0;
    #define TUI_PHASE1_TICK() do {                                              \
        if ((++tick_ctr & TUI_TICK_MASK) != 0) break;                           \
        time_t now = time(NULL);                                                \
        if (now - t_last_log >= LOG_EVERY_SEC) {                                \
            int64_t elapsed = (int64_t)(now - t_phase1_start);                  \
            st_logInfo("tui: phase 1 progress: %" PRIi64 " blocks, "            \
                       "%" PRIi64 " universal columns, "                        \
                       "%" PRIi64 " genomes seen, %" PRIi64 "s elapsed "        \
                       "(%.1f kblocks/s)\n",                                    \
                       block_count, p1.T,                                       \
                       (int64_t) stHash_size(p1.spill_ents),                    \
                       elapsed,                                                 \
                       elapsed > 0 ? (block_count / 1000.0) / elapsed : 0.0);   \
            t_last_log = now;                                                   \
        }                                                                       \
    } while (0)
    if (input_format == 1) {
        while (1) {
            idx_anchor(&p1, li, 1, false);
            aln = maf_read_block(li);
            if (aln == NULL) break;
            tui_phase1_block(&p1, aln);
            alignment_destruct(aln, 1);
            block_count++;
            TUI_PHASE1_TICK();
        }
    } else {
        while (1) {
            idx_anchor(&p1, li, 0, rle);
            aln = taf_read_block(p_aln, rle, li);
            if (aln == NULL) break;
            tui_phase1_block(&p1, aln);
            if (p_aln != NULL) alignment_destruct(p_aln, 1);
            p_aln = aln;
            block_count++;
            TUI_PHASE1_TICK();
        }
        if (p_aln != NULL) alignment_destruct(p_aln, 1);
    }
    #undef TUI_PHASE1_TICK
    // fclose surfaces buffered short writes (e.g. disk full) -> fail loudly
    // here rather than emit a silently-truncated spill that Phase 2 would
    // later misread as "corrupt spill line".
    if (p1.idx_fh != NULL) {
        if (fclose(p1.idx_fh) != 0) st_errAbort("tui: idx spill close failed "
            "(disk full / write error) -- %s", p1.idx_path);
        p1.idx_fh = NULL;
    }
    // close any currently-open spill FILE*s (paths kept on each SpillEnt
    // for phase 2).  Under the LRU pool only n_open ents are open at this
    // point; previously-evicted ones already had their write buffers
    // flushed at eviction-fclose time, so we only need to close survivors.
    stHashIterator *hit = stHash_getIterator(p1.spill_ents);
    char *gk;
    while ((gk = stHash_getNext(hit)) != NULL) {
        SpillEnt *e = stHash_search(p1.spill_ents, gk);
        if (e->fh != NULL) {
            if (fclose(e->fh) != 0) {
                st_errAbort("tui: spill close failed for genome '%s' "
                            "(disk full / write error) -- %s", gk, strerror(errno));
            }
            e->fh = NULL;
        }
    }
    stHash_destructIterator(hit);
    int64_t n_genomes = (int64_t) stHash_size(p1.spill_ents);
    st_logInfo("tui: phase 1 done in %" PRIi64 " seconds "
               "(T=%" PRIi64 " columns, %" PRIi64 " blocks, "
               "%" PRIi64 " genomes, %" PRIi64 " idx anchors)\n",
               (int64_t)(time(NULL) - t_phase1_start),
               p1.T, block_count, n_genomes, p1.idx_n);

    // deterministic global order: genome-major (true resolved genome), then
    // sequence.  Resolve each genome once here (memoized in SeqKey.genome).
    int64_t n_seqs = stList_length(p1.seq_keys);
    SeqKey *seqks = st_malloc((n_seqs ? n_seqs : 1) * sizeof(SeqKey));
    for (int64_t i = 0; i < n_seqs; i++) {
        seqks[i].seq = stList_get(p1.seq_keys, i);          // owned by seq_len
        seqks[i].genome = genome_of(seqks[i].seq, eff_map); // we own this
    }
    qsort(seqks, n_seqs, sizeof(SeqKey), seqkey_cmp);

    OneSchema *schema = NULL;
    OneFile *of = tui_open_write(out_path, "taffy", "tui",
                                 "universal column index", &schema);
    if (of == NULL) {
        tui_cleanup(&p1, seqks, n_seqs, eff_tmp, tree_map);
        tui_atexit_disarm();
        return 1;
    }

    // t: total columns
    oneInt(of, 0) = p1.T;
    oneWriteLine(of, 't', 0, NULL);

    // X: universal-column -> file_pos index (the .tai-equivalent).  Read the
    // column-ordered anchor spill, encode (deflated SoA), write one X line.
    time_t t_idxX_start = time(NULL);
    st_logInfo("tui: Index X: encoding universal-column->file-offset anchors\n");
    {
        int64_t *ic = NULL, *iff = NULL, in = 0, icap = 0;
        if (p1.idx_path != NULL) {
            FILE *xf = fopen(p1.idx_path, "r");
            if (xf == NULL) st_errAbort("tui: cannot reopen idx spill %s", p1.idx_path);
            char line[4096]; int64_t c, f, prevc = -1;
            while (fgets(line, sizeof(line), xf) != NULL) {
                if (sscanf(line, "%" SCNi64 "\t%" SCNi64, &c, &f) != 2) {
                    st_errAbort("tui: corrupt idx spill line: %s", line);
                }
                if (c <= prevc) st_errAbort("tui: idx anchors not strictly "
                    "increasing (%" PRIi64 " after %" PRIi64 ")", c, prevc);
                prevc = c;
                if (in == icap) { icap = icap ? icap*2 : 4096;
                    ic = st_realloc(ic, icap*sizeof(int64_t));
                    iff = st_realloc(iff, icap*sizeof(int64_t)); }
                ic[in] = c; iff[in] = f; in++;
            }
            fclose(xf);
        }
        int64_t x_raw = 0, x_def = 0;
        uint8_t *xdef = tui_encode_idx(ic, iff, in, &x_raw, &x_def);
        int64_t *cc = st_malloc((in?in:1)*sizeof(int64_t));
        int64_t *cf = st_malloc((in?in:1)*sizeof(int64_t));
        int64_t dn = decode_idx(xdef, x_def, x_raw, cc, cf);
        assert(dn == in);
        for (int64_t k = 0; k < in; k++) assert(cc[k]==ic[k] && cf[k]==iff[k]);
        free(cc); free(cf);
        oneInt(of, 0) = x_raw;
        oneInt(of, 1) = in;
        oneWriteLine(of, 'X', x_def, xdef);
        st_logInfo("tui: Index X done in %" PRIi64 " seconds "
                   "(%" PRIi64 " anchors, %" PRIi64 " deflated bytes)\n",
                   (int64_t)(time(NULL) - t_idxX_start), in, x_def);
        free(xdef); free(ic); free(iff);
    }

    // d: directory.  Written in NAME-sorted order (independent of the seqks
    // genome-major order used for S/R), so the reader can binary-search the
    // d-lines via oneGoto-by-index without preloading the directory.  The
    // S-ordinal field still points back to the matching S object (== position
    // in seqks).
    st_logInfo("tui: writing directory (%" PRIi64 " sequences across %" PRIi64
               " genomes)\n", n_seqs, n_genomes);
    {
        DKey *dks = st_malloc((n_seqs ? n_seqs : 1) * sizeof(DKey));
        for (int64_t i = 0; i < n_seqs; i++) {
            dks[i].idx = i; dks[i].name = seqks[i].seq;
        }
        qsort(dks, n_seqs, sizeof(DKey), dkey_cmp);
        for (int64_t k = 0; k < n_seqs; k++) {
            int64_t i = dks[k].idx;
            char *sk = seqks[i].seq;
            int64_t slen = *(int64_t *)stHash_search(p1.seq_len, sk);
            oneInt(of, 1) = i;          // S-ord = position in seqks
            oneInt(of, 2) = slen;
            oneWriteLine(of, 'd', strlen(sk), (void *)sk);
        }
        free(dks);
    }

    // per genome: g + S/R objects.  seqks is sorted by the TRUE resolved
    // genome, so one genome's sequences are ONE contiguous block by
    // construction (no first-dot-collision hazard, no re-resolve).
    //
    // Phase 2 is per-genome-independent CPU work; the actual writes into
    // OneFile must stay serial (the c_ord_emit counter is the file-position
    // C ordinal stamped INTO each C record, and ONElib isn't thread-safe).
    // We split into worker (sort/merge/chunk/encode -> Phase2Genome buffer)
    // and writer (oneWriteLine S/C/R) so the parallel loop below can use
    // schedule(dynamic) + ordered: workers in parallel, writes in seqks
    // order.  See the Phase2Genome / phase2_genome_work / phase2_genome_write
    // block above tui_create for the data shapes.
    time_t t_phase2_start = time(NULL);
    st_logInfo("tui: starting phase 2 (per-genome sort + encode + write, "
               "%" PRIi64 " genomes)\n", n_genomes);

    // Per-genome scheduling metadata: seqks window + cached spill_path.
    // Building this upfront (single thread) avoids stHash_search in workers
    // (sonLib's hash isn't reentrant for concurrent reads).
    P2GenomeRange *gr = st_malloc((n_genomes ? n_genomes : 1) * sizeof(P2GenomeRange));
    {
        int64_t gi = 0, j = 0;
        while (j < n_seqs) {
            char *gname = seqks[j].genome;
            gr[gi].seqks_start = j;
            gr[gi].gname       = gname;
            SpillEnt *ent      = stHash_search(p1.spill_ents, gname);
            gr[gi].spill_path  = (ent != NULL) ? ent->path : NULL;
            while (j < n_seqs && strcmp(seqks[j].genome, gname) == 0) j++;
            gr[gi].seqks_end = j;
            gi++;
        }
        assert(gi == n_genomes);
    }
    // Flat per-seqks slen lookup -- workers read this directly instead of
    // p1.seq_len's stHash.
    int64_t *slen_by_idx = st_malloc((n_seqs ? n_seqs : 1) * sizeof(int64_t));
    for (int64_t k = 0; k < n_seqs; k++) {
        slen_by_idx[k] = *(int64_t *)stHash_search(p1.seq_len, seqks[k].seq);
    }

    // File-order C ordinal counter; assigned inside the writer.
    int64_t c_ord_emit = 0;

    // Parallel-for over genomes.  schedule(dynamic, 1) gives one genome
    // per task -- well-matched to the heavy-tailed work (median 7s,
    // p95 407s on 577-way; the rare giants would clog a static partition).
    // ordered serialises the write region in iteration order = seqks order
    // so the c_ord_emit counter and OneFile state stay un-contended and
    // S record positions match the s_ord values already committed to the
    // d directory above.  n_threads<=1 falls back to a serial loop (no
    // OpenMP runtime cost).
    //
    // MEMORY: each in-flight worker holds the FULL runs[] for its genome
    // (loaded by load_genome_runs, ~Run-struct + seq-string-copy per run)
    // PLUS its Phase2Genome with all encoded chunks.  For the 577-way
    // root-ish ancestors that's tens of GB per worker; with n_threads=64
    // we observed an OOM near the start of phase 2 on a 1.5 TB box.
    // Caller (taf_index.c) defaults phase2-threads to 1 for safety;
    // user opts into more with --phase2Threads after sizing their
    // giants.  ordered does NOT bound in-flight memory: with dynamic
    // schedule the worker fetches the next iter immediately on yield,
    // so n_threads workers ARE running simultaneously even though the
    // ordered region serialises their writes.
    int nt = (n_threads > 1) ? n_threads : 1;
    #pragma omp parallel for schedule(dynamic, 1) num_threads(nt) ordered
    for (int64_t gi = 0; gi < n_genomes; gi++) {
        time_t t_g_start = time(NULL);
        Phase2Genome *gres = phase2_genome_work(&gr[gi], seqks, slen_by_idx);
        time_t t_g_done = time(NULL);
        #pragma omp ordered
        {
            int64_t nr_saved = gres->nr;
            phase2_genome_write(of, gres, &c_ord_emit);
            st_logInfo("tui: phase 2 genome %" PRIi64 "/%" PRIi64 " '%s' "
                       "worked in %" PRIi64 "s, written in %" PRIi64 "s "
                       "(%" PRIi64 " runs)\n",
                       gi + 1, n_genomes, gr[gi].gname,
                       (int64_t)(t_g_done - t_g_start),
                       (int64_t)(time(NULL) - t_g_done),
                       nr_saved);
        }
    }
    // --- Genome roster (g records) -----------------------------------------
    // One record per resolved genome with (name, total_bp, n_chroms).
    // Written here -- after every per-seq d/S/C/R is on disk -- so the
    // reader has a stable, deterministic genome list without having to
    // guess where the "<genome>.<sequence>" split is in d-line keys
    // (matters for NCBI-style versioned accessions like
    // "GCA_028858775.2" where the genome name itself contains a dot).
    for (int64_t gi = 0; gi < n_genomes; gi++) {
        int64_t total_bp = 0;
        int64_t n_chroms = gr[gi].seqks_end - gr[gi].seqks_start;
        for (int64_t i = gr[gi].seqks_start; i < gr[gi].seqks_end; i++)
            total_bp += slen_by_idx[i];
        oneInt(of, 1) = total_bp;
        oneInt(of, 2) = n_chroms;
        oneWriteLine(of, 'g', strlen(gr[gi].gname), (void *)gr[gi].gname);
    }

    free(gr);
    free(slen_by_idx);
    st_logInfo("tui: phase 2 done in %" PRIi64 " seconds\n",
               (int64_t)(time(NULL) - t_phase2_start));

    oneFileClose(of);
    oneSchemaDestroy(schema);
    st_logInfo("tui: build complete in %" PRIi64 " seconds total\n",
               (int64_t)(time(NULL) - t_total_start));
    // We made it; release the crash-cleanup safety net (the atexit list owns
    // its own strdup'd copies, so order vs tui_cleanup is not load-bearing,
    // but disarming first avoids any chance of remove()ing a finalized file).
    tui_atexit_disarm();
    tui_cleanup(&p1, seqks, n_seqs, eff_tmp, tree_map);
    return 0;
}

/////////////////////////////////////////////////////////////////////////////
// Reader / query side.  Open keeps only the small global state in RAM (T +
// the X-anchor table); the directory and per-sequence runs stay on disk and
// are resolved per query: binary-search the name-sorted d-lines via
// oneGoto-by-index for the S-ordinal + seq length, then oneGoto on S to fetch
// the matching R blob.  No directory hashes are built at open time.
/////////////////////////////////////////////////////////////////////////////

struct _Tui {
    char   *path;       // .tui path
    int64_t T;          // total universal columns
    int64_t n_d;        // number of d-lines (binary-search upper bound)
    // Universal-column -> file_pos index (X track); both strictly increasing.
    int64_t *idxCol, *idxFpos;
    int64_t  idxN;
    // Cached OneFile handle, reused across tui_query / tui_load_seq_runs.
    // oneFileOpenRead on a multi-GB .tui parses the embedded schema + reads
    // the full footer index for every object type (S / C / R / d / X / t) --
    // ~1-2 s on a 92 GB cluster .tui.  Previously we re-opened on every call,
    // so an N-bed-line lift paid N+3 opens.  Cached, the cost is paid once at
    // tui_load.  Callers must serialize on the same Tui * (see tui.h) -- the
    // cursor and read buffer are mutated by every oneGoto / oneReadLine.
    OneFile *of;
};

Tui *tui_load(const char *tui_path) {
    OneFile *of = oneFileOpenRead(tui_path, NULL, "tui", 1);
    if (of == NULL) return NULL;
    Tui *tui = st_calloc(1, sizeof(Tui));
    tui->path = stString_copy(tui_path);
    // Read only `t` and `X` (small, used by every query).  Stop before the
    // d-lines.  The directory is left on disk -- tui_query binary-searches
    // it per call (O(log n_d) seeks).
    int seen_t = 0, seen_x = 0;
    char c;
    while (!(seen_t && seen_x) && (c = oneReadLine(of)) != 0) {
        if (c == 't') {
            tui->T = oneInt(of, 0);
            seen_t = 1;
        } else if (c == 'X') {
            int64_t x_raw = oneInt(of, 0);
            int64_t xn    = oneInt(of, 1);
            int64_t x_def = oneLen(of);
            const uint8_t *xdef = (const uint8_t *)oneString(of);
            int64_t cap = xn ? xn : 1;
            tui->idxCol  = st_malloc(cap * sizeof(int64_t));
            tui->idxFpos = st_malloc(cap * sizeof(int64_t));
            tui->idxN = decode_idx(xdef, x_def, x_raw, tui->idxCol, tui->idxFpos);
            assert(tui->idxN == xn);
            seen_x = 1;
        }
    }
    // Directory size for the binary-search bound; populated from the file
    // footer so it works regardless of how far we read above.
    I64 nd = 0;
    oneStats(of, 'd', &nd, NULL, NULL);
    tui->n_d = nd;
    // Keep the handle around for the lift hot path -- closing here and
    // re-opening per query was costing ~1-2 s per call on multi-GB .tui's.
    tui->of = of;
    return tui;
}

void tui_destruct(Tui *tui) {
    if (tui == NULL) return;
    if (tui->of != NULL) oneFileClose(tui->of);
    free(tui->idxCol); free(tui->idxFpos);
    free(tui->path);
    free(tui);
}

int64_t tui_total_columns(const Tui *tui) { return tui->T; }

int64_t        tui_idx_n   (const Tui *tui) { return tui->idxN; }
const int64_t *tui_idx_cols(const Tui *tui) { return tui->idxCol; }
const int64_t *tui_idx_fpos(const Tui *tui) { return tui->idxFpos; }

// Binary-search the d-lines (name-sorted by the writer) for `seq_name`.
// Returns the S-ordinal (0-indexed; position in the genome-major S/R order),
// or -1 if not found.  Optionally fills *seqlen.
// `of` must be a freshly opened (or otherwise idle) tui OneFile -- we move
// its position via oneGoto.
static int64_t tui_find_d(OneFile *of, int64_t n_d,
                          const char *seq_name, int64_t *seqlen_out) {
    int64_t lo = 1, hi = n_d;
    // 8192 mirrors the writer's spill scanner "%8191s" upper bound -- the two
    // sides need to move together if either gets bumped.
    char buf[8192];
    while (lo <= hi) {
        int64_t mid = lo + (hi - lo) / 2;        // overflow-safe shape
        if (!oneGoto(of, 'd', mid)) return -1;
        if (oneReadLine(of) != 'd') return -1;
        int64_t n = oneLen(of);
        if (n < 0 || n >= (int64_t)sizeof(buf)) return -1;  // see size note above
        memcpy(buf, oneString(of), (size_t)n);
        buf[n] = '\0';
        int cmp = strcmp(buf, seq_name);
        if (cmp == 0) {
            if (seqlen_out != NULL) *seqlen_out = oneInt(of, 2);
            return oneInt(of, 1);
        } else if (cmp < 0) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }
    return -1;
}

static int tui_iv_cmp(const void *a, const void *b) {
    const TuiInterval *x = a, *y = b;
    return (x->start < y->start) ? -1 : (x->start > y->start) ? 1 : 0;
}

TuiInterval *tui_query(Tui *tui, const char *seq_name,
                       int64_t start, int64_t end, int64_t *n_out) {
    *n_out = 0;
    if (start >= end) return NULL;

    OneFile *of = tui->of;
    if (of == NULL) return NULL;

    // Resolve name -> S-ordinal by binary-searching the (name-sorted) d-lines
    // via oneGoto; O(log n_d) seeks, no preloaded directory hashes.
    int64_t ord = tui_find_d(of, tui->n_d, seq_name, NULL);
    if (ord < 0) return NULL;

    // Jump to the (ord+1)-th S object via the ONElib footer index; after the
    // goto, oneReadLine() returns S then (C, R)+ chunk pairs.
    if (!oneGoto(of, 'S', ord + 1)) return NULL;
    int64_t *runs = NULL;
    char c = oneReadLine(of);
    if (c != 'S') return NULL;
    int64_t sn = oneLen(of);
    if (sn != (int64_t)strlen(seq_name) ||
        memcmp(oneString(of), seq_name, sn) != 0) {
        return NULL;
    }
    // Per-chunk t-range skip.  C records added t_min / t_max in field 4/5 to
    // let the reader skip chunks whose source-coord range can't overlap
    // [start, end) without paying the zlib decompress of the chunk's R blob
    // (decode_runs+inflate is ~77% of taffy lift wall on rodents-scale .tui).
    // Backward-compat: old .tui without the extra fields has nField == 4 on
    // C; we treat that as "no skip info, decode every chunk" -- same wall
    // as before this change.
    int c_has_t_range = (of->info[(int)'C'] != NULL &&
                         of->info[(int)'C']->nField >= 6);
    int64_t total = 0, cap = 0;
    while ((c = oneReadLine(of)) == 'C') {
        int skip = 0;
        if (c_has_t_range) {
            int64_t t_min = oneInt(of, 4);
            int64_t t_max = oneInt(of, 5);
            if (t_max <= start || t_min >= end) skip = 1;
        }
        if (oneReadLine(of) != 'R') break;
        if (skip) continue;
        int64_t raw_len = oneInt(of, 0);
        int64_t def_len = oneLen(of);
        const uint8_t *def = (const uint8_t *)oneString(of);
        int64_t need = total + raw_len + 3;
        if (need > cap) {
            cap = need;
            runs = st_realloc(runs, (cap ? cap : 1) * sizeof(int64_t));
        }
        int64_t cm = decode_runs(def, def_len, raw_len, runs + total, raw_len + 3);
        total += 3 * cm;
    }
    int64_t m = total / 3;
    if (m == 0) { free(runs); return NULL; }

    // Clip each overlapping run to [start,end), map to a column interval,
    // and remember the source-genome position at iv.start + the rev bit.
    // For rev=0 the column-ascending direction matches source-ascending,
    // so iv.start corresponds to the LOW source coord `a`.  For rev=1
    // the directions are opposite, so iv.start corresponds to the HIGH
    // source coord (b - 1; we store b - 1 explicitly so the visitor can
    // recover source_pos = iv.t_start - (c - iv.start)).
    TuiInterval *iv = st_malloc((size_t)m * sizeof(TuiInterval));
    int64_t k = 0;
    for (int64_t r = 0; r < m; r++) {
        int64_t t = runs[3*r+0], g = runs[3*r+1], lenc = runs[3*r+2];
        int64_t len = lenc >> 1, rev = lenc & 1, te = t + len;
        int64_t a = start > t ? start : t;
        int64_t b = end < te ? end : te;
        if (a >= b) continue;
        if (!rev) {
            iv[k].start   = g + (a - t);
            iv[k].end     = g + (b - t);
            iv[k].t_start = a;
            iv[k].rev     = 0;
        } else {
            iv[k].start   = g + (t + len - b);
            iv[k].end     = g + (t + len - a);
            iv[k].t_start = b - 1;
            iv[k].rev     = 1;
        }
        k++;
    }
    free(runs);
    if (k == 0) { free(iv); return NULL; }

    // Sort by (start, rev) so same-rev adjacent intervals end up next to
    // each other; the existing tui_iv_cmp orders by start which is fine
    // -- we add the rev / t_start safety check in the merge step.
    qsort(iv, k, sizeof(TuiInterval), tui_iv_cmp);

    // Merge only if (a) the two intervals overlap or abut on the column
    // axis AND (b) they share the same rev AND (c) the joined t_start
    // mapping is linearly consistent.  Without (b)+(c), merging would
    // erase rev / t_start information the visitor needs to compute the
    // correct source pos + relative strand.
    int64_t w = 0;
    for (int64_t i = 1; i < k; i++) {
        int merge_ok = (iv[i].start <= iv[w].end) && (iv[i].rev == iv[w].rev);
        if (merge_ok) {
            // Check the t_start mapping is contiguous through the abut.
            int64_t cw_extent = iv[w].end - iv[w].start;
            int64_t expected_t_at_i_start = (iv[w].rev == 0)
                ? iv[w].t_start + cw_extent
                : iv[w].t_start - cw_extent;
            if (iv[i].t_start != expected_t_at_i_start) merge_ok = 0;
        }
        if (merge_ok) {
            if (iv[i].end > iv[w].end) iv[w].end = iv[i].end;
        } else {
            iv[++w] = iv[i];
        }
    }
    *n_out = w + 1;
    return iv;
}

int64_t *tui_load_seq_runs(Tui *tui, const char *seq_name, int64_t *n_out) {
    *n_out = 0;
    if (tui == NULL || seq_name == NULL) return NULL;

    OneFile *of = tui->of;
    if (of == NULL) return NULL;

    int64_t ord = tui_find_d(of, tui->n_d, seq_name, NULL);
    if (ord < 0) return NULL;
    if (!oneGoto(of, 'S', ord + 1)) return NULL;
    char c = oneReadLine(of);
    if (c != 'S') return NULL;
    int64_t sn = oneLen(of);
    if (sn != (int64_t)strlen(seq_name) ||
        memcmp(oneString(of), seq_name, sn) != 0) {
        return NULL;
    }
    // Walk (C, R)+ chunk pairs and concat decoded triples into one flat
    // array.  NB: tui_query has a per-chunk t-range skip (see C-record
    // fields 4/5 + the c_has_t_range gate), this function deliberately
    // does NOT skip -- its contract is "return ALL runs of the seq" so
    // batch callers can binary-search the result by t_start in RAM.
    // Skipping by t-range would defeat that contract.
    int64_t *runs = NULL;
    int64_t total = 0, cap = 0;
    while ((c = oneReadLine(of)) == 'C') {
        if (oneReadLine(of) != 'R') break;
        int64_t raw_len = oneInt(of, 0);
        int64_t def_len = oneLen(of);
        const uint8_t *def = (const uint8_t *)oneString(of);
        int64_t need = total + raw_len + 3;
        if (need > cap) {
            cap = need;
            runs = st_realloc(runs, (cap ? cap : 1) * sizeof(int64_t));
        }
        int64_t cm = decode_runs(def, def_len, raw_len, runs + total, raw_len + 3);
        total += 3 * cm;
    }
    if (total == 0) { free(runs); return NULL; }
    int64_t n = total / 3;
    // Writer pre-sorts each chunk's runs by g_start for the column-keyed
    // target index.  Source-side callers (taffy lift's per-record coord
    // lookup) binary-search by t_start instead, so restore the t-sorted
    // view here.  O(n log n) once per source seq -- much smaller than the
    // O(N) loop we feed the result into.
    if (n > 1) qsort(runs, n, 3 * sizeof(int64_t), triple_cmp_by_t);
    *n_out = n;
    return runs;
}

/////////////////////////////////////////////////////////////////////////////
// Reverse lookup: universal column -> target genome coord.
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    int64_t g_start;     // universal-column start of the run
    int64_t t_start;     // forward coord in the genome sequence
    int64_t length;      // run length in bases
    int64_t seq_idx;     // index into TuiGenomeLift::seq_names
    int     strand;      // 1 = '+', 0 = '-'
} GLRun;

// One chunk of a sequence's runs.  The .tui's C metadata is read at open
// (cheap) into these structs; the GLRun array is built lazily on first
// query hit via chunk_decode().  g_min / g_max enable a "skip if column
// outside" test before paying the decode cost.
typedef struct {
    int64_t g_min, g_max;        // column range from the C record
    int64_t c_ord;               // 1-based C ordinal for oneGoto lazy load
    int64_t seq_idx;             // index into TuiGenomeLift::seq_names
    GLRun  *runs;                // NULL until decoded; sorted by g_start
    int64_t n_runs;
    // max_end_prefix[i] = max(runs[k].g_start + runs[k].length for k in
    // [0..i]).  Used as an early-exit during the backward scan in
    // chunk_collect: if max_end_prefix[j-1] <= column, no earlier run can
    // cover.  Bounds the scan to O(log n + paralog_count).
    int64_t *max_end_prefix;
} TGLChunk;

struct _TuiGenomeLift {
    char     **seq_names;        // owned; one per sequence loaded
    int64_t    n_seq;
    TGLChunk  *chunks;           // owned; sorted by g_min
    int64_t    n_chunks;
    // Running max(chunk.g_max) over chunks[0..i] -- same shape as
    // TGLChunk::max_end_prefix, but over chunks for the outer scan.
    int64_t   *chunk_max_end;
    // OneFile cursor for lazy R decoding.  Single-threaded; callers
    // serialize (blockViz holds g_mutex; CLI tools are single-threaded).
    OneFile   *of;
};

static int glrun_cmp(const void *a, const void *b) {
    const GLRun *x = a, *y = b;
    return (x->g_start < y->g_start) ? -1 : (x->g_start > y->g_start) ? 1 : 0;
}

static int tglchunk_cmp_gmin(const void *a, const void *b) {
    const TGLChunk *x = a, *y = b;
    return (x->g_min < y->g_min) ? -1 : (x->g_min > y->g_min) ? 1 : 0;
}

// Smallest d-line ordinal whose name >= prefix (1-based, like tui_find_d),
// or n_d+1 if all names sort before prefix.  Uses oneGoto-by-index just like
// tui_find_d; same 8192 buf size convention.
static int64_t tui_find_d_lower_bound(OneFile *of, int64_t n_d, const char *prefix) {
    size_t plen = strlen(prefix);
    int64_t lo = 1, hi = n_d + 1;
    char buf[8192];
    while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        if (!oneGoto(of, 'd', mid)) return n_d + 1;
        if (oneReadLine(of) != 'd')  return n_d + 1;
        int64_t n = oneLen(of);
        if (n < 0 || n >= (int64_t)sizeof(buf)) return n_d + 1;
        memcpy(buf, oneString(of), (size_t)n);
        buf[n] = '\0';
        int cmp = strncmp(buf, prefix, plen);
        // First plen chars compared.  If equal so far, buf >= prefix iff
        // either buf is strictly longer or strcmp == 0 (== prefix means
        // buf == prefix exactly, also >= prefix).
        int ge = (cmp > 0) || (cmp == 0);
        if (ge) hi = mid; else lo = mid + 1;
    }
    return lo;
}

stHash *tui_sequence_lengths(const char *tui_path) {
    if (tui_path == NULL) return NULL;
    OneFile *of = oneFileOpenRead(tui_path, NULL, "tui", 1);
    if (of == NULL) return NULL;

    // n_d is the d-line count; populated from the footer regardless of
    // where the cursor is.
    I64 n_d = 0;
    oneStats(of, 'd', &n_d, NULL, NULL);

    // Value field is the seq length, stored as a void* (cast via intptr_t).
    // Same convention tai_sequence_lengths uses.
    stHash *out = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
                                    free, NULL);
    // 8192 mirrors tui_find_d's buffer (and the writer's "%8191s" upper
    // bound) -- the two sides need to move together if either gets bumped.
    char buf[8192];
    for (int64_t i = 1; i <= (int64_t)n_d; i++) {
        if (!oneGoto(of, 'd', i)) { stHash_destruct(out); oneFileClose(of); return NULL; }
        if (oneReadLine(of) != 'd') { stHash_destruct(out); oneFileClose(of); return NULL; }
        int64_t n = oneLen(of);
        if (n < 0 || n >= (int64_t)sizeof(buf)) {
            stHash_destruct(out); oneFileClose(of); return NULL;
        }
        memcpy(buf, oneString(of), (size_t)n);
        buf[n] = '\0';
        int64_t slen = oneInt(of, 2);
        // d-records are name-sorted and unique, so no need to check for
        // existing entry -- and the void* value path here means stHash
        // wouldn't free the old value anyway.
        stHash_insert(out, stString_copy(buf), (void *)(intptr_t)slen);
    }
    oneFileClose(of);
    return out;
}

TuiGenomeInfo *tui_genome_names(const char *tui_path, int64_t *n_out) {
    if (n_out) *n_out = 0;
    if (tui_path == NULL) return NULL;
    OneFile *of = oneFileOpenRead(tui_path, NULL, "tui", 1);
    if (of == NULL) return NULL;
    // Backward-compat: old .tui files predating the g-record schema lack
    // the 'g' line type entirely; ONElib leaves of->info['g']==NULL in
    // that case.  Return NULL so callers can fall back to a heuristic.
    if (of->info[(int)'g'] == NULL) { oneFileClose(of); return NULL; }
    I64 n_g = 0;
    oneStats(of, 'g', &n_g, NULL, NULL);
    if (n_g <= 0) { oneFileClose(of); return NULL; }
    TuiGenomeInfo *arr = (TuiGenomeInfo *) st_calloc((size_t)n_g, sizeof(TuiGenomeInfo));
    char buf[8192];
    for (int64_t i = 1; i <= (int64_t)n_g; i++) {
        if (!oneGoto(of, 'g', i)) { tui_genome_info_free(arr, n_g); oneFileClose(of); return NULL; }
        if (oneReadLine(of) != 'g') { tui_genome_info_free(arr, n_g); oneFileClose(of); return NULL; }
        int64_t nl = oneLen(of);
        if (nl < 0 || nl >= (int64_t)sizeof(buf)) {
            tui_genome_info_free(arr, n_g); oneFileClose(of); return NULL;
        }
        memcpy(buf, oneString(of), (size_t)nl);
        buf[nl] = '\0';
        arr[i - 1].name     = stString_copy(buf);
        arr[i - 1].total_bp = oneInt(of, 1);
        arr[i - 1].n_chroms = oneInt(of, 2);
    }
    oneFileClose(of);
    if (n_out) *n_out = n_g;
    return arr;
}

void tui_genome_info_free(TuiGenomeInfo *info, int64_t n) {
    if (info == NULL) return;
    for (int64_t i = 0; i < n; i++) free(info[i].name);
    free(info);
}

TuiGenomeLift *tui_genome_lift_load(Tui *tui, const char *genome_name) {
    if (tui == NULL || genome_name == NULL || *genome_name == 0) return NULL;

    // Use a dedicated OneFile handle (not tui->of): this one persists on the
    // returned gl for lazy chunk_decode on first column-query hit, and would
    // otherwise have its file position clobbered by interleaved tui_query
    // calls in the source-side path.  Closed in tui_genome_lift_destruct.
    OneFile *of = oneFileOpenRead(tui->path, NULL, "tui", 1);
    if (of == NULL) return NULL;

    // Match d-lines on prefix "<genome_name>." -- the writer's d-line names
    // are always "<genome>.<sequence>".
    size_t gn = strlen(genome_name);
    char *prefix = st_malloc(gn + 2);
    memcpy(prefix, genome_name, gn);
    prefix[gn] = '.';
    prefix[gn + 1] = 0;
    size_t plen = gn + 1;

    int64_t first = tui_find_d_lower_bound(of, tui->n_d, prefix);
    if (first > tui->n_d) { free(prefix); oneFileClose(of); return NULL; }

    // Pass 1: collect (S-ord, name, seqlen) for all matching d-lines.
    // d-lines are name-sorted so a prefix match is a contiguous range; stop
    // at the first non-match.
    typedef struct { int64_t ord; char *name; } DEnt;
    int64_t cap = 16, n = 0;
    DEnt *ents = st_malloc(cap * sizeof(DEnt));
    char buf[8192];
    for (int64_t i = first; i <= tui->n_d; i++) {
        if (!oneGoto(of, 'd', i)) break;
        if (oneReadLine(of) != 'd') break;
        int64_t bn = oneLen(of);
        if (bn < 0 || bn >= (int64_t)sizeof(buf)) break;
        memcpy(buf, oneString(of), (size_t)bn);
        buf[bn] = 0;
        if (strncmp(buf, prefix, plen) != 0) break;
        if (n == cap) { cap *= 2; ents = st_realloc(ents, cap * sizeof(DEnt)); }
        ents[n].ord = oneInt(of, 1);
        // d-line field 2 (seqlen) is unused by this lookup path; we only
        // need the S-ordinal + sequence name to walk the chunk metadata.
        // Store just the sequence part (skip the "<genome>." prefix) so the
        // returned wig records carry the bare contig name.
        ents[n].name = stString_copy(buf + plen);
        n++;
    }
    free(prefix);

    if (n == 0) { free(ents); oneFileClose(of); return NULL; }

    // Pass 2: per-S walk -- for each of our target's sequences, oneGoto its
    // S, then iterate (C, R)+ pairs reading only the C (oneGoto past R via
    // its next-C ordinal).  Bypasses the global C list (which can be
    // ~hundreds of thousands at vertebrate scale -- this targets only our
    // few-tens to few-thousand seqs).  R blobs stay on disk; chunk_decode
    // pulls them lazily on first column-query hit.
    int64_t cap_c = 1024, nc = 0;
    TGLChunk *chunks = st_malloc(cap_c * sizeof(TGLChunk));
    for (int64_t k_seq = 0; k_seq < n; k_seq++) {
        int64_t s_ord = ents[k_seq].ord + 1;
        if (!oneGoto(of, 'S', s_ord)) continue;
        if (oneReadLine(of) != 'S') continue;
        char c;
        while ((c = oneReadLine(of)) == 'C') {
            // The writer always pairs each C with the same seq's R+ chunks,
            // so the parent S-ord on a C should match the seq we jumped to.
            // A mismatch means we've walked into the next seq's chunks --
            // bail to the next outer iteration.  c_ord (the 4th field) is
            // the writer-stamped file-order ordinal; we use it instead of
            // ONElib's oneObject() since the accumulator is unreliable
            // after a cross-type oneGoto in this lib version.
            if (oneInt(of, 2) != s_ord) break;
            int64_t cur_c_ord = oneInt(of, 3);
            if (nc == cap_c) { cap_c *= 2; chunks = st_realloc(chunks, cap_c * sizeof(TGLChunk)); }
            chunks[nc].g_min          = oneInt(of, 0);
            chunks[nc].g_max          = oneInt(of, 1);
            chunks[nc].c_ord          = cur_c_ord;
            chunks[nc].seq_idx        = k_seq;
            chunks[nc].runs           = NULL;
            chunks[nc].n_runs         = 0;
            chunks[nc].max_end_prefix = NULL;
            nc++;
            // Skip the paired R by oneGoto-ing to the next C in the file.
            // Avoids reading R's deflated bytes -- only the C metadata
            // crosses the disk at open; R waits for the first column hit.
            if (!oneGoto(of, 'C', cur_c_ord + 1)) break;
        }
        (void)c;
    }

    if (nc == 0) {
        // Target has no chunks (empty sequences or unknown genome).
        free(chunks);
        for (int64_t i = 0; i < n; i++) free(ents[i].name);
        free(ents);
        oneFileClose(of);
        return NULL;
    }
    qsort(chunks, nc, sizeof(TGLChunk), tglchunk_cmp_gmin);

    // Running max(g_max) over chunks[0..i] for the outer scan's early exit.
    int64_t *cme = st_malloc(nc * sizeof(int64_t));
    int64_t running = chunks[0].g_max;
    cme[0] = running;
    for (int64_t i = 1; i < nc; i++) {
        if (chunks[i].g_max > running) running = chunks[i].g_max;
        cme[i] = running;
    }

    TuiGenomeLift *gl = st_calloc(1, sizeof(TuiGenomeLift));
    gl->n_seq = n;
    gl->seq_names = st_malloc(n * sizeof(char*));
    for (int64_t i = 0; i < n; i++) gl->seq_names[i] = ents[i].name;
    free(ents);
    gl->chunks        = chunks;
    gl->n_chunks      = nc;
    gl->chunk_max_end = cme;
    gl->of            = of;     // kept open; destruct closes
    return gl;
}

void tui_genome_lift_destruct(TuiGenomeLift *gl) {
    if (gl == NULL) return;
    for (int64_t i = 0; i < gl->n_seq; i++) free(gl->seq_names[i]);
    free(gl->seq_names);
    for (int64_t i = 0; i < gl->n_chunks; i++) {
        free(gl->chunks[i].runs);
        free(gl->chunks[i].max_end_prefix);
    }
    free(gl->chunks);
    free(gl->chunk_max_end);
    if (gl->of != NULL) oneFileClose(gl->of);
    free(gl);
}

// Lazy-decode one chunk: oneGoto its C, read past C into the paired R, decode
// the run triples, build GLRun array, build per-chunk max_end_prefix.
// No-op if already decoded.  The writer pre-sorts each chunk's runs by
// g_start, so the decoded triples come out g-sorted; no reader-side sort
// needed.  Caller must serialize calls on the same gl (the OneFile handle
// is not thread-safe).
static void chunk_decode(TGLChunk *ch, OneFile *of) {
    if (ch->runs != NULL) return;
    if (!oneGoto(of, 'C', ch->c_ord)) return;
    char c = oneReadLine(of);                  // C
    if (c != 'C') return;
    c = oneReadLine(of);                       // R (paired)
    if (c != 'R') return;
    int64_t raw_len = oneInt(of, 0);
    int64_t def_len = oneLen(of);
    const uint8_t *def = (const uint8_t *)oneString(of);
    int64_t bcap = raw_len + 3;
    int64_t *tmp = st_malloc((bcap ? bcap : 1) * sizeof(int64_t));
    int64_t m = decode_runs(def, def_len, raw_len, tmp, bcap);
    if (m == 0) { free(tmp); ch->n_runs = 0; ch->runs = (GLRun *)st_calloc(1, 1); return; }
    GLRun *rs = st_malloc((size_t)m * sizeof(GLRun));
    for (int64_t i = 0; i < m; i++) {
        int64_t t = tmp[3*i+0], g = tmp[3*i+1], lenc = tmp[3*i+2];
        rs[i].g_start = g;
        rs[i].t_start = t;
        rs[i].length  = lenc >> 1;
        rs[i].strand  = (int)(1 - (lenc & 1));
        rs[i].seq_idx = ch->seq_idx;
        // Sanity: the writer pre-sorts each chunk's runs by g_start.  If
        // this assert fires we're either reading an old .tui (built by
        // the t-sorted writer) or the format invariant got broken some
        // other way -- chunk_collect's binary search + max_end_prefix
        // both require g_start to be non-decreasing, so a violation
        // would silently miss matches downstream.
        assert(i == 0 || g >= rs[i - 1].g_start);
    }
    free(tmp);
    int64_t *mep = st_malloc((size_t)m * sizeof(int64_t));
    int64_t running = rs[0].g_start + rs[0].length;
    mep[0] = running;
    for (int64_t i = 1; i < m; i++) {
        int64_t e = rs[i].g_start + rs[i].length;
        if (e > running) running = e;
        mep[i] = running;
    }
    ch->runs           = rs;
    ch->n_runs         = m;
    ch->max_end_prefix = mep;
}

// Search within one decoded chunk for runs covering `column`; append to out[].
// Mirrors the legacy backward-scan-with-prefix-bound; bounded by per-chunk
// paralog density.
static int chunk_collect(const TGLChunk *ch, int64_t column,
                         TuiGenomeMatch *out, int cap, int already, const char *seq_name) {
    if (ch->n_runs == 0) return already;
    int64_t lo = 0, hi = ch->n_runs;
    while (lo < hi) {
        int64_t m = lo + (hi - lo) / 2;
        if (ch->runs[m].g_start <= column) lo = m + 1; else hi = m;
    }
    int64_t j = lo - 1;
    int count = already;
    while (j >= 0) {
        const GLRun *r = &ch->runs[j];
        if (column < r->g_start + r->length) {
            if (count < cap && out != NULL) {
                int64_t offset = column - r->g_start;
                out[count].seq    = seq_name;
                out[count].pos    = r->strand
                    ? r->t_start + offset
                    : r->t_start + r->length - 1 - offset;
                out[count].strand = (r->strand != 0);
            }
            count++;
        }
        if (j == 0) break;
        if (ch->max_end_prefix[j - 1] <= column) break;
        j--;
    }
    return count;
}

// Find chunks intersecting `column`, lazy-decode them on first hit, search
// within.  The outer scan terminates via chunk_max_end[i-1] <= column (same
// early-exit shape as the per-chunk max_end_prefix used in chunk_collect,
// just lifted to chunks).
int tui_genome_lift_column(const TuiGenomeLift *gl, int64_t column,
                           TuiGenomeMatch *out, int cap) {
    if (gl == NULL || gl->n_chunks == 0) return 0;
    // chunks are logically mutable for lazy decode -- the gl pointer is
    // const-by-API but the underlying entries get filled in on first hit.
    TGLChunk *cs = gl->chunks;
    int64_t lo = 0, hi = gl->n_chunks;
    while (lo < hi) {
        int64_t m = lo + (hi - lo) / 2;
        if (cs[m].g_min <= column) lo = m + 1; else hi = m;
    }
    int64_t j = lo - 1;
    int count = 0;
    while (j >= 0) {
        TGLChunk *ch = &cs[j];
        if (column < ch->g_max) {
            if (ch->runs == NULL) chunk_decode(ch, gl->of);
            count = chunk_collect(ch, column, out, cap, count,
                                  gl->seq_names[ch->seq_idx]);
        }
        if (j == 0) break;
        if (gl->chunk_max_end[j - 1] <= column) break;
        j--;
    }
    return count;
}

int64_t tui_genome_lift_n_chunks(const TuiGenomeLift *gl) {
    return gl == NULL ? 0 : gl->n_chunks;
}

// Visit runs in `ch` intersecting [c_lo, c_hi).  Binary-search for the first
// run with g_start >= c_lo, then walk backward (catching earlier runs whose
// length extends past c_lo via max_end_prefix) and forward (until g_start
// >= c_hi).  Caller must have decoded ch->runs.
static void visit_chunk_runs(const TGLChunk *ch, const TuiGenomeLift *gl,
                             int64_t c_lo, int64_t c_hi,
                             void (*cb)(const TuiRun *run, void *user),
                             void *user) {
    if (ch->n_runs == 0) return;
    int64_t lo = 0, hi = ch->n_runs;
    while (lo < hi) {
        int64_t m = lo + (hi - lo) / 2;
        if (ch->runs[m].g_start < c_lo) lo = m + 1; else hi = m;
    }
    // Backward pass: earlier runs that extend into [c_lo, c_hi).
    for (int64_t j = lo - 1; j >= 0; j--) {
        const GLRun *r = &ch->runs[j];
        int64_t r_end = r->g_start + r->length;
        if (r_end > c_lo && r->g_start < c_hi) {
            TuiRun out = { gl->seq_names[r->seq_idx],
                           r->g_start, r->length, r->t_start, r->strand };
            cb(&out, user);
        }
        if (j == 0) break;
        if (ch->max_end_prefix[j - 1] <= c_lo) break;
    }
    // Forward pass: g_start >= c_lo by binary-search result.
    for (int64_t j = lo; j < ch->n_runs; j++) {
        const GLRun *r = &ch->runs[j];
        if (r->g_start >= c_hi) break;
        TuiRun out = { gl->seq_names[r->seq_idx],
                       r->g_start, r->length, r->t_start, r->strand };
        cb(&out, user);
    }
}

void tui_genome_lift_visit_runs(TuiGenomeLift *gl, int64_t c_lo, int64_t c_hi,
                                void (*cb)(const TuiRun *run, void *user),
                                void *user) {
    if (gl == NULL || cb == NULL || c_lo >= c_hi) return;
    TGLChunk *cs = gl->chunks;
    // Chunks are sorted by g_min.  Binary-search the lower_bound on
    // g_min >= c_lo, then collect straddlers (g_min < c_lo but g_max
    // > c_lo) via a BACKWARD pass terminated by chunk_max_end[j-1]
    // <= c_lo, and in-range chunks via a FORWARD pass terminated by
    // g_min >= c_hi.  This mirrors tui_genome_lift_column's per-column
    // outer scan.
    //
    // Pass 1 (cheap, serial): collect visited chunk indices into
    // visited[].  Pass 2 (parallel): chunk_decode each chunk that
    // hasn't been decoded yet -- the zlib-inflate of a ~1.6 MB R blob
    // dominates wall on cold cache, and each chunk is decoded by
    // exactly one thread thanks to OpenMP's iteration partitioning
    // (no race on ch->runs).  Pass 3 (serial): walk runs and call cb
    // -- pending_push / bin accumulators are not thread-safe, and
    // the per-run cost is microseconds anyway so this is cheap.
    int64_t lo = 0, hi = gl->n_chunks;
    while (lo < hi) {
        int64_t m = lo + (hi - lo) / 2;
        if (cs[m].g_min < c_lo) lo = m + 1; else hi = m;
    }
    int64_t start = lo;

    int64_t v_cap = 64, v_n = 0;
    int64_t *visited = st_malloc((size_t)v_cap * sizeof(int64_t));
    // Backward pass collects in DESCENDING ci order; reverse below so
    // visited[] is g_min-ascending (preserves cb call order vs the
    // previous in-place implementation, which mattered to downstream
    // pending_push merging).
    int64_t back_start = v_n;
    for (int64_t ci = start - 1; ci >= 0; ci--) {
        if (cs[ci].g_max > c_lo) {
            if (v_n == v_cap) {
                v_cap *= 2;
                visited = st_realloc(visited, (size_t)v_cap * sizeof(int64_t));
            }
            visited[v_n++] = ci;
        }
        if (ci == 0) break;
        if (gl->chunk_max_end[ci - 1] <= c_lo) break;
    }
    for (int64_t i = back_start, j = v_n - 1; i < j; i++, j--) {
        int64_t t = visited[i]; visited[i] = visited[j]; visited[j] = t;
    }
    for (int64_t ci = start; ci < gl->n_chunks; ci++) {
        if (cs[ci].g_min >= c_hi) break;
        if (v_n == v_cap) {
            v_cap *= 2;
            visited = st_realloc(visited, (size_t)v_cap * sizeof(int64_t));
        }
        visited[v_n++] = ci;
    }

    // Walk visited chunks in g_min-ascending order, lazy-decoding any
    // that haven't been hit by a prior call.  Cached chunks (ch->runs
    // populated by a previous visit) skip the inflate.
    for (int64_t i = 0; i < v_n; i++) {
        TGLChunk *ch = &cs[visited[i]];
        if (ch->runs == NULL) chunk_decode(ch, gl->of);
        visit_chunk_runs(ch, gl, c_lo, c_hi, cb, user);
    }
    free(visited);
}

void tui_genome_lift_stream_runs(TuiGenomeLift *gl,
                                 void (*cb)(const TuiRun *run, void *user),
                                 void *user) {
    if (gl == NULL || cb == NULL) return;
    /* One-shot scan: visit each chunk in g_min order, decode-visit-free.
     * Unlike tui_genome_lift_visit_runs (which caches decoded chunks
     * forever on the assumption that browser pan/zoom will hit them
     * again), this trades cache reuse for bounded memory -- only one
     * chunk's decoded runs are resident at a time (~3 MB peak).  For
     * whole-genome scanning workloads (taf_coarsen) where each chunk
     * is touched exactly once and the cache would only inflate peak
     * RSS to many tens of GB on giant ancestors. */
    for (int64_t ci = 0; ci < gl->n_chunks; ci++) {
        TGLChunk *ch = &gl->chunks[ci];
        int was_cached = (ch->runs != NULL);
        if (!was_cached) chunk_decode(ch, gl->of);
        visit_chunk_runs(ch, gl, ch->g_min, ch->g_max + 1, cb, user);
        if (!was_cached) {
            free(ch->runs);          ch->runs          = NULL;
            free(ch->max_end_prefix); ch->max_end_prefix = NULL;
            ch->n_runs = 0;
        }
        /* If the chunk was already cached before this call, leave it
         * alone -- a concurrent / earlier caller (e.g. a browser query)
         * may rely on it. */
    }
}

/////////////////////////////////////////////////////////////////////////////
// Universal-column block extractor.  Single forward scan in column order from
// the X-index anchor; reuses tai's seek/resync/read primitives.  No per-
// contig key => the C1 (non-monotone row-0 ancestor) failure is impossible.
/////////////////////////////////////////////////////////////////////////////

struct _TuiExtractIt {
    Tui *tui;               // borrowed; owns idxCol/idxFpos used for per-iv seeks
    int is_maf;
    bool rle;
    Alignment *(*readblk)(Alignment *, bool, LI *);
    TuiInterval *iv;        // owned copy, sorted+merged universal-col intervals
    int64_t n_iv, iv_cur;
    int64_t c_hi;           // last interval end (stop bound)
    int64_t scan_col;       // universal column of the next block to read
    Alignment *phys;        // current physical block overlapping iv (owned)
    int64_t phys_col;       // universal column of phys's first column
    int64_t *runs;          // covered local-column runs of phys: [j0,j1) pairs
    int64_t runs_cap, n_runs, run_i;
    Alignment *to_free;     // sub-block returned last call (freed next/destruct)
    int64_t last_col_start; // universal col of the first col of `to_free`
};

// Build a standalone Alignment = columns [j0,j1) of `src`, every row clipped.
// Reuses .tai clip_alignment's convention: a row's start advances by the
// non-gap bases removed on the left, length = the non-gap bases kept (works
// for '+' and '-' identically, as clip_alignment does).  src is unchanged so
// multiple runs of one physical block stay independent.
//
// Note: emitted rows have l_row/r_row implicitly NULL (st_calloc).  This is
// what makes the link-skip in tai_maf_read_block_nolink safe: nothing the
// .tui extractor emits depends on adjacent-block linkage.  If you ever start
// using l_row/r_row downstream of this function, you'll need to call
// alignment_link_adjacent on the emitted chain yourself (or wire the
// extractor back onto tai_maf_read_block).
static Alignment *tui_subblock(const Alignment *src, int64_t j0, int64_t j1) {
    int64_t newcol = j1 - j0;
    Alignment *aln = st_calloc(1, sizeof(Alignment));
    aln->column_number = newcol;
    aln->column_tags = st_calloc(newcol ? newcol : 1, sizeof(Tag *));
    aln->row_number = src->row_number;
    Alignment_Row **pp = &aln->row;
    for (Alignment_Row *r = src->row; r != NULL; r = r->n_row) {
        int64_t pre = 0, mid = 0;
        for (int64_t c = 0; c < j0; c++) if (r->bases[c] != '-') pre++;
        for (int64_t c = j0; c < j1; c++) if (r->bases[c] != '-') mid++;
        Alignment_Row *nr = st_calloc(1, sizeof(Alignment_Row));
        nr->sequence_name = stString_copy(r->sequence_name);
        nr->start = r->start + pre;
        nr->length = mid;
        nr->strand = r->strand;
        nr->sequence_length = r->sequence_length;
        nr->bases = stString_getSubString(r->bases, j0, newcol);
        *pp = nr;
        pp = &nr->n_row;
    }
    return aln;
}

// Covered local-column runs of the physical block at universal [pc, pc+cn).
// Fills it->runs; advances it->iv_cur only past intervals fully consumed by
// this block (an interval extending beyond pc+cn stays for the next block).
static void tui_compute_runs(TuiExtractIt *it, int64_t pc, int64_t cn) {
    it->n_runs = 0;
    int64_t pend = pc + cn;
    while (it->iv_cur < it->n_iv && it->iv[it->iv_cur].end <= pc) it->iv_cur++;
    int64_t k = it->iv_cur;
    while (k < it->n_iv && it->iv[k].start < pend) {
        int64_t s = it->iv[k].start > pc ? it->iv[k].start : pc;
        int64_t e = it->iv[k].end < pend ? it->iv[k].end : pend;
        if (s < e) {
            if (2 * (it->n_runs + 1) > it->runs_cap) {
                it->runs_cap = it->runs_cap ? it->runs_cap * 2 : 16;
                it->runs = st_realloc(it->runs, it->runs_cap * sizeof(int64_t));
            }
            it->runs[2*it->n_runs]   = s - pc;
            it->runs[2*it->n_runs+1] = e - pc;
            it->n_runs++;
        }
        if (it->iv[k].end <= pend) { k++; it->iv_cur = k; }   // fully consumed
        else break;                                           // continues next block
    }
}

// Greatest anchor index `a` with idxCol[a] <= c (idxCol strictly increasing,
// always has idxCol[0]=0).  O(log N).
static int64_t tui_anchor_for(const Tui *tui, int64_t c) {
    int64_t lo = 0, hi = tui->idxN - 1, a = 0;
    while (lo <= hi) {
        int64_t mid = (lo + hi) / 2;
        if (tui->idxCol[mid] <= c) { a = mid; lo = mid + 1; } else hi = mid - 1;
    }
    return a;
}

// Seek the underlying stream to anchor `a` and mirror the (TAF) resync the
// existing tai_iterator uses, so the next readblk() works regardless of
// format.  Returns the universal column we're now positioned at.
static int64_t tui_seek_to_anchor(LI *li, const Tui *tui, int64_t a, int is_maf) {
    LI_seek(li, tui->idxFpos[a]);
    LI_get_next_line(li);                       // mirror tai_iterator
    if (!is_maf) tai_resync_taf_line(LI_peek_at_next_line(li));
    return tui->idxCol[a];
}

// Forward-scan (ctx = taf context, not freed here) skipping non-overlapping
// physical blocks, until one overlaps iv -> it->phys/phys_col/runs(run_i=0),
// or past the last interval / EOF -> it->phys = NULL.
//
// Per-interval seek: at the top of each iteration we check whether the X
// anchor before iv[iv_cur].start is meaningfully ahead of scan_col -- if so,
// we LI_seek there before reading the next block, bounding the forward scan
// per emitted block to at most one X-anchor gap (~10000 cols).  Eliminates
// the O(c_hi - c_lo) blowup when lifted intervals are scattered in column
// space (was costing minutes for kilobase-scale leaf-coordinate queries).
static void tui_extract_advance(TuiExtractIt *it, LI *li, Alignment *ctx) {
    Alignment *prev = ctx;
    it->phys = NULL;
    while (1) {
        if (it->iv_cur < it->n_iv) {
            int64_t a = tui_anchor_for(it->tui, it->iv[it->iv_cur].start);
            if (it->tui->idxCol[a] > it->scan_col) {
                if (prev != NULL && prev != ctx) alignment_destruct(prev, 1);
                it->scan_col = tui_seek_to_anchor(li, it->tui, a, it->is_maf);
                prev = NULL;                            // TAF chain broken by seek
                continue;
            }
        }
        Alignment *b = it->readblk(prev, it->rle, li);
        if (b == NULL) break;
        int64_t cn = b->column_number;
        int64_t bstart = it->scan_col, bend = bstart + cn;
        while (it->iv_cur < it->n_iv && it->iv[it->iv_cur].end <= bstart) it->iv_cur++;
        if (it->iv_cur < it->n_iv && it->iv[it->iv_cur].start < bend) {
            it->phys = b; it->phys_col = bstart; it->run_i = 0;
            tui_compute_runs(it, bstart, cn);
            break;
        }
        it->scan_col = bend;
        if (it->iv_cur >= it->n_iv || it->scan_col >= it->c_hi) {
            alignment_destruct(b, 1); break;            // past the last interval
        }
        // prev can be NULL here (we just seeked); the original two-way check
        // didn't have to consider that.  Guard explicitly.
        if (prev != NULL && prev != ctx) alignment_destruct(prev, 1);
        prev = b;
    }
    if (prev != NULL && prev != ctx && prev != it->phys) alignment_destruct(prev, 1);
}

TuiExtractIt *tui_extract_iterator(Tui *tui, LI *li, int is_maf, bool rle,
                                   const TuiInterval *iv, int64_t n_iv) {
    TuiExtractIt *it = st_calloc(1, sizeof(TuiExtractIt));
    it->tui = tui;
    it->is_maf = is_maf;
    it->rle = rle;
    // _nolink: the extractor emits whole MAF blocks and never reads the
    // per-row l_row/r_row neighbour pointers that alignment_link_adjacent
    // would populate.  Skipping that step removes ~19% of CPU on the fish
    // 1 Mb profile (the WFA chain underneath).  taf_read_block doesn't run
    // alignment_link_adjacent itself, so the TAF side stays as-is.
    it->readblk = is_maf ? tai_maf_read_block_nolink : taf_read_block;
    if (n_iv <= 0 || tui->idxN == 0) return it;        // empty iterator
    it->iv = st_malloc((size_t)n_iv * sizeof(TuiInterval));
    memcpy(it->iv, iv, (size_t)n_iv * sizeof(TuiInterval));
    it->n_iv = n_iv;
    it->c_hi = iv[n_iv - 1].end;
    // Position scan_col before the first interval; tui_extract_advance will
    // do its first per-interval seek and start reading from the right anchor.
    it->scan_col = -1;
    it->iv_cur = 0;
    tui_extract_advance(it, li, NULL);
    return it;
}

Alignment *tui_extract_next(TuiExtractIt *it, LI *li) {
    if (it->to_free != NULL) { alignment_destruct(it->to_free, 1); it->to_free = NULL; }
    if (it->phys == NULL) return NULL;
    int64_t j0 = it->runs[2*it->run_i], j1 = it->runs[2*it->run_i + 1];
    Alignment *sub = tui_subblock(it->phys, j0, j1);
    it->last_col_start = it->phys_col + j0;             // captured BEFORE advance
    it->run_i++;
    if (it->run_i >= it->n_runs) {                      // this physical block done
        Alignment *old = it->phys;
        it->scan_col = it->phys_col + old->column_number;
        tui_extract_advance(it, li, old);               // old = taf context
        alignment_destruct(old, 1);                     // context consumed; subs copied
    }
    it->to_free = sub;                                  // caller uses before next call
    return sub;
}

int64_t tui_extract_col_start(const TuiExtractIt *it) { return it->last_col_start; }

bool tui_extract_has_next(TuiExtractIt *it) { return it->phys != NULL; }

void tui_extract_take_ownership(TuiExtractIt *it) {
    // Disown the most recent yield so the next _next() / _destruct doesn't
    // free it.  Caller takes responsibility for `alignment_destruct`.
    it->to_free = NULL;
}

void tui_extract_iterator_destruct(TuiExtractIt *it) {
    if (it->to_free != NULL) alignment_destruct(it->to_free, 1);
    if (it->phys != NULL) alignment_destruct(it->phys, 1);
    free(it->runs);
    free(it->iv);
    free(it);
}
