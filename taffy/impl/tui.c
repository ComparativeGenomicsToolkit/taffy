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
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <signal.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
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
    "D R 2 3 INT 6 STRING\n";                // runs: inflatedLen, deflate(blob)

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

static void tui_atexit_cleanup(void) {
    if (atexit_spill_paths != NULL) {
        int64_t n = stList_length(atexit_spill_paths);
        for (int64_t i = 0; i < n; i++) {
            (void) remove((char *) stList_get(atexit_spill_paths, i));
        }
        stList_destruct(atexit_spill_paths);
        atexit_spill_paths = NULL;
    }
    if (atexit_tui_path != NULL) {
        (void) remove(atexit_tui_path);
        free(atexit_tui_path);
        atexit_tui_path = NULL;
    }
}

static void tui_signal_handler(int sig) {
    tui_atexit_cleanup();
    // restore default disposition and re-raise so the parent sees the right
    // exit (don't swallow SIGINT into exit-code-0)
    signal(sig, SIG_DFL);
    raise(sig);
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
static void tui_atexit_track_tui(const char *path) {
    if (atexit_tui_path != NULL) { free(atexit_tui_path); }
    atexit_tui_path = stString_copy(path);
}

// Disarm the crash cleanup on the normal success path, BEFORE tui_cleanup()
// frees its path strings.  Idempotent.
static void tui_atexit_disarm(void) {
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

// One spilled run, text line: "seqName\tt_start\tg_start\tlength\tstrand\n"
static void spill_run(FILE *fh, const char *seq_name, int64_t t_start,
                      int64_t g_start, int64_t length, bool strand) {
    if (fprintf(fh, "%s\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%c\n",
                seq_name, t_start, g_start, length, strand ? '+' : '-') < 0) {
        st_errAbort("tui: failed writing run spill (disk full / write error)");
    }
}

// Emit the maximal gap-free stretches of one row of one block.  `g` is the
// global column of this block's first column.
static void emit_row_runs(FILE *spill, Alignment_Row *row, int64_t g,
                          int64_t column_number) {
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
        spill_run(spill, row->sequence_name, t_start, g_start, slen, row->strand);
        pre += slen;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Phase 1: one streaming scan, spill runs per-genome (column order, O(1) RAM)
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    stHash *spill_fh;     // genomeName -> FILE*  (kept open through phase 1)
    stHash *spill_path;   // genomeName -> char*  (we own/free these)
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

static FILE *spill_for(Phase1 *p1, const char *genome) {
    FILE *fh = stHash_search(p1->spill_fh, (void *)genome);
    if (fh == NULL) {
        char *path = stString_print("%s/%s.tuiSpill.%ld.%d", p1->tmp_dir,
                                    p1->out_base, (long)getpid(),
                                    p1->next_spill_id++);
        tui_atexit_track_spill(path);
        fh = fopen(path, "w");
        if (fh == NULL) { fprintf(stderr, "tui: cannot open spill %s\n", path); exit(1); }
        stHash_insert(p1->spill_fh, stString_copy(genome), fh);
        stHash_insert(p1->spill_path, stString_copy(genome), path);
    }
    return fh;
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
        emit_row_runs(spill_for(p1, gname), row, p1->T, cn);
        free(gname);
    }
    p1->T += cn;
}

/////////////////////////////////////////////////////////////////////////////
// Phase 2: per-genome in-RAM sort by (seq, t_start), write ONEcode container
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

// Load every run line of one genome spill into a Run array (caller frees).
static Run *load_genome_runs(const char *path, int64_t *n_out) {
    FILE *fh = fopen(path, "r");
    if (fh == NULL) { *n_out = 0; return NULL; }
    int64_t cap = 1024, n = 0;
    Run *runs = st_malloc(cap * sizeof(Run));
    char line[16384];
    while (fgets(line, sizeof(line), fh) != NULL) {
        if (n == cap) { cap *= 2; runs = st_realloc(runs, cap * sizeof(Run)); }
        Run *r = &runs[n];
        char sbuf[8192], st;
        if (sscanf(line, "%8191s\t%" SCNi64 "\t%" SCNi64 "\t%" SCNi64 "\t%c",
                   sbuf, &r->t, &r->g, &r->len, &st) == 5) {
            r->seq = stString_copy(sbuf);
            r->strand = st;
            n++;
        } else {
            st_errAbort("tui: corrupt spill line in %s: %s", path, line);
        }
    }
    fclose(fh);
    *n_out = n;
    return runs;
}

/////////////////////////////////////////////////////////////////////////////
// Per-sequence run codec: zigzag-delta + LEB128 varint, then zlib deflate.
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
static uint8_t *encode_runs(const int64_t *buf, int64_t m,
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
static uint8_t *encode_idx(const int64_t *col, const int64_t *fpos, int64_t n,
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
    stHashIterator *hit = stHash_getIterator(p1->spill_path);
    char *gk;
    while ((gk = stHash_getNext(hit)) != NULL) {
        remove((char *)stHash_search(p1->spill_path, gk));
    }
    stHash_destructIterator(hit);
    if (p1->idx_fh != NULL) fclose(p1->idx_fh);
    if (p1->idx_path != NULL) { remove(p1->idx_path); free(p1->idx_path); }
    for (int64_t i = 0; i < n_seqs; i++) free(seqks[i].genome);  // .seq owned by seq_len
    free(seqks);
    stHash_destruct(p1->spill_fh);
    stHash_destruct(p1->spill_path);
    stHash_destruct(p1->seq_len);
    stList_destruct(p1->seq_keys);
    if (tree_map != NULL) stHash_destruct(tree_map);
    free(eff_tmp);
}

int tui_create(LI *li, const char *out_path, const char *tmp_dir,
                stHash *genome_name_map) {
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

    Phase1 p1;
    p1.tmp_dir = eff_tmp;
    p1.out_base = base_of(out_path);
    p1.next_spill_id = 0;
    p1.gmap = eff_map;
    p1.spill_fh = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    // we own the spill paths -> free them via the hash on destruct
    p1.spill_path = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    p1.seq_len = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    p1.seq_keys = stList_construct();
    p1.T = 0;
    p1.idx_fh = NULL; p1.idx_path = NULL; p1.idx_last_col = 0; p1.idx_n = 0;

    // idx_anchor captures the block's anchor offset itself (LI_get_position =
    // start of the peeked block-first line; the exact tai_create_taf
    // convention so LI_seek + the shared readers behave identically at query).
    time_t t_total_start = time(NULL);
    time_t t_phase1_start = t_total_start;
    st_logInfo("tui: starting phase 1 (streaming %s scan, spilling per-genome runs)\n",
               input_format == 0 ? "TAF" : "MAF");
    int64_t block_count = 0;
    Alignment *aln, *p_aln = NULL;
    if (input_format == 1) {
        while (1) {
            idx_anchor(&p1, li, 1, false);
            aln = maf_read_block(li);
            if (aln == NULL) break;
            tui_phase1_block(&p1, aln);
            alignment_destruct(aln, 1);
            block_count++;
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
        }
        if (p_aln != NULL) alignment_destruct(p_aln, 1);
    }
    // fclose surfaces buffered short writes (e.g. disk full) -> fail loudly
    // here rather than emit a silently-truncated spill that Phase 2 would
    // later misread as "corrupt spill line".
    if (p1.idx_fh != NULL) {
        if (fclose(p1.idx_fh) != 0) st_errAbort("tui: idx spill close failed "
            "(disk full / write error) -- %s", p1.idx_path);
        p1.idx_fh = NULL;
    }
    // close all spill FILE*s (paths kept for phase 2)
    stHashIterator *hit = stHash_getIterator(p1.spill_fh);
    char *gk;
    while ((gk = stHash_getNext(hit)) != NULL) {
        if (fclose(stHash_search(p1.spill_fh, gk)) != 0) {
            st_errAbort("tui: spill close failed for genome '%s' "
                        "(disk full / write error)", gk);
        }
    }
    stHash_destructIterator(hit);
    int64_t n_genomes = (int64_t) stHash_size(p1.spill_path);
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

    OneSchema *schema = oneSchemaCreateFromText(TUI_SCHEMA);
    OneFile *of = oneFileOpenWriteNew(out_path, schema, "tui", true, 1);
    if (of == NULL) {
        fprintf(stderr, "tui: cannot write %s\n", out_path);
        oneSchemaDestroy(schema);
        tui_cleanup(&p1, seqks, n_seqs, eff_tmp, tree_map);
        tui_atexit_disarm();
        return 1;
    }
    tui_atexit_track_tui(out_path);
    oneAddProvenance(of, "taffy", "tui", "universal column index", 0);

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
        uint8_t *xdef = encode_idx(ic, iff, in, &x_raw, &x_def);
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
    time_t t_phase2_start = time(NULL);
    st_logInfo("tui: starting phase 2 (per-genome sort + encode + write, "
               "%" PRIi64 " genomes)\n", n_genomes);
    int64_t genome_idx = 0;
    int64_t i = 0;
    while (i < n_seqs) {
        char *gname = seqks[i].genome;  // groups spill files; not emitted
        time_t t_g_start = time(NULL);

        int64_t nr = 0;
        Run *runs = NULL;
        char *path = stHash_search(p1.spill_path, gname);
        if (path != NULL) {
            runs = load_genome_runs(path, &nr);
            if (runs != NULL) qsort(runs, nr, sizeof(Run), run_cmp);
        }
        // walk this genome's seqs (contiguous block in seqks[]); runs[] are
        // seq-sorted consistently (strcmp), so advance a single forward cursor.
        int64_t rc = 0;
        while (i < n_seqs && strcmp(seqks[i].genome, gname) == 0) {
            char *sk = seqks[i].seq;
            int64_t slen = *(int64_t *)stHash_search(p1.seq_len, sk);
            oneInt(of, 1) = slen;
            oneWriteLine(of, 'S', strlen(sk), (void *)sk);
            // R: runs of this seq (contiguous in sorted runs[]); same strcmp
            // order as seqks[] so the forward cursor can't skip/misattribute
            while (rc < nr && strcmp(runs[rc].seq, sk) < 0) rc++;  // forward only
            int64_t a = rc;
            int64_t bnd = a;
            while (bnd < nr && strcmp(runs[bnd].seq, sk) == 0) bnd++;
            int64_t cnt = bnd - a;
            int64_t *buf = st_malloc((cnt ? 3 * cnt : 1) * sizeof(int64_t));
            // Colinear merge: runs[a..bnd) are this seq's runs sorted by t
            // ascending and partition its forward coords disjointly.  Two
            // consecutive runs describe the same linear p->col map and so fold
            // into one iff same strand, forward-adjacent (prev.t+prev.len ==
            // cur.t), AND column-contiguous.  col(p): '+' = g+(p-t) (slope +1,
            // so cur.g == prev.g+prev.len); '-' = g+(t+len-1-p) (slope -1, so
            // cur.g == prev.g-cur.len, and the merged run keeps cur.g).  An
            // in-row alignment gap leaves t adjacent but g discontinuous, so
            // the column test correctly refuses to merge across it.
            int64_t m = 0;  // merged-run count (triples written into buf)
            for (int64_t k = a; k < bnd; k++) {
                Run *r = &runs[k];
                int64_t rt = r->t, rg = r->g, rl = r->len;
                char rs = r->strand;
                if (m > 0) {
                    int64_t *pp = &buf[3 * (m - 1)];
                    int64_t pt = pp[0], pg = pp[1], pl = pp[2] >> 1;
                    char ps = (pp[2] & 1) ? '-' : '+';
                    if (ps == rs && rt == pt + pl) {
                        if (rs == '+' && rg == pg + pl) {
                            pp[2] = (pl + rl) << 1;            // extend, '+'
                            continue;
                        }
                        if (rs == '-' && rg == pg - rl) {
                            pp[1] = rg;                        // merged g = cur.g
                            pp[2] = ((pl + rl) << 1) | 1;      // extend, '-'
                            continue;
                        }
                    }
                }
                buf[3*m+0] = rt;
                buf[3*m+1] = rg;
                buf[3*m+2] = (rl << 1) | (rs == '-' ? 1 : 0);
                m++;
            }
            // Encode merged triples -> zigzag-delta varint -> deflate.
            int64_t raw_len = 0, def_len = 0;
            uint8_t *def = encode_runs(buf, m, &raw_len, &def_len);
            // Self-check: decode must reproduce the triples exactly.  taffy
            // builds with asserts on, so this validates the codec on every
            // file (incl. the real vertebrate-scale ones), losslessly.
            int64_t *chk = st_malloc((m ? 3 * m : 1) * sizeof(int64_t));
            int64_t dm = decode_runs(def, def_len, raw_len, chk, 3 * m);
            assert(dm == m && memcmp(chk, buf, (size_t)(3 * m) * sizeof(int64_t)) == 0);
            free(chk);
            oneInt(of, 0) = raw_len;
            oneWriteLine(of, 'R', def_len, def);
            free(def);
            free(buf);
            rc = bnd;   // advance cursor past this seq's runs
            i++;
        }
        st_logInfo("tui: phase 2 genome %" PRIi64 "/%" PRIi64 " '%s' done in "
                   "%" PRIi64 " seconds (%" PRIi64 " runs)\n",
                   ++genome_idx, n_genomes, gname,
                   (int64_t)(time(NULL) - t_g_start), nr);
        for (int64_t k = 0; k < nr; k++) free(runs[k].seq);
        free(runs);
        // gname is seqks[i].genome -- owned by seqks, freed in tui_cleanup
    }
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
    char   *path;       // .tui path (re-opened per query)
    int64_t T;          // total universal columns
    int64_t n_d;        // number of d-lines (binary-search upper bound)
    // Universal-column -> file_pos index (X track); both strictly increasing.
    int64_t *idxCol, *idxFpos;
    int64_t  idxN;
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
    oneFileClose(of);
    return tui;
}

void tui_destruct(Tui *tui) {
    if (tui == NULL) return;
    free(tui->idxCol); free(tui->idxFpos);
    free(tui->path);
    free(tui);
}

int64_t tui_total_columns(const Tui *tui) { return tui->T; }

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

int tui_has_sequence(const Tui *tui, const char *seq_name) {
    OneFile *of = oneFileOpenRead(tui->path, NULL, "tui", 1);
    if (of == NULL) return 0;
    int64_t ord = tui_find_d(of, tui->n_d, seq_name, NULL);
    oneFileClose(of);
    return ord >= 0 ? 1 : 0;
}

static int tui_iv_cmp(const void *a, const void *b) {
    const TuiInterval *x = a, *y = b;
    return (x->start < y->start) ? -1 : (x->start > y->start) ? 1 : 0;
}

TuiInterval *tui_query(Tui *tui, const char *seq_name,
                       int64_t start, int64_t end, int64_t *n_out) {
    *n_out = 0;
    if (start >= end) return NULL;

    OneFile *of = oneFileOpenRead(tui->path, NULL, "tui", 1);
    if (of == NULL) return NULL;

    // Resolve name -> S-ordinal by binary-searching the (name-sorted) d-lines
    // via oneGoto; O(log n_d) seeks, no preloaded directory hashes.
    int64_t ord = tui_find_d(of, tui->n_d, seq_name, NULL);
    if (ord < 0) { oneFileClose(of); return NULL; }

    // Jump straight to the (ord+1)-th S object via the ONElib footer index.
    // After the goto, oneReadLine() returns the S line, then the matching R
    // line follows immediately (writer emits S; R; S; R; ... per genome).
    if (!oneGoto(of, 'S', ord + 1)) { oneFileClose(of); return NULL; }
    int64_t *runs = NULL, m = 0;
    char c = oneReadLine(of);
    if (c != 'S') { oneFileClose(of); return NULL; }
    // Sanity check the S name matches what the directory said (would only
    // fail on a corrupted / wrongly-built .tui).
    int64_t sn = oneLen(of);
    if (sn != (int64_t)strlen(seq_name) ||
        memcmp(oneString(of), seq_name, sn) != 0) {
        oneFileClose(of);
        return NULL;
    }
    c = oneReadLine(of);
    if (c != 'R') { oneFileClose(of); return NULL; }
    int64_t raw_len = oneInt(of, 0);
    int64_t def_len = oneLen(of);
    const uint8_t *def = (const uint8_t *)oneString(of);
    int64_t cap = raw_len + 3;              // m <= raw_len/3 -> 3*m <= raw_len
    runs = st_malloc((cap ? cap : 1) * sizeof(int64_t));
    m = decode_runs(def, def_len, raw_len, runs, cap);
    oneFileClose(of);
    if (m == 0) { free(runs); return NULL; }

    // Clip each overlapping run to [start,end), map to a column interval.
    TuiInterval *iv = st_malloc((size_t)m * sizeof(TuiInterval));
    int64_t k = 0;
    for (int64_t r = 0; r < m; r++) {
        int64_t t = runs[3*r+0], g = runs[3*r+1], lenc = runs[3*r+2];
        int64_t len = lenc >> 1, rev = lenc & 1, te = t + len;
        int64_t a = start > t ? start : t;
        int64_t b = end < te ? end : te;
        if (a >= b) continue;
        if (!rev) {                             // col = g + (p - t), increasing
            iv[k].start = g + (a - t);
            iv[k].end   = g + (b - t);
        } else {                                // col = g + (t+len-1 - p), decreasing
            iv[k].start = g + (t + len - b);
            iv[k].end   = g + (t + len - a);
        }
        k++;
    }
    free(runs);
    if (k == 0) { free(iv); return NULL; }

    qsort(iv, k, sizeof(TuiInterval), tui_iv_cmp);
    int64_t w = 0;                              // merge adjacent / overlapping
    for (int64_t i = 1; i < k; i++) {
        if (iv[i].start <= iv[w].end) {
            if (iv[i].end > iv[w].end) iv[w].end = iv[i].end;
        } else {
            iv[++w] = iv[i];
        }
    }
    *n_out = w + 1;
    return iv;
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

void tui_extract_iterator_destruct(TuiExtractIt *it) {
    if (it->to_free != NULL) alignment_destruct(it->to_free, 1);
    if (it->phys != NULL) alignment_destruct(it->phys, 1);
    free(it->runs);
    free(it->iv);
    free(it);
}
