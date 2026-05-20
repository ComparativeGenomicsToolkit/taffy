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
#include <time.h>
#include <unistd.h>
#include <zlib.h>

// ONEcode schema. Exactly one list field per line (STRING / INT_LIST); other
// fields are scalar INT.  The directory keys on the full genome.sequence name
// (genome derived identically at build/load via genome_of()).
//
// R holds one sequence's merged runs as a per-sequence structure-of-arrays
// delta+varint blob, zlib-deflated (see encode_runs).  The scalar INT is the
// inflated blob length (sizes the inflate buffer); the STRING list is the
// deflated bytes.  Absolute (t,g,len) triples defeat ONElib's Huffman; the
// SoA gap|gsk|lenc form is ~5x smaller than absolute end to end.
//
// A is the explicit Index-A reference track (the materialized canonical BED):
// the column-ordered row-0 (block-reference) segments tiling [0,T).  Stored
// like R -- (inflatedLen INT, deflated STRING) -- as a SoA delta+varint blob
// (see encode_refseg); colStart is dropped (= prefix sum of len, since the
// segments tile [0,T)).  Step-2 query uses it to turn universal columns into
// row-0 coords for the existing .tai.  Written once, right after `t`, before
// `d`.  Strand is always '+' (asserted at build; tai rejects a '-' row-0).
static const char *TUI_SCHEMA =
    "P 3 tui\n"
    "D t 1 3 INT\n"                          // total columns T (global)
    "D A 3 3 INT 3 INT 6 STRING\n"           // Index-A: inflatedLen, nSeg, deflate(SoA blob)
    "D d 4 6 STRING 3 INT 3 INT 3 INT\n"     // dir: seqName, S-ordinal, seqLen, isRef
    "O S 2 6 STRING 3 INT\n"                 // sequence object: seqName, seqLen
    "D R 2 3 INT 6 STRING\n"                 // runs: inflatedLen, deflate(blob)
    "D X 3 3 INT 3 INT 6 STRING\n";          // univ-col index: inflatedLen, nRec, deflate(SoA)

char *tui_path(const char *maf_path) {
    char *p = st_malloc(strlen(maf_path) + 5);
    sprintf(p, "%s.tui", maf_path);
    return p;
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
    stSet  *ref_genomes;  // genomes that appear as row-0 (ancestral)
    int64_t T;            // total columns
    const char *tmp_dir;  // where spill files go
    const char *out_base; // basename of out path (for unique spill names)
    int next_spill_id;
    stHash *gmap;         // optional name map for >1-dot genome resolution
    // Index-A reference track: one open row-0 segment, colinear-merged on the
    // fly, flushed (column-ordered) to a single spill file.
    FILE   *ref_fh;       // ref-track spill (lazy open); NULL until first flush
    char   *ref_path;     // we own/free this
    char   *ro_seq;       // open segment's row-0 sequence name (strdup)
    int64_t ro_col0, ro_row0, ro_len;
    int     ro_open;
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

// Flush the open row-0 segment (if any) to the ref-track spill.
static void ref_flush(Phase1 *p1) {
    if (!p1->ro_open) return;
    if (p1->ref_fh == NULL) {
        p1->ref_path = stString_print("%s/%s.tuiRef.%ld", p1->tmp_dir,
                                      p1->out_base, (long)getpid());
        p1->ref_fh = fopen(p1->ref_path, "w");
        if (p1->ref_fh == NULL) {
            fprintf(stderr, "tui: cannot open ref spill %s\n", p1->ref_path);
            exit(1);
        }
    }
    if (fprintf(p1->ref_fh, "%" PRIi64 "\t%s\t%" PRIi64 "\t%" PRIi64 "\n",
                p1->ro_col0, p1->ro_seq, p1->ro_row0, p1->ro_len) < 0) {
        st_errAbort("tui: failed writing ref spill (disk full / write error)");
    }
    free(p1->ro_seq);
    p1->ro_seq = NULL;
    p1->ro_open = 0;
}

// Record this block's row-0 segment into the on-the-fly-merged ref track.
// row0 = aln->row (the block reference), col0 = first universal column of the
// block, cn = column_number (== row-0 base count, since row-0 is gap-free).
static void ref_track_block(Phase1 *p1, Alignment_Row *row0,
                            int64_t col0, int64_t cn) {
    if (!row0->strand) {
        st_errAbort("tui: row-0 sequence '%s' is on '-' strand; a universal "
                    "MAF row-0 must be '+' (the .tai requires it)",
                    row0->sequence_name);
    }
    if (row0->length != cn) {
        st_errAbort("tui: row-0 '%s' is not gap-free (length %" PRIi64
                    " != column_number %" PRIi64 "); --universal / "
                    "maxRefGap==0 expected", row0->sequence_name,
                    row0->length, cn);
    }
    if (p1->ro_open && p1->ro_row0 + p1->ro_len == row0->start &&
        strcmp(p1->ro_seq, row0->sequence_name) == 0) {
        assert(p1->ro_col0 + p1->ro_len == col0);  // columns globally sequential
        p1->ro_len += cn;
        return;
    }
    ref_flush(p1);
    p1->ro_seq  = stString_copy(row0->sequence_name);
    p1->ro_col0 = col0;
    p1->ro_row0 = row0->start;
    p1->ro_len  = cn;
    p1->ro_open = 1;
}

static FILE *spill_for(Phase1 *p1, const char *genome) {
    FILE *fh = stHash_search(p1->spill_fh, (void *)genome);
    if (fh == NULL) {
        char *path = stString_print("%s/%s.tuiSpill.%ld.%d", p1->tmp_dir,
                                    p1->out_base, (long)getpid(),
                                    p1->next_spill_id++);
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
    char *ref_g = genome_of(aln->row->sequence_name, p1->gmap);
    if (stSet_search(p1->ref_genomes, ref_g) == NULL) stSet_insert(p1->ref_genomes, ref_g);
    else free(ref_g);
    ref_track_block(p1, aln->row, p1->T, cn);   // Index A: this block's row-0
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

// Index-A reference-track codec.  Same shape as encode_runs (SoA delta+varint,
// one zlib deflate) but tuned for the row-0 segments.  `colStart` is NOT
// stored: the segments tile [0,T) exactly (builder asserts contiguity), so
// colStart is the prefix sum of `len` -- reconstructed at decode.  Streams:
//   O = zigzag(sOrd - prevOrd)     (often 0: same row-0 seq across a jump)
//   P = zigzag(row0Start - prevRow0)
//   L = uvarint(len)
// header = uvarint(nSeg), uvarint(|O|), uvarint(|P|); raw = hdr+O+P+L.
static uint8_t *encode_refseg(const int64_t *sOrd, const int64_t *row0,
                              const int64_t *len, int64_t nSeg,
                              int64_t *raw_len, int64_t *def_len) {
    uint8_t *O = st_malloc((size_t)(nSeg * 10 + 1));
    uint8_t *P = st_malloc((size_t)(nSeg * 10 + 1));
    uint8_t *L = st_malloc((size_t)(nSeg * 10 + 1));
    size_t on = 0, pn = 0, ln = 0;
    int64_t po = 0, pr = 0;
    for (int64_t k = 0; k < nSeg; k++) {
        on += put_uvarint(O + on, zigzag(sOrd[k] - po));
        pn += put_uvarint(P + pn, zigzag(row0[k] - pr));
        ln += put_uvarint(L + ln, (uint64_t)len[k]);
        po = sOrd[k]; pr = row0[k];
    }
    uint8_t hdr[30];
    size_t hn = 0;
    hn += put_uvarint(hdr + hn, (uint64_t)nSeg);
    hn += put_uvarint(hdr + hn, (uint64_t)on);
    hn += put_uvarint(hdr + hn, (uint64_t)pn);
    size_t rn = hn + on + pn + ln;
    uint8_t *raw = st_malloc(rn ? rn : 1);
    memcpy(raw, hdr, hn);
    memcpy(raw + hn, O, on);
    memcpy(raw + hn + on, P, pn);
    memcpy(raw + hn + on + pn, L, ln);
    free(O); free(P); free(L);
    uLongf cap = compressBound((uLong)rn);
    uint8_t *def = st_malloc(cap ? cap : 1);
    uLongf dl = cap;
    if (compress2(def, &dl, raw, (uLong)rn, 9) != Z_OK) {
        st_errAbort("tui: zlib compress2 failed (A)");
    }
    free(raw);
    *raw_len = (int64_t)rn;
    *def_len = (int64_t)dl;
    return def;
}

// Inverse of encode_refseg.  Fills colStart/sOrd/row0/len[0..nSeg) (colStart
// reconstructed as the running prefix sum of len, from 0).  Returns nSeg.
static int64_t decode_refseg(const uint8_t *def, int64_t def_len,
                             int64_t raw_len, int64_t *colStart,
                             int64_t *sOrd, int64_t *row0, int64_t *len) {
    uint8_t *raw = st_malloc((size_t)(raw_len ? raw_len : 1));
    uLongf rl = (uLongf)raw_len;
    if (uncompress(raw, &rl, def, (uLong)def_len) != Z_OK || (int64_t)rl != raw_len) {
        st_errAbort("tui: zlib uncompress failed (A)");
    }
    const uint8_t *h = raw;
    int64_t nSeg = (int64_t)get_uvarint(&h);
    int64_t on = (int64_t)get_uvarint(&h);
    int64_t pn = (int64_t)get_uvarint(&h);
    const uint8_t *op = h, *pp = h + on, *lp = h + on + pn;
    int64_t po = 0, pr = 0, col = 0;
    for (int64_t k = 0; k < nSeg; k++) {
        int64_t o = po + unzigzag(get_uvarint(&op));
        int64_t r = pr + unzigzag(get_uvarint(&pp));
        int64_t l = (int64_t)get_uvarint(&lp);
        colStart[k] = col; sOrd[k] = o; row0[k] = r; len[k] = l;
        col += l; po = o; pr = r;
    }
    free(raw);
    return nSeg;
}

// Universal-column index codec.  n anchors (col, file_pos), both STRICTLY
// increasing in anchor order, so plain non-negative uvarint deltas (no
// zigzag).  Two SoA streams C|F; header = uvarint(n), uvarint(|C bytes|);
// one zlib deflate.  Mirrors encode_refseg.
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
    if (p1->ref_fh != NULL) fclose(p1->ref_fh);          // error-path: still open
    if (p1->ref_path != NULL) { remove(p1->ref_path); free(p1->ref_path); }
    if (p1->idx_fh != NULL) fclose(p1->idx_fh);
    if (p1->idx_path != NULL) { remove(p1->idx_path); free(p1->idx_path); }
    if (p1->ro_open) free(p1->ro_seq);                   // error-path: open segment
    for (int64_t i = 0; i < n_seqs; i++) free(seqks[i].genome);  // .seq owned by seq_len
    free(seqks);
    stHash_destruct(p1->spill_fh);
    stHash_destruct(p1->spill_path);
    stHash_destruct(p1->seq_len);
    stList_destruct(p1->seq_keys);
    stSet_destruct(p1->ref_genomes);
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
    p1.ref_genomes = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    p1.T = 0;
    p1.ref_fh = NULL; p1.ref_path = NULL; p1.ro_seq = NULL; p1.ro_open = 0;
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
    ref_flush(&p1);                       // flush the last open row-0 segment
    // fclose surfaces buffered short writes (e.g. disk full) -> fail loudly
    // here rather than emit a silently-truncated spill that Phase 2 would
    // later misread as "not contiguous" / "corrupt spill line".
    if (p1.ref_fh != NULL) {
        if (fclose(p1.ref_fh) != 0) st_errAbort("tui: ref spill close failed "
            "(disk full / write error) -- %s", p1.ref_path);
        p1.ref_fh = NULL;
    }
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
        return 1;
    }
    oneAddProvenance(of, "taffy", "tui", "universal column index", 0);

    // t: total columns
    oneInt(of, 0) = p1.T;
    oneWriteLine(of, 't', 0, NULL);

    // A: Index-A reference track.  Resolve each spilled row-0 segment's
    // sequence name to its S-ordinal (== position in the sorted global order),
    // then SoA delta+varint+deflate (colStart dropped = prefix sum of len).
    time_t t_idxA_start = time(NULL);
    st_logInfo("tui: Index A: encoding row-0 reference track\n");
    stHash *name2ord = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
                                         NULL, NULL);
    for (int64_t i = 0; i < n_seqs; i++) {
        stHash_insert(name2ord, seqks[i].seq, (void *)(intptr_t)(i + 1));
    }
    int64_t *segOrd = NULL, *segRow0 = NULL, *segLen = NULL, *segCol0 = NULL;
    int64_t nSeg = 0, segCap = 0;
    if (p1.ref_path != NULL) {
        FILE *rf = fopen(p1.ref_path, "r");
        if (rf == NULL) st_errAbort("tui: cannot reopen ref spill %s", p1.ref_path);
        char line[16384], sbuf[8192];
        int64_t c0, r0, ln, runcol = 0;
        while (fgets(line, sizeof(line), rf) != NULL) {
            if (sscanf(line, "%" SCNi64 "\t%8191s\t%" SCNi64 "\t%" SCNi64,
                       &c0, sbuf, &r0, &ln) != 4) {
                st_errAbort("tui: corrupt ref spill line: %s", line);
            }
            void *ov = stHash_search(name2ord, sbuf);
            if (ov == NULL) st_errAbort("tui: ref seq '%s' not in directory", sbuf);
            // segments must tile [0,T): colStart == running prefix sum of len
            if (c0 != runcol) {
                st_errAbort("tui: ref track not contiguous (seg colStart %"
                            PRIi64 " != expected %" PRIi64 ")", c0, runcol);
            }
            runcol += ln;
            if (nSeg == segCap) {
                segCap = segCap ? segCap * 2 : 4096;
                segOrd  = st_realloc(segOrd,  segCap * sizeof(int64_t));
                segRow0 = st_realloc(segRow0, segCap * sizeof(int64_t));
                segLen  = st_realloc(segLen,  segCap * sizeof(int64_t));
                segCol0 = st_realloc(segCol0, segCap * sizeof(int64_t));
            }
            segCol0[nSeg] = c0;
            segOrd[nSeg]  = (int64_t)(intptr_t)ov - 1;
            segRow0[nSeg] = r0;
            segLen[nSeg]  = ln;
            nSeg++;
        }
        fclose(rf);
    }
    {
        int64_t a_raw = 0, a_def = 0;
        uint8_t *adef = encode_refseg(segOrd, segRow0, segLen, nSeg,
                                      &a_raw, &a_def);
        // self-check: decode must reproduce the segments (asserts on in taffy)
        int64_t *cC = st_malloc((nSeg ? nSeg : 1) * sizeof(int64_t));
        int64_t *cO = st_malloc((nSeg ? nSeg : 1) * sizeof(int64_t));
        int64_t *cR = st_malloc((nSeg ? nSeg : 1) * sizeof(int64_t));
        int64_t *cL = st_malloc((nSeg ? nSeg : 1) * sizeof(int64_t));
        int64_t dn = decode_refseg(adef, a_def, a_raw, cC, cO, cR, cL);
        assert(dn == nSeg);
        for (int64_t k = 0; k < nSeg; k++) {
            assert(cC[k] == segCol0[k] && cO[k] == segOrd[k] &&
                   cR[k] == segRow0[k] && cL[k] == segLen[k]);
        }
        free(cC); free(cO); free(cR); free(cL);
        oneInt(of, 0) = a_raw;
        oneInt(of, 1) = nSeg;
        oneWriteLine(of, 'A', a_def, adef);
        st_logInfo("tui: Index A done in %" PRIi64 " seconds "
                   "(%" PRIi64 " segments, %" PRIi64 " raw bytes, "
                   "%" PRIi64 " deflated bytes)\n",
                   (int64_t)(time(NULL) - t_idxA_start), nSeg, a_raw, a_def);
        free(adef);
    }
    free(segOrd); free(segRow0); free(segLen); free(segCol0);
    stHash_destruct(name2ord);

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

    // d: directory (front of file), S-ordinal == position in global order
    st_logInfo("tui: writing directory (%" PRIi64 " sequences across %" PRIi64
               " genomes)\n", n_seqs, n_genomes);
    for (int64_t i = 0; i < n_seqs; i++) {
        char *sk = seqks[i].seq;
        int64_t slen = *(int64_t *)stHash_search(p1.seq_len, sk);
        int is_ref = stSet_search(p1.ref_genomes, seqks[i].genome) != NULL ? 1 : 0;
        oneInt(of, 1) = i;
        oneInt(of, 2) = slen;
        oneInt(of, 3) = is_ref;
        oneWriteLine(of, 'd', strlen(sk), (void *)sk);
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
    tui_cleanup(&p1, seqks, n_seqs, eff_tmp, tree_map);
    return 0;
}

/////////////////////////////////////////////////////////////////////////////
// Reader / query side (Index B).  Baseline: re-scan to find a sequence's R
// blob, inflate+decode it, clip to the query interval.  Correctness first.
/////////////////////////////////////////////////////////////////////////////

struct _Tui {
    char   *path;       // .tui path (re-opened per query for now)
    int64_t T;          // total universal columns
    stHash *seq_len;    // full "genome.sequence" -> int64_t* sequence length
    // Index-A reference track (column-ordered, tiles [0,T)).
    int64_t *segCol0, *segOrd, *segRow0, *segLen;
    int64_t  nSeg;
    char   **ord2name;  // S-ordinal -> sequence name (we own the strings)
    int64_t  nOrd;
    // Universal-column -> file_pos index (X track); both strictly increasing.
    int64_t *idxCol, *idxFpos;
    int64_t  idxN;
};

Tui *tui_load(const char *tui_path) {
    OneFile *of = oneFileOpenRead(tui_path, NULL, "tui", 1);
    if (of == NULL) return NULL;
    Tui *tui = st_calloc(1, sizeof(Tui));
    tui->path = stString_copy(tui_path);
    tui->seq_len = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    int64_t ord_cap = 0;
    char c;
    while ((c = oneReadLine(of)) != 0) {
        if (c == 't') {
            tui->T = oneInt(of, 0);
        } else if (c == 'A') {                 // Index-A: (rawLen, nSeg, deflate blob)
            int64_t a_raw = oneInt(of, 0);
            int64_t nSeg  = oneInt(of, 1);
            int64_t a_def = oneLen(of);
            const uint8_t *adef = (const uint8_t *)oneString(of);
            int64_t cap = nSeg ? nSeg : 1;
            tui->segCol0 = st_malloc(cap * sizeof(int64_t));
            tui->segOrd  = st_malloc(cap * sizeof(int64_t));
            tui->segRow0 = st_malloc(cap * sizeof(int64_t));
            tui->segLen  = st_malloc(cap * sizeof(int64_t));
            tui->nSeg = decode_refseg(adef, a_def, a_raw, tui->segCol0,
                                      tui->segOrd, tui->segRow0, tui->segLen);
            assert(tui->nSeg == nSeg);
        } else if (c == 'X') {                 // univ-col index: (rawLen,n,blob)
            int64_t x_raw = oneInt(of, 0);
            int64_t xn    = oneInt(of, 1);
            int64_t x_def = oneLen(of);
            const uint8_t *xdef = (const uint8_t *)oneString(of);
            int64_t cap = xn ? xn : 1;
            tui->idxCol  = st_malloc(cap * sizeof(int64_t));
            tui->idxFpos = st_malloc(cap * sizeof(int64_t));
            tui->idxN = decode_idx(xdef, x_def, x_raw, tui->idxCol, tui->idxFpos);
            assert(tui->idxN == xn);
        } else if (c == 'd') {                 // dir: name(STRING), ord, seqLen, isRef
            int64_t n = oneLen(of);
            char *nm = st_malloc(n + 1);
            memcpy(nm, oneString(of), n);
            nm[n] = '\0';
            int64_t ord = oneInt(of, 1);
            if (ord >= ord_cap) {
                int64_t nc = ord_cap ? ord_cap * 2 : 1024;
                while (nc <= ord) nc *= 2;
                tui->ord2name = st_realloc(tui->ord2name, nc * sizeof(char *));
                for (int64_t z = ord_cap; z < nc; z++) tui->ord2name[z] = NULL;
                ord_cap = nc;
            }
            tui->ord2name[ord] = stString_copy(nm);
            if (ord + 1 > tui->nOrd) tui->nOrd = ord + 1;
            if (stHash_search(tui->seq_len, nm) == NULL) {
                int64_t *v = st_malloc(sizeof(int64_t));
                *v = oneInt(of, 2);
                stHash_insert(tui->seq_len, nm, v);
            } else {
                free(nm);
            }
        } else if (c == 'S') {
            break;                              // directory precedes the objects
        }
    }
    oneFileClose(of);
    return tui;
}

void tui_destruct(Tui *tui) {
    if (tui == NULL) return;
    stHash_destruct(tui->seq_len);
    for (int64_t i = 0; i < tui->nOrd; i++) free(tui->ord2name[i]);
    free(tui->ord2name);
    free(tui->segCol0); free(tui->segOrd); free(tui->segRow0); free(tui->segLen);
    free(tui->idxCol); free(tui->idxFpos);
    free(tui->path);
    free(tui);
}

int64_t tui_total_columns(const Tui *tui) { return tui->T; }

int tui_has_sequence(const Tui *tui, const char *seq_name) {
    return stHash_search(tui->seq_len, (void *)seq_name) != NULL;
}

static int tui_iv_cmp(const void *a, const void *b) {
    const TuiInterval *x = a, *y = b;
    return (x->start < y->start) ? -1 : (x->start > y->start) ? 1 : 0;
}

TuiInterval *tui_query(Tui *tui, const char *seq_name,
                       int64_t start, int64_t end, int64_t *n_out) {
    *n_out = 0;
    int64_t *expect = stHash_search(tui->seq_len, (void *)seq_name);
    if (expect == NULL || start >= end) return NULL;

    OneFile *of = oneFileOpenRead(tui->path, NULL, "tui", 1);
    if (of == NULL) return NULL;

    int64_t *runs = NULL, m = 0;
    char *cur = NULL;
    int found = 0;
    char c;
    while (!found && (c = oneReadLine(of)) != 0) {
        if (c == 'S') {
            int64_t n = oneLen(of);
            free(cur);
            cur = st_malloc(n + 1);
            memcpy(cur, oneString(of), n);
            cur[n] = '\0';
        } else if (c == 'R' && cur != NULL && strcmp(cur, seq_name) == 0) {
            int64_t raw_len = oneInt(of, 0);
            int64_t def_len = oneLen(of);
            const uint8_t *def = (const uint8_t *)oneString(of);
            int64_t cap = raw_len + 3;          // m <= raw_len/3 -> 3*m <= raw_len
            runs = st_malloc((cap ? cap : 1) * sizeof(int64_t));
            m = decode_runs(def, def_len, raw_len, runs, cap);
            found = 1;
        }
    }
    free(cur);
    oneFileClose(of);
    if (!found || m == 0) { free(runs); return NULL; }

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

// Index A: map universal-column intervals -> row-0 (ancestor) coordinate
// pieces, using the explicit reference track.  segCol0 is sorted ascending
// and the segments tile [0,T), so a column has exactly one segment (binary
// search), and an interval spans a contiguous run of segments.  One input
// interval can yield several pieces on DIFFERENT row-0 sequences/ancestors
// (e.g. hg38.chr1:1-100 -> AncR..  + AncC..).  Returned array is in column
// order; .seq points into tui->ord2name (do NOT free per piece -- free the
// array only).  Strand is always '+' so the map is the affine
// row0 = segRow0 + (col - segCol0).
static int64_t tui_seg_find(const Tui *tui, int64_t col) {
    int64_t lo = 0, hi = tui->nSeg - 1, ans = -1;
    while (lo <= hi) {
        int64_t mid = (lo + hi) / 2;
        if (tui->segCol0[mid] <= col) { ans = mid; lo = mid + 1; }
        else hi = mid - 1;
    }
    return ans;  // greatest seg with segCol0 <= col, or -1
}

TuiRef *tui_col_range_to_ref(const Tui *tui, const TuiInterval *uiv,
                             int64_t n_uiv, int64_t *n_out) {
    *n_out = 0;
    if (tui->nSeg == 0 || n_uiv == 0) return NULL;
    int64_t cap = 16, n = 0;
    TuiRef *out = st_malloc(cap * sizeof(TuiRef));
    for (int64_t q = 0; q < n_uiv; q++) {
        int64_t a = uiv[q].start, b = uiv[q].end;     // [a, b) universal columns
        int64_t s = tui_seg_find(tui, a);
        if (s < 0) s = 0;                              // before first seg (shouldn't happen)
        for (; s < tui->nSeg && tui->segCol0[s] < b; s++) {
            int64_t segA = tui->segCol0[s];
            int64_t segB = segA + tui->segLen[s];
            int64_t ca = a > segA ? a : segA;
            int64_t cb = b < segB ? b : segB;
            if (ca >= cb) continue;                    // no overlap with this seg
            if (n == cap) { cap *= 2; out = st_realloc(out, cap * sizeof(TuiRef)); }
            out[n].seq   = tui->ord2name[tui->segOrd[s]];
            out[n].start = tui->segRow0[s] + (ca - segA);
            out[n].len   = cb - ca;
            n++;
        }
    }
    *n_out = n;
    if (n == 0) { free(out); return NULL; }
    return out;
}

/////////////////////////////////////////////////////////////////////////////
// Universal-column block extractor.  Single forward scan in column order from
// the X-index anchor; reuses tai's seek/resync/read primitives.  No per-
// contig key => the C1 (non-monotone row-0 ancestor) failure is impossible.
/////////////////////////////////////////////////////////////////////////////

struct _TuiExtractIt {
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

// Forward-scan (ctx = taf context, not freed here) skipping non-overlapping
// physical blocks, until one overlaps iv -> it->phys/phys_col/runs(run_i=0),
// or past the last interval / EOF -> it->phys = NULL.
static void tui_extract_advance(TuiExtractIt *it, LI *li, Alignment *ctx) {
    Alignment *prev = ctx;
    it->phys = NULL;
    while (1) {
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
        if (prev != ctx) alignment_destruct(prev, 1);   // intermediate skipped
        prev = b;
    }
    if (prev != ctx && prev != it->phys) alignment_destruct(prev, 1);
}

TuiExtractIt *tui_extract_iterator(Tui *tui, LI *li, int is_maf, bool rle,
                                   const TuiInterval *iv, int64_t n_iv) {
    TuiExtractIt *it = st_calloc(1, sizeof(TuiExtractIt));
    it->is_maf = is_maf;
    it->rle = rle;
    it->readblk = is_maf ? tai_maf_read_block : taf_read_block;
    if (n_iv <= 0 || tui->idxN == 0) return it;        // empty iterator
    it->iv = st_malloc((size_t)n_iv * sizeof(TuiInterval));
    memcpy(it->iv, iv, (size_t)n_iv * sizeof(TuiInterval));
    it->n_iv = n_iv;
    it->c_hi = iv[n_iv - 1].end;
    int64_t c_lo = iv[0].start;
    // anchor = greatest idxCol <= c_lo (idxCol strictly increasing, [0]=0)
    int64_t lo = 0, hi = tui->idxN - 1, a = 0;
    while (lo <= hi) {
        int64_t mid = (lo + hi) / 2;
        if (tui->idxCol[mid] <= c_lo) { a = mid; lo = mid + 1; } else hi = mid - 1;
    }
    LI_seek(li, tui->idxFpos[a]);
    LI_get_next_line(li);                              // mirror tai_iterator
    if (!is_maf) tai_resync_taf_line(LI_peek_at_next_line(li));
    it->scan_col = tui->idxCol[a];
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
