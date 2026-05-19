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
static const char *TUI_SCHEMA =
    "P 3 tui\n"
    "D t 1 3 INT\n"                          // total columns T (global)
    "D d 4 6 STRING 3 INT 3 INT 3 INT\n"     // dir: seqName, S-ordinal, seqLen, isRef
    "O S 2 6 STRING 3 INT\n"                 // sequence object: seqName, seqLen
    "D R 2 3 INT 6 STRING\n";                // runs: inflatedLen, deflate(blob)

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
    fprintf(fh, "%s\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%c\n",
            seq_name, t_start, g_start, length, strand ? '+' : '-');
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
} Phase1;

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

    Alignment *aln, *p_aln = NULL;
    if (input_format == 1) {
        while ((aln = maf_read_block(li)) != NULL) {
            tui_phase1_block(&p1, aln);
            alignment_destruct(aln, 1);
        }
    } else {
        while ((aln = taf_read_block(p_aln, rle, li)) != NULL) {
            tui_phase1_block(&p1, aln);
            if (p_aln != NULL) alignment_destruct(p_aln, 1);
            p_aln = aln;
        }
        if (p_aln != NULL) alignment_destruct(p_aln, 1);
    }
    // close all spill FILE*s (paths kept for phase 2)
    stHashIterator *hit = stHash_getIterator(p1.spill_fh);
    char *gk;
    while ((gk = stHash_getNext(hit)) != NULL) fclose(stHash_search(p1.spill_fh, gk));
    stHash_destructIterator(hit);

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

    // d: directory (front of file), S-ordinal == position in global order
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
    int64_t i = 0;
    while (i < n_seqs) {
        char *gname = seqks[i].genome;  // groups spill files; not emitted

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
        for (int64_t k = 0; k < nr; k++) free(runs[k].seq);
        free(runs);
        // gname is seqks[i].genome -- owned by seqks, freed in tui_cleanup
    }

    oneFileClose(of);
    oneSchemaDestroy(schema);
    tui_cleanup(&p1, seqks, n_seqs, eff_tmp, tree_map);
    return 0;
}
