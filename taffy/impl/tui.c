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

// ONEcode schema. Exactly one list field per line (STRING / INT_LIST); other
// fields are scalar INT.  The directory keys on the full genome.sequence name
// (genome = prefix before the first '.', derived identically at build/load).
static const char *TUI_SCHEMA =
    "P 3 tui\n"
    "D t 1 3 INT\n"                          // total columns T (global)
    "D d 4 6 STRING 3 INT 3 INT 3 INT\n"     // dir: seqName, S-ordinal, seqLen, isRef
    "O S 2 6 STRING 3 INT\n"                 // sequence object: seqName, seqLen
    "D R 1 8 INT_LIST\n";                    // runs: [t, g, (len<<1|revStrand)]*

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

// genome = prefix of seq name before the first '.'  (dots-in-names deferred)
static char *genome_of(const char *seq_name) {
    const char *dot = strchr(seq_name, '.');
    int64_t n = dot ? (dot - seq_name) : (int64_t)strlen(seq_name);
    char *g = st_malloc(n + 1);
    memcpy(g, seq_name, n);
    g[n] = '\0';
    return g;
}

// Order full "genome.sequence" names genome-major (split at first '.') so each
// genome's sequences form one contiguous block.  THE single source of truth for
// sequence ordering: both seqs[] (qsort) and runs[] (run_cmp) must use this, or
// the Phase-2 forward cursor over runs[] can skip/misattribute a seq's runs.
static int seqname_cmp(const char *x, const char *y) {
    const char *xd = strchr(x, '.'), *yd = strchr(y, '.');
    size_t xg = xd ? (size_t)(xd - x) : strlen(x);
    size_t yg = yd ? (size_t)(yd - y) : strlen(y);
    int c = strncmp(x, y, xg < yg ? xg : yg);
    if (c != 0) return c;
    if (xg != yg) return xg < yg ? -1 : 1;   // shorter genome prefix first
    return strcmp(x, y);                      // same genome -> by full seq name
}

static int seqkey_cmp(const void *a, const void *b) {
    return seqname_cmp(*(const char **)a, *(const char **)b);
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
} Phase1;

static FILE *spill_for(Phase1 *p1, const char *genome) {
    FILE *fh = stHash_search(p1->spill_fh, (void *)genome);
    if (fh == NULL) {
        char *path = stString_print("%s/%s.tuiSpill.%d", p1->tmp_dir,
                                    p1->out_base, p1->next_spill_id++);
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
    char *ref_g = genome_of(aln->row->sequence_name);
    if (stSet_search(p1->ref_genomes, ref_g) == NULL) stSet_insert(p1->ref_genomes, ref_g);
    else free(ref_g);
    for (Alignment_Row *row = aln->row; row != NULL; row = row->n_row) {
        char *gname = genome_of(row->sequence_name);
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

static int run_cmp(const void *a, const void *b) {
    const Run *x = a, *y = b;
    int s = seqname_cmp(x->seq, y->seq);   // same ordering as seqs[] (seqkey_cmp)
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
        }
    }
    fclose(fh);
    *n_out = n;
    return runs;
}

// Remove spill files and free phase-1 state.  Shared by the success and the
// error-return paths so a failed run never leaks spill files on disk.
static void tui_cleanup(Phase1 *p1, char **seqs, char *eff_tmp) {
    stHashIterator *hit = stHash_getIterator(p1->spill_path);
    char *gk;
    while ((gk = stHash_getNext(hit)) != NULL) {
        remove((char *)stHash_search(p1->spill_path, gk));
    }
    stHash_destructIterator(hit);
    free(seqs);
    stHash_destruct(p1->spill_fh);
    stHash_destruct(p1->spill_path);
    stHash_destruct(p1->seq_len);
    stList_destruct(p1->seq_keys);
    stSet_destruct(p1->ref_genomes);
    free(eff_tmp);
}

int tui_create(LI *li, const char *out_path, const char *tmp_dir) {
    int input_format = check_input_format(LI_peek_at_next_line(li));  // 0=TAF 1=MAF
    if (input_format != 0 && input_format != 1) {
        fprintf(stderr, "tui: input must be MAF or TAF\n");
        return 1;
    }
    bool rle = 0;
    Tag *tag = (input_format == 0) ? taf_read_header_2(li, &rle) : maf_read_header(li);
    tag_destruct(tag);

    // spill dir: explicit --tmpDir, else the output file's own directory
    char *eff_tmp = (tmp_dir != NULL && tmp_dir[0] != '\0')
                        ? stString_copy(tmp_dir) : dir_of(out_path);

    Phase1 p1;
    p1.tmp_dir = eff_tmp;
    p1.out_base = base_of(out_path);
    p1.next_spill_id = 0;
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

    // deterministic global order: genome-major, then sequence
    int64_t n_seqs = stList_length(p1.seq_keys);
    char **seqs = st_malloc((n_seqs ? n_seqs : 1) * sizeof(char *));
    for (int64_t i = 0; i < n_seqs; i++) seqs[i] = stList_get(p1.seq_keys, i);
    qsort(seqs, n_seqs, sizeof(char *), seqkey_cmp);

    OneSchema *schema = oneSchemaCreateFromText(TUI_SCHEMA);
    OneFile *of = oneFileOpenWriteNew(out_path, schema, "tui", true, 1);
    if (of == NULL) {
        fprintf(stderr, "tui: cannot write %s\n", out_path);
        oneSchemaDestroy(schema);
        tui_cleanup(&p1, seqs, eff_tmp);
        return 1;
    }
    oneAddProvenance(of, "taffy", "tui", "universal column index", 0);

    // t: total columns
    oneInt(of, 0) = p1.T;
    oneWriteLine(of, 't', 0, NULL);

    // d: directory (front of file), S-ordinal == position in global order
    for (int64_t i = 0; i < n_seqs; i++) {
        char *sk = seqs[i];
        int64_t slen = *(int64_t *)stHash_search(p1.seq_len, sk);
        char *gg = genome_of(sk);
        int is_ref = stSet_search(p1.ref_genomes, gg) != NULL ? 1 : 0;
        free(gg);
        oneInt(of, 1) = i;
        oneInt(of, 2) = slen;
        oneInt(of, 3) = is_ref;
        oneWriteLine(of, 'd', strlen(sk), (void *)sk);
    }

    // per genome (contiguous in seqs[] thanks to seqkey_cmp): g + S/R objects
    int64_t i = 0;
    while (i < n_seqs) {
        char *gname = genome_of(seqs[i]);   // groups spill files; not emitted

        int64_t nr = 0;
        Run *runs = NULL;
        char *path = stHash_search(p1.spill_path, gname);
        if (path != NULL) {
            runs = load_genome_runs(path, &nr);
            if (runs != NULL) qsort(runs, nr, sizeof(Run), run_cmp);
        }
        // walk this genome's seqs (contiguous block in seqs[]); runs[] are
        // seq-sorted consistently, so advance a single forward cursor.
        int64_t rc = 0;
        while (i < n_seqs) {
            char *gg = genome_of(seqs[i]);
            int same = strcmp(gg, gname) == 0;
            free(gg);
            if (!same) break;
            char *sk = seqs[i];
            int64_t slen = *(int64_t *)stHash_search(p1.seq_len, sk);
            oneInt(of, 1) = slen;
            oneWriteLine(of, 'S', strlen(sk), (void *)sk);
            // R: runs of this seq (contiguous in sorted runs[]); same ordering
            // (seqname_cmp) as seqs[] so the forward cursor can't skip/misattribute
            while (rc < nr && seqname_cmp(runs[rc].seq, sk) < 0) rc++;  // forward only
            int64_t a = rc;
            int64_t bnd = a;
            while (bnd < nr && seqname_cmp(runs[bnd].seq, sk) == 0) bnd++;
            int64_t cnt = bnd - a;
            int64_t *buf = st_malloc((cnt ? 3 * cnt : 1) * sizeof(int64_t));
            for (int64_t k = 0; k < cnt; k++) {
                Run *r = &runs[a + k];
                buf[3*k+0] = r->t;
                buf[3*k+1] = r->g;
                buf[3*k+2] = (r->len << 1) | (r->strand == '-' ? 1 : 0);
            }
            oneWriteLine(of, 'R', 3 * cnt, buf);
            free(buf);
            rc = bnd;   // advance cursor past this seq's runs
            i++;
        }
        for (int64_t k = 0; k < nr; k++) free(runs[k].seq);
        free(runs);
        free(gname);
    }

    oneFileClose(of);
    oneSchemaDestroy(schema);
    tui_cleanup(&p1, seqs, eff_tmp);
    return 0;
}
