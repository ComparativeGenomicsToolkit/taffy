/*
 * taf_summary: bigMafSummary-equivalent zoom-out coverage from the universal
 * .tui.  For a reference window, emit one BED-like row per (query species,
 * alignment chain): the reference interval the species covers, a coverage-
 * fraction score (aligned bp / interval span), and gap-derived
 * leftStatus/rightStatus synteny codes.  This replaces UCSC bigMafSummary's
 * role (the >1Mb browser zoom-out track) by querying the .tui at runtime
 * instead of reading a precomputed .summary.bb.
 *
 * It is a thin wrapper over the blockViz COVERAGE path (taffyGetBlocksInTarget
 * Range): we loop the query species, run the existing coverage engine over the
 * window for each, and reshape its per-chain blocks into bigMafSummary rows.
 * No score is stored in the .tui -- coverage fraction is derived from the chain
 * geometry the engine already computes (see docs/tui_summary_design.md).
 *
 * Output (tab-separated, bed3+4 like mafSummary.as; NOT sorted -- pipe through
 * `sort -k1,1 -k2,2n` if a sorted stream is needed):
 *   chrom  start  end  src  score  leftStatus  rightStatus
 *
 *  Released under the MIT license, see LICENSE.txt
 */

extern "C" {
#include "taf.h"
#include "sonLib.h"
#include "tui.h"
}
#include "taffyBlockViz.h"
#include <getopt.h>
#include <time.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// Status-code gap threshold: a gap to the adjacent same-species chain larger
// than this is a break ('N'); within it is an insert ('I'); zero is contig
// ('C').  Mirrors UCSC's 100kb NEW threshold (kent maf.h).
static const int64_t STATUS_BREAK_GAP = 100000;

static void usage() {
    fprintf(stderr, "taffy summary [options]\n");
    fprintf(stderr, "bigMafSummary-equivalent zoom-out coverage from a universal .tui: one BED row\n");
    fprintf(stderr, "per (query species, chain) over a reference window, scored by coverage fraction.\n\n");
    fprintf(stderr, "-i --inputFile FILE   Universal MAF/TAF (.tui resolved alongside) or a .tui directly\n");
    fprintf(stderr, "                      (incl. a chained .tui -- recommended for zoom-out).\n");
    fprintf(stderr, "-R --refGenome NAME   Reference genome (the window's coordinate system; tSpecies).\n");
    fprintf(stderr, "-r --region REGION    chrom:start-end (bare chrom, 0-based half-open) on the reference.\n");
    fprintf(stderr, "-q --query CSV        Comma-list of query species to summarize. Default: every\n");
    fprintf(stderr, "                      genome in the .tui except the reference.\n");
    fprintf(stderr, "-m --maxOutputBlocks N  Coverage bin budget per species (default engine value).\n");
    fprintf(stderr, "-o --outputFile FILE  Output (default stdout).\n");
    fprintf(stderr, "-t --time             Print per-species + total query wall time to stderr (de-risk).\n");
    fprintf(stderr, "-S --sweep            All-species coverage in ONE column-sweep (low-level tui API,\n");
    fprintf(stderr, "                      no per-species chaining); emits chrom start end src frac.\n");
    fprintf(stderr, "-T --threads N        (--sweep) parallel genome loads, one Tui handle per thread.\n");
    fprintf(stderr, "-P --sample N         Forward block-sampling: read the whole block (all species) at\n");
    fprintf(stderr, "                      up to N X-index-anchor columns; cost scales with N, not species.\n");
    fprintf(stderr, "-b --bins N           Per-pixel COVERAGE (bigMafSummary-matching): N reference bins,\n");
    fprintf(stderr, "                      coverage fraction per (species,bin); emits bed3+4 with status.\n");
    fprintf(stderr, "-C --sampleCols S     (--bins) estimate coverage from S columns/anchor (direct seeks)\n");
    fprintf(stderr, "                      instead of the full window -- cheap; 0 = exact full-window.\n");
    fprintf(stderr, "-l --logLevel LEVEL   Set the log level.\n");
    fprintf(stderr, "-h --help             Print this help message.\n");
}

// A reshaped chain: reference span + covered bp + strand, accumulated from the
// engine's mappedBlocks grouped by chainId.
struct ChainAgg {
    int64_t start = INT64_MAX;
    int64_t end = 0;
    int64_t covered_bp = 0;
    char strand = '+';
    bool seen = false;
};

static double now_seconds() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + (double) ts.tv_nsec * 1e-9;
}

// ---------------------------------------------------------------------------
// --sweep: all-species coverage in ONE pass over the reference window's
// universal columns, via the low-level tui API (no per-species tui_query, no
// chaining, no blockViz mutex).  Prototype to test whether amortizing those
// costs collapses the serial per-species latency (the design risk).
//
// For coverage we only need, per species, how many of the window's universal
// columns it aligns to: the reference is row-0 so each window column == one
// reference bp, hence covered-columns == covered-reference-bp.  (This v1 SUMS
// run overlaps, so paralogs that revisit a column over-count slightly; the
// fraction is clamped to 1.0.  A faithful version would union the columns.)
// ---------------------------------------------------------------------------
struct SweepCtx { int64_t c_lo; int64_t c_hi; int64_t covered; };

static void sweep_cb(const TuiRun *run, void *user) {
    SweepCtx *cx = (SweepCtx *) user;
    int64_t a = run->g_start, b = run->g_start + run->length;
    int64_t lo = a > cx->c_lo ? a : cx->c_lo;
    int64_t hi = b < cx->c_hi ? b : cx->c_hi;
    if (hi > lo) cx->covered += (hi - lo);
}

static int run_sweep(const char *input, const char *ref_genome, const string &chrom,
                     int64_t t_start, int64_t t_end, const vector<string> *query_subset,
                     FILE *out, bool do_time, int n_threads) {
    // tui_load wants the .tui path itself (taffyOpen resolves the sibling, this
    // doesn't): accept either a .tui or its companion MAF/TAF.
    size_t ilen = strlen(input);
    char *tp = (ilen >= 4 && strcmp(input + ilen - 4, ".tui") == 0) ? NULL : tui_path(input);
    string tui_file = tp ? tp : input;
    free(tp);

    Tui *tui = tui_load(tui_file.c_str());
    if (tui == NULL) { fprintf(stderr, "ERROR: tui_load(%s) failed\n", tui_file.c_str()); return 1; }

    string ref_seq = string(ref_genome) + "." + chrom;
    int64_t ref_len = tui_seq_length(tui, ref_seq.c_str());
    if (ref_len < 0) { fprintf(stderr, "ERROR: ref seq %s not in .tui\n", ref_seq.c_str()); tui_destruct(tui); return 1; }
    if (t_end == 0) t_end = ref_len;

    // Reference window -> universal-column intervals (ONCE, shared by all species;
    // plain ints, so it outlives this Tui and is read-only across threads).
    int64_t n_iv = 0;
    TuiInterval *iv = tui_query(tui, ref_seq.c_str(), t_start, t_end, &n_iv);
    int64_t window_cols = 0;
    for (int64_t i = 0; i < n_iv; i++) window_cols += iv[i].end - iv[i].start;

    // Genome list: explicit subset, or the whole roster minus the reference.
    vector<string> genomes;
    if (query_subset != NULL) {
        genomes = *query_subset;
    } else {
        int64_t ngg = 0;
        TuiGenomeInfo *gi = tui_genome_names(tui, &ngg);
        for (int64_t i = 0; i < ngg; i++) if (strcmp(gi[i].name, ref_genome) != 0) genomes.push_back(gi[i].name);
        tui_genome_info_free(gi, ngg);
    }
    tui_destruct(tui);   // done with the shared handle; threads each load their own
    if (window_cols == 0) {
        if (do_time) fprintf(stderr, "sweep: reference window has no columns\n");
        free(iv); return 0;
    }

    const int64_t ng = (int64_t) genomes.size();
    vector<double> frac(ng, -1.0);   // -1 = no coverage (skipped on emit)

    // Per-genome lift-table loading is the bottleneck (random I/O into the big
    // .tui), NOT query/chaining -- so the win is parallelism.  The blockViz API
    // serializes on a per-handle mutex; this low-level sweep doesn't, so each
    // thread opens its OWN Tui* (cursors are not shareable) and the per-genome
    // loads run concurrently.
#ifdef _OPENMP
    if (n_threads < 1) n_threads = 1;
#else
    n_threads = 1;
#endif
    double t0 = now_seconds();
#ifdef _OPENMP
    #pragma omp parallel num_threads(n_threads)
#endif
    {
        Tui *lt = tui_load(tui_file.c_str());
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 8)
#endif
        for (int64_t k = 0; k < ng; k++) {
            if (lt == NULL) continue;
            TuiGenomeLift *gl = tui_genome_lift_load(lt, genomes[k].c_str());
            if (gl == NULL) continue;
            int64_t covered = 0;
            for (int64_t i = 0; i < n_iv; i++) {
                SweepCtx cx = { iv[i].start, iv[i].end, 0 };
                tui_genome_lift_visit_runs(gl, iv[i].start, iv[i].end, sweep_cb, &cx);
                covered += cx.covered;
            }
            tui_genome_lift_destruct(gl);
            if (covered > 0) {
                double f = (double) covered / (double) window_cols;
                frac[k] = f > 1.0 ? 1.0 : f;
            }
        }
        if (lt != NULL) tui_destruct(lt);
    }
    double total = now_seconds() - t0;

    int64_t n_rows = 0, n_empty = 0;
    for (int64_t k = 0; k < ng; k++) {
        if (frac[k] < 0.0) { n_empty++; continue; }
        fprintf(out, "%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%.4g\n",
                chrom.c_str(), t_start, t_end, genomes[k].c_str(), frac[k]);
        n_rows++;
    }
    if (do_time) {
        fprintf(stderr, "sweep: %" PRIi64 " species, %" PRIi64 " rows (%" PRIi64 " empty), %" PRIi64 " window-cols"
                        " in %.3fs (%.1f ms/species, %d threads)\n",
                ng, n_rows, n_empty, window_cols, total,
                ng ? 1000.0 * total / (double) ng : 0.0, n_threads);
    }
    free(iv);
    return 0;
}

// ---------------------------------------------------------------------------
// --sample: FORWARD block-sampling.  The universal MAF/TAF is in column order
// (== file order) and every block carries ALL species, so we sample ~N columns
// across the window and read the whole block at each -- per-pixel coverage in
// ~N block reads, INDEPENDENT of species count (vs the reverse-lift sweep's
// per-genome floor).  We sample at the .tui X-index anchors (column->file-offset,
// ~every 10k cols) so each read is a direct seek, not a scan.
// ---------------------------------------------------------------------------
static int run_sample(const char *input, const char *ref_genome, const string &chrom,
                      int64_t t_start, int64_t t_end, int64_t n_want, const char *maf_override,
                      const vector<string> *query_subset, FILE *out, bool do_time) {
    const char *extract_path = maf_override ? maf_override : input;
    size_t ilen = strlen(input);
    char *tp = (ilen >= 4 && strcmp(input + ilen - 4, ".tui") == 0) ? NULL : tui_path(input);
    string tui_file = tp ? tp : input;
    free(tp);
    Tui *tui = tui_load(tui_file.c_str());
    if (tui == NULL) { fprintf(stderr, "ERROR: tui_load(%s) failed\n", tui_file.c_str()); return 1; }

    string ref_seq = string(ref_genome) + "." + chrom;
    int64_t ref_len = tui_seq_length(tui, ref_seq.c_str());
    if (ref_len < 0) { fprintf(stderr, "ERROR: ref seq %s not in .tui\n", ref_seq.c_str()); tui_destruct(tui); return 1; }
    if (t_end == 0) t_end = ref_len;
    int64_t n_uiv = 0;
    TuiInterval *uiv = tui_query(tui, ref_seq.c_str(), t_start, t_end, &n_uiv);
    if (uiv == NULL || n_uiv == 0) { fprintf(stderr, "sample: window maps to no columns\n"); tui_destruct(tui); return 0; }
    int64_t win_lo = uiv[0].start, win_hi = uiv[n_uiv - 1].end;

    // Pick sample columns at X-index anchors inside the window (direct seeks);
    // stride down to <= n_want.  A small helper: a column is "in window" if it
    // lands in one of the reference's uiv column intervals.
    int64_t na = tui_idx_n(tui);
    const int64_t *acol = tui_idx_cols(tui);
    vector<int64_t> anchors;
    for (int64_t i = 0; i < na; i++) {
        int64_t c = acol[i];
        if (c < win_lo || c >= win_hi) continue;
        for (int64_t j = 0; j < n_uiv; j++) if (c >= uiv[j].start && c < uiv[j].end) { anchors.push_back(c); break; }
    }
    // If anchors are sparser than requested, just use them; if denser, stride.
    vector<TuiInterval> sample_iv;
    int64_t have = (int64_t) anchors.size();
    int64_t stride = (n_want > 0 && have > n_want) ? (have / n_want) : 1;
    for (int64_t i = 0; i < have; i += stride) {
        TuiInterval s = { anchors[i], anchors[i] + 1, 0, 0 };
        sample_iv.push_back(s);
    }

    // genome_name_map for extract_genome_name (resolves dotted accession names; a map
    // returns NULL on a miss, whereas a hal_species SET st_errAbort()s the process).
    int64_t ng = 0;
    TuiGenomeInfo *gi = tui_genome_names(tui, &ng);
    stHash *gmap = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    for (int64_t i = 0; i < ng; i++) { char *nm = stString_copy(gi[i].name); stHash_insert(gmap, nm, nm); }
    tui_genome_info_free(gi, ng);

    // Open the underlying MAF/TAF for forward extraction (may differ from -i:
    // e.g. local .tui for the load + a remote MAF for the block reads).
    FILE *fh = fopen(extract_path, "r");
    if (fh == NULL) { fprintf(stderr, "ERROR: cannot open %s for extraction (needs the real MAF/TAF, not a stub)\n", extract_path); stHash_destruct(gmap); free(uiv); tui_destruct(tui); return 1; }
    LI *li = LI_construct(fh);
    int fmt = check_input_format(LI_peek_at_next_line(li));
    bool maf_input = (fmt == 1);
    bool rle = false;
    Tag *hdr = maf_input ? maf_read_header(li) : taf_read_header(li);
    if (!maf_input) { Tag *t = tag_find(hdr, "run_length_encode_bases"); if (t && strcmp(t->value, "1") == 0) rle = true; }
    if (hdr != NULL) tag_destruct(hdr);   // only the rle flag is needed

    map<string, int64_t> present;   // genome -> #sampled blocks it appears in
    int64_t n_blocks = 0;

    double t0 = now_seconds();
    TuiExtractIt *xit = tui_extract_iterator(tui, li, maf_input, rle, sample_iv.data(), (int64_t) sample_iv.size());
    Alignment *aln = NULL;
    while ((aln = tui_extract_next(xit, li)) != NULL) {
        n_blocks++;
        // Dedup genomes within the block (paralogs share a genome).
        set<string> here;
        for (Alignment_Row *r = aln->row; r != NULL; r = r->n_row) {
            char *g = extract_genome_name(r->sequence_name, NULL, gmap);
            if (g) { here.insert(g); free(g); }
        }
        for (const string &g : here) present[g]++;
    }
    tui_extract_iterator_destruct(xit);
    double total = now_seconds() - t0;

    // Emit per-genome presence rate (= sampled coverage fraction estimate).
    set<string> qset;
    if (query_subset) for (const string &s : *query_subset) qset.insert(s);
    int64_t n_rows = 0;
    for (auto &kv : present) {
        if (kv.first == ref_genome) continue;
        if (query_subset && qset.count(kv.first) == 0) continue;   // honor -q/--query
        double frac = n_blocks ? (double) kv.second / (double) n_blocks : 0.0;
        fprintf(out, "%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%.4g\n", chrom.c_str(), t_start, t_end, kv.first.c_str(), frac);
        n_rows++;
    }
    if (do_time) {
        fprintf(stderr, "sample: %" PRIi64 " sample-cols (of %" PRIi64 " anchors in window), %" PRIi64 " blocks read, "
                        "%" PRIi64 " species-rows in %.3fs (%.2f ms/block)\n",
                (int64_t) sample_iv.size(), have, n_blocks, n_rows, total,
                n_blocks ? 1000.0 * total / (double) n_blocks : 0.0);
    }
    LI_destruct(li);
    fclose(fh);
    stHash_destruct(gmap);
    free(uiv);
    tui_destruct(tui);
    return 0;
}

// ---------------------------------------------------------------------------
// --bins: the real bigMafSummary-matching output -- per-pixel COVERAGE (not
// presence).  Forward-extract every block in the window; for each block, find
// the reference row (its [start,start+length) is the block's reference span)
// and attribute that span to every non-reference species present, accumulated
// into N reference bins.  coverage[species][bin] = covered_ref_bp / bin_ref_bp.
// Emit bed3+4 per (species, bin) with coverage score + adjacent-bin status.
//
// Block-extent attribution (species present in a block "covers" the block's
// reference span) is a coarse approximation -- within-block bp gaps are below
// the ~bin resolution, fine for a zoom-out overview.  Reads the whole window
// (exact, not sampled); the sampled-cheap variant is future work.
// ---------------------------------------------------------------------------
static int run_coverage(const char *input, const char *ref_genome, const string &chrom,
                        int64_t t_start, int64_t t_end, int64_t n_bins, int64_t sample_cols,
                        int n_threads, const char *maf_override, const vector<string> *query_subset,
                        FILE *out, bool do_time) {
    const char *extract_path = maf_override ? maf_override : input;
    size_t ilen = strlen(input);
    char *tp = (ilen >= 4 && strcmp(input + ilen - 4, ".tui") == 0) ? NULL : tui_path(input);
    string tui_file = tp ? tp : input;
    free(tp);
    Tui *tui = tui_load(tui_file.c_str());
    if (tui == NULL) { fprintf(stderr, "ERROR: tui_load(%s) failed\n", tui_file.c_str()); return 1; }

    string ref_seq = string(ref_genome) + "." + chrom;
    // Validate the reference seq ALWAYS (not only when t_end==0): with an explicit
    // end a bogus -R/chrom otherwise slips through to an empty tui_query and exit 0.
    int64_t ref_len = tui_seq_length(tui, ref_seq.c_str());
    if (ref_len < 0) { fprintf(stderr, "ERROR: ref seq %s not in .tui\n", ref_seq.c_str()); tui_destruct(tui); return 1; }
    if (t_end == 0) t_end = ref_len;
    int64_t span = t_end - t_start;
    if (n_bins < 1) n_bins = 1;
    if (n_bins > span) n_bins = span;   // no more bins than reference bp

    int64_t n_uiv = 0;
    TuiInterval *uiv = tui_query(tui, ref_seq.c_str(), t_start, t_end, &n_uiv);
    if (uiv == NULL || n_uiv == 0) { fprintf(stderr, "coverage: window maps to no columns\n"); tui_destruct(tui); return 0; }

    int64_t ng = 0;
    TuiGenomeInfo *gi = tui_genome_names(tui, &ng);
    // genome_name_map (not a hal_species SET) so extract_genome_name returns NULL on
    // a miss instead of st_errAbort()ing the process; value=key (non-NULL), freed once.
    stHash *gmap = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    for (int64_t i = 0; i < ng; i++) { char *nm = stString_copy(gi[i].name); stHash_insert(gmap, nm, nm); }
    tui_genome_info_free(gi, ng);

    FILE *fh = fopen(extract_path, "r");
    if (fh == NULL) { fprintf(stderr, "ERROR: cannot open %s for extraction\n", extract_path); stHash_destruct(gmap); free(uiv); tui_destruct(tui); return 1; }
    LI *li = LI_construct(fh);
    int fmt = check_input_format(LI_peek_at_next_line(li));
    bool maf_input = (fmt == 1);
    bool rle = false;
    Tag *hdr = maf_input ? maf_read_header(li) : taf_read_header(li);
    if (!maf_input) { Tag *t = tag_find(hdr, "run_length_encode_bases"); if (t && strcmp(t->value, "1") == 0) rle = true; }
    if (hdr != NULL) tag_destruct(hdr);   // only the rle flag is needed

    // Extraction intervals: the full window (exact) or thin slices at the X-index
    // anchors (sampled).  CAVEAT on the sampled estimate: the number of sampling
    // LOCATIONS is the fixed X-index anchor set (~1 per TUI_IDX_BLOCK=10k columns),
    // INDEPENDENT of sample_cols -- so convergence to exact is governed by
    // sample_cols / bin-width, NOT by sample_cols.  With ~1 anchor per bin the per-bin
    // estimate is noisy (~0.09 stdev observed) and bins with no in-window anchor are
    // dropped.  Good for a cheap >>1Mb preview; use sample_cols=0 (exact) when per-bin
    // fidelity matters.  NOTE: tui_query above (line ~362) and this anchor scan are
    // WINDOW-bound, not pixel-bound -- only the block reads scale with the sample count.
    const TuiInterval *xiv = uiv;
    int64_t n_xiv = n_uiv;
    vector<TuiInterval> slices;
    if (sample_cols > 0) {
        int64_t win_lo = uiv[0].start, win_hi = uiv[n_uiv - 1].end;
        int64_t na = tui_idx_n(tui);
        const int64_t *acol = tui_idx_cols(tui);   // sorted ascending
        // Binary-search the window's anchor range [win_lo, win_hi) (the anchor
        // array is genome-wide -- ~7M entries at 577 -- so a linear scan per
        // query would dominate).
        const int64_t *lo = std::lower_bound(acol, acol + na, win_lo);
        const int64_t *hi = std::lower_bound(acol, acol + na, win_hi);
        // Keep only anchors that fall inside a REFERENCE column-interval (uiv):
        // in a universal MAF most columns gap the reference (insertions in other
        // genomes), and a gap-column slice samples no reference position -- so
        // sampling there both wastes the read and starves real bins.  Merge the
        // two sorted lists (anchors, uiv) -- O(anchors_in_window + n_uiv).
        vector<int64_t> kept;
        int64_t j = 0;
        for (const int64_t *p = lo; p < hi; ++p) {
            int64_t a = *p;
            while (j < n_uiv && uiv[j].end <= a) j++;
            if (j < n_uiv && a >= uiv[j].start && a < uiv[j].end) kept.push_back(a);
        }
        // Stride anchors down to ~n_bins (the pixel count): the display is fixed
        // at ~N px regardless of window, so cost stays ∝ pixels, not window size.
        int64_t have = (int64_t) kept.size();
        int64_t stride = (have > n_bins && n_bins > 0) ? (have / n_bins) : 1;
        for (int64_t i = 0; i < have; i += stride) {
            int64_t a = kept[i];
            int64_t e = a + sample_cols; if (e > win_hi) e = win_hi;
            TuiInterval s = { a, e, 0, 0 };
            slices.push_back(s);
        }
        // tui_extract_iterator requires sorted + MERGED intervals; slices are
        // sorted (anchors ascending) but can overlap when a strided gap is
        // smaller than sample_cols, so coalesce abutting/overlapping ones.
        if (!slices.empty()) {
            vector<TuiInterval> merged;
            merged.push_back(slices[0]);
            for (size_t i = 1; i < slices.size(); i++) {
                if (slices[i].start <= merged.back().end) {
                    if (slices[i].end > merged.back().end) merged.back().end = slices[i].end;
                } else {
                    merged.push_back(slices[i]);
                }
            }
            slices.swap(merged);
        }
        xiv = slices.data();
        n_xiv = (int64_t) slices.size();
        if (n_xiv == 0)
            fprintf(stderr, "WARNING: no X-index anchors fell in the window; --sampleCols produced"
                            " nothing (window span < ~anchor gap) -- widen the window or use"
                            " --sampleCols 0 (exact).\n");
    }

    vector<int64_t> bin_ref(n_bins, 0);                  // reference bp per bin
    map<string, vector<int64_t>> cov;                    // species -> covered ref bp per bin
    // Add ref interval [a,b) to a per-bin accumulator (overlap-weighted).  Start at
    // floor((a-t_start)*n_bins/span), which is provably <= the bin actually
    // containing `a`, then WALK forward (clipping per bin) and stop once a bin starts
    // at/after b.  The old upper bound b1=floor((b-1-t_start)*n_bins/span) is NOT the
    // inverse of the lo=t_start+bin*span/n_bins boundaries when span%n_bins!=0, so it
    // under-shot the top bin and dropped a span's high-end tail into no bin at all.
    auto add_span = [&](vector<int64_t> &arr, int64_t a, int64_t b) {
        if (a < t_start) a = t_start;
        if (b > t_end) b = t_end;
        if (b <= a) return;
        int64_t bin = (a - t_start) * n_bins / span;
        if (bin < 0) bin = 0;
        for (; bin < n_bins; bin++) {
            int64_t lo = t_start + bin * span / n_bins;
            if (lo >= b) break;
            int64_t hi = t_start + (bin + 1) * span / n_bins;
            int64_t ov = (b < hi ? b : hi) - (a > lo ? a : lo);
            if (ov > 0) arr[bin] += ov;
        }
    };

    // Per-extractor accumulation into the given (bin_ref, cov) pair -- reused by
    // the serial path and by each parallel worker (with its own local arrays).
    // t_ex/t_pr accumulate the time spent in block-DECODE (tui_extract_next) vs
    // ROW-PROCESSING (extract_genome_name + add_span) so -t can show which side
    // dominates at scale -- the bottleneck isn't measurable on the small local set.
    // thr_gmap is a genome_name_map (NOT a hal_species SET): with a map,
    // extract_genome_name RETURNS NULL on an unresolvable name; with a set it
    // st_errAbort()s -- i.e. kills the whole process mid-stream on one odd row.
    auto accumulate = [&](TuiExtractIt *xit, LI *lli, stHash *thr_gmap, vector<int64_t> &lbin,
                          map<string, vector<int64_t>> &lcov, double &t_ex, double &t_pr) -> int64_t {
        int64_t nb = 0;
        vector<pair<int64_t, int64_t>> ref_spans;   // forward [s,e) of the block's reference row(s)
        set<string> here;
        while (true) {
            double a = now_seconds();
            Alignment *aln = tui_extract_next(xit, lli);
            t_ex += now_seconds() - a;
            if (aln == NULL) break;
            double b = now_seconds();
            nb++;
            ref_spans.clear();
            here.clear();
            for (Alignment_Row *r = aln->row; r != NULL; r = r->n_row) {
                char *g = extract_genome_name(r->sequence_name, NULL, thr_gmap);
                if (!g) continue;   // sequence whose genome prefix isn't in the .tui roster
                if (strcmp(g, ref_genome) == 0) {
                    // r->start is MAF-native: forward on '+', but srcSize-relative on
                    // '-'.  Convert to the FORWARD reference span so it bins correctly --
                    // a leaf reference is '-' in ~half its rows, and those were being
                    // silently clipped to nothing (read as forward, they fall outside
                    // the window).  Union (not last-wins) if the reference dups.
                    int64_t s, e;
                    if (r->strand) { s = r->start; e = r->start + r->length; }
                    else { s = r->sequence_length - r->start - r->length; e = r->sequence_length - r->start; }
                    if (e > s) ref_spans.push_back({s, e});
                } else {
                    here.insert(g);
                }
                free(g);
            }
            for (const auto &sp : ref_spans) {   // attribute each reference span (usually exactly one)
                add_span(lbin, sp.first, sp.second);
                for (const string &g : here) {
                    auto it = lcov.find(g);
                    if (it == lcov.end()) it = lcov.emplace(g, vector<int64_t>(n_bins, 0)).first;
                    add_span(it->second, sp.first, sp.second);
                }
            }
            t_pr += now_seconds() - b;
        }
        return nb;
    };

    int64_t n_blocks = 0;
    double t_extract = 0, t_process = 0;   // decode vs row-processing (thread-summed in parallel)
    double t0 = now_seconds();
#ifdef _OPENMP
    // Parallelize the SAMPLED path: the anchor slices are independent, so split
    // them across threads -- each opens its OWN Tui + LI + extractor (the .tui
    // cursor is not shareable) and accumulates locally, then merges.
    if (sample_cols > 0 && n_threads > 1 && n_xiv > 1) {
        #pragma omp parallel num_threads(n_threads) reduction(+:n_blocks)
        {
            int tid = omp_get_thread_num(), nth = omp_get_num_threads();
            int64_t per = (n_xiv + nth - 1) / nth;
            int64_t s0 = (int64_t) tid * per, s1 = s0 + per;
            if (s1 > n_xiv) s1 = n_xiv;
            if (s0 < s1) {
                Tui *lt = tui_load(tui_file.c_str());
                FILE *lfh = fopen(extract_path, "r");
                LI *lli = lfh ? LI_construct(lfh) : NULL;
                if (lt != NULL && lli != NULL) {
                    int lf = check_input_format(LI_peek_at_next_line(lli));
                    Tag *lhdr = (lf == 1) ? maf_read_header(lli) : taf_read_header(lli);
                    if (lhdr != NULL) tag_destruct(lhdr);
                    // Per-thread genome set: extract_genome_name does a (read-only
                    // but unsynchronized) stSet_search, so each thread gets its own
                    // set rather than racing on a shared one.
                    int64_t lng = 0;
                    TuiGenomeInfo *lgi = tui_genome_names(lt, &lng);
                    stHash *lgmap = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
                    for (int64_t i = 0; i < lng; i++) { char *nm = stString_copy(lgi[i].name); stHash_insert(lgmap, nm, nm); }
                    tui_genome_info_free(lgi, lng);
                    vector<int64_t> lbin(n_bins, 0);
                    map<string, vector<int64_t>> lcov;
                    double ltex = 0, ltpr = 0;
                    TuiExtractIt *lxit = tui_extract_iterator(lt, lli, maf_input, rle, &xiv[s0], s1 - s0);
                    n_blocks += accumulate(lxit, lli, lgmap, lbin, lcov, ltex, ltpr);
                    tui_extract_iterator_destruct(lxit);
                    stHash_destruct(lgmap);
                    #pragma omp critical
                    {
                        t_extract += ltex; t_process += ltpr;
                        for (int64_t i = 0; i < n_bins; i++) bin_ref[i] += lbin[i];
                        for (auto &kv : lcov) {
                            auto &dst = cov[kv.first];
                            if ((int64_t) dst.size() != n_bins) dst.assign(n_bins, 0);
                            for (int64_t i = 0; i < n_bins; i++) dst[i] += kv.second[i];
                        }
                    }
                }
                if (lli != NULL) LI_destruct(lli);
                if (lfh != NULL) fclose(lfh);
                if (lt != NULL) tui_destruct(lt);
            }
        }
    } else
#endif
    {
        TuiExtractIt *xit = tui_extract_iterator(tui, li, maf_input, rle, xiv, n_xiv);
        n_blocks = accumulate(xit, li, gmap, bin_ref, cov, t_extract, t_process);
        tui_extract_iterator_destruct(xit);
    }
    double total = now_seconds() - t0;

    // Emit per (species, bin) coverage where the reference bin has data; status
    // codes from the adjacent same-species bin (C contiguous / N break; M omitted).
    set<string> qset;
    if (query_subset) for (const string &s : *query_subset) qset.insert(s);
    int64_t n_rows = 0;
    for (auto &kv : cov) {
        if (kv.first == ref_genome) continue;                       // never the reference
        if (query_subset && qset.count(kv.first) == 0) continue;    // honor -q/--query
        const vector<int64_t> &c = kv.second;
        for (int64_t bin = 0; bin < n_bins; bin++) {
            if (bin_ref[bin] <= 0 || c[bin] <= 0) continue;
            double frac = (double) c[bin] / (double) bin_ref[bin];
            if (frac > 1.0) frac = 1.0;
            int64_t lo = t_start + bin * span / n_bins;
            int64_t hi = t_start + (bin + 1) * span / n_bins;
            // 'C' only if the adjacent bin is ALSO emitted (has both ref bp AND
            // coverage); a skipped (bin_ref==0) neighbour is a break, not contig.
            char left  = (bin > 0 && bin_ref[bin - 1] > 0 && c[bin - 1] > 0) ? 'C' : 'N';
            char right = (bin + 1 < n_bins && bin_ref[bin + 1] > 0 && c[bin + 1] > 0) ? 'C' : 'N';
            fprintf(out, "%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%.4g\t%c\t%c\n",
                    chrom.c_str(), lo, hi, kv.first.c_str(), frac, left, right);
            n_rows++;
        }
    }
    if (do_time) {
        int eff_threads = (sample_cols > 0 ? n_threads : 1);
        fprintf(stderr, "coverage: %" PRIi64 " bins, %" PRIi64 " blocks (%s%" PRIi64 " intervals), %zu species, %" PRIi64 " rows, %d threads in %.3fs\n",
                n_bins, n_blocks, sample_cols > 0 ? "sampled " : "full ", n_xiv, cov.size(), n_rows,
                eff_threads, total);
        // Thread-summed: divide by threads for the wall-clock share.  Shows whether
        // block-DECODE (taf column walk) or ROW-PROCESSING (extract_genome_name +
        // add_span) is the bottleneck -- the thing to optimize.
        fprintf(stderr, "  profile: decode %.2fs  process %.2fs  (thread-summed; /%d for wall share)\n",
                t_extract, t_process, eff_threads);
    }
    LI_destruct(li);
    fclose(fh);
    stHash_destruct(gmap);
    free(uiv);
    tui_destruct(tui);
    return 0;
}

int taf_summary_main(int argc, char *argv[]) {
    char *input = NULL;
    char *ref_genome = NULL;
    char *region = NULL;
    char *query_csv = NULL;
    char *output = NULL;
    int64_t max_output_blocks = -1;   // -1 = leave engine default
    bool do_time = false;
    bool sweep = false;
    int n_threads = 1;
    int64_t n_sample = 0;   // >0 = forward block-sampling mode (N samples)
    int64_t n_bins = 0;     // >0 = per-bin coverage mode (the real bigMafSummary output)
    int64_t sample_cols = 0;   // --sampleCols S: sample S cols/anchor instead of full window
    char *maf_override = NULL;   // --maf: extraction source != -i (local .tui + remote MAF)
    char *log_level = NULL;

    while (true) {
        static struct option opts[] = {
            {"inputFile", required_argument, 0, 'i'},
            {"refGenome", required_argument, 0, 'R'},
            {"region", required_argument, 0, 'r'},
            {"query", required_argument, 0, 'q'},
            {"maxOutputBlocks", required_argument, 0, 'm'},
            {"outputFile", required_argument, 0, 'o'},
            {"time", no_argument, 0, 't'},
            {"sweep", no_argument, 0, 'S'},
            {"threads", required_argument, 0, 'T'},
            {"sample", required_argument, 0, 'P'},
            {"bins", required_argument, 0, 'b'},
            {"sampleCols", required_argument, 0, 'C'},
            {"maf", required_argument, 0, 'M'},
            {"logLevel", required_argument, 0, 'l'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        int oi = 0;
        int c = getopt_long(argc, argv, "i:R:r:q:m:o:tST:P:b:C:M:l:h", opts, &oi);
        if (c == -1) break;
        switch (c) {
            case 'i': input = optarg; break;
            case 'R': ref_genome = optarg; break;
            case 'r': region = optarg; break;
            case 'q': query_csv = optarg; break;
            case 'm': max_output_blocks = atoll(optarg); break;
            case 'o': output = optarg; break;
            case 't': do_time = true; break;
            case 'S': sweep = true; break;
            case 'T': n_threads = atoi(optarg); break;
            case 'P': n_sample = atoll(optarg); break;
            case 'b': n_bins = atoll(optarg); break;
            case 'C': sample_cols = atoll(optarg); break;
            case 'M': maf_override = optarg; break;
            case 'l': log_level = optarg; break;
            case 'h': usage(); return 0;
            default: usage(); return 1;
        }
    }
    if (log_level != NULL) st_setLogLevelFromString(log_level);

    if (input == NULL)      { fprintf(stderr, "ERROR: -i/--inputFile required\n"); usage(); return 1; }
    if (ref_genome == NULL) { fprintf(stderr, "ERROR: -R/--refGenome required\n"); usage(); return 1; }
    if (region == NULL)     { fprintf(stderr, "ERROR: -r/--region required\n"); usage(); return 1; }

    // Parse region "chrom:start-end" (chrom may itself contain no ':'; split on
    // the LAST ':' so it tolerates odd names, then the single '-' in the range).
    string reg(region);
    size_t colon = reg.rfind(':');
    if (colon == string::npos) { fprintf(stderr, "ERROR: --region must be chrom:start-end\n"); return 1; }
    string chrom = reg.substr(0, colon);
    string range = reg.substr(colon + 1);
    size_t dash = range.find('-');
    if (dash == string::npos) { fprintf(stderr, "ERROR: --region must be chrom:start-end\n"); return 1; }
    int64_t t_start = atoll(range.substr(0, dash).c_str());
    int64_t t_end = atoll(range.substr(dash + 1).c_str());
    if (t_start < 0 || t_end < 0) { fprintf(stderr, "ERROR: region coordinates must be >= 0\n"); return 1; }
    if (t_end != 0 && t_end <= t_start) { fprintf(stderr, "ERROR: region end must be > start (or 0 for chrom end)\n"); return 1; }
    if (n_threads < 1) n_threads = 1;
    if (n_bins < 0 || n_sample < 0 || sample_cols < 0) { fprintf(stderr, "ERROR: --bins/--sample/--sampleCols must be >= 0\n"); return 1; }
    if (((n_bins > 0) + (n_sample > 0) + (sweep ? 1 : 0)) > 1)
        fprintf(stderr, "WARNING: --bins/--sample/--sweep are mutually exclusive; precedence is bins>sample>sweep\n");

    // Optional explicit query subset (shared by both paths).
    vector<string> qsubset;
    bool have_subset = false;
    if (query_csv != NULL) {
        have_subset = true;
        string s(query_csv), cur;
        for (char ch : s) {
            if (ch == ',') { if (!cur.empty()) qsubset.push_back(cur); cur.clear(); }
            else cur.push_back(ch);
        }
        if (!cur.empty()) qsubset.push_back(cur);
    }

    FILE *out = output ? fopen(output, "w") : stdout;
    if (out == NULL) { fprintf(stderr, "ERROR: cannot open output %s\n", output); return 1; }

    // --bins: per-pixel coverage (the real bigMafSummary output).
    if (n_bins > 0) {
        int rc = run_coverage(input, ref_genome, chrom, t_start, t_end, n_bins, sample_cols, n_threads,
                              maf_override, have_subset ? &qsubset : NULL, out, do_time);
        if (out != stdout) fclose(out);
        return rc;
    }

    // --sample: forward block-sampling (reads the real MAF/TAF, all species per block).
    if (n_sample > 0) {
        int rc = run_sample(input, ref_genome, chrom, t_start, t_end, n_sample, maf_override,
                            have_subset ? &qsubset : NULL, out, do_time);
        if (out != stdout) fclose(out);
        return rc;
    }

    // --sweep: low-level all-species column-sweep (no blockViz handle).
    if (sweep) {
        int rc = run_sweep(input, ref_genome, chrom, t_start, t_end,
                           have_subset ? &qsubset : NULL, out, do_time, n_threads);
        if (out != stdout) fclose(out);
        return rc;
    }

    char *err = NULL;
    int h = taffyOpen(input, &err);
    if (h < 0) { fprintf(stderr, "ERROR: taffyOpen(%s) failed: %s\n", input, err ? err : "?"); free(err); if (out != stdout) fclose(out); return 1; }
    if (max_output_blocks > 0 && taffySetMaxOutputBlocks(h, max_output_blocks, &err) != 0) {
        fprintf(stderr, "ERROR: taffySetMaxOutputBlocks: %s\n", err ? err : "?"); free(err); taffyClose(h, NULL); if (out != stdout) fclose(out); return 1;
    }

    // Resolve the query-species list for the blockViz path.
    vector<string> species;
    if (have_subset) {
        species = qsubset;
    } else {
        struct taffy_species_t *sp = taffyGetSpecies(h, &err);
        if (sp == NULL && err != NULL) { fprintf(stderr, "ERROR: taffyGetSpecies: %s\n", err); free(err); taffyClose(h, NULL); if (out != stdout) fclose(out); return 1; }
        bool ref_found = false;
        for (struct taffy_species_t *p = sp; p != NULL; p = p->next) {
            if (strcmp(p->name, ref_genome) == 0) ref_found = true;
            else species.push_back(p->name);
        }
        taffyFreeSpeciesList(sp);
        if (!ref_found) { fprintf(stderr, "ERROR: reference genome %s not found in the .tui\n", ref_genome); taffyClose(h, NULL); if (out != stdout) fclose(out); return 1; }
    }

    double t0 = now_seconds();
    double slowest = 0.0; string slowest_sp;
    int64_t n_rows = 0, n_empty = 0, n_err = 0;

    for (const string &q : species) {
        double s0 = now_seconds();
        char *qerr = NULL;
        struct taffy_block_results_t *res = taffyGetBlocksInTargetRange(
            h, q.c_str(), ref_genome, chrom.c_str(), t_start, t_end,
            0, TAFFY_NO_SEQUENCES, TAFFY_QUERY_AND_TARGET_DUPS, 0, NULL, &qerr);
        double dt = now_seconds() - s0;
        if (dt > slowest) { slowest = dt; slowest_sp = q; }
        if (res == NULL) {
            // NULL = error (an empty list -- species simply not aligned in the
            // window -- comes back as a non-NULL result with mappedBlocks=NULL).
            if (qerr) { st_logInfo("summary: %s skipped: %s\n", q.c_str(), qerr); free(qerr); }
            n_err++;
            continue;
        }

        // Group the engine's mappedBlocks by chainId into ChainAgg.
        map<taffy_int_t, ChainAgg> chains;
        for (struct taffy_block_t *b = res->mappedBlocks; b != NULL; b = b->next) {
            ChainAgg &c = chains[b->chainId];
            c.seen = true;
            if (b->tStart < c.start) c.start = b->tStart;
            if (b->tStart + b->size > c.end) c.end = b->tStart + b->size;
            c.covered_bp += b->size;
            c.strand = b->strand;
        }
        if (chains.empty()) { n_empty++; taffyFreeBlockResults(res); continue; }

        // Order chains by reference start so we can derive left/right status
        // from the gap to the adjacent same-species chain.
        vector<ChainAgg> cv;
        cv.reserve(chains.size());
        for (auto &kv : chains) cv.push_back(kv.second);
        sort(cv.begin(), cv.end(), [](const ChainAgg &a, const ChainAgg &b) { return a.start < b.start; });

        for (size_t i = 0; i < cv.size(); i++) {
            const ChainAgg &c = cv[i];
            int64_t span = c.end - c.start;
            double score = span > 0 ? (double) c.covered_bp / (double) span : 0.0;
            if (score > 1.0) score = 1.0;   // paralog overlap safety clamp
            // Status vs the adjacent same-species chain (M = assembly-gap is not
            // derivable from the .tui -- no N-content -- so it is never emitted).
            char left = (i == 0) ? 'N' : 'C';
            char right = (i + 1 == cv.size()) ? 'N' : 'C';
            if (i > 0) {
                int64_t gap = c.start - cv[i - 1].end;
                left = gap < 0 ? 'T' : (gap == 0 ? 'C' : (gap > STATUS_BREAK_GAP ? 'N' : 'I'));
            }
            if (i + 1 < cv.size()) {
                int64_t gap = cv[i + 1].start - c.end;
                right = gap < 0 ? 'T' : (gap == 0 ? 'C' : (gap > STATUS_BREAK_GAP ? 'N' : 'I'));
            }
            fprintf(out, "%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%.4g\t%c\t%c\n",
                    chrom.c_str(), (int64_t) c.start, (int64_t) c.end, q.c_str(), score, left, right);
            n_rows++;
        }
        taffyFreeBlockResults(res);
    }

    if (out != stdout) fclose(out);
    double total = now_seconds() - t0;
    if (do_time) {
        fprintf(stderr, "summary: %zu species, %" PRIi64 " rows (%" PRIi64 " empty, %" PRIi64 " err) in %.3fs"
                        " (%.1f ms/species; slowest %s %.3fs)\n",
                species.size(), n_rows, n_empty, n_err, total,
                species.empty() ? 0.0 : 1000.0 * total / (double) species.size(),
                slowest_sp.c_str(), slowest);
    }
    taffyClose(h, NULL);
    return 0;
}
