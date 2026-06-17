/*
 * taffyBlockViz implementation -- thin C++ shim over taffy's existing
 * Tui / TuiGenomeLift primitives, exposing a C API shaped like HAL's
 * halBlockViz so a browser snake-track client can swap data sources
 * with minimal code change.  See taffyBlockViz.h for the public
 * contract and the list of unsupported features in this initial cut.
 */

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <vector>

#include "taffyBlockViz.h"

extern "C" {
#include "chain.h"
#include "sonLib.h"
#include "taf.h"
#include "tui.h"
}

/* ------------------------------------------------------------------ */
/* Handle state                                                        */
/* ------------------------------------------------------------------ */

// One row of a per-(genome, seq) runs table cached on the handle to
// drive mapBackAdjacencies' qSpecies-coord neighbor lookup.  Decoded
// out of tui_load_seq_runs's flat (t_start, g_start, lenc) triples;
// kept sorted by t_start (the load order) so a lower_bound bsearch
// over t_start finds the run owning a given qSpecies position.
struct SeqRun {
    int64_t t_start;   // qSpecies forward position
    int64_t g_start;   // universal column
    int64_t length;
    int     strand;    // +1 / -1
};

// Default cap on per-query mappedBlocks count.  Per-handle, tunable at
// runtime via taffySetMaxOutputBlocks.  Set to HAL's halSnakeTrack.c
// NUM_LEVELS=1000 -- the conservative 500 left dense viewports a few %
// short of full coverage (the chained .tui fragments even a near-
// identical genome into ~900 blocks per Mb), and snake levels track
// distinct chains, not raw block count, so 1000 blocks stays safe.
#define TAFFY_DEFAULT_MAX_OUTPUT_BLOCKS ((int64_t) 1000)

struct TaffyHandle {
    // Per-handle mutex.  All non-trivial public entries lock this for
    // the duration of their work.  Two threads on DIFFERENT handles run
    // fully concurrent (they only meet at the small g_table_mutex when
    // looking the handle up).  Two threads on the SAME handle still
    // serialize -- correct, since Tui*, TuiGenomeLift*, lift_cache and
    // qseq_runs_cache are NOT thread-safe across calls on one Tui*
    // (tui.h:97-104 and tui.h:200-205).  The kent halSnakeTrack flow
    // gives each track-thread its own taffyOpen handle, so per-handle
    // locking captures the full parallel-load win.
    std::mutex mu;
    Tui *tui = nullptr;
    std::string tui_path_str;   // resolved .tui path (for tui_sequence_lengths)
    std::map<std::string, TuiGenomeLift *> lift_cache;
    // Chain tuning -- defaults match TAFFY_CHAIN_DEFAULT_{OPEN,EXTEND,MAX_GAP}.
    // Tunable per-handle via taffySetChainParams.
    int64_t chain_open      = TAFFY_CHAIN_DEFAULT_OPEN;
    int64_t chain_extend    = TAFFY_CHAIN_DEFAULT_EXTEND;
    int64_t max_gap_length  = TAFFY_CHAIN_DEFAULT_MAX_GAP;
    // Overlap-frac paralogy filter on the chain output (see chain.h):
    // walk chains in score-desc order, accept iff their q-coverage
    // overlaps the running union of kept chains by at most this fraction
    // of the candidate's q-bp.  Default 0.5: a paralog whose q-coverage
    // is >50% redundant with a kept chain is dropped, while real
    // ortholog coverage survives.  Strict 0.0 is WRONG for the universal
    // .tui: it drops a whole chain for ANY (even 1-bp) overlap, and the
    // alignment is fragmented into many chains with tiny boundary
    // overlaps, so 0.0 culls ~17% of real hs1->hg38 coverage (measured).
    // 0.5 keeps that coverage AND still cuts SD-region paralog redundancy
    // ~8x.  Set to -1.0 to disable the filter entirely.
    double  chain_overlap_frac = 0.5;
    // Hard cap on mappedBlocks length per query.  Default is
    // TAFFY_DEFAULT_MAX_OUTPUT_BLOCKS (1000); tunable at runtime via
    // taffySetMaxOutputBlocks.
    int64_t max_output_blocks = TAFFY_DEFAULT_MAX_OUTPUT_BLOCKS;
    // Noise filter.  Drop an output block iff its size is BOTH <
    // min_block_span_frac of the query window AND < min_block_rel_frac of the
    // largest output block -- i.e. < min(window*span_frac, max*rel_frac)
    // ("below A and below B" == "below min(A,B)").  The window span stands in
    // for screen width (span_frac ~ 1/track-pixels), so this is a sub-pixel
    // test without a pixel count; the relative term is self-protecting (a
    // uniformly-small region keeps -- nothing is small relative to its own max
    // -- while slivers beside a real feature drop).  ON by default with
    // conservative fractions (~1/1000 px AND < 10% of the largest), so it
    // declutters out of the box; the browser can override per-track via
    // taffySetMinBlockFilter (e.g. exact 1/pixel-width).  Set BOTH to 0 to
    // disable; span_frac=1 makes it relative-only, rel_frac=1 window-only.
    double  min_block_span_frac = 0.001;
    double  min_block_rel_frac  = 0.1;
    // Adaptive run floor: drop input runs shorter than this (bp) BEFORE
    // chaining, on wide overview queries where they're sub-pixel.  -1 =
    // auto (span/500000, ~0 for narrow detail queries); >= 0 overrides.
    int64_t min_run_size = -1;
    // mapBackAdjacencies: per-(qSpecies, qChrom) sorted run table.
    // Key is the fully-qualified "<genome>.<chrom>" string the .tui
    // uses for d-line lookups.  Lazily populated on first flank scan
    // for each qChrom; browser pan/zoom keeps hitting the same set.
    std::map<std::string, std::vector<SeqRun>> qseq_runs_cache;
};

static std::map<int, TaffyHandle *> g_handles;
static int  g_next_handle = 1;
// Guards ONLY the g_handles map + g_next_handle.  Per-query work locks
// the handle's own H->mu instead, so two threads on different handles
// don't contend here except for the brief lookup.
static std::mutex g_table_mutex;

static void set_err(char **errStr, const char *msg) {
    if (errStr) *errStr = strdup(msg);
}

// Look up a handle's TaffyHandle* under g_table_mutex and return it.
// The caller then locks H->mu for the duration of its work.  Safe
// because TaffyHandle pointers are stable for the handle's lifetime
// (only taffyClose removes them, and it drains H->mu before delete).
static TaffyHandle *lookup_handle(int h, char **errStr) {
    std::lock_guard<std::mutex> lock(g_table_mutex);
    auto it = g_handles.find(h);
    if (it == g_handles.end()) {
        set_err(errStr, "taffyBlockViz: invalid handle");
        return nullptr;
    }
    return it->second;
}

// Cache the per-genome lift table across calls (browser pans + zooms
// hit the same target/query pair repeatedly).  Caller holds H->mu.
static TuiGenomeLift *get_or_load_gl(TaffyHandle *H, const std::string &genome,
                                     char **errStr) {
    auto it = H->lift_cache.find(genome);
    if (it != H->lift_cache.end()) return it->second;
    TuiGenomeLift *gl = tui_genome_lift_load(H->tui, genome.c_str());
    if (!gl) {
        std::string msg = "taffyBlockViz: genome '" + genome + "' not found in .tui";
        set_err(errStr, msg.c_str());
        return nullptr;
    }
    H->lift_cache[genome] = gl;
    return gl;
}

/* ------------------------------------------------------------------ */
/* Open / Close                                                        */
/* ------------------------------------------------------------------ */

extern "C" int taffyOpen(const char *path, char **errStr) {
    if (!path) { set_err(errStr, "taffyBlockViz: NULL path"); return -1; }
    // Resolve <foo>.maf.gz -> <foo>.maf.gz.tui, OR accept an already-
    // fully-qualified .tui path as-is.  tui_path() unconditionally appends
    // ".tui", so guard against a path that already ends in it (otherwise
    // an explicit ".tui" arg loads ".tui.tui" and fails).
    size_t plen = strlen(path);
    char *p = (plen >= 4 && strcmp(path + plen - 4, ".tui") == 0)
                  ? strdup(path)
                  : tui_path(path);
    Tui *tui = tui_load(p);
    if (!tui) {
        std::string msg = std::string("taffyBlockViz: tui_load failed for ") + (p ? p : path);
        set_err(errStr, msg.c_str());
        free(p);
        return -1;
    }

    TaffyHandle *H = new TaffyHandle();
    H->tui = tui;
    H->tui_path_str = p;
    free(p);
    std::lock_guard<std::mutex> lock(g_table_mutex);
    int h = g_next_handle++;
    g_handles[h] = H;
    return h;
}

extern "C" int taffyClose(int h, char **errStr) {
    // Erase from the handle table FIRST so no new query can find this
    // handle, then briefly take H->mu to drain any concurrent in-flight
    // query, then destruct.  This is the standard "remove-then-drain"
    // pattern; correct as long as callers don't call taffyClose racing
    // with their own ongoing query on the same handle.
    TaffyHandle *H = nullptr;
    {
        std::lock_guard<std::mutex> lock(g_table_mutex);
        auto it = g_handles.find(h);
        if (it == g_handles.end()) {
            set_err(errStr, "taffyBlockViz: invalid handle");
            return -1;
        }
        H = it->second;
        g_handles.erase(it);
    }
    { std::lock_guard<std::mutex> drain(H->mu); }   // wait for in-flight
    for (auto &kv : H->lift_cache) tui_genome_lift_destruct(kv.second);
    tui_destruct(H->tui);
    delete H;
    return 0;
}

extern "C" int taffyCloseGenome(int h, const char *genome, char **errStr) {
    if (!genome) { set_err(errStr, "taffyBlockViz: NULL genome"); return -1; }
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return -1;
    std::lock_guard<std::mutex> lock(H->mu);
    auto it = H->lift_cache.find(genome);
    if (it != H->lift_cache.end()) {
        tui_genome_lift_destruct(it->second);
        H->lift_cache.erase(it);
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Chain tuning setter/getter                                          */
/* ------------------------------------------------------------------ */

extern "C" int taffySetChainParams(int h,
                                   int64_t chain_open,
                                   int64_t chain_extend,
                                   int64_t max_gap_length,
                                   char **errStr) {
    // Validate: -1 means "leave unchanged"; any other negative is bogus.
    if (chain_open      < -1) { set_err(errStr, "taffyBlockViz: chain_open must be >= 0 or -1"); return -1; }
    if (chain_extend    < -1) { set_err(errStr, "taffyBlockViz: chain_extend must be >= 0 or -1"); return -1; }
    if (max_gap_length  < -1) { set_err(errStr, "taffyBlockViz: max_gap_length must be >= 0 or -1"); return -1; }
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return -1;
    std::lock_guard<std::mutex> lock(H->mu);
    if (chain_open     != -1) H->chain_open     = chain_open;
    if (chain_extend   != -1) H->chain_extend   = chain_extend;
    if (max_gap_length != -1) H->max_gap_length = max_gap_length;
    return 0;
}

extern "C" int taffyGetChainParams(int h,
                                   int64_t *chain_open,
                                   int64_t *chain_extend,
                                   int64_t *max_gap_length,
                                   char **errStr) {
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return -1;
    std::lock_guard<std::mutex> lock(H->mu);
    if (chain_open)     *chain_open     = H->chain_open;
    if (chain_extend)   *chain_extend   = H->chain_extend;
    if (max_gap_length) *max_gap_length = H->max_gap_length;
    return 0;
}

extern "C" int taffySetChainOverlapFrac(int h, double frac, char **errStr) {
    // -1.0 disables the filter; otherwise must be in [0.0, 1.0].
    if (!(frac == -1.0 || (frac >= 0.0 && frac <= 1.0))) {
        set_err(errStr, "taffyBlockViz: chain_overlap_frac must be -1 (off) or in [0.0, 1.0]");
        return -1;
    }
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return -1;
    std::lock_guard<std::mutex> lock(H->mu);
    H->chain_overlap_frac = frac;
    return 0;
}

extern "C" int taffyGetChainOverlapFrac(int h, double *frac, char **errStr) {
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return -1;
    std::lock_guard<std::mutex> lock(H->mu);
    if (frac) *frac = H->chain_overlap_frac;
    return 0;
}

extern "C" int taffySetMaxOutputBlocks(int h, int64_t n, char **errStr) {
    if (n < 1) {
        set_err(errStr, "taffyBlockViz: max_output_blocks must be >= 1");
        return -1;
    }
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return -1;
    std::lock_guard<std::mutex> lock(H->mu);
    H->max_output_blocks = n;
    return 0;
}

extern "C" int taffyGetMaxOutputBlocks(int h, int64_t *n, char **errStr) {
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return -1;
    std::lock_guard<std::mutex> lock(H->mu);
    if (n) *n = H->max_output_blocks;
    return 0;
}

extern "C" int taffySetMinBlockFilter(int h, double spanFrac, double relFrac,
                                      char **errStr) {
    if (spanFrac < 0.0 || spanFrac > 1.0 || relFrac < 0.0 || relFrac > 1.0) {
        set_err(errStr, "taffyBlockViz: min-block filter fractions must be in [0.0, 1.0]");
        return -1;
    }
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return -1;
    std::lock_guard<std::mutex> lock(H->mu);
    H->min_block_span_frac = spanFrac;
    H->min_block_rel_frac  = relFrac;
    return 0;
}

extern "C" int taffyGetMinBlockFilter(int h, double *spanFrac, double *relFrac,
                                      char **errStr) {
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return -1;
    std::lock_guard<std::mutex> lock(H->mu);
    if (spanFrac) *spanFrac = H->min_block_span_frac;
    if (relFrac)  *relFrac  = H->min_block_rel_frac;
    return 0;
}

// Run floor (bp): drop input runs shorter than this before chaining.
// -1 (default) = auto per query span; 0 = off; > 0 = explicit floor.
extern "C" int taffySetMinRunSize(int h, int64_t minRun, char **errStr) {
    if (minRun < -1) { set_err(errStr, "taffyBlockViz: minRun must be >= -1 (-1 = auto)"); return -1; }
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return -1;
    std::lock_guard<std::mutex> lock(H->mu);
    H->min_run_size = minRun;
    return 0;
}

extern "C" int taffyGetMinRunSize(int h, int64_t *minRun, char **errStr) {
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return -1;
    std::lock_guard<std::mutex> lock(H->mu);
    if (minRun) *minRun = H->min_run_size;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Free helpers                                                        */
/* ------------------------------------------------------------------ */

extern "C" void taffyFreeBlocks(struct taffy_block_t *b) {
    while (b) {
        struct taffy_block_t *nx = b->next;
        free(b->qChrom);
        free(b->qSequence);
        free(b->tSequence);
        free(b);
        b = nx;
    }
}

extern "C" void taffyFreeTargetDupeLists(struct taffy_target_dupe_list_t *d) {
    while (d) {
        struct taffy_target_dupe_list_t *nx = d->next;
        struct taffy_target_range_t *r = d->tRange;
        while (r) { struct taffy_target_range_t *rn = r->next; free(r); r = rn; }
        free(d->qChrom);
        free(d);
        d = nx;
    }
}

static void free_chain_summaries(struct taffy_chain_summary_t *c) {
    while (c) {
        struct taffy_chain_summary_t *nx = c->next;
        free(c);
        c = nx;
    }
}

extern "C" void taffyFreeBlockResults(struct taffy_block_results_t *res) {
    if (!res) return;
    taffyFreeBlocks(res->mappedBlocks);
    taffyFreeTargetDupeLists(res->targetDupeBlocks);
    free_chain_summaries(res->chainSummaries);
    free(res);
}

extern "C" void taffyFreeSpeciesList(struct taffy_species_t *s) {
    while (s) {
        struct taffy_species_t *nx = s->next;
        free(s->name);
        free(s->parentName);
        free(s);
        s = nx;
    }
}

extern "C" void taffyFreeChromList(struct taffy_chromosome_t *c) {
    while (c) {
        struct taffy_chromosome_t *nx = c->next;
        free(c->name);
        free(c);
        c = nx;
    }
}

extern "C" void taffyFreeMetadataList(struct taffy_metadata_t *m) {
    while (m) {
        struct taffy_metadata_t *nx = m->next;
        free(m->key);
        free(m->value);
        free(m);
        m = nx;
    }
}

/* ------------------------------------------------------------------ */
/* Species / Chroms                                                    */
/* ------------------------------------------------------------------ */

// Per-genome row.  Source of truth is the .tui's g-record roster
// (tui_genome_names); every .tui carries it, so the full genome name
// -- including any dotted version suffix (NCBI accessions like
// "GCA_028858775.2") -- comes straight from the roster.
struct GenomeRow { std::string genome; taffy_int_t total_bp = 0; taffy_int_t n_chroms = 0; };

static std::map<std::string, GenomeRow> enumerate_genomes(TaffyHandle *H) {
    std::map<std::string, GenomeRow> out;

    int64_t n_g = 0;
    TuiGenomeInfo *roster = tui_genome_names(H->tui, &n_g);
    if (roster == nullptr || n_g <= 0) return out;
    for (int64_t i = 0; i < n_g; i++) {
        GenomeRow row;
        row.genome   = roster[i].name;
        row.total_bp = roster[i].total_bp;
        row.n_chroms = roster[i].n_chroms;
        out[row.genome] = row;
    }
    tui_genome_info_free(roster, n_g);
    return out;
}

extern "C" struct taffy_species_t *taffyGetSpecies(int h, char **errStr) {
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return nullptr;
    std::lock_guard<std::mutex> lock(H->mu);
    auto rows = enumerate_genomes(H);

    struct taffy_species_t *head = nullptr, *tail = nullptr;
    for (auto &kv : rows) {
        struct taffy_species_t *s = (struct taffy_species_t *) st_calloc(1, sizeof(*s));
        s->name = strdup(kv.second.genome.c_str());
        s->length = kv.second.total_bp;
        s->numChroms = kv.second.n_chroms;
        s->parentName = strdup("");      // tree-walk not wired in initial cut
        s->parentBranchLength = 0.0;
        if (!head) head = s;
        if (tail) tail->next = s;
        tail = s;
    }
    return head;
}

extern "C" struct taffy_chromosome_t *taffyGetChroms(int h, const char *species, char **errStr) {
    if (!species) { set_err(errStr, "taffyBlockViz: NULL speciesName"); return nullptr; }
    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return nullptr;
    std::lock_guard<std::mutex> lock(H->mu);

    // Only this species' sequences, via a lower-bound seek + contiguous prefix
    // scan of the directory (O(log n_d + k)) instead of enumerating every
    // genome's sequences (the whole-directory slurp that dominated query time
    // on the 1155-genome .tui).  Keys are still full "<species>.<seq>" names,
    // so the prefix filter below is a cheap pass-through.
    stHash *seqs = tui_genome_seq_lengths(H->tui, species);
    if (!seqs) return nullptr;

    // Filter to entries whose key starts with "<species>." -- collect
    // into a vector + sort for stable output.
    std::string prefix = std::string(species) + ".";
    std::vector<std::pair<std::string, taffy_int_t>> hits;
    stHashIterator *it = stHash_getIterator(seqs);
    char *k;
    while ((k = (char *) stHash_getNext(it)) != NULL) {
        if (strncmp(k, prefix.c_str(), prefix.size()) != 0) continue;
        int64_t len = (int64_t)(intptr_t) stHash_search(seqs, k);
        hits.emplace_back(k + prefix.size(), len);
    }
    stHash_destructIterator(it);
    stHash_destruct(seqs);

    std::sort(hits.begin(), hits.end());

    struct taffy_chromosome_t *head = nullptr, *tail = nullptr;
    for (auto &h2 : hits) {
        struct taffy_chromosome_t *c = (struct taffy_chromosome_t *) st_calloc(1, sizeof(*c));
        c->name = strdup(h2.first.c_str());
        c->length = h2.second;
        if (!head) head = c;
        if (tail) tail->next = c;
        tail = c;
    }
    return head;
}

/* ------------------------------------------------------------------ */
/* mapBackAdjacencies helpers                                          */
/* ------------------------------------------------------------------ */

// HAL caps its mapAdjacencies inner loop at ~10 steps per side per
// mapped segment.  Match that here so the off-screen scan doesn't run
// away on regions with many short paralogous fragments.
static const int MAX_ADJ_SCAN = 10;

// Lazily load + cache the qSpecies sequence's runs table.  Key is
// "<qSpecies>.<qChrom>" (the .tui's d-line name).  Caller holds H->mu.
static const std::vector<SeqRun> *get_seq_runs(TaffyHandle *H,
                                               const std::string &full_seq_name) {
    auto it = H->qseq_runs_cache.find(full_seq_name);
    if (it != H->qseq_runs_cache.end()) return &it->second;

    int64_t n = 0;
    int64_t *flat = tui_load_seq_runs(H->tui, full_seq_name.c_str(), &n);
    std::vector<SeqRun> out;
    if (flat && n > 0) {
        out.reserve((size_t) n);
        for (int64_t i = 0; i < n; i++) {
            // The .tui stores lenc = (length << 1) | (1 - strand).
            // Decoded by tui_load_seq_runs into the third field of each
            // triple as-is; reverse-strand bit is in the LSB.
            int64_t lenc = flat[3*i + 2];
            SeqRun r;
            r.t_start = flat[3*i + 0];
            r.g_start = flat[3*i + 1];
            r.length  = lenc >> 1;
            r.strand  = (lenc & 1) ? -1 : +1;
            out.push_back(r);
        }
    }
    free(flat);
    auto ins = H->qseq_runs_cache.emplace(full_seq_name, std::move(out));
    return &ins.first->second;
}

// Build a taffy_block_t for an off-screen flank candidate.
// Forward-strand half-open coords; strand reflects the candidate run's
// strand (the alignment direction, NOT post-XOR with input strand).
static struct taffy_block_t *make_flank_block(const SeqRun &cand,
                                              const char *qChrom,
                                              int64_t tStart_fwd) {
    struct taffy_block_t *b = (struct taffy_block_t *) st_calloc(1, sizeof(*b));
    b->qChrom    = strdup(qChrom);
    b->tStart    = tStart_fwd;
    b->qStart    = cand.t_start;
    b->size      = cand.length;
    b->strand    = (cand.strand > 0) ? '+' : '-';
    b->qSequence = nullptr;
    b->tSequence = nullptr;
    return b;
}

// Scan up to MAX_ADJ_SCAN qSpecies-coord neighbors of the chain edge,
// back-map each candidate to tSpecies via gl_t, and return the first
// off-screen flank (i.e. whose tSpecies range falls outside [tStartUser,
// tEndUser) on tChrom_filter).  dir = -1 walks left from `start_idx`,
// dir = +1 walks right.  Returns NULL if no off-screen flank is found
// within the scan budget.
static struct taffy_block_t *find_off_screen_flank(
        const std::vector<SeqRun> *runs, int64_t start_idx, int dir,
        const TuiGenomeLift *gl_t,
        const char *qChrom,
        const char *tChrom_filter,
        int64_t tStartUser, int64_t tEndUser) {
    const int64_t n = (int64_t) runs->size();
    for (int step = 1; step <= MAX_ADJ_SCAN; step++) {
        int64_t idx = start_idx + dir * step;
        if (idx < 0 || idx >= n) break;
        const SeqRun &cand = (*runs)[idx];
        // Back-map BOTH ends of the candidate to tSpecies via gl_t.  If
        // either end's seq differs from tChrom_filter, the run lifts to
        // a different tChrom or is unmapped on tSpecies -- not a usable
        // flank on the same browser track.
        TuiGenomeMatch m_lo[1], m_hi[1];
        int nlo = tui_genome_lift_column((TuiGenomeLift *) gl_t,
                                         cand.g_start, m_lo, 1);
        int nhi = tui_genome_lift_column((TuiGenomeLift *) gl_t,
                                         cand.g_start + cand.length - 1, m_hi, 1);
        if (nlo == 0 || nhi == 0) continue;
        if (strcmp(m_lo[0].seq, tChrom_filter) != 0) continue;
        if (strcmp(m_hi[0].seq, tChrom_filter) != 0) continue;
        int64_t tLo = (m_lo[0].pos < m_hi[0].pos) ? m_lo[0].pos : m_hi[0].pos;
        int64_t tHi = ((m_lo[0].pos > m_hi[0].pos) ? m_lo[0].pos : m_hi[0].pos) + 1;
        // Off-screen test: the candidate's tSpecies range must NOT
        // overlap the user's [tStartUser, tEndUser).
        if (tLo < tEndUser && tHi > tStartUser) continue;   // on-screen
        return make_flank_block(cand, qChrom, tLo);
    }
    return nullptr;
}

/* ------------------------------------------------------------------ */
/* Block query (the hot path)                                          */
/* ------------------------------------------------------------------ */

// One buffered, clipped qSpecies run awaiting the windowed chain pass.
// POD, appended by the visitor with NO per-run heap traffic -- the old
// per-run st_calloc + strdup + std::vector<ptr> + std::set::emplace was
// the chromosome-scale bottleneck (100x slower than `taffy lift`, whose
// architecture this mirrors).  Same shape as taf_lift.c's LiftRun, plus
// the (t_ref, q_qsp) block geometry blockViz emits.  `seq` is borrowed
// from the qSpecies TuiGenomeLift -- pointer stable until that gl is
// destructed (tui.h TuiRun) -- so no interning is needed: taffy_chain
// compares t_name by strcmp and the gl hands out one stable pointer per
// chrom.
struct BlockRun {
    const char *seq;            // qChrom (gl-borrowed, stable)
    int64_t     c_start, c_end; // universal-column range (window decision)
    int64_t     t_ref;          // tSpecies block coord (block.tStart)
    int64_t     q_qsp;          // qSpecies forward coord (block.qStart)
    int64_t     size;           // bp
    int         strand_sign;    // +1 / -1 relative strand (-> TaffyAln.strand)
};

// One survivor of the windowed chain+merge pass: a collapsed, gap-free
// 2-D block.  Accumulated flat across all windows/ivs; a single GLOBAL
// chain pass over these (cheap -- they're already merged, so there are
// few) does the paralogy filter, chain-id assignment, primary pick,
// dupMode routing, and budget.  Filtering CANNOT be done per window: a
// browser query's paralogy lives on the tSpecies axis, which spans ivs
// (the same tSpecies bp reached via two universal-column regions is a
// dupe), so the overlap-frac select must see every iv at once -- unlike
// `taffy lift`, whose chain q-axis is the column and is disjoint across
// ivs by construction.  qChrom is the gl-borrowed pointer (strdup'd
// only when the final taffy_block_t is built).
struct SurvBlock {
    const char *qChrom;
    int64_t     tStart, qStart, size;
    char        strand;
};

struct BlockCtx {
    // ---- per-iv source-axis state (set before each visit) -------------
    int64_t     c_lo = 0, c_hi = 0;     // current iv's column clip range
    int64_t     tpos_at_c_lo = 0;       // tSpecies pos at column c_lo (iv.t_start)
    int         iv_rev = 0;             // iv.rev: 0 col-asc->t-asc, 1 col-asc->t-desc.
                                        // XOR'd with each run's strand for the true
                                        // tSpecies<->qSpecies relative strand.
    const char *qChromFilter = nullptr; // optional qChrom restriction
    std::string tFullName;              // "<tSpecies>.<tChrom>" -- the chain q_name

    // ---- chain tuning (copied from the handle) ------------------------
    int64_t     chain_open = TAFFY_CHAIN_DEFAULT_OPEN;
    int64_t     chain_extend = TAFFY_CHAIN_DEFAULT_EXTEND;
    int64_t     chain_max_gap = TAFFY_CHAIN_DEFAULT_MAX_GAP;
    double      chain_overlap_frac = 0.0;   // <0 disables the per-window filter
    int64_t     min_run = 0;                // per-query run floor (bp); set from span in get_blocks_impl

    // ---- windowed chain processing (mirrors taf_lift.c) ---------------
    int64_t     chain_window = 0;        // W (column bp); 0 = no within-iv windowing
    int64_t     chain_window_overlap = 0;// K (column bp) carried across boundaries
    int64_t     window_start = 0;        // current window left column edge (INT64_MIN=unset)

    // ---- flat run buffer (lift-style; no per-run heap) ----------------
    BlockRun   *run_buf = nullptr;
    size_t      n_run_buf = 0, cap_run_buf = 0;

    // ---- survivor accumulator -----------------------------------------
    std::vector<SurvBlock> survivors;
};

// blockViz windowing constants, mirroring taf_lift.c's W/K defaults.  A
// browser query is one viewport (<= a chromosome), so within-iv
// windowing only engages on multi-hundred-Mb ranges; there it bounds
// the chain pass + buffer memory.  K >= chain_max_gap so a chain
// crossing a boundary re-forms in the carried tail.
static const int64_t BLOCK_CHAIN_WINDOW         = 100LL * 1000 * 1000; // W
static const int64_t BLOCK_CHAIN_WINDOW_OVERLAP =  10LL * 1000 * 1000; // K

// Max-gap bound for the GLOBAL (level-2) greedy regroup of survivors.
// Level 1 already chained each window's runs; the global greedy pass
// only needs to BRIDGE the reference-indel gaps that tui_query split a
// diagonal across.  Two consecutive survivors join one chain only if
// their gap on each axis is within this bound -- generous enough to span
// real indels, tight enough that an off-diagonal paralog (whose qSpecies
// coord jumps far) starts a fresh chain.  BLOCK_L2_GROUP_SLACK is the
// small reverse tolerance allowed for boundary round-off / tiny overlaps.
static const int64_t BLOCK_L2_MAX_GAP      = 1000 * 1000;   // 1 Mb
// Max disagreement between the reference- and qSpecies-axis gaps for two
// consecutive blocks to count as collinear (same diagonal) and join one
// chain.  The g1000 index emits redundant nested/overlapping fragments
// on the SAME diagonal (a small block inside a big one -> NEGATIVE gaps
// that nonetheless agree, t_gap == q_gap); grouping on proximity alone
// breaks the chain at each, shattering one orthologous arm into dozens.
// Collinearity joins them while a true off-diagonal paralog (gaps
// disagree) or the centromere step (huge one-axis jump) still splits.
// 50 kb (not 10 kb): a real inverted/orthologous arm carries internal
// indels up to tens of kb -- e.g. the 8p23.1 inversion shatters into 5
// reverse chains at 10 kb but reunites into ONE at >= 25 kb.  Validated
// not to over-merge paralogs (they sit Mb+ off-diagonal): chr9:62-63M SD
// chain count is unchanged 10kb->100kb and chr9:72-73M is identical.
static const int64_t BLOCK_L2_COLLINEAR_TOL = 50000;        // bp

// Coarsen (bin) when the per-block budget would make each output block
// span at least this many reference bp -- i.e. when (queried span /
// max_output_blocks) >= this.  This makes the detail<->coverage switch a
// pure function of the ZOOM (query span), NOT of how many blocks happen
// to fall in the window.  A count-based switch flips as you pan a fixed
// zoom across a dense SD locus -- the window total crosses the cap, the
// whole result re-renders coverage-vs-detail, and dupe bars blink.  A
// span-based switch is panning-stable: the same locus renders the same
// way at the same zoom regardless of the window's edges.
static const int64_t BLOCK_MIN_COARSE_BIN = 10000;         // bp

// taffy_chain_merge_collinear callback.  chain.c has already grown the
// KEPT aln's q/t extents to the union; here we additionally (1) grow the
// kept run's COLUMN extent so the window "decided vs live" test sees the
// merged span, and (2) mark the absorbed run dead (user=NULL) so the
// flush's decide loop skips it (and it isn't carried forward).
static void bv_on_merge(TaffyAln *kept, TaffyAln *absorbed, void *user) {
    (void) user;
    BlockRun *K = (BlockRun *) kept->user;
    BlockRun *A = (BlockRun *) absorbed->user;
    if (K && A) {
        if (A->c_start < K->c_start) K->c_start = A->c_start;
        if (A->c_end   > K->c_end)   K->c_end   = A->c_end;
    }
    absorbed->user = nullptr;
}

/* Chain + paralogy-filter + collinear-merge one window's worth of
 * buffered runs, emit the decided (kept, merged) survivors, and carry
 * live runs forward.  This is the SPEED tier and is structurally
 * taf_lift.c::chain_flush_window: the per-window overlap-frac filter is
 * what keeps the chain pass small in paralog-dense (SD / pericentromeric)
 * regions -- without it the survivor set explodes and the global pass in
 * get_blocks_impl blows up (125x slower, measured).  A residual GLOBAL
 * filter still runs there to catch CROSS-iv dupes the per-window pass
 * can't see (same tSpecies bp via two column regions); see SurvBlock.
 *
 * force_all != 0 (end of iv): every buffered run is decided and the
 * buffer cleared.  Otherwise runs whose (possibly merged) column extent
 * ends at/below boundary_c are decided; the rest stay live for the next
 * window's pass so cross-boundary chains re-form.  Live runs always
 * carry (even if currently filtered out) -- their chain may re-form with
 * more context; only DECIDED runs consult the keep set.  Survivors
 * become 2-D SurvBlock records (via taffy_chain_merge_collinear) rather
 * than being replayed through lift's 1-D pending_push BED merge. */
static void block_flush_window(BlockCtx *cx, int64_t boundary_c, int force_all) {
    int64_t n = (int64_t) cx->n_run_buf;
    if (n == 0) return;

    TaffyAln *alns = (TaffyAln *) st_malloc((size_t) n * sizeof(TaffyAln));
    for (int64_t i = 0; i < n; i++) {
        BlockRun *L = &cx->run_buf[i];
        alns[i].q_name  = cx->tFullName.c_str();   // single tSpecies.tChrom partition
        alns[i].q_start = L->t_ref;
        alns[i].q_end   = L->t_ref + L->size;
        alns[i].t_name  = L->seq;
        alns[i].t_start = L->q_qsp;
        alns[i].t_end   = L->q_qsp + L->size;
        alns[i].strand  = L->strand_sign;
        alns[i].score   = L->size;
        alns[i].user    = (void *) L;              // run_buf is stable during a flush
    }

    TaffyChainCostParams cost = { cx->chain_open, cx->chain_extend };
    int64_t *chain_id = (int64_t *) st_calloc((size_t) n, sizeof(int64_t));
    TaffyChainInfo *chains_out = nullptr;
    int64_t n_chains_out = 0;
    taffy_chain(alns, n, taffy_chain_default_gap_cost, &cost,
                cx->chain_max_gap, chain_id, &chains_out, &n_chains_out);

    int64_t max_id = 0;
    for (int64_t k = 0; k < n_chains_out; k++)
        if (chains_out[k].id > max_id) max_id = chains_out[k].id;

    // Per-window overlap-frac paralogy filter (disabled => keep all).
    char *keep_chain = (char *) st_calloc((size_t)(max_id + 1), sizeof(char));
    if (cx->chain_overlap_frac >= 0) {
        taffy_chain_overlap_frac_select(alns, n, chain_id, chains_out,
                                        n_chains_out, max_id,
                                        cx->chain_overlap_frac, /*cap=*/0,
                                        keep_chain);
    } else {
        for (int64_t k = 0; k <= max_id; k++) keep_chain[k] = 1;
    }

    // Collapse adjacent collinear alns within a chain into one block.
    // bv_on_merge grows the kept run's column extent + kills absorbed.
    taffy_chain_merge_collinear(alns, n, chain_id, bv_on_merge, cx);

    // Decide each surviving (non-absorbed) run: emit it if past the
    // window boundary AND its chain is kept, else carry it live.
    // Geometry is read from the merged alns[i]; live runs sync it back
    // into their BlockRun so the next window re-chains the merged extent.
    char *is_live = (char *) st_calloc((size_t) n, sizeof(char));
    for (int64_t i = 0; i < n; i++) {
        BlockRun *L = (BlockRun *) alns[i].user;
        if (L == nullptr) continue;                 // absorbed by merge
        int decided = force_all || L->c_end <= boundary_c;
        if (decided) {
            if (keep_chain[chain_id[i]]) {
                SurvBlock s;
                s.qChrom = alns[i].t_name;
                s.tStart = alns[i].q_start;
                s.size   = alns[i].q_end - alns[i].q_start;
                s.qStart = alns[i].t_start;         // min qSpecies (both strands post-merge)
                s.strand = (alns[i].strand > 0) ? '+' : '-';
                cx->survivors.push_back(s);
            }
            // filtered-out decided runs simply vanish
        } else {
            L->t_ref = alns[i].q_start;
            L->size  = alns[i].q_end - alns[i].q_start;
            L->q_qsp = alns[i].t_start;
            is_live[L - cx->run_buf] = 1;
        }
    }

    // Compact live runs to the head of run_buf for the next window.
    int64_t w = 0;
    for (int64_t bi = 0; bi < n; bi++) {
        if (!is_live[bi]) continue;
        if (w != bi) cx->run_buf[w] = cx->run_buf[bi];
        w++;
    }
    cx->n_run_buf = (size_t) w;

    free(is_live);
    free(keep_chain);
    free(chain_id);
    free(chains_out);
    free(alns);
}

// Visitor: clip one qSpecies run to the current iv's column window,
// compute its block geometry + relative strand (the iv_rev XOR run
// strand fix), and append a flat BlockRun -- NO heap traffic.  When the
// buffered column span crosses the chain window, flush + slide.
static void block_visit_cb(const TuiRun *r, void *user) {
    BlockCtx *cx = (BlockCtx *) user;
    if (cx->qChromFilter && strcmp(r->seq, cx->qChromFilter) != 0) return;
    int64_t r_end = r->g_start + r->length;
    int64_t cs = r->g_start > cx->c_lo ? r->g_start : cx->c_lo;
    int64_t ce = r_end < cx->c_hi      ? r_end      : cx->c_hi;
    if (cs >= ce) return;
    int64_t size = ce - cs;

    // tSpecies position at the clipped run [cs, ce).  iv_rev=0: column-asc
    // matches tSpecies-asc so the low end is at column cs.  iv_rev=1:
    // column-asc is tSpecies-desc so the low end is at column ce-1.
    int64_t tStart = cx->iv_rev
        ? cx->tpos_at_c_lo - (ce - 1 - cx->c_lo)
        : cx->tpos_at_c_lo + (cs - cx->c_lo);

    // Relative strand tSpecies<->qSpecies = iv_rev XOR run strand.  The
    // XOR is essential: an ancestor block reverse-mapped to BOTH genomes
    // is biologically forward, and without it shows a spurious '-'.
    // r->strand is 1 for forward, 0 for reverse.
    int actual_fwd = cx->iv_rev ? !r->strand : r->strand;

    // qSpecies forward coord at the run's low tSpecies end.  Uses
    // r->strand (column->q direction), independent of the iv_rev factor
    // that flips the block's strand label.
    int64_t qStart = r->strand
        ? r->t_start + (cs - r->g_start)
        : r->t_start + r->length - (ce - r->g_start);

    if (size < cx->min_run) return;   // sub-floor run: skip before it enters the chain buffer
    if (cx->n_run_buf == cx->cap_run_buf) {
        cx->cap_run_buf = cx->cap_run_buf ? cx->cap_run_buf * 2 : 4096;
        cx->run_buf = (BlockRun *) st_realloc(cx->run_buf,
                                              cx->cap_run_buf * sizeof(*cx->run_buf));
    }
    BlockRun *L    = &cx->run_buf[cx->n_run_buf++];
    L->seq         = r->seq;
    L->c_start     = cs;
    L->c_end       = ce;
    L->t_ref       = tStart;
    L->q_qsp       = qStart;
    L->size        = size;
    L->strand_sign = actual_fwd ? +1 : -1;

    // Windowed flush: anchor the window on this iv's first run, then
    // flush + slide whenever the buffered column span reaches W.
    if (cx->window_start == INT64_MIN) cx->window_start = cs;
    if (cx->chain_window > 0 && ce - cx->window_start >= cx->chain_window) {
        int64_t window_end = cx->window_start + cx->chain_window;
        int64_t boundary_c = window_end - cx->chain_window_overlap;
        block_flush_window(cx, boundary_c, /*force_all=*/0);
        cx->window_start = boundary_c;
    }
}

static struct taffy_block_results_t *
get_blocks_impl(int h, const char *qSpecies, const char *tSpecies,
                const char *tChrom, taffy_int_t tStart, taffy_int_t tEnd,
                taffy_int_t tReversed, taffy_seqmode_type_t seqMode,
                taffy_dup_type_t dupMode, int mapBackAdjacencies,
                const char *qChrom, const char *coalescenceLimitName,
                char **errStr) {
    // Reject unsupported parameters loudly so the caller knows.
    if (tReversed != 0) { set_err(errStr, "taffyBlockViz: tReversed not supported"); return nullptr; }
    // mapBackAdjacencies: implemented below; honour the flag.
    if (coalescenceLimitName) { set_err(errStr, "taffyBlockViz: coalescenceLimitName not supported"); return nullptr; }
    if (seqMode != TAFFY_NO_SEQUENCES) { set_err(errStr, "taffyBlockViz: seqMode != NO_SEQUENCES not supported"); return nullptr; }
    if (dupMode != TAFFY_NO_DUPS && dupMode != TAFFY_QUERY_DUPS && dupMode != TAFFY_QUERY_AND_TARGET_DUPS) {
        set_err(errStr, "taffyBlockViz: unknown dupMode");
        return nullptr;
    }
    if (!qSpecies || !tSpecies || !tChrom) { set_err(errStr, "taffyBlockViz: missing required arg"); return nullptr; }
    if (tEnd < 0 || (tEnd > 0 && tEnd <= tStart)) { set_err(errStr, "taffyBlockViz: bad tStart/tEnd"); return nullptr; }

    TaffyHandle *H = lookup_handle(h, errStr);
    if (!H) return nullptr;
    std::lock_guard<std::mutex> lock(H->mu);

    // Resolve the .tui d-line key for tChrom: "<tSpecies>.<tChrom>".
    std::string tFullName = std::string(tSpecies) + "." + tChrom;

    // tEnd == 0 -> use full chrom length.  Resolve it with a single directory
    // binary search (O(log n_d)) instead of building a hash of every genome's
    // sequence table -- the whole-directory slurp cost ~2.4 s / 794 MB per
    // query on the 1155-genome .tui, just to read one length.
    if (tEnd == 0) {
        int64_t fullLen = tui_seq_length(H->tui, tFullName.c_str());
        if (fullLen > 0) tEnd = fullLen;
        if (tEnd == 0) { set_err(errStr, "taffyBlockViz: unknown tChrom"); return nullptr; }
    }

    // Get the qSpecies lift table (cached).
    TuiGenomeLift *gl_q = get_or_load_gl(H, qSpecies, errStr);
    if (!gl_q) return nullptr;

    // tSpecies range -> column intervals via tui_query.  The intervals
    // are sorted + merged in column order; the cumulative tSpecies
    // bp within them equals the input range covered by alignment.
    int64_t n_iv = 0;
    TuiInterval *iv = tui_query(H->tui, tFullName.c_str(), tStart, tEnd, &n_iv);
    if (!iv || n_iv == 0) {
        struct taffy_block_results_t *res = (struct taffy_block_results_t *) st_calloc(1, sizeof(*res));
        free(iv);
        return res;   // empty result, not an error
    }

    // Visit each iv, buffering clipped runs into a flat POD array and
    // chaining + merging them in column windows (block_flush_window).
    // This is `taffy lift`'s architecture -- the source of its
    // chromosome-scale speed -- and replaces the old per-run
    // taffy_block_t + std::vector<ptr> + std::set bottleneck.
    BlockCtx cx;
    cx.qChromFilter         = qChrom;
    cx.tFullName            = tFullName;          // single chain q_name partition
    cx.chain_open           = H->chain_open;
    cx.chain_extend         = H->chain_extend;
    cx.chain_max_gap        = H->max_gap_length;
    cx.chain_overlap_frac   = H->chain_overlap_frac;
    cx.chain_window         = BLOCK_CHAIN_WINDOW;
    cx.chain_window_overlap = BLOCK_CHAIN_WINDOW_OVERLAP;
    // Adaptive run floor: on a wide overview, drop sub-pixel runs before
    // they reach the chain DP.  A whole-chromosome view is ~90% sub-100bp
    // runs (<1% of bp), and the O(n log n) sweep chokes on them in segdup-
    // dense regions.  Auto = span/500000 (deeply sub-pixel: ~100-500 bp at
    // whole-chrom, 0 below ~500 kb so detail views stay exact).  A handle
    // override (>= 0) wins.
    cx.min_run = (H->min_run_size >= 0) ? H->min_run_size
                                        : (int64_t)(tEnd - tStart) / 500000;

    // Per-iv loop: tui_query hands us iv.t_start (tSpecies pos at
    // iv.start) and iv.rev (source-to-column orientation) per chunk.
    // Each iv is a disjoint, internally-monotone universal-column range,
    // so the chain window anchors on the iv's first run and the buffer
    // is force-drained at the iv boundary (chains never span ivs --
    // they're column-disjoint, far past chain_max_gap).
    for (int64_t k = 0; k < n_iv; k++) {
        cx.c_lo         = iv[k].start;
        cx.c_hi         = iv[k].end;
        cx.tpos_at_c_lo = iv[k].t_start;
        cx.iv_rev       = iv[k].rev;
        cx.window_start = INT64_MIN;
        tui_genome_lift_visit_runs(gl_q, iv[k].start, iv[k].end,
                                   block_visit_cb, &cx);
        if (cx.n_run_buf > 0)
            block_flush_window(&cx, /*boundary_c=*/0, /*force_all=*/1);
    }
    free(iv);
    free(cx.run_buf);

    struct taffy_block_results_t *res = (struct taffy_block_results_t *) st_calloc(1, sizeof(*res));
    if (cx.survivors.empty()) return res;  // nothing aligned in range

    // Materialize one taffy_block_t per survivor, each carried as the
    // TaffyAln.user back-pointer for the global grouping pass below.  The
    // routing loop later frees each block or hands it to the output list.
    // The old code's mistake was this per-block alloc per RAW RUN; here
    // it's once per merged survivor.  t_name is the gl-borrowed qChrom
    // (stable through the grouping + overlap-frac calls).
    std::vector<TaffyAln> alns;
    alns.reserve(cx.survivors.size());
    for (const SurvBlock &s : cx.survivors) {
        struct taffy_block_t *b = (struct taffy_block_t *) st_calloc(1, sizeof(*b));
        b->qChrom = strdup(s.qChrom);
        b->tStart = s.tStart;
        b->qStart = s.qStart;
        b->size   = s.size;
        b->strand = s.strand;

        TaffyAln a = {0};
        a.q_name  = cx.tFullName.c_str();    // single tSpecies.tChrom partition
        a.q_start = s.tStart;
        a.q_end   = s.tStart + s.size;
        a.t_name  = s.qChrom;                // gl-borrowed, stable through the chain call
        a.t_start = s.qStart;
        a.t_end   = s.qStart + s.size;
        a.strand  = (s.strand == '+') ? +1 : -1;
        a.score   = s.size;
        a.user    = (void *) b;
        alns.push_back(a);
    }
    int64_t n = (int64_t) alns.size();

    // Group survivors into chains, then paralogy-filter + budget them.
    //
    // A genome-browser reference is usually a LEAF (e.g. hg38), whose
    // alignment to the universal columns is broken by every indel into
    // tens of thousands of tiny ivs -- so tui_query hands back a heavily
    // fragmented, paralog-laden survivor set that the per-window filter
    // (it sees one iv at a time) cannot dedup.  The full taffy_chain DP
    // would reunite the orthologous diagonal correctly, but its
    // O(n*active) sweep -- and the O(kept^2) of an overlap-frac pass run
    // over an ungrouped survivor set -- both go quadratic at chromosome
    // scale (measured ~20 s on 138-248 Mb chromosomes).  Instead a single
    // GREEDY collinear sweep (O(n log n)) collapses the diagonal into one
    // chain and peels off-diagonal paralogs into their own; overlap-frac
    // then runs over that SMALL chain set (cheap), and the budget caps
    // the output.  (lift stays fast on the same data only because it is
    // queried on the ANCESTOR, which maps contiguously -> few big ivs.)
    int64_t primary_id = 0;
    std::vector<TaffyChainInfo> chains;
    int64_t n_chains = 0;
    std::vector<int64_t> chain_id((size_t) n, 0);

    // Sort by (qChrom, strand, tStart): contiguous, tStart-ascending
    // partitions -- the order overlap-frac expects and the greedy needs.
    std::sort(alns.begin(), alns.end(),
              [](const TaffyAln &a, const TaffyAln &b) {
                  int c = strcmp(a.t_name ? a.t_name : "", b.t_name ? b.t_name : "");
                  if (c != 0) return c < 0;
                  if (a.strand  != b.strand)  return a.strand  < b.strand;
                  if (a.q_start != b.q_start) return a.q_start < b.q_start;
                  // TOTAL order: same (t_name,strand,q_start) alns must sort
                  // build-independently (taffy_chain's upstream sort breaks
                  // such ties by POINTER), else the greedy groups them by heap
                  // layout -> non-deterministic coverage across rebuilds.
                  if (a.q_end   != b.q_end)   return a.q_end   < b.q_end;
                  if (a.t_start != b.t_start) return a.t_start < b.t_start;
                  return a.t_end < b.t_end;
              });

    // Greedy collinear grouping within a (qChrom, strand) partition.
    // Extend the current chain with the next aln iff it's the same qChrom
    // + strand and lies on the SAME diagonal: the reference-axis and
    // qSpecies-axis gaps between the two blocks agree within
    // BLOCK_L2_COLLINEAR_TOL, and neither axis jumps more than
    // BLOCK_L2_MAX_GAP.  Both gaps are SIGNED -- the g1000 index's nested
    // fragments give negative (overlapping) gaps that still agree, so
    // collinearity (not mere proximity) is what keeps one orthologous arm
    // a single chain instead of shattering it at every nested fragment.
    // An off-diagonal paralog (gaps disagree) or the centromere step (a
    // huge q-jump) opens a new chain.
    for (int64_t i = 0; i < n; i++) {
        bool extend = false;
        if (i > 0) {
            const TaffyAln &A = alns[i - 1];
            const TaffyAln &B = alns[i];
            if (A.strand == B.strand &&
                strcmp(A.t_name ? A.t_name : "", B.t_name ? B.t_name : "") == 0) {
                int64_t t_gap = B.q_start - A.q_end;          // reference axis
                int64_t q_gap = (A.strand > 0) ? (B.t_start - A.t_end)
                                               : (A.t_start - B.t_end);
                int64_t at = t_gap < 0 ? -t_gap : t_gap;
                int64_t aq = q_gap < 0 ? -q_gap : q_gap;
                int64_t dd = (t_gap - q_gap) < 0 ? (q_gap - t_gap) : (t_gap - q_gap);
                if (at <= BLOCK_L2_MAX_GAP && aq <= BLOCK_L2_MAX_GAP &&
                    dd <= BLOCK_L2_COLLINEAR_TOL)
                    extend = true;
            }
        }
        if (!extend) {
            TaffyChainInfo ci;
            ci.id          = ++n_chains;    // 1-based, monotonic
            ci.total_score = 0;
            ci.total_bp    = 0;
            ci.n_alns      = 0;
            chains.push_back(ci);
        }
        TaffyChainInfo &cur = chains.back();
        chain_id[i]      = cur.id;
        cur.total_score += alns[i].score;
        cur.total_bp    += (alns[i].q_end - alns[i].q_start);
        cur.n_alns      += 1;
    }

    // Rank chains by score desc (primary first), matching taffy_chain's
    // output contract so the budget + summaries logic stays unchanged.
    std::sort(chains.begin(), chains.end(),
              [](const TaffyChainInfo &a, const TaffyChainInfo &b) {
                  if (a.total_score != b.total_score) return a.total_score > b.total_score;
                  return a.id < b.id;
              });
    primary_id = (n_chains > 0) ? chains[0].id : 0;

    // Overlap-frac paralogy filter on the tSpecies (q) axis over the
    // SMALL greedy chain set: drops chains whose q-coverage is already
    // held by a higher-scoring kept chain (default frac 0.5).  Cheap now
    // that the diagonal is grouped, not 10^5 disjoint pieces.  The output
    // cap is then applied POSITIONALLY in the emit step (below), not as a
    // by-score chain budget here.
    std::set<int64_t> kept_chains;
    std::vector<char> ovr_keep;
    if (n_chains > 0) {
        if (H->chain_overlap_frac >= 0) {
            int64_t max_id = 0;
            for (int64_t k = 0; k < n_chains; k++)
                if (chains[k].id > max_id) max_id = chains[k].id;
            ovr_keep.assign((size_t)(max_id + 1), 0);
            taffy_chain_overlap_frac_select(alns.data(), n,
                                            chain_id.data(),
                                            chains.data(), n_chains, max_id,
                                            H->chain_overlap_frac,
                                            /*cap=*/0,
                                            ovr_keep.data());
        }

        // Keep EVERY overlap-frac survivor chain -- NO chain-count
        // budget.  A by-score block budget collapses the output onto
        // whichever chain has the most bp, which on a chromosome with a
        // dense pericentromeric satellite (chr9 9q12, the acrocentrics,
        // chr16) is that satellite -- it swallows the whole cap and the
        // arms vanish.  Instead the cap is enforced POSITIONALLY in the
        // emit step, so the budget spreads across [tStart, tEnd).
        kept_chains.insert(primary_id);
        for (int64_t k = 0; k < n_chains; k++) {
            int64_t cid = chains[k].id;
            if (ovr_keep.empty() || ovr_keep[cid]) kept_chains.insert(cid);
        }
    }

    // Collect the kept candidate blocks (overlap-frac survivors, dupMode-
    // filtered) into a flat vector; free the rest.  Then either emit them
    // all at full detail (narrow query, count <= cap) or POSITIONALLY
    // COARSEN to the cap (wide query) so the budget spreads across
    // [tStart, tEnd) instead of piling onto the densest locus.
    struct taffy_block_t *mapped_head = nullptr, *mapped_tail = nullptr;
    int64_t mapped_count = 0;
    // Output ceiling.  Coverage (wide/binned) output is capped at
    // max_output_blocks (one per bin).  Detail (narrow) output emits every
    // block in the window -- the window is span-bounded, so the count is
    // bounded by the LOCAL alignment density, not the cap; a high backstop
    // only guards a pathologically dense locus.  Capping detail at
    // max_output_blocks instead would truncate (window-dependently) and
    // re-introduce the very instability the span-based switch removes.
    int64_t emit_cap = H->max_output_blocks;
    auto append_mapped = [&](struct taffy_block_t *b) {
        if (mapped_count >= emit_cap) {
            free(b->qChrom);
            free(b);
            return;
        }
        b->next = nullptr;
        if (!mapped_head) mapped_head = b;
        if (mapped_tail) mapped_tail->next = b;
        mapped_tail = b;
        mapped_count++;
    };
    // Drop sub-resolution noise from a to-be-emitted block vector: a block is
    // dropped iff size < min(window*span_frac, maxSize*rel_frac), where maxSize
    // is taken over the SAME vector (so the largest block always survives and a
    // uniformly-small region is preserved).  Frees the dropped blocks.  No-op
    // unless BOTH fractions are > 0.  See the handle fields for the rationale.
    auto filter_noise = [&](std::vector<struct taffy_block_t *> &v) {
        if (H->min_block_span_frac <= 0.0 || H->min_block_rel_frac <= 0.0 || v.empty())
            return;
        int64_t span  = (int64_t) tEnd - (int64_t) tStart;
        int64_t maxsz = 0;
        for (struct taffy_block_t *b : v) if (b->size > maxsz) maxsz = b->size;
        int64_t win_thr = (int64_t) ((double) span  * H->min_block_span_frac);
        int64_t rel_thr = (int64_t) ((double) maxsz * H->min_block_rel_frac);
        int64_t thr = win_thr < rel_thr ? win_thr : rel_thr;   // below BOTH == below min
        v.erase(std::remove_if(v.begin(), v.end(),
                    [&](struct taffy_block_t *b) {
                        if (b->size >= thr) return false;
                        free(b->qChrom); free(b);
                        return true;
                    }), v.end());
    };
    // For QUERY_AND_TARGET_DUPS: collect non-primary alns grouped by
    // chain_id so each non-primary chain becomes one dupe-list entry.
    std::map<int64_t, taffy_target_dupe_list_t *> dupe_by_chain;
    // Per-qChrom chain summaries, populated only on the coarse path.
    std::vector<TaffyChainInfo> coarse_summaries;

    std::vector<struct taffy_block_t *> cand;
    cand.reserve((size_t) n);
    for (int64_t i = 0; i < n; i++) {
        struct taffy_block_t *b = (struct taffy_block_t *) alns[i].user;
        if (b == nullptr) continue;   // defensive (no L2 merge nulls user now)
        int64_t cid = chain_id[i];
        if (!kept_chains.count(cid) ||
            (dupMode == TAFFY_NO_DUPS && cid != primary_id)) {
            free(b->qChrom);
            free(b);
            continue;
        }
        b->chainId = cid;   // browser groups blocks of a chain for snake-trace
        cand.push_back(b);
    }

    // Within each chain, fold redundant nested/overlapping blocks into
    // one tile.  The g1000 index emits the SAME diagonal at two
    // granularities -- a merged big block plus the fine ungapped fragments
    // nested inside it -- and a snake renderer then loops from the big
    // block's edge back to each nested sliver (a jarring fake
    // "rearrangement").  Two blocks are merged ONLY when they overlap on
    // BOTH the reference AND the query axis (i.e. they're the same diagonal
    // region, redundant).  Blocks that merely abut on the reference axis
    // but step on the query axis are DIFFERENT alignments (a real indel /
    // dup boundary) and stay separate, as do real gaps -- so genuine fine
    // structure survives; only redundant nesting collapses.
    std::sort(cand.begin(), cand.end(),
              [](const struct taffy_block_t *a, const struct taffy_block_t *b) {
                  if (a->chainId != b->chainId) return a->chainId < b->chainId;
                  if (a->tStart  != b->tStart)  return a->tStart  < b->tStart;
                  if (a->size    != b->size)    return a->size    < b->size;
                  return a->qStart < b->qStart;   // TOTAL order: build-stable
              });
    {
        std::vector<struct taffy_block_t *> merged;
        merged.reserve(cand.size());
        for (struct taffy_block_t *b : cand) {
            struct taffy_block_t *last = merged.empty() ? nullptr : merged.back();
            bool do_merge = false;
            if (last && last->chainId == b->chainId) {
                int64_t l_te = last->tStart + last->size, b_te = b->tStart + b->size;
                int64_t l_qe = last->qStart + last->size, b_qe = b->qStart + b->size;
                bool ref_ov = b->tStart < l_te && last->tStart < b_te;
                bool q_ov   = b->qStart < l_qe && last->qStart < b_qe;
                do_merge = ref_ov && q_ov;     // overlap on BOTH axes
            }
            if (do_merge) {
                int64_t last_end = last->tStart + last->size;
                int64_t b_end    = b->tStart + b->size;
                int64_t new_end  = b_end > last_end ? b_end : last_end;
                // '-' strand: qStart is the LOW query coord, held by the
                // block with the highest tStart, so it can only shrink.
                if (last->strand == '-' && b->qStart < last->qStart)
                    last->qStart = b->qStart;
                last->size = new_end - last->tStart;
                free(b->qChrom);
                free(b);
            } else {
                merged.push_back(b);
            }
        }
        cand.swap(merged);
    }

    // De-duplicate the residual redundant layer: drop any block whose
    // reference interval is fully contained in a LARGER block of the same
    // chain.  The big (merged) block already represents that reference
    // span; the nested fragment -- whether on the same diagonal or a few
    // bp off it -- is the second granularity that makes the snake loop
    // back.  Removing it leaves one consistent granularity per query with
    // NO size/pixel threshold, and coverage is unchanged (the container
    // covers the span).  Re-sort largest-first within (chain, tStart) so a
    // container is always seen before the blocks it contains.
    std::sort(cand.begin(), cand.end(),
              [](const struct taffy_block_t *a, const struct taffy_block_t *b) {
                  if (a->chainId != b->chainId) return a->chainId < b->chainId;
                  if (a->tStart  != b->tStart)  return a->tStart  < b->tStart;
                  if (a->size    != b->size)    return a->size > b->size;
                  return a->qStart < b->qStart;   // TOTAL order: build-stable
              });
    {
        std::vector<struct taffy_block_t *> tiled;
        tiled.reserve(cand.size());
        int64_t cur_chain = -1, max_end = 0;
        for (struct taffy_block_t *b : cand) {
            int64_t b_end = b->tStart + b->size;
            if (b->chainId == cur_chain) {
                if (b_end <= max_end) {            // reference-contained -> redundant
                    free(b->qChrom);
                    free(b);
                    continue;
                }
                if (b_end > max_end) max_end = b_end;
            } else {
                cur_chain = b->chainId;
                max_end   = b_end;
            }
            tiled.push_back(b);
        }
        cand.swap(tiled);
    }

    // Trim the residual partial overlaps so adjacent blocks within a chain
    // ABUT instead of sharing a few reference bases.  A snake leveler's
    // rule is binary -- a block whose target-start falls before the
    // running edge of an earlier block on its row gets bumped to a new
    // row, no matter how small the overlap or how big the block -- so a
    // 47 bp overlap can exile an entire 49 kb block to row 2.  The shared
    // bases are the chainer extending two blocks (flanking a small query
    // insertion) into the same reference; we keep them on the EARLIER
    // block and clip the later block's start up to the running edge.
    // cand is already in (chainId, tStart) order from the de-dup sort, and
    // post-de-dup every block extends past the edge, so the clip is always
    // a proper sub-trim (d < size).
    {
        int64_t cur_chain = -1, edge = 0;
        for (struct taffy_block_t *b : cand) {
            if (b->chainId != cur_chain) {
                cur_chain = b->chainId;
                edge = b->tStart + b->size;
                continue;
            }
            if (b->tStart < edge) {
                int64_t d = edge - b->tStart;        // shared reference bp
                if (d < b->size) {
                    b->tStart += d;
                    b->size   -= d;
                    if (b->strand == '+') b->qStart += d;  // '-' keeps qStart
                }
            }
            int64_t b_end = b->tStart + b->size;
            if (b_end > edge) edge = b_end;
        }
    }

    // Coarsen by ZOOM, not by how many blocks landed in the window: bin
    // iff the per-block budget would make each block span >= a coverage
    // bin.  This is panning-stable (a fixed locus renders the same at the
    // same zoom).  Detail output then emits EVERY block (span-bounded
    // count), so raise its ceiling off max_output_blocks.
    int64_t coarse_cap = H->max_output_blocks > 0 ? H->max_output_blocks : 1;
    bool coarsened = ((int64_t)(tEnd - tStart) / coarse_cap) >= BLOCK_MIN_COARSE_BIN;
    if (!coarsened) emit_cap = 1000000;     // detail: emit all (high backstop)
    if (!coarsened) {
        // Narrow query: emit every candidate at full detail, plus (for
        // QUERY_AND_TARGET_DUPS) the per-dupe-chain tRange lists.
        filter_noise(cand);   // drop sub-resolution slivers (no-op unless enabled)
        for (struct taffy_block_t *b : cand) {
            bool is_primary = (b->chainId == primary_id);
            // Build the dupe-list entry BEFORE append_mapped(b): at the
            // emit_cap backstop append_mapped frees b, so b->chainId /
            // qChrom / tStart / size must be read first.  append_mapped(b)
            // must be the LAST use of b.
            if (dupMode == TAFFY_QUERY_AND_TARGET_DUPS && !is_primary) {
                auto it = dupe_by_chain.find(b->chainId);
                if (it == dupe_by_chain.end()) {
                    taffy_target_dupe_list_t *d = (taffy_target_dupe_list_t *)
                        st_calloc(1, sizeof(*d));
                    d->id     = b->chainId;
                    d->qChrom = strdup(b->qChrom);
                    d->tRange = nullptr;
                    dupe_by_chain[b->chainId] = d;
                    it = dupe_by_chain.find(b->chainId);
                }
                taffy_target_range_t *r = (taffy_target_range_t *)
                    st_calloc(1, sizeof(*r));
                r->tStart = b->tStart;
                r->size   = b->size;
                r->next   = it->second->tRange;
                it->second->tRange = r;
            }
            append_mapped(b);
        }
    } else {
        // Wide query: bin [tStart, tEnd) into max_output_blocks bins, take
        // ONE coverage block per occupied bin (the bin's dominant/largest
        // block, sized to that bin's aligned coverage so it stays visible
        // at chromosome zoom and anchored on a real diagonal), then GROUP
        // those per-bin blocks BY DIAGONAL with the same collinear rule the
        // narrow path uses.  Grouping by diagonal -- not by qChrom -- keeps
        // a chromosome's orthologous backbone as one clean diagonal per arm
        // and gives each paralog copy (a different diagonal that dominates
        // some bins) its OWN row, instead of every copy of a qChrom winding
        // through a single row.  Per-dupe / mapBack detail is dropped here.
        int64_t span     = (int64_t) tEnd - (int64_t) tStart;
        int64_t cap      = H->max_output_blocks;
        int64_t bin_size = (span + cap - 1) / cap;
        if (bin_size < 1) bin_size = 1;
        struct BinAgg { struct taffy_block_t *dom; int64_t cov; };
        std::map<int64_t, BinAgg> bins;
        for (struct taffy_block_t *b : cand) {
            // Bin index is ABSOLUTE (b->tStart / bin_size), not relative to
            // the query start -- so a locus lands in the same bin no matter
            // where the window's edge is, keeping the coverage view stable
            // as you pan a fixed zoom (bin_size is fixed by the span).
            int64_t bin = b->tStart / bin_size;
            auto it = bins.find(bin);
            if (it == bins.end()) {
                BinAgg g; g.dom = b; g.cov = b->size; bins[bin] = g;
            } else {
                it->second.cov += b->size;
                // Deterministic dominant: larger wins; ties break by lower
                // tStart then qStart, so the bin's anchor never depends on
                // block iteration / heap order.
                struct taffy_block_t *d = it->second.dom;
                if (b->size > d->size ||
                    (b->size == d->size && b->tStart < d->tStart) ||
                    (b->size == d->size && b->tStart == d->tStart && b->qStart < d->qStart))
                    it->second.dom = b;
            }
        }
        // Materialize one coverage block per occupied bin (chainId TBD).
        std::vector<struct taffy_block_t *> coarse;
        coarse.reserve(bins.size());
        for (auto &kv : bins) {
            struct taffy_block_t *dom = kv.second.dom;
            int64_t sz = kv.second.cov < bin_size ? kv.second.cov : bin_size;
            if (sz < dom->size) sz = dom->size;
            struct taffy_block_t *nb = (struct taffy_block_t *) st_calloc(1, sizeof(*nb));
            nb->qChrom  = strdup(dom->qChrom);
            nb->tStart  = dom->tStart;
            nb->qStart  = dom->qStart;
            nb->size    = sz;
            nb->strand  = dom->strand;
            nb->chainId = 0;
            coarse.push_back(nb);
        }
        for (struct taffy_block_t *b : cand) { free(b->qChrom); free(b); }
        filter_noise(coarse);  // drop sparse / sub-resolution bins (no-op unless enabled)

        // Emit every coverage bin.  Absolute binning over [tStart,tEnd) can
        // yield cap+1 occupied bins; each is counted into coarse_summaries
        // below, so the emit ceiling must cover them all -- otherwise the
        // last bin's block is dropped by append_mapped while its summary
        // survives (an orphan chainSummary + a silently lost bin).  Bounded:
        // coarse.size() <= max_output_blocks + 1.
        emit_cap = (int64_t) coarse.size();

        // Diagonal grouping: sort by (qChrom, strand, tStart), then chain
        // consecutive bins that stay on one diagonal.  dd here is exactly
        // the offset (qStart-tStart) jump between adjacent bins, so a
        // paralog copy or the centromere step breaks to a new chain while
        // the backbone (slow drift) stays one chain.  Bins are bin_size bp
        // apart, so allow that much reference/query gap.
        std::sort(coarse.begin(), coarse.end(),
                  [](const struct taffy_block_t *a, const struct taffy_block_t *b) {
                      int c = strcmp(a->qChrom ? a->qChrom : "", b->qChrom ? b->qChrom : "");
                      if (c != 0) return c < 0;
                      if (a->strand != b->strand) return a->strand < b->strand;
                      return a->tStart < b->tStart;
                  });
        int64_t cc_gap = bin_size > BLOCK_L2_MAX_GAP ? bin_size : BLOCK_L2_MAX_GAP;
        // Per-bin offset drift along ONE diagonal can reach bin scale (a
        // single bin may hide a multi-kb indel), so allow that much dd --
        // still orders of magnitude below the Mb+ jump to a paralog copy.
        int64_t cc_tol = bin_size > BLOCK_L2_COLLINEAR_TOL ? bin_size : BLOCK_L2_COLLINEAR_TOL;
        int64_t cc_n = 0;
        std::map<int64_t, TaffyChainInfo> cc_agg;
        for (size_t i = 0; i < coarse.size(); i++) {
            bool extend = false;
            if (i > 0) {
                struct taffy_block_t *A = coarse[i - 1];
                struct taffy_block_t *B = coarse[i];
                if (A->strand == B->strand &&
                    strcmp(A->qChrom ? A->qChrom : "", B->qChrom ? B->qChrom : "") == 0) {
                    int64_t t_gap = B->tStart - (A->tStart + A->size);
                    int64_t q_gap = (A->strand == '+')
                                  ? (B->qStart - (A->qStart + A->size))
                                  : (A->qStart - (B->qStart + B->size));
                    int64_t at = t_gap < 0 ? -t_gap : t_gap;
                    int64_t aq = q_gap < 0 ? -q_gap : q_gap;
                    int64_t dd = (t_gap - q_gap) < 0 ? (q_gap - t_gap) : (t_gap - q_gap);
                    if (at <= cc_gap && aq <= cc_gap && dd <= cc_tol)
                        extend = true;
                }
            }
            if (!extend) {
                cc_n++;
                TaffyChainInfo ci; ci.id = cc_n; ci.total_score = 0;
                ci.total_bp = 0; ci.n_alns = 0; cc_agg[cc_n] = ci;
            }
            coarse[i]->chainId = cc_n;
            cc_agg[cc_n].total_bp    += coarse[i]->size;
            cc_agg[cc_n].total_score += coarse[i]->size;
            cc_agg[cc_n].n_alns      += 1;
            append_mapped(coarse[i]);
        }
        for (auto &kv : cc_agg) coarse_summaries.push_back(kv.second);
    }

    // ---- mapBackAdjacencies ---------------------------------------------
    // For each emitted (chain, qChrom) span, find the immediate qSpecies-
    // forward-coord neighbors on either side and emit any that back-map
    // to a tSpecies region OUTSIDE [tStart, tEnd) on tChrom.  Matches
    // HAL's mapAdjacencies (halBlockMapper.cpp::mapAdjacencies) but at
    // chain granularity (one flank per side per chain) rather than per
    // mapped segment -- a chain is the snake-track unit.  Skipped when the
    // output was positionally coarsened (the per-chain edge geometry the
    // flank scan needs no longer matches the emitted blocks).
    if (mapBackAdjacencies && !coarsened) {
        TuiGenomeLift *gl_t = get_or_load_gl(H, tSpecies, errStr);
        if (gl_t) {
            // Span per (chain_id, qChrom): min q_start aln-index and max
            // q_end aln-index.  Walking emitted alns suffices because the
            // chainer sorted them by q_start within each (q_name, strand)
            // partition; the first/last aln of a chain on a given t_name
            // is what we want to scan from.
            struct ChainEdge {
                int64_t left_idx;   // aln-buffer index of leftmost run
                int64_t right_idx;  // aln-buffer index of rightmost run
                std::string qChrom;
            };
            std::map<int64_t, ChainEdge> edges;
            for (int64_t i = 0; i < n; i++) {
                int64_t cid = chain_id[i];
                // For NO_DUPS we already dropped non-primary blocks;
                // restrict flanks to chains we actually emitted.
                if (dupMode == TAFFY_NO_DUPS && cid != primary_id) continue;
                // Budget cap: only emit flanks for chains we kept.
                if (!kept_chains.count(cid)) continue;
                auto it = edges.find(cid);
                if (it == edges.end()) {
                    ChainEdge e;
                    e.left_idx  = i;
                    e.right_idx = i;
                    e.qChrom    = alns[i].t_name;
                    edges[cid] = e;
                } else {
                    if (alns[i].q_start < alns[it->second.left_idx].q_start)
                        it->second.left_idx = i;
                    if (alns[i].q_end > alns[it->second.right_idx].q_end)
                        it->second.right_idx = i;
                }
            }

            for (auto &kv : edges) {
                const ChainEdge &E = kv.second;
                std::string full = std::string(qSpecies) + "." + E.qChrom;
                const std::vector<SeqRun> *runs = get_seq_runs(H, full);
                if (!runs || runs->empty()) continue;

                // Binary-search runs[] for the chain's leftmost / rightmost
                // qStart -- those are the SeqRun indices we'll walk outward
                // from.  Aln.t_start IS qSpecies forward-coord (qStart in
                // block-output terms), so lower_bound on t_start finds the
                // owning run.
                auto bsearch_t_start = [&](int64_t target) -> int64_t {
                    int64_t lo = 0, hi = (int64_t) runs->size();
                    while (lo < hi) {
                        int64_t m = lo + (hi - lo) / 2;
                        if ((*runs)[m].t_start < target) lo = m + 1;
                        else hi = m;
                    }
                    return lo;
                };
                int64_t left_qStart  = alns[E.left_idx].t_start;
                int64_t right_qEnd   = alns[E.right_idx].t_end;
                int64_t left_idx_runs  = bsearch_t_start(left_qStart);
                int64_t right_idx_runs = bsearch_t_start(right_qEnd) - 1;
                if (right_idx_runs < 0) right_idx_runs = 0;

                struct taffy_block_t *L =
                    find_off_screen_flank(runs, left_idx_runs,  -1,
                                          gl_t, E.qChrom.c_str(),
                                          tChrom, tStart, tEnd);
                if (L) append_mapped(L);
                struct taffy_block_t *R =
                    find_off_screen_flank(runs, right_idx_runs, +1,
                                          gl_t, E.qChrom.c_str(),
                                          tChrom, tStart, tEnd);
                if (R) append_mapped(R);
            }
        }
        // Note: gl_t isn't released here -- it stays in lift_cache for
        // later browser pans against the same tSpecies.
    }

    // Link dupe entries into a single linked list.
    struct taffy_target_dupe_list_t *dupe_head = nullptr, *dupe_tail = nullptr;
    for (auto &kv : dupe_by_chain) {
        if (!dupe_head) dupe_head = kv.second;
        if (dupe_tail) dupe_tail->next = kv.second;
        dupe_tail = kv.second;
    }

    res->mappedBlocks     = mapped_head;
    res->targetDupeBlocks = dupe_head;

    // Build chainSummaries, score-desc.  Coarse path: one per query
    // chromosome (the per-qChrom groups built above).  Detail path: one
    // per real chain that actually appears in mappedBlocks (not merely
    // kept_chains).
    struct taffy_chain_summary_t *cs_head = nullptr, *cs_tail = nullptr;
    auto push_summary = [&](int64_t id, int64_t score, int64_t bp, int64_t na) {
        struct taffy_chain_summary_t *cs = (struct taffy_chain_summary_t *)
            st_calloc(1, sizeof(*cs));
        cs->id = id; cs->totalScore = score; cs->totalBp = bp; cs->nAlns = na;
        cs->next = nullptr;
        if (!cs_head) cs_head = cs;
        if (cs_tail) cs_tail->next = cs;
        cs_tail = cs;
    };
    if (coarsened) {
        std::sort(coarse_summaries.begin(), coarse_summaries.end(),
                  [](const TaffyChainInfo &a, const TaffyChainInfo &b) {
                      if (a.total_score != b.total_score) return a.total_score > b.total_score;
                      return a.id < b.id;
                  });
        for (const TaffyChainInfo &ci : coarse_summaries)
            push_summary(ci.id, ci.total_score, ci.total_bp, ci.n_alns);
    } else {
        std::set<int64_t> emitted_chains;
        for (struct taffy_block_t *b = mapped_head; b; b = b->next)
            if (b->chainId != 0) emitted_chains.insert(b->chainId);
        for (int64_t k = 0; k < n_chains; k++) {
            if (!emitted_chains.count(chains[k].id)) continue;
            push_summary(chains[k].id, chains[k].total_score,
                         chains[k].total_bp, chains[k].n_alns);
        }
    }
    res->chainSummaries = cs_head;

    return res;
}

extern "C" struct taffy_block_results_t *
taffyGetBlocksInTargetRange(int h, const char *qSpecies, const char *tSpecies,
                            const char *tChrom, taffy_int_t tStart, taffy_int_t tEnd,
                            taffy_int_t tReversed, taffy_seqmode_type_t seqMode,
                            taffy_dup_type_t dupMode, int mapBackAdjacencies,
                            const char *coalescenceLimitName, char **errStr) {
    return get_blocks_impl(h, qSpecies, tSpecies, tChrom, tStart, tEnd,
                           tReversed, seqMode, dupMode, mapBackAdjacencies,
                           /*qChrom=*/ nullptr, coalescenceLimitName, errStr);
}

extern "C" struct taffy_block_results_t *
taffyGetBlocksInTargetRange_filterByChrom(int h, const char *qSpecies, const char *tSpecies,
                                          const char *tChrom, taffy_int_t tStart, taffy_int_t tEnd,
                                          taffy_int_t tReversed, taffy_seqmode_type_t seqMode,
                                          taffy_dup_type_t dupMode, int mapBackAdjacencies,
                                          const char *qChrom, const char *coalescenceLimitName,
                                          char **errStr) {
    return get_blocks_impl(h, qSpecies, tSpecies, tChrom, tStart, tEnd,
                           tReversed, seqMode, dupMode, mapBackAdjacencies,
                           qChrom, coalescenceLimitName, errStr);
}

/* ------------------------------------------------------------------ */
/* DNA stub                                                            */
/* ------------------------------------------------------------------ */

extern "C" char *taffyGetDna(int h, const char *species, const char *chrom,
                             taffy_int_t start, taffy_int_t end, char **errStr) {
    (void) h; (void) species; (void) chrom;
    if (end < start) { set_err(errStr, "taffyBlockViz: bad start/end"); return nullptr; }
    int64_t len = (int64_t)(end - start);
    char *s = (char *) st_malloc((size_t) len + 1);
    memset(s, 'N', (size_t) len);
    s[len] = 0;
    return s;
}

/* ------------------------------------------------------------------ */
/* Genome metadata stub                                                */
/* ------------------------------------------------------------------ */

extern "C" struct taffy_metadata_t *taffyGetGenomeMetadata(int h, const char *genome, char **errStr) {
    (void) h; (void) genome; (void) errStr;
    return nullptr;  // no .tui-side metadata in initial cut
}
