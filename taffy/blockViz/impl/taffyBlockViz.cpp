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

struct TaffyHandle {
    Tui *tui = nullptr;
    std::string tui_path_str;   // resolved .tui path (for tui_sequence_lengths)
    std::map<std::string, TuiGenomeLift *> lift_cache;
    // Chain tuning -- defaults match TAFFY_CHAIN_DEFAULT_{OPEN,EXTEND,MAX_GAP}.
    // Tunable per-handle via taffySetChainParams.
    int64_t chain_open      = TAFFY_CHAIN_DEFAULT_OPEN;
    int64_t chain_extend    = TAFFY_CHAIN_DEFAULT_EXTEND;
    int64_t max_gap_length  = TAFFY_CHAIN_DEFAULT_MAX_GAP;
};

static std::map<int, TaffyHandle *> g_handles;
static int  g_next_handle = 1;
static std::mutex g_mutex;

static void set_err(char **errStr, const char *msg) {
    if (errStr) *errStr = strdup(msg);
}

static TaffyHandle *get_handle(int h, char **errStr) {
    auto it = g_handles.find(h);
    if (it == g_handles.end()) {
        set_err(errStr, "taffyBlockViz: invalid handle");
        return nullptr;
    }
    return it->second;
}

// Cache the per-genome lift table across calls (browser pans + zooms
// hit the same target/query pair repeatedly).  Caller must hold g_mutex.
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
    // tui_path() resolves either <foo>.maf.gz -> <foo>.maf.gz.tui or
    // accepts an already-fully-qualified .tui path.
    char *p = tui_path(path);
    Tui *tui = tui_load(p);
    if (!tui) {
        std::string msg = std::string("taffyBlockViz: tui_load failed for ") + (p ? p : path);
        set_err(errStr, msg.c_str());
        free(p);
        return -1;
    }

    std::lock_guard<std::mutex> lock(g_mutex);
    int h = g_next_handle++;
    TaffyHandle *H = new TaffyHandle();
    H->tui = tui;
    H->tui_path_str = p;
    free(p);
    g_handles[h] = H;
    return h;
}

extern "C" int taffyClose(int h, char **errStr) {
    std::lock_guard<std::mutex> lock(g_mutex);
    auto it = g_handles.find(h);
    if (it == g_handles.end()) {
        set_err(errStr, "taffyBlockViz: invalid handle");
        return -1;
    }
    TaffyHandle *H = it->second;
    for (auto &kv : H->lift_cache) tui_genome_lift_destruct(kv.second);
    tui_destruct(H->tui);
    delete H;
    g_handles.erase(it);
    return 0;
}

extern "C" int taffyCloseGenome(int h, const char *genome, char **errStr) {
    if (!genome) { set_err(errStr, "taffyBlockViz: NULL genome"); return -1; }
    std::lock_guard<std::mutex> lock(g_mutex);
    TaffyHandle *H = get_handle(h, errStr);
    if (!H) return -1;
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
    std::lock_guard<std::mutex> lock(g_mutex);
    TaffyHandle *H = get_handle(h, errStr);
    if (!H) return -1;
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
    std::lock_guard<std::mutex> lock(g_mutex);
    TaffyHandle *H = get_handle(h, errStr);
    if (!H) return -1;
    if (chain_open)     *chain_open     = H->chain_open;
    if (chain_extend)   *chain_extend   = H->chain_extend;
    if (max_gap_length) *max_gap_length = H->max_gap_length;
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

extern "C" void taffyFreeBlockResults(struct taffy_block_results_t *res) {
    if (!res) return;
    taffyFreeBlocks(res->mappedBlocks);
    taffyFreeTargetDupeLists(res->targetDupeBlocks);
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
// (tui_genome_names) when present.  Older .tui files without g-records
// fall back to a first-dot split of "<genome>.<seq>" d-line keys, which
// is wrong for genome names containing dots (NCBI accessions like
// "GCA_028858775.2") but the only option short of a rebuild.
struct GenomeRow { std::string genome; taffy_int_t total_bp = 0; taffy_int_t n_chroms = 0; };

static std::map<std::string, GenomeRow> enumerate_genomes(TaffyHandle *H) {
    std::map<std::string, GenomeRow> out;

    int64_t n_g = 0;
    TuiGenomeInfo *roster = tui_genome_names(H->tui_path_str.c_str(), &n_g);
    if (roster != nullptr && n_g > 0) {
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

    // Fallback for pre-roster .tui: split d-line keys on the FIRST dot.
    // Misclassifies dotted genome names; reported as a known limitation
    // in the taffyGetSpecies header doc.
    stHash *seqs = tui_sequence_lengths(H->tui_path_str.c_str());
    if (!seqs) return out;
    stHashIterator *it = stHash_getIterator(seqs);
    char *k;
    while ((k = (char *) stHash_getNext(it)) != NULL) {
        char *dot = strchr(k, '.');
        if (!dot) continue;
        std::string genome(k, dot - k);
        int64_t len = (int64_t)(intptr_t) stHash_search(seqs, k);
        GenomeRow &row = out[genome];
        if (row.genome.empty()) row.genome = genome;
        row.total_bp += len;
        row.n_chroms++;
    }
    stHash_destructIterator(it);
    stHash_destruct(seqs);
    return out;
}

extern "C" struct taffy_species_t *taffyGetSpecies(int h, char **errStr) {
    std::lock_guard<std::mutex> lock(g_mutex);
    TaffyHandle *H = get_handle(h, errStr);
    if (!H) return nullptr;
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
    std::lock_guard<std::mutex> lock(g_mutex);
    TaffyHandle *H = get_handle(h, errStr);
    if (!H) return nullptr;

    stHash *seqs = tui_sequence_lengths(H->tui_path_str.c_str());
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
/* Block query (the hot path)                                          */
/* ------------------------------------------------------------------ */

struct BlockCtx {
    // Buffer one taffy_block_t plus a parallel TaffyAln per visited run.
    // Linking into the output list happens after taffy_chain has assigned
    // chain IDs and we know how to route each block per the dupMode.
    std::vector<struct taffy_block_t *> blocks;
    std::vector<TaffyAln>                alns;
    // Intern owned strings (qChroms) so TaffyAln const char * stays
    // valid through the chain call.  std::set guarantees that
    // references to elements aren't invalidated by future insertions,
    // so the c_str() pointer we stash in TaffyAln.t_name is stable.
    std::set<std::string>                qchrom_interns;
    std::string                          tFullName;
    int64_t c_lo = 0;             // column-start of current interval
    int64_t tpos_at_c_lo = 0;     // tSpecies position at c_lo
    const char *qChromFilter = nullptr;  // optional filter
};

// Visitor: clip the q run, compute (tStart, qStart, size, strand),
// and buffer both a taffy_block_t (for the eventual output linked list)
// and a TaffyAln (for the chain call after the visit loop).
static void block_visit_cb(const TuiRun *r, void *user) {
    BlockCtx *cx = (BlockCtx *) user;
    if (cx->qChromFilter && strcmp(r->seq, cx->qChromFilter) != 0) return;

    int64_t r_end = r->g_start + r->length;
    int64_t cs = r->g_start > cx->c_lo ? r->g_start : cx->c_lo;
    int64_t ce = r_end;

    // tSpecies position at column `cs`: linear within an interval, since
    // an interval is a contiguous run of tSpecies bases.
    int64_t tStart = cx->tpos_at_c_lo + (cs - cx->c_lo);

    // Compute qStart honouring strand (mirrors taf_lift.c chunk_lift_visit_cb).
    int64_t qStart;
    if (r->strand) {
        qStart = r->t_start + (cs - r->g_start);
    } else {
        qStart = r->t_start + r->length - (ce - r->g_start);
    }
    int64_t size = ce - cs;

    struct taffy_block_t *b = (struct taffy_block_t *) st_calloc(1, sizeof(*b));
    b->qChrom = strdup(r->seq);
    b->tStart = tStart;
    b->qStart = qStart;
    b->size   = size;
    b->strand = r->strand ? '+' : '-';
    cx->blocks.push_back(b);

    // Intern qChrom; set::insert returns iterator stable across future
    // inserts so .c_str() remains valid for the rest of the visit.
    auto ins = cx->qchrom_interns.emplace(r->seq);
    const char *qc_interned = ins.first->c_str();

    TaffyAln a = {0};
    a.q_name  = cx->tFullName.c_str();        // single tSpecies.tChrom for whole query
    a.q_start = tStart;
    a.q_end   = tStart + size;
    a.t_name  = qc_interned;
    a.t_start = qStart;
    a.t_end   = qStart + size;
    a.strand  = r->strand ? +1 : -1;
    a.score   = size;
    a.user    = (void *) b;
    cx->alns.push_back(a);
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
    if (mapBackAdjacencies != 0) { set_err(errStr, "taffyBlockViz: mapBackAdjacencies not supported"); return nullptr; }
    if (coalescenceLimitName) { set_err(errStr, "taffyBlockViz: coalescenceLimitName not supported"); return nullptr; }
    if (seqMode != TAFFY_NO_SEQUENCES) { set_err(errStr, "taffyBlockViz: seqMode != NO_SEQUENCES not supported"); return nullptr; }
    if (dupMode != TAFFY_NO_DUPS && dupMode != TAFFY_QUERY_DUPS && dupMode != TAFFY_QUERY_AND_TARGET_DUPS) {
        set_err(errStr, "taffyBlockViz: unknown dupMode");
        return nullptr;
    }
    if (!qSpecies || !tSpecies || !tChrom) { set_err(errStr, "taffyBlockViz: missing required arg"); return nullptr; }
    if (tEnd < 0 || (tEnd > 0 && tEnd <= tStart)) { set_err(errStr, "taffyBlockViz: bad tStart/tEnd"); return nullptr; }

    std::lock_guard<std::mutex> lock(g_mutex);
    TaffyHandle *H = get_handle(h, errStr);
    if (!H) return nullptr;

    // Resolve the .tui d-line key for tChrom: "<tSpecies>.<tChrom>".
    std::string tFullName = std::string(tSpecies) + "." + tChrom;

    // tEnd == 0 -> use full chrom length.
    if (tEnd == 0) {
        stHash *seqs = tui_sequence_lengths(H->tui_path_str.c_str());
        if (seqs) {
            tEnd = (int64_t)(intptr_t) stHash_search(seqs, (void *) tFullName.c_str());
            stHash_destruct(seqs);
        }
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

    // Walk each interval; buffer one taffy_block_t + TaffyAln per visited
    // qSpecies run.  Output linking happens AFTER chaining.
    BlockCtx cx;
    cx.qChromFilter = qChrom;
    cx.tFullName    = tFullName;  // stable interned string for chain partition
    int64_t cum_tpos = tStart;
    for (int64_t k = 0; k < n_iv; k++) {
        cx.c_lo = iv[k].start;
        cx.tpos_at_c_lo = cum_tpos;
        tui_genome_lift_visit_runs(gl_q, iv[k].start, iv[k].end,
                                   block_visit_cb, &cx);
        cum_tpos += (iv[k].end - iv[k].start);
    }
    free(iv);

    struct taffy_block_results_t *res = (struct taffy_block_results_t *) st_calloc(1, sizeof(*res));
    int64_t n = (int64_t) cx.alns.size();
    if (n == 0) return res;  // empty query, nothing to chain

    // Run the chainer.  Tuning rationale: browser snake tracks want to
    // join collinear runs aggressively (one snake per real alignment),
    // and the universal MAF's runs are gap-free by construction so the
    // q_gap + t_gap between adjacent runs of the same alignment is
    // usually small.  Defaults (TAFFY_CHAIN_DEFAULT_OPEN/EXTEND/MAX_GAP)
    // chain anything in the same syntenic block while keeping truly
    // distant paralogs in their own chain; callers can override via
    // taffySetChainParams.
    TaffyChainCostParams cost = { H->chain_open, H->chain_extend };
    std::vector<int64_t> chain_id((size_t) n, 0);
    TaffyChainInfo *chains = nullptr;
    int64_t n_chains = 0;
    taffy_chain(cx.alns.data(), n,
                taffy_chain_default_gap_cost, &cost,
                H->max_gap_length,
                chain_id.data(), &chains, &n_chains);

    // Primary chain = chains[0] (sorted desc by score).
    int64_t primary_id = (n_chains > 0) ? chains[0].id : 0;

    // Walk alns (post-chain order) and route blocks per dupMode.  Each
    // aln's `user` field is the taffy_block_t* we created in the
    // visitor; chain_id[i] is its chain.
    struct taffy_block_t *mapped_head = nullptr, *mapped_tail = nullptr;
    auto append_mapped = [&](struct taffy_block_t *b) {
        b->next = nullptr;
        if (!mapped_head) mapped_head = b;
        if (mapped_tail) mapped_tail->next = b;
        mapped_tail = b;
    };
    // For QUERY_AND_TARGET_DUPS: collect non-primary alns grouped by
    // chain_id so each non-primary chain becomes one dupe-list entry.
    std::map<int64_t, taffy_target_dupe_list_t *> dupe_by_chain;

    for (int64_t i = 0; i < n; i++) {
        struct taffy_block_t *b = (struct taffy_block_t *) cx.alns[i].user;
        int64_t cid = chain_id[i];
        bool is_primary = (cid == primary_id);

        if (dupMode == TAFFY_NO_DUPS && !is_primary) {
            // Drop non-primary blocks entirely.
            free(b->qChrom);
            free(b);
            continue;
        }

        // Always append to mappedBlocks under QUERY_DUPS / QUERY_AND_TARGET_DUPS.
        append_mapped(b);

        if (dupMode == TAFFY_QUERY_AND_TARGET_DUPS && !is_primary) {
            // Also build a dupe-list entry for this non-primary chain.
            auto it = dupe_by_chain.find(cid);
            if (it == dupe_by_chain.end()) {
                taffy_target_dupe_list_t *d = (taffy_target_dupe_list_t *)
                    st_calloc(1, sizeof(*d));
                d->id     = cid;
                d->qChrom = strdup(b->qChrom);  // own a copy; mappedBlocks owns its own
                d->tRange = nullptr;
                dupe_by_chain[cid] = d;
                it = dupe_by_chain.find(cid);
            }
            // Append tRange node (in source-coord order; aln list is
            // already q_start-sorted within the chain by the sweep).
            taffy_target_range_t *r = (taffy_target_range_t *)
                st_calloc(1, sizeof(*r));
            r->tStart = b->tStart;
            r->size   = b->size;
            r->next   = it->second->tRange;   // prepend (we'll reverse below if needed)
            it->second->tRange = r;
        }
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
    free(chains);
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
