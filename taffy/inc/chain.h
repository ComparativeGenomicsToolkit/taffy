/*
 * Generic alignment chaining for taffy.
 *
 * Ports the DP sweep algorithm from `paffy/impl/chaining.c` into a
 * taffy-internal module so both taffyBlockViz (snake-track dupe
 * filtering) and `taffy lift` (the future --chainFilter option) can
 * share one principled chainer instead of inventing their own
 * heuristics or copying HAL's hacky greedy stack + 2.5%-of-total
 * fixed-threshold filter.
 *
 * INPUT: an array of TaffyAln intervals (forward-strand half-open
 * coords on both axes; strand carried separately).  Each aln carries
 * a back-pointer `user` to whatever the caller's per-aln record is
 * (TuiRun, Paf, etc.) so we don't impose a representation.
 *
 * OUTPUT: a 1-based chain ID per aln + an array of per-chain stats
 * (score, bp, n_alns) sorted by score descending -- so chains_out[0]
 * is the "primary" chain.  Callers ranking dupes can pick the top N
 * chains, route the rest into dupe lists, etc.
 *
 * Caller arranges:
 *   - score per aln (in TaffyAln.score) -- typically the aln length.
 *     Chaining maximizes (sum of scores) - (sum of gap_costs).
 *   - gap_cost callback (lastz-style chain_open + chain_extend gap
 *     is provided as taffy_chain_default_gap_cost).
 *   - max_gap_length: any pair of alns whose gap exceeds this on
 *     either axis breaks the chain.  Set to INT64_MAX to disable.
 *
 * COMPLEXITY: O(n log n) sort + O(n * k) sweep where k is the average
 * number of predecessor chains scanned per aln.  In practice k stays
 * small because the sweep prunes chains that exceed max_gap_length on
 * the query axis.
 */

#ifndef TAF_CHAIN_H_
#define TAF_CHAIN_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* One alignment interval to feed the chainer.  All coords are
 * forward-strand half-open ([start, end)).  `strand` is +1 or -1;
 * chains never join alns of opposite strands or differing seq names. */
typedef struct {
    const char *q_name;
    int64_t     q_start, q_end;
    const char *t_name;
    int64_t     t_start, t_end;
    int         strand;        // +1 or -1
    int64_t     score;          // per-aln score; typically aln length
    void       *user;           // back-pointer to caller's record
} TaffyAln;

/* Per-chain summary returned by taffy_chain.  `id` is 1-based and
 * matches the chain_id[] output array.  Sorted by total_score
 * descending so chains[0] is the primary. */
typedef struct {
    int64_t id;
    int64_t total_score;       // sum(aln.score) - sum(gap_cost)
    int64_t total_bp;          // sum(aln.q_end - aln.q_start)
    int64_t n_alns;            // # alns assigned to this chain
} TaffyChainInfo;

/* Default lastz-style gap cost:
 *   chain_open + chain_extend * (q_gap + t_gap)   when (q_gap + t_gap) > 0
 *   0                                              otherwise
 * Tunable via TaffyChainCostParams passed as gap_cost_params; pass
 * NULL for defaults (open=5000, extend=1, matching paffy's defaults). */
typedef struct {
    int64_t chain_open;
    int64_t chain_extend;
} TaffyChainCostParams;

/* Browser/blockViz tuning preset, shared by taffyBlockViz's
 * get_blocks_impl and taffy lift --chainFilter so the constants live
 * in one place.  open=0 + extend=1 chains aggressively (any forward-
 * collinear runs join; nothing is dropped on cost alone); max_gap of
 * 10 Mb breaks chains only across truly large rearrangements.  These
 * differ from the lastz-style 5000/1 defaults baked into
 * taffy_chain_default_gap_cost() -- those are tuned for whole-genome
 * PAF chaining where alignment scores are kbs, not tens of bp.
 */
#define TAFFY_CHAIN_DEFAULT_OPEN     ((int64_t) 0)
#define TAFFY_CHAIN_DEFAULT_EXTEND   ((int64_t) 1)
#define TAFFY_CHAIN_DEFAULT_MAX_GAP  ((int64_t) 10 * 1000 * 1000)

int64_t taffy_chain_default_gap_cost(int64_t q_gap, int64_t t_gap, void *params);

/* Chain alignments and assign chain IDs.
 *
 * `alns` is sorted in place by (q_name, strand, q_start) -- so the
 * caller's original ordering is destroyed but the back-pointer
 * (alns[i].user) lets the caller recover its own record.
 *
 * `chain_id[i]` is the 1-based ID for alns[i] (one chain may contain
 * multiple alns; conversely a single-aln chain is still a chain with
 * its own id).
 *
 * `*chains_out` is a malloc'd array of TaffyChainInfo, length
 * `*n_chains_out`, sorted by total_score descending.  Caller frees.
 *
 * Implementation notes:
 *   - alns are partitioned by (q_name, strand) into independent
 *     sub-problems.  Reverse-strand alns get q_start/q_end mirrored
 *     internally (so the sweep stays forward-only) and restored to
 *     forward coords on return.
 *   - chain scoring: total_score = sum(aln.score) -
 *     sum_{chain links} gap_cost(q_gap, t_gap).  An aln joins a
 *     predecessor chain only if doing so increases score AND the
 *     gap_cost is less than the new aln's own score (otherwise
 *     starting a new chain is better).
 */
void taffy_chain(TaffyAln *alns, int64_t n,
                 int64_t (*gap_cost)(int64_t, int64_t, void *), void *gap_cost_params,
                 int64_t max_gap_length,
                 int64_t *chain_id,
                 TaffyChainInfo **chains_out, int64_t *n_chains_out);

#ifdef __cplusplus
}
#endif

#endif /* TAF_CHAIN_H_ */
