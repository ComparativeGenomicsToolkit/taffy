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
    int64_t     q_id, t_id;     // interned q_name/t_name ids; set by taffy_chain (callers leave unset)
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

/* Default constant sets in taffy chaining:
 *
 * 1. TAFFY_CHAIN_DEFAULT_OPEN/EXTEND/MAX_GAP (below): the CANONICAL
 *    taffy preset.  open=0 + extend=1 chains aggressively (any forward-
 *    collinear runs join; nothing is dropped on cost alone); max_gap
 *    of 10 Mb breaks chains only across truly large rearrangements.
 *    Used by taffyBlockViz get_blocks_impl, taffy lift --chainFilter,
 *    view_chain_dup_filter.  Tuned for short, dense, gap-free runs as
 *    emitted by the .tui visitor.
 *
 * 2. taffy_chain_default_gap_cost()'s NULL-params fallback: open=5000,
 *    extend=1.  Legacy: matches paffy's whole-genome-PAF chaining
 *    defaults where per-aln scores are kbs.  Only fires when callers
 *    pass NULL for gap_cost_params -- callers using
 *    TaffyChainCostParams get the (1) preset via TAFFY_CHAIN_DEFAULT_*.
 *
 * 3. taf_chain.c's --maxGap=1000 default: sidecar-specific (axtChain-
 *    equivalent).  taf_chain.c does NOT use this gap-cost machinery --
 *    it's a streaming greedy merge with a symmetric maxGap on both
 *    axes, no scoring, no top-N.  The constant lives in that file's
 *    CLI parser and intentionally matches axtChain so output is
 *    comparable to UCSC chain files.
 */
#define TAFFY_CHAIN_DEFAULT_OPEN     ((int64_t) 0)
#define TAFFY_CHAIN_DEFAULT_EXTEND   ((int64_t) 1)
#define TAFFY_CHAIN_DEFAULT_MAX_GAP  ((int64_t) 10 * 1000 * 1000)

int64_t taffy_chain_default_gap_cost(int64_t q_gap, int64_t t_gap, void *params);

/* Post-chain collinear merge.  After taffy_chain has assigned chain
 * IDs, this folds adjacent alns within the same chain whose extents
 * abut on BOTH axes (forward-collinear runs) into a single merged
 * entry.  Without this, the browser snake renderer / view block
 * emitter can produce thousands of one-block-per-run outputs for
 * what is logically a single gap-free alignment.
 *
 * MERGE RULE: two alns A (earlier in q) and B (later) merge iff
 *   - chain_id[A] == chain_id[B]
 *   - A.q_end == B.q_start (q-axis abut, strict)
 *   - + strand: A.t_end   == B.t_start  (t-axis abut, forward)
 *   - - strand: B.t_end   == A.t_start  (t-axis abut, reverse)
 *   - same (q_name, t_name, strand) -- guaranteed within a chain by
 *     taffy_chain's invariant, so checked only via chain_id match.
 *
 * On merge: the EARLIER aln (A) is kept and its extents grow to cover
 * both; its `score` is accumulated.  The LATER aln (B) is logically
 * removed.  Callers that maintain per-aln auxiliary records (eg
 * blockViz's taffy_block_t* via TaffyAln.user) supply an `on_merge`
 * callback to fold or free those records.  Pass NULL to skip.
 *
 * Returns the count of surviving alns (n - n_merges).  alns[] is
 * mutated in place; the removed entries' contents are undefined --
 * downstream code must use the returned count, OR walk all n entries
 * and detect removal via the user-supplied marker (eg user==NULL).
 */
int64_t taffy_chain_merge_collinear(TaffyAln *alns, int64_t n,
                                    const int64_t *chain_id,
                                    void (*on_merge)(TaffyAln *kept,
                                                     TaffyAln *absorbed,
                                                     void *user),
                                    void *user);

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

/* Overlap-aware paralogy filter: pick chain survivors based on q-axis
 * overlap with already-kept chains.  Walks chains_out in (score desc)
 * order; for each candidate chain, computes its union-of-aln q-coverage
 * and accepts iff the q-bp intersection with the running union of kept
 * chains' q-coverage is at most overlap_frac of the candidate's own
 * q-bp.  Paralogs (same q bp mapped to multiple t loci) hit ~100% and
 * drop; inversions / disjoint-q runs hit 0% and stay.
 *
 * Inputs:
 *   alns/n           result of a prior taffy_chain() call.  alns are
 *                    expected to be (q_name, strand, q_start)-sorted
 *                    -- if you just returned from taffy_chain that is
 *                    already true.
 *   chain_id         length-n; the chain id per aln (1-based).
 *   chains_out       sorted by score desc (as taffy_chain returns).
 *   n_chains_out     length of chains_out.
 *   max_id           max chain id; keep_chain must have (max_id + 1)
 *                    bytes.
 *   overlap_frac     threshold in [0, 1].  0 = drop on ANY overlap.
 *                    Negative is treated as "off" by the caller; this
 *                    function doesn't check.
 *   cap              optional safety cap on survivor count.  0 = no cap.
 *
 * Output:
 *   keep_chain[cid]  set to 1 for accepted chain ids, otherwise 0.
 *                    Caller must zero-init (st_calloc) before calling.
 *
 * Complexity: O(n log n) total in the typical case (CSR bucketing in
 * O(n), merged kept-union maintained in O(m) per merge with proper
 * 2-list union).  Tractable at chromosome scale.
 */
void taffy_chain_overlap_frac_select(const TaffyAln *alns, int64_t n,
                                     const int64_t *chain_id,
                                     const TaffyChainInfo *chains_out,
                                     int64_t n_chains_out,
                                     int64_t max_id,
                                     double overlap_frac,
                                     int64_t cap,
                                     char *keep_chain);

#ifdef __cplusplus
}
#endif

#endif /* TAF_CHAIN_H_ */
