/*
 * Generic alignment chaining for taffy.
 *
 * Faithful port of paffy/impl/chaining.c's DP sweep algorithm,
 * generalized away from the Paf-specific record type so the same
 * routine drives taffyBlockViz's dupe filter AND a future
 * `taffy lift --chainFilter`.  See taffy/inc/chain.h for the API.
 *
 * Algorithm (per (q_name, strand) partition):
 *   - Sort input alns by q_start.
 *   - Maintain a sorted set of "active chains" keyed on
 *     (t_name, t_end, q_end).
 *   - For each new aln in q_start order:
 *       * walk active chains backward (predecessor lookup);
 *       * drop predecessors whose q_gap exceeds max_gap_length
 *         (they can never be extended further);
 *       * for each viable predecessor, compute
 *           new_score = aln.score + chain.score
 *                       - gap_cost(q_gap, t_gap)
 *         keep best new_score IF the gap_cost is less than the new
 *         aln's own score (otherwise starting a new chain is better);
 *       * insert the new chain (joined or standalone) into the active
 *         set and the global all-chains set.
 *   - At the end, pop chains in score order high -> low; for each
 *     popped chain follow pChain links back to assign all its alns
 *     the chain ID; if a link was already claimed by a higher-scoring
 *     chain, cut the chain there.
 *
 * +/- strand handling: reverse-strand alns get q_start/q_end
 * "mirrored" (negated) so the sweep stays forward-monotone; we
 * restore forward coords before returning.
 */

#include "chain.h"
#include "sonLib.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Default gap cost                                                    */
/* ------------------------------------------------------------------ */

int64_t taffy_chain_default_gap_cost(int64_t q_gap, int64_t t_gap, void *params) {
    int64_t open = 5000, ext = 1;
    if (params != NULL) {
        TaffyChainCostParams *p = (TaffyChainCostParams *) params;
        open = p->chain_open;
        ext  = p->chain_extend;
    }
    if (q_gap + t_gap == 0) return 0;
    return open + ext * (q_gap + t_gap);
}

/* ------------------------------------------------------------------ */
/* Internal types                                                      */
/* ------------------------------------------------------------------ */

/* One node in the chain DP graph.  Multiple Chain nodes can share an
 * `aln` pointer in degenerate cases (none reachable from the public
 * API), but in practice each aln gets exactly one node.  Each node
 * has a back-link `pChain` to its best predecessor (or NULL if it's
 * a chain root). */
typedef struct _ChainNode {
    TaffyAln          *aln;
    int64_t            score;        // best score reaching this aln
    struct _ChainNode *pChain;       // predecessor link, or NULL
} ChainNode;

/* ------------------------------------------------------------------ */
/* Comparators (sonLib stSortedSet)                                   */
/* ------------------------------------------------------------------ */

static int intcmp(int64_t a, int64_t b) { return a > b ? 1 : (a < b ? -1 : 0); }

/* Sort input alns by (q_name, strand, q_start).  Stable on equal keys
 * via pointer tie-break -- matters because the sweep iterates in this
 * order. */
static int aln_cmp_by_q_loc(const void *a, const void *b) {
    const TaffyAln *x = a, *y = b;
    int c = strcmp(x->q_name ? x->q_name : "", y->q_name ? y->q_name : "");
    if (c) return c;
    if (x->strand != y->strand) return x->strand < y->strand ? -1 : 1;
    c = intcmp(x->q_start, y->q_start);
    if (c) return c;
    c = intcmp((int64_t)(intptr_t) x, (int64_t)(intptr_t) y);
    return c;
}

/* Active-chains key: (t_name, t_end, q_end), with pointer tiebreak. */
static int chain_cmp_by_t_loc(const void *a, const void *b) {
    const ChainNode *c1 = a, *c2 = b;
    const TaffyAln *p1 = c1->aln, *p2 = c2->aln;
    int c = strcmp(p1->t_name ? p1->t_name : "", p2->t_name ? p2->t_name : "");
    if (c) return c;
    c = intcmp(p1->t_end, p2->t_end);
    if (c) return c;
    c = intcmp(p1->q_end, p2->q_end);
    if (c) return c;
    return intcmp((int64_t)(intptr_t) c1, (int64_t)(intptr_t) c2);
}

/* All-chains key: by score (so popping last == highest). */
static int chain_cmp_by_score(const void *a, const void *b) {
    const ChainNode *c1 = a, *c2 = b;
    int r = intcmp(c1->score, c2->score);
    if (r) return r;
    return intcmp((int64_t)(intptr_t) c1, (int64_t)(intptr_t) c2);
}

/* ------------------------------------------------------------------ */
/* Predecessor iterator                                                */
/* ------------------------------------------------------------------ */

/* Find the position in active_set that's <= `chain` if we treat
 * chain's (t_end, q_end) as (t_start, q_start) for the purpose of
 * lookup.  Returns an iterator positioned so stSortedSet_getPrevious
 * yields the closest predecessor, then walks backward through
 * progressively earlier ones (paffy's same trick of temporarily
 * mutating the lookup key in the search comparator). */
static stSortedSetIterator *predecessor_iter(stSortedSet *active, ChainNode *cn) {
    TaffyAln *p = cn->aln;
    int64_t q_end_save = p->q_end, t_end_save = p->t_end;
    p->q_end = p->q_start;
    p->t_end = p->t_start;
    ChainNode *cn2 = stSortedSet_searchLessThanOrEqual(active, cn);
    p->q_end = q_end_save;
    p->t_end = t_end_save;
    if (cn2 == NULL) return stSortedSet_getIterator(active);
    stSortedSetIterator *it = stSortedSet_getIteratorFrom(active, cn2);
    /* same nudge as paffy: advance past so the very next
     * stSortedSet_getPrevious returns cn2 itself. */
    stSortedSet_getNext(it);
    stSortedSet_getNext(it);
    return it;
}

/* ------------------------------------------------------------------ */
/* Per-(q_name, strand) sweep                                          */
/* ------------------------------------------------------------------ */

/* Run the DP sweep on a contiguous run of alns sharing (q_name,
 * strand).  Writes per-aln chain assignments into the array starting
 * at chain_id_base.  Returns the new starting chain_id for the next
 * call. */
static int64_t chain_partition(TaffyAln *alns, int64_t n,
                               int64_t (*gap_cost)(int64_t, int64_t, void *),
                               void *gap_cost_params,
                               int64_t max_gap_length,
                               int64_t *chain_id_out,
                               TaffyChainInfo **chains_out, int64_t *n_chains_out,
                               int64_t next_chain_id) {
    if (n == 0) return next_chain_id;

    /* alns is already sorted by q_start by the caller. */

    stSortedSet *active = stSortedSet_construct3(chain_cmp_by_t_loc, NULL);
    /* NULL element-destructor: ID assignment below stSortedSet_remove()s
     * every node from `chains` (remove does NOT invoke the destructor), so
     * the set never frees them -- we free them explicitly via nodes_by_idx. */
    stSortedSet *chains = stSortedSet_construct3(chain_cmp_by_score, NULL);
    stList *to_remove   = stList_construct();
    ChainNode **nodes_by_idx = st_calloc((size_t) n, sizeof(ChainNode *));

    for (int64_t i = 0; i < n; i++) {
        TaffyAln *aln = &alns[i];
        ChainNode *cn = st_calloc(1, sizeof(ChainNode));
        cn->aln    = aln;
        cn->score  = aln->score;
        cn->pChain = NULL;
        nodes_by_idx[i] = cn;

        stSortedSetIterator *it = predecessor_iter(active, cn);
        ChainNode *pcn;
        while ((pcn = stSortedSet_getPrevious(it)) != NULL) {
            TaffyAln *p = pcn->aln;
            /* names + strand must match -- but in this partition the
             * q_name + strand are fixed, so just check t_name. */
            if (strcmp(aln->t_name ? aln->t_name : "",
                       p->t_name   ? p->t_name   : "") != 0) {
                /* sorted-set key starts with t_name; once we hit a
                 * different t_name, no earlier predecessor matches. */
                break;
            }
            if (aln->q_start < p->q_end) {
                /* query overlap on this predecessor; skip but keep scanning */
                continue;
            }
            int64_t q_gap = aln->q_start - p->q_end;
            if (q_gap > max_gap_length) {
                /* this predecessor can never be extended further along q;
                 * retire it from the active set. */
                stList_append(to_remove, pcn);
                continue;
            }
            if (aln->t_start < p->t_end) {
                continue;
            }
            int64_t t_gap = aln->t_start - p->t_end;
            if (t_gap > max_gap_length) {
                /* sorted-set key is t_end-ascending; subsequent (earlier
                 * t_end) predecessors will have even bigger t_gap. */
                break;
            }
            int64_t g = gap_cost(q_gap, t_gap, gap_cost_params);
            int64_t cand = aln->score + pcn->score - g;
            /* Match paffy's two guards: gap_cost must be less than the
             * new aln's own score (else starting fresh is at least as
             * good), AND the candidate must beat the current best. */
            if (g < aln->score && cand > cn->score) {
                cn->score  = cand;
                cn->pChain = pcn;
            }
        }
        stSortedSet_destructIterator(it);

        stSortedSet_insert(active, cn);
        stSortedSet_insert(chains, cn);

        while (stList_length(to_remove) > 0) {
            stSortedSet_remove(active, stList_pop(to_remove));
        }
    }

    /* Pop chains highest-score-first, assigning chain IDs.  When a
     * link is already claimed (its pChain is gone from `chains`), cut
     * the current chain there. */
    while (stSortedSet_size(chains) > 0) {
        ChainNode *root = stSortedSet_remove(chains, stSortedSet_getLast(chains));
        int64_t cid = next_chain_id++;
        int64_t bp = 0, n_alns_in_chain = 0;
        int64_t total_score = root->score;

        ChainNode *walker = root;
        while (1) {
            /* assign id to this aln; recover index from the aln pointer */
            int64_t idx = walker->aln - alns;
            assert(idx >= 0 && idx < n);
            chain_id_out[idx] = cid;
            bp += walker->aln->q_end - walker->aln->q_start;
            n_alns_in_chain++;
            if (walker->pChain == NULL) break;
            if (stSortedSet_search(chains, walker->pChain) == NULL) {
                /* predecessor already taken by a higher-scoring chain;
                 * cut here. */
                walker->pChain = NULL;
                break;
            }
            walker = walker->pChain;
            ChainNode *taken = stSortedSet_remove(chains, walker);
            assert(taken == walker);
        }

        /* grow chains_out */
        TaffyChainInfo *arr = *chains_out;
        int64_t k = (*n_chains_out)++;
        if ((k & (k + 1)) == 0) {  /* k+1 is a power of 2 -> realloc */
            int64_t new_cap = (k + 1) * 2;
            arr = st_realloc(arr, (size_t) new_cap * sizeof(TaffyChainInfo));
            *chains_out = arr;
        }
        arr[k].id          = cid;
        arr[k].total_score = total_score;
        arr[k].total_bp    = bp;
        arr[k].n_alns      = n_alns_in_chain;
    }

    /* Free every ChainNode via the full index.  ID assignment removed them
     * all from `chains` (and `active` holds borrowed pointers with a NULL
     * destructor), so neither set owns them -- nodes_by_idx[i] is the only
     * surviving handle to each. */
    for (int64_t i = 0; i < n; i++) free(nodes_by_idx[i]);
    free(nodes_by_idx);
    stList_setDestructor(to_remove, NULL);
    stList_destruct(to_remove);
    stSortedSet_destruct(active);
    stSortedSet_destruct(chains);
    return next_chain_id;
}

/* ------------------------------------------------------------------ */
/* Public entry                                                        */
/* ------------------------------------------------------------------ */

/* Comparator for the final per-chain ranking. */
static int chain_info_cmp_by_score_desc(const void *a, const void *b) {
    const TaffyChainInfo *x = a, *y = b;
    if (x->total_score != y->total_score)
        return x->total_score > y->total_score ? -1 : 1;
    return intcmp(x->id, y->id);
}

/* Mirror q coords for reverse-strand alns: q_start' = -q_end,
 * q_end' = -q_start.  Sweeping in q_start order then visits reverse-
 * strand alns in the natural "right-to-left on the forward axis"
 * order, which matches paffy's invert_query_strand trick. */
static void mirror_q(TaffyAln *a) {
    int64_t s = a->q_start;
    a->q_start = -a->q_end;
    a->q_end   = -s;
}

void taffy_chain(TaffyAln *alns, int64_t n,
                 int64_t (*gap_cost)(int64_t, int64_t, void *), void *gap_cost_params,
                 int64_t max_gap_length,
                 int64_t *chain_id,
                 TaffyChainInfo **chains_out, int64_t *n_chains_out) {
    *chains_out = NULL;
    *n_chains_out = 0;
    if (n == 0) return;

    /* Mirror reverse-strand q coords. */
    for (int64_t i = 0; i < n; i++)
        if (alns[i].strand < 0) mirror_q(&alns[i]);

    /* Sort all alns by (q_name, strand, q_start) -- so same-partition
     * alns are contiguous and each sub-sweep is one pass. */
    qsort(alns, (size_t) n, sizeof(TaffyAln), aln_cmp_by_q_loc);

    /* Sweep each partition. */
    int64_t next_id = 1;
    int64_t i = 0;
    while (i < n) {
        int64_t j = i + 1;
        const char *qn = alns[i].q_name ? alns[i].q_name : "";
        int strand = alns[i].strand;
        while (j < n
               && strcmp(alns[j].q_name ? alns[j].q_name : "", qn) == 0
               && alns[j].strand == strand) j++;
        next_id = chain_partition(alns + i, j - i, gap_cost, gap_cost_params,
                                  max_gap_length,
                                  chain_id + i,
                                  chains_out, n_chains_out, next_id);
        i = j;
    }

    /* Restore forward-strand q coords on reverse-strand alns. */
    for (int64_t i = 0; i < n; i++)
        if (alns[i].strand < 0) mirror_q(&alns[i]);

    /* Sort chains_out by total_score descending. */
    qsort(*chains_out, (size_t) *n_chains_out, sizeof(TaffyChainInfo),
          chain_info_cmp_by_score_desc);
}

/* ------------------------------------------------------------------ */
/* Post-chain collinear merge.  See chain.h for the rule + contract.  */
/* ------------------------------------------------------------------ */

typedef struct {
    int64_t chain_id;
    int64_t q_start;
    int64_t idx;
} TaffyMergeKey;

static int taffy_merge_key_cmp(const void *a, const void *b) {
    const TaffyMergeKey *p = (const TaffyMergeKey *) a;
    const TaffyMergeKey *q = (const TaffyMergeKey *) b;
    if (p->chain_id != q->chain_id) return p->chain_id < q->chain_id ? -1 : 1;
    if (p->q_start  != q->q_start)  return p->q_start  < q->q_start  ? -1 : 1;
    return 0;
}

int64_t taffy_chain_merge_collinear(TaffyAln *alns, int64_t n,
                                    const int64_t *chain_id,
                                    void (*on_merge)(TaffyAln *kept,
                                                     TaffyAln *absorbed,
                                                     void *user),
                                    void *user) {
    if (n <= 1) return n;

    /* Group by chain_id, then walk each chain's alns in q_start order.
     * Within a chain, taffy_chain's invariant already guarantees
     * q_start order -- so we could skip the per-chain re-sort -- but
     * we don't rely on that here: alns[] may have been mutated by an
     * earlier caller pass.  The combined sort is O(n log n). */
    TaffyMergeKey *keys = (TaffyMergeKey *) st_malloc((size_t) n * sizeof(TaffyMergeKey));
    for (int64_t i = 0; i < n; i++) {
        keys[i].chain_id = chain_id[i];
        keys[i].q_start  = alns[i].q_start;
        keys[i].idx      = i;
    }
    qsort(keys, (size_t) n, sizeof(TaffyMergeKey), taffy_merge_key_cmp);

    int64_t n_merges = 0;
    int64_t i = 0;
    while (i < n) {
        /* Find the end of this chain_id group. */
        int64_t j = i + 1;
        while (j < n && keys[j].chain_id == keys[i].chain_id) j++;

        /* Walk [i, j) in q_start order, extending `kept` through
         * consecutive abutting neighbors. */
        int64_t kept_idx = keys[i].idx;
        for (int64_t p = i + 1; p < j; p++) {
            int64_t cur_idx = keys[p].idx;
            TaffyAln *prev = &alns[kept_idx];
            TaffyAln *cur  = &alns[cur_idx];
            /* Within one chain, q_name / t_name / strand are
             * guaranteed identical by taffy_chain (chain.h invariant);
             * only the coord abut needs checking. */
            int abut_q = (prev->q_end == cur->q_start);
            int abut_t = (prev->strand > 0)
                       ? (prev->t_end  == cur->t_start)
                       : (cur->t_end   == prev->t_start);
            if (abut_q && abut_t) {
                /* Grow kept's extents to cover both.  + strand: t_end
                 * advances; - strand: t_start slides down (the target
                 * range still spans [t_start, t_end) post-merge). */
                prev->q_end = cur->q_end;
                if (prev->strand > 0) prev->t_end   = cur->t_end;
                else                  prev->t_start = cur->t_start;
                prev->score += cur->score;
                if (on_merge) on_merge(prev, cur, user);
                n_merges++;
                /* kept_idx stays the same: the merged prev may chain
                 * with the next p, too (consecutive merges). */
            } else {
                /* No merge: cur becomes the new kept for subsequent
                 * abut checks within this chain. */
                kept_idx = cur_idx;
            }
        }
        i = j;
    }
    free(keys);
    return n - n_merges;
}

typedef struct { int64_t lo, hi; } QIv;
static int qiv_cmp_lo(const void *a, const void *b) {
    int64_t la = ((const QIv *) a)->lo, lb = ((const QIv *) b)->lo;
    return (la > lb) - (la < lb);
}

/* Overlap-aware paralogy filter -- see chain.h.
 *
 * Implementation:
 *   - CSR-style chain bucket built in one O(n) pass: bucket_off[cid..
 *     cid+1) is the contiguous index range in bucket_aln for chain cid.
 *     The global sort is ASCENDING forward-q for forward chains but
 *     DESCENDING forward-q for reverse chains (taffy_chain mirrors
 *     reverse-strand q-coords), so each chain's q-intervals are re-sorted
 *     ascending before the union -- see the per-chain merge below.
 *   - kept_iv is maintained sorted+disjoint via proper 2-list union
 *     (O(m+n) per accept) into a ping-pong buffer; both buffers grow
 *     geometrically on demand.
 *   - O(n log n) total in the common case.  Earlier "rescan all alns
 *     per chain" was O(n_chains * n) and intractable for chromosome
 *     lifts. */
void taffy_chain_overlap_frac_select(const TaffyAln *alns, int64_t n,
                                     const int64_t *chain_id,
                                     const TaffyChainInfo *chains_out,
                                     int64_t n_chains_out,
                                     int64_t max_id,
                                     double overlap_frac,
                                     int64_t cap,
                                     char *keep_chain) {
    if (n == 0 || n_chains_out == 0) return;

    /* CSR bucket of alns by chain id. */
    int64_t *bucket_off = st_calloc((size_t)(max_id + 2), sizeof(int64_t));
    int64_t *bucket_aln = st_malloc((size_t) n * sizeof(int64_t));
    for (int64_t i = 0; i < n; i++) bucket_off[chain_id[i] + 1]++;
    for (int64_t k = 1; k <= max_id + 1; k++) bucket_off[k] += bucket_off[k-1];
    {
        int64_t *cursor = st_malloc((size_t)(max_id + 2) * sizeof(int64_t));
        memcpy(cursor, bucket_off, (size_t)(max_id + 2) * sizeof(int64_t));
        for (int64_t i = 0; i < n; i++)
            bucket_aln[cursor[chain_id[i]]++] = i;
        free(cursor);
    }

    /* Ping-pong kept_iv buffers + growable candidate buffer. */
    int64_t cap_kept = 4096, cap_cand = 4096;
    QIv *kept_a   = st_malloc((size_t)cap_kept * sizeof(QIv));
    QIv *kept_b   = st_malloc((size_t)cap_kept * sizeof(QIv));
    QIv *cand_iv  = st_malloc((size_t)cap_cand * sizeof(QIv));
    QIv *kept_iv  = kept_a;
    QIv *kept_nxt = kept_b;
    int64_t n_kept_iv = 0;
    int64_t survivor_cap = (cap > 0 && cap < n_chains_out) ? cap : n_chains_out;
    int64_t n_kept_chains = 0;

    for (int64_t k = 0; k < n_chains_out; k++) {
        int64_t cid = chains_out[k].id;
        int64_t b_lo = bucket_off[cid];
        int64_t b_hi = bucket_off[cid + 1];
        if (b_lo == b_hi) continue;
        int64_t n_in = b_hi - b_lo;
        if (n_in > cap_cand) {
            while (cap_cand < n_in) cap_cand *= 2;
            cand_iv = st_realloc(cand_iv, (size_t)cap_cand * sizeof(QIv));
        }
        /* Gather this chain's q-intervals, then SORT ascending by lo
         * before merging.  A reverse-strand chain's alns sit in DESCENDING
         * forward-q order, so a single forward pass over raw array order
         * would fold every interval but the rightmost into one (extending
         * only .hi) -- undercounting cand_bp AND leaving cand_iv unsorted,
         * which corrupts the kept_iv 2-list union below (it requires both
         * inputs sorted+disjoint).  Forward chains are already sorted, so
         * the sort is a no-op for them. */
        int64_t nm = 0;
        for (int64_t z = b_lo; z < b_hi; z++) {
            cand_iv[nm].lo = alns[bucket_aln[z]].q_start;
            cand_iv[nm].hi = alns[bucket_aln[z]].q_end;
            nm++;
        }
        qsort(cand_iv, (size_t) nm, sizeof(QIv), qiv_cmp_lo);
        {   /* collapse overlapping/abutting intervals in place */
            int64_t w = 0;
            for (int64_t r = 0; r < nm; r++) {
                if (w > 0 && cand_iv[r].lo <= cand_iv[w-1].hi) {
                    if (cand_iv[r].hi > cand_iv[w-1].hi)
                        cand_iv[w-1].hi = cand_iv[r].hi;
                } else {
                    cand_iv[w++] = cand_iv[r];
                }
            }
            nm = w;
        }
        int64_t cand_bp = 0;
        for (int64_t i = 0; i < nm; i++) cand_bp += cand_iv[i].hi - cand_iv[i].lo;
        if (cand_bp == 0) continue;
        /* Overlap with kept_iv (both sorted+disjoint). */
        int64_t ovr_bp = 0;
        {
            int64_t i = 0, j = 0;
            while (i < n_kept_iv && j < nm) {
                int64_t lo = kept_iv[i].lo > cand_iv[j].lo
                                ? kept_iv[i].lo : cand_iv[j].lo;
                int64_t hi = kept_iv[i].hi < cand_iv[j].hi
                                ? kept_iv[i].hi : cand_iv[j].hi;
                if (hi > lo) ovr_bp += hi - lo;
                if (kept_iv[i].hi < cand_iv[j].hi) i++;
                else j++;
            }
        }
        if ((double) ovr_bp > overlap_frac * (double) cand_bp)
            continue;  /* drop */
        /* Accept: merge cand_iv into kept_iv via 2-list union. */
        int64_t need = n_kept_iv + nm;
        if (need > cap_kept) {
            int kept_is_a = (kept_iv == kept_a);
            while (cap_kept < need) cap_kept *= 2;
            kept_a = st_realloc(kept_a, (size_t)cap_kept * sizeof(QIv));
            kept_b = st_realloc(kept_b, (size_t)cap_kept * sizeof(QIv));
            kept_iv  = kept_is_a ? kept_a : kept_b;
            kept_nxt = kept_is_a ? kept_b : kept_a;
        }
        int64_t i = 0, j = 0, k_out = 0;
        while (i < n_kept_iv || j < nm) {
            QIv pick;
            if (j >= nm || (i < n_kept_iv && kept_iv[i].lo <= cand_iv[j].lo)) {
                pick = kept_iv[i++];
            } else {
                pick = cand_iv[j++];
            }
            if (k_out > 0 && pick.lo <= kept_nxt[k_out-1].hi) {
                if (pick.hi > kept_nxt[k_out-1].hi)
                    kept_nxt[k_out-1].hi = pick.hi;
            } else {
                kept_nxt[k_out++] = pick;
            }
        }
        n_kept_iv = k_out;
        QIv *swap = kept_iv; kept_iv = kept_nxt; kept_nxt = swap;

        keep_chain[cid] = 1;
        n_kept_chains++;
        if (n_kept_chains >= survivor_cap) break;
    }
    free(kept_a); free(kept_b); free(cand_iv);
    free(bucket_off); free(bucket_aln);
}
