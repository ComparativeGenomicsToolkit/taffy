/*
 * view_chain_dup_filter: see view_chain_dup_filter.h.
 *
 *  Released under the MIT license, see LICENSE.txt
 */

#include "view_chain_dup_filter.h"
#include "chain.h"
#include <string.h>

/* First-dot genome partition.  Matches the .tui's heuristic when no
 * g-record roster / hal phylogeny is available -- internally consistent
 * across all rows in the buffer is what matters here, not whether the
 * resulting "genome" exactly equals the cactus species name. */
static char *first_dot_genome(const char *seq_name) {
    const char *dot = strchr(seq_name, '.');
    if (dot == NULL) return stString_copy(seq_name);
    return stString_getSubString(seq_name, 0, dot - seq_name);
}

typedef struct {
    int64_t        block_idx;
    Alignment_Row *row;
} RowRef;

void view_chain_dup_filter(stList *blocks, double overlap_frac, int64_t cap) {
    if (overlap_frac < 0) return;            /* disabled */
    if (overlap_frac > 1) overlap_frac = 1;
    int64_t n_blocks = stList_length(blocks);
    if (n_blocks == 0) return;

    /* 1. Group all non-row-0 rows by genome (first-dot partition).
     *    Separately track WHICH genomes have any block with 2+ rows of
     *    that genome -- only those genomes are dup candidates.  A genome
     *    that appears <= once per block (e.g. a continuous syntenic
     *    target across many blocks) is NOT a dup candidate even if it
     *    spans the whole buffer. */
    stHash *genome_to_rows = stHash_construct3(
        stHash_stringKey, stHash_stringEqualKey, free,
        (void(*)(void*)) stList_destruct);
    stSet  *dup_genomes = stSet_construct3(
        stHash_stringKey, stHash_stringEqualKey, free);

    for (int64_t bi = 0; bi < n_blocks; bi++) {
        Alignment *aln = (Alignment*) stList_get(blocks, bi);
        if (aln == NULL || aln->row == NULL) continue;

        /* Per-block count for the dup-genome detector.  We use a fresh
         * temp hash so a genome appearing in many blocks-once doesn't
         * get marked as a dup. */
        stHash *block_counts = stHash_construct3(
            stHash_stringKey, stHash_stringEqualKey, free, free);

        Alignment_Row *r = aln->row->n_row;          /* skip row-0 */
        while (r != NULL) {
            char *gn = first_dot_genome(r->sequence_name);

            /* Add the row to the per-genome list.  If gn is already in
             * the hash, free our copy and reuse. */
            stList *L = stHash_search(genome_to_rows, gn);
            if (L == NULL) {
                L = stList_construct3(0, free);
                stHash_insert(genome_to_rows, gn, L);  /* hash owns gn */
            } else {
                free(gn);
            }
            RowRef *rr = st_malloc(sizeof(*rr));
            rr->block_idx = bi;
            rr->row       = r;
            stList_append(L, rr);

            /* Per-block count -- mark gn as a dup-candidate once it
             * hits 2 in this block.  String-keyed; dup_genomes owns
             * its own copy of the name (independent of genome_to_rows
             * key ownership). */
            char *gn_count_key = first_dot_genome(r->sequence_name);
            int64_t *cnt = stHash_search(block_counts, gn_count_key);
            if (cnt == NULL) {
                cnt = st_malloc(sizeof(int64_t));
                *cnt = 1;
                stHash_insert(block_counts, gn_count_key, cnt);
            } else {
                (*cnt)++;
                free(gn_count_key);
                if (*cnt == 2) {
                    char *dg = first_dot_genome(r->sequence_name);
                    if (stSet_search(dup_genomes, dg) == NULL) {
                        stSet_insert(dup_genomes, dg);   /* set owns */
                    } else {
                        free(dg);
                    }
                }
            }
            r = r->n_row;
        }
        stHash_destruct(block_counts);
    }

    /* If no genome ever had a per-block dup, there's nothing to filter. */
    if (stSet_size(dup_genomes) == 0) {
        stHash_destruct(genome_to_rows);
        stSet_destruct(dup_genomes);
        return;
    }

    /* 2. Per dup-candidate genome: chain its alns; top_n chain ids
     *    flag survivors.  Survivors are tracked by Alignment_Row pointer. */
    stSet *survivors = stSet_construct();
    stHashIterator *it = stHash_getIterator(genome_to_rows);
    char *gn;
    while ((gn = (char*)stHash_getNext(it)) != NULL) {
        stList *rows = stHash_search(genome_to_rows, gn);
        int64_t n = stList_length(rows);
        /* If gn isn't a dup-candidate, keep all its rows unchanged. */
        if (stSet_search(dup_genomes, gn) == NULL) {
            for (int64_t i = 0; i < n; i++) {
                RowRef *rr = stList_get(rows, i);
                stSet_insert(survivors, rr->row);
            }
            continue;
        }

        TaffyAln *alns = st_malloc(n * sizeof(TaffyAln));
        for (int64_t i = 0; i < n; i++) {
            RowRef *rr = stList_get(rows, i);
            Alignment_Row *rt = rr->row;        /* target row */
            int strand = rt->strand ? +1 : -1;
            /* Forward-strand coords on both axes; mirror target on '-'. */
            int64_t t_start, t_end;
            if (rt->strand) {
                t_start = rt->start;
                t_end   = rt->start + rt->length;
            } else {
                t_end   = rt->sequence_length - rt->start;
                t_start = t_end - rt->length;
            }
            /* q axis: block-order index so all alns of this target
             * genome land in ONE partition.  Using row-0's sequence
             * name would partition by ancestor, breaking chains across
             * blocks where the universal-MAF row-0 ancestor switches
             * (Anc1 -> Anc0 etc.) -- a chain-break artifact unrelated
             * to target synteny.  +1 q_end ensures non-zero q-extent. */
            alns[i].q_name  = "block_axis";
            alns[i].q_start = rr->block_idx;
            alns[i].q_end   = rr->block_idx + 1;
            alns[i].t_name  = rt->sequence_name;
            alns[i].t_start = t_start;
            alns[i].t_end   = t_end;
            alns[i].strand  = strand;
            alns[i].score   = rt->length;
            alns[i].user    = rr;
        }

        int64_t *chain_id = st_malloc(n * sizeof(int64_t));
        TaffyChainInfo *chains = NULL;
        int64_t n_chains = 0;
        TaffyChainCostParams cp = {
            TAFFY_CHAIN_DEFAULT_OPEN, TAFFY_CHAIN_DEFAULT_EXTEND };
        taffy_chain(alns, n,
                    taffy_chain_default_gap_cost, &cp,
                    TAFFY_CHAIN_DEFAULT_MAX_GAP,
                    chain_id, &chains, &n_chains);

        /* Run the shared overlap-frac selector: keep_chain[cid] flags
         * surviving chain ids (1-based, sized max_id + 1). */
        int64_t max_id = 0;
        for (int64_t k = 0; k < n_chains; k++)
            if (chains[k].id > max_id) max_id = chains[k].id;
        char *keep_chain = st_calloc((size_t)(max_id + 1), sizeof(char));
        taffy_chain_overlap_frac_select(alns, n, chain_id,
                                        chains, n_chains, max_id,
                                        overlap_frac, cap,
                                        keep_chain);
        for (int64_t i = 0; i < n; i++) {
            if (keep_chain[chain_id[i]]) {
                RowRef *rr = (RowRef*) alns[i].user;
                stSet_insert(survivors, rr->row);
            }
        }
        free(keep_chain);
        free(chain_id);
        free(chains);
        free(alns);
    }
    stHash_destructIterator(it);

    /* 3. Drop any non-row-0 row whose genome is a dup candidate AND row
     *    isn't a survivor.  Non-dup-candidate genomes are untouched.
     *    Mirrors the in-place unlink template in taf_view.c's
     *    filter_out_ancestor_rows. */
    for (int64_t bi = 0; bi < n_blocks; bi++) {
        Alignment *aln = (Alignment*) stList_get(blocks, bi);
        if (aln == NULL || aln->row == NULL || aln->row->n_row == NULL) continue;
        Alignment_Row **link = &aln->row->n_row;
        while (*link != NULL) {
            Alignment_Row *r = *link;
            char *gn = first_dot_genome(r->sequence_name);
            int is_dup_genome = (stSet_search(dup_genomes, gn) != NULL);
            free(gn);
            if (is_dup_genome && stSet_search(survivors, r) == NULL) {
                *link = r->n_row;
                r->n_row = NULL;
                if (r->l_row != NULL) r->l_row->r_row = NULL;
                r->l_row = NULL;
                alignment_row_destruct(r);
                aln->row_number--;
            } else {
                link = &r->n_row;
            }
        }
    }

    stSet_destruct(survivors);
    stHash_destruct(genome_to_rows);
    stSet_destruct(dup_genomes);
}
