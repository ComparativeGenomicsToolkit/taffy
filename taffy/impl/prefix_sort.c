#include "taf.h"
#include "tai.h"
#include "sonLib.h"
#include <stdio.h>
#include <ctype.h>

/*
 * Functions for sorting alignment rows by prefixes of their sequence name
 */

void sequence_prefix_destruct(Sequence_Prefix *sequence_prefix) {
    free(sequence_prefix->prefix);
    free(sequence_prefix);
}

Sequence_Prefix *sequence_prefix_construct(char *prefix, int64_t index) {
    Sequence_Prefix *sequence_prefix = st_calloc(1, sizeof(Sequence_Prefix));
    sequence_prefix->prefix = prefix;
    sequence_prefix->prefix_length = strlen(prefix);
    if(sequence_prefix->prefix_length == 0) {
        st_errAbort("Found an empty sequence prefix");
    }
    sequence_prefix->index = index;
    return sequence_prefix;
}

int sequence_prefix_cmp_fn(Sequence_Prefix *p1, Sequence_Prefix *p2) {
    return strcmp(p1->prefix, p2->prefix);
}

stList *sequence_prefix_load(FILE *sort_fh) {
    stList *prefixes_to_sort_by = stList_construct3(0, (void (*)(void *))sequence_prefix_destruct);
    int64_t index = 0;
    char *line;
    while((line = stFile_getLineFromFile(sort_fh)) != NULL) {
        stList *tokens = stString_split(line);
        if(stList_length(tokens) != 1) {
            st_errAbort("Expected exactly one string in sort file on line: %s", line);
        }
        Sequence_Prefix *sequence_prefix = sequence_prefix_construct(stList_pop(tokens), index++);
        stList_append(prefixes_to_sort_by, sequence_prefix);
        // Clean up
        free(line);
        stList_destruct(tokens);
    }
    stList_sort(prefixes_to_sort_by, (int (*)(const void *, const void *))sequence_prefix_cmp_fn);
    return prefixes_to_sort_by;
}

static int get_closest_prefix_cmp_fn(char *sequence_name, Sequence_Prefix *sp) {
    int64_t i = strcmp(sequence_name, sp->prefix);
    if(i > 0) { // If sequence_name is lexicographically larger than sequence_prefix could
        // be a prefix (can not be a prefix is i < 0)
        for(int64_t j=0; j<sp->prefix_length; j++) {
            if(sequence_name[j] != sp->prefix[j]) {
                return 1;
            }
        }
        return 0;
    }
    return i;
}

int64_t alignment_row_get_closest_sequence_prefix(Alignment_Row *row, stList *prefixes_to_sort_by) {
    // Binary search the sequence name
    Sequence_Prefix *sp = stList_binarySearch(prefixes_to_sort_by, row->sequence_name,
                                              (int (*)(const void *a, const void *b))get_closest_prefix_cmp_fn);
    if(sp == NULL) {
        st_logDebug("Did not find a valid prefix to match: %s\n", row->sequence_name);
    }
    return sp != NULL ? sp->index : -1; // Sequences that don't have a match will appear first in the sort
}

int alignment_sequence_prefix_cmp_fn(Alignment_Row *a1, Alignment_Row *a2,
                                     stList *prefixes_to_sort_by) {
    int i = alignment_row_get_closest_sequence_prefix(a1, prefixes_to_sort_by);
    int j = alignment_row_get_closest_sequence_prefix(a2, prefixes_to_sort_by);
    return i < j ? -1 : (i > j ? 1 : strcmp(a1->sequence_name, a2->sequence_name));
}

void alignment_sort_the_rows(Alignment *p_alignment, Alignment *alignment, stList *prefixes_to_sort_by, bool ignore_first_row) {
    // Get the rows
    stList *rows = alignment_get_rows_in_a_list(ignore_first_row && alignment->row ? alignment->row->n_row : alignment->row);
    assert(stList_length(rows) == (ignore_first_row ? alignment->row_number -1 : alignment->row_number)); // Quick sanity check

    // Sort the rows by the prefix ordering
    stList_sort2(rows, (int (*)(const void *, const void *, void *))alignment_sequence_prefix_cmp_fn, prefixes_to_sort_by);

    // Add back the first row if ignored
    if(ignore_first_row && alignment->row) {
        // hack for now cos no insert method in stList!
        stList_reverse(rows);
        stList_append(rows, alignment->row);
        stList_reverse(rows);
        assert(stList_get(rows, 0) == alignment->row);
    }
    assert(stList_length(rows) == alignment->row_number); // One more sanity check

    // Re-connect the rows
    alignment_set_rows(alignment, rows);

    // Reset the alignment of the rows with the prior row
    if(p_alignment != NULL) {
        alignment_link_adjacent(p_alignment, alignment, 1);
    }
}

static void remove_rows(Alignment *alignment, int (*delete_row)(Alignment_Row *, void *),
                        void *extra_arg, bool ignore_first_row) {
    Alignment_Row *row = alignment->row, **p_row = &(alignment->row);
    if(ignore_first_row && row != NULL) {  // Keep the first row
        p_row = &(row->n_row);  // Update pointers
        row = row->n_row;
    }
    while(row != NULL) {
        if(delete_row(row, extra_arg)) { //stSet_search(rows_to_delete, row) != -1) { // Filter the row
            alignment->row_number--; // Reduce the number of row
            assert(alignment->row_number >= 0);
            *p_row = row->n_row; // Update the previous link to point at the row after the current row
            alignment_row_destruct(row); // Clean up the row
            row = (*p_row); // Set row to point at the next row
        }
        else { // Keep the row, do nothing
            p_row = &(row->n_row);  // Update pointers
            row = row->n_row;
        }
    }
}

static int alignment_filter_fn(Alignment_Row *row, stList *prefixes_to_filter_by) {
    return alignment_row_get_closest_sequence_prefix(row, prefixes_to_filter_by) != -1;
}

void alignment_filter_the_rows(Alignment *alignment, stList *prefixes_to_filter_by, bool ignore_first_row) {
    remove_rows(alignment, (int (*)(Alignment_Row *, void *))alignment_filter_fn,
                       prefixes_to_filter_by, ignore_first_row);
}

void alignment_show_only_lineage_differences(Alignment *alignment, char mask_char, stList *sequence_prefixes, stList *tree_nodes) {
    // First create map of tree nodes to bases
    stHash *tree_nodes_to_bases = stHash_construct2(NULL, (void (*)(void *))stList_destruct);
    Alignment_Row *row = alignment->row;
    while(row != NULL) { // For each row
        // Get corresponding tree node using sequence prefixes
        int64_t i = alignment_row_get_closest_sequence_prefix(row, sequence_prefixes);
        if(i != -1) { // We found it in the tree
            // Add to the tree_nodes_to_bases map
            stTree *node = stList_get(tree_nodes, i);
            stList *ancestor_sequences = stHash_search(tree_nodes_to_bases, node);
            if(ancestor_sequences == NULL) {
                ancestor_sequences = stList_construct3(0, free);
                stHash_insert(tree_nodes_to_bases, node, ancestor_sequences);
            }
            stList_append(ancestor_sequences, stString_copy(row->bases));
        }
        else {
            st_logDebug("Alignment row sequence not found in tree: %s\n", row->sequence_name);
        }
        row = row->n_row;
    }

    // Now identify mutations
    row = alignment->row;
    while(row != NULL) { // For each row
        //  Identify node in tree using prefix search
        int64_t i = alignment_row_get_closest_sequence_prefix(row, sequence_prefixes);
        if(i != -1) { // We found it in the tree
            stTree *node = stList_get(tree_nodes, i);
            stTree *ancestor = stTree_getParent(node);
            if(ancestor != NULL) { // If we're not at the root of the tree - otherwise we must report the base
                stList *ancestor_sequences = stHash_search(tree_nodes_to_bases, ancestor);
                for(int64_t j=0; j<alignment->column_number; j++) { // For each alignment column
                    char base = row->bases[j];
                    if(base != '-') { // If not a gap base
                        for (int64_t k = 0; k < stList_length(ancestor_sequences); k++) { // For each ancestor base
                            char *ancestor_sequence = stList_get(ancestor_sequences, k);
                            if(toupper(base) == toupper(ancestor_sequence[j])) { // If identical to ancestor base
                                row->bases[j] = mask_char; // Switch to a star character
                                break;
                            }
                        }
                    }
                }
            }
        } // If not found do nothing, as we log missing sequences in the loop above
        row = row->n_row;
    }

    // Clean up
    stHash_destruct(tree_nodes_to_bases);
}

/*
 * Functions to pad an alignment block with an extra dummy row for each missing sequences - helpful for
 * making normalized alignments
 */

static int get_closest_row_cmp_fn(Sequence_Prefix *p, Alignment_Row *r) {
    int64_t i = strcmp(p->prefix, r->sequence_name);
    if(i < 0) { // If prefix_name is lexicographically smaller than row's sequence could
        // be a prefix (can not be a prefix is i > 0)
        for(int64_t j=0; j<p->prefix_length; j++) {
            if(r->sequence_name[j] != p->prefix[j]) {
                return -1;
            }
        }
        return 0;
    }
    return i;
}

static int alignment_row_cmp_fn(Alignment_Row *r1, Alignment_Row *r2) {
    return strcmp(r1->sequence_name, r2->sequence_name);
}

void alignment_pad_the_rows(Alignment *p_alignment, Alignment *alignment, stList *sequence_prefixes) {
    // Get the rows in a list
    stList *rows = alignment_get_rows_in_a_list(alignment->row);
    assert(stList_length(rows) == (alignment->row_number)); // Quick sanity check

    // Sort the rows by name
    stList_sort(rows, (int (*)(const void *, const void *))alignment_row_cmp_fn);

    // Get the pointer to the last row of the alignment so we can add rows
    Alignment_Row **p_r = &(alignment->row);
    while(*p_r != NULL) {
        p_r = &((*p_r)->n_row);
    }

    // For each sequence prefix
    for(int64_t i=0; i<stList_length(sequence_prefixes); i++) {
        Sequence_Prefix *sp = stList_get(sequence_prefixes, i);

        // Check if there is a corresponding sequence name using binary search
        Alignment_Row *r = stList_binarySearch(rows, sp,(int (*)(const void *a, const void *b))get_closest_row_cmp_fn);

        if(r == NULL) { // If there isn't a corresponding row, add one to the alignment at the end setting the coordinates to zero
            r = st_calloc(1, sizeof(Alignment_Row));
            alignment->row_number++; // Increment the row number
            r->sequence_name = stString_copy(sp->prefix);
            r->bases = st_calloc(alignment->column_number+1, sizeof(char));
            for(int64_t j=0; j<alignment->column_number; j++) {
                r->bases[j] = '-';
            }
            r->bases[alignment->column_number] = '\0';
            r->strand = 1;
            *p_r = r;
            p_r = &(r->n_row);
        }
    }

    // Clean up
    stList_destruct(rows);

    // Reset the alignment of the rows with the prior row
    if(p_alignment != NULL) {
        alignment_link_adjacent(p_alignment, alignment, 1);
    }
}

/*
 * Functions to filter duplicate sequence rows
 */

int alignment_filter_rows_fn(Alignment_Row *row, stSet *rows_to_delete) {
    return stSet_search(rows_to_delete, row) != NULL;
}

void alignment_filter_duplicate_rows(Alignment *alignment, stList *prefixes_to_match_on, bool ignore_first_row) {
    // Create a map from prefixes to rows
    stHash *prefixes_to_rows = stHash_construct2(NULL, (void (*)(void *))stList_destruct);
    stSet *rows_to_delete = stSet_construct();
    Alignment_Row *r = alignment->row;
    while(r != NULL) {
        Sequence_Prefix *sp = stList_binarySearch(prefixes_to_match_on, r->sequence_name,
                                                  (int (*)(const void *a, const void *b))get_closest_prefix_cmp_fn);
        if(sp != NULL) {
            stList *l = stHash_search(prefixes_to_rows, sp);
            if (!l) {
                l = stList_construct();
                stHash_insert(prefixes_to_rows, sp, l);
            }
            stList_append(l, r);
        }
        r = r->n_row;
    }

    // For each sequence prefix
    for(int64_t i=0; i<stList_length(prefixes_to_match_on); i++) {
        Sequence_Prefix *sp = stList_get(prefixes_to_match_on, i);
        stList *l = stHash_search(prefixes_to_rows, sp);

        // Where there is a sequence prefix with multiple rows
        if(l != NULL && stList_length(l) > 1) {

            // TODO: Add option to more intelligently rank rows....

            // Add weakest rows to the set to delete
            for(int64_t j=1; j<stList_length(l); j++) {
                stSet_insert(rows_to_delete, stList_get(l, j));
            }
        }
    }

    // Now delete the rows
    remove_rows(alignment, (int (*)(Alignment_Row *, void *))alignment_filter_rows_fn,
                rows_to_delete, ignore_first_row);

    // Cleanup
    stSet_destruct(rows_to_delete);
    stHash_destruct(prefixes_to_rows);
}
