#include "taf.h"
#include "ond.h"
#include "sonLib.h"

#define ANSI_COLOR_RED     "\x1b[41m"
#define ANSI_COLOR_GREEN   "\x1b[42m"
#define ANSI_COLOR_YELLOW  "\x1b[43m"
#define ANSI_COLOR_BLUE    "\x1b[44m"
#define ANSI_COLOR_MAGENTA "\x1b[45m"
#define ANSI_COLOR_CYAN    "\x1b[46m"
#define ANSI_COLOR_RESET   "\x1b[0m"

char *color_base_char(char base) {
    switch(base) {
        case 'A':
        case 'a':
            return stString_print(ANSI_COLOR_RED "%c" ANSI_COLOR_RESET, base);
        case 'C':
        case 'c':
            return stString_print(ANSI_COLOR_BLUE "%c" ANSI_COLOR_RESET, base);
        case 'G':
        case 'g':
            return stString_print(ANSI_COLOR_YELLOW "%c" ANSI_COLOR_RESET, base);
        case 'T':
        case 't':
            return stString_print(ANSI_COLOR_GREEN "%c" ANSI_COLOR_RESET, base);
        case '*':
            return stString_print(ANSI_COLOR_MAGENTA "%c" ANSI_COLOR_RESET, base);
        case '-':
            return stString_print(ANSI_COLOR_CYAN "%c" ANSI_COLOR_RESET, base);
        default:
            return stString_print("%c", base);
    }
}

char *color_base_string(char *bases, int64_t length) {
    stList *colored_base_strings = stList_construct3(0, free);
    for(int64_t i=0; i<length; i++) {
        char *base_string = color_base_char(bases[i]);
        stList_append(colored_base_strings, base_string);
    }
    char *colored_bases = stString_join2("", colored_base_strings);
    stList_destruct(colored_base_strings);
    return colored_bases;
}

void tag_destruct(Tag *tag) {
    while(tag != NULL) {
        Tag *p_tag = tag;
        tag = tag->n_tag;
        free(p_tag);
    }
}

Tag *tag_find(Tag *tag, char *key) {
    while(tag != NULL) {
        if(strcmp(tag->key, key) == 0) {
            return tag;
        }
        tag = tag->n_tag;
    }
    return NULL;
}

Tag *tag_remove(Tag *first_tag, char *key) {
    if(strcmp(first_tag->key, key) == 0) { // If first tag is one to remove
        return first_tag->n_tag;
    }
    // Tag to remove is not the first link
    Tag *tag = first_tag;
    while(tag->n_tag != NULL) {
        if(strcmp(tag->n_tag->key, key) == 0) {
            Tag *t = tag->n_tag;
            tag->n_tag = tag->n_tag->n_tag; // Remove the link
            free(t);
            break;
        }
        tag = tag->n_tag;
    }
    return first_tag;
}

Tag *tag_construct(char *key, char *value, Tag *n_tag) {
    Tag *tag = st_calloc(1, sizeof(Tag));
    tag->key = stString_copy(key);
    tag->value = stString_copy(value);
    tag->n_tag = n_tag;
    return tag;
}

Tag *tag_parse(char *tag_string, char *delimiter, Tag *p_tag) {
    stList *tag_tokens = stString_splitByString(tag_string, delimiter);
    if (stList_length(tag_tokens) != 2) {
        st_errAbort("Tag not separated by '%s' character: %s\n", delimiter, tag_string);
    }
    Tag *tag = st_calloc(1, sizeof(Tag));
    if(p_tag != NULL) { // If the p_tag is not null, set its n_tag to point at tag
        p_tag->n_tag = tag;
    }
    tag->key = stList_get(tag_tokens, 0);
    tag->value = stList_get(tag_tokens, 1);
    stList_setDestructor(tag_tokens, NULL);
    stList_destruct(tag_tokens);
    return tag;
}

void alignment_row_destruct(Alignment_Row *row) {
    if(row->l_row != NULL) { // If there is a preceding, left row then unlink it
        assert(row->l_row->r_row == row);
        row->l_row->r_row = NULL;
    }
    if(row->r_row != NULL) { // If there is a proceeding, right row then unlink it
        assert(row->r_row->l_row == row);
        row->r_row->l_row = NULL;
    }
    if(row->bases != NULL) {
        free(row->bases);
    }
    if(row->sequence_name != NULL) {
        free(row->sequence_name);
    }
    if(row->left_gap_sequence != NULL) {
        free(row->left_gap_sequence);
    }
    free(row);
}

void alignment_destruct(Alignment *alignment, bool cleanup_rows) {
    Alignment_Row *row = alignment->row;
    if(cleanup_rows) {
        while (row != NULL) {
            Alignment_Row *r = row;
            row = row->n_row;
            alignment_row_destruct(r);
        }
    }
    assert(alignment->column_tags != NULL);  // Clean up column tags
    for(int64_t i=0; i<alignment->column_number; i++) {
        tag_destruct(alignment->column_tags[i]);
    }
    free(alignment->column_tags);
    free(alignment);
}

stList *alignment_get_rows_in_a_list(Alignment_Row *row) {
    stList *l = stList_construct();
    while(row != NULL) {
        stList_append(l, row);
        row = row->n_row;
    }
    return l;
}

void alignment_set_rows(Alignment *alignment, stList *rows) {
    alignment->row_number = stList_length(rows);
    if(stList_length(rows) > 0) {
        alignment->row = stList_get(rows, 0);
        for (int64_t i = 1; i < stList_length(rows); i++) {
            Alignment_Row *p_row = stList_get(rows, i-1);
            Alignment_Row *row = stList_get(rows, i);
            p_row->n_row = row;
            row->n_row = NULL;
        }
    }
}

bool alignment_row_is_predecessor(Alignment_Row *left_row, Alignment_Row *right_row) {
    // Do the rows match
    return strcmp(left_row->sequence_name, right_row->sequence_name) == 0 && left_row->strand == right_row->strand &&
            left_row->start + left_row->length <= right_row->start;
}

bool alignment_row_is_predecessor_2(Alignment_Row **left_row, Alignment_Row **right_row) {
    // Do the rows match - this one is needed to work with the OND aligner which compares pointers to the objects being
    // compared
    return alignment_row_is_predecessor(left_row[0], right_row[0]);
}

char *alignment_row_to_string(Alignment_Row *row) {
    return stString_print("%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%" PRIi64 "\t%s",
                          row->sequence_name, row->start, row->length,
                          row->strand ? "+" : "-", row->sequence_length, row->bases);
}

void alignment_link_adjacent(Alignment *left_alignment, Alignment *right_alignment, bool allow_row_substitutions) {
    stList *left_rows = alignment_get_rows_in_a_list(left_alignment->row);
    stList *right_rows = alignment_get_rows_in_a_list(right_alignment->row);
    // get the alignment of the rows
    WFA *wfa = WFA_construct(stList_getBackingArray(left_rows), stList_getBackingArray(right_rows),
                             stList_length(left_rows), stList_length(right_rows),
                             sizeof(void *), (bool (*)(void *, void *))alignment_row_is_predecessor_2, 1,
                             allow_row_substitutions ? 1 : 100000000); // Use unit gap and mismatch costs for the diff
                             // unless we disallow substitutions, in which case use an arbitrarily large mismatch cost
    int64_t aligned_rows[stList_length(left_rows)];
    WFA_get_alignment(wfa, aligned_rows);
    // Remove any previous links
    Alignment_Row *row = left_alignment->row;
    while(row != NULL) {
        row->r_row = NULL;
        // Clean up any gap sequences because once unlinked we can not guarantee continuity
        if(row->left_gap_sequence != NULL) {
            free(row->left_gap_sequence);
            row->left_gap_sequence = NULL;
        }
        row = row->n_row;
    }
    row = right_alignment->row;
    while(row != NULL) {
        row->l_row = NULL;
        row = row->n_row;
    }
    // connect up the rows according to the alignment
    for(int64_t i=0; i<stList_length(left_rows); i++) {
        if(aligned_rows[i] != -1) {
            Alignment_Row *left_row = stList_get(left_rows, i);
            Alignment_Row *right_row = stList_get(right_rows, aligned_rows[i]);
            left_row->r_row = right_row;
            right_row->l_row = left_row;
            if(!allow_row_substitutions) {
                assert(alignment_row_is_predecessor(left_row, right_row));
            }
        }
    }
    // clean up
    stList_destruct(left_rows);
    stList_destruct(right_rows);
    WFA_destruct(wfa);
}

int64_t alignment_length(Alignment *alignment) {
    return alignment->column_number;
}

int64_t alignment_total_gap_length(Alignment *left_alignment) {
    Alignment_Row *l_row = left_alignment->row;
    int64_t total_interstitial_gap_length = 0;
    while(l_row != NULL) {
        if(l_row->r_row != NULL && alignment_row_is_predecessor(l_row, l_row->r_row)) {
            int64_t i = l_row->r_row->start - (l_row->start + l_row->length);
            if (i > total_interstitial_gap_length) {
                total_interstitial_gap_length = i;
            }
        }
        l_row = l_row->n_row; // Move to the next left alignment row
    }
    return total_interstitial_gap_length;
}

char *alignment_to_string(Alignment *alignment) {
    stList *strings = stList_construct3(0, free);
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        stList_append(strings,alignment_row_to_string(row));
        row = row->n_row;
    }
    char *block_string = stString_join2("\n", strings);
    stList_destruct(strings);
    return block_string;
}

void alignment_row_mask_identical_bases(Alignment *alignment, Alignment_Row *ref, Alignment_Row *non_ref, char mask_char) {
    for(int64_t i=0; i<alignment->column_number; i++) {
        if(ref->bases[i] == non_ref->bases[i]) {
            non_ref->bases[i] = mask_char;
        }
    }
}

void alignment_mask_reference_bases(Alignment *alignment, char mask_char) {
    Alignment_Row *ref_row = alignment->row;
    if(ref_row) {
        Alignment_Row *non_ref_row = ref_row->n_row;
        while (non_ref_row != NULL) {
            alignment_row_mask_identical_bases(alignment, ref_row, non_ref_row, mask_char);
            non_ref_row = non_ref_row->n_row;
        }
    }
}

/*
 * Returns a sequence of tags from the tokens, starting at starting_token
 */
Tag *parse_tags(stList *tokens, int64_t starting_token, char *delimiter) {
    Tag *first_tag = NULL, *tag = NULL;
    if(starting_token < stList_length(tokens)) { // parse first tag pair
        tag = tag_parse(stList_get(tokens, starting_token), delimiter, NULL);
        first_tag = tag;
    }
    for(int64_t i=starting_token+1; i<stList_length(tokens); i++) {
        tag = tag_parse(stList_get(tokens, i), delimiter, tag);
    }
    return first_tag;
}

Tag *parse_header(stList *tokens, char *header_prefix, char *delimiter) {
    if(stList_length(tokens) == 0 || strcmp((char *)stList_get(tokens, 0), header_prefix) != 0) {
        st_errAbort("Header line does not start with %s\n", header_prefix);
    }
    return parse_tags(tokens, 1, delimiter);
}

void write_header(Tag *tag, LW *lw, char *header_prefix, char *delimiter, char *end) {
    LW_write(lw, "%s", header_prefix);
    while(tag != NULL) {
        LW_write(lw, " %s%s%s", tag->key, delimiter, tag->value);
        tag = tag->n_tag;
    }
    LW_write(lw, "%s", end);
}

int64_t alignment_number_of_common_rows(Alignment *left_alignment, Alignment *right_alignment) {
    // First un-link any rows that are substitutions as these can't be merged
    Alignment_Row *r_row = right_alignment->row;
    int64_t shared_rows = 0;
    while (r_row != NULL) {
        if (r_row->l_row != NULL && alignment_row_is_predecessor(r_row->l_row, r_row)) {
            shared_rows++;
        }
        r_row = r_row->n_row;
    }
    return shared_rows;
}

void alignment_get_column_in_buffer(Alignment *alignment, int64_t column_index, char *buffer) {
    assert(column_index >= 0);
    assert(column_index < alignment->column_number);
    Alignment_Row *row = alignment->row;
    for(int64_t i=0; i<alignment->row_number; i++) {
        buffer[i] = row->bases[column_index];
        row = row->n_row;
    }
    assert(row == NULL);
}

char *alignment_get_column(Alignment *alignment, int64_t column_index) {
    char *column_string = st_malloc(sizeof(char) * (alignment->row_number+1));
    column_string[alignment->row_number] = '\0';
    assert(column_index >= 0);
    assert(column_index < alignment->column_number);
    Alignment_Row *row = alignment->row;
    for(int64_t i=0; i<alignment->row_number; i++) {
        column_string[i] = row->bases[column_index];
        row = row->n_row;
    }
    assert(row == NULL);
    return column_string;
}
