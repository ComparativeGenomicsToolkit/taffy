//
// Created by Benedict Paten on 8/25/22.

#include "taf.h"
#include "sonLib.h"
#include "line_iterator.h"

Alignment *maf_read_block(LI *li) {
    while(1) {
        char *line = LI_get_next_line(li);
        if(line == NULL) {
            return NULL;
        }
        stList *tokens = stString_split(line);
        free(line);
        if(stList_length(tokens) == 0) {
            stList_destruct(tokens);
            continue;
        }
        if(strcmp(stList_get(tokens, 0), "a") == 0) { // If is an "a" line
            stList_destruct(tokens);
            Alignment *alignment = st_calloc(1, sizeof(Alignment));
            Alignment_Row **p_row = &(alignment->row);
            while(1) {
                line = LI_get_next_line(li);
                if(line == NULL) {
                    return alignment;
                }
                tokens = stString_split(line);
                free(line);
                if(stList_length(tokens) == 0) {
                    stList_destruct(tokens);
                    return alignment;
                }
                if(strcmp(stList_get(tokens, 0), "s") != 0) {
                    assert(strcmp(stList_get(tokens, 0), "i") == 0 || strcmp(stList_get(tokens, 0), "e") == 0); // Must be an "i" or "e" line, which we ignore
                    stList_destruct(tokens);
                    continue;
                }
                assert(strcmp(stList_get(tokens, 0), "s") == 0); // Must be an "s" line
                Alignment_Row *row = st_calloc(1, sizeof(Alignment_Row));
                alignment->row_number++;
                *p_row = row;
                p_row = &(row->n_row);
                row->sequence_name = stString_copy(stList_get(tokens, 1));
                row->start = atol(stList_get(tokens, 2));
                row->length = atol(stList_get(tokens, 3));
                assert(strcmp(stList_get(tokens, 4), "+") == 0 || strcmp(stList_get(tokens, 4), "-") == 0);
                row->strand = strcmp(stList_get(tokens, 4), "+") == 0;
                row->sequence_length = atol(stList_get(tokens, 5));
                row->bases = stString_copy(stList_get(tokens, 6));
                stList_destruct(tokens);
                if(alignment->row_number == 1) {
                    alignment->column_number = strlen(row->bases);
                    alignment->column_tags = st_calloc(alignment->column_number, sizeof(Tag *));
                }
                else {
                    assert(alignment->column_number == strlen(row->bases));
                }
            }
            return alignment;
        }
        else {
            assert(strcmp(stList_get(tokens, 0), "s") != 0); // Can not be an s line without a prior a line - we will ignore this line
            stList_destruct(tokens);
        }
    }
}

Tag *parse_header(stList *tokens, char *header_prefix, char *delimiter);

Tag *maf_read_header(LI *li) {
    char *line = LI_get_next_line(li);
    stList *tokens = stString_split(line);
    free(line);
    Tag *tag = parse_header(tokens, "##maf", "=");
    stList_destruct(tokens);

    return tag;
}

void maf_write_block2(Alignment *alignment, LW *lw, bool color_bases) {
    LW_write(lw, "a\n");
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        char *bases = color_bases ? color_base_string(row->bases, alignment->column_number) : row->bases;
        LW_write(lw, "s\t%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%" PRIi64 "\t%s\n", row->sequence_name, row->start, row->length,
                row->strand ? "+" : "-", row->sequence_length, bases);
        row = row->n_row;
        if(color_bases) {
            free(bases);
        }
    }
    LW_write(lw, "\n"); // Add a blank line at the end of the block
}

void maf_write_block(Alignment *alignment, LW *lw) {
    maf_write_block2(alignment, lw, 0);
}

void write_header(Tag *tag, LW *lw, char *header_prefix, char *delimiter, char *end);

void maf_write_header(Tag *tag, LW *lw) {
    write_header(tag, lw, "##maf", "=", "\n\n");
}




