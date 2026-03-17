//
// Created by Benedict Paten on 8/25/22.

#include "taf.h"
#include "sonLib.h"
#include "line_iterator.h"

Alignment *maf_read_block(LI *li) {
    char *tokens_buf[MAX_TOKENS];
    while(1) {
        char *line = LI_get_next_line(li);
        if(line == NULL) {
            return NULL;
        }
        int n_tokens = tokenize_line(line, tokens_buf, MAX_TOKENS);
        if(n_tokens == 0) {
            free(line);
            continue;
        }
        if(strcmp(tokens_buf[0], "a") == 0) { // If is an "a" line
            free(line);
            Alignment *alignment = st_calloc(1, sizeof(Alignment));
            Alignment_Row **p_row = &(alignment->row);
            while(1) {
                line = LI_get_next_line(li);
                if(line == NULL) {
                    return alignment;
                }
                n_tokens = tokenize_line(line, tokens_buf, MAX_TOKENS);
                if(n_tokens == 0) {
                    free(line);
                    return alignment;
                }
                if(strcmp(tokens_buf[0], "s") != 0) {
                    assert(strcmp(tokens_buf[0], "i") == 0 || strcmp(tokens_buf[0], "e") == 0); // Must be an "i" or "e" line, which we ignore
                    free(line);
                    continue;
                }
                assert(strcmp(tokens_buf[0], "s") == 0); // Must be an "s" line
                Alignment_Row *row = st_calloc(1, sizeof(Alignment_Row));
                alignment->row_number++;
                *p_row = row;
                p_row = &(row->n_row);
                row->sequence_name = stString_copy(tokens_buf[1]);
                row->start = atol(tokens_buf[2]);
                row->length = atol(tokens_buf[3]);
                assert(strcmp(tokens_buf[4], "+") == 0 || strcmp(tokens_buf[4], "-") == 0);
                row->strand = strcmp(tokens_buf[4], "+") == 0;
                row->sequence_length = atol(tokens_buf[5]);
                row->bases = stString_copy(tokens_buf[6]);
                free(line);
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
            assert(strcmp(tokens_buf[0], "s") != 0); // Can not be an s line without a prior a line - we will ignore this line
            free(line);
        }
    }
}

Tag *parse_header(char **tokens, int n_tokens, char *header_prefix, char *delimiter);

Tag *maf_read_header(LI *li) {
    char *line = LI_get_next_line(li);
    char *tokens_buf[MAX_TOKENS];
    int n_tokens = tokenize_line(line, tokens_buf, MAX_TOKENS);
    Tag *tag = parse_header(tokens_buf, n_tokens, "##maf", "=");
    free(line);
    return tag;
}

void maf_write_block2(Alignment *alignment, LW *lw, bool color_bases) {
    LW_write(lw, "a\n");
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        char *bases = color_bases ? color_base_string(row->bases, alignment->column_number) : row->bases;
        int bases_len = color_bases ? (int)strlen(bases) : (int)alignment->column_number;
        LW_write(lw, "s\t%s\t%" PRIi64 "\t%" PRIi64 "\t%c\t%" PRIi64 "\t",
                 row->sequence_name, row->start, row->length,
                 row->strand ? '+' : '-', row->sequence_length);
        LW_write_bytes(lw, bases, bases_len);
        LW_write(lw, "\n");
        if(color_bases) {
            free(bases);
        }
        row = row->n_row;
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




