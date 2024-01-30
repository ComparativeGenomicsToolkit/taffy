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

void maf_write_block(Alignment *alignment, LW *lw) {
    LW_write(lw, "a\n");

    // check for base quality ('q' tag). Just looking at position[0] to not be too inefficient
    // so the assumption is either you have a quality for every base in every column, or nothing
    // will fail on an error if this doesn't hold (but we can relax if needed later)
    bool has_qualities = false;
    Tag **col_qualities = NULL;
    char *qual_buffer = NULL;
    if (alignment->column_number > 0) { 
        for (Tag *col_tag = alignment->column_tags[0]; col_tag && !has_qualities; col_tag = col_tag->n_tag) {
            if (strcmp(col_tag->key, TAF_BASE_QUALITY_TAG_KEY) == 0) {
                has_qualities = true;
            }
        }
        if (has_qualities) {
            col_qualities = (Tag**)st_calloc(alignment->column_number, sizeof(Tag*));
            qual_buffer = (char*)st_calloc(alignment->column_number + 1, sizeof(char));
            // if we have qualites, fill in col_qualities so we don't have to fish for the tags again
            for (int64_t col = 0; col < alignment->column_number; ++col) {
                for (Tag *col_tag = alignment->column_tags[col]; col_tag && !col_qualities[col]; col_tag = col_tag->n_tag) {
                    if (strcmp(col_tag->key, TAF_BASE_QUALITY_TAG_KEY) == 0) {
                        col_qualities[col] = col_tag;
                    }
                }
                if (col_qualities[col] == NULL) {
                    // see comment above
                    fprintf(stderr, "Error: missing base quality at column in block with base qualities\n");
                    exit(1);
                }
            }
        }
    }
    
    Alignment_Row *row = alignment->row;
    int64_t row_idx = 0;
    while(row != NULL) {
        LW_write(lw, "s\t%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%" PRIi64 "\t%s\n", row->sequence_name, row->start, row->length,
                row->strand ? "+" : "-", row->sequence_length, row->bases);

        if (has_qualities && row->length > 0) {
            int64_t i = 0;
            for (int64_t col = 0; col < alignment->column_number; ++col) {
                if (row->bases[col] != '-') {
                    // this is an ascii-shifted phred score
                    unsigned char qual = col_qualities[col]->value[row_idx] - (unsigned char)33;
                    // do the transformation shown here
                    // https://genome.ucsc.edu/FAQ/FAQformat.html#format5
                    // MAF quality value = min( floor(actual quality value/5), 9 )
                    qual_buffer[i++] = qual >= 99 ? 'F' : (qual >= 45 ? '9' : '0' + qual/5);
                } else {
                    qual_buffer[i++] = '-';
                }
            }
            assert(i == alignment->column_number);
            qual_buffer[i] = '\0';
            LW_write(lw, "q\t%s\t%s\n", row->sequence_name, qual_buffer);
        }
        
        row = row->n_row;
        ++row_idx;
    }
    LW_write(lw, "\n"); // Add a blank line at the end of the block
    free(qual_buffer);
    free(col_qualities);
}

void write_header(Tag *tag, LW *lw, char *header_prefix, char *delimiter, char *end);

void maf_write_header(Tag *tag, LW *lw) {
    write_header(tag, lw, "##maf", "=", "\n\n");
}




