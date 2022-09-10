//
// Created by Benedict Paten on 8/25/22.

#include "taf.h"

Alignment *maf_read_block(FILE *fh) {
    while(1) {
        char *line = stFile_getLineFromFile(fh);
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
                line = stFile_getLineFromFile(fh);
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
            }
            return alignment;
        }
        else {
            assert(strcmp(stList_get(tokens, 0), "s") != 0); // Can not be an s line without a prior a line - we will ignore this line
            stList_destruct(tokens);
        }
    }
}

stList *maf_read_header(FILE *fh) {
    char *line = stFile_getLineFromFile(fh);
    stList *tokens = stString_split(line);
    free(line);
    stList *tags = parse_header(tokens, "##maf", "=");
    stList_destruct(tokens);

    return tags;
}

void maf_write_block(Alignment *alignment, FILE *fh) {
    fprintf(fh, "a\n");
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        fprintf(fh, "s\t%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%" PRIi64 "\t%s\n", row->sequence_name, row->start, row->length,
                row->strand ? "+" : "-", row->sequence_length, row->bases);
        row = row->n_row;
    }
    fprintf(fh, "\n"); // Add a blank line at the end of the block
}

void maf_write_header(stList *tags, FILE *fh) {
    assert(stList_length(tags) % 2 == 0); // list must be a sequence of alternative key:value pairs
    fprintf(fh, "##maf");
    for(int64_t i=0; i<stList_length(tags); i+=2) {
        fprintf(fh, " %s=%s", (char *)stList_get(tags, i), (char *)stList_get(tags, i+1));
    }
    fprintf(fh, "\n\n"); // Add an extra new line to space between the header and the a lines
}



