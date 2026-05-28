//
// Created by Benedict Paten on 8/25/22.

#include "taf.h"
#include "sonLib.h"
#include "line_iterator.h"

// Hand-rolled "s name start length strand srcSize bases" parser.  Mutates
// `line` in place by NUL-terminating field boundaries; returned pointers
// reference its bytes (the row owns its own strdup'd name + bases though).
// Returns true on a well-formed s line.
//
// Why: stString_split allocated ~7 strings + an stList per s line and
// freed five immediately after; the per-line malloc/free traffic
// dominated MAF reads.  This parser walks the line once with no
// allocations beyond the two strdups we actually keep.
static bool maf_parse_s_line(char *line, Alignment_Row *row) {
    char *p = line;
    if (*p++ != 's') return false;
    if (*p != '\t' && *p != ' ') return false;
    while (*p == '\t' || *p == ' ') p++;

    // Field 1: sequence name (run to next whitespace).
    char *name = p;
    while (*p && *p != '\t' && *p != ' ') p++;
    if (!*p) return false;
    *p++ = '\0';
    while (*p == '\t' || *p == ' ') p++;
    row->sequence_name = stString_copy(name);

    // Field 2: start.  strtoll consumes leading whitespace and stops at
    // the first non-digit, advancing p to that boundary.
    row->start = (int64_t)strtoll(p, &p, 10);
    while (*p == '\t' || *p == ' ') p++;

    // Field 3: alignment length.
    row->length = (int64_t)strtoll(p, &p, 10);
    while (*p == '\t' || *p == ' ') p++;

    // Field 4: strand: '+' or '-' (single char, then whitespace).
    if (*p != '+' && *p != '-') return false;
    row->strand = (*p == '+');
    p++;
    while (*p == '\t' || *p == ' ') p++;

    // Field 5: source sequence length.
    row->sequence_length = (int64_t)strtoll(p, &p, 10);
    while (*p == '\t' || *p == ' ') p++;

    // Field 6: bases (to end-of-line / NUL).  No internal whitespace.
    char *bases = p;
    while (*p && *p != '\n' && *p != '\r') p++;
    *p = '\0';
    row->bases = stString_copy(bases);

    return true;
}

Alignment *maf_read_block(LI *li) {
    while(1) {
        char *line = LI_get_next_line(li);
        if(line == NULL) {
            return NULL;
        }
        // Cheap byte-level dispatch on the leading char (and that the
        // following byte is whitespace or end-of-line, so we don't match
        // a sequence name that happens to start with 'a').  The old code
        // tokenised every line via stString_split before testing this.
        char c = line[0];
        char c1 = line[1];
        bool is_a_line = (c == 'a' && (c1 == '\t' || c1 == ' ' ||
                                       c1 == '\n' || c1 == '\r' || c1 == '\0'));
        if (!is_a_line) {
            // Skip blank lines, comments, stray e/i/q/s lines outside a
            // block.  Out-of-block "s" lines were assertion errors in the
            // previous code path; preserve the same hard failure.
            assert(!(c == 's' && (c1 == '\t' || c1 == ' ')));
            free(line);
            continue;
        }
        free(line);

        Alignment *alignment = st_calloc(1, sizeof(Alignment));
        Alignment_Row **p_row = &(alignment->row);
        while(1) {
            line = LI_get_next_line(li);
            if(line == NULL) {
                return alignment;
            }
            c = line[0];
            c1 = line[1];
            if (c == '\0' || c == '\n' || c == '\r') { // blank line = block end
                free(line);
                return alignment;
            }
            if (c == 's' && (c1 == '\t' || c1 == ' ')) {
                Alignment_Row *row = st_calloc(1, sizeof(Alignment_Row));
                bool ok = maf_parse_s_line(line, row);
                assert(ok);
                alignment->row_number++;
                *p_row = row;
                p_row = &(row->n_row);
                if(alignment->row_number == 1) {
                    alignment->column_number = strlen(row->bases);
                    alignment->column_tags = st_calloc(alignment->column_number, sizeof(Tag *));
                }
                else {
                    assert(alignment->column_number == strlen(row->bases));
                }
                free(line);
                continue;
            }
            // i/e/q lines (or anything else) -- ignore, matching prior behaviour.
            assert(c == 'i' || c == 'e' || c == 'q');
            free(line);
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
    LW_putn(lw, "a\n", 2);
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        char *bases = color_bases ? color_base_string(row->bases, alignment->column_number) : row->bases;
        // "s\t%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%" PRIi64 "\t%s\n"
        LW_putn(lw, "s\t", 2);
        LW_puts(lw, row->sequence_name);
        LW_putc(lw, '\t');
        LW_puti64(lw, row->start);
        LW_putc(lw, '\t');
        LW_puti64(lw, row->length);
        LW_putc(lw, '\t');
        LW_putc(lw, row->strand ? '+' : '-');
        LW_putc(lw, '\t');
        LW_puti64(lw, row->sequence_length);
        LW_putc(lw, '\t');
        LW_puts(lw, bases);
        LW_putc(lw, '\n');
        row = row->n_row;
        if(color_bases) {
            free(bases);
        }
    }
    LW_putc(lw, '\n'); // Add a blank line at the end of the block
}

void maf_write_block(Alignment *alignment, LW *lw) {
    maf_write_block2(alignment, lw, 0);
}

void write_header(Tag *tag, LW *lw, char *header_prefix, char *delimiter, char *end);

void maf_write_header(Tag *tag, LW *lw) {
    write_header(tag, lw, "##maf", "=", "\n\n");
}




