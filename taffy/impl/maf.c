//
// Created by Benedict Paten on 8/25/22.

#include "taf.h"
#include "sonLib.h"
#include "line_iterator.h"

// Skip ' ' and '\t' (the only whitespace MAF uses inside a line).  LI strips
// the trailing newline already so we don't need to worry about \r/\n here.
static inline char *maf_skip_ws(char *p) {
    while (*p == ' ' || *p == '\t') ++p;
    return p;
}

// Find the end of the current token (first whitespace or NUL) and return a
// pointer to it.  Caller is responsible for advancing past it.
static inline char *maf_tok_end(char *p) {
    while (*p && *p != ' ' && *p != '\t') ++p;
    return p;
}

Alignment *maf_read_block(LI *li) {
    // Outer loop: skip blank/comment lines until we find an 'a' block start.
    // Inner loop (entered on 'a'): consume successive 's' rows in place, then
    // return the alignment when we hit a blank line or EOF.
    //
    // The fast path tokenises each row in place on the LI-owned buffer
    // (no stString_split / stList / per-token malloc/free pairs).  Only
    // sequence_name and bases are heap-copied -- they're the only fields
    // we keep beyond the line's lifetime.  See the older stList-based
    // version in git history if you need to compare semantics.
    while (1) {
        char *line = LI_get_next_line(li);
        if (line == NULL) return NULL;
        char *p = maf_skip_ws(line);
        if (*p == '\0') { free(line); continue; }            // blank line
        if (*p != 'a') { free(line); continue; }             // header / comment / non-block
        free(line);

        Alignment *alignment = st_calloc(1, sizeof(Alignment));
        Alignment_Row **p_row = &(alignment->row);

        while (1) {
            line = LI_get_next_line(li);
            if (line == NULL) return alignment;              // EOF mid-block
            p = maf_skip_ws(line);
            if (*p == '\0') { free(line); return alignment; }  // blank line ends block
            if (*p == 'i' || *p == 'e') { free(line); continue; }
            if (*p != 's') {
                // Unrecognised line (could be a '##' comment in the middle of
                // a block).  Original parser ignored these.  Do the same.
                free(line); continue;
            }

            // 's' line: in-place tokenise.  Format:
            //   s  <name>  <start>  <length>  <strand>  <seqLen>  <bases>
            ++p;                                              // skip 's'
            p = maf_skip_ws(p);

            char *name = p;
            p = maf_tok_end(p);
            if (*p == '\0') { free(line); continue; }
            *p++ = '\0';                                      // terminate name

            p = maf_skip_ws(p);
            int64_t start = atoll(p);
            p = maf_tok_end(p);
            if (*p == '\0') { free(line); continue; }
            *p++ = '\0';

            p = maf_skip_ws(p);
            int64_t length = atoll(p);
            p = maf_tok_end(p);
            if (*p == '\0') { free(line); continue; }
            *p++ = '\0';

            p = maf_skip_ws(p);
            assert(*p == '+' || *p == '-');
            bool strand = (*p == '+');
            ++p;
            if (*p != ' ' && *p != '\t') { free(line); continue; }  // malformed
            p = maf_skip_ws(p);

            int64_t seq_length = atoll(p);
            p = maf_tok_end(p);
            if (*p == '\0') { free(line); continue; }
            *p++ = '\0';

            p = maf_skip_ws(p);
            char *bases_start = p;
            // Bases run to end-of-line.  LI strips '\n' but be defensive about
            // trailing whitespace (some MAFs have stray spaces post-bases).
            char *bases_end = p + strlen(p);
            while (bases_end > bases_start && (bases_end[-1] == ' ' || bases_end[-1] == '\t'))
                --bases_end;
            *bases_end = '\0';
            int64_t cn = bases_end - bases_start;

            Alignment_Row *row = st_calloc(1, sizeof(Alignment_Row));
            alignment->row_number++;
            *p_row = row;
            p_row = &(row->n_row);
            row->sequence_name   = stString_copy(name);
            row->start           = start;
            row->length          = length;
            row->strand          = strand;
            row->sequence_length = seq_length;
            row->bases           = stString_copy(bases_start);
            free(line);

            if (alignment->row_number == 1) {
                alignment->column_number = cn;
                alignment->column_tags = st_calloc(alignment->column_number, sizeof(Tag *));
            } else {
                assert(alignment->column_number == cn);
            }
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
    tag = header_capture_hal_comment(li, tag); // preserve the hal2maf tree

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
    write_header(tag, lw, "##maf", "=", "\n"); // tag list (hal tag skipped within)
    Tag *hal = tag_find(tag, (char *) TAF_HAL_TREE_KEY);
    if(hal != NULL) {
        LW_write(lw, "%s %s\n", TAF_HAL_TREE_KEY, hal->value);
    }
    LW_write(lw, "\n"); // blank line separating the header from the blocks
}




