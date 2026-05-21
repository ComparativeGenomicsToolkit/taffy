//
// Created by Benedict Paten on 8/25/22.

#include "taf.h"
#include "sonLib.h"
#include "line_iterator.h"
#ifdef USE_HTSLIB
#include "htslib/bgzf.h"
#endif

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

// Parse a non-negative decimal integer from `p` into `*out`.  Returns the
// pointer to the first non-digit character (typically whitespace or NUL).
// Faster than atoll because it skips strtoll's locale, sign, and base
// detection; MAF s-line integer fields are unsigned in practice and the
// caller has already done maf_skip_ws.  On the 50 Mb fish profile, atoll +
// strtoll were ~13% of CPU before this replacement.
static inline char *maf_parse_u64(char *p, int64_t *out) {
    int64_t v = 0;
    while ((unsigned)(*p - '0') < 10u) {
        v = v * 10 + (*p - '0');
        ++p;
    }
    *out = v;
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
            int64_t start;
            p = maf_parse_u64(p, &start);
            if (*p == '\0') { free(line); continue; }
            ++p;

            p = maf_skip_ws(p);
            int64_t length;
            p = maf_parse_u64(p, &length);
            if (*p == '\0') { free(line); continue; }
            ++p;

            p = maf_skip_ws(p);
            assert(*p == '+' || *p == '-');
            bool strand = (*p == '+');
            ++p;
            if (*p != ' ' && *p != '\t') { free(line); continue; }  // malformed
            p = maf_skip_ws(p);

            int64_t seq_length;
            p = maf_parse_u64(p, &seq_length);
            if (*p == '\0') { free(line); continue; }
            ++p;

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

// Format an int64 into `buf` (no leading-zero handling, no NUL terminator).
// Returns the number of bytes written.  Faster than snprintf because no
// format-string parsing is involved.
static inline int maf_i64_str(int64_t v, char *buf) {
    if (v == 0) { *buf = '0'; return 1; }
    int neg = v < 0;
    uint64_t u = neg ? -(uint64_t)v : (uint64_t)v;
    char tmp[24]; int n = 0;
    do { tmp[n++] = '0' + (int)(u % 10); u /= 10; } while (u > 0);
    int o = 0;
    if (neg) buf[o++] = '-';
    while (n > 0) buf[o++] = tmp[--n];
    return o;
}

// Write `n` bytes to the LW, routing to bgzf or fh as the LW was opened with.
// Return values are discarded with the same fail-loud convention as LW_write
// (which also doesn't check), since downstream errors surface as a truncated
// output file -- cheap to detect with md5sum, easier than threading errors.
static inline void maf_lw_putn(LW *lw, const char *buf, int n) {
#ifdef USE_HTSLIB
    if (lw->bgzf) {
        ssize_t r = bgzf_write(lw->bgzf, buf, n);
        (void)r;   // suppress -Wunused-result
        return;
    }
#endif
    size_t r = fwrite(buf, 1, n, lw->fh);
    (void)r;
}

void maf_write_block2(Alignment *alignment, LW *lw, bool color_bases) {
    // Hot path: avoid vfprintf's format-string parsing on every row.  vfprintf
    // walks "s\t%s\t%lld\t%lld\t%s\t%lld\t%s\n" once per row; at fish 50 Mb
    // scale that was ~25% of CPU.  Build each row's preamble in a stack buffer
    // with hand-rolled itoa + memcpy and emit it as a single write, then write
    // the bases buffer (already a contiguous heap string) directly, then "\n".
    //
    // color_bases is a viewer-only path that mutates the base characters; keep
    // it on the original LW_write/vfprintf code path -- it isn't perf-critical.
    maf_lw_putn(lw, "a\n", 2);
    int64_t cn = alignment->column_number;
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        if (color_bases) {
            char *cb = color_base_string(row->bases, cn);
            LW_write(lw, "s\t%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%" PRIi64 "\t%s\n",
                     row->sequence_name, row->start, row->length,
                     row->strand ? "+" : "-", row->sequence_length, cb);
            free(cb);
            row = row->n_row;
            continue;
        }
        // Compose the preamble: "s\t<name>\t<start>\t<length>\t<strand>\t<seqLen>\t"
        // 1024 covers sequence_name plus six small ints; we fall back to the
        // slow path defensively if a name ever exceeds that (it shouldn't).
        char pre[1024];
        const char *name = row->sequence_name;
        size_t name_n = strlen(name);
        if (name_n + 64 >= sizeof(pre)) {
            LW_write(lw, "s\t%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%" PRIi64 "\t%s\n",
                     name, row->start, row->length,
                     row->strand ? "+" : "-", row->sequence_length, row->bases);
            row = row->n_row;
            continue;
        }
        int k = 0;
        pre[k++] = 's'; pre[k++] = '\t';
        memcpy(pre + k, name, name_n); k += (int)name_n; pre[k++] = '\t';
        k += maf_i64_str(row->start, pre + k); pre[k++] = '\t';
        k += maf_i64_str(row->length, pre + k); pre[k++] = '\t';
        pre[k++] = row->strand ? '+' : '-'; pre[k++] = '\t';
        k += maf_i64_str(row->sequence_length, pre + k); pre[k++] = '\t';

        maf_lw_putn(lw, pre, k);
        maf_lw_putn(lw, row->bases, (int)cn);
        maf_lw_putn(lw, "\n", 1);
        row = row->n_row;
    }
    maf_lw_putn(lw, "\n", 1);   // Add a blank line at the end of the block
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




