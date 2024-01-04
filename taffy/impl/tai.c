#include "taf.h"
#include "tai.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include <ctype.h>
#include <time.h>

char *tai_path(const char *taf_path) {
    char *ret = (char*)st_calloc(strlen(taf_path) + 5, sizeof(char));
    sprintf(ret, "%s.tai", taf_path);
    return ret;    
}

char *tai_parse_region(const char *region, int64_t *start, int64_t *length) {
    int64_t n = strlen(region);
    char *colon = strrchr(region, ':');
    int64_t contig_length = colon ? colon - region : n;
    *start = 0;
    *length = LONG_MAX;
    if (colon) {
        char *dash = strchr(colon + 1, '-');
        int64_t start_length = dash ? dash - colon - 1: strlen(colon + 1);
        if (start_length > 0) {
            char *start_string = stString_getSubString(colon + 1, 0, start_length);
            for (int64_t i = 0; i < strlen(start_string); ++i) {
                if (!isdigit(start_string[i])) {
                    contig_length = 0;
                    break;
                }
            }
            *start = contig_length > 0 ? atol(start_string) : -1;
            free(start_string);
            *length = 1;
        } else {
            return NULL;
        }
        int64_t end_length = n - (contig_length + start_length + 2);
        if (end_length > 0) {
            for (int64_t i = contig_length + 2 + start_length; i < n; ++i) {
                if (!isdigit(region[i])) {
                    contig_length = 0;
                    break;
                }
            }
            if (contig_length > 0) {
                int64_t end = atol(colon + 2 + start_length);
                if (end < *start) {
                    contig_length = 0;
                } else {
                    *length = end - *start;
                }
            }
        }
    }
    return contig_length > 0 ? stString_getSubString(region, 0, contig_length) : NULL;
}

// gets the "reference" (first row) coordinate information
// but only returns everything if there's a coordinate for every row in the column
// this can happen when everything is an "i" like on the first line
// or when everything is an "s" like on a repeat-coordinates-every-n-columns line
static char *parse_coordinates_line(stList *tokens, int64_t *start, bool *strand,
                                    bool run_length_encode_bases) {
    int64_t j = -1;
    if (!has_coordinates(tokens, &j)) {
        return NULL;
    }

    int64_t num_bases = 0;
    // todo: really should go through api for this one
    if (run_length_encode_bases) {
        assert(j > 1);
        for (int64_t i = 0; i < j; ++i) {
            char *token = (char*)stList_get(tokens, i);
            if (isdigit(token[0])) {
                num_bases += atol(token);
            }
        }
    } else {
        assert(j == 1);
        num_bases = strlen(stList_get(tokens, 0));
    }
    int64_t num_coordinates = 0;
    char *seq = NULL;
    int64_t sequence_length;
    bool dummy;
    int64_t n = stList_length(tokens);

    ++j;
    while (j < n) {
        // copied from taf.c
        char *op_type = stList_get(tokens, j++); // This is the operation
        assert(strlen(op_type) == 1); // Must be a single character in length
        int64_t row_index = atol(stList_get(tokens, j++)); // Get the index of the affected row
        if(op_type[0] == 'i' || op_type[0] == 's') { // We have coordinates!
            num_coordinates++;
            if (row_index == 0) {
                seq = parse_coordinates(&j, tokens, start, strand, &sequence_length);
            } else {
                // we parse but don't use
                // todo: smoother api
                char *s = parse_coordinates(&j, tokens, &sequence_length, &dummy, &sequence_length);
                free(s);
            }
        } else if (op_type[0] == 'd') {
        } else if (op_type[0] == 'g') {
            j++;
        } else {
            assert(op_type[0] == 'G');
            j++;
        }
    }

    if (num_coordinates == num_bases) {
        assert(seq != NULL);
        return seq;
    }

    *start = -1;
    return NULL;
}

// this is rather hacky: because our index points to
// taf lines where all coordinates are listed as "s", but we
// want to start new block parsers on these lines, we have to
// convert the s's to i's to pretend we're starting new files
static void change_s_coordinates_to_i(char *line) {
    stList* tokens = stString_split(line);

    int64_t j = -1;
    bool hc = has_coordinates(tokens, &j);
    int64_t dummy_int;
    bool dummy_bool;
    bool found_s = false;

    if (hc) {
        int64_t n = stList_length(tokens);
        bool *mask = st_calloc(n, sizeof(bool));        
        ++j;
        while (j < n) {
            char *op_type = stList_get(tokens, j++); // This is the operation
            j++; // the row index;
            assert(strlen(op_type) == 1); // Must be a single character in length
            if(op_type[0] == 'i' || op_type[0] == 's') { // We have coordinates!
                found_s = found_s || op_type[0] == 's';
                op_type[0] = 'i';
                // only parse to increment j (would rather use api than just add 4)
                parse_coordinates(&j, tokens, &dummy_int, &dummy_bool, &dummy_int);
            } else if (op_type[0] == 'd') {
                mask[j-2] = true;
                mask[j-1] = true;
            } else if (op_type[0] == 'g') {
                mask[j-2] = true;
                mask[j-1] = true;
                mask[j] = true;
                j++;
            } else {
                assert(op_type[0] == 'G');
                mask[j-2] = true;
                mask[j-1] = true;
                mask[j] = true;
                j++;
            }
        }
        
        if (found_s) {
            // overwrite our line
            int64_t line_len = strlen(line);
            int64_t k = 0;
            int64_t n = stList_length(tokens);
            assert(mask[0] == false);
            for (j = 0; j < n; ++j) {
                char *token = (char*)stList_get(tokens, j);
                if (mask[j] == false) {
                    int64_t token_len = strlen(token);
                    char *spacer = j == 0 ? "" : " ";
                    assert(k + token_len + strlen(spacer) <= line_len);
                    sprintf(line + k, "%s%s", spacer, token);
                    k += strlen(token) + strlen(spacer);
                }
            }
        }
    } else {
        fprintf(stderr, "Error loading coordinates from indexed taf line: %s\n", line);
        exit(1);
    }

    stList_destruct(tokens);
}

static int tai_create_taf(LI *li, FILE *idx_fh, int64_t index_block_size, bool run_length_encode_bases) {
    char *prev_ref = NULL;
    int64_t prev_pos = 0;
    int64_t prev_file_pos = 0;

    // scan the taf line by line
    for (char *line = LI_get_next_line(li); line != NULL; line = LI_get_next_line(li)) {
        stList* tokens = stString_split(line);
        int64_t pos;
        assert(sizeof(int64_t) == sizeof(off_t));
        bool strand;
        char *ref = parse_coordinates_line(tokens, &pos, &strand, run_length_encode_bases);
        if (ref != NULL) {
            // shouldn't need to handle negative strand on reference, right?
            assert(strand == true);

            // we need to update our index if we're on a new reference contig
            // or we're on the same contig but >= index_block_size bases away
            bool same_ref = prev_ref && strcmp(ref, prev_ref) == 0;
            if (!same_ref || pos - prev_pos >= index_block_size) {
                int64_t file_pos = LI_tell(li);
                if (same_ref) {
                    // save a little space by writing relative coordinates
                    fprintf(idx_fh, "*\t%" PRIi64 "\t%" PRIi64 "\n", pos-prev_pos, file_pos-prev_file_pos);
                } else {
                    fprintf(idx_fh, "%s\t%" PRIi64 "\t%" PRIi64 "\n", ref, pos, file_pos);
                }
                free(prev_ref);
                prev_ref = ref;
                prev_pos = pos;
                prev_file_pos = file_pos;
            }
        }
        stList_destruct(tokens);
        free(line);
    }
    free(prev_ref);
    return 0;
}

static int tai_create_maf(LI *li, FILE *idx_fh, int64_t index_block_size) {
    char *prev_ref = NULL;
    int64_t prev_pos = 0;
    int64_t prev_file_pos = 0;

    // scan the maf block by block line by line
    Alignment *alignment, *p_alignment = NULL;
    int64_t file_pos = LI_tell(li);
    while((alignment = maf_read_block(li)) != NULL) {
        if(p_alignment != NULL) {
            alignment_link_adjacent(p_alignment, alignment, 1);
        }
        if (alignment->row->strand != 1) {
            fprintf(stderr, "Can't index maf because reference (row 0) sequence found on negative strand\n");
            exit(1);
        }
        // todo: error message when out of order
        bool same_ref = prev_ref && strcmp(alignment->row->sequence_name, prev_ref) == 0;
        int64_t pos = alignment->row->start;
        char *ref = alignment->row->sequence_name;
        if (!same_ref || pos - prev_pos >= index_block_size) {
            if (same_ref) {
                // save a little space by writing relative coordinates
                fprintf(idx_fh, "*\t%" PRIi64 "\t%" PRIi64 "\n", pos-prev_pos, file_pos-prev_file_pos);
            } else {
                fprintf(idx_fh, "%s\t%" PRIi64 "\t%" PRIi64 "\n", ref, pos, file_pos);
            }
            free(prev_ref);
            prev_ref = stString_copy(ref);
            prev_pos = pos;
            prev_file_pos = file_pos;
        }
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment, 1);
        }
        p_alignment = alignment;
        file_pos = LI_tell(li);                    
    }
    if(p_alignment != NULL) {
        alignment_destruct(p_alignment, 1);
    }
    free(prev_ref);
    return 0;
}
    
int tai_create(LI *li, FILE* idx_fh, int64_t index_block_size) {

    int input_format = check_input_format(LI_peek_at_next_line(li));
    assert(input_format == 0 || input_format == 1);
    
    // read the TAF header
    bool run_length_encode_bases = 0;
    Tag *tag = NULL;
    if (input_format == 0) {
        tag = taf_read_header(li);
        Tag *t = tag_find(tag, "run_length_encode_bases");
        if(t != NULL && strcmp(t->value, "1") == 0) {
            run_length_encode_bases = 1;
        }
    } else {
        tag = maf_read_header(li);
    }
    tag_destruct(tag);

    if (input_format == 0) {
        return tai_create_taf(li, idx_fh, index_block_size, run_length_encode_bases);
    } else {
        return tai_create_maf(li, idx_fh, index_block_size);
    }

    return -1;
}

// basically all the information in the index file
typedef struct _TaiRec {
    char *name;
    int64_t seq_pos;
    int64_t file_pos;
} TaiRec;

static int tai_record_cmp(const void *v1, const void *v2) {
    TaiRec *tr1 = (TaiRec*)v1;
    TaiRec *tr2 = (TaiRec*)v2;
    int ret = strcmp(tr1->name, tr2->name);
    if (ret == 0) {
        if (tr1->seq_pos < tr2->seq_pos) {
            ret = -1;
        } else if (tr1->seq_pos > tr2->seq_pos) {
            ret = 1;
        }
    }
    return ret;
}

static Tai *tai_construct() {
    Tai *tai = st_calloc(1, sizeof(Tai));
    tai->idx = stSortedSet_construct3(tai_record_cmp, free);
    tai->names = stList_construct3(0, free);
    return tai;
}

void tai_destruct(Tai* tai) {
    stSortedSet_destruct(tai->idx);
    stList_destruct(tai->names);
    free(tai);
}

Tai *tai_load(FILE* idx_fh, bool maf) {
    time_t start_time = time(NULL);
    Tai *tai = tai_construct();
    tai->maf = maf;
    LI* li = LI_construct(idx_fh);
    char *line;
    TaiRec *prev_rec = NULL;
    while ((line = LI_get_next_line(li)) != NULL) {
        stList* tokens = stString_splitByString(line, "\t");
        if (stList_length(tokens) != 3) {
            fprintf(stderr, "Skipping tai line that does not have 3 columns: %s\n", line);
            continue;
        }
        TaiRec* rec = (TaiRec*)st_calloc(1, sizeof(TaiRec));
        rec->seq_pos = atol(stList_get(tokens, 1));
        rec->file_pos = atol(stList_get(tokens, 2));
        char *name = (char*)stList_get(tokens, 0);
        if (strcmp(name, "*") == 0) {
            if (prev_rec == NULL) {
                fprintf(stderr, "Unable to deduce name from tai line: %s\n", line);
                exit(1);
            }
            rec->name = prev_rec->name;
            rec->seq_pos += prev_rec->seq_pos;
            rec->file_pos += prev_rec->file_pos;
        } else {
            rec->name = stString_copy(stList_get(tokens, 0));
            stList_append(tai->names, rec->name);
        }
        stSortedSet_insert(tai->idx, rec);
        stList_destruct(tokens);
        prev_rec = rec;
    }
    LI_destruct(li);
    st_logInfo("Loaded .tai index in %" PRIi64 " seconds\n", time(NULL) - start_time);
    return tai;
}

// dummy function to let us toggle between maf/taf reading at runtime (by providing a
// maf reader with same interface as taf reader)
static Alignment *maf_read_block_3(Alignment *p_block, bool run_length_encode_bases, LI *li) {
    Alignment *alignment = maf_read_block(li);
    if(p_block != NULL) {
        alignment_link_adjacent(p_block, alignment, 1);
    }
    return alignment;
}

TaiIt *tai_iterator(Tai* tai, LI *li, bool run_length_encode_bases, const char *contig, int64_t start, int64_t length) {
    time_t start_time = time(NULL);
    TaiIt *tai_it = st_calloc(1, sizeof(TaiIt));
    tai_it->maf = tai->maf;

    // parse the region into the iterator
    tai_it->name = stString_copy(contig);
    tai_it->start = start;
    tai_it->end = length < 0 ? LONG_MAX : start + length;
    tai_it->run_length_encode_bases = run_length_encode_bases;

    // look up the region in the taf index, which takes a dummy record
    TaiRec qr;
    qr.name = tai_it->name;
    qr.seq_pos = tai_it->start;

    TaiRec *tair_1 = stSortedSet_searchLessThanOrEqual(tai->idx, &qr);
    if (tair_1 == NULL) {
        tai_iterator_destruct(tai_it);
        // there's no chance of finding the region as its contig isn't
        // in the index or its start position is too low
        // up to caller to spit out an error
        return NULL;
    }

    // sorted set doesn't let us iterate, so we use a second lookup to get
    // the next record
    TaiRec qr2;
    qr2.name = tai_it->name;
    qr2.seq_pos = tai_it->end;
    TaiRec *tair_2 = stSortedSet_searchGreaterThanOrEqual(tai->idx, &qr2);
    st_logInfo("Queried the in-memory .tai index in %" PRIi64 " seconds\n", time(NULL) - start_time);

    // now we know that the start of our region is somewhere in [tair_1, tair_2)
    // (with the possibility of tair_2 not existing)

    // move to the first record in our file
    start_time = time(NULL);
    LI_seek(li, tair_1->file_pos);
    st_logInfo("Seeked to the queried anchor position with taf file in %" PRIi64 " seconds\n", time(NULL) - start_time);
    LI_get_next_line(li);

    // all maf / taf logic toggling is handled right here
    Alignment* (*maftaf_read_block)(Alignment*, bool, LI*) = maf_read_block_3;
    if (!tai_it->maf) {
        maftaf_read_block = taf_read_block;
        // force taf to start a new alignment at our current file position by making
        // sure all coordinates are expressed as insertions
        change_s_coordinates_to_i(LI_peek_at_next_line(li));
    }
    
    // now we have to scan forward until we overlap actually the region
    // TODO: this will surely need speeding up with binary search for giant files...
    //       but i think that's best done once tests are set up based on the simpler version
    start_time = time(NULL);
    size_t scan_block_count = 0;
    tai_it->alignment = NULL;
    tai_it->p_alignment = NULL;
    Alignment *alignment = NULL;
    Alignment *p_alignment = NULL;
    int64_t file_pos = LI_tell(li);
    while((alignment = maftaf_read_block(p_alignment, tai_it->run_length_encode_bases, li)) != NULL) {
        ++scan_block_count;
        if (tair_2 && file_pos >= tair_2->file_pos) {
            // we've gone past our query region: there's no hope
            alignment_destruct(alignment, true);
            if (p_alignment) {
                alignment_destruct(p_alignment, true);
            }
            alignment = NULL;
            p_alignment = NULL;
            break;
        } else {
            if (p_alignment != NULL) {
                alignment_destruct(p_alignment, true);
            }
            p_alignment = alignment;
            if (strcmp(alignment->row->sequence_name, qr.name) == 0 &&
                alignment->row->start < tai_it->end &&
                (alignment->row->start + alignment->row->length) > tai_it->start) {
                // important: need to cut off p_alignment to get our absolute coordinates
                for (Alignment_Row *row = alignment->row; row != NULL; row = row->n_row) {
                    row->l_row = NULL;
                }
                // we've found an intersection at "alignment"
                tai_it->alignment = alignment;
                p_alignment = NULL;
                break;
            }
        }
        file_pos = LI_tell(li);
    }
    if (p_alignment != NULL) {
        alignment_destruct(p_alignment, true);
        p_alignment = NULL;
    }

    // we scanned past our region, which could happen if your taf doesn't contain the whole
    // sequence -- jsut return nothing
    if (tai_it->alignment == NULL) {
        if (tai_it->p_alignment) {
            alignment_destruct(tai_it->p_alignment, true);
        }
        tai_iterator_destruct(tai_it);
        st_logInfo("Scanned %" PRIi64 " blocks to NOT find region start in %" PRIi64 " seconds\n", scan_block_count, time(NULL) - start_time);
        return NULL;
    }

    st_logInfo("Scanned %" PRIi64 " blocks to find region start in %" PRIi64 " seconds\n", scan_block_count, time(NULL) - start_time);
    
    return tai_it;
}

/** 
 * clip an alignment, assumes start/end overlap aln and would leave a non-empty remainder
 * returns 0: no cut
 *         1: right side cut 
 *         2: left side cut
 *         3: left and right side cut
 */ 
static unsigned int clip_alignment(Alignment *aln, Alignment *p_aln, int64_t start, int64_t end) {

    unsigned int ret = 0;
    // clip the left side
    int64_t left_trim = start - aln->row->start;
    if (left_trim > 0) {
        ret = ret | 2;
        // we assume that the current alignment overlaps our range
        assert(aln->column_number > left_trim);

        // we need to find the cut point by counting off left_trim non-gap bases from the end of the reference row
        int64_t cut_point = 0;
        for(int64_t cut_count = 0; cut_count < left_trim && cut_point < aln->row->length; ++cut_point) {
            if (aln->row->bases[cut_point] != '-') {
                ++cut_count;
            }
        }
        // then we trim out every row, making sure that we adjust the start/length fields only
        // by non-gap bases we removed
        for (Alignment_Row *row = aln->row; row != NULL; row = row->n_row) {
            for (int64_t col = 0; col < cut_point; ++col) {
                if (row->bases[col] != '-') {
                    ++row->start;
                    --row->length;
                }
            }
            if (row->length == 0) {
                row->bases[0] = '\0';
            } else {
                char *bases = row->bases;
                row->bases = stString_getSubString(bases, cut_point, strlen(row->bases) - cut_point);
                free(bases);
            }
            assert(strlen(row->bases) >= row->length);            
        }
        aln->column_number -= left_trim;
    }

    //clip the right side
    int64_t right_trim = (aln->row->start + aln->row->length) - end;
    if (right_trim > 0) {
        ret = ret | 1;
        assert(aln->column_number > right_trim);

        // we need to find the cut point by counting off right_trim non-gap bases from the end of the reference row
        int64_t cut_point = strlen(aln->row->bases) - 1;
        for (int64_t cut_count = 0; cut_count < right_trim && cut_point >= 0; --cut_point) {
            if (aln->row->bases[cut_point] != '-') {
                ++cut_count;
            }
        }
        // then we trim out every row, making sure that we adjust the start/length fields only
        // by non-gap bases we removed
        for (Alignment_Row *row = aln->row; row != NULL; row = row->n_row) {
            int64_t rowlen = strlen(row->bases);
            for (int64_t col = rowlen - 1; col > cut_point; --col) {
                if (row->bases[col] != '-') {
                    --row->length;
                }
            }
            if (row->length == 0) {
                row->bases[0] = '\0';
            } else {
                char *bases = row->bases;
                row->bases = stString_getSubString(bases, 0, cut_point + 1);
                free(bases);
            }
            assert(strlen(row->bases) >= row->length);
        }
        aln->column_number -= right_trim;                
    }

    // now we make sure any deleted rows are unlinked from prev alignment.
    if (p_aln) {
        for (Alignment_Row *row = p_aln->row; row != NULL; row = row->n_row) {
            if (row->r_row && row->r_row->length == 0) {
                assert(strlen(row->r_row->bases) == 0);
                row->r_row = NULL;
            }
        }
    }
    
    // then get rid of empty rows as they api doesn't handle them
    Alignment_Row *prev = NULL;
    Alignment_Row *next = NULL;
    for (Alignment_Row *row = aln->row; row != NULL; row = next) {
        next = row->n_row;
        if (row->length == 0) {
            assert(strlen(row->bases) == 0);
            assert(prev != NULL);
            prev->n_row = next;
            alignment_row_destruct(row);
            --aln->row_number;
        } else {
            prev = row;
        }
    }
    assert(aln->column_number > 0);
    
    return ret;
}

Alignment *tai_next(TaiIt *tai_it, LI *li) {    
    if (tai_it->alignment == NULL) {
        return NULL;
    }

    // start by clamping the alignment block to the region
    assert(strcmp(tai_it->alignment->row->sequence_name, tai_it->name) == 0);
    unsigned int ret = clip_alignment(tai_it->alignment, tai_it->p_alignment, tai_it->start, tai_it->end);

    // save this alignment, it's what we're gonna return
    tai_it->p_alignment = tai_it->alignment;

    Alignment* (*maftaf_read_block)(Alignment*, bool, LI*) = maf_read_block_3;
    if (!tai_it->maf) {
        maftaf_read_block = taf_read_block;
    }
    
    // scan forward
    if (ret & 1) {
        tai_it->alignment = NULL;
    } else {
        tai_it->alignment = maftaf_read_block(tai_it->p_alignment, tai_it->run_length_encode_bases, li);
        if (tai_it->alignment != NULL &&
            (strcmp(tai_it->alignment->row->sequence_name, tai_it->name) != 0 ||
             tai_it->alignment->row->start >= tai_it->end)) {
            alignment_destruct(tai_it->alignment, true);
            tai_it->alignment = NULL;
        }
    }
    
    return tai_it->p_alignment;
}

void tai_iterator_destruct(TaiIt *tai_it) {
    free(tai_it->name);
    free(tai_it);
}

stHash *tai_sequence_lengths(Tai *tai, LI *li) {
    // read the header
    LI_seek(li, 0);
    LI_get_next_line(li);
    int input_format = check_input_format(LI_peek_at_next_line(li));
    assert(input_format == 0 || input_format == 1);
    assert((input_format == 1) == tai->maf);
    bool run_length_encode_bases = 0;
    Tag *tag = NULL;
    if (tai->maf == 0) {
        tag = taf_read_header(li);
        Tag *t = tag_find(tag, "run_length_encode_bases");
        if(t != NULL && strcmp(t->value, "1") == 0) {
            run_length_encode_bases = 1;
        }
    } else {
        tag = maf_read_header(li);
    }
    tag_destruct(tag);

    Alignment* (*maftaf_read_block)(Alignment*, bool, LI*) = maf_read_block_3;
    if (!tai->maf) {
        maftaf_read_block = taf_read_block;
    }
    
    stHash *seq_to_len = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);

    // dummy record for querying the set
    TaiRec qr;
    qr.seq_pos = 0;

    // the index keeps a set of sequence names, which is handy here
    // we iterate that, using the position index to hook into the first occurrence of
    // each name.  but we add every sequence name we can find for each block, not just reference
    for (int64_t i = 0; i < stList_length(tai->names); ++i) {
        char *seq = (char*)stList_get(tai->names, i);
        if (stHash_search(seq_to_len, seq) == NULL) {
            qr.name = seq;
            TaiRec *tair = stSortedSet_searchGreaterThanOrEqual(tai->idx, &qr);
            assert(tair != NULL);
            LI_seek(li, tair->file_pos);
            LI_get_next_line(li);

            if (!tai->maf) {
                // force taf to start a new alignment at our current file position by making
                // sure all coordinates are expressed as insertions
                change_s_coordinates_to_i(LI_peek_at_next_line(li));
            }
            
            Alignment *alignment = maftaf_read_block(NULL, run_length_encode_bases, li);
            assert(alignment != NULL);

            Alignment_Row *row = alignment->row;
            assert(row != NULL);
            assert(strcmp(row->sequence_name, seq) == 0);
            stHash_insert(seq_to_len, stString_copy(row->sequence_name), (void*)row->sequence_length);
        }
    }

    assert(stHash_size(seq_to_len) == stList_length(tai->names));
    return seq_to_len;
}
