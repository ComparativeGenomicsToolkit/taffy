#include "taf.h"
#include "tai.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include <ctype.h>

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
                parse_coordinates(&j, tokens, &sequence_length, &dummy, &sequence_length);
            }
        } else {
            seq = NULL;
            break;
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
        ++j;
        while (j < n) {
            char *op_type = stList_get(tokens, j++); // This is the operation
            assert(strlen(op_type) == 1); // Must be a single character in length
            if(op_type[0] == 'i' || op_type[0] == 's') { // We have coordinates!
                found_s = op_type[0] == 's';
                op_type[0] = 'i';
                j++; // the row index;
                // only parse to increment j (would rather use api than just add 4)
                parse_coordinates(&j, tokens, &dummy_int, &dummy_bool, &dummy_int);
            } else {
                hc = false;
                break;
            }
        }
    }

    if (hc) {
        if (found_s) {
            // overwrite our line
            int64_t line_len = strlen(line);
            int64_t k = 0;
            int64_t n = stList_length(tokens);
            for (j = 0; j < n; ++j) {
                char *token = (char*)stList_get(tokens, j);
                int64_t token_len = strlen(token);
                char *spacer = j == 0 ? "" : " ";
                assert(k + token_len + strlen(spacer) <= line_len);
                sprintf(line + k, "%s%s", spacer, token);
                k += strlen(token) + strlen(spacer);
            }
        }
    } else {
        fprintf(stderr, "Error loading coordinates from indexed taf line: %s\n", line);
        exit(1);
    }

    stList_destruct(tokens);
}

int tai_index(LI *li, FILE* idx_fh, int64_t index_block_size){
    char *prev_ref = NULL;
    int64_t prev_pos = 0;
    
    // read the TAF header
    stList *tags = taf_read_header(li);
    assert(stList_length(tags) % 2 == 0);
    bool run_length_encode_bases = 0;
    for(int64_t i=0; i<stList_length(tags); i+=2) {
        char *key = (char*)stList_get(tags, i);
        char *value = (char*)stList_get(tags, i+1);
        if(strcmp(key, "run_length_encode_bases") == 0 && strcmp(value, "1") == 0) {
            run_length_encode_bases = 1;
        }
    }

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
            if (!prev_ref || strcmp(ref, prev_ref) != 0 ||
                pos - prev_pos >= index_block_size) {
                // we write contig - position - file_offset and that's it
                fprintf(idx_fh, "%s\t%" PRIi64 "\t%" PRIi64 "\n", ref, pos, LI_tell(li));
            }
            prev_ref = ref;
            prev_pos = pos;
        }
        stList_destruct(tokens);
        free(line);
    }
    return 0;
}

// basically all the information in the index file
typedef struct _TaiRec {
    char *name;
    int64_t seq_pos;
    int64_t file_pos;
} TaiRec;

static int taf_record_cmp(const void *v1, const void *v2) {
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

Tai *tai_load(FILE* idx_fh) {
    stSortedSet *taf_index = stSortedSet_construct3(taf_record_cmp, free);
    LI* li = LI_construct(idx_fh);
    char *line;
    while ((line = LI_get_next_line(li)) != NULL) {
        stList* tokens = stString_splitByString(line, "\t");
        if (stList_length(tokens) != 3) {
            fprintf(stderr, "Skipping tai line that does not have 3 columns: %s\n", line);
            continue;
        }
        TaiRec* tr = (TaiRec*)st_calloc(1, sizeof(TaiRec));
        tr->name = stString_copy(stList_get(tokens, 0));
        tr->seq_pos = atol(stList_get(tokens, 1));
        tr->file_pos = atol(stList_get(tokens, 2));
        stSortedSet_insert(taf_index, tr);
        stList_destruct(tokens);
    }
    LI_destruct(li);
    return taf_index;
}

void tai_destruct(Tai* idx) {
    stSortedSet_destruct(idx);
}

struct _TaiIt {
    char *name;
    // these are bed-like 0-based half-open
    int64_t start;
    int64_t end;
    Alignment *alignment;
    Alignment *p_alignment;
    bool run_length_encode_bases;
};

TaiIt *tai_iterator(Tai* tai, LI *li, bool run_length_encode_bases, const char *contig, int64_t start, int64_t length) {

    TaiIt *tai_it = st_calloc(1, sizeof(TaiIt));

    // parse the region into the iterator
    tai_it->name = stString_copy(contig);
    tai_it->start = start;
    tai_it->end = length < 0 ? LONG_MAX : start + length;
    tai_it->run_length_encode_bases = run_length_encode_bases;

    // look up the region in the taf index, which takes a dummy record
    TaiRec qr;
    qr.name = tai_it->name;
    qr.seq_pos = tai_it->start;

    TaiRec *tair_1 = stSortedSet_searchLessThanOrEqual(tai, &qr);
    if (tair_1 == NULL) {
        tai_iterator_destruct(tai_it);
        // there's no chance of finding the region as its contig isn't
        // in the index or its start position is too low
        // up to caller to spit out an error
        return NULL;
    }

    // sorted set doesn't let us iterate, so we use a second lookup to get
    // the next record
    TaiRec *tair_2 = stSortedSet_searchGreaterThan(tai, &qr);

    // now we know that the start of our region is somewhere in [tair_1, tair_2)
    // (with the possibility of tair_2 not existing)

    // move to the first record in our file
    LI_seek(li, tair_1->file_pos);
    LI_get_next_line(li);

    // force taf to start a new alignment at our current file position by making
    // sure all coordinates are expressed as insertions
    change_s_coordinates_to_i(li->line);
    
    // now we have to scan forward until we overlap actually the region
    // TODO: this will surely need speeding up with binary search for giant files...
    //       but i think that's best done once tests are set up based on the simpler version
    tai_it->alignment = NULL;
    tai_it->p_alignment = NULL;
    Alignment *alignment = NULL;
    Alignment *p_alignment = NULL;
    int64_t file_pos = LI_tell(li);
    while((alignment = taf_read_block(p_alignment, tai_it->run_length_encode_bases, li)) != NULL) {
        if (tair_2 && file_pos >= tair_2->file_pos) {
            // we've gone past our query region: there's no hope
            alignment_destruct(alignment);
            if (p_alignment) {
                alignment_destruct(p_alignment);
            }
            alignment = NULL;
            p_alignment = NULL;
            break;
        } else {
            if (p_alignment != NULL) {
                alignment_destruct(p_alignment);
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
        alignment_destruct(p_alignment);
        p_alignment = NULL;
    }

    // we scanned past our region, which could happen if your taf doesn't contain the whole
    // sequence -- jsut return nothing
    if (tai_it->alignment == NULL) {
        if (tai_it->p_alignment) {
            alignment_destruct(tai_it->p_alignment);
        }
        tai_iterator_destruct(tai_it);
        return NULL;
    }
    
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
            Alignment_Row_destruct(row);
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

    // scan forward
    if (ret & 1) {
        tai_it->alignment = NULL;
    } else {
        tai_it->alignment = taf_read_block(tai_it->p_alignment, tai_it->run_length_encode_bases, li);
        if (tai_it->alignment != NULL &&
            (strcmp(tai_it->alignment->row->sequence_name, tai_it->name) != 0 ||
             tai_it->alignment->row->start >= tai_it->end)) {
            alignment_destruct(tai_it->alignment);
            tai_it->alignment = NULL;
        }
    }
    
    return tai_it->p_alignment;
}
void tai_iterator_destruct(TaiIt *tai_it) {
    free(tai_it->name);
    free(tai_it);
}
