#include "taf.h"
#include "taf_index.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h" // just for hts_parse_reg()
#include <ctype.h>

char *tai_path(const char *taf_path) {
    char *ret = (char*)st_calloc(strlen(taf_path) + 5, sizeof(char));
    sprintf(ret, "%s.tai", taf_path);
    return ret;    
}

// todo: it would be nice to be able to use something from the taf.h API for this!!
// but it's currently inside write_coordinates()
static void write_tai_coorindates(Alignment *alignment, FILE* fh) {
    int64_t i = 0;
    for (Alignment_Row *row = alignment->row; row; row = row->n_row, ++i) {
        // the first row is already present in our index as a "reference" column
        if (i > 0) {
            if (i > 1) {
                fprintf(fh, " ");
            }
            fprintf(fh, "i %" PRIi64 " %s %" PRIi64 " %c %" PRIi64 "",
                    i, row->sequence_name, row->start, row->strand ? '+' : '-', row->sequence_length);
        }
    }
}

int tai_index(LI *li, FILE* idx_fh, int64_t index_block_size){
    Alignment *alignment, *p_alignment = NULL;
    char *prev_ref = NULL;
    int64_t prev_pos = 0;
    int64_t prev_file_offset = bgzf_utell(li->bgzf);
    
    // read the TAF header
    stList *tags = taf_read_header(li);    
    assert(stList_length(tags) % 2 == 0);

    while((alignment = taf_read_block(p_alignment, false, li)) != NULL) {

        char *cur_ref = alignment->row->sequence_name;
        int64_t cur_offset = alignment->row->start;

        // shouldn't need to handle negative strand on reference, right?
        assert(alignment->row->strand != 0);

        // we need to update our index if we're on a new reference contig
        // or we're on the same contig but >= index_block_size bases away
        if (!prev_ref || strcmp(cur_ref, prev_ref) != 0 ||
            cur_offset - prev_pos >= index_block_size) {

            fprintf(idx_fh, "%s\t%ld\t%ld\t%ld\t", cur_ref, cur_offset, prev_file_offset, alignment->row->sequence_length);
            write_tai_coorindates(alignment, idx_fh);
            fprintf(idx_fh, "\n");
        }

        prev_ref = alignment->row->sequence_name;
        prev_pos = alignment->row->start;
        prev_file_offset = bgzf_utell(li->bgzf);
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment);
        }
        p_alignment = alignment;

    }
    if(p_alignment != NULL) {
        alignment_destruct(p_alignment);
    }

    stList_destruct(tags);
    
    return 0;
}

// basically all the information in the index file
typedef struct _TaiRec {
    char *name;
    int64_t seq_pos;
    int64_t file_pos;
    int64_t seq_len;
    char *coords;    
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
        if (stList_length(tokens) != 5) {
            fprintf(stderr, "Skipping tai line that does not have 4 columns: %s\n", line);
            continue;
        }
        TaiRec* tr = (TaiRec*)st_calloc(1, sizeof(TaiRec));
        tr->name = stString_copy(stList_get(tokens, 0));
        tr->seq_pos = atol(stList_get(tokens, 1));
        tr->file_pos = atol(stList_get(tokens, 2));
        tr->seq_len = atol(stList_get(tokens, 3));
        tr->coords = stString_copy(stList_get(tokens, 4));
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
};

/*
 * Overwrite (!) LI's line buffer with a version that has full coordinate 
 * specification as obtained from the TaiRec.
 *
 * TODO: this is rather dirty, but simpler (and faster) than the alternative
 *       of scanning to the nearest full coordinate specification.  But...
 *       if/when we move to a binary-search mode that relies on scanning to full
 *       coordinates anyway, it will no longer be necessary
 */
static void inject_full_coorindates(TaiRec *tr, LI *li) {
    assert(li->line);
    size_t end_of_bases = 0;
    while (li->line[end_of_bases] != '\0' && !isspace(li->line[end_of_bases])) {
        ++end_of_bases;
    }
    char *old_line = li->line;
    li->line = st_calloc(1, sizeof(char) * (strlen(old_line) + strlen(tr->name) + strlen(tr->coords) + 256));
    strncpy(li->line, old_line, end_of_bases);
    free(old_line);
    sprintf(li->line + end_of_bases, " ; i 0 %s %ld + %ld %s", tr->name, tr->seq_pos, tr->seq_len, tr->coords);
}

TaiIt *tai_iterator(Tai* tai, LI *li, const char *region) {

    TaiIt *tai_it = st_calloc(1, sizeof(TaiIt));

    // parse the region into the iterator
    const char* ce = hts_parse_reg64(region, &tai_it->start, &tai_it->end);
    tai_it->name = stString_getSubString(region, 0, ce - region);

    // look up the region in the taf index, which takes a dummy record
    TaiRec qr;
    qr.name = tai_it->name;
    qr.seq_pos = tai_it->start;

    TaiRec *tair_1 = stSortedSet_searchLessThanOrEqual(tai, &qr);
    if (tair_1 == NULL ||
        strcmp(tair_1->name, qr.name) != 0) {
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
        
    // we know that the start of our region is somewhere between TaiRec and the next one

    // move to the first record in our file
    LI_seek(li, tair_1->file_pos);
    LI_get_next_line(li);

    fprintf(stderr, "query n=%s sp=%ld\n", qr.name, qr.seq_pos);
    fprintf(stderr, "tair_1 n=%s sp=%ld fp=%ld tags=%s\n", tair_1->name, tair_1->seq_pos, tair_1->file_pos, tair_1->coords);

    if (tair_2) {
        fprintf(stderr, "tair_2 n=%s sp=%ld fp=%ld tags=%s\n", tair_2->name, tair_2->seq_pos, tair_2->file_pos, tair_2->coords);
    } else {
        fprintf(stderr, "tair_2 NULL\n");
    }

    // force taf to start a new alignment at our current file position by making
    // a full coordinates line and reading it
    inject_full_coorindates(tair_1, li);    
    
    // now we have to scan forward until we overlap actually the region
    // TODO: this will surely need speeding up with binary search for giant files...
    //       but i think that's best done once tests are set up based on the simpler version
    Alignment *alignment = NULL;
    Alignment *p_alignment = NULL;
    Alignment *pp_alignment = NULL;
    while((alignment = taf_read_block(p_alignment, false, li)) != NULL) {
        bool in_contig = strcmp(alignment->row->sequence_name, qr.name) == 0;
        fprintf(stderr, "visit alignmenet %s %ld %ld\n", alignment->row->sequence_name, alignment->row->start,
                alignment->row->start + alignment->row->length);
        if (!in_contig ||
            alignment->row->start >= tai_it->end) {
            // we've gone past our query region: there's no hope
            alignment_destruct(alignment);
            if (p_alignment) {
                alignment_destruct(p_alignment);
            }
            if (pp_alignment) {
                alignment_destruct(pp_alignment);
            }
            alignment = NULL;
            p_alignment = NULL;
            pp_alignment = NULL;
            break;
        } else {
            if (pp_alignment != NULL) {
                alignment_destruct(pp_alignment);
            }
            pp_alignment = p_alignment;
            p_alignment = alignment;
            if (in_contig && alignment->row->start < tai_it->end &&
                (alignment->row->start + alignment->row->length) > tai_it->start) {
                // we've found an intersection at "alignment"
                break;
            }
        }
    }
    if (pp_alignment != NULL) {
        alignment_destruct(pp_alignment);
    }
    tai_it->alignment = alignment;
    tai_it->p_alignment = p_alignment;

    if (tai_it->alignment) {
        fprintf(stderr, "settle alignmenet %s %ld %ld\n", alignment->row->sequence_name, alignment->row->start,
                alignment->row->start + alignment->row->length);
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

    fprintf(stderr, "query found at %ld is: %s\n\n", tair_1->file_pos, li->line);

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
        fprintf(stderr, "found left trim %ld of block start=%ld len=%ld vs query start %ld end %ld\n",
                left_trim, aln->row->start, aln->row->length, start, end);
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
            fprintf(stderr, "before left trim cp %ld we have %s start %ld len %ld bases %s\n", cut_point,
                    row->sequence_name, row->start, row->length, row->bases);            
            for (int64_t col = 0; col < cut_point; ++col) {
                if (row->bases[col] != '-') {
                    ++row->start;
                    --row->length;
                }
            }
            char *bases = row->bases;
            row->bases = row->length > 0 ? stString_getSubString(bases, cut_point, strlen(row->bases) - cut_point) : "";
            free(bases);
            fprintf(stderr, "after left trim we have %s start %ld len %ld bases %s\n", row->sequence_name, row->start, row->length, row->bases);

        }
        aln->column_number -= left_trim;
    }

    //clip the right side
    int64_t right_trim = (aln->row->start + aln->row->length) - end;
    if (right_trim > 0) {
        ret = ret | 1;
        fprintf(stderr, "found right trim %ld of block start=%ld len=%ld vs query start %ld end %ld\n",
                right_trim, aln->row->start, aln->row->length, start, end);

        assert(aln->column_number > right_trim);

        // we need to find the cut point by counting off right_trim non-gap bases from the end of the reference row
        int64_t cut_point = strlen(aln->row->bases) - 1;
        for (int64_t cut_count = 0; cut_count < right_trim && cut_point >= 0; --cut_point) {
            fprintf(stderr, "cut count %ld cut point %ld\n", cut_count, cut_point);
            if (aln->row->bases[cut_point] != '-') {
                ++cut_count;
            }
        }
        // then we trim out every row, making sure that we adjust the start/length fields only
        // by non-gap bases we removed
        for (Alignment_Row *row = aln->row; row != NULL; row = row->n_row) {
            fprintf(stderr, "before right trim cp %ld we have %s start %ld len %ld bases %s\n", cut_point,
                    row->sequence_name, row->start, row->length, row->bases);
            int64_t rowlen = strlen(row->bases);
            for (int64_t col = rowlen - 1; col > cut_point; --col) {
                if (row->bases[col] != '-') {
                    --row->length;
                }
            }
            char *bases = row->bases;
            row->bases = row->length > 0 ? stString_getSubString(bases, 0, cut_point) : "";
            free(bases);
            fprintf(stderr, "after right trim we have %s start %ld len %ld bases %s\n", row->sequence_name, row->start, row->length, row->bases);
        }
        aln->column_number -= right_trim;                
    }

    // break the links from the p_alignment to delete rows
    for (Alignment_Row *row = p_aln->row; row != NULL; row = row->n_row) {
        if (row->r_row && row->r_row->length == 0) {
            row->r_row = NULL;
        }
    }
    
    // need to get rid of empty rows as they api doesn't handle them
    Alignment_Row *prev = NULL;
    Alignment_Row *next = NULL;
    for (Alignment_Row *row = aln->row; row != NULL; row = next) {
        next = row->n_row;
        if (row->length == 0) {
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
        tai_it->alignment = taf_read_block(tai_it->p_alignment, false, li);
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
