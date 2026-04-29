#include "block_reader.h"
#include "sonLib.h"

struct _BlockReader {
    LI *li;
    bool is_maf;
    bool run_length_encode_bases;
    Tag *header;  // owned until block_reader_take_header is called
};

BlockReader *block_reader_open(LI *li) {
    int input_format = check_input_format(LI_peek_at_next_line(li));
    if (input_format != 0 && input_format != 1) {
        fprintf(stderr, "Input not supported: unable to detect ##maf or #taf header\n");
        return NULL;
    }
    BlockReader *r = (BlockReader *) st_calloc(1, sizeof(BlockReader));
    r->li = li;
    r->is_maf = (input_format == 1);
    if (r->is_maf) {
        r->header = maf_read_header(li);
        r->run_length_encode_bases = false;
    } else {
        r->header = taf_read_header_2(li, &r->run_length_encode_bases);
    }
    return r;
}

Tag *block_reader_take_header(BlockReader *r) {
    Tag *tag = r->header;
    r->header = NULL;
    return tag;
}

bool block_reader_is_maf(const BlockReader *r) {
    return r->is_maf;
}

bool block_reader_run_length_encoded(const BlockReader *r) {
    return r->run_length_encode_bases;
}

Alignment *block_reader_next(BlockReader *r, Alignment *prev_block) {
    Alignment *a;
    if (r->is_maf) {
        a = maf_read_block(r->li);
        if (a != NULL && prev_block != NULL) {
            alignment_link_adjacent(prev_block, a, 1);
        }
    } else {
        a = taf_read_block(prev_block, r->run_length_encode_bases, r->li);
    }
    return a;
}

void block_reader_destruct(BlockReader *r) {
    if (r == NULL) return;
    if (r->header != NULL) {
        tag_destruct(r->header);
    }
    free(r);
}
