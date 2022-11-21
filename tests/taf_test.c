#include "CuTest.h"
#include "taf.h"
#include "sonLib.h"

static char *get_random_tags() {
    stList *tags = stList_construct3(0, free);
    while(st_random() > 0.5) {
        stList_append(tags, stString_print("%f:%f", st_random(), st_random()));
    }
    char *tag_string = stString_join2(" ", tags);
    stList_destruct(tags);
    return tag_string;
}

Tag *parse_tags(stList *tokens, int64_t starting_token, char *delimiter);

void check_tags(CuTest *testCase, Tag *t, Tag *j) {
    while(t != NULL) {
        CuAssertTrue(testCase, j != NULL);
        CuAssertStrEquals(testCase, t->key, j->key);
        CuAssertStrEquals(testCase, t->value, j->value);
        t = t->n_tag;
        j = j->n_tag;
    }
    CuAssertTrue(testCase, j == NULL);
}

static void test_taf(CuTest *testCase) {
    // Example maf file
    char *example_file = "./tests/evolverMammals.maf"; // "./tests/chr2_KI270776v1_alt.maf.1"; // "./tests/chr2_KI270893v1_alt.maf"; //"./tests/chr2_KI270776v1_alt.maf";
    char *temp_copy = "./tests/evolverMammals.taf"; //"./tests/chr2_KI270776v1_alt.taf.1"; // "./tests/chr2_KI270893v1_alt.taf"; //"./tests/chr2_KI270776v1_alt.taf";
    bool run_length_encode_bases = 0;
    stList *column_tags = stList_construct3(0, (void (*)(void *))tag_destruct);

    // Write out the taf file
    FILE *file = fopen(example_file, "r");
    FILE *out_file = fopen(temp_copy, "w");
    Alignment *alignment, *p_alignment = NULL;
    while((alignment = maf_read_block(file)) != NULL) {
        if(p_alignment != NULL) {
            alignment_link_adjacent(p_alignment, alignment, 1);
        }

        // Make random tags for each column
        alignment->column_tags = st_malloc(sizeof(Tag *) * alignment->column_number);
        for(int64_t i=0; i<alignment->column_number; i++) {
            char *tag_string = get_random_tags();
            stList *tokens = stString_split(tag_string);
            alignment->column_tags[i] = parse_tags(tokens, 0, ":");
            stList_append(column_tags, parse_tags(tokens, 0, ":"));
            stList_destruct(tokens);
        }

        taf_write_block(p_alignment, alignment, run_length_encode_bases, 1000, out_file);
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment);
        }
        p_alignment = alignment;
    }
    if(p_alignment != NULL) {
        alignment_destruct(p_alignment);
    }
    fclose(file);
    fclose(out_file);

    // Now parse the taf
    file = fopen(example_file, "r");
    FILE *file_copy = fopen(temp_copy, "r");
    LI *li = LI_construct(file_copy);
    Alignment *alignment2 = NULL;
    Alignment *p_alignment2 = NULL;
    int64_t column_index = 0;
    while((alignment = maf_read_block(file)) != NULL) {
        // Alignment *taf_read_block(Alignment *p_block, bool run_length_encode_bases, LI *li)
        alignment2 = taf_read_block(p_alignment2, run_length_encode_bases, li);
        CuAssertTrue(testCase, alignment2 != NULL);
        // Check that the blocks are the same
        CuAssertIntEquals(testCase, alignment->row_number, alignment2->row_number);
        CuAssertIntEquals(testCase, alignment->column_number, alignment2->column_number);

        // Check the rows
        Alignment_Row *row = alignment->row, *row2 = alignment2->row;
        while(row != NULL) {
            CuAssertTrue(testCase, row2 != NULL);
            CuAssertStrEquals(testCase, row->sequence_name, row2->sequence_name);
            CuAssertIntEquals(testCase, row->start, row2->start);
            CuAssertIntEquals(testCase, row->length, row2->length);
            CuAssertIntEquals(testCase, row->sequence_length, row2->sequence_length);
            CuAssertIntEquals(testCase, row->strand, row2->strand);
            CuAssertStrEquals(testCase, row->bases, row2->bases);
            row = row->n_row; row2 = row2->n_row;
        }
        CuAssertTrue(testCase, row2 == NULL);

        // Check the tags are the same
        for(int64_t i=0; i<alignment2->column_number; i++) {
            check_tags(testCase, alignment2->column_tags[i], stList_get(column_tags, column_index++));
        }

        alignment_destruct(alignment);
        if(p_alignment2 != NULL) {
            alignment_destruct(p_alignment2);
        }
        p_alignment2 = alignment2;
    }
    if(p_alignment2 != NULL) {
        alignment_destruct(p_alignment2);
    }
    CuAssertTrue(testCase, taf_read_block(p_alignment2, run_length_encode_bases, li) == NULL);
    CuAssertIntEquals(testCase, stList_length(column_tags), column_index);
    // clean up
    stList_destruct(column_tags);
    LI_destruct(li);
    fclose(file);
}

CuSuite* taf_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_taf);
    return suite;
}
