#include "CuTest.h"
#include "taf.h"

static char *make_row_string(Alignment_Row *row) {
    int64_t length = row->length;
    if(row->r_row != NULL) {
        length += row->r_row->length;
        length += row->r_row->start - (row->start + row->length);
    }
    char *row_string = stString_print("%s %i %s %i %i %i", row->sequence_name, (int)row->start,
                    row->strand ? "+" : "-", (int)row->sequence_length, (int)length);
    return row_string;
}

static void test_normalize(CuTest *testCase) {
    // Example maf file
    char *example_file = "./tests/evolverMammals.maf";
    // Read a maf and write a copy of it
    FILE *file = fopen(example_file, "r");
    Alignment *alignment, *p_alignment = NULL;
    bool align_gap_sequences = 1;
    while((alignment = maf_read_block(file)) != NULL) {
        if(p_alignment != NULL) {
            // Align the rows of the block
            alignment_link_adjacent(p_alignment, alignment, 0);

            // And the expected length of the alignment
            int64_t combined_alignment_length = alignment_length(p_alignment) + alignment_length(alignment) + alignment_total_gap_length(p_alignment, align_gap_sequences);

            if((alignment_length(p_alignment) < 50 || alignment_length(alignment) < 50) &&
                alignment_total_gap_length(p_alignment, align_gap_sequences) < 50) {

                // Calculate the number of expected rows and get list of row coordinates
                uint64_t combined_alignment_rows=0;
                Alignment_Row *row = p_alignment->row;
                stList *row_strings = stList_construct3(0, free); // list of strings representing the expected row coordinates
                while(row != NULL) {
                    combined_alignment_rows++;
                    stList_append(row_strings, make_row_string(row));
                    row = row->n_row;
                }
                row = alignment->row;
                while(row != NULL) {
                    if(row->l_row == NULL) {
                        combined_alignment_rows++;
                        stList_append(row_strings, make_row_string(row));
                    }
                    row = row->n_row;
                }

                p_alignment = alignment_merge_adjacent(p_alignment, alignment, align_gap_sequences);

                // Check we have the expected number of rows
                CuAssertIntEquals(testCase, combined_alignment_rows, p_alignment->row_number);

                // Check we have the expected length
                row = p_alignment->row;
                while(row != NULL) {
                    CuAssertIntEquals(testCase, combined_alignment_length, strlen(row->bases));
                    char *row_string = make_row_string(row);
                    bool seen = 0;
                    for(int64_t i=0; i<stList_length(row_strings); i++) {
                        char *row_string_2 = stList_get(row_strings, i);
                        if(strcmp(row_string, row_string_2) == 0) {
                            seen = 1;
                            break;
                        }
                    }
                    CuAssertTrue(testCase, seen);
                    free(row_string);
                    row = row->n_row;
                }

                // Clean up
                stList_destruct(row_strings);
            }
            else {
                alignment_destruct(p_alignment);
                p_alignment = alignment;
            }
        }
        else {
            p_alignment = alignment;
        }
    }
    if(p_alignment != NULL) {
        alignment_destruct(p_alignment);
    }
    fclose(file);
}

static void test_maf_norm_to_maf(CuTest *testCase) {
    /*
     * Run taf_norm with the evolver mammals to output a maf and check command succeeds.
     */
    // Example maf file
    char *example_file = "./tests/evolverMammals.maf";
    char *output_file = "./tests/evolverMammals.maf.norm";
    int i = st_system("./bin/maf_to_taf -i %s | ./bin/taf_add_gap_bases ./tests/seqs/* | ./bin/taf_norm -k -p > %s",
              example_file, output_file);
    CuAssertIntEquals(testCase, 0, i); // return value should be zero
}

static void test_maf_norm(CuTest *testCase) {
    /*
     * Run taf_norm with the evolver mammals files and check command succeeds.
     */
    // Example maf file
    char *example_file = "./tests/evolverMammals.maf";
    char *output_file = "./tests/evolverMammals.taf.norm";
    int i = st_system("./bin/maf_to_taf -i %s | ./bin/taf_add_gap_bases ./tests/seqs/* | ./bin/taf_norm -p > %s",
                      example_file, output_file);
    CuAssertIntEquals(testCase, 0, i); // return value should be zero
}

CuSuite* normalize_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_normalize);
    SUITE_ADD_TEST(suite, test_maf_norm);
    SUITE_ADD_TEST(suite, test_maf_norm_to_maf);
    return suite;
}
