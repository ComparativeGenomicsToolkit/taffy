#include "CuTest.h"
#include "taf.h"
#include "sonLib.h"
#include "bioioC.h"

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
    // Read a maf and normalize it
    FILE *file = fopen(example_file, "r");
    LI *li = LI_construct(file);
    Alignment *alignment, *p_alignment = NULL;
    while((alignment = maf_read_block(li)) != NULL) {
        if(p_alignment != NULL) {
            // Align the rows of the block
            alignment_link_adjacent(p_alignment, alignment, 0);

            // And the expected length of the alignment
            int64_t combined_alignment_length = alignment_length(p_alignment) + alignment_length(alignment) + alignment_total_gap_length(p_alignment);

            if((alignment_length(p_alignment) < 50 || alignment_length(alignment) < 50) &&
                alignment_total_gap_length(p_alignment) < 50) {

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

                p_alignment = alignment_merge_adjacent(p_alignment, alignment);

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
                alignment_destruct(p_alignment, 1);
                p_alignment = alignment;
            }
        }
        else {
            p_alignment = alignment;
        }
    }
    if(p_alignment != NULL) {
        alignment_destruct(p_alignment, 1);
    }
    fclose(file);
    LI_destruct(li);
}

void add_to_hash(void *fastas, const char *fasta_header, const char *sequence, int64_t length) {
    if(stHash_search((stHash *)fastas, (void *)fasta_header) != NULL) {
        // c++ gives an angry warning if we try to send our string literal directly to st_errAbort, so we do this
        char msg[8192];
        sprintf(msg, "Found duplicate sequence header: %s\n", fasta_header);
        st_errAbort(msg);
    }
    stHash_insert((stHash *)fastas, stString_copy(fasta_header), stString_copy(sequence));
}

static void test_maf_norm_to_maf(CuTest *testCase) {
    /*
     * Run taf_norm with the evolver mammals to output a maf and check command succeeds, then
     * check bases in alignment correspond to actual sequence in input string files
     */
    // Example maf file
    char *example_file = "./tests/evolverMammals.maf";
    char *output_file = "./tests/evolverMammals.maf.norm";
    int i = st_system("./bin/taffy view -i %s | ./bin/taffy add-gap-bases ./tests/seqs/* | ./bin/taffy norm -k > %s",
              example_file, output_file);
    CuAssertIntEquals(testCase, 0, i); // return value should be zero

    // Parse the sequence files
    stHash *fastas = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    char *seqs[5] = { "./tests/seqs/simCow.chr6", "./tests/seqs/simDog.chr6",
                     "./tests/seqs/simHuman.chr6", "./tests/seqs/simMouse.chr6", "./tests/seqs/simRat.chr6" };
    for(int64_t i=0; i<5; i++) {
        st_logInfo("Parsing sequence file : %s\n", seqs[i]);
        FILE *fh = fopen(seqs[i], "r");
        fastaReadToFunction(fh, fastas, add_to_hash);
        fclose(fh);
    }

    // Now load the maf and check that the bases match the sequences
    FILE *file = fopen(example_file, "r");
    Alignment *alignment;
    LI *li = LI_construct(file);
    while((alignment = maf_read_block(li)) != NULL) {
        Alignment_Row *row = alignment->row;
        while(row != NULL) {
            // Get the sequence
            char *string = stHash_search(fastas, row->sequence_name);
            if(string != NULL) {
                // Get the sequence of bases for each row
                int64_t j = row->start;
                for (int64_t i = 0; i < row->length; i++) {
                    if (row->bases[i] != '-') {  // If not a gap
                        if(row->strand) { // Case on positive strand
                            CuAssertIntEquals(testCase, row->bases[i], string[j++]);
                        }
                        else { // Case on reverse strand
                            CuAssertIntEquals(testCase, row->bases[i], stString_reverseComplementChar(string[row->sequence_length - 1 - j++]));
                        }
                    }
                }
            }
            else { // Check is an ancestor sequence
                st_logDebug("Didn't find sequence for: %s\n", row->sequence_name);
            }
            row = row->n_row;
        }
        alignment_destruct(alignment, 1);
    }

    // Cleanup
    stHash_destruct(fastas);
    fclose(file);
    LI_destruct(li);
}

static void test_maf_norm(CuTest *testCase) {
    /*
     * Run taf_norm with the evolver mammals files and check command succeeds.
     */
    // Example maf file
    char *example_file = "./tests/evolverMammals.maf";
    char *output_file = "./tests/evolverMammals.taf.norm";
    int i = st_system("./bin/taffy view -i %s | ./bin/taffy add-gap-bases ./tests/seqs/* | ./bin/taffy norm > %s",
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
