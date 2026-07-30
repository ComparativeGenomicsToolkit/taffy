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
            int64_t combined_alignment_length = alignment_length(p_alignment) + alignment_length(alignment) + alignment_max_gap_length(p_alignment);

            if((alignment_length(p_alignment) < 50 || alignment_length(alignment) < 50) &&
                alignment_max_gap_length(p_alignment) < 50) {

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
    int i = st_system("./bin/taffy view -i %s | ./bin/taffy norm -k -b ./tests/seqs/* > %s",
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
    st_system("rm %s", output_file);
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
    st_system("rm %s", output_file);
}

static void test_dupe_filter(CuTest *testCase) {
    /*
     * Run taf_norm -d on really simple test maf to make sure that it drops a gap-causing dupe
     */
    char *input_file = "./tests/dupe_test.maf";
    char *output_file = "./tests/dupe_test_out.maf";
    int i = st_system("./bin/taffy view -i %s | ./bin/taffy norm -d | ./bin/taffy view -m > %s",
                      input_file, output_file);
    CuAssertIntEquals(testCase, 0, i); // return value should be zero        
    char *truth_file = "./tests/dupe_test_truth.maf";
    int diff_ret = st_system("diff %s %s", output_file, truth_file);
    CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames
    st_system("rm %s", output_file);
}

static void test_norm_pipeline(CuTest *testCase) {
    /*
     * Run taf_norm -d on really simple test maf to make sure that it drops a gap-causing dupe
     */
    char *input_file = "./tests/evolverMammals.maf.mini";
    char *filter_file = "./tests/filter_file.txt";
    char *output_file = "./tests/evolverMammals.maf.norm_pipeline";
    int i = st_system("./bin/taffy view -i %s | ./bin/taffy sort -f %s -r | ./bin/taffy norm -b ./tests/seqs/* -k > %s",
                      input_file, filter_file, output_file);
    CuAssertIntEquals(testCase, 0, i); // return value should be zero
    char *truth_file = "./tests/evolverMammals.maf.mini.filtered.normalized";
    int diff_ret = st_system("diff %s %s", output_file, truth_file);
    CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames
    st_system("rm %s", output_file);
}

// Verifies that taffy add-gap-bases -i x.maf produces output identical to
// taffy view -i x.maf | taffy add-gap-bases. add-gap-bases relies on adjacent
// blocks being linked (alignment_add_gap_strings reads l_row/r_row coordinates),
// so this is the most sensitive of the dual-input migration tests.
static void test_add_gap_bases_maf_input(CuTest *testCase) {
    char *example_file = "./tests/evolverMammals.maf";
    char *piped_out = "./tests/agb.piped.taf";
    char *direct_out = "./tests/agb.direct.taf";
    int i = st_system("./bin/taffy view -i %s | ./bin/taffy add-gap-bases ./tests/seqs/* > %s",
                      example_file, piped_out);
    CuAssertIntEquals(testCase, 0, i);
    int j = st_system("./bin/taffy add-gap-bases -i %s -o %s ./tests/seqs/*",
                      example_file, direct_out);
    CuAssertIntEquals(testCase, 0, j);
    int diff_ret = st_system("diff %s %s", piped_out, direct_out);
    CuAssertIntEquals(testCase, 0, diff_ret);
    st_system("rm -f %s %s", piped_out, direct_out);
}

// Verifies that taffy norm -i x.maf produces output identical to
// taffy view -i x.maf | taffy norm. This is the per-tool dual-input
// regression test for the BlockReader migration.
static void test_norm_maf_input(CuTest *testCase) {
    char *example_file = "./tests/evolverMammals.maf";
    char *piped_out = "./tests/norm.piped.taf";
    char *direct_out = "./tests/norm.direct.taf";
    int i = st_system("./bin/taffy view -i %s | ./bin/taffy norm > %s", example_file, piped_out);
    CuAssertIntEquals(testCase, 0, i);
    int j = st_system("./bin/taffy norm -i %s -o %s", example_file, direct_out);
    CuAssertIntEquals(testCase, 0, j);
    int diff_ret = st_system("diff %s %s", piped_out, direct_out);
    CuAssertIntEquals(testCase, 0, diff_ret);
    st_system("rm -f %s %s", piped_out, direct_out);
}

// Checks alignment_remove_all_gap_columns directly: the gap-only columns go, the remaining bases
// keep their order, the column tags follow the columns they belong to, and the row coordinates are
// left alone.
static void test_remove_all_gap_columns(CuTest *testCase) {
    char *example_file = "./tests/gap_column_test.maf";
    FILE *file = fopen(example_file, "r");
    LI *li = LI_construct(file);
    Tag *header = maf_read_header(li);
    tag_destruct(header);

    // The first block of the toy file is 8 columns of which 4 (the first, one in the middle and the
    // last two) are entirely gaps
    Alignment *alignment = maf_read_block(li);
    CuAssertTrue(testCase, alignment != NULL);
    CuAssertIntEquals(testCase, 8, alignment->column_number);
    CuAssertIntEquals(testCase, 3, alignment->row_number);

    // Tag a column that survives and a column that does not, so we can check the tags are remapped
    alignment->column_tags[1] = tag_construct("keep", "1", alignment->column_tags[1]);
    alignment->column_tags[3] = tag_construct("drop", "1", alignment->column_tags[3]);

    // Record the row coordinates so we can check they are untouched
    int64_t starts[3], lengths[3], sequence_lengths[3];
    int64_t i = 0;
    for(Alignment_Row *row = alignment->row; row != NULL; row = row->n_row, i++) {
        starts[i] = row->start;
        lengths[i] = row->length;
        sequence_lengths[i] = row->sequence_length;
    }
    CuAssertIntEquals(testCase, 3, i);

    CuAssertIntEquals(testCase, 4, alignment_remove_all_gap_columns(alignment));
    CuAssertIntEquals(testCase, 4, alignment->column_number);
    CuAssertIntEquals(testCase, 3, alignment->row_number); // No rows are removed

    i = 0;
    for(Alignment_Row *row = alignment->row; row != NULL; row = row->n_row, i++) {
        CuAssertStrEquals(testCase, "ACGT", row->bases);
        CuAssertIntEquals(testCase, starts[i], row->start);
        CuAssertIntEquals(testCase, lengths[i], row->length);
        CuAssertIntEquals(testCase, sequence_lengths[i], row->sequence_length);
    }
    CuAssertIntEquals(testCase, 3, i);

    // The tag of the kept column moved with it from index 1 to index 0, and the tag of the removed
    // column is gone
    CuAssertTrue(testCase, tag_find(alignment->column_tags[0], "keep") != NULL);
    for(int64_t j=0; j<alignment->column_number; j++) {
        CuAssertTrue(testCase, tag_find(alignment->column_tags[j], "drop") == NULL);
    }

    // Running it again is a no-op now that no column is all gaps
    CuAssertIntEquals(testCase, 0, alignment_remove_all_gap_columns(alignment));
    CuAssertIntEquals(testCase, 4, alignment->column_number);

    alignment_destruct(alignment, 1);
    LI_destruct(li);
    fclose(file);
}

// A block in which every column is a gap keeps one column rather than being reduced to nothing. A
// zero column block can not be written, and dropping such a block from the alignment instead would
// strand the interstitial gap sequence that the following block records relative to it.
static void test_remove_all_gap_columns_all_gap_block(CuTest *testCase) {
    Alignment *alignment = st_calloc(1, sizeof(Alignment));
    alignment->column_number = 3;
    alignment->column_tags = st_calloc(3, sizeof(Tag *));
    alignment->row_number = 2;
    Alignment_Row *row_1 = st_calloc(1, sizeof(Alignment_Row));
    row_1->sequence_name = stString_copy("simCow.chr1");
    row_1->bases = stString_copy("---");
    row_1->strand = 1;
    row_1->sequence_length = 100;
    Alignment_Row *row_2 = st_calloc(1, sizeof(Alignment_Row));
    row_2->sequence_name = stString_copy("simDog.chr1");
    row_2->bases = stString_copy("---");
    row_2->strand = 1;
    row_2->sequence_length = 100;
    alignment->row = row_1;
    row_1->n_row = row_2;

    CuAssertIntEquals(testCase, 2, alignment_remove_all_gap_columns(alignment));
    CuAssertIntEquals(testCase, 1, alignment->column_number);
    CuAssertStrEquals(testCase, "-", row_1->bases);
    CuAssertStrEquals(testCase, "-", row_2->bases);

    // And it stays at one column if run again
    CuAssertIntEquals(testCase, 0, alignment_remove_all_gap_columns(alignment));
    CuAssertIntEquals(testCase, 1, alignment->column_number);

    alignment_destruct(alignment, 1);
}

/*
 * A block that is nothing but gaps has to survive as a block, because the block after it records
 * its interstitial gap sequence (a TAF "G" record) relative to it. Removing the all-gap block
 * would leave that gap sequence describing a gap that no longer exists, which silently drops the
 * unaligned bases and leaves each row's length field disagreeing with its aligned bases.
 */
static void test_all_gap_block_keeps_gap_sequence(CuTest *testCase) {
    char *input_file = "./tests/gap_column_gap_seq_test.taf";
    char *output_file = "./tests/gap_column_gap_seq_out.maf";
    int i = st_system("./bin/taffy norm -i %s -k -o %s", input_file, output_file);
    CuAssertIntEquals(testCase, 0, i); // must not abort

    // The two blocks either side of the all-gap block merge, so the 5 bases of gap sequence before
    // the all-gap block and the 3 before the final block are both part of the merged row
    FILE *file = fopen(output_file, "r");
    CuAssertTrue(testCase, file != NULL);
    LI *li = LI_construct(file);
    Tag *header = maf_read_header(li);
    tag_destruct(header);
    Alignment *alignment = maf_read_block(li);
    CuAssertTrue(testCase, alignment != NULL);
    for(Alignment_Row *row = alignment->row; row != NULL; row = row->n_row) {
        int64_t bases = 0;
        for(int64_t j=0; j<alignment->column_number; j++) {
            if(row->bases[j] != '-') {
                bases++;
            }
        }
        // 2 aligned + 5 unaligned + 3 unaligned + 2 aligned
        CuAssertIntEquals(testCase, 12, row->length);
        CuAssertIntEquals(testCase, row->length, bases); // the gap bases must still be there
    }
    alignment_destruct(alignment, 1);
    LI_destruct(li);
    fclose(file);
    st_system("rm %s", output_file);
}

static void test_gap_column_filter(CuTest *testCase) {
    /*
     * Run taffy norm on a small maf containing all-gap columns and check they are removed. The toy
     * file also exercises the interaction with block merging: its second and third blocks are
     * adjacent on the same contig and get merged, and the merged block must not carry over the gap
     * columns of either input block.
     */
    char *input_file = "./tests/gap_column_test.maf";
    char *taf_file = "./tests/gap_column_test_out.taf";
    char *output_file = "./tests/gap_column_test_out.maf";
    // Run norm as its own command rather than in the middle of a pipeline, so that a norm failure
    // is not masked by the exit status of the last stage
    int i = st_system("./bin/taffy view -i %s -o %s", input_file, taf_file);
    CuAssertIntEquals(testCase, 0, i);
    int j = st_system("./bin/taffy norm -i %s -k -o %s", taf_file, output_file);
    CuAssertIntEquals(testCase, 0, j); // return value should be zero
    char *truth_file = "./tests/gap_column_test_truth.maf";
    int diff_ret = st_system("diff %s %s", output_file, truth_file);
    CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files same
    st_system("rm %s %s", output_file, taf_file);
}

// Counts the columns of a maf file that consist entirely of gaps, additionally checking that each
// row's stated length agrees with its number of non-gap bases. The latter is the invariant that
// removing columns has to preserve: it breaks if we ever drop a column that held a base.
static int64_t count_all_gap_columns(CuTest *testCase, char *maf_file) {
    FILE *file = fopen(maf_file, "r");
    CuAssertTrue(testCase, file != NULL);
    LI *li = LI_construct(file);
    Tag *header = maf_read_header(li);
    tag_destruct(header);
    int64_t all_gap_columns = 0;
    Alignment *alignment;
    while((alignment = maf_read_block(li)) != NULL) {
        for(int64_t i=0; i<alignment->column_number; i++) {
            bool all_gaps = true;
            for(Alignment_Row *row = alignment->row; row != NULL; row = row->n_row) {
                if(row->bases[i] != '-') {
                    all_gaps = false;
                    break;
                }
            }
            if(all_gaps) {
                all_gap_columns++;
            }
        }
        for(Alignment_Row *row = alignment->row; row != NULL; row = row->n_row) {
            int64_t bases = 0;
            for(int64_t i=0; i<alignment->column_number; i++) {
                if(row->bases[i] != '-') {
                    bases++;
                }
            }
            CuAssertIntEquals(testCase, row->length, bases);
        }
        alignment_destruct(alignment, 1);
    }
    LI_destruct(li);
    fclose(file);
    return all_gap_columns;
}

// The motivating case: filtering out a species leaves behind the columns that only it had a base
// in, and norm should clean them up. The toy input gives simMouse two insertions relative to the
// other two rows, so filtering simMouse out strands four all-gap columns.
static void test_gap_columns_after_filtering_a_species(CuTest *testCase) {
    char *input_file = "./tests/gap_column_filter_test.maf";
    char *filter_file = "./tests/gap_column_filter.txt";
    char *filtered_file = "./tests/gap_column_filtered.maf";
    char *normalized_file = "./tests/gap_column_filtered.norm.maf";

    int i = st_system("./bin/taffy view -i %s | ./bin/taffy sort -f %s | ./bin/taffy view -m > %s",
                      input_file, filter_file, filtered_file);
    CuAssertIntEquals(testCase, 0, i);
    int j = st_system("./bin/taffy view -i %s | ./bin/taffy sort -f %s | ./bin/taffy norm -k > %s",
                      input_file, filter_file, normalized_file);
    CuAssertIntEquals(testCase, 0, j);

    // Filtering on its own strands the four columns simMouse was the only row with a base in
    CuAssertIntEquals(testCase, 4, count_all_gap_columns(testCase, filtered_file));
    // Normalizing afterwards removes all of them
    CuAssertIntEquals(testCase, 0, count_all_gap_columns(testCase, normalized_file));

    char *truth_file = "./tests/gap_column_filter_test_truth.maf";
    int diff_ret = st_system("diff %s %s", normalized_file, truth_file);
    CuAssertIntEquals(testCase, 0, diff_ret);

    st_system("rm -f %s %s", filtered_file, normalized_file);
}

/*
 * taffy add-gap-bases records the unaligned sequence between two blocks as a TAF "G" record, whose
 * length is tied to the gap between a row and the row it continues from. taffy norm re-links rows,
 * which can pair a row with a different, or a differently distant, predecessor, so a gap sequence
 * that no longer fits its gap must not be spliced into a merged row: doing so puts more bases in the
 * row than its length field accounts for. This is the pipeline cactus uses, via the HAL backend.
 */
static void test_norm_with_gap_sequences(CuTest *testCase) {
    char *gap_file = "./tests/evolverMammals.gapseq.taf";
    char *output_file = "./tests/evolverMammals.gapseq.norm.maf";
    int i = st_system("./bin/taffy view -i ./tests/evolverMammals.maf | "
                      "./bin/taffy add-gap-bases ./tests/seqs/* > %s", gap_file);
    CuAssertIntEquals(testCase, 0, i);
    int j = st_system("./bin/taffy norm -i %s -k -o %s", gap_file, output_file);
    CuAssertIntEquals(testCase, 0, j);
    // count_all_gap_columns checks each row's length field against its non-gap base count, which is
    // the invariant a stale gap sequence breaks
    count_all_gap_columns(testCase, output_file);

    // And again with the sequences available, where the gap sequences that no longer fit should be
    // refilled with real bases rather than left as Ns
    int k = st_system("./bin/taffy norm -i %s -k -b ./tests/seqs/* -o %s", gap_file, output_file);
    CuAssertIntEquals(testCase, 0, k);
    count_all_gap_columns(testCase, output_file);

    st_system("rm -f %s %s", gap_file, output_file);
}

/*
 * Merging two blocks moves the right block's column tags into the merged block, so the right block
 * must not then free them: reading them back when the merged block is written is a use after free,
 * which crashed taffy norm on any taf carrying column tags on a mergeable block. The input also has
 * an all-gap column with a tag on it, so this covers the tags being remapped by both the merge and
 * the gap column removal at once.
 */
static void test_column_tags_through_merge(CuTest *testCase) {
    char *input_file = "./tests/column_tag_merge_test.taf";
    char *output_file = "./tests/column_tag_merge_out.taf";
    int i = st_system("./bin/taffy norm -i %s -o %s", input_file, output_file);
    CuAssertIntEquals(testCase, 0, i); // must not crash
    char *truth_file = "./tests/column_tag_merge_test_truth.taf";
    int diff_ret = st_system("diff %s %s", output_file, truth_file);
    CuAssertIntEquals(testCase, 0, diff_ret);
    st_system("rm -f %s", output_file);
}

CuSuite* normalize_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_remove_all_gap_columns);
    SUITE_ADD_TEST(suite, test_remove_all_gap_columns_all_gap_block);
    SUITE_ADD_TEST(suite, test_all_gap_block_keeps_gap_sequence);
    SUITE_ADD_TEST(suite, test_gap_column_filter);
    SUITE_ADD_TEST(suite, test_gap_columns_after_filtering_a_species);
    SUITE_ADD_TEST(suite, test_norm_with_gap_sequences);
    SUITE_ADD_TEST(suite, test_column_tags_through_merge);
    SUITE_ADD_TEST(suite, test_normalize);
    SUITE_ADD_TEST(suite, test_maf_norm);
    SUITE_ADD_TEST(suite, test_maf_norm_to_maf);
    SUITE_ADD_TEST(suite, test_dupe_filter);
    SUITE_ADD_TEST(suite, test_norm_pipeline);
    SUITE_ADD_TEST(suite, test_norm_maf_input);
    SUITE_ADD_TEST(suite, test_add_gap_bases_maf_input);
    return suite;
}
