extern "C" {
#include "taf.h"
#include "sonLib.h"
}
#include "bioioC.h"
#include <getopt.h>
#include <time.h>
#ifdef USE_HAL
#include "halBlockViz.h"
#endif

static void add_to_hash(void *fastas, const char *fasta_header, const char *sequence, int64_t length) {
    if(stHash_search((stHash *)fastas, (void *)fasta_header) != NULL) {
        // c++ gives an angry warning if we try to send our string literal directly to st_errAbort, so we do this
        char msg[8192];
        snprintf(msg, 8192, "Found duplicate sequence header: %s\n", fasta_header);
        st_errAbort(msg);
    }
    stHash_insert((stHash *)fastas, stString_copy(fasta_header), stString_copy(sequence));
}

// get a dna interval either from the fastas hash file or from the hal_handle
// note the string returned needs to be freed
static char *get_sequence_fragment(const char* sequence_name, int64_t start, int64_t length, stHash *fastas, int hal_handle, stSet *hal_species) {
    char *fragment = NULL;
    if (fastas) {
        assert(hal_handle == -1);
        char *seq = (char*)stHash_search(fastas, (void*)sequence_name);
        if(seq != NULL) {
            fragment = stString_getSubString(seq, start, length);
        }
    } else {
#ifdef USE_HAL
        assert(fastas == NULL);
        char* species_name = extract_genome_name(sequence_name, hal_species, NULL);
        char* chrom_name = (char*)sequence_name + strlen(species_name) + 1;
        fragment = halGetDna(hal_handle, species_name, chrom_name, start, start + length, NULL);
        free(species_name);
#else
        assert(false);
#endif
    }

    return fragment;
}

void alignment_add_gap_strings(Alignment *p_alignment, Alignment *alignment, stHash *fastas, int hal_handle, stSet *hal_species,
                     int64_t maximum_gap_string_length) {
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        if(row->l_row != NULL && alignment_row_is_predecessor(row->l_row, row)) {
            int64_t gap_length = row->start - (row->l_row->start + row->l_row->length);
            if((maximum_gap_string_length < 0 || gap_length <= maximum_gap_string_length) && row->left_gap_sequence == NULL) {
                char* seq_interval = NULL;
                int64_t i = row->l_row->start + row->l_row->length;
                assert(i >= 0 && i < row->sequence_length);
                if(row->strand) {
                    seq_interval = get_sequence_fragment(row->sequence_name, i, gap_length, fastas, hal_handle, hal_species);
                }
                else { // Case sequence is on the negative strand
                    assert(row->sequence_length - i - gap_length >= 0);
                    char *s = get_sequence_fragment(row->sequence_name,
                                                    row->sequence_length - i - gap_length, gap_length, fastas, hal_handle, hal_species);
                    if(s != NULL) {
                        seq_interval = stString_reverseComplementString(s);
                        free(s);
                    }
                }
                if(seq_interval == NULL) {
                    st_logDebug("[taf] Missing sequence for gap, seq name: %s, skipping!\n", row->sequence_name);
                }
                else {
                    row->left_gap_sequence = seq_interval;
                    assert(strlen(row->left_gap_sequence) == (size_t)gap_length);
                }
            }
        }
        row = row->n_row;
    }
}

stHash *load_sequences_from_fasta_files(char **seq_file_names, int64_t seq_file_number) {
    stHash *fastas = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    for (int64_t i = 0; i < seq_file_number; i++) { //optind < argc) {
        char *seq_file = seq_file_names[i];
        st_logInfo("Parsing sequence file : %s\n", seq_file);
        FILE *fh = fopen(seq_file, "r");
        fastaReadToFunction(fh, fastas, add_to_hash);
        fclose(fh);
    }
    st_logInfo("Finished parsing sequence fasta files\n");
    return fastas;
}

stSet *load_sequences_from_hal_file(char *hal_file, int *hal_handle) {
    stSet *hal_species = NULL;
    st_logInfo("Parsing hal file : %s\n", hal_file);
#ifdef USE_HAL
    *hal_handle = halOpen(hal_file, NULL);
    // make a table of every species in the hal. this will be used to try to parse
    // contig names with lots of dots in them where we don't know which dot separates genome from chrom...
    hal_species = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    struct hal_species_t *species_list = halGetSpecies(*hal_handle, NULL);
    for (struct hal_species_t *species = species_list; species != NULL; species = species->next) {
        stSet_insert(hal_species, stString_copy(species->name));
    }
    halFreeSpeciesList(species_list);
#endif
    st_logInfo("Finished parsing hal file\n");
    return hal_species;
}
