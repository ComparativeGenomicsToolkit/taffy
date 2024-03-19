/*
 * taf_coverage: Compute some basic coverage statistics from a TAF/MAF file
 *
 *  Released under the MIT license, see LICENSE.txt
*/

extern "C" {
#include "taf.h"
#include "sonLib.h"
}
#include <getopt.h>
#include <time.h>
#include <unordered_map>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <functional>
#include <iostream>
#include <iomanip>
#include <limits>

using namespace std;

// keep track of very basic coverage stats, broken down
// into total and regions with 1-1 alignments (single)
struct CoverageCounts {
    int64_t tot_aligned = 0;
    int64_t tot_identical = 0;
    int64_t single_aligned = 0;
    int64_t single_identical = 0;
    int64_t prev_ref_pos = 0;
    map<int64_t, int64_t> gap_hist;
};
// coverage stats for a given reference contig
struct CoverageMap {
    int64_t ref_length = -1;
    map<string, CoverageCounts> genome_map;
    CoverageCounts& operator[](const string& key) { return genome_map[key]; }    
};
// map ref contig name to query genome contig coverage counts
typedef map<string, CoverageMap> ContigCoverageMap;

// update the coverage map for a given block
static void update_block_coverage(Alignment* aln, Alignment* prev_aln, const string& ref_name,
                                  stHash* genome_names, ContigCoverageMap& contig_cov_map);
// sum up all the coverages and add a total coverage entry in the map
static void update_total_coverage(ContigCoverageMap& contig_cov_map, const string& key = "_Total_");
// add the final gap in each ref contig and 
static void add_final_gap(ContigCoverageMap& contig_cov_map);
// transform so gaps counts are cumulative
static void postprocess_gap_hist(ContigCoverageMap& contig_cov_map);
// print the coverage tsv
static void print_coverage_tsv(const ContigCoverageMap& contig_cov_map, const set<int64_t>& gap_thresholds, ostream& os);

static void usage() {
    fprintf(stderr, "taffy coverage [options]\n");    
    fprintf(stderr, "Compute very basic pairwise coverage stats as fraction and bp for a TAF file\n");
    fprintf(stderr, "-i --inputFile : Input taf file to normalize. If not specified reads from stdin\n");
    fprintf(stderr, "-r --reference : Name of reference genome. If note specified used first row in block\n");
    fprintf(stderr, "-g --genomeNames : List of genome names (quoted, space-separated), ex from \"$(halStats --genomes aln.hal)\". This can help contig name parsing which otherwise uses everything up to first . as genome name\n");
    fprintf(stderr, "-a, --gapThreshold : Breakdown rows using given gap threshold, to restrict aligned bp to exclude gaps>threshold. Multiple allowed. \n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int taf_coverage_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    string reference;
    char *genomeNames = NULL;
    set<int64_t> gap_thresholds = {-1};

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "reference", required_argument, 0, 'r' },
                                                { "genomeNames", required_argument, 0, 'g' },
                                                { "gapThreshold", required_argument, 0, 'a' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:r:g:a:", long_options, &option_index);
        if (key == -1) {
            break;
        }
        
        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'i':
                inputFile = optarg;
                break;
            case 'r':
                reference = optarg;
                break;
            case 'g':
                genomeNames = optarg;
                break;
            case 'a':
                gap_thresholds.insert(atol(optarg));
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    if (optind != argc) {
        usage();
        return 1;
    }
    
    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Input file string : %s\n", inputFile);
    if (!reference.empty()) {
        st_logInfo("Reference : %s\n", reference.c_str());
    }
    if (genomeNames) {
        st_logInfo("Genome names : %s\n", genomeNames);
    }

    // per-genome results collected here
    ContigCoverageMap contig_coverage_map;
    
    // load the given genome names into a stHash (since that's what the existing name parser machinery wants)
    // values don't matter, just keys...
    stHash* genome_names_hash = NULL;
    if (genomeNames != NULL) {
        stList* tokens = stString_splitByString(genomeNames, " ");
        genome_names_hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
        for (int64_t i = 0; i < stList_length(tokens); ++i) {
            stHash_insert(genome_names_hash, stString_copy((char*)stList_get(tokens, i)),
                          stString_copy((char*)stList_get(tokens, i)));
        }
        stList_destruct(tokens);
    }

    // Open TAF    
    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    LI *li = LI_construct(input);

    // Parse the header
    bool run_length_encode_bases;
    Tag *tag = taf_read_header_2(li, &run_length_encode_bases);
    tag_destruct(tag);
    
    Alignment *alignment, *p_alignment = NULL;
    while((alignment = taf_read_block(p_alignment, run_length_encode_bases, li)) != NULL) {
        // update the coverage
        update_block_coverage(alignment, p_alignment, reference, genome_names_hash, contig_coverage_map);

        // Clean up the previous alignment
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment, 1);
        }
        p_alignment = alignment; // Update the previous alignment
    }
    if(p_alignment != NULL) { // Clean up the final alignment
        alignment_destruct(p_alignment, 1);
    }

    // add gaps from last covered base to ends of contigs
    add_final_gap(contig_coverage_map);
    
    // total up and print coverage
    update_total_coverage(contig_coverage_map);

    // finalize the gap coverage, making it cumulative in bp
    postprocess_gap_hist(contig_coverage_map);

    // write the table to stdout
    print_coverage_tsv(contig_coverage_map, gap_thresholds, cout);    

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    LI_destruct(li);
    if(inputFile != NULL) {
        fclose(input);
    }
    if (genome_names_hash != NULL) {
        stHash_destruct(genome_names_hash);
    }

    st_logInfo("taffy coverage is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}


void update_block_coverage(Alignment* aln, Alignment* prev_aln, const string& ref_name, stHash* genome_names,
                           ContigCoverageMap& contig_cov_map) {
    
    // random access rows
    vector<Alignment_Row*> rows(aln->row_number, NULL);

    // remember parsed names
    vector<string> row_to_name(aln->row_number);

    // group rows by genome
    unordered_map<string, int64_t> name_to_group;
    vector<vector<int64_t>> groups;
    vector<int64_t> row_to_group(aln->row_number, -1);

    // reference row
    int64_t ref_row_idx = -1;

    // index rows and group them by genome name, also find the reference row
    int64_t row_idx = 0;
    for (Alignment_Row* row = aln->row; row != NULL; row = row->n_row, ++row_idx) {
        // resolve the genome name from the full sequence name
        char* name = NULL;
        // check the input list if given
        if (genome_names != NULL) {
            name = extract_genome_name(row->sequence_name, NULL, genome_names);
        }
        row_to_name[row_idx] = name != NULL ? name : row->sequence_name;
        // if the name wasn't in the list, try parsing on first .
        if (name == NULL) {
            auto dotpos = row_to_name[row_idx].find('.');
            if (dotpos > 0 && dotpos != string::npos) {
                row_to_name[row_idx] = row_to_name[row_idx].substr(0, dotpos);
            }                
        }
     
        const string& name_str = row_to_name[row_idx];        
        free(name);

        rows[row_idx] = row;
        
        // update the ref row
        if (ref_row_idx == -1 && (ref_name.empty() || name_str == ref_name)) {
            ref_row_idx = row_idx;
        }

        // update the groups
        int64_t group_idx = -1;
        if (name_to_group.count(name_str)) {
            group_idx = name_to_group[name_str];
        } else {
            group_idx = name_to_group.size();
            name_to_group[name_str] = group_idx;
            groups.push_back({});            
        }
        groups[group_idx].push_back(row_idx);
        row_to_group[row_idx] = group_idx;
    }

    // we ignore blocks with no reference.  todo: should there be a warning?
    if (ref_row_idx == -1) {
        return;
    }
    
    // find / initialize the coverage data structure (can only be done after find ref row)
    CoverageMap& cov_map = contig_cov_map[rows[ref_row_idx]->sequence_name];
    if (cov_map.ref_length < 0) {
        cov_map.ref_length = rows[ref_row_idx]->sequence_length;
    }

    //link rows to their coverage counters
    vector<CoverageCounts*> row_to_cov(aln->row_number, NULL);
    for (row_idx = 0; row_idx < aln->row_number; ++row_idx) {
        row_to_cov[row_idx] = &cov_map[row_to_name[row_idx]];
    }

    int64_t ref_count = groups[row_to_group[ref_row_idx]].size();
    int64_t ref_pos = rows[ref_row_idx]->start;
    
    // update the coverage column by column
    for (int64_t col = 0; col < aln->column_number; ++col) {
        char ref_base = toupper(rows[ref_row_idx]->bases[col]);
        if (ref_base != '-' && ref_base != 'N') {
            for (const auto& group : groups) {
                CoverageCounts& coverage = *row_to_cov[group.front()];
                bool found_aligned = false;
                bool found_identical = false;
                for (int64_t row_idx : group) {
                    char alt_base = toupper(rows[row_idx]->bases[col]);
                    if (alt_base != '-' && alt_base != 'N') {
                        if (!found_aligned) {
                            ++coverage.tot_aligned;
                            if (ref_count == 1 && group.size() == 1) {
                                ++coverage.single_aligned;
                            }
                            found_aligned = true;                            
                        }
                        if (!found_identical && ref_base == alt_base) {
                            ++coverage.tot_identical;
                            if (ref_count == 1 && group.size() == 1) {
                                ++coverage.single_identical;
                            };
                            found_identical = true;
                        }
                        // update gap information for given species
                        int64_t gap_len = ref_pos - coverage.prev_ref_pos - 1;
                        if (gap_len > 0) {
                            coverage.gap_hist[gap_len] += 1;
                        }
                        coverage.prev_ref_pos = ref_pos;
                    }
                    if (found_aligned && found_identical) {
                        break;
                    }
                }
            }
        }
        if (ref_base != '-') {
            ++ref_pos;
        }
    }
}

void update_total_coverage(ContigCoverageMap& contig_cov_map, const string& key) {
    string fixed_key = key;
    if (contig_cov_map.count(key)) {
        string new_key = key + "_";
        while (contig_cov_map.count(new_key)) {
            new_key += "_";            
        }
        cerr << "[taffy coverage] Warning: Total coverage stored as \"" << new_key << "\" because \"" << key << "\" was in map" << endl;
        fixed_key = new_key;
    }

    CoverageMap& tot_cov = contig_cov_map[fixed_key];
    tot_cov.ref_length = 0;
    for (const auto& contig_covmap : contig_cov_map) {
        if (contig_covmap.first == fixed_key) {
            continue;
        }
        assert(contig_covmap.second.ref_length >= 0);
        tot_cov.ref_length += contig_covmap.second.ref_length;
        for (const auto& genome_counts : contig_covmap.second.genome_map) {
            CoverageCounts& tot_counts = tot_cov[genome_counts.first];
            tot_counts.tot_aligned += genome_counts.second.tot_aligned;
            tot_counts.tot_identical += genome_counts.second.tot_identical;
            tot_counts.single_aligned += genome_counts.second.single_aligned;
            tot_counts.single_identical += genome_counts.second.single_identical;
            tot_counts.prev_ref_pos = numeric_limits<int64_t>::max();
            for (const auto& gc : genome_counts.second.gap_hist) {
                tot_counts.gap_hist[gc.first] += gc.second;
            }
        }
    }
}

void add_final_gap(ContigCoverageMap& contig_cov_map) {
    for (auto& contig_covmap : contig_cov_map) {
        for (auto& genome_cov : contig_covmap.second.genome_map) {
            // add in the final gap
            int64_t gap_len = contig_covmap.second.ref_length - genome_cov.second.prev_ref_pos - 1;
            if (gap_len > 0) {
                genome_cov.second.gap_hist[gap_len] += 1;
            }
        }
    }
}

void postprocess_gap_hist(ContigCoverageMap& contig_cov_map) {
    for (auto& contig_covmap : contig_cov_map) {
        for (auto& genome_cov : contig_covmap.second.genome_map) {
            // make gap_hist cumulative
            // before: gap_hist[i] is the number of gaps with length == i
            // after: gap_hist[i] is the number of gap BASES with length >=i 
            int64_t running_total = 0;
            for (auto gci = genome_cov.second.gap_hist.rbegin();
                 gci != genome_cov.second.gap_hist.rend();
                 ++gci) {
                running_total += gci->second * gci->first;
                gci->second = running_total;
            }

            // add in a point for whole contig
            genome_cov.second.gap_hist[numeric_limits<int64_t>::max()] = 0;
        }        
    }
}


void print_coverage_tsv(const ContigCoverageMap& contig_cov_map, const set<int64_t>& gap_thresholds, ostream& os) {
    os << "contig" << "\t"
       << "max-gap" << "\t"
       << "len" << "\t"
       << "query" << "\t"
       << "aln" << "\t"
       << "ident" << "\t"
       << "1:1-aln" << "\t"
       << "1:1-ident" << "\t"        
       << "aln-bp" << "\t"
       << "ident-bp" << "\t"
       << "1:1-aln-bp" << "\t"
       << "1:1-ident-bp" << endl;

    for (const auto& contig_cov : contig_cov_map) {
        set<int64_t> contig_gap_thresholds;
        // replace -1 with the contig length (for prettier output)
        for (int64_t gt : gap_thresholds) {
            if (gt >= 0) {
                contig_gap_thresholds.insert(gt);
            } else {
                assert(gt == -1);
                contig_gap_thresholds.insert(contig_cov.second.ref_length);
            }
        }        
        for (const auto& genome_counts : contig_cov.second.genome_map) {
            for (int64_t max_gap : contig_gap_thresholds) {
                int64_t ref_length = contig_cov.second.ref_length;
                int64_t gap_length = genome_counts.second.gap_hist.upper_bound(max_gap)->second;
                ref_length -= gap_length;
                double tot_aligned_pct = 0;
                double tot_identical_pct = 0;
                double single_aligned_pct = 0;
                double single_identical_pct = 0;
                if (ref_length > 0) {
                    tot_aligned_pct = (double)genome_counts.second.tot_aligned / ref_length;
                    single_aligned_pct = (double)genome_counts.second.single_aligned / ref_length;                
                }
                if (genome_counts.second.tot_aligned > 0) {
                    tot_identical_pct = (double)genome_counts.second.tot_identical / genome_counts.second.tot_aligned;
                }
                if (genome_counts.second.single_aligned > 0) {
                    single_identical_pct = (double)genome_counts.second.single_identical / genome_counts.second.single_aligned;
                }
                os << contig_cov.first << "\t"
                   << max_gap << "\t"
                   << ref_length << "\t"
                   << genome_counts.first << "\t"
                   << std::setprecision(4) << std::fixed
                   << tot_aligned_pct << "\t"
                   << tot_identical_pct << "\t"
                   << single_aligned_pct << "\t"
                   << single_identical_pct << "\t"
                   << genome_counts.second.tot_aligned << "\t"
                   << genome_counts.second.tot_identical << "\t"
                   << genome_counts.second.single_aligned << "\t"
                   << genome_counts.second.single_identical << endl;
            }
        }        
    }
}
