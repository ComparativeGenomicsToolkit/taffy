/*
 * taf_index: Make a .tai index from a TAF or MAF file (which can be bgzipped or uncompressed)
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "tai.h"
#include "tui.h"
#include "sonLib.h"
#include <getopt.h>
#include <time.h>

// Load a whitespace-separated (space/tab/newline) list of genome names into an
// stHash (keys = names, dummy non-NULL values) -- the shape extract_genome_name
// wants for its prefix-walk.  This is a plain membership list: no renaming.
static stHash *load_genome_name_list(const char *path) {
    FILE *fh = fopen(path, "r");
    if (fh == NULL) {
        fprintf(stderr, "Unable to open genome name list: %s\n", path);
        exit(1);
    }
    stHash *names = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    char tok[8192];
    int c, n = 0;
    do {
        c = fgetc(fh);
        if (c == EOF || c == ' ' || c == '\t' || c == '\n' ||
            c == '\r' || c == '\f' || c == '\v') {
            if (n > 0) {
                tok[n] = '\0';
                if (stHash_search(names, tok) == NULL) {
                    stHash_insert(names, stString_copy(tok), (void *)1);
                }
                n = 0;
            }
        } else if (n < (int)sizeof(tok) - 1) {
            tok[n++] = (char)c;
        }
    } while (c != EOF);
    fclose(fh);
    return names;
}

static void usage(void) {
    fprintf(stderr, "taffy index [options]\n");
    fprintf(stderr, "Index a TAF or MAF file, output goes in <file>.tai\n");
    fprintf(stderr, "-i --inputFile : Input taf or maf file [REQUIRED]\n");
    fprintf(stderr, "-b --blockSize : Write an index line for intervals of this many bp [default:10000]\n");
    fprintf(stderr, "-u --universal : Build a universal column index in <file>.tui (ONEcode) "
                     "instead of the .tai. For a `cactus-hal2maf --universal` MAF/TAF.\n");
    fprintf(stderr, "-d --tmpDir : Directory for --universal spill files "
                     "[default: directory of the output .tui]\n");
    fprintf(stderr, "-n --genomeNames : File listing genome names "
                     "(whitespace-separated), ex from `halStats --genomes`.  "
                     "Used only to split genome.sequence names; no renaming.  "
                     "For --universal when a name has >1 '.': if omitted, the "
                     "genome set is taken from the `# hal` tree in the header "
                     "(hal2maf writes one); this option overrides it.\n");
    fprintf(stderr, "-T --threads N : Use N threads for bgzf I/O (default 1, only effective on bgzipped streams)\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int taf_index_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *taf_fn = NULL;
    int64_t block_size = 10000;
    bool universal = false;
    char *tmp_dir = NULL;
    char *genome_names_file = NULL;
    int bgzf_threads = 1;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "blockSize", required_argument, 0, 'b' },
                                                { "universal", no_argument, 0, 'u' },
                                                { "tmpDir", required_argument, 0, 'd' },
                                                { "genomeNames", required_argument, 0, 'n' },
                                                { "threads", required_argument, 0, 'T' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:b:ud:n:hT:", long_options, &option_index);
        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'i':
                taf_fn = optarg;
                break;
            case 'b':
                block_size = atoi(optarg);
                break;
            case 'u':
                universal = true;
                break;
            case 'd':
                tmp_dir = optarg;
                break;
            case 'n':
                genome_names_file = optarg;
                break;
            case 'T':
                bgzf_threads = atoi(optarg);
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    LI_set_bgzf_threads(bgzf_threads);
    st_logInfo("Input file string : %s\n", taf_fn);
    st_logInfo("Block size : %" PRIi64 "\n", block_size);
    
    //////////////////////////////////////////////
    // Make the .tai index
    //////////////////////////////////////////////

    if (taf_fn == NULL) {
        fprintf(stderr, "Input file must be specified with -i\n");
        return 1;        
    }
    FILE *taf_fh = fopen(taf_fn, "r");
    if (taf_fh == NULL) {
        fprintf(stderr, "Unable to open input file: %s\n", taf_fn);
        return 1;
    }
    LI *li = LI_construct(taf_fh);
    if (!LI_indexable(li)) {
        fprintf(stderr, "Input file must be either uncompressed or bgzipped: gzip not supported: %s\n", taf_fn);
        return 1;
    }

    char *tai_fn = NULL;
    FILE *tai_fh = NULL;
    int rv = 0;
    if (universal) {
        stHash *genome_name_map = NULL;
        if (genome_names_file != NULL) {
            st_logInfo("Genome names file string : %s\n", genome_names_file);
            genome_name_map = load_genome_name_list(genome_names_file);
        }
        char *tui_fn = tui_path(taf_fn);
        st_logInfo("Output index file : %s\n", tui_fn);
        rv = tui_create(li, tui_fn, tmp_dir, genome_name_map, bgzf_threads);
        free(tui_fn);
        if (genome_name_map != NULL) {
            stHash_destruct(genome_name_map);
        }
    } else {
        tai_fn = tai_path(taf_fn);
        st_logInfo("Output index file : %s\n", tai_fn);
        tai_fh = fopen(tai_fn, "w");
        tai_create(li, tai_fh, block_size);
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if(taf_fh != NULL) {
        fclose(taf_fh);
    }
    if(tai_fh != NULL) {
        fclose(tai_fh);
    }
    free(tai_fn);

    LI_destruct(li);
    
    st_logInfo("taffy index is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return rv;
}


