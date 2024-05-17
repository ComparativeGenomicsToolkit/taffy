/*
 * taffy: The TAFFY command line interface, whose subcommands run different tools
 *
 *  Released under the MIT license, see LICENSE.txt
*/
extern "C" {
#include "taf.h"
#include "sonLib.h"
extern int taf_norm_main(int argc, char *argv[]);
extern int taf_index_main(int argc, char *argv[]);
extern int taf_view_main(int argc, char *argv[]);
extern int taf_sort_main(int argc, char *argv[]);
extern int taf_stats_main(int argc, char *argv[]);
extern int taf_annotate_main(int argc, char *argv[]);
}

extern int taf_add_gap_bases_main(int argc, char *argv[]);
extern int taf_coverage_main(int argc, char *argv[]);

void usage() {
    fprintf(stderr, "taffy: toolkit for working with TAF and MAF multiple alignment files\n\n");
    fprintf(stderr, "usage: taffy <command> [options]\n\n");
    fprintf(stderr, "available commands:\n");
    fprintf(stderr, "    view           MAF / TAF conversion and region extraction\n");
    fprintf(stderr, "    norm           normalize TAF blocks\n");
    fprintf(stderr, "    add-gap-bases  add sequences from HAL or FASTA files into TAF gaps\n");
    fprintf(stderr, "    index          create a .tai index (required for region extraction)\n");
    fprintf(stderr, "    sort           sort the rows of a TAF file to a desired order\n");
    fprintf(stderr, "    stats          print statistics of a TAF file\n");
    fprintf(stderr, "    coverage       print coverage statistics of a given genome in a TAF file\n");
    fprintf(stderr, "    annotate       annotate a TAF file with labels from a wiggle file\n");
    fprintf(stderr, "\n");

#ifdef USE_HTSLIB
    fprintf(stderr, "all commands accept uncompressed or bgzipped TAF input\n");
#else
    fprintf(stderr, "taffy was compiled without bgzip support: only uncompressed inputs accepted\n");
#endif
    fprintf(stderr, "\nrun taffy <command> -h to show the given command's interface\n\n");
}

int main(int argc, char *argv[]) {

    if (argc < 2) {
        usage();
        return 0;
    }
    if (strcmp(argv[1], "view") == 0) {
        return taf_view_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "norm") == 0) {
        return taf_norm_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "add-gap-bases") == 0) {
        return taf_add_gap_bases_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "index") == 0) {
        return taf_index_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "sort") == 0) {
        return taf_sort_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "stats") == 0) {
        return taf_stats_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "coverage") == 0) {
        return taf_coverage_main(argc - 1, argv + 1);        
    } else if (strcmp(argv[1], "annotate") == 0) {
        return taf_annotate_main(argc - 1, argv + 1);
    } else {
        fprintf(stderr, "%s is not a valid taffy command\n", argv[1]);
        usage();
        return 1;
    }
    return 1;
}
