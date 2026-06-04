/*
 * taffyBlockVizTest -- CLI mirror of hal's blockVizTest.cpp.
 *
 * Output format matches blockVizTest as closely as the data shape
 * allows so the two tools can be diffed on the same alignment.
 *
 * Args (positional):
 *   <tuiPath>  <qSpecies>  <tSpecies>  <tChrom>  <tStart>  <tEnd>
 *
 * Options:
 *   --doSeq             request sequence (currently ignored -- our
 *                       blockViz returns no sequence in the initial cut)
 *   --doDupes           include target-dupe list (currently a no-op
 *                       since the initial cut routes paralogs into
 *                       mappedBlocks)
 *   --verbose           verbose tracing
 *   --filterByChrom S   pass S as qChrom to filter the output
 */

#include "taffyBlockViz.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

struct Args {
    std::string path;
    std::string qSpecies;
    std::string tSpecies;
    std::string tChrom;
    taffy_int_t tStart = 0;
    taffy_int_t tEnd   = 0;
    int doSeq = 0;
    int doDupes = 0;
    int verbose = 0;
    int mapBack = 0;
    std::string qChrom;       // "" = unfiltered
    taffy_dup_type_t dupMode = TAFFY_QUERY_AND_TARGET_DUPS;
};

static void usage(const char *prog) {
    fprintf(stderr,
        "usage: %s [options] <tuiPath> <qSpecies> <tSpecies> <tChrom> <tStart> <tEnd>\n"
        "options:\n"
        "  --doSeq                request sequence (no-op in initial cut)\n"
        "  --doDupes              include target-dupe list (no-op in initial cut)\n"
        "  --verbose              verbose tracing on stderr\n"
        "  --filterByChrom <chr>  restrict output to a single qChrom\n"
        "  -h / --help            this help\n", prog);
}

static int parse_args(int argc, char **argv, Args *a) {
    int i = 1;
    while (i < argc && argv[i][0] == '-') {
        std::string opt = argv[i];
        if (opt == "--doSeq")        { a->doSeq = 1; i++; }
        else if (opt == "--doDupes") { a->doDupes = 1; i++; }
        else if (opt == "--verbose") { a->verbose = 1; i++; }
        else if (opt == "--mapBack") { a->mapBack = 1; i++; }
        else if (opt == "--filterByChrom") {
            if (i + 1 >= argc) { fprintf(stderr, "--filterByChrom needs an argument\n"); return -1; }
            a->qChrom = argv[i + 1]; i += 2;
        }
        else if (opt == "--dupMode") {
            if (i + 1 >= argc) { fprintf(stderr, "--dupMode needs an argument (noDups|queryDups|all)\n"); return -1; }
            std::string m = argv[i + 1];
            if (m == "noDups")      a->dupMode = TAFFY_NO_DUPS;
            else if (m == "queryDups") a->dupMode = TAFFY_QUERY_DUPS;
            else if (m == "all")    a->dupMode = TAFFY_QUERY_AND_TARGET_DUPS;
            else { fprintf(stderr, "--dupMode must be noDups|queryDups|all\n"); return -1; }
            i += 2;
        }
        else if (opt == "-h" || opt == "--help") { usage(argv[0]); return 1; }
        else { fprintf(stderr, "unknown option: %s\n", opt.c_str()); return -1; }
    }
    if (argc - i != 6) { usage(argv[0]); return -1; }
    a->path     = argv[i++];
    a->qSpecies = argv[i++];
    a->tSpecies = argv[i++];
    a->tChrom   = argv[i++];
    a->tStart   = (taffy_int_t) atoll(argv[i++]);
    a->tEnd     = (taffy_int_t) atoll(argv[i++]);
    return 0;
}

static void die(const char *what, char *errStr) {
    fprintf(stderr, "ERROR: %s: %s\n", what, errStr ? errStr : "(no message)");
    free(errStr);
    exit(1);
}

int main(int argc, char **argv) {
    Args args;
    int rc = parse_args(argc, argv, &args);
    if (rc != 0) return rc < 0 ? 1 : 0;

    char *errStr = nullptr;
    int h = taffyOpen(args.path.c_str(), &errStr);
    if (h < 0) die("taffyOpen", errStr);

    if (args.verbose) {
        fprintf(stderr, ">> opened %s (handle=%d)\n", args.path.c_str(), h);

        // List species (small alignment) -- HAL's blockVizTest does
        // this too in --verbose.
        struct taffy_species_t *sp = taffyGetSpecies(h, &errStr);
        if (!sp) die("taffyGetSpecies", errStr);
        fprintf(stderr, ">> species in alignment:\n");
        for (struct taffy_species_t *s = sp; s; s = s->next) {
            fprintf(stderr, "      %s\tlen=%ld\tnumChroms=%ld\n",
                    s->name, (long) s->length, (long) s->numChroms);
        }
        taffyFreeSpeciesList(sp);

        // List chroms for tSpecies.
        struct taffy_chromosome_t *cl = taffyGetChroms(h, args.tSpecies.c_str(), &errStr);
        if (!cl) die("taffyGetChroms", errStr);
        fprintf(stderr, ">> %s chroms:\n", args.tSpecies.c_str());
        int n_show = 0;
        for (struct taffy_chromosome_t *c = cl; c; c = c->next) {
            fprintf(stderr, "      %s\t%ld\n", c->name, (long) c->length);
            if (++n_show >= 5) { fprintf(stderr, "      ...\n"); break; }
        }
        taffyFreeChromList(cl);

        fprintf(stderr, ">> querying blocks for %s:%ld-%ld in %s (qSpecies=%s%s%s)\n",
                args.tChrom.c_str(), (long) args.tStart, (long) args.tEnd,
                args.tSpecies.c_str(), args.qSpecies.c_str(),
                args.qChrom.empty() ? "" : ", qChrom=",
                args.qChrom.empty() ? "" : args.qChrom.c_str());
    }

    // The hot call.
    struct taffy_block_results_t *res = nullptr;
    if (args.qChrom.empty()) {
        res = taffyGetBlocksInTargetRange(
            h, args.qSpecies.c_str(), args.tSpecies.c_str(),
            args.tChrom.c_str(), args.tStart, args.tEnd,
            /*tReversed=*/0,
            /*seqMode=*/TAFFY_NO_SEQUENCES,
            /*dupMode=*/args.dupMode,
            /*mapBackAdjacencies=*/args.mapBack,
            /*coalescenceLimitName=*/nullptr,
            &errStr);
    } else {
        res = taffyGetBlocksInTargetRange_filterByChrom(
            h, args.qSpecies.c_str(), args.tSpecies.c_str(),
            args.tChrom.c_str(), args.tStart, args.tEnd,
            /*tReversed=*/0,
            /*seqMode=*/TAFFY_NO_SEQUENCES,
            /*dupMode=*/args.dupMode,
            /*mapBackAdjacencies=*/args.mapBack,
            args.qChrom.c_str(),
            /*coalescenceLimitName=*/nullptr,
            &errStr);
    }
    if (!res) die("taffyGetBlocksInTargetRange", errStr);

    // Output: one line per block, tab-separated:
    //   qChrom\ttStart\tqStart\tsize\tstrand
    // (HAL's blockVizTest prints similar fields; the strand char +
    // numeric columns line up for direct diff after sorting.)
    int64_t n = 0;
    for (struct taffy_block_t *b = res->mappedBlocks; b; b = b->next) {
        printf("%s\t%ld\t%ld\t%ld\t%c\n",
               b->qChrom, (long) b->tStart, (long) b->qStart, (long) b->size, b->strand);
        n++;
    }
    if (args.verbose) fprintf(stderr, ">> emitted %ld mapped blocks\n", (long) n);

    if (args.doDupes) {
        int64_t nd = 0, nr = 0;
        printf("--- targetDupeBlocks ---\n");
        for (struct taffy_target_dupe_list_t *d = res->targetDupeBlocks; d; d = d->next) {
            nd++;
            for (struct taffy_target_range_t *r = d->tRange; r; r = r->next) {
                printf("dupe id=%ld\tqChrom=%s\ttStart=%ld\tsize=%ld\n",
                       (long) d->id, d->qChrom, (long) r->tStart, (long) r->size);
                nr++;
            }
        }
        if (args.verbose) fprintf(stderr, ">> emitted %ld dupe entries (%ld ranges)\n", (long) nd, (long) nr);
    }

    taffyFreeBlockResults(res);
    if (taffyClose(h, &errStr) != 0) die("taffyClose", errStr);
    return 0;
}
