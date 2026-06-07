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
    int64_t chainOpen     = -1;     // -1 = leave at handle default
    int64_t chainExtend   = -1;
    int64_t chainMaxGap   = -1;
    double  chainOverlapFrac = -2.0;  // -2 = leave at handle default (0.0)
    int64_t maxOutputBlocks = -1;     // -1 = leave at handle default (500)
};

static void usage(const char *prog) {
    fprintf(stderr,
        "usage: %s [options] <tuiPath> <qSpecies> <tSpecies> <tChrom> <tStart> <tEnd>\n"
        "options:\n"
        "  --doSeq                request sequence (no-op in initial cut)\n"
        "  --doDupes              include target-dupe list\n"
        "  --verbose              verbose tracing on stderr\n"
        "  --mapBack              include mapBackAdjacencies off-screen flanks\n"
        "  --filterByChrom <chr>  restrict output to a single qChrom\n"
        "  --dupMode MODE         noDups | queryDups | all  (default: all)\n"
        "  --maxGap N             override per-handle chain max_gap_length\n"
        "  --chainOpen N          override per-handle chain_open\n"
        "  --chainExtend N        override per-handle chain_extend\n"
        "  --chainOverlapFrac F   override per-handle chain_overlap_frac (-1 = filter off, [0,1] = active)\n"
        "  --maxOutputBlocks N    override per-handle max_output_blocks cap (default 500)\n"
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
        else if (opt == "--maxGap" || opt == "--chainMaxGap") {
            if (i + 1 >= argc) { fprintf(stderr, "%s needs an integer\n", opt.c_str()); return -1; }
            a->chainMaxGap = (int64_t) atoll(argv[i + 1]); i += 2;
        }
        else if (opt == "--chainOpen") {
            if (i + 1 >= argc) { fprintf(stderr, "--chainOpen needs an integer\n"); return -1; }
            a->chainOpen = (int64_t) atoll(argv[i + 1]); i += 2;
        }
        else if (opt == "--chainExtend") {
            if (i + 1 >= argc) { fprintf(stderr, "--chainExtend needs an integer\n"); return -1; }
            a->chainExtend = (int64_t) atoll(argv[i + 1]); i += 2;
        }
        else if (opt == "--chainOverlapFrac") {
            if (i + 1 >= argc) { fprintf(stderr, "--chainOverlapFrac needs a float\n"); return -1; }
            a->chainOverlapFrac = atof(argv[i + 1]); i += 2;
        }
        else if (opt == "--maxOutputBlocks") {
            if (i + 1 >= argc) { fprintf(stderr, "--maxOutputBlocks needs an integer\n"); return -1; }
            a->maxOutputBlocks = (int64_t) atoll(argv[i + 1]); i += 2;
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

    if (args.chainOpen >= 0 || args.chainExtend >= 0 || args.chainMaxGap >= 0) {
        if (taffySetChainParams(h, args.chainOpen, args.chainExtend,
                                args.chainMaxGap, &errStr) != 0) {
            die("taffySetChainParams", errStr);
        }
        if (args.verbose) {
            int64_t co = 0, ce = 0, mg = 0;
            taffyGetChainParams(h, &co, &ce, &mg, nullptr);
            fprintf(stderr, ">> chain params: open=%ld extend=%ld max_gap=%ld\n",
                    (long) co, (long) ce, (long) mg);
        }
    }
    if (args.chainOverlapFrac > -1.5) {
        if (taffySetChainOverlapFrac(h, args.chainOverlapFrac, &errStr) != 0) {
            die("taffySetChainOverlapFrac", errStr);
        }
    }
    if (args.maxOutputBlocks > 0) {
        if (taffySetMaxOutputBlocks(h, args.maxOutputBlocks, &errStr) != 0) {
            die("taffySetMaxOutputBlocks", errStr);
        }
    }
    if (args.verbose) {
        double f = 0;
        taffyGetChainOverlapFrac(h, &f, nullptr);
        fprintf(stderr, ">> chain_overlap_frac = %.3f\n", f);
    }

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
        printf("%s\t%ld\t%ld\t%ld\t%c\t%ld\n",
               b->qChrom, (long) b->tStart, (long) b->qStart, (long) b->size,
               b->strand, (long) b->chainId);
        n++;
    }
    if (args.verbose) fprintf(stderr, ">> emitted %ld mapped blocks "
                                       "(last column = chainId; 0 = bin or flank)\n",
                                       (long) n);

    if (args.verbose && res->chainSummaries) {
        fprintf(stderr, ">> chainSummaries (score-desc; primary first):\n");
        int64_t cn = 0;
        for (struct taffy_chain_summary_t *c = res->chainSummaries; c; c = c->next) {
            fprintf(stderr, "      id=%ld score=%ld bp=%ld n_alns=%ld\n",
                    (long) c->id, (long) c->totalScore,
                    (long) c->totalBp, (long) c->nAlns);
            if (++cn >= 10) { fprintf(stderr, "      ...\n"); break; }
        }
    }

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
