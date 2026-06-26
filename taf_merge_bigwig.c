/*
 * taffy merge-bigwig: merge per-shard 64-bit per-species VECTOR bigWigs into one.
 *
 * Each input is a valid vector bigWig on a single chrom uni0 of size T, carrying
 * data only in its own disjoint column slice; the shards are given in COLUMN
 * order (the SLURM build's zero-padded shard glob is already column-ordered).
 * bigWig is binary, so the shards can't just be cat'd -- this re-streams every
 * shard's intervals into one output bigWig (rebuilding the R-tree index + zoom),
 * and copies the <shard0>.names sidecar to <out>.names.
 *
 *  Released under the MIT license, see LICENSE.txt
 */

#include "taf.h"
#include "sonLib.h"
#include "bigWig.h"
#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void usage(void) {
    fprintf(stderr, "usage: taffy merge-bigwig -o OUT.bw SHARD0.bw SHARD1.bw ...\n\n");
    fprintf(stderr, "  Merge per-shard per-species VECTOR bigWigs (disjoint uni0 column slices of the\n");
    fprintf(stderr, "  same axis, given in COLUMN ORDER) into a single bigWig.  Also writes OUT.bw.names\n");
    fprintf(stderr, "  from SHARD0.bw.names (the .names sidecars are identical across shards).\n\n");
    fprintf(stderr, "  -o FILE   output bigWig (required)\n");
    fprintf(stderr, "  -h        show this help\n");
}

// copy <src_bw>.names -> <dst_bw>.names; returns 0 on success (or if src missing)
static int copy_names(const char *src_bw, const char *dst_bw) {
    char *src = stString_print("%s.names", src_bw);
    char *dst = stString_print("%s.names", dst_bw);
    int rv = 0;
    FILE *in = fopen(src, "rb");
    if (!in) {
        fprintf(stderr, "taffy merge-bigwig: WARNING: %s missing -- no .names sidecar written\n", src);
    } else {
        FILE *out = fopen(dst, "wb");
        if (!out) {
            fprintf(stderr, "taffy merge-bigwig: cannot write %s\n", dst);
            rv = 1;
        } else {
            char buf[65536];
            size_t n;
            while ((n = fread(buf, 1, sizeof buf, in)) > 0) {
                if (fwrite(buf, 1, n, out) != n) { rv = 1; break; }
            }
            if (fclose(out) != 0) rv = 1;
            if (!rv) fprintf(stderr, "taffy merge-bigwig: wrote %s\n", dst);
        }
        fclose(in);
    }
    free(src);
    free(dst);
    return rv;
}

int taf_merge_bigwig_main(int argc, char *argv[]) {
    char *out_path = NULL;
    int c;
    while ((c = getopt(argc, argv, "o:h")) != -1) {
        switch (c) {
            case 'o': out_path = optarg; break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }
    if (out_path == NULL) { fprintf(stderr, "taffy merge-bigwig: -o OUT.bw is required\n"); usage(); return 1; }
    int n_shards = argc - optind;
    if (n_shards < 1) { fprintf(stderr, "taffy merge-bigwig: need at least one shard bigWig\n"); usage(); return 1; }
    char **shards = argv + optind;

    if (bwInit(1 << 17) != 0) { fprintf(stderr, "taffy merge-bigwig: bwInit failed\n"); return 1; }

    // ---- read N + T (+ chrom name) from the first shard's header ----
    bigWigFile_t *s0 = bwOpen(shards[0], NULL, "r");
    if (s0 == NULL) { fprintf(stderr, "taffy merge-bigwig: cannot open %s\n", shards[0]); bwCleanup(); return 1; }
    uint32_t N = s0->vecN;
    if (N < 1) {
        fprintf(stderr, "taffy merge-bigwig: %s is not a vector bigWig (vecN=0) -- expected `depth --perSpecies` output\n", shards[0]);
        bwClose(s0); bwCleanup(); return 1;
    }
    if (s0->cl == NULL || s0->cl->nKeys != 1) {
        fprintf(stderr, "taffy merge-bigwig: %s must be a single-chrom (uni0) bigWig\n", shards[0]);
        bwClose(s0); bwCleanup(); return 1;
    }
    uint64_t T = s0->cl->len[0];
    char *chrom = stString_copy(s0->cl->chrom[0]);   // "uni0"
    bwClose(s0);
    fprintf(stderr, "taffy merge-bigwig: N=%u components, T=%" PRIu64 " columns (chrom=%s), %d shards\n",
            N, T, chrom, n_shards);

    // ---- open the output ----
    bigWigFile_t *out = bwOpen(out_path, NULL, "w");
    if (out == NULL) { fprintf(stderr, "taffy merge-bigwig: cannot open output %s\n", out_path); free(chrom); bwCleanup(); return 1; }
    if (bwCreateHdrVec(out, 10, N)) {
        fprintf(stderr, "taffy merge-bigwig: bwCreateHdrVec failed\n");
        bwClose(out); free(chrom); bwCleanup(); return 1;
    }
    const char *cn[1] = { chrom };
    uint64_t    cl[1] = { T };
    out->cl = bwCreateChromList(cn, cl, 1);
    if (out->cl == NULL || bwWriteHdr(out)) {
        fprintf(stderr, "taffy merge-bigwig: bwWriteHdr failed\n");
        bwClose(out); free(chrom); bwCleanup(); return 1;
    }

    // reusable chrom[] for the (single) first bwAddIntervalsVec batch
    const uint32_t BATCH = 1u << 20;   // bound the chrom[] array + per-call size
    const char **chroms = st_malloc((size_t) BATCH * sizeof(char *));
    for (uint32_t i = 0; i < BATCH; i++) chroms[i] = chrom;

    int wrote = 0;            // first non-empty batch uses bwAddIntervalsVec, rest bwAppendIntervalsVec
    int fail = 0;
    uint64_t total_intervals = 0;
    uint64_t prev_end = 0;    // shards must be in non-overlapping column order

    for (int s = 0; s < n_shards && !fail; s++) {
        bigWigFile_t *sh = bwOpen(shards[s], NULL, "r");
        if (sh == NULL) { fprintf(stderr, "taffy merge-bigwig: cannot open shard %s\n", shards[s]); fail = 1; break; }
        if (sh->vecN != N || sh->cl == NULL || sh->cl->nKeys != 1 || sh->cl->len[0] != T) {
            fprintf(stderr, "taffy merge-bigwig: shard %s inconsistent (vecN=%u T=%" PRIu64 "; expected vecN=%u T=%" PRIu64 ")\n",
                    shards[s], sh->vecN, (sh->cl ? sh->cl->len[0] : (uint64_t) 0), N, T);
            bwClose(sh); fail = 1; break;
        }
        bwOverlappingIntervalsVec_t *iv = bwGetOverlappingIntervalsVec(sh, chrom, 0, T);
        if (iv == NULL) { fprintf(stderr, "taffy merge-bigwig: read failed on shard %s\n", shards[s]); bwClose(sh); fail = 1; break; }
        if (iv->l > 0) {
            if (iv->start[0] < prev_end) {
                fprintf(stderr, "taffy merge-bigwig: shards out of order / overlapping at %s (start %" PRIu64 " < prev end %" PRIu64 ")\n",
                        shards[s], iv->start[0], prev_end);
                bwDestroyOverlappingIntervalsVec(iv); bwClose(sh); fail = 1; break;
            }
            for (uint64_t off = 0; off < iv->l && !fail; off += BATCH) {
                uint32_t n = (uint32_t) ((iv->l - off < (uint64_t) BATCH) ? (iv->l - off) : (uint64_t) BATCH);
                int rv;
                if (!wrote) {
                    rv = bwAddIntervalsVec(out, chroms, iv->start + off, iv->end + off, iv->value + (size_t) off * N, n);
                    if (!rv) wrote = 1;
                } else {
                    rv = bwAppendIntervalsVec(out, iv->start + off, iv->end + off, iv->value + (size_t) off * N, n);
                }
                if (rv) { fprintf(stderr, "taffy merge-bigwig: write failed (shard %s)\n", shards[s]); fail = 1; break; }
            }
            total_intervals += iv->l;
            prev_end = iv->end[iv->l - 1];
        }
        bwDestroyOverlappingIntervalsVec(iv);
        bwClose(sh);
    }

    bwClose(out);   // finalize: builds the index + zoom
    bwCleanup();
    free(chroms);

    if (fail) {
        remove(out_path);   // never leave a partial/corrupt .bw
        free(chrom);
        return 1;
    }
    fprintf(stderr, "taffy merge-bigwig: merged %" PRIu64 " intervals -> %s\n", total_intervals, out_path);
    int names_rv = copy_names(shards[0], out_path);
    free(chrom);
    return names_rv;
}
