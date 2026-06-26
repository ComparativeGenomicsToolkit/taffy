#include "bigWig.h"
#include "bwCommon.h"
#include <errno.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <string.h>

//Returns -1 if there are no applicable levels, otherwise an integer indicating the most appropriate level.
//Like Kent's library, this divides the desired bin size by 2 to minimize the effect of blocks overlapping multiple bins
static int32_t determineZoomLevel(const bigWigFile_t *fp, int basesPerBin) {
    int32_t out = -1;
    int64_t diff;
    uint32_t bestDiff = -1;
    uint16_t i;

    basesPerBin/=2;
    for(i=0; i<fp->hdr->nLevels; i++) {
        diff = basesPerBin - (int64_t) fp->hdr->zoomHdrs->level[i];
        if(diff >= 0 && diff < bestDiff) {
            bestDiff = diff;
            out = i;
        }
    }
    return out;
}

/// @cond SKIP
struct val_t {
    uint32_t nBases;
    float min, max, sum, sumsq;
    double scalar;
};

struct vals_t {
    uint32_t n;
    struct val_t **vals;
};
/// @endcond

void destroyVals_t(struct vals_t *v) {
    uint32_t i;
    if(!v) return;
    for(i=0; i<v->n; i++) free(v->vals[i]);
    if(v->vals) free(v->vals);
    free(v);
}

//Determine the base-pair overlap between an interval and a block
double getScalar(uint64_t i_start, uint64_t i_end, uint64_t b_start, uint64_t b_end) {  //64-bit fork: positions u64
    double rv = 0.0;
    if(b_start <= i_start) {
        //64-bit fork (H4): cap the numerator at i_end when the block fully spans the query
        if(b_end > i_start) rv = ((double)((b_end < i_end ? b_end : i_end) - i_start))/(b_end-b_start);
    } else if(b_start < i_end) {
        if(b_end < i_end) rv = ((double)(b_end - b_start))/(b_end-b_start);
        else rv = ((double)(i_end - b_start))/(b_end-b_start);
    }

    return rv;
}

//Returns NULL on error
static struct vals_t *getVals(bigWigFile_t *fp, bwOverlapBlock_t *o, int i, uint32_t tid, uint64_t start, uint64_t end) {
    void *buf = NULL, *compBuf = NULL;
    uLongf sz = fp->hdr->bufSize;
    int compressed = 0, rv;
    uint8_t *p; uint32_t vtid; uint64_t vstart, vend;  //64-bit fork: positions u64, byte-ptr for 40B records
    struct vals_t *vals = NULL;
    struct val_t *v = NULL;

    if(sz) {
        compressed = 1;
        buf = malloc(sz);
    }
    sz = 0; //This is now the size of the compressed buffer

    if(bwSetPos(fp, o->offset[i])) goto error;

    vals = calloc(1,sizeof(struct vals_t));
    if(!vals) goto error;

    v = malloc(sizeof(struct val_t));
    if(!v) goto error;

    if(sz < o->size[i]) compBuf = malloc(o->size[i]);
    if(!compBuf) goto error;

    if(bwRead(compBuf, o->size[i], 1, fp) != 1) goto error;
    if(compressed) {
        sz = fp->hdr->bufSize;
        rv = uncompress(buf, &sz, compBuf, o->size[i]);
        if(rv != Z_OK) goto error;
    } else {
        buf = compBuf;
        sz = o->size[i];
    }

    p = (uint8_t*)buf;
    while(((uLongf)(p - (uint8_t*)buf)) < sz) {
        //64-bit fork: 40B zoom record -- chromId u32@0, start u64@4, end u64@12, nBases u32@20, min/max/sum/sumsq @24/28/32/36
        vtid = *(uint32_t*)(p+0);
        vstart = *(uint64_t*)(p+4);
        vend = *(uint64_t*)(p+12);
        v->nBases = *(uint32_t*)(p+20);
        v->min = *(float*)(p+24);
        v->max = *(float*)(p+28);
        v->sum = *(float*)(p+32);
        v->sumsq = *(float*)(p+36);
        v->scalar = getScalar(start, end, vstart, vend);

        if(tid == vtid) {
            if((start <= vstart && end > vstart) || (start < vend && start >= vstart)) {
                vals->vals = realloc(vals->vals, sizeof(struct val_t*)*(vals->n+1));
                if(!vals->vals) goto error;
                vals->vals[vals->n++] = v;
                v = malloc(sizeof(struct val_t));
                if(!v) goto error;
            }
            if(vstart > end) break;
        } else if(vtid > tid) {
            break;
        }
        p += 40;  //64-bit fork: advance one 40B zoom record
    }

    free(v);
    free(buf);
    if(compressed) free(compBuf);
    return vals;

error:
    if(buf) free(buf);
    if(compBuf && compressed) free(compBuf);
    if(v) free(v);
    destroyVals_t(vals);
    return NULL;
}

//On error, errno is set to ENOMEM and NaN is returned (though NaN can be returned normally)
static double blockMean(bigWigFile_t *fp, bwOverlapBlock_t *blocks, uint32_t tid, uint64_t start, uint64_t end) {
    uint32_t i, j;
    double output = 0.0, coverage = 0.0;
    struct vals_t *v = NULL;

    if(!blocks->n) return strtod("NaN", NULL);

    //Iterate over the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            output += v->vals[j]->sum * v->vals[j]->scalar;
            coverage += v->vals[j]->nBases * v->vals[j]->scalar;
        }
        destroyVals_t(v);
    }


    if(!coverage) return strtod("NaN", NULL);

    return output/coverage;

error:
    if(v) free(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

static double intMean(bwOverlappingIntervals_t* ints, uint64_t start, uint64_t end) {
    double sum = 0.0;
    uint32_t i; uint64_t nBases = 0, start_use, end_use;  //64-bit fork: positions u64

    if(!ints->l) return strtod("NaN", NULL);

    for(i=0; i<ints->l; i++) {
        start_use = ints->start[i];
        end_use = ints->end[i];
        if(ints->start[i] < start) start_use = start;
        if(ints->end[i] > end) end_use = end;
        nBases += end_use-start_use;
        sum += (end_use-start_use)*((double) ints->value[i]);
    }

    return sum/nBases;
}

//Does UCSC compensate for partial block/range overlap?
static double blockDev(bigWigFile_t *fp, bwOverlapBlock_t *blocks, uint32_t tid, uint64_t start, uint64_t end) {
    uint32_t i, j;
    double mean = 0.0, ssq = 0.0, coverage = 0.0, diff;
    struct vals_t *v = NULL;

    if(!blocks->n) return strtod("NaN", NULL);

    //Iterate over the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            coverage += v->vals[j]->nBases * v->vals[j]->scalar;
            mean += v->vals[j]->sum * v->vals[j]->scalar;
            ssq += v->vals[j]->sumsq * v->vals[j]->scalar;
        }
        destroyVals_t(v);
        v = NULL;
    }

    if(coverage<=1.0) return strtod("NaN", NULL);
    diff = ssq-mean*mean/coverage;
    if(coverage > 1.0) diff /= coverage-1;
    if(fabs(diff) > 1e-8) { //Ignore floating point differences
        return sqrt(diff);
    } else {
        return 0.0;
    }

error:
    if(v) destroyVals_t(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

//This uses compensated summation to account for finite precision math
static double intDev(bwOverlappingIntervals_t* ints, uint64_t start, uint64_t end) {
    double v1 = 0.0, mean, rv;
    uint32_t i; uint64_t nBases = 0, start_use, end_use;  //64-bit fork: positions u64

    if(!ints->l) return strtod("NaN", NULL);
    mean = intMean(ints, start, end);

    for(i=0; i<ints->l; i++) {
        start_use = ints->start[i];
        end_use = ints->end[i];
        if(ints->start[i] < start) start_use = start;
        if(ints->end[i] > end) end_use = end;
        nBases += end_use-start_use;
        v1 += (end_use-start_use) * pow(ints->value[i]-mean, 2.0); //running sum of squared difference
    }

    if(nBases>=2) rv = sqrt(v1/(nBases-1));
    else if(nBases==1) rv = sqrt(v1);
    else rv = strtod("NaN", NULL);

    return rv;
}

static double blockMax(bigWigFile_t *fp, bwOverlapBlock_t *blocks, uint32_t tid, uint64_t start, uint64_t end) {
    uint32_t i, j, isNA = 1;
    double o = strtod("NaN", NULL);
    struct vals_t *v = NULL;

    if(!blocks->n) return o;

    //Iterate the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            if(isNA) {
                o = v->vals[j]->max;
                isNA = 0;
            } else if(v->vals[j]->max > o) {
                o = v->vals[j]->max;
            }
        }
        destroyVals_t(v);
    }

    return o;

error:
    destroyVals_t(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

static double intMax(bwOverlappingIntervals_t* ints) {
    uint32_t i;
    double o;

    if(ints->l < 1) return strtod("NaN", NULL);

    o = ints->value[0];
    for(i=1; i<ints->l; i++) {
        if(ints->value[i] > o) o = ints->value[i];
    }

    return o;
}

static double blockMin(bigWigFile_t *fp, bwOverlapBlock_t *blocks, uint32_t tid, uint64_t start, uint64_t end) {
    uint32_t i, j, isNA = 1;
    double o = strtod("NaN", NULL);
    struct vals_t *v = NULL;

    if(!blocks->n) return o;

    //Iterate the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            if(isNA) {
                o = v->vals[j]->min;
                isNA = 0;
            } else if(v->vals[j]->min < o) o = v->vals[j]->min;
        }
        destroyVals_t(v);
    }

    return o;

error:
    destroyVals_t(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

static double intMin(bwOverlappingIntervals_t* ints) {
    uint32_t i;
    double o;

    if(ints->l < 1) return strtod("NaN", NULL);

    o = ints->value[0];
    for(i=1; i<ints->l; i++) {
        if(ints->value[i] < o) o = ints->value[i];
    }

    return o;
}

//Does UCSC compensate for only partial block/interval overlap?
static double blockCoverage(bigWigFile_t *fp, bwOverlapBlock_t *blocks, uint32_t tid, uint64_t start, uint64_t end) {
    uint32_t i, j;
    double o = 0.0;
    struct vals_t *v = NULL;

    if(!blocks->n) return strtod("NaN", NULL);

    //Iterate over the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            o+= v->vals[j]->nBases * v->vals[j]->scalar;
        }
        destroyVals_t(v);
    }

    if(o == 0.0) return strtod("NaN", NULL);
    return o;

error:
    destroyVals_t(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

static double intCoverage(bwOverlappingIntervals_t* ints, uint64_t start, uint64_t end) {
    uint32_t i; uint64_t start_use, end_use;  //64-bit fork: positions u64
    double o = 0.0;

    if(!ints->l) return strtod("NaN", NULL);

    for(i=0; i<ints->l; i++) {
        start_use = ints->start[i];
        end_use = ints->end[i];
        if(start_use < start) start_use = start;
        if(end_use > end) end_use = end;
        o += end_use - start_use;
    }

    return o/(end-start);
}

static double blockSum(bigWigFile_t *fp, bwOverlapBlock_t *blocks, uint32_t tid, uint64_t start, uint64_t end) {
    uint32_t i, j;
    double o = 0.0, coverage = 0.0;  //64-bit fork (H3): track coverage to tell no-data from sum==0
    struct vals_t *v = NULL;

    if(!blocks->n) return strtod("NaN", NULL);

    //Iterate over the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            //64-bit fork: getScalar is fractional (commit ad0bbad), so the record's stored
            //sum contributes scaled by its fractional overlap with the query bin. (The old
            //sizeUse=scalar truncated the fraction to a uint32_t 0/1 -> wrong sum at zoom levels.)
            o += v->vals[j]->sum * v->vals[j]->scalar;
            coverage += v->vals[j]->nBases * v->vals[j]->scalar;  //64-bit fork (H3)
        }
        destroyVals_t(v);
    }

    //64-bit fork (H3): NaN only when there is genuinely no data; a covered bin that
    //sums to exactly 0 (e.g. an uncovered depth column) must return 0.0, like intSum.
    if(!coverage) return strtod("NaN", NULL);
    return o;

error:
    destroyVals_t(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

static double intSum(bwOverlappingIntervals_t* ints, uint64_t start, uint64_t end) {
    uint32_t i; uint64_t start_use, end_use;  //64-bit fork: positions u64
    double o = 0.0;

    if(!ints->l) return strtod("NaN", NULL);

    for(i=0; i<ints->l; i++) {
        start_use = ints->start[i];
        end_use = ints->end[i];
        if(start_use < start) start_use = start;
        if(end_use > end) end_use = end;
        o += (end_use - start_use) * ints->value[i];
    }

    return o;
}

//Returns NULL on error, otherwise a double* that needs to be free()d
static double *bwStatsFromZoom(bigWigFile_t *fp, int32_t level, uint32_t tid, uint64_t start, uint64_t end, uint32_t nBins, enum bwStatsType type) {
    bwOverlapBlock_t *blocks = NULL;
    double *output = NULL;
    uint64_t pos = start, end2; uint32_t i;  //64-bit fork: pos/end2 are positions

    if(!fp->hdr->zoomHdrs->idx[level]) {
        fp->hdr->zoomHdrs->idx[level] = bwReadIndex(fp, fp->hdr->zoomHdrs->indexOffset[level]);
        if(!fp->hdr->zoomHdrs->idx[level]) return NULL;
    }
    errno = 0; //Sometimes libCurls sets and then doesn't unset errno on errors

    output = malloc(sizeof(double)*nBins);
    if(!output) return NULL;

    for(i=0, pos=start; i<nBins; i++) {
        end2 = start + ((double)(end-start)*(i+1))/((int) nBins);
        blocks = walkRTreeNodes(fp, fp->hdr->zoomHdrs->idx[level]->root, tid, pos, end2);
        if(!blocks) goto error;

        switch(type) {
        case 0:
            //mean
            output[i] = blockMean(fp, blocks, tid, pos, end2);
            break;
        case 1:
            //stdev
            output[i] = blockDev(fp, blocks, tid, pos, end2);
            break;
        case 2:
            //max
            output[i] = blockMax(fp, blocks, tid, pos, end2);
            break;
        case 3:
            //min
            output[i] = blockMin(fp, blocks, tid, pos, end2);
            break;
        case 4:
            //cov
            output[i] = blockCoverage(fp, blocks, tid, pos, end2)/(end2-pos);
            break;
        case 5:
            //sum
            output[i] = blockSum(fp, blocks, tid, pos, end2);
            break;
        default:
            goto error;
            break;
        }
        if(errno) goto error;
        destroyBWOverlapBlock(blocks);
        pos = end2;
    }

    return output;

error:
    fprintf(stderr, "got an error in bwStatsFromZoom in the range %"PRIu64"-%"PRIu64": %s\n", pos, end2, strerror(errno));
    if(blocks) destroyBWOverlapBlock(blocks);
    if(output) free(output);
    return NULL;
}

double *bwStatsFromFull(bigWigFile_t *fp, const char *chrom, uint64_t start, uint64_t end, uint32_t nBins, enum bwStatsType type) {
    bwOverlappingIntervals_t *ints = NULL;
    double *output = malloc(sizeof(double)*nBins);
    uint32_t i;
    uint64_t pos = start, end2;  //64-bit fork: pos/end2 are positions, must not truncate
    if(!output) return NULL;

    for(i=0; i<nBins; i++) {
        end2 = start + ((double)(end-start)*(i+1))/((int) nBins);
        ints = bwGetOverlappingIntervals(fp, chrom, pos, end2);

        if(!ints) {
            output[i] = strtod("NaN", NULL);
            continue;
        }

        switch(type) {
        default :
        case 0:
            output[i] = intMean(ints, pos, end2);
            break;
        case 1:
            output[i] = intDev(ints, pos, end2);
            break;
        case 2:
            output[i] = intMax(ints);
            break;
        case 3:
            output[i] = intMin(ints);
            break;
        case 4:
            output[i] = intCoverage(ints, pos, end2);
            break;
        case 5:
            output[i] = intSum(ints, pos, end2);
            break;
        }
        bwDestroyOverlappingIntervals(ints);
        pos = end2;
    }

    return output;
}

//Returns a list of floats of length nBins that must be free()d
//On error, NULL is returned
double *bwStats(bigWigFile_t *fp, const char *chrom, uint64_t start, uint64_t end, uint32_t nBins, enum bwStatsType type) {
    int32_t level = determineZoomLevel(fp, ((double)(end-start))/((int) nBins));
    uint32_t tid = bwGetTid(fp, chrom);
    if(tid == (uint32_t) -1) return NULL;

    if(level == -1) return bwStatsFromFull(fp, chrom, start, end, nBins, type);
    return bwStatsFromZoom(fp, level, tid, start, end, nBins, type);
}

//================ 64-bit vector fork: per-component stats (VECTOR_FORMAT.md) ================
//Data-level only (zoom deferred).  out is nBins*N (bin b, component c at out[b*N+c]).
//type sum:  out[b*N+c] = sum over overlapping intervals of value_c * overlap_bases.
//type mean: that / covered bases.  Uncovered bin -> 0 (sum) / NaN (mean).
//vector fork: accumulate per-component sums (+ covered bases) over [start,end) from a zoom level's
//overlapping blocks.  out[c] += rec_sum_c * scalar; *cov += rec_nBases * scalar (scalar = getScalar
//fractional overlap).  Mirrors getVals for the 24+4N record.  Returns 0 on success, 1 on error.
static int accumZoomVec(bigWigFile_t *fp, bwOverlapBlock_t *o, uint32_t tid, uint64_t start, uint64_t end,
                        uint32_t N, double *out, double *cov) {
    size_t zrec = 24 + 4*(size_t)N;
    void *buf = NULL, *compBuf = NULL;
    int compressed = (fp->hdr->bufSize > 0), ret = 1;
    uint64_t i;
    if(compressed) { buf = malloc(fp->hdr->bufSize); if(!buf) return 1; }
    for(i=0; i<o->n; i++) {
        uint8_t *p, *pend;
        uLongf usz;
        if(bwSetPos(fp, o->offset[i])) goto done;
        compBuf = malloc(o->size[i]);
        if(!compBuf) goto done;
        if(bwRead(compBuf, o->size[i], 1, fp) != 1) goto done;
        if(compressed) {
            usz = fp->hdr->bufSize;
            if(uncompress(buf, &usz, compBuf, o->size[i]) != Z_OK) goto done;
        } else { buf = compBuf; usz = o->size[i]; }
        for(p = (uint8_t*)buf, pend = p + usz; p + zrec <= pend; p += zrec) {
            uint32_t vtid = *(uint32_t*)(p+0);
            uint64_t vs = *(uint64_t*)(p+4), ve = *(uint64_t*)(p+12);
            if(vtid == tid) {
                if((start <= vs && end > vs) || (start < ve && start >= vs)) {
                    double scalar = getScalar(start, end, vs, ve);
                    uint32_t c;
                    for(c=0;c<N;c++) out[c] += scalar * (double)(*(float*)(p+24+4*(size_t)c));
                    *cov += scalar * (double)(*(uint32_t*)(p+20));
                }
                if(vs > end) break;
            } else if(vtid > tid) break;
        }
        free(compBuf); compBuf = NULL;
        if(!compressed) buf = NULL;
    }
    ret = 0;
done:
    if(compBuf) free(compBuf);
    if(compressed && buf) free(buf);
    return ret;
}

int bwStatsVec(bigWigFile_t *fp, const char *chrom, uint64_t start, uint64_t end, uint32_t nBins, enum bwStatsType type, double *out) {
    uint32_t b, c, N;
    uint64_t i, pos, end2;
    bwOverlappingIntervalsVec_t *ints;
    if(!fp->vecN || !out || !nBins) return 1;
    if(type != sum && type != mean) return 2;
    N = fp->vecN;

    //vector fork: coarse query -> use the per-component-SUM zoom pyramid (mirrors bwStatsFromZoom)
    if(fp->hdr->nLevels && end > start) {
        int basesPerBin = (int)((end - start)/nBins);
        int32_t level = determineZoomLevel(fp, basesPerBin);
        if(level >= 0) {
            uint32_t tid = bwGetTid(fp, chrom);
            if(tid != (uint32_t)-1) {
                if(!fp->hdr->zoomHdrs->idx[level])
                    fp->hdr->zoomHdrs->idx[level] = bwReadIndex(fp, fp->hdr->zoomHdrs->indexOffset[level]);
                if(fp->hdr->zoomHdrs->idx[level]) {
                    for(b=0, pos=start; b<nBins; b++) {
                        bwOverlapBlock_t *blocks;
                        double cov = 0.0;
                        end2 = start + (uint64_t)(((double)(end-start)*(b+1))/((double)nBins));
                        for(c=0;c<N;c++) out[(size_t)b*N+c] = 0.0;
                        blocks = walkRTreeNodes(fp, fp->hdr->zoomHdrs->idx[level]->root, tid, pos, end2);
                        if(blocks) {
                            int arv = accumZoomVec(fp, blocks, tid, pos, end2, N, out + (size_t)b*N, &cov);
                            destroyBWOverlapBlock(blocks);
                            if(arv) return 3;   //vector fork: I/O/OOM mid-zoom-stat -> fail, don't return partial sums
                        }
                        if(type == mean) {
                            if(cov > 0.0) for(c=0;c<N;c++) out[(size_t)b*N+c] /= cov;
                            else          for(c=0;c<N;c++) out[(size_t)b*N+c] = strtod("NaN", NULL);
                        }
                        pos = end2;
                    }
                    return 0;
                }
            }
        }
    }

    pos = start;
    for(b=0; b<nBins; b++) {
        end2 = start + (uint64_t)(((double)(end-start)*(b+1))/((double) nBins));
        for(c=0; c<N; c++) out[(size_t)b*N + c] = 0.0;
        ints = bwGetOverlappingIntervalsVec(fp, chrom, pos, end2);
        if(ints && ints->l) {
            double cov = 0.0;
            for(i=0; i<ints->l; i++) {
                uint64_t s = ints->start[i], e = ints->end[i];
                if(s < pos) s = pos;
                if(e > end2) e = end2;
                if(e <= s) continue;
                double ov = (double)(e - s);
                cov += ov;
                for(c=0; c<N; c++) out[(size_t)b*N + c] += ov * (double) ints->value[(size_t)i*N + c];
            }
            if(type == mean) {
                if(cov > 0.0) for(c=0; c<N; c++) out[(size_t)b*N + c] /= cov;
                else          for(c=0; c<N; c++) out[(size_t)b*N + c] = strtod("NaN", NULL);
            }
        } else if(type == mean) {
            for(c=0; c<N; c++) out[(size_t)b*N + c] = strtod("NaN", NULL);
        }
        if(ints) bwDestroyOverlappingIntervalsVec(ints);
        pos = end2;
    }
    return 0;
}
