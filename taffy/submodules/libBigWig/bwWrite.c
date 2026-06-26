#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bigWig.h"
#include "bwCommon.h"

/// @cond SKIP
struct val_t {
    uint32_t tid;
    uint32_t start;
    uint32_t nBases;
    float min, max, sum, sumsq;
    double scalar;
    struct val_t *next;
};
/// @endcond

//Create a chromList_t and attach it to a bigWigFile_t *. Returns NULL on error
//Note that chroms and lengths are duplicated, so you MUST free the input
chromList_t *bwCreateChromList(const char* const* chroms, const uint64_t *lengths, int64_t n) {  //64-bit fork: chrom lengths u64
    int64_t i = 0;
    chromList_t *cl = calloc(1, sizeof(chromList_t));
    if(!cl) return NULL;

    cl->nKeys = n;
    cl->chrom = malloc(sizeof(char*)*n);
    cl->len = malloc(sizeof(uint64_t)*n);  //64-bit fork
    if(!cl->chrom) goto error;
    if(!cl->len) goto error;

    for(i=0; i<n; i++) {
        cl->len[i] = lengths[i];
        cl->chrom[i] = bwStrdup(chroms[i]);
        if(!cl->chrom[i]) goto error;
    }

    return cl;

error:
    if(i) {
        int64_t j;
        for(j=0; j<i; j++) free(cl->chrom[j]);
    }
    if(cl) {
        if(cl->chrom) free(cl->chrom);
        if(cl->len) free(cl->len);
        free(cl);
    }
    return NULL;
}

//If maxZooms == 0, then 0 is used (i.e., there are no zoom levels). If maxZooms < 0 or > 65535 then 10 is used.
//TODO allow changing bufSize and blockSize
int bwCreateHdr(bigWigFile_t *fp, int32_t maxZooms) {
    if(!fp->isWrite) return 1;
    bigWigHdr_t *hdr = calloc(1, sizeof(bigWigHdr_t));
    if(!hdr) return 2;

    hdr->version = 4;
    if(maxZooms < 0 || maxZooms > 65535) {
        hdr->nLevels = 10;
    } else {
        hdr->nLevels = maxZooms;
    }

    hdr->bufSize = 32768; //When the file is finalized this is reset if fp->writeBuffer->compressPsz is 0!
    hdr->minVal = DBL_MAX;
    hdr->maxVal = DBL_MIN;
    fp->hdr = hdr;
    fp->writeBuffer->blockSize = 64;

    //Allocate the writeBuffer buffers
    fp->writeBuffer->compressPsz = compressBound(hdr->bufSize);
    fp->writeBuffer->compressP = malloc(fp->writeBuffer->compressPsz);
    if(!fp->writeBuffer->compressP) return 3;
    fp->writeBuffer->p = calloc(1,hdr->bufSize);
    if(!fp->writeBuffer->p) return 4;

    return 0;
}

//64-bit vector fork: like bwCreateHdr but marks the file vector (N float components/record).
//bwWriteHdr then writes BIGWIG64VEC_MAGIC + fieldCount=N.
int bwCreateHdrVec(bigWigFile_t *fp, int32_t maxZooms, uint32_t N) {
    int rv;
    if(N < 1) return 10;
    rv = bwCreateHdr(fp, maxZooms);
    if(rv) return rv;
    //A record is 16 + 4N bytes; it must fit in a data block (after the 32-byte block header).
    if(32 + 16 + 4*(uint64_t)N > fp->hdr->bufSize) return 11;
    fp->vecN = N;
    fp->hdr->fieldCount = (uint16_t) N;
    return 0;
}

//return 0 on success
static int writeAtPos(void *ptr, size_t sz, size_t nmemb, size_t pos, FILE *fp) {
    size_t curpos = ftell(fp);
    if(fseek(fp, pos, SEEK_SET)) return 1;
    if(fwrite(ptr, sz, nmemb, fp) != nmemb) return 2;
    if(fseek(fp, curpos, SEEK_SET)) return 3;
    return 0;
}

//We lose keySize bytes on error
static int writeChromList(FILE *fp, chromList_t *cl) {
    uint16_t k;
    uint32_t j, magic = CIRTREE_MAGIC;
    uint32_t nperblock = (cl->nKeys > 0x7FFF) ? 0x7FFF : cl->nKeys; //Items per leaf/non-leaf, there are no unsigned ints in java :(
    uint32_t nblocks, keySize = 0, valSize = 12; //64-bit fork: value is chromId(u32)+chromSize(u64)=12
    uint64_t i, nonLeafEnd, leafSize, nextLeaf;
    uint8_t eight;
    int64_t i64;
    char *chrom;
    size_t l;

    if(cl->nKeys > 1073676289) {
        fprintf(stderr, "[writeChromList] Error: Currently only 1,073,676,289 contigs are supported. If you really need more then please post a request on github.\n");
        return 1;
    }
    nblocks = cl->nKeys/nperblock;
    nblocks += ((cl->nKeys % nperblock) > 0)?1:0;

    for(i64=0; i64<cl->nKeys; i64++) {
        l = strlen(cl->chrom[i64]);
        if(l>keySize) keySize = l;
    }
    l--; //We don't null terminate strings, because schiess mich tot
    chrom = calloc(keySize, sizeof(char));

    //Write the root node of a largely pointless tree
    if(fwrite(&magic, sizeof(uint32_t), 1, fp) != 1) return 1;
    if(fwrite(&nperblock, sizeof(uint32_t), 1, fp) != 1) return 2;
    if(fwrite(&keySize, sizeof(uint32_t), 1, fp) != 1) return 3;
    if(fwrite(&valSize, sizeof(uint32_t), 1, fp) != 1) return 4;
    if(fwrite(&(cl->nKeys), sizeof(uint64_t), 1, fp) != 1) return 5;

    //Padding?
    i=0;
    if(fwrite(&i, sizeof(uint64_t), 1, fp) != 1) return 6;

    //Do we need a non-leaf node?
    if(nblocks > 1) {
        eight = 0;
        if(fwrite(&eight, sizeof(uint8_t), 1, fp) != 1) return 7;
        if(fwrite(&eight, sizeof(uint8_t), 1, fp) != 1) return 8; //padding
        if(fwrite(&nblocks, sizeof(uint16_t), 1, fp) != 1) return 8;
        nonLeafEnd = ftell(fp) + nperblock * (keySize + 8);
        leafSize = nperblock * (keySize + 8) + 4;
        for(i=0; i<nblocks; i++) { //Why yes, this is pointless
            chrom = strncpy(chrom, cl->chrom[i * nperblock], keySize);
            nextLeaf = nonLeafEnd + i * leafSize;
            if(fwrite(chrom, keySize, 1, fp) != 1) return 9;
            if(fwrite(&nextLeaf, sizeof(uint64_t), 1, fp) != 1) return 10;
        }
        for(i=0; i<keySize; i++) chrom[i] = '\0';
        nextLeaf = 0;
        for(i=nblocks; i<nperblock; i++) {
            if(fwrite(chrom, keySize, 1, fp) != 1) return 9;
            if(fwrite(&nextLeaf, sizeof(uint64_t), 1, fp) != 1) return 10;
        }
    }

    //Write the leaves
    nextLeaf = 0;
    for(i=0, j=0; i<nblocks; i++) {
        eight = 1;
        if(fwrite(&eight, sizeof(uint8_t), 1, fp) != 1) return 11;
        eight = 0;
        if(fwrite(&eight, sizeof(uint8_t), 1, fp) != 1) return 12;
        if(cl->nKeys - j < nperblock) {
            k = cl->nKeys - j;
            if(fwrite(&k, sizeof(uint16_t), 1, fp) != 1) return 13;
        } else {
            if(fwrite(&nperblock, sizeof(uint16_t), 1, fp) != 1) return 13;
        }
        for(k=0; k<nperblock; k++) {
            if(j>=cl->nKeys) {
                if(chrom[0]) {
                    for(l=0; l<keySize; l++) chrom[l] = '\0';
                }
                if(fwrite(chrom, keySize, 1, fp) != 1) return 15;
                if(fwrite(&nextLeaf, sizeof(uint64_t), 1, fp) != 1) return 16;
            } else {
                chrom = strncpy(chrom, cl->chrom[j], keySize);
                if(fwrite(chrom, keySize, 1, fp) != 1) return 15;
                if(fwrite(&j, sizeof(uint32_t), 1, fp) != 1) return 16;
                if(fwrite(&(cl->len[j++]), sizeof(uint64_t), 1, fp) != 1) return 17; //64-bit fork: chromSize u64
            }
        }
    }

    free(chrom);
    return 0;
}

//returns 0 on success
//Still need to fill in indexOffset
int bwWriteHdr(bigWigFile_t *bw) {
    uint32_t magic = bw->vecN ? BIGWIG64VEC_MAGIC : BIGWIG64_MAGIC;  //64-bit fork: non-standard magic (vec variant if vecN)
    uint16_t two = 4;
    FILE *fp;
    const uint8_t pbuff[58] = {0}; // 58 bytes of nothing
    const void *p = (const void *)&pbuff;
    if(!bw->isWrite) return 1;

    //The header itself, largely just reserving space...
    fp = bw->URL->x.fp;
    if(!fp) return 2;
    if(fseek(fp, 0, SEEK_SET)) return 3;
    if(fwrite(&magic, sizeof(uint32_t), 1, fp) != 1) return 4;
    if(fwrite(&two, sizeof(uint16_t), 1, fp) != 1) return 5;
    if(fwrite(p, sizeof(uint8_t), 58, fp) != 58) return 6;
    if(bw->vecN) {  //64-bit vector fork: store N in fieldCount (offset 0x20)
        uint16_t fc = (uint16_t) bw->vecN;
        if(writeAtPos(&fc, sizeof(uint16_t), 1, 0x20, fp)) return 12;
    }

    //Empty zoom headers
    if(bw->hdr->nLevels) {
        for(two=0; two<bw->hdr->nLevels; two++) {
            if(fwrite(p, sizeof(uint8_t), 24, fp) != 24) return 9;
        }
    }

    //Update summaryOffset and write an empty summary block
    bw->hdr->summaryOffset = ftell(fp);
    if(fwrite(p, sizeof(uint8_t), 40, fp) != 40) return 10;
    if(writeAtPos(&(bw->hdr->summaryOffset), sizeof(uint64_t), 1, 0x2c, fp)) return 11;

    //Write the chromosome list as a stupid freaking tree (because let's TREE ALL THE THINGS!!!)
    bw->hdr->ctOffset = ftell(fp);
    if(writeChromList(fp, bw->cl)) return 7;
    if(writeAtPos(&(bw->hdr->ctOffset), sizeof(uint64_t), 1, 0x8, fp)) return 8;

    //Update the dataOffset
    bw->hdr->dataOffset = ftell(fp);
    if(writeAtPos(&bw->hdr->dataOffset, sizeof(uint64_t), 1, 0x10, fp)) return 12;

    //Save space for the number of blocks
    if(fwrite(p, sizeof(uint8_t), 8, fp) != 8) return 13;

    return 0;
}

static int insertIndexNode(bigWigFile_t *fp, bwRTreeNode_t *leaf) {
    bwLL *l = malloc(sizeof(bwLL));
    if(!l) return 1;
    l->node = leaf;
    l->next = NULL;

    if(!fp->writeBuffer->firstIndexNode) {
        fp->writeBuffer->firstIndexNode = l;
    } else {
        fp->writeBuffer->currentIndexNode->next = l;
    }
    fp->writeBuffer->currentIndexNode = l;
    return 0;
}

//0 on success
static int appendIndexNodeEntry(bigWigFile_t *fp, uint32_t tid0, uint32_t tid1, uint64_t start, uint64_t end, uint64_t offset, uint64_t size) {  //64-bit fork: start/end u64 (else entries 2+ in a node truncate baseStart/baseEnd to u32)
    bwLL *n = fp->writeBuffer->currentIndexNode;
    if(!n) return 1;
    if(n->node->nChildren >= fp->writeBuffer->blockSize) return 2;

    n->node->chrIdxStart[n->node->nChildren] = tid0;
    n->node->baseStart[n->node->nChildren] = start;
    n->node->chrIdxEnd[n->node->nChildren] = tid1;
    n->node->baseEnd[n->node->nChildren] = end;
    n->node->dataOffset[n->node->nChildren] = offset;
    n->node->x.size[n->node->nChildren] = size;
    n->node->nChildren++;
    return 0;
}

//Returns 0 on success
static int addIndexEntry(bigWigFile_t *fp, uint32_t tid0, uint32_t tid1, uint64_t start, uint64_t end, uint64_t offset, uint64_t size) {  //64-bit fork: start/end u64 (don't truncate the R-tree feed)
    bwRTreeNode_t *node;

    if(appendIndexNodeEntry(fp, tid0, tid1, start, end, offset, size)) {
        //The last index node is full, we need to add a new one
        node = calloc(1, sizeof(bwRTreeNode_t));
        if(!node) return 1;

        //Allocate and set the fields
        node->isLeaf = 1;
        node->nChildren = 1;
        node->chrIdxStart = malloc(sizeof(uint32_t)*fp->writeBuffer->blockSize);
        if(!node->chrIdxStart) goto error;
        node->baseStart = malloc(sizeof(uint64_t)*fp->writeBuffer->blockSize);  //64-bit fork
        if(!node->baseStart) goto error;
        node->chrIdxEnd = malloc(sizeof(uint32_t)*fp->writeBuffer->blockSize);
        if(!node->chrIdxEnd) goto error;
        node->baseEnd = malloc(sizeof(uint64_t)*fp->writeBuffer->blockSize);  //64-bit fork
        if(!node->baseEnd) goto error;
        node->dataOffset = malloc(sizeof(uint64_t)*fp->writeBuffer->blockSize);
        if(!node->dataOffset) goto error;
        node->x.size = malloc(sizeof(uint64_t)*fp->writeBuffer->blockSize);
        if(!node->x.size) goto error;

        node->chrIdxStart[0] = tid0;
        node->baseStart[0] = start;
        node->chrIdxEnd[0] = tid1;
        node->baseEnd[0] = end;
        node->dataOffset[0] = offset;
        node->x.size[0] = size;

        if(insertIndexNode(fp, node)) goto error;
    }

    return 0;

error:
    if(node->chrIdxStart) free(node->chrIdxStart);
    if(node->baseStart) free(node->baseStart);
    if(node->chrIdxEnd) free(node->chrIdxEnd);
    if(node->baseEnd) free(node->baseEnd);
    if(node->dataOffset) free(node->dataOffset);
    if(node->x.size) free(node->x.size);
    return 2;
}

/*
 * TODO:
 *     The buffer size and compression sz need to be determined elsewhere (and p and compressP filled in!)
 */
static int flushBuffer(bigWigFile_t *fp) {
    bwWriteBuffer_t *wb = fp->writeBuffer;
    uLongf sz = wb->compressPsz;
    uint16_t nItems;
    if(!fp->writeBuffer->l) return 0;
    if(!wb->ltype) return 0;

    //Fill in the header
    //64-bit fork: data-block header is now 32 bytes (start/end widened to u64):
    //tid u32@0, start u64@4, end u64@12, step u32@20, span u32@24, ltype u8@28, pad@29, nItems u16@30
    if(!memcpy((char*)wb->p, &(wb->tid), sizeof(uint32_t))) return 1;
    if(!memcpy((char*)wb->p+4, &(wb->start), sizeof(uint64_t))) return 2;
    if(!memcpy((char*)wb->p+12, &(wb->end), sizeof(uint64_t))) return 3;
    if(!memcpy((char*)wb->p+20, &(wb->step), sizeof(uint32_t))) return 4;
    if(!memcpy((char*)wb->p+24, &(wb->span), sizeof(uint32_t))) return 5;
    if(!memcpy((char*)wb->p+28, &(wb->ltype), sizeof(uint8_t))) return 6;
    //1 byte padding (offset 29)
    //Determine the number of items (32-byte header; records: bedGraph 20, varStep 12, fixedStep 4)
    switch(wb->ltype) {
    case 1:
        nItems = (wb->l-32) / (fp->vecN ? (16 + 4*fp->vecN) : 20);  //64-bit vector fork: stride 16+4N
        break;
    case 2:
        nItems = (wb->l-32)/12;
        break;
    case 3:
        nItems = (wb->l-32)/4;
        break;
    default:
        return 7;
    }
    if(!memcpy((char*)wb->p+30, &nItems, sizeof(uint16_t))) return 8;

    if(sz) {
        //compress
        if(compress(wb->compressP, &sz, wb->p, wb->l) != Z_OK) return 9;

        //write the data to disk
        if(fwrite(wb->compressP, sizeof(uint8_t), sz, fp->URL->x.fp) != sz) return 10;
    } else {
        sz = wb->l;
        if(fwrite(wb->p, sizeof(uint8_t), wb->l, fp->URL->x.fp) != wb->l) return 10;
    }

    //Add an entry into the index
    if(addIndexEntry(fp, wb->tid, wb->tid, wb->start, wb->end, bwTell(fp)-sz, sz)) return 11;

    wb->nBlocks++;
    wb->l = 32;   //64-bit fork: header is 32 bytes
    return 0;
}

static void updateStats(bigWigFile_t *fp, uint32_t span, float val) {
    if(val < fp->hdr->minVal) fp->hdr->minVal = val;
    else if(val > fp->hdr->maxVal) fp->hdr->maxVal = val;
    fp->hdr->nBasesCovered += span;
    fp->hdr->sumData += span*val;
    fp->hdr->sumSquared += span*pow(val,2);

    fp->writeBuffer->nEntries++;
    fp->writeBuffer->runningWidthSum += span;
}

//12 bytes per entry
int bwAddIntervals(bigWigFile_t *fp, const char* const* chrom, const uint64_t *start, const uint64_t *end, const float *values, uint32_t n) {
    uint32_t tid = 0, i;
    const char *lastChrom = NULL;
    bwWriteBuffer_t *wb = fp->writeBuffer;
    if(!n) return 0; //Not an error per se
    if(!fp->isWrite) return 1;
    if(!wb) return 2;

    //Flush if needed
    if(wb->ltype != 1) if(flushBuffer(fp)) return 3;
    if(wb->l+52 > fp->hdr->bufSize) if(flushBuffer(fp)) return 4;  //64-bit fork: 32 hdr + 20 rec
    lastChrom = chrom[0];
    tid = bwGetTid(fp, chrom[0]);
    if(tid == (uint32_t) -1) return 5;
    if(tid != wb->tid) {
        if(flushBuffer(fp)) return 6;
        wb->tid = tid;
        wb->start = start[0];
        wb->end = end[0];
    }

    //Ensure that everything is set correctly
    wb->ltype = 1;
    if(wb->l <= 32) {
        wb->l = 32;   //64-bit fork: reserve the 32-byte data header for a fresh block
        wb->start = start[0];
        wb->span = 0;
        wb->step = 0;
    }
    if(!memcpy((char*)wb->p+wb->l, start, sizeof(uint64_t))) return 7;
    if(!memcpy((char*)wb->p+wb->l+8, end, sizeof(uint64_t))) return 8;
    if(!memcpy((char*)wb->p+wb->l+16, values, sizeof(float))) return 9;
    updateStats(fp, end[0]-start[0], values[0]);
    wb->l += 20;

    for(i=1; i<n; i++) {
        if(strcmp(chrom[i],lastChrom) != 0) {
            wb->end = end[i-1];
            flushBuffer(fp);
            lastChrom = chrom[i];
            tid = bwGetTid(fp, chrom[i]);
            if(tid == (uint32_t) -1) return 10;
            wb->tid = tid;
            wb->start = start[i];
        }
        if(wb->l+20 > fp->hdr->bufSize) { //20 bytes/entry (64-bit fork)
            wb->end = end[i-1];
            flushBuffer(fp);
            wb->start = start[i];
        }
        if(!memcpy((char*)wb->p+wb->l, &(start[i]), sizeof(uint64_t))) return 11;
        if(!memcpy((char*)wb->p+wb->l+8, &(end[i]), sizeof(uint64_t))) return 12;
        if(!memcpy((char*)wb->p+wb->l+16, &(values[i]), sizeof(float))) return 13;
        updateStats(fp, end[i]-start[i], values[i]);
        wb->l += 20;
    }
    wb->end = end[i-1];

    return 0;
}

int bwAppendIntervals(bigWigFile_t *fp, const uint64_t *start, const uint64_t *end, const float *values, uint32_t n) {
    uint32_t i;
    bwWriteBuffer_t *wb = fp->writeBuffer;
    if(!n) return 0;
    if(!fp->isWrite) return 1;
    if(!wb) return 2;
    if(wb->ltype != 1) return 3;

    for(i=0; i<n; i++) {
        if(wb->l+20 > fp->hdr->bufSize) {  //64-bit fork (B3): bedGraph record is 20B
            if(i>0) { //otherwise it's already set
                wb->end = end[i-1];
            }
            flushBuffer(fp);
            wb->start = start[i];
        }
        if(!memcpy((char*)wb->p+wb->l, &(start[i]), sizeof(uint64_t))) return 4;
        if(!memcpy((char*)wb->p+wb->l+8, &(end[i]), sizeof(uint64_t))) return 5;
        if(!memcpy((char*)wb->p+wb->l+16, &(values[i]), sizeof(float))) return 6;
        updateStats(fp, end[i]-start[i], values[i]);
        wb->l += 20;
    }
    wb->end = end[i-1];

    return 0;
}

//64-bit vector fork: type-1 bedGraph packer with N float components/record (stride 16+4N).
//values is n*N flattened (interval i, component c at values[i*N+c]); sorted, non-overlapping.
int bwAddIntervalsVec(bigWigFile_t *fp, const char* const* chrom, const uint64_t *start, const uint64_t *end, const float *values, uint32_t n) {
    uint32_t tid = 0, i, N, rec;
    const char *lastChrom = NULL;
    bwWriteBuffer_t *wb = fp->writeBuffer;
    if(!n) return 0;
    if(!fp->isWrite) return 1;
    if(!wb) return 2;
    if(!fp->vecN) return 14;  //not a vector file
    N = fp->vecN; rec = 16 + 4*N;

    if(wb->ltype != 1) if(flushBuffer(fp)) return 3;
    if(wb->l + 32 + rec > fp->hdr->bufSize) if(flushBuffer(fp)) return 4;
    lastChrom = chrom[0];
    tid = bwGetTid(fp, chrom[0]);
    if(tid == (uint32_t) -1) return 5;
    if(tid != wb->tid) {
        if(flushBuffer(fp)) return 6;
        wb->tid = tid;
        wb->start = start[0];
        wb->end = end[0];
    }

    wb->ltype = 1;
    if(wb->l <= 32) {
        wb->l = 32;   //reserve the 32-byte block header
        wb->start = start[0];
        wb->span = 0;
        wb->step = 0;
    }
    if(!memcpy((char*)wb->p+wb->l, start, sizeof(uint64_t))) return 7;
    if(!memcpy((char*)wb->p+wb->l+8, end, sizeof(uint64_t))) return 8;
    if(!memcpy((char*)wb->p+wb->l+16, values, sizeof(float)*N)) return 9;
    updateStats(fp, end[0]-start[0], values[0]);
    wb->l += rec;

    for(i=1; i<n; i++) {
        if(strcmp(chrom[i],lastChrom) != 0) {
            wb->end = end[i-1];
            flushBuffer(fp);
            lastChrom = chrom[i];
            tid = bwGetTid(fp, chrom[i]);
            if(tid == (uint32_t) -1) return 10;
            wb->tid = tid;
            wb->start = start[i];
        }
        if(wb->l+rec > fp->hdr->bufSize) {
            wb->end = end[i-1];
            flushBuffer(fp);
            wb->start = start[i];
        }
        if(!memcpy((char*)wb->p+wb->l, &(start[i]), sizeof(uint64_t))) return 11;
        if(!memcpy((char*)wb->p+wb->l+8, &(end[i]), sizeof(uint64_t))) return 12;
        if(!memcpy((char*)wb->p+wb->l+16, values + (size_t)i*N, sizeof(float)*N)) return 13;
        updateStats(fp, end[i]-start[i], values[(size_t)i*N]);
        wb->l += rec;
    }
    wb->end = end[i-1];
    return 0;
}

int bwAppendIntervalsVec(bigWigFile_t *fp, const uint64_t *start, const uint64_t *end, const float *values, uint32_t n) {
    uint32_t i, N, rec;
    bwWriteBuffer_t *wb = fp->writeBuffer;
    if(!n) return 0;
    if(!fp->isWrite) return 1;
    if(!wb) return 2;
    if(wb->ltype != 1) return 3;
    if(!fp->vecN) return 7;
    N = fp->vecN; rec = 16 + 4*N;

    for(i=0; i<n; i++) {
        if(wb->l+rec > fp->hdr->bufSize) {
            if(i>0) wb->end = end[i-1];
            flushBuffer(fp);
            wb->start = start[i];
        }
        if(!memcpy((char*)wb->p+wb->l, &(start[i]), sizeof(uint64_t))) return 4;
        if(!memcpy((char*)wb->p+wb->l+8, &(end[i]), sizeof(uint64_t))) return 5;
        if(!memcpy((char*)wb->p+wb->l+16, values + (size_t)i*N, sizeof(float)*N)) return 6;
        updateStats(fp, end[i]-start[i], values[(size_t)i*N]);
        wb->l += rec;
    }
    wb->end = end[i-1];
    return 0;
}

//8 bytes per entry
int bwAddIntervalSpans(bigWigFile_t *fp, const char *chrom, const uint64_t *start, uint32_t span, const float *values, uint32_t n) {
    uint32_t i, tid;
    bwWriteBuffer_t *wb = fp->writeBuffer;
    if(!n) return 0;
    if(!fp->isWrite) return 1;
    if(!wb) return 2;
    if(wb->ltype != 2) if(flushBuffer(fp)) return 3;
    if(flushBuffer(fp)) return 4;

    tid = bwGetTid(fp, chrom);
    if(tid == (uint32_t) -1) return 5;
    wb->tid = tid;
    wb->start = start[0];
    wb->step = 0;
    wb->span = span;
    wb->ltype = 2;

    for(i=0; i<n; i++) {
        if(wb->l + 12 >= fp->hdr->bufSize) { //64-bit fork (B2): varStep record is 12B
            if(i) wb->end = start[i-1]+span;
            flushBuffer(fp);
            wb->start = start[i];
        }
        if(!memcpy((char*)wb->p+wb->l, &(start[i]), sizeof(uint64_t))) return 5;
        if(!memcpy((char*)wb->p+wb->l+8, &(values[i]), sizeof(float))) return 6;
        updateStats(fp, span, values[i]);
        wb->l += 12;
    }
    wb->end = start[n-1] + span;

    return 0;
}

int bwAppendIntervalSpans(bigWigFile_t *fp, const uint64_t *start, const float *values, uint32_t n) {
    uint32_t i;
    bwWriteBuffer_t *wb = fp->writeBuffer;
    if(!n) return 0;
    if(!fp->isWrite) return 1;
    if(!wb) return 2;
    if(wb->ltype != 2) return 3;

    for(i=0; i<n; i++) {
        if(wb->l + 12 >= fp->hdr->bufSize) {  //64-bit fork (B2): varStep record is 12B
            if(i) wb->end = start[i-1]+wb->span;
            flushBuffer(fp);
            wb->start = start[i];
        }
        if(!memcpy((char*)wb->p+wb->l, &(start[i]), sizeof(uint64_t))) return 4;
        if(!memcpy((char*)wb->p+wb->l+8, &(values[i]), sizeof(float))) return 5;
        updateStats(fp, wb->span, values[i]);
        wb->l += 12;
    }
    wb->end = start[n-1] + wb->span;

    return 0;
}

//4 bytes per entry
int bwAddIntervalSpanSteps(bigWigFile_t *fp, const char *chrom, uint64_t start, uint32_t span, uint32_t step, const float *values, uint32_t n) {
    uint32_t i, tid;
    bwWriteBuffer_t *wb = fp->writeBuffer;
    if(!n) return 0;
    if(!fp->isWrite) return 1;
    if(!wb) return 2;
    if(wb->ltype != 3) flushBuffer(fp);
    if(flushBuffer(fp)) return 3;

    tid = bwGetTid(fp, chrom);
    if(tid == (uint32_t) -1) return 4;
    wb->tid = tid;
    wb->start = start;
    wb->step = step;
    wb->span = span;
    wb->ltype = 3;

    for(i=0; i<n; i++) {
        if(wb->l + 4 >= fp->hdr->bufSize) {
            wb->end = wb->start + ((wb->l-32)>>2) * step;  //64-bit fork (B4): 32B header
            flushBuffer(fp);
            wb->start = wb->end;
        }
        if(!memcpy((char*)wb->p+wb->l, &(values[i]), sizeof(float))) return 5;
        updateStats(fp, wb->span, values[i]);
        wb->l += 4;
    }
    wb->end = wb->start + ((wb->l-32)>>2) * step;  //64-bit fork (B4): 32B header

    return 0;
}

int bwAppendIntervalSpanSteps(bigWigFile_t *fp, const float *values, uint32_t n) {
    uint32_t i;
    bwWriteBuffer_t *wb = fp->writeBuffer;
    if(!n) return 0;
    if(!fp->isWrite) return 1;
    if(!wb) return 2;
    if(wb->ltype != 3) return 3;

    for(i=0; i<n; i++) {
        if(wb->l + 4 >= fp->hdr->bufSize) {
            wb->end = wb->start + ((wb->l-32)>>2) * wb->step;  //64-bit fork (B4): 32B header
            flushBuffer(fp);
            wb->start = wb->end;
        }
        if(!memcpy((char*)wb->p+wb->l, &(values[i]), sizeof(float))) return 4;
        updateStats(fp, wb->span, values[i]);
        wb->l += 4;
    }
    wb->end = wb->start + ((wb->l-32)>>2) * wb->step;  //64-bit fork (B4): 32B header

    return 0;
}

//0 on success
int writeSummary(bigWigFile_t *fp) {
    if(writeAtPos(&(fp->hdr->nBasesCovered), sizeof(uint64_t), 1, fp->hdr->summaryOffset, fp->URL->x.fp)) return 1;
    if(writeAtPos(&(fp->hdr->minVal), sizeof(double), 1, fp->hdr->summaryOffset+8, fp->URL->x.fp)) return 2;
    if(writeAtPos(&(fp->hdr->maxVal), sizeof(double), 1, fp->hdr->summaryOffset+16, fp->URL->x.fp)) return 3;
    if(writeAtPos(&(fp->hdr->sumData), sizeof(double), 1, fp->hdr->summaryOffset+24, fp->URL->x.fp)) return 4;
    if(writeAtPos(&(fp->hdr->sumSquared), sizeof(double), 1, fp->hdr->summaryOffset+32, fp->URL->x.fp)) return 5;
    return 0;
}

static bwRTreeNode_t *makeEmptyNode(uint32_t blockSize) {
    bwRTreeNode_t *n = calloc(1, sizeof(bwRTreeNode_t));
    if(!n) return NULL;

    n->chrIdxStart = malloc(blockSize*sizeof(uint32_t));
    if(!n->chrIdxStart) goto error;
    n->baseStart = malloc(blockSize*sizeof(uint64_t));  //64-bit fork: positions are u64
    if(!n->baseStart) goto error;
    n->chrIdxEnd = malloc(blockSize*sizeof(uint32_t));
    if(!n->chrIdxEnd) goto error;
    n->baseEnd = malloc(blockSize*sizeof(uint64_t));  //64-bit fork: positions are u64
    if(!n->baseEnd) goto error;
    n->dataOffset = calloc(blockSize,sizeof(uint64_t)); //This MUST be 0 for node writing!
    if(!n->dataOffset) goto error;
    n->x.child = malloc(blockSize*sizeof(uint64_t));
    if(!n->x.child) goto error;

    return n;

error:
    if(n->chrIdxStart) free(n->chrIdxStart);
    if(n->baseStart) free(n->baseStart);
    if(n->chrIdxEnd) free(n->chrIdxEnd);
    if(n->baseEnd) free(n->baseEnd);
    if(n->dataOffset) free(n->dataOffset);
    if(n->x.child) free(n->x.child);
    free(n);
    return NULL;
}

//Returns 0 on success. This doesn't attempt to clean up!
static bwRTreeNode_t *addLeaves(bwLL **ll, uint64_t *sz, uint64_t toProcess, uint32_t blockSize) {
    uint32_t i;
    uint64_t foo;
    bwRTreeNode_t *n = makeEmptyNode(blockSize);
    if(!n) return NULL;

    if(toProcess <= blockSize) {
        for(i=0; i<toProcess; i++) {
            n->chrIdxStart[i] = (*ll)->node->chrIdxStart[0];
            n->baseStart[i] = (*ll)->node->baseStart[0];
            n->chrIdxEnd[i] = (*ll)->node->chrIdxEnd[(*ll)->node->nChildren-1];
            n->baseEnd[i] = (*ll)->node->baseEnd[(*ll)->node->nChildren-1];
            n->x.child[i] = (*ll)->node;
            *sz += 4 + 32*(*ll)->node->nChildren;
            *ll = (*ll)->next;
            n->nChildren++;
        }
    } else {
        for(i=0; i<blockSize; i++) {
            foo = ceil(((double) toProcess)/((double) blockSize-i));
            if(!ll) break;
            n->x.child[i] = addLeaves(ll, sz, foo, blockSize);
            if(!n->x.child[i]) goto error;
            n->chrIdxStart[i] = n->x.child[i]->chrIdxStart[0];
            n->baseStart[i] = n->x.child[i]->baseStart[0];
            n->chrIdxEnd[i] = n->x.child[i]->chrIdxEnd[n->x.child[i]->nChildren-1];
            n->baseEnd[i] = n->x.child[i]->baseEnd[n->x.child[i]->nChildren-1];
            n->nChildren++;
            toProcess -= foo;
        }
    }

    *sz += 4 + 32*n->nChildren;  //64-bit fork: twig child now 32B (leaf 40B; idxSize is an estimate, as upstream)
    return n;

error:
    bwDestroyIndexNode(n);
    return NULL;
}

//Returns 1 on error
int writeIndexTreeNode(FILE *fp, bwRTreeNode_t *n, uint8_t *wrote, int level) {
    uint8_t one = 0;
    uint32_t i, j;            //64-bit fork: positions written individually (baseStart/baseEnd are u64)
    uint64_t zero64 = 0;

    if(n->isLeaf) return 0;

    for(i=0; i<n->nChildren; i++) {
        if(n->dataOffset[i]) { //traverse into child
            if(n->isLeaf) return 0; //Only write leaves once!
            if(writeIndexTreeNode(fp, n->x.child[i], wrote, level+1)) return 1;
        } else {
            n->dataOffset[i] = ftell(fp);
            if(fwrite(&(n->x.child[i]->isLeaf), sizeof(uint8_t), 1, fp) != 1) return 1;
            if(fwrite(&one, sizeof(uint8_t), 1, fp) != 1) return 1; //one byte of padding
            if(fwrite(&(n->x.child[i]->nChildren), sizeof(uint16_t), 1, fp) != 1) return 1;
            for(j=0; j<n->x.child[i]->nChildren; j++) {
                //64-bit fork: chrIdxStart(u32) baseStart(u64) chrIdxEnd(u32) baseEnd(u64) = 24 bytes
                if(fwrite(&(n->x.child[i]->chrIdxStart[j]), sizeof(uint32_t), 1, fp) != 1) return 1;
                if(fwrite(&(n->x.child[i]->baseStart[j]), sizeof(uint64_t), 1, fp) != 1) return 1;
                if(fwrite(&(n->x.child[i]->chrIdxEnd[j]), sizeof(uint32_t), 1, fp) != 1) return 1;
                if(fwrite(&(n->x.child[i]->baseEnd[j]), sizeof(uint64_t), 1, fp) != 1) return 1;
                if(n->x.child[i]->isLeaf) {
                    //leaf child = 24 + dataOffset(u64) + size(u64) = 40 bytes
                    if(fwrite(&(n->x.child[i]->dataOffset[j]), sizeof(uint64_t), 1, fp) != 1) return 1;
                    if(fwrite(&(n->x.child[i]->x.size[j]), sizeof(uint64_t), 1, fp) != 1) return 1;
                } else {
                    //twig child = 24 + child-offset placeholder(u64) = 32 bytes
                    if(fwrite(&zero64, sizeof(uint64_t), 1, fp) != 1) return 1;
                }
            }
            *wrote = 1;
        }
    }

    return 0;
}

//returns 1 on success
int writeIndexOffsets(FILE *fp, bwRTreeNode_t *n, uint64_t offset) {
    uint32_t i;

    if(n->isLeaf) return 0;
    for(i=0; i<n->nChildren; i++) {
        if(writeIndexOffsets(fp, n->x.child[i], n->dataOffset[i])) return 1;
        if(writeAtPos(&(n->dataOffset[i]), sizeof(uint64_t), 1, offset+28+32*i, fp)) return 2; //64-bit fork: twig child 32B, child-offset@24 within child
    }
    return 0;
}

//Returns 0 on success
int writeIndexTree(bigWigFile_t *fp) {
    uint64_t offset;
    uint8_t wrote = 0;
    int rv;

    while((rv = writeIndexTreeNode(fp->URL->x.fp, fp->idx->root, &wrote, 0)) == 0) {
        if(!wrote) break;
        wrote = 0;
    }

    if(rv || wrote) return 1;

    //Save the file position
    offset = bwTell(fp);

    //Write the offsets
    if(writeIndexOffsets(fp->URL->x.fp, fp->idx->root, fp->idx->rootOffset)) return 2;

    //Move the file pointer back to the end
    bwSetPos(fp, offset);

    return 0;
}

//Returns 0 on success. The original state SHOULD be preserved on error
int writeIndex(bigWigFile_t *fp) {
    uint32_t four = IDX_MAGIC;
    uint64_t idxSize = 0, foo;
    uint8_t one = 0;
    uint32_t i;  //64-bit fork: root node written field-by-field (no u32 vector)
    bwLL *ll = fp->writeBuffer->firstIndexNode, *p;
    bwRTreeNode_t *root = NULL;

    if(!fp->writeBuffer->nBlocks) return 0;
    fp->idx = malloc(sizeof(bwRTree_t));
    if(!fp->idx) return 2;
    fp->idx->root = root;

    //Update the file header to indicate the proper index position
    foo = bwTell(fp);
    if(writeAtPos(&foo, sizeof(uint64_t), 1, 0x18, fp->URL->x.fp)) return 3;

    //Make the tree
    if(ll == fp->writeBuffer->currentIndexNode) {
        root = ll->node;
        idxSize = 4 + 32*root->nChildren;  //64-bit fork: twig child 32B
    } else {
        root = addLeaves(&ll, &idxSize, ceil(((double)fp->writeBuffer->nBlocks)/fp->writeBuffer->blockSize), fp->writeBuffer->blockSize);
    }
    if(!root) return 4;
    fp->idx->root = root;

    ll = fp->writeBuffer->firstIndexNode;
    while(ll) {
        p = ll->next;
        free(ll);
        ll=p;
    }

    //write the header
    if(fwrite(&four, sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 5;
    if(fwrite(&(fp->writeBuffer->blockSize), sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 6;
    if(fwrite(&(fp->writeBuffer->nBlocks), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 7;
    if(fwrite(&(root->chrIdxStart[0]), sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 8;
    if(fwrite(&(root->baseStart[0]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 9;  //64-bit fork (L1)
    if(fwrite(&(root->chrIdxEnd[root->nChildren-1]), sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 10;
    if(fwrite(&(root->baseEnd[root->nChildren-1]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 11;  //64-bit fork (L1)
    if(fwrite(&idxSize, sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 12;
    four = 1;
    if(fwrite(&four, sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 13;
    four = 0;
    if(fwrite(&four, sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 14; //padding
    fp->idx->rootOffset = bwTell(fp);

    //Write the root node, since writeIndexTree writes the children and fills in the offset
    if(fwrite(&(root->isLeaf), sizeof(uint8_t), 1, fp->URL->x.fp) != 1) return 16;
    if(fwrite(&one, sizeof(uint8_t), 1, fp->URL->x.fp) != 1) return 17; //one byte of padding
    if(fwrite(&(root->nChildren), sizeof(uint16_t), 1, fp->URL->x.fp) != 1) return 18;
    for(i=0; i<root->nChildren; i++) {
        //64-bit fork: per child = chrIdxStart(u32) baseStart(u64) chrIdxEnd(u32) baseEnd(u64) = 24B,
        //then leaf: dataOffset(u64)+size(u64) [40B]; twig: child-offset placeholder(u64) [32B]
        if(fwrite(&(root->chrIdxStart[i]), sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 19;
        if(fwrite(&(root->baseStart[i]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 19;
        if(fwrite(&(root->chrIdxEnd[i]), sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 19;
        if(fwrite(&(root->baseEnd[i]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 19;
        if(root->isLeaf) {
            if(fwrite(&(root->dataOffset[i]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 20;
            if(fwrite(&(root->x.size[i]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 21;
        } else {
            uint64_t zero64 = 0;
            root->dataOffset[i] = 0;
            if(fwrite(&zero64, sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 22;
        }
    }

    //Write each level
    if(writeIndexTree(fp)) return 23;

    return 0;
}

//The first zoom level has a resolution of 4x mean entry size
//This may or may not produce the requested number of zoom levels
int makeZoomLevels(bigWigFile_t *fp) {
    uint32_t meanBinSize, i;
    uint32_t multiplier = 4;
    uint64_t zoom = 10, maxZoom = 0;  //64-bit fork (H1): chrom lengths are u64
    uint16_t nLevels = 0;
    size_t zrec = fp->vecN ? (24 + 4*(size_t)fp->vecN) : 40;  //vector fork: zoom-summary record size (24+4N vs scalar 40)

    meanBinSize = ((double) fp->writeBuffer->runningWidthSum)/(fp->writeBuffer->nEntries);
    //In reality, one level is skipped
    meanBinSize *= 4;
    //N.B., we must ALWAYS check that the zoom doesn't overflow a uint32_t!
    if(((uint32_t)-1)>>2 < meanBinSize) return 0; //No zoom levels!
    if(meanBinSize*4 > zoom) zoom = multiplier*meanBinSize;

    fp->hdr->zoomHdrs = calloc(1, sizeof(bwZoomHdr_t));
    if(!fp->hdr->zoomHdrs) return 1;
    fp->hdr->zoomHdrs->level = malloc(fp->hdr->nLevels * sizeof(uint32_t));
    fp->hdr->zoomHdrs->dataOffset = calloc(fp->hdr->nLevels, sizeof(uint64_t));
    fp->hdr->zoomHdrs->indexOffset = calloc(fp->hdr->nLevels, sizeof(uint64_t));
    fp->hdr->zoomHdrs->idx = calloc(fp->hdr->nLevels, sizeof(bwRTree_t*));
    if(!fp->hdr->zoomHdrs->level) return 2;
    if(!fp->hdr->zoomHdrs->dataOffset) return 3;
    if(!fp->hdr->zoomHdrs->indexOffset) return 4;
    if(!fp->hdr->zoomHdrs->idx) return 5;

    //There's no point in having a zoom level larger than the largest chromosome
    //This will none the less allow at least one zoom level, which is generally needed for IGV et al.
    for(i=0; i<fp->cl->nKeys; i++) {
        if(fp->cl->len[i] > maxZoom) maxZoom = fp->cl->len[i];
    }
    if(zoom > maxZoom) zoom = maxZoom;

    for(i=0; i<fp->hdr->nLevels; i++) {
        if(zoom > maxZoom) break; //prevent absurdly large zoom levels
        fp->hdr->zoomHdrs->level[i] = zoom;
        nLevels++;
        if(((uint32_t)-1)/multiplier < zoom) break;
        zoom *= multiplier;
    }
    fp->hdr->nLevels = nLevels;

    fp->writeBuffer->firstZoomBuffer = calloc(nLevels,sizeof(bwZoomBuffer_t*));
    if(!fp->writeBuffer->firstZoomBuffer) goto error;
    fp->writeBuffer->lastZoomBuffer = calloc(nLevels,sizeof(bwZoomBuffer_t*));
    if(!fp->writeBuffer->lastZoomBuffer) goto error;
    fp->writeBuffer->nNodes = calloc(nLevels, sizeof(uint64_t));

    for(i=0; i<fp->hdr->nLevels; i++) {
        fp->writeBuffer->firstZoomBuffer[i] = calloc(1, sizeof(bwZoomBuffer_t));
        if(!fp->writeBuffer->firstZoomBuffer[i]) goto error;
        //64-bit fork: zoom records are 40 bytes scalar (chromId u32@0, start u64@4, end u64@12, ...)
        //or 24+4N bytes vector (...nBases u32@20, N float sums @24); zrec selects.
        fp->writeBuffer->firstZoomBuffer[i]->p = calloc(fp->hdr->bufSize/zrec, zrec);
        if(!fp->writeBuffer->firstZoomBuffer[i]->p) goto error;
        fp->writeBuffer->firstZoomBuffer[i]->m = (fp->hdr->bufSize/zrec)*zrec;
        {
            uint8_t *zp = (uint8_t*)fp->writeBuffer->firstZoomBuffer[i]->p;
            *(uint32_t*)(zp+0) = 0;   //chromId
            *(uint64_t*)(zp+4) = 0;   //start
            uint64_t e0 = fp->hdr->zoomHdrs->level[i];
            if(e0 > fp->cl->len[0]) e0 = fp->cl->len[0];
            *(uint64_t*)(zp+12) = e0; //end
        }
        fp->writeBuffer->lastZoomBuffer[i] =  fp->writeBuffer->firstZoomBuffer[i];
    }

    return 0;

error:
    if(fp->writeBuffer->firstZoomBuffer) {
        for(i=0; i<fp->hdr->nLevels; i++) {
            if(fp->writeBuffer->firstZoomBuffer[i]) {
                if(fp->writeBuffer->firstZoomBuffer[i]->p) free(fp->writeBuffer->firstZoomBuffer[i]->p);
                free(fp->writeBuffer->firstZoomBuffer[i]);
            }
        }
        free(fp->writeBuffer->firstZoomBuffer);
    }
    if(fp->writeBuffer->lastZoomBuffer) free(fp->writeBuffer->lastZoomBuffer);
    if(fp->writeBuffer->nNodes) free(fp->writeBuffer->lastZoomBuffer);
    return 6;
}

//Given an interval start, calculate the next one at a zoom level
void nextPos(bigWigFile_t *fp, uint32_t size, uint32_t *pos, uint32_t desiredTid) {
    uint32_t *tid = pos;
    uint32_t *start = pos+1;
    uint32_t *end = pos+2;
    *start += size;
    if(*start >= fp->cl->len[*tid]) {
        (*start) = 0;
        (*tid)++;
    }

    //prevent needless iteration when changing chromosomes
    if(*tid < desiredTid) {
        *tid = desiredTid;
        *start = 0;
    }

    (*end) = *start+size;
    if(*end > fp->cl->len[*tid]) (*end) = fp->cl->len[*tid];
}

//Return the number of bases two intervals overlap
uint32_t overlapsInterval(uint32_t tid0, uint64_t start0, uint64_t end0, uint32_t tid1, uint64_t start1, uint64_t end1) {  //64-bit fork: positions u64
    if(tid0 != tid1) return 0;
    if(end0 <= start1) return 0;
    if(end1 <= start0) return 0;
    if(end0 <= end1) {
        if(start1 > start0) return end0-start1;
        return end0-start0;
    } else {
        if(start1 > start0) return end1-start1;
        return end1-start0;
    }
}

//Returns the number of bases of the interval written
//64-bit fork: zoom-summary record is 40 bytes -- chromId u32@0, start u64@4,
//end u64@12, nBases u32@20, min/max/sum/sumsq float@24/28/32/36. Positions are
//u64; the bin size + the per-call overlap (rv) stay u32 (bounded for sane levels).
uint32_t updateInterval(bigWigFile_t *fp, bwZoomBuffer_t *buffer, double *sum, double *sumsq, uint32_t size, uint32_t tid, uint64_t start, uint64_t end, float value) {
    uint8_t *b = (uint8_t*) buffer->p;
    uint32_t rv = 0, offset = 0;
    if(!buffer) return 0;
    if(buffer->l+40 >= buffer->m) return 0;
    #define ZC(o)   (*(uint32_t*)(b + (size_t)(o)*40 + 0))
    #define ZS(o)   (*(uint64_t*)(b + (size_t)(o)*40 + 4))
    #define ZE(o)   (*(uint64_t*)(b + (size_t)(o)*40 + 12))
    #define ZN(o)   (*(uint32_t*)(b + (size_t)(o)*40 + 20))
    #define ZMIN(o) (*(float*)(b + (size_t)(o)*40 + 24))
    #define ZMAX(o) (*(float*)(b + (size_t)(o)*40 + 28))
    #define ZSUM(o) (*(float*)(b + (size_t)(o)*40 + 32))
    #define ZSQ(o)  (*(float*)(b + (size_t)(o)*40 + 36))

    if(buffer->l) {
        offset = buffer->l/40;
    } else {
        ZC(0) = tid;
        ZS(0) = start;
        ZE(0) = (start+size < end) ? start+size : end;
    }

    //Do we have any overlap with the previously added interval?
    if(offset) {
        rv = overlapsInterval(ZC(offset-1), ZS(offset-1), ZS(offset-1) + size, tid, start, end);
        if(rv) {
            ZE(offset-1) = start + rv;
            ZN(offset-1) += rv;
            if(ZMIN(offset-1) > value) ZMIN(offset-1) = value;
            if(ZMAX(offset-1) < value) ZMAX(offset-1) = value;
            *sum += rv*value;
            *sumsq += rv*pow(value, 2.0);
            return rv;
        } else {
            ZSUM(offset-1) = *sum;
            ZSQ(offset-1) = *sumsq;
            *sum = 0.0;
            *sumsq = 0.0;
        }
    }

    //If we move to a new interval then skip obviously non-overlapping intervals
    if(offset && ZE(offset) == 0) {
        ZC(offset) = tid;
        ZS(offset) = start;
        ZE(offset) = (start+size < end) ? start+size : end;
    }

    //Add a new entry
    while(!(rv = overlapsInterval(ZC(offset), ZS(offset), ZS(offset) + size, tid, start, end))) {
        ZC(offset) = tid;
        ZS(offset) = start;
        ZE(offset) = (start+size < end) ? start+size : end;
    }
    ZN(offset) = rv;
    ZMIN(offset) = value; //min
    ZMAX(offset) = value; //max
    *sum += rv * value;
    *sumsq += rv * pow(value,2.0);
    buffer->l += 40;
    return rv;
    #undef ZC
    #undef ZS
    #undef ZE
    #undef ZN
    #undef ZMIN
    #undef ZMAX
    #undef ZSUM
    #undef ZSQ
}

//Returns 0 on success
int addIntervalValue(bigWigFile_t *fp, uint64_t *nEntries, double *sum, double *sumsq, bwZoomBuffer_t *buffer, uint32_t itemsPerSlot, uint32_t zoom, uint32_t tid, uint64_t start, uint64_t end, float value) {  //64-bit fork: start/end u64
    bwZoomBuffer_t *newBuffer = NULL;
    uint32_t rv;

    while(start < end) {
        rv = updateInterval(fp, buffer, sum, sumsq, zoom, tid, start, end, value);
        if(!rv) {
            //Allocate a new buffer
            newBuffer = calloc(1, sizeof(bwZoomBuffer_t));
            if(!newBuffer) return 1;
            newBuffer->p = calloc(itemsPerSlot, 40);  //64-bit fork: 40B records
            if(!newBuffer->p) goto error;
            newBuffer->m = itemsPerSlot*40;
            //carry the last record's chromId(@0,4B) + start(@4,8B u64); end(@12) = start + zoom
            memcpy(newBuffer->p, (unsigned char*)buffer->p+buffer->l-40, 4);
            memcpy((unsigned char*)newBuffer->p+4, (unsigned char*)buffer->p + buffer->l-36, 8);
            *(uint64_t*)((unsigned char*)newBuffer->p+12) = *(uint64_t*)((unsigned char*)newBuffer->p+4) + zoom;
            *sum = *sumsq = 0.0;
            rv = updateInterval(fp, newBuffer, sum, sumsq, zoom, tid, start, end, value);
            if(!rv) goto error;
            buffer->next = newBuffer;
            buffer = buffer->next;
            *nEntries += 1;
        }
        start += rv;
    }

    return 0;

error:
    if(newBuffer) {
        if(newBuffer->m) free(newBuffer->p);
        free(newBuffer);
    }
    return 2;
}

//vector fork: append one grid-aligned per-component-SUM zoom record (24+4N) to level k's buffer chain.
static int appendZoomRecVec(bigWigFile_t *fp, uint32_t k, uint32_t N, uint32_t tid,
                            uint64_t binStart, uint64_t binEnd, uint32_t nBases, const double *sumv) {
    size_t zrec = 24 + 4*(size_t)N, c;
    bwZoomBuffer_t *buf = fp->writeBuffer->lastZoomBuffer[k];
    if(buf->l + zrec > buf->m) {
        bwZoomBuffer_t *nb = calloc(1, sizeof(bwZoomBuffer_t));
        size_t cap;
        if(!nb) return 1;
        cap = fp->hdr->bufSize/zrec;
        if(cap < 1) cap = 1;
        nb->p = calloc(cap, zrec);
        if(!nb->p) { free(nb); return 1; }
        nb->m = cap*zrec;
        buf->next = nb;
        fp->writeBuffer->lastZoomBuffer[k] = nb;
        fp->writeBuffer->nNodes[k] += 1;
        buf = nb;
    }
    {
        uint8_t *r = (uint8_t*)buf->p + buf->l;
        *(uint32_t*)(r+0)  = tid;
        *(uint64_t*)(r+4)  = binStart;
        *(uint64_t*)(r+12) = binEnd;
        *(uint32_t*)(r+20) = nBases;
        for(c=0;c<N;c++) *(float*)(r+24+4*(size_t)c) = (float) sumv[c];
        buf->l += zrec;
    }
    return 0;
}

//vector fork: build the per-component-SUM zoom pyramid by re-reading the data (grid-aligned bins).
//Each interval contributes value_c*overlap to each grid bin [b*zoom,(b+1)*zoom) it overlaps; one
//record per non-empty bin per level.  Re-read in windows; an interval is processed once, in the
//window holding its start.  Mirrors constructZoomLevels but with N per-component sums.
int constructZoomLevelsVec(bigWigFile_t *fp) {
    uint32_t N = fp->vecN, k, c, ci;
    uint16_t nL = fp->hdr->nLevels;
    int rv = 1;
    int64_t *curBin = NULL;
    uint64_t *curStart = NULL, *curEnd = NULL;
    uint32_t *curN = NULL, *curTid = NULL;
    double *curSum = NULL;   //nL*N running per-component sums for the open bin
    bwOverlappingIntervalsVec_t *iv = NULL;
    const uint64_t W = 50000000ULL;  //re-read window (bounded memory)

    curBin   = malloc((size_t)nL*sizeof(int64_t));
    curStart = malloc((size_t)nL*sizeof(uint64_t));
    curEnd   = malloc((size_t)nL*sizeof(uint64_t));
    curN     = calloc(nL, sizeof(uint32_t));
    curTid   = calloc(nL, sizeof(uint32_t));
    curSum   = calloc((size_t)nL*N, sizeof(double));
    if(!curBin || !curStart || !curEnd || !curN || !curTid || !curSum) goto done;
    for(k=0;k<nL;k++) curBin[k] = -1;

    for(ci=0; ci<fp->cl->nKeys; ci++) {
        uint64_t len = fp->cl->len[ci], w;
        for(w=0; w<len; w += W) {
            uint64_t we = (w+W < len) ? w+W : len, j;
            iv = bwGetOverlappingIntervalsVec(fp, fp->cl->chrom[ci], w, we);
            if(!iv) goto done;   //vector fork: NULL is error-only here (empty window -> non-NULL n=0); propagate
            for(j=0; j<iv->l; j++) {
                uint64_t s = iv->start[j], e = iv->end[j], st;
                const float *v = iv->value + (size_t)j*N;
                if(s < w || s >= we) continue;   //process each interval once, in its start's window
                for(k=0;k<nL;k++) {
                    uint64_t zoom = fp->hdr->zoomHdrs->level[k];
                    st = s;
                    while(st < e) {
                        int64_t b = (int64_t)(st / zoom);
                        uint64_t binEnd = (uint64_t)(b+1)*zoom, ov;
                        if(binEnd > len) binEnd = len;
                        ov = ((e < binEnd) ? e : binEnd) - st;
                        if(curBin[k] != b || curTid[k] != ci) {
                            if(curBin[k] >= 0) {
                                if(appendZoomRecVec(fp, k, N, curTid[k], curStart[k], curEnd[k], curN[k], curSum + (size_t)k*N)) goto done;
                            }
                            curBin[k] = b; curTid[k] = ci;
                            curStart[k] = (uint64_t)b*zoom;
                            curEnd[k] = binEnd;
                            curN[k] = 0;
                            for(c=0;c<N;c++) curSum[(size_t)k*N+c] = 0.0;
                        }
                        for(c=0;c<N;c++) curSum[(size_t)k*N+c] += (double)ov * (double)v[c];
                        curN[k] += (uint32_t) ov;
                        st += ov;
                    }
                }
            }
            bwDestroyOverlappingIntervalsVec(iv);
            iv = NULL;
        }
    }
    //flush the final open bin per level (the grid-aligned build has no successor to trigger it)
    for(k=0;k<nL;k++) {
        if(curBin[k] >= 0) {
            if(appendZoomRecVec(fp, k, N, curTid[k], curStart[k], curEnd[k], curN[k], curSum + (size_t)k*N)) goto done;
        }
    }
    //make an index for each zoom level (mirrors constructZoomLevels)
    for(k=0;k<nL;k++) {
        fp->hdr->zoomHdrs->idx[k] = calloc(1, sizeof(bwRTree_t));
        if(!fp->hdr->zoomHdrs->idx[k]) goto done;
        fp->hdr->zoomHdrs->idx[k]->blockSize = fp->writeBuffer->blockSize;
    }
    rv = 0;
done:
    if(iv) bwDestroyOverlappingIntervalsVec(iv);
    free(curBin); free(curStart); free(curEnd); free(curN); free(curTid); free(curSum);
    return rv;
}

//Get all of the intervals and add them to the appropriate zoomBuffer
int constructZoomLevels(bigWigFile_t *fp) {
    bwOverlapIterator_t *it = NULL;
    if(fp->vecN) return constructZoomLevelsVec(fp);
    double *sum = NULL, *sumsq = NULL;
    uint32_t i, j, k;

    sum = calloc(fp->hdr->nLevels, sizeof(double));
    sumsq = calloc(fp->hdr->nLevels, sizeof(double));
    if(!sum || !sumsq) goto error;

    for(i=0; i<fp->cl->nKeys; i++) {
        it = bwOverlappingIntervalsIterator(fp, fp->cl->chrom[i], 0, fp->cl->len[i], 100000);
        if(!it) goto error;
	while(it->data != NULL){
	  for(j=0;j<it->intervals->l;j++){
		for(k=0;k<fp->hdr->nLevels;k++){
			if(addIntervalValue(fp, &(fp->writeBuffer->nNodes[k]), sum+k, sumsq+k, fp->writeBuffer->lastZoomBuffer[k], fp->hdr->bufSize/40, fp->hdr->zoomHdrs->level[k], i, it->intervals->start[j], it->intervals->end[j], it->intervals->value[j])) goto error;
			while(fp->writeBuffer->lastZoomBuffer[k]->next) fp->writeBuffer->lastZoomBuffer[k] = fp->writeBuffer->lastZoomBuffer[k]->next;
		}
	  }
	  it = bwIteratorNext(it);
	}
	bwIteratorDestroy(it);

    }

    //Make an index for each zoom level
    for(i=0; i<fp->hdr->nLevels; i++) {
        fp->hdr->zoomHdrs->idx[i] = calloc(1, sizeof(bwRTree_t));
        if(!fp->hdr->zoomHdrs->idx[i]) return 1;
        fp->hdr->zoomHdrs->idx[i]->blockSize = fp->writeBuffer->blockSize;
    }


    free(sum);
    free(sumsq);

    return 0;

error:
    if(it) bwIteratorDestroy(it);
    if(sum) free(sum);
    if(sumsq) free(sumsq);
    return 1;
}

int writeZoomLevels(bigWigFile_t *fp) {
    uint64_t offset1, offset2, idxSize = 0;
    uint32_t i, j, four = 0;  //64-bit fork: root node written field-by-field (no u32 vector)
    uint8_t wrote, one = 0;
    uint16_t actualNLevels = 0;
    int rv;
    bwLL *ll, *p;
    bwRTreeNode_t *root;
    bwZoomBuffer_t *zb, *zb2;
    bwWriteBuffer_t *wb = fp->writeBuffer;
    uLongf sz;
    size_t zrec = fp->vecN ? (24 + 4*(size_t)fp->vecN) : 40;  //vector fork: zoom record size

    for(i=0; i<fp->hdr->nLevels; i++) {
        if(i) {
            //Is this a duplicate level?
            if(fp->writeBuffer->nNodes[i] == fp->writeBuffer->nNodes[i-1]) break;
        }
        actualNLevels++;

        //reserve a uint32_t for the number of blocks
        fp->hdr->zoomHdrs->dataOffset[i] = bwTell(fp);
        fp->writeBuffer->nBlocks = 0;
        fp->writeBuffer->l = 24;
        if(fwrite(&four, sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 1;
        zb = fp->writeBuffer->firstZoomBuffer[i];
        fp->writeBuffer->firstIndexNode = NULL;
        fp->writeBuffer->currentIndexNode = NULL;
        while(zb) {
            sz = fp->hdr->bufSize;
            if(compress(wb->compressP, &sz, zb->p, zb->l) != Z_OK) return 2;

            //write the data to disk
            if(fwrite(wb->compressP, sizeof(uint8_t), sz, fp->URL->x.fp) != sz) return 3;

            //Add an entry into the index
            //64-bit fork: 40B zoom records -- index entry = {first chromId@0, last chromId,
            //first start u64@4, last end u64@+12}. Last record begins at zb->l-40.
            {
                uint8_t *zp = (uint8_t*)zb->p;
                uint8_t *lastr = zp + zb->l - zrec;
                if(addIndexEntry(fp, *(uint32_t*)(zp+0), *(uint32_t*)(lastr+0),
                                     *(uint64_t*)(zp+4), *(uint64_t*)(lastr+12),
                                     bwTell(fp)-sz, sz)) return 4;
            }

            wb->nBlocks++;
            wb->l = 24;
            zb = zb->next;
        }
        if(writeAtPos(&(wb->nBlocks), sizeof(uint32_t), 1, fp->hdr->zoomHdrs->dataOffset[i], fp->URL->x.fp)) return 5;

        //Make the tree
        ll = fp->writeBuffer->firstIndexNode;
        if(ll == fp->writeBuffer->currentIndexNode) {
            root = ll->node;
            idxSize = 4 + 32*root->nChildren;  //64-bit fork: twig child 32B
        } else {
            root = addLeaves(&ll, &idxSize, ceil(((double)fp->writeBuffer->nBlocks)/fp->writeBuffer->blockSize), fp->writeBuffer->blockSize);
        }
        if(!root) return 4;
        fp->hdr->zoomHdrs->idx[i]->root = root;

        ll = fp->writeBuffer->firstIndexNode;
        while(ll) {
            p = ll->next;
            free(ll);
            ll=p;
        }


        //write the index
        wrote = 0;
        fp->hdr->zoomHdrs->indexOffset[i] = bwTell(fp);
        four = IDX_MAGIC;
        if(fwrite(&four, sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 1;
        root = fp->hdr->zoomHdrs->idx[i]->root;
        if(fwrite(&(fp->writeBuffer->blockSize), sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 6;
        if(fwrite(&(fp->writeBuffer->nBlocks), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 7;
        if(fwrite(&(root->chrIdxStart[0]), sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 8;
        if(fwrite(&(root->baseStart[0]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 9;  //64-bit fork (L1)
        if(fwrite(&(root->chrIdxEnd[root->nChildren-1]), sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 10;
        if(fwrite(&(root->baseEnd[root->nChildren-1]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 11;  //64-bit fork (L1)
        if(fwrite(&idxSize, sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 12;
        four = fp->hdr->bufSize/zrec;
        if(fwrite(&four, sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 13;
        four = 0;
        if(fwrite(&four, sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 14; //padding
        fp->hdr->zoomHdrs->idx[i]->rootOffset = bwTell(fp);

        //Write the root node, since writeIndexTree writes the children and fills in the offset
        offset1 = bwTell(fp);
        if(fwrite(&(root->isLeaf), sizeof(uint8_t), 1, fp->URL->x.fp) != 1) return 16;
        if(fwrite(&one, sizeof(uint8_t), 1, fp->URL->x.fp) != 1) return 17; //one byte of padding
        if(fwrite(&(root->nChildren), sizeof(uint16_t), 1, fp->URL->x.fp) != 1) return 18;
        for(j=0; j<root->nChildren; j++) {
            //64-bit fork: per child = chrIdxStart(u32) baseStart(u64) chrIdxEnd(u32) baseEnd(u64) = 24B,
            //then leaf: dataOffset(u64)+size(u64) [40B]; twig: child-offset placeholder(u64) [32B]
            if(fwrite(&(root->chrIdxStart[j]), sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 19;
            if(fwrite(&(root->baseStart[j]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 19;
            if(fwrite(&(root->chrIdxEnd[j]), sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 19;
            if(fwrite(&(root->baseEnd[j]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 19;
            if(root->isLeaf) {
                if(fwrite(&(root->dataOffset[j]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 20;
                if(fwrite(&(root->x.size[j]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 21;
            } else {
                uint64_t zero64 = 0;
                root->dataOffset[j] = 0;
                if(fwrite(&zero64, sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 22;
            }
        }

        while((rv = writeIndexTreeNode(fp->URL->x.fp, fp->hdr->zoomHdrs->idx[i]->root, &wrote, 0)) == 0) {
            if(!wrote) break;
            wrote = 0;
        }

        if(rv || wrote) return 6;

        //Save the file position
        offset2 = bwTell(fp);

        //Write the offsets
        if(writeIndexOffsets(fp->URL->x.fp, root, offset1)) return 2;

        //Move the file pointer back to the end
        bwSetPos(fp, offset2);


        //Free the linked list
        zb = fp->writeBuffer->firstZoomBuffer[i];
        while(zb) {
            if(zb->p) free(zb->p);
            zb2 = zb->next;
            free(zb);
            zb = zb2;
        }
        fp->writeBuffer->firstZoomBuffer[i] = NULL;
    }

    //Free unused zoom levels
    for(i=actualNLevels; i<fp->hdr->nLevels; i++) {
        zb = fp->writeBuffer->firstZoomBuffer[i];
        while(zb) {
            if(zb->p) free(zb->p);
            zb2 = zb->next;
            free(zb);
            zb = zb2;
        }
        fp->writeBuffer->firstZoomBuffer[i] = NULL;
    }

    //Write the zoom headers to disk
    offset1 = bwTell(fp);
    if(bwSetPos(fp, 0x40)) return 7;
    four = 0;
    for(i=0; i<actualNLevels; i++) {
        if(fwrite(&(fp->hdr->zoomHdrs->level[i]), sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 8;
        if(fwrite(&four, sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 9;
        if(fwrite(&(fp->hdr->zoomHdrs->dataOffset[i]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 10;
        if(fwrite(&(fp->hdr->zoomHdrs->indexOffset[i]), sizeof(uint64_t), 1, fp->URL->x.fp) != 1) return 11;
    }

    //Write the number of levels if needed
    if(bwSetPos(fp, 0x6)) return 12;
    if(fwrite(&actualNLevels, sizeof(uint16_t), 1, fp->URL->x.fp) != 1) return 13;

    if(bwSetPos(fp, offset1)) return 14;

    return 0;
}

//0 on success
int bwFinalize(bigWigFile_t *fp) {
    uint32_t four;
    uint64_t offset;
    if(!fp->isWrite) return 0;

    //Flush the buffer
    if(flushBuffer(fp)) return 1; //Valgrind reports a problem here!

    //Update the data section with the number of blocks written
    if(fp->hdr) {
        if(writeAtPos(&(fp->writeBuffer->nBlocks), sizeof(uint64_t), 1, fp->hdr->dataOffset, fp->URL->x.fp)) return 2;
    } else {
        //The header wasn't written!
        return 1;
    }

    //write the bufferSize
    if(fp->hdr->bufSize) {
        if(writeAtPos(&(fp->hdr->bufSize), sizeof(uint32_t), 1, 0x34, fp->URL->x.fp)) return 2;
    }

    //write the summary information
    if(writeSummary(fp)) return 3;

    //Convert the linked-list to a tree and write to disk
    if(writeIndex(fp)) return 4;

    //Zoom level stuff here?
    if(fp->hdr->nLevels && fp->writeBuffer->nBlocks) {
        offset = bwTell(fp);
        if(makeZoomLevels(fp)) return 5;
        if(constructZoomLevels(fp)) return 6;
        bwSetPos(fp, offset);
        if(writeZoomLevels(fp)) return 7; //This write nLevels as well
    }

    //write magic at the end of the file
    four = BIGWIG_MAGIC;
    if(fwrite(&four, sizeof(uint32_t), 1, fp->URL->x.fp) != 1) return 9;

    return 0;
}

/*
data chunk:
uint64_t number of blocks (2 / 110851)
some blocks

an uncompressed data block (24 byte header)
uint32_t Tid	0-4
uint32_t start	4-8
uint32_t end	8-12
uint32_t step	12-16
uint32_t span	16-20
uint8_t type	20
uint8_t padding
uint16_t nItems	22
nItems of:
    type 1: //12 bytes
        uint32_t start
        uint32_t end
        float value
    type 2: //8 bytes
        uint32_t start
        float value
    type 3: //4 bytes
        float value

data block index header
uint32_t magic
uint32_t blockSize (256 in the example) maximum number of children
uint64_t number of blocks (2 / 110851)
uint32_t startTid
uint32_t startPos
uint32_t endTid
uint32_t endPos
uint64_t index size? (0x1E7 / 0x1AF0401F) index address?
uint32_t itemsPerBlock (1 / 1) 1024 for zoom headers 1024 for zoom headers
uint32_t padding

data block index node non-leaf (4 bytes + 24*nChildren)
uint8_t isLeaf
uint8_t padding
uint16_t nChildren (2, 256)
uint32_t startTid
uint32_t startPos
uint32_t endTid
uint32_t endPos
uint64_t dataOffset (0x1AF05853, 0x1AF07057)

data block index node leaf (4 bytes + 32*nChildren)
uint8_t isLeaf
uint8_t padding
uint16_t nChildren (2)
uint32_t startTid
uint32_t startPos
uint32_t endTid
uint32_t endPos
uint64_t dataOffset (0x198, 0x1CF)
uint64_t dataSize (55, 24)

zoom data block
uint32_t number of blocks (10519766)
some data blocks
*/
