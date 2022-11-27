/*
 * Quick implementation of O(ND) algorithm from Myers, using terminology
 * from S. Marco-Sola et al., Fast gap-affine pairwise alignment using the wavefront algorithm.
 * Is useful for quickly computing the diff between two similar lists.
 */

#include "ond.h"
#include "sonLib.h"

typedef struct _WF {
    /*
     * Represents a "wavefront", a series of points along the x+y diagonal that represent "furthest points".
    */
    int64_t min_diag, max_diag; // Min and max diag are the bounds (inclusive) on the diagonal
    int64_t original_min_diag; // We keep track of the first min_diag value, because this is used to reference points
    // in the fpa array even if we subsequently trim the min and max diags
    int64_t *fpa; // the array of furthest points
} WF;

WF *WF_construct(int64_t min_diag, int64_t max_diag) {
    WF *wf = st_calloc(1, sizeof(WF));
    wf->min_diag = min_diag;
    wf->max_diag = max_diag;
    wf->original_min_diag = min_diag;
    assert(max_diag >= min_diag);
    wf->fpa = st_calloc(1+max_diag-min_diag, sizeof(int64_t));
    return wf;
}

void WF_destruct(WF *wf) {
    if(wf != NULL) {
        free(wf->fpa);
        free(wf);
    }
}

int64_t WF_get_fp(WF *wf, int64_t k) {
    /*
     * Returns the further point (an x coordinate) on the x - y = k antidiagonal
    */
    if(k < wf->min_diag || k > wf->max_diag) {
        return -1000000;  //if the point is not on the anti-diagonal then return a very small number, indicating
        // it is unreachable
    }
    return wf->fpa[k - wf->original_min_diag];
}

void WF_set_fp(WF *wf, int64_t k, int64_t h) {
    /*
     * Set the further point (an x coordinate) on the x - y = k antidiagonal
    */
    assert(wf->min_diag <= k);
    assert(k <= wf->max_diag);  // Otherwise we're trying to set a point not on the wavefront
    wf->fpa[k - wf->original_min_diag] = h;
}

void WF_trim_diagonal(WF *wf, int64_t new_min_diag, int64_t new_max_diag) {
    /*
     * Trim the wavefront (e.g. to remove outermost points).
    */
    // We can only trim to make the diagonal smaller or the same size
    assert(wf->min_diag <= new_min_diag);
    assert(new_min_diag <= new_max_diag);
    assert(new_max_diag <= wf->max_diag);

    // No need to change the fpa array itself
    wf->min_diag = new_min_diag;
    wf->max_diag = new_max_diag;
}

typedef struct _WFS {
    /*
     * Represents a wavefront set, i.e. a series of wavefronts, one for each score.
    */
    stList *wfl; // list of wavefronts indexed by score
} WFS;

WFS *WFS_construct() {
    WFS *wfs = st_calloc(1, sizeof(WFS));
    wfs->wfl = stList_construct3(0, (void (*)(void *))WF_destruct);
    stList_append(wfs->wfl, WF_construct(0, 0));
    return wfs;
}

void WFS_destruct(WFS *wfs) {
    stList_destruct(wfs->wfl);
    free(wfs);
}

WF *WFS_get_wf(WFS *wfs, int64_t s) {
    /*
    * Get the wavefront for score s
    */
    return (s >= 0 && s < stList_length(wfs->wfl)) ? stList_get(wfs->wfl, s) : NULL;
}

int64_t WFS_get_fp(WFS *wfs, int64_t s, int64_t k) {
    /*
    * Get the furthest point for a score s and antidiaonal k = x - y
    */
    WF *wf = WFS_get_wf(wfs, s);
    if (wf == NULL) {
        return -100000;  // If the furthest point is not defined return a very small number
    }
    return WF_get_fp(wf, k);
}

void WFS_set_fp(WFS *wfs, int64_t s, int64_t k, int64_t h) {
    /*
    * Set the furthest point for a score s and antidiaonal k = x - y
    */
    WF_set_fp(stList_get(wfs->wfl, s), k, h);
}

WF *WFS_add_wf(WFS *wfs, int64_t min_diag, int64_t max_diag, int64_t s) {
    /*
     * Adds a wavefront to the set.
    */
    WF *wf = WF_construct(min_diag, max_diag);
    assert(s >= stList_length(wfs->wfl));
    while(s > stList_length(wfs->wfl)) { // pad out any intermediate points
        stList_append(wfs->wfl, NULL);
    }
    stList_append(wfs->wfl, wf);
    assert(stList_get(wfs->wfl, s) == wf);
    return wf;
}

int64_t WFS_get_min_diag(WFS *wfs, int64_t s) {
    /*
     * Get the minimum k=x-y for the wavefront for score s
    */
    WF *wf = WFS_get_wf(wfs, s);
    if (wf == NULL) {
        return 1000000000;  // If the wfs is not defined return a very large number
    }
    return wf->min_diag;
}

int64_t WFS_get_max_diag(WFS *wfs, int64_t s) {
    /*
    * Get the maximum k=x-y for the wavefront for score s
    */
    WF *wf = WFS_get_wf(wfs, s);
    if (wf == NULL) {
        return -1000000000;  // If the wfs is not defined return a very small number
    }
    return wf->max_diag;
}

struct _WFA {
    void *string1;
    void *string2;
    int64_t string1_length;
    int64_t string2_length;
    size_t element_size; // the size in bytes of each element in string1 and string2
    int64_t gap_score, mismatch_score;
    bool (*elements_equal)(void *, void *);
    int64_t s; // The starting alignment score
    WFS *wfs; // The wavefront set
};

void WFA_destruct(WFA *wfa) {
    WFS_destruct(wfa->wfs);
    free(wfa);
}

void *get_element(void *string, size_t element_size, int64_t i) {
    return &(((char *)string)[i * element_size]);
}

void WFA_extend(WFA *wfa) {
    /*
    * Extends each point on the current wavefront by alignment matches.
    */
    // Get the current wavefront, whose points are to be extended
    WF *wf = WFS_get_wf(wfa->wfs, wfa->s);
    assert(wf != NULL);
    // For each diagonal on the wf extend it by the maximum number of matches from the current furthest point
    for(int64_t k=wf->min_diag; k<=wf->max_diag; k++) {
        int64_t h = WF_get_fp(wf, k);
        if(h >= 0 && h - k >= 0) {  // If h = x-y such that x >= 0 and y >= 0
            while(h < wfa->string1_length && h - k < wfa->string2_length &&
            wfa->elements_equal(get_element(wfa->string1, wfa->element_size, h),
                                get_element(wfa->string2, wfa->element_size, h - k))) {
                // Extend the furthest point
                h += 1;
                WF_set_fp(wf, k, h);
            }
        }
    }
}

bool WFA_done(WFA *wfa) {
    /*
     * Are we at the end of the dp matrix?
    */
    return WFS_get_fp(wfa->wfs, wfa->s, wfa->string1_length - wfa->string2_length) == wfa->string1_length;
}

static int64_t max(int64_t i, int64_t j) {
    return i > j ? i : j;
}

static int64_t min(int64_t i, int64_t j) {
    return i < j ? i : j;
}

void WFA_next(WFA *wfa) {
    /*
     * Adds the next score wavefront to the set.
    */
    while(1) { // Get the next score by increasing s until we find s minus mismatch or gap score has a
        // wavefront
        wfa->s++; // Increment s
        if (WFS_get_wf(wfa->wfs, wfa->s - wfa->gap_score) != NULL ||
            WFS_get_wf(wfa->wfs, wfa->s - wfa->mismatch_score) != NULL) {
            break;  // There is a prior wavefront to connect to
        }
    }

    // Update min and max diag
    int64_t min_diag = min(WFS_get_min_diag(wfa->wfs, wfa->s - wfa->gap_score),
                           WFS_get_min_diag(wfa->wfs, wfa->s - wfa->mismatch_score)) - 1;
    int64_t max_diag = max(WFS_get_max_diag(wfa->wfs, wfa->s - wfa->gap_score),
                           WFS_get_max_diag(wfa->wfs, wfa->s - wfa->mismatch_score)) + 1;

    // Add the next WFS line
    WF *wf = WFS_add_wf(wfa->wfs, min_diag, max_diag, wfa->s);

    // Do dp calcs
    for(int64_t k=wf->min_diag; k<=wf->max_diag; k++) {
        WF_set_fp(wf, k, max(max(WFS_get_fp(wfa->wfs, wfa->s - wfa->gap_score, k - 1) + 1,  // insert in string1
                                 WFS_get_fp(wfa->wfs, wfa->s - wfa->gap_score, k + 1)),  // insert in string2
                             WFS_get_fp(wfa->wfs, wfa->s - wfa->mismatch_score, k) + 1));  // mismatch
    }
}

WFA *WFA_construct(void *string1, void *string2, int64_t string1_length, int64_t string2_length,
                   size_t element_size, bool (*elements_equal)(void *, void *),
                   int64_t gap_score, int64_t mismatch_score) {
    /* Finds an optimal global alignment of two strings using WFS algorithm.
    * The algorithm is as described in https://doi.org/10.1093/bioinformatics/btaa777
    * The notation/language somewhat follows the paper, but is otherwise as follows:
    * The two input strings string1 (x) and string2 (y)
    * In the dp matrix we have (x, y) row,column coordinates, thus string1 is along the rows and string2
    *         is along the columns.
    * The anti-diagonal is k = x-y
    * The diagonal is x+y
    *         The coordinates can be visualized as follows:
    * y-1 y+0 y+1
    * x-1 k+0 k-1 k-2
    * x+0 k+1 k+0 k-1
    * x+1 k+2 k+1 k+0
    * As in the paper, the further points, "fp", are represented as x coordinates along the anti-diagonal.
    */
    WFA *wfa = st_calloc(1, sizeof(WFA));
    wfa->string1 = string1;
    wfa->string2 = string2;
    wfa->string1_length = string1_length;
    wfa->string2_length = string2_length;
    wfa->element_size = element_size;
    wfa->gap_score = gap_score;
    wfa->mismatch_score = mismatch_score;
    wfa->elements_equal = elements_equal;
    wfa->wfs = WFS_construct();  // The wavefront set
    wfa->s = 0;  // The starting alignment score

    // Run the wavefront dynamic programming process to find the optimal alignment
    while(1) {
        WFA_extend(wfa);  // Extend the wavefront
        if (WFA_done(wfa)) {  // We're done if we reach the end of the dp matrix
            break;
        }
        WFA_next(wfa);  // Set up the next wavefront
    }
    return wfa;
}

int64_t WFA_get_alignment_score(WFA *wfa) {
    /*
     * Return the alignment score
    */
    return wfa->s;
}

void WFA_get_alignment(WFA *wfa, int64_t *elements_aligned_to_string1) {
    /*
    * Returns an alignment of the two string. Implements the traceback algorithm. The alignment is represented
     * as an array of integers of length equal to the length of string1, the entries represent the indices of positions
     * in string2, with -1 if the position in string1 is aligned to a gap.
    */
    int64_t t = wfa->s;  // The score of the sub-alignment that we're tracing back
    int64_t k = wfa->string1_length - wfa->string2_length;  // The diagonal we're tracing back on
    int64_t f = wfa->string1_length;  // The furthest point
    //stList *alignment = stList_construct();  // The alignment, represented as a sequence of (x, y) pairs
    for(int64_t i=0; i<wfa->string1_length; i++) { // Initialize the alignment positions to be all gaps
        elements_aligned_to_string1[i] = -1;
    }
    assert(WFS_get_fp(wfa->wfs, t, k) == f);  // This is the condition that must be true at the beginning of trace back
    while (k != 0 || f != 0) {  // While we haven't gotten to the first cell in the dp matrix
        // Do backtrace dp calcs
        int64_t a = WFS_get_fp(wfa->wfs, t - wfa->mismatch_score, k);  // match
        int64_t b = WFS_get_fp(wfa->wfs, t - wfa->gap_score, k - 1);  // insert in string1 (x)
        int64_t c = WFS_get_fp(wfa->wfs, t - wfa->gap_score, k + 1);  // insert in string2 (y)

        while (f > max(max(a, b + 1), max(c, 0))) {  // The plus one for an insert in string1 is necessary
            // k = x - y, f = x
            int64_t x = f;
            int64_t y = -(k - f);
            elements_aligned_to_string1[x-1] = y-1; // subtract one to get sequence coordinates
            f -= 1;
        }

        if (a >= b && a >= c) {  // we must take a mis-match
            t -= wfa->mismatch_score;
        } else if (b >= c) {  // alignment has insert in string1
            assert(b >= a);
            k -= 1;
            f -= 1;
            t -= wfa->gap_score;
        } else {  // alignment has insert in string2
            assert(c >= a && c >= b);
            k += 1;
            t -= wfa->gap_score;
        }
    }
}
