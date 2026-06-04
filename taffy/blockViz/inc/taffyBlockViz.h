/*
 * taffyBlockViz: C API mirroring HAL's halBlockViz for use by genome
 * browsers (or other consumers) that need block-level alignment data
 * out of a universal .tui index.  All structs / functions follow the
 * hal_X / halX shapes so a thin compat shim (typedefs + #defines)
 * gets existing hal blockViz callers running against taffy.
 *
 * SCOPE -- THIS IS THE INITIAL BASIC-INTEGRATION SUBSET:
 *   Supported:
 *     - taffyOpen / taffyClose / taffyCloseGenome
 *     - taffyGetSpecies                       (list genomes)
 *     - taffyGetChroms                        (list chroms in a genome)
 *     - taffyGetBlocksInTargetRange[_filterByChrom]
 *     - taffyGetDna                           (STUB: returns "NNN..."
 *                                              of requested length)
 *     - taffyGetGenomeMetadata                (returns NULL today)
 *     - taffyFree* for every returned list
 *
 *   NOT supported (will set errStr and return NULL / -1):
 *     - tReversed != 0
 *     - mapBackAdjacencies != 0
 *     - coalescenceLimitName != NULL          (taffy has no analog)
 *     - dupMode != TAFFY_QUERY_AND_TARGET_DUPS (only the default
 *                                              "include all paralogs"
 *                                              shape is wired today)
 *     - LOD-related calls                     (no LOD concept in .tui)
 *     - MAF emit                              (not in initial cut)
 *
 * `errStr`:  if non-NULL, *errStr is set to a malloc'd error message
 * (caller frees) on failure.  If NULL, the call still fails the same
 * way but the message is dropped on the floor (matching halBlockViz).
 *
 * COORDINATES: BED-like, 0-based, half-open, forward-strand-relative
 * for both qStart and tStart.  `strand` in taffy_block_t encodes the
 * alignment strand of qChrom relative to the target reference.
 */

#ifndef _TAFFY_BLOCK_VIZ_H
#define _TAFFY_BLOCK_VIZ_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* keep integer type definition in one place; matches HAL exactly. */
typedef long taffy_int_t;

/** range of coordinates in the target */
struct taffy_target_range_t {
    struct taffy_target_range_t *next;
    taffy_int_t tStart;
    taffy_int_t size;
};

/** paralogous ranges in the *target* genome (e.g. self-alignment of
 * the target that landed in the same column).  Browsers render these
 * as a separate track ("blue lines") so they're split out from the
 * main mappedBlocks list. */
struct taffy_target_dupe_list_t {
    struct taffy_target_dupe_list_t *next;
    taffy_int_t id;
    struct taffy_target_range_t *tRange;
    char *qChrom;
};

/** Mapped + target-dupe block results from a single
 * taffyGetBlocksInTargetRange() call. */
struct taffy_block_results_t {
    struct taffy_block_t *mappedBlocks;
    struct taffy_target_dupe_list_t *targetDupeBlocks;
};

/** One contiguous aligned block (one gap-free run of qChrom that
 * lifts to a contiguous run on the target reference).
 * Coordinates are forward-strand-relative.
 */
struct taffy_block_t {
    struct taffy_block_t *next;
    char *qChrom;
    taffy_int_t tStart;
    taffy_int_t qStart;
    taffy_int_t size;
    char strand;                /* '+' or '-' */
    char *qSequence;            /* query DNA, NULL unless seqMode != NO_SEQUENCES */
    char *tSequence;            /* target DNA, NULL unless seqMode != NO_SEQUENCES */
};

/** Genome (species) metadata. */
struct taffy_species_t {
    struct taffy_species_t *next;
    char *name;
    taffy_int_t length;         /* sum of chrom lengths */
    taffy_int_t numChroms;
    char *parentName;           /* "" (no tree-walk needed for browser blocks) */
    double parentBranchLength;  /* 0.0 */
};

/** Chromosome / sequence within a genome. */
struct taffy_chromosome_t {
    struct taffy_chromosome_t *next;
    char *name;
    taffy_int_t length;
};

/** Optional key/value metadata; today taffyGetGenomeMetadata returns NULL. */
struct taffy_metadata_t {
    struct taffy_metadata_t *next;
    char *key;
    char *value;
};

/** seqMode -- whether to fetch sequence into qSequence / tSequence on
 * each block.  Same values as HAL's hal_seqmode_type_t.
 * Initial cut: only HAL_NO_SEQUENCES is honoured; the others set errStr. */
typedef enum {
    TAFFY_NO_SEQUENCES      = 0,
    TAFFY_HARD_MASKED_SEQ   = 1,
    TAFFY_SOFT_MASKED_SEQ   = 2,
    TAFFY_LOD0_SEQ          = 3
} taffy_seqmode_type_t;

/** dupMode -- which paralog edges to include.  Initial cut: only
 * TAFFY_QUERY_AND_TARGET_DUPS is honoured. */
typedef enum {
    TAFFY_NO_DUPS               = 0,
    TAFFY_QUERY_DUPS            = 1,
    TAFFY_QUERY_AND_TARGET_DUPS = 2
} taffy_dup_type_t;

/* ------------------------------------------------------------------ */
/* Lifecycle                                                           */
/* ------------------------------------------------------------------ */

/** Open a universal .tui file; returns an opaque non-negative handle
 * on success or -1 on failure (with errStr set if non-NULL).
 *
 * @param tuiPath path to the .tui sidecar OR to its companion .maf.gz
 *  (we accept either and resolve the .tui via the same convention as
 *  `taffy lift -i`).
 */
int taffyOpen(const char *tuiPath, char **errStr);

/** Close a handle returned by taffyOpen().  Returns 0 on success, -1
 * on failure. */
int taffyClose(int taffyHandle, char **errStr);

/** Release any cached per-genome lift table held by `genomeName` on
 * the open handle.  Useful to drop memory after a browser pan moves
 * to a new species.  Returns 0 on success / -1 on failure. */
int taffyCloseGenome(int taffyHandle, const char *genomeName, char **errStr);

/* ------------------------------------------------------------------ */
/* Free functions for the returned linked lists.                       */
/* ------------------------------------------------------------------ */

void taffyFreeBlockResults(struct taffy_block_results_t *results);
void taffyFreeBlocks(struct taffy_block_t *block);
void taffyFreeTargetDupeLists(struct taffy_target_dupe_list_t *dupes);
void taffyFreeSpeciesList(struct taffy_species_t *species);
void taffyFreeChromList(struct taffy_chromosome_t *chromosome);
void taffyFreeMetadataList(struct taffy_metadata_t *metadata);

/* ------------------------------------------------------------------ */
/* Block query (the hot path).                                         */
/* ------------------------------------------------------------------ */

/** Return blocks of `qSpecies` aligned within `tChrom`:[tStart, tEnd)
 * on `tSpecies`.  Output is a linked list of taffy_block_t.
 *
 * @param taffyHandle      handle from taffyOpen
 * @param qSpecies         genome to fetch blocks FROM
 * @param tSpecies         genome whose coordinates [tStart, tEnd) frame
 *                         the query range (any genome the .tui knows --
 *                         leaf or ancestor)
 * @param tChrom           chromosome name in tSpecies (bare, no genome prefix)
 * @param tStart           half-open range start
 * @param tEnd             half-open range end (0 = chrom length)
 * @param tReversed        UNSUPPORTED -- must be 0
 * @param seqMode          UNSUPPORTED for !=TAFFY_NO_SEQUENCES (stub today)
 * @param dupMode          UNSUPPORTED for !=TAFFY_QUERY_AND_TARGET_DUPS
 * @param mapBackAdjacencies UNSUPPORTED -- must be 0
 * @param coalescenceLimitName UNSUPPORTED -- must be NULL
 * @param errStr           error message destination (caller-frees) or NULL
 *
 * @return malloc'd taffy_block_results_t (free with taffyFreeBlockResults),
 *         or NULL on error.
 */
struct taffy_block_results_t *taffyGetBlocksInTargetRange(
    int taffyHandle,
    const char *qSpecies, const char *tSpecies,
    const char *tChrom, taffy_int_t tStart, taffy_int_t tEnd,
    taffy_int_t tReversed,
    taffy_seqmode_type_t seqMode,
    taffy_dup_type_t dupMode,
    int mapBackAdjacencies,
    const char *coalescenceLimitName,
    char **errStr);

/** Same as taffyGetBlocksInTargetRange but only emits blocks whose
 * qChrom matches `qChrom` (still on `qSpecies`).  qChrom == NULL is
 * equivalent to the unfiltered call. */
struct taffy_block_results_t *taffyGetBlocksInTargetRange_filterByChrom(
    int taffyHandle,
    const char *qSpecies, const char *tSpecies,
    const char *tChrom, taffy_int_t tStart, taffy_int_t tEnd,
    taffy_int_t tReversed,
    taffy_seqmode_type_t seqMode,
    taffy_dup_type_t dupMode,
    int mapBackAdjacencies,
    const char *qChrom,
    const char *coalescenceLimitName,
    char **errStr);

/* ------------------------------------------------------------------ */
/* Introspection                                                       */
/* ------------------------------------------------------------------ */

/** List of all genomes in the .tui (sorted alphabetically).
 * `parentName` and `parentBranchLength` are stubbed to "" / 0.0 in
 * this initial cut.
 *
 * KNOWN LIMITATION: the implementation infers genome boundaries by
 * splitting d-line keys on the FIRST dot, which misclassifies genome
 * names that contain dots themselves (e.g. NCBI versioned accessions
 * "GCA_028858775.2" come out as "GCA_028858775").  The actual query
 * functions (taffyGetBlocksInTargetRange, taffyGetChroms) take the
 * full genome string from the caller and work correctly regardless;
 * this only affects the cosmetic listing. */
struct taffy_species_t *taffyGetSpecies(int taffyHandle, char **errStr);

/** List of chromosomes in the named genome (sorted alphabetically). */
struct taffy_chromosome_t *taffyGetChroms(int taffyHandle, const char *speciesName, char **errStr);

/** Fetch DNA in [start, end) of `chromName` on `speciesName`.
 *
 * STUB IN THIS INITIAL CUT: always returns a malloc'd string of
 * (end - start) 'N' characters and sets no error.  Browser code that
 * pulls sequence to render won't crash; output is just opaque.  Will
 * be replaced with a real fetch once we wire up DNA access. */
char *taffyGetDna(int taffyHandle, const char *speciesName, const char *chromName,
                  taffy_int_t start, taffy_int_t end, char **errStr);

/** Per-genome metadata as a linked list.  Returns NULL today (.tui
 * does not yet carry per-genome metadata). */
struct taffy_metadata_t *taffyGetGenomeMetadata(int taffyHandle, const char *genomeName, char **errStr);

#ifdef __cplusplus
}
#endif

#endif /* _TAFFY_BLOCK_VIZ_H */
