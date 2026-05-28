#ifndef TAF_TUI_H_
#define TAF_TUI_H_

/*
 * tui: Universal Column Index for a `cactus-hal2maf --universal` MAF/TAF.
 *
 * The universal column id is the k-th alignment column in file order (blocks
 * are in canonical order and maxRefGap==0 so every block has exactly
 * column_number columns; the count is a stable global integer in [0,T)).
 *
 * This builds, in one streaming scan over all rows, a per-(genome,sequence)
 * piecewise-linear map  genome.seq:pos -> universal column.  Ancestral
 * (row-0/reference) genomes are recorded like any other genome.
 *
 * Output is a single provenanced ONEcode container (schema `P 3 tui`),
 * written to <maf>.tui, holding every genome's run lists plus a directory,
 * the total column count T, and a sparse universal-column -> file-offset
 * anchor index for random access into the underlying MAF/TAF stream.
 *
 * A "run" maps a contiguous, gap-free, colinear stretch of a sequence to a
 * contiguous column range.  Stored fields (all 0-based, half-open in the
 * sequence's FORWARD coordinates):
 *
 *   t_start  forward start of the stretch in the sequence
 *   g_start  universal column of the forward base at t_start
 *   length   number of bases in the stretch
 *   strand   '+' or '-'
 *
 * Query  seq:p  (p in forward coords), for the run with t_start <= p < t_start+length:
 *   '+' : col = g_start + (p - t_start)
 *   '-' : col = g_start + (t_start + length - 1 - p)
 * p in no run  ->  unmapped (leaf/lineage-specific insertion; by design).
 *
 * Consecutive colinear runs of a sequence are merged (one linear map); the
 * remaining runs are stored per sequence as a structure-of-arrays blob: three
 * concatenated LEB128-varint streams gap|gsk|lenc (gap = forward gap to the
 * prev run end ~always 0; gsk = intervening universal columns; lenc =
 * len<<1|strand), zlib-deflated.  Absolute (t,g,len) triples defeat ONElib's
 * Huffman; this is ~5x smaller end to end.  The R line is (inflatedLen INT,
 * deflated bytes STRING).  encode_runs/decode_runs in tui.c are exact
 * inverses and the builder self-checks every sequence's round-trip (asserts
 * on in taffy).
 *
 *  Released under the MIT license, see LICENSE.txt
 */

#include "line_iterator.h"
#include "sonLib.h"

/* Return maf_path + ".tui" (caller frees). */
char *tui_path(const char *maf_path);

/*
 * Build the universal column index for the TAF/MAF behind `li` and write the
 * ONEcode container to `out_path`.  Reads the header and dispatches MAF/TAF
 * internally (mirrors tai_create).  Returns 0 on success.
 *
 * Phase 1: one streaming scan, runs spilled per-genome (column order, O(1) RAM).
 * Phase 2: per-genome in-RAM sort by (sequence, t_start) and ONEcode write.
 *
 * Per-genome spill files are written under `tmp_dir`.  If `tmp_dir` is NULL or
 * empty, the directory of `out_path` is used (roomy, same filesystem as the
 * output, and avoids a possibly tmpfs/small /tmp at vertebrate scale).
 *
 * Genome name resolution: a name with <=1 '.' is split at the (only) '.' (or
 * is itself the genome if no '.').  A name with MORE than one '.' is ambiguous
 * (e.g. assembly accessions `GCF_947179515.1.NC_067472.1`) and needs the set
 * of known genome names.  Precedence: an explicit `genome_name_map` (stHash,
 * keys = names, dummy values; a membership set, NOT a renaming map) wins; else
 * the genome set is derived from the `# hal` Newick tree comment in the header
 * if present (hal2maf emits one, preserved through maf<->taf).  Resolution
 * reuses taf.c's extract_genome_name(); an unresolvable >1-dot name with no
 * source for the genome set is a fatal error.  Pass NULL for genome_name_map
 * to use the header tree (or when the input has no >1-dot names).
 */
/*
 * `n_threads` controls the phase-2 worker pool (per-genome sort + encode);
 * phase 1 stays serial since it's the streaming MAF/TAF scan.  Values <= 1
 * disable OpenMP and execute the loop sequentially.
 */
int tui_create(LI *li, const char *out_path, const char *tmp_dir,
                stHash *genome_name_map, int n_threads);

/////////////////////////////////////////////////////////////////////////////
// Reader / query side (genome.seq:pos -> universal columns).
/////////////////////////////////////////////////////////////////////////////

/* A half-open universal-column interval [start, end). */
typedef struct { int64_t start; int64_t end; } TuiInterval;

typedef struct _Tui Tui;

/*
 * Open a .tui and read its directory (sequence names, lengths, T) into RAM.
 * Returns NULL if the file can't be opened / isn't a tui container.
 */
/*
 * Thread-safety: a `Tui` is immutable after `tui_load` returns (it holds only
 * the .tui path and a small in-memory X index).  Each query call re-opens the
 * file privately, so multiple threads can concurrently issue `tui_query`,
 * `tui_load_seq_runs`, and `tui_genome_lift_load` on the SAME `Tui`.  The
 * `TuiGenomeLift` returned by `tui_genome_lift_load`, however, is NOT
 * thread-safe -- see that function's docs.
 */
Tui *tui_load(const char *tui_path);

void tui_destruct(Tui *tui);

/* Total number of universal columns T (the global column count). */
int64_t tui_total_columns(const Tui *tui);

/*
 * Map the forward-coordinate interval [start, end) of sequence `seq_name` to
 * universal columns.  Loads that one sequence's runs (via oneGoto on the
 * directory-resolved S-ordinal), clips, and returns a malloc'd array (caller
 * frees) of sorted, merged half-open column intervals; sets *n_out to the
 * count.  Returns NULL with *n_out==0 if the sequence is absent or nothing
 * overlaps (e.g. a lineage-specific insertion, by design).
 */
TuiInterval *tui_query(Tui *tui, const char *seq_name,
                       int64_t start, int64_t end, int64_t *n_out);

/*
 * Decode and return ALL of `seq_name`'s runs as a flat array of
 * (t_start, g_start, lenc) triples (where lenc = len << 1 | strand, same
 * encoding the on-disk format uses).  Returned array has 3 * (*n_out)
 * int64 elements; caller frees.  *n_out == 0 and NULL return iff the
 * sequence is absent from the index.
 *
 * Intended for batch callers that need to do many lookups against the
 * SAME sequence's runs: load once, then binary-search the array directly
 * by t_start.  Each call opens + closes the .tui internally, so this is
 * exactly as expensive as one tui_query call -- but avoids the per-lookup
 * open/close that makes tui_query unsuitable in a tight loop.  (See
 * taf_lift.c for the cache-by-source-seq usage pattern.)
 */
int64_t *tui_load_seq_runs(Tui *tui, const char *seq_name, int64_t *n_out);

/////////////////////////////////////////////////////////////////////////////
// Reverse lookup: universal column -> a target genome's coordinate.
//
// The on-disk lift table is sorted by (sequence, t_start), i.e. genome
// coord -- the primary `taffy view` use case is "genome.seq:p -> columns".
// For column-keyed lookups (e.g. annotation liftover from row-0 to a leaf)
// we need a column-sorted view.  The .tui groups each sequence's runs into
// column-window chunks (writer's TUI_CHUNK_RUNS); load time pulls only each
// chunk's (g_min, g_max) plus a pointer to its R blob.  The R blob's runs
// are decoded LAZILY on the first column query that hits that chunk -- so
// narrow queries (annotation lift to one chromosome) pay only a handful of
// decodes, not the full per-genome decode cost.
//
// Memory: once decoded, a chunk's runs stay resident for the lifetime of
// the TuiGenomeLift (no eviction; the lift access pattern is a sequential
// sweep through universal columns, where LRU would thrash).  Worst case
// (full-genome lift touches every chunk) holds the entire target genome's
// lift table in memory, same memory ceiling as a single-shot eager load
// would have used, just amortized as queries fire.
//
// Thread-safety: the TuiGenomeLift holds an open OneFile cursor for lazy R
// decoding AND lazily mutates chunk state on first hit.  All calls on the
// same TuiGenomeLift must be serialized.  (The parent Tui is a separate
// object and can be shared across threads as long as each thread builds
// its own TuiGenomeLift via tui_genome_lift_load.)
/////////////////////////////////////////////////////////////////////////////

typedef struct _TuiGenomeLift TuiGenomeLift;

/*
 * Open the per-genome lift table for "<genome_name>.*" sequences.  Reads
 * each sequence's chunk metadata (small) but does NOT decode the per-chunk
 * R blobs -- those are loaded on first query hit (see header comment).
 * Returns NULL on I/O error or if no d-line matches the prefix.
 */
TuiGenomeLift *tui_genome_lift_load(Tui *tui, const char *genome_name);

void tui_genome_lift_destruct(TuiGenomeLift *gl);

/*
 * One (sequence, position, strand) match returned by tui_genome_lift_column.
 * `seq` is BORROWED from the gl -- valid until tui_genome_lift_destruct.
 */
typedef struct {
    const char *seq;
    int64_t     pos;
    bool        strand;
} TuiGenomeMatch;

/*
 * Look up ALL of the target genome's bases mapped to universal column
 * `column`.  The one-to-many direction: every base of a genome maps to at
 * most one column (cactus-hal2maf --universal --noRefDupes invariant), but
 * one column can carry several bases of the same genome -- paralogs aligned
 * to the same ancestral position.
 *
 * Fills the caller-supplied `out[]` buffer with up to `cap` matches (in
 * unspecified order) and returns the ACTUAL match count.  If the return
 * value exceeds `cap`, only the first `cap` entries are filled; the caller
 * may retry with a larger buffer to get the remainder.  Pass cap=0 (out can
 * be NULL) to just count.
 */
int tui_genome_lift_column(const TuiGenomeLift *gl, int64_t column,
                           TuiGenomeMatch *out, int cap);

/* Number of column chunks indexed for this genome (set at load time;
 * does NOT grow with lazy decodes).  Use for the startup diagnostic. */
int64_t tui_genome_lift_n_chunks(const TuiGenomeLift *gl);

/////////////////////////////////////////////////////////////////////////////
// Universal-column block extractor (replaces the .tai for the -U path).
//
// Drives the SHARED random-access primitives (LI_seek + tai_resync_taf_line +
// the maf/taf block readers, from tai.h) off the .tui's universal-column
// index.  The universal column is a single globally-monotone key (file order
// == column order), so the per-contig non-monotonicity that breaks
// tai_iterator on a universal MAF cannot occur here.  Emits WHOLE blocks (all
// rows) whose universal-column span overlaps any of `iv`, in column order,
// each exactly once.  `iv` must be sorted+merged (tui_query output).
/////////////////////////////////////////////////////////////////////////////

typedef struct _TuiExtractIt TuiExtractIt;

TuiExtractIt *tui_extract_iterator(Tui *tui, LI *li, int is_maf, bool rle,
                                   const TuiInterval *iv, int64_t n_iv);
Alignment *tui_extract_next(TuiExtractIt *it, LI *li);
bool tui_extract_has_next(TuiExtractIt *it);
void tui_extract_iterator_destruct(TuiExtractIt *it);

/* Universal column of the FIRST column of the sub-block just returned by
 * tui_extract_next() (i.e. tcol of its row-0).  Call right after _next(). */
int64_t tui_extract_col_start(const TuiExtractIt *it);

#endif /* TAF_TUI_H_ */

// Local Variables:
// mode: c
// End:
