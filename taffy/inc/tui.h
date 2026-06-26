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

/* A half-open universal-column interval [start, end), carrying enough
 * info for callers to recover the source genome's position and strand
 * orientation at any column inside it without a separate lookup.
 *
 *   start, end  half-open universal column range
 *   t_start     source genome position at column == start.  For rev=0
 *               this is the LOW source coord (start of the queried sub-
 *               range that this iv covers).  For rev=1 it is the
 *               highest source coord, because column-ascending traverses
 *               source-descending under reverse mapping.
 *   rev         0 = source-to-column mapping is forward (column c maps
 *               to source pos t_start + (c - start)).
 *               1 = reverse (column c maps to source pos
 *               t_start - (c - start)).
 *
 * The rev bit is the source genome's strand at this universal column
 * range.  Combining it with the strand of a run from another genome's
 * lift (via XOR) yields the actual relative strand between source and
 * target -- the only correct way to label strand for paralog / SD
 * mappings where ancestor blocks are inverted relative to both genomes.
 */
typedef struct {
    int64_t start;
    int64_t end;
    int64_t t_start;
    int     rev;
} TuiInterval;

typedef struct _Tui Tui;

/*
 * Open a .tui and read its directory (sequence names, lengths, T) into RAM.
 * Returns NULL if the file can't be opened / isn't a tui container.
 */
/*
 * Thread-safety: tui_query, tui_load_seq_runs, tui_sequence_lengths,
 * tui_genome_names, tui_seq_length, tui_genome_seq_lengths and
 * tui_genome_lift_load all share the cached OneFile cursor on tui->of (opened
 * once in tui_load), so they are NOT thread-safe on the same Tui *.  Callers
 * serialize externally (blockViz holds g_mutex; CLI callers are single-
 * threaded).  The `TuiGenomeLift` returned by tui_genome_lift_load is NOT
 * thread-safe -- see that function's docs.
 */
Tui *tui_load(const char *tui_path);

void tui_destruct(Tui *tui);

/* Total number of universal columns T (the global column count). */
int64_t tui_total_columns(const Tui *tui);

/* X-record (universal-column -> source MAF file-position) accessors.
 * Used by sidecar writers (tui-chain) that need to copy the X-track from
 * a source .tui forward into a derived .tui so that view -r still works
 * against the original .maf.gz.  Returned arrays alias internal storage;
 * do not free or modify. */
int64_t        tui_idx_n   (const Tui *tui);
const int64_t *tui_idx_cols(const Tui *tui);
const int64_t *tui_idx_fpos(const Tui *tui);

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

/*
 * Enumerate every sequence in the .tui's directory: name -> length (bp).
 * Walks the d-records (already name-sorted by the builder) via oneGoto and
 * returns an stHash mapping a freshly-allocated char* (the seq name) to
 * the seq length stored directly as the value pointer (cast via intptr_t,
 * mirroring tai_sequence_lengths's convention).
 *
 * Returns NULL on I/O error.  Reuses tui's cached cursor (no re-open);
 * serialize on `tui` (same note as tui_query).
 *
 * Includes ALL sequences in the universal alignment -- both anchor (row-0)
 * and leaf-genome.  The .tui doesn't tag them; callers that need anchor-
 * only (e.g. for sharding without double-counting universal columns) must
 * filter by name prefix against the # hal tree's internal labels, OR shard
 * by universal-column range instead (see tui_total_columns).
 */
stHash *tui_sequence_lengths(Tui *tui);

/*
 * Length (bp) of ONE sequence by its full "genome.sequence" d-line name, via
 * a single O(log n_d) binary search of the name-sorted directory -- the
 * targeted alternative to tui_sequence_lengths when you need just one length
 * (e.g. resolving a query's full-chrom end).  Returns -1 if absent.  Uses the
 * cached cursor on `tui`, so the caller must serialize on it (same note as
 * tui_query / tui_load_seq_runs).
 */
int64_t tui_seq_length(Tui *tui, const char *seq_name);

/*
 * Sequences of ONE genome: full "genome.sequence" name -> length (bp), same
 * stHash value convention as tui_sequence_lengths.  Lower-bound seeks to the
 * genome's "<genome>." prefix then scans its contiguous d-records (the
 * directory is name-sorted, so a genome's sequences are adjacent) -- O(log
 * n_d + k) for k the genome's sequence count, NOT the whole directory.
 * Returns NULL if the genome has no sequences.  Uses the cached cursor on
 * `tui`; serialize on it.
 */
stHash *tui_genome_seq_lengths(Tui *tui, const char *genome);

/*
 * Per-resolved-genome metadata, returned by tui_genome_names().  `name` is
 * malloc'd; free via tui_genome_info_free.
 */
typedef struct {
    char    *name;
    int64_t  total_bp;
    int64_t  n_chroms;
} TuiGenomeInfo;

/*
 * Enumerate the genome roster (g records) embedded in the .tui at index
 * time, in writer order.  Returns NULL with *n_out == 0 when the .tui
 * predates the g-record schema (those callers should fall back to a
 * heuristic on `tui_sequence_lengths` output).  Caller frees via
 * tui_genome_info_free.
 *
 * The roster is deterministic and resolves the "<genome>.<sequence>"
 * boundary without ambiguity for genome names that contain dots
 * themselves (NCBI versioned accessions like "GCA_028858775.2").
 */
TuiGenomeInfo *tui_genome_names(Tui *tui, int64_t *n_out);
void tui_genome_info_free(TuiGenomeInfo *info, int64_t n);

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

/* Smallest universal column this genome maps to (chunks[0].g_min), or -1 if it
 * has no chunks.  Cheap (chunk metadata only, no R decode).  For row-0 ancestors
 * in a universal MAF -- column-contiguous in pre-order -- this is the genome's
 * row-0 backbone start, the correct key for ordering its row-0 column range. */
int64_t tui_genome_lift_min_col(const TuiGenomeLift *gl);

/*
 * One gap-free run of source-to-target alignment, as visited by
 * tui_genome_lift_visit_runs.  `seq` is borrowed from gl (valid until
 * tui_genome_lift_destruct).  Source columns covered by the run are
 * [g_start, g_start + length); target coords are computed by the
 * caller based on `strand` (1 = '+': target advances with source;
 * 0 = '-': target advances inversely, target_pos(c) =
 * t_start + length - 1 - (c - g_start)).  Range-precise target
 * range covered by the full run is [t_start, t_start + length)
 * regardless of strand.
 */
typedef struct {
    const char *seq;
    int64_t     g_start;
    int64_t     length;
    int64_t     t_start;
    int         strand;
} TuiRun;

/*
 * Visit every run in `gl` whose source-column range
 * [g_start, g_start + length) intersects [c_lo, c_hi).  For each
 * such run, calls cb(run, user).  Chunks are lazily decoded on
 * demand (single-threaded).  Order: chunks visited in g_min order;
 * within each chunk, runs in g_start order.
 *
 * Use this for the "lift a range" idiom (whole-chromosome browser
 * views, bulk annotation transfer); it's O(n_runs_in_range) versus
 * the O(n_columns_in_range) cost of repeated tui_genome_lift_column
 * calls.
 */
void tui_genome_lift_visit_runs(TuiGenomeLift *gl, int64_t c_lo, int64_t c_hi,
                                void (*cb)(const TuiRun *run, void *user),
                                void *user);

/*
 * Stream every run of every chunk in this genome's lift table in g_min
 * order.  Unlike tui_genome_lift_visit_runs, decoded chunks are FREED
 * immediately after the visit (no per-chunk caching) -- peak memory is
 * one decoded chunk (~3 MB) rather than the full lift table (tens of GB
 * on giant ancestors).  Chunks that were ALREADY decoded by an earlier
 * caller are left intact -- safe to interleave with visit_runs callers.
 *
 * For one-shot whole-genome scans (tui-chain) where each chunk is
 * touched exactly once and the visit_runs cache would only inflate
 * peak RSS.  The cb callback contract matches visit_runs.
 */
void tui_genome_lift_stream_runs(TuiGenomeLift *gl,
                                 void (*cb)(const TuiRun *run, void *user),
                                 void *user);

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

/*
 * Disown the alignment from the iterator's "auto-free on next call" slot.
 * Use when the caller needs to keep the most recent yield alive across the
 * next tui_extract_next() call -- e.g. to feed it as the p_alignment to
 * taf_write_block2 so consecutive blocks' shared rows get properly
 * delta-encoded.  After this call, the iterator no longer tracks the
 * yield; the caller is responsible for `alignment_destruct`-ing it.
 *
 * Safe to call when there's no current yield (no-op).
 */
void tui_extract_take_ownership(TuiExtractIt *it);

/* Universal column of the FIRST column of the sub-block just returned by
 * tui_extract_next() (i.e. tcol of its row-0).  Call right after _next(). */
int64_t tui_extract_col_start(const TuiExtractIt *it);

/////////////////////////////////////////////////////////////////////////////
// internal: shared writer-side primitives.  Exported so secondary writers
// (currently tui-chain, which builds chained .tui files from existing
// .tui inputs) reuse the same on-disk encoding without copy-paste
// of the encode/cleanup logic.  NOT part of the stable public API; format
// changes here must update BOTH tui_create AND every other writer.
//
// Pulls in ONElib.h for OneFile / OneSchema -- these types ride the writer
// API.  Public callers (taffy view / lift / blockViz) already see ONElib
// transitively; this just makes the inclusion explicit.
/////////////////////////////////////////////////////////////////////////////

#include "ONElib.h"

/* Per-chunk caps used by tui_create's phase-2 chunk-boundary loop.  Match
 * these in coarsened writers so the reader's per-chunk t-range skip + the
 * "TuiRun is at most 65k per chunk" invariant hold.  Comment block at
 * tui.c:81-108 explains the rationale. */
#define TUI_CHUNK_RUNS   65536
#define TUI_CHUNK_G_MAX  1000000   /* universal-column span per chunk */

/* Open a new .tui in WRITE mode, register a crash-cleanup atexit hook,
 * stamp the OneCode provenance record.  Returns the OneFile* on success
 * and writes *schema_out (caller frees via oneSchemaDestroy AFTER
 * oneFileClose).  Returns NULL on open failure; *schema_out is set NULL.
 * Caller is responsible for any pre-open state cleanup + calling
 * tui_atexit_disarm() on the error path.  prog/what/blurb populate the
 * provenance record (e.g. "taffy", "tui-coarsen", "LOD B=1000000"). */
OneFile *tui_open_write(const char *out_path, const char *prog,
                        const char *what, const char *blurb,
                        OneSchema **schema_out);

/* Write the `t` header record: total universal columns + the current on-disk
 * format version (major, minor).  Call once, right after tui_open_write. */
void tui_write_header(OneFile *of, int64_t T);

/* One output chunk for tui_write_sequence: the chunk's universal-column range
 * [g_min, g_max) and source-coord range [t_min, t_max), plus its run blob
 * (raw_len = inflated length; def/def_len = the deflate from tui_encode_runs).
 * The caller owns `def`. */
typedef struct {
    int64_t  g_min, g_max;
    int64_t  t_min, t_max;
    int64_t  raw_len;
    uint8_t *def;
    int64_t  def_len;
} TuiWriteChunk;

/* Write one sequence -- the `S` record (seqName, seqLen, n_chunks) followed by
 * its (C, R)+ chunk pairs (C is a data line, not indexed).  THE single source
 * of truth for the on-disk S/C/R layout: every .tui producer (taffy index,
 * tui-chain) emits through here, so a format change can't be applied to one
 * writer and missed in another.  Does NOT free the chunks' `def` blobs. */
void tui_write_sequence(OneFile *of, const char *seq_name, int64_t seq_len,
                        const TuiWriteChunk *chunks, int64_t n_chunks);

/* Encode `m` (t, g, lenc) triples in `buf` (m*3 int64) as the standard
 * .tui R-record payload: header + three SoA varint streams (gap | gsk |
 * lenc) + zlib deflate.  Caller frees the returned buffer.  *raw_len and
 * *def_len are set to the uncompressed and compressed sizes. */
uint8_t *tui_encode_runs(const int64_t *buf, int64_t m,
                         int64_t *raw_len, int64_t *def_len);

/* Encode `n` (col, fpos) anchor pairs as the standard .tui X-record
 * payload (delta-varint per stream + zlib).  Caller frees. */
uint8_t *tui_encode_idx(const int64_t *col, const int64_t *fpos, int64_t n,
                        int64_t *raw_len, int64_t *def_len);

/* Register/clear the crash-cleanup hook that removes a half-written .tui
 * on SIGINT / st_errAbort.  Idempotent.  Both should bracket the
 * tui_open_write / oneFileClose pair. */
void tui_atexit_track_tui(const char *path);
void tui_atexit_disarm(void);

#endif /* TAF_TUI_H_ */

// Local Variables:
// mode: c
// End:
