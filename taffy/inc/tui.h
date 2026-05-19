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
 * piecewise-linear map  genome.seq:pos -> universal column  (Index B).  The
 * ancestral (row-0/reference) genomes are recorded like any other genome;
 * their runs ARE Index A (column <-> ancestral), reconstructed at query time.
 *
 * Output is a single provenanced ONEcode container (schema `P 3 tui`),
 * written to <maf>.tui, holding every genome's run lists plus a directory,
 * the total column count T, and the set of reference (ancestral) genomes.
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
int tui_create(LI *li, const char *out_path, const char *tmp_dir,
                stHash *genome_name_map);

/////////////////////////////////////////////////////////////////////////////
// Reader / query side (Index B: genome.seq:pos -> universal columns).
/////////////////////////////////////////////////////////////////////////////

/* A half-open universal-column interval [start, end). */
typedef struct { int64_t start; int64_t end; } TuiInterval;

/* A row-0 (ancestor) coordinate piece: [start, start+len) of `seq`. */
typedef struct { char *seq; int64_t start; int64_t len; } TuiRef;

typedef struct _Tui Tui;

/*
 * Open a .tui and read its directory (sequence names, lengths, T) into RAM.
 * Returns NULL if the file can't be opened / isn't a tui container.
 */
Tui *tui_load(const char *tui_path);

void tui_destruct(Tui *tui);

/* Total number of universal columns T (the global column count). */
int64_t tui_total_columns(const Tui *tui);

/* 1 if the full "genome.sequence" name is present in the index, else 0. */
int tui_has_sequence(const Tui *tui, const char *seq_name);

/*
 * Map the forward-coordinate interval [start, end) of sequence `seq_name` to
 * universal columns.  Loads that one sequence's runs, clips, and returns a
 * malloc'd array (caller frees) of sorted, merged half-open column intervals;
 * sets *n_out to the count.  Returns NULL with *n_out==0 if the sequence is
 * absent or nothing overlaps (e.g. a lineage-specific insertion, by design).
 *
 * Baseline implementation: re-scans the container to find the sequence's R
 * blob (no oneGoto yet) -- correctness first; expected to be slow.
 */
TuiInterval *tui_query(Tui *tui, const char *seq_name,
                       int64_t start, int64_t end, int64_t *n_out);

/*
 * Index A: map universal-column intervals to row-0 (ancestor) coordinate
 * pieces via the explicit reference track.  `uiv` is the column-ordered,
 * merged output of tui_query.  Returns a malloc'd array (caller frees the
 * array; the `.seq` pointers are owned by `tui`, do NOT free them) in column
 * order; *n_out = count.  One input interval can yield several pieces on
 * different row-0 sequences/ancestors.  Returns NULL with *n_out==0 if empty.
 */
TuiRef *tui_col_range_to_ref(const Tui *tui, const TuiInterval *uiv,
                             int64_t n_uiv, int64_t *n_out);

#endif /* TAF_TUI_H_ */

// Local Variables:
// mode: c
// End:
