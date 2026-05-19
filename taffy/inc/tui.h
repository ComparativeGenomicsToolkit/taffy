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

#endif /* TAF_TUI_H_ */

// Local Variables:
// mode: c
// End:
