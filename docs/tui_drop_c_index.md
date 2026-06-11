# Design: drop the `.tui` C-index (cheap open) — demote `C` to a data line

**Status:** Proposal / design — 2026-06-11
**Format version:** 0.1 → **0.2** (breaking; regenerate all `.tui`)
**Scope:** option A/C from the open-overhead investigation — eliminate the per-open
cost that blocks remote access and adds 2–3 s to every 577-way browser query.
Remote access (htslib range I/O) is out of scope here but this is its prerequisite.

---

## 1. Problem

The `.tui` universal-MAF index is built on ONElib/OneCode. ONElib gives random
access (`oneGoto(type, i)`) by materializing, **per *object* line type, a full
`I64[count+1]` array of file offsets in RAM on open** — allocated at the `#`
count line (ONElib.c:1546-1550) and filled from a list-codec-**compressed** blob
in the footer at the `&` line (ONElib.c:1625-1631).

The `.tui` object types are `d` (directory), `S` (sequence), `C` (chunk),
`g` (genome roster). **The `C` (chunk) index dominates** because `#chunks ≫ #seqs`:

| `.tui` | file | open wall | open RSS | dominant index |
|---|---|---|---|---|
| apes (15 genomes) | 357 MB | 0.15 s | 23 MB | C |
| 577-way (measured earlier) | ~92 GB | ~2.4 s | ~0.8–1 GB | C |

A single query touches only `O(log n_seqs)` directory entries plus one genome's
handful of chunks, so **~90 % of what open materializes is the `C` array, which
any one query barely reads.** Worse, a browser query opens the `.tui` **3–4
times** (the shared handle, the lift handle, the roster open for `GetSpecies`,
the seq-length open for `GetChroms`) — each paying the full `C`-index cost.

ONElib is stdio-only, so this also makes remote access (download-a-small-index +
range reads, the `.tai` model) impossible.

## 2. Why `C` is an object today, and why it no longer needs to be

`C` was made an ONElib *object* (`O C 4 …`, tui.c:71) for exactly one reason:
to support `oneGoto(of,'C',ord)` for random chunk access. A 13-agent audit
(2026-06-11) confirmed there are **exactly two** `oneGoto('C')` sites in the whole
tree, both inside the lift:

1. **`tui_genome_lift_load` pass-2** (tui.c:2164) — the eager load-time metadata
   scan: per target seq, `oneGoto('S')` then per chunk `oneGoto('C', first_c_ord+k)`
   reading only `g_min/g_span`, **jumping C→C and skipping the paired R blob.**
2. **`chunk_decode`** (tui.c:2232) — lazy decode on a query hit:
   `oneGoto('C', ch->c_ord)`, read C then its paired R, decode the runs.

Audit findings (all with `file:line` evidence in §10):

- There is **no** `oneStats('C')`, no reliance on the `C` object *count*, no
  `oneCountUntilNext`, no `oneObject('C')` anywhere.
- Every other reader of C records — `tui_query` (tui.c:1770), `tui_load_seq_runs`
  (tui.c:1876) — reaches them by **sequential `oneReadLine` after `oneGoto('S')`**,
  consuming only the `S` index, never the `C` index. They keep working verbatim.
- `taf_lift`, `taf_view`, `taf_chain`, blockViz, and the tests touch chunks only
  through the public `tui_genome_lift_*` / `tui_load_seq_runs` API — **zero**
  direct `oneGoto`/`c_ord` use.
- The `X` column-anchor record stores `(col, fpos)` = **byte offsets into the
  source MAF** (tui.c:40, `tui_idx_fpos`) — independent of `C`. Unaffected.
- **Safety invariant (verified, high confidence, survived all refuters):**
  `chunk_decode` is only ever reached via `gl->chunks`, which is built **once**
  in pass-2 and never reallocated. So an offset recorded during pass-2 is
  *always* available at decode time — there is no cold chunk access.

**Conclusion:** the only thing keeping `C` an object is `oneGoto('C')`, and both
call sites can navigate by a byte offset recorded during pass-2 instead. So `C`
can become a plain data line and the index can be dropped from the format.

## 3. Decision: drop, don't skip

Two ways to kill the open cost were considered:

- **(A) Keep the C index on disk, skip *loading* it** (an `oneFileOpenReadSkipIndex`
  flag that suppresses the alloc + the `&` read for `C`). Read-side only;
  preserves compatibility. **Rejected:** it keeps a now-dead structure in every
  file, still requires non-trivial ONElib surgery in the footer parser
  (the `&` handler asserts/`memcpy`s through the index), and has to be applied at
  4 separate open sites.

- **(B) Demote `C` from object (`O`) to data line (`D`) — drop the index from the
  format.** **Chosen.** No index is written (smaller files, smaller footer),
  nothing to skip (cheap open *automatically*, at *every* open site, no flag),
  and the format is simpler. Breaks on-disk compatibility — acceptable: the
  format is intentionally young (0.x), `tui_load` already hard-rejects any version
  ≠ the current one, and the policy is "every format change bumps the version and
  regenerates everything."

## 4. Format changes (version 0.2)

```
  before (0.1)                              after (0.2)
  ----------------------------------------  ----------------------------------------
  O S 4 6 STRING 3 INT 3 INT 3 INT          O S 3 6 STRING 3 INT 3 INT
    seqName, seqLen, first_c_ord, n_chunks    seqName, seqLen, n_chunks
  O C 4 3 INT 3 INT 3 INT 3 INT             D C 4 3 INT 3 INT 3 INT 3 INT
    g_min, g_span, t_min, t_span              g_min, g_span, t_min, t_span  (fields unchanged)
  t = {T, 0, 1}                             t = {T, 0, 2}
```

- `C`: **`O` → `D`**. Fields unchanged. The footer no longer carries a `C` `&`
  index (and ONElib never allocates/loads one).
- `S`: **drop `first_c_ord`** (4 → 3 fields). It only existed to seed
  `oneGoto('C', first_c_ord+k)`. `n_chunks` is kept (pre-allocation + bounds; the
  scan could instead stop at the next `S`, see §8 open decision).
- `d`, `g` stay objects (the directory binary-search and roster lookup still need
  their small indexes). `X` stays a data line, unchanged.
- Net file size: the footer loses 8 bytes × #chunks of index; nothing is added.
  **Files get smaller.**

## 5. Reader / lift changes (tui.c)

- **`TGLChunk`** (tui.c:1917): replace `int64_t c_ord` with `int64_t c_byte`
  (the file offset of the chunk's `C` line — i.e. `ftello(of->f)` taken *before*
  the `oneReadLine` that reads the `C`; this is byte-identical to what `oneGoto`
  would have seeked to, ONElib.c:2424/1771).

- **`tui_genome_lift_load` pass-2** (tui.c:2156-2177): per target seq,
  `oneGoto('S', ord)` (still an object — works), then **walk its `(C,R)` data
  lines sequentially**, recording `c_byte = ftello(of->f)` before each `C`,
  reading `g_min/g_span`, and **skipping the R blob** to reach the next `C`.

- **`chunk_decode`** (tui.c:2230): `fseeko(of->f, ch->c_byte, SEEK_SET)` then
  `oneReadLine` (C), `oneReadLine` (R), decode — replacing `oneGoto('C', ch->c_ord)`.
  The decode itself is unchanged; the R list-codec is navigation-order-independent
  (deserialized once from the footer `;` line), so decoded runs are **byte-identical**.
  Manual `fseeko` + `oneReadLine` is stdio-coherent — `vfGetc` is plain `getc`
  with no second buffer layer (ONElib.c:752-757, verified).

### 5.1 Skipping R during the scan (the one real mechanism choice)

Pass-2 must advance past each ~1.6 MB R blob **without reading it** (today's C→C
index jump reads *zero* R bytes; a naïve sequential scan that `oneReadLine`s R
would `fread` the entire blob — reading a whole genome's R data during load,
a severe regression). The on-disk skip distance is fully derivable from bytes
already present (no format change), but ONElib has no header-only-skip API today.

**Recommended:** add a small (~15-line) ONElib primitive
`oneReadLineSkipList(vf)` that reads the line's fixed fields + the list field's
length header and `fseeko`s past the list payload (branching on the per-line
`x&0x1` compression flag: compressed ⇒ skip `(nBits+7)>>3` bytes; uncompressed ⇒
`listLen * listEltSize`). It is general, reusable, and is also exactly what a
future remote reader needs (skip without fetching). This is the **only** ONElib
change required by this design.

Alternatives (see §8): de-interleave C/R so the scan reads a contiguous metadata
block (cleaner scan, heavier writer change); or read R fully (simple, but the
load-time regression makes it a non-starter at 577-way scale).

## 6. Writer changes

- **`tui_write_sequence`** (tui.c:1147): write `C` with the data-line opcode;
  drop the `c_ord_emit` counter and the `first_c_ord` field on `S` (3 fields).
- **`tui_write_header`**: bump the version ints to `0, 2`.
- `taf_chain.c write_genome` already emits through `tui_write_sequence`, so it
  inherits the change for free (no hand-rolled layout — see the centralized-writer
  invariant).

## 7. What does **not** change

- The `d`/`S`/`g` indexes and the directory binary search; `tui_find_d`,
  `tui_seq_length`, `tui_genome_seq_lengths`, `tui_genome_names`.
- `tui_query` and `tui_load_seq_runs` — already sequential after `oneGoto('S')`.
- The `X` anchor path and `tui_extract_*`.
- The public lift API (`tui_genome_lift_load/column/visit_runs/stream_runs`),
  blockViz, and all CLI tools — the change is entirely under the API.
- Decode logic, the run encoding, chaining, everything downstream of a decoded chunk.

**No `oneFileOpenReadSkipIndex` flag, and no per-open-site changes** — every open
(shared handle, lift handle, the `GetSpecies` roster open, the `GetChroms`
seq-length open) becomes cheap automatically because the C index no longer exists.

## 8. Open decisions

1. **R-skip mechanism** — `oneReadLineSkipList` primitive (recommended) vs.
   de-interleaving C and R (cleaner scan + better for remote, but a two-pass /
   buffered writer). Recommend the primitive now; revisit de-interleave when we
   build remote access.
2. **Keep `n_chunks` on `S`?** It's now optional (the scan can stop at the next
   `S`/`d`). Recommended to keep it for pre-allocation and as a sanity bound.
3. **Version bump** — `0.2` (minor) is proposed; the gate is exact-match so the
   distinction is cosmetic. Stay 0.x.

## 9. Migration & validation

- Bump to **0.2**; `tui_load` already rejects any other version, so old `.tui`
  fail loudly with the regenerate message. Regenerate all `.tui` (apes, rodents,
  577-way, chains) — same drill as the 0.1 cutover.
- **Tests:** `bin/stTafTests` (94 C tests) and `python3 -m pytest tests/tui/test_tui.py`
  (70) — build per the htslib note. Plus the **byte-identical lift guard**:
  `g0`-chained lift == base lift (known baseline 2814 blocks / 173128 cov on the
  evolver fixture). **Choose a lift target whose chunks span both the
  uncompressed-R region (first ~100 KB, before the codec trains) and a later
  compressed-R region**, so both `fseeko`/skip paths (`x&0x1` = 0 and 1) are
  exercised.
- **Measure:** open wall-time + RSS on the apes `.tui` before/after; confirm the
  footer shrinks and RSS drops to the d/S/g level. Extrapolate to 577-way.

## 10. Risks & mitigations (from adversarial review)

| Risk | Mitigation |
|---|---|
| `oneReadLineSkipList` must branch on per-line compression (`x&0x1`); early R blobs are uncompressed, later ones compressed (codec trains after 100 KB). | Branch explicitly; hard-code the STRING list path; assert `fieldType[listField]==oneSTRING` so a future schema change trips loudly. Test across the codec boundary (§9). |
| Manual `fseeko` desyncs ONElib's internal accumulators. | tui never reads those accumulators (navigates via S/n_chunks); `oneReadLine` correctness doesn't consult them. Verified. |
| A future caller reintroduces `oneGoto('C')`. | It will fail to compile against the new schema (C is not an object) — fails loud, not silent. |
| Scan reads more I/O than the old C→C jump. | The skip primitive reads only each R's ~12-byte header (vs the full blob), and only for the *target* genome — far below the ~1 GB index it eliminates. |

## 11. Why this is the right foundation for remote access

With `C` demoted, the only resident indexes are `d` + `S` (+ `g`), all `O(n_seqs)`
— small enough to download. Chunk navigation is already "seek to a byte offset
recorded from a sequential scan, fetch C+R" — **which is exactly the range-request
shape a remote reader needs.** The `oneReadLineSkipList` primitive doubles as the
remote "skip without fetching" operation. So option A/C here is a clean stepping
stone to option B (htslib `hFILE` range I/O), not a detour.

## 12. Appendix — key evidence (file:line)

```
ONElib.c:1546-1550   per-type index alloc on '#' count line (the cost source)
ONElib.c:1625-1631   '&' footer fill (memcpy from compressed blob)
ONElib.c:1765-1799   oneGoto = index[i] + fseek; cross-type accum loop tolerates NULL index (1782)
ONElib.c:2424,1771   index[i] is the line-start offset captured before the type byte == ftello target
ONElib.c:752-757     vfGetc is plain getc (manual fseeko + oneReadLine is coherent)
tui.c:40,1770,1876   X = (col, source-MAF fpos); tui_query/tui_load_seq_runs are sequential (no C index)
tui.c:2156-2177      pass-2 scan (oneGoto('C'), skips R)  -> sequential + record c_byte
tui.c:2230-2236      chunk_decode (oneGoto('C'))          -> fseeko(c_byte)
tui.c:1147-1164,1151 writer: tui_write_sequence + first_c_ord seed
tui.c:68-71          schema: D X, O d, O S 4..., O C 4... (C is currently an object)
```
