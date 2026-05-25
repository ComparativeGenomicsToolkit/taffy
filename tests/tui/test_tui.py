"""
Tests for taffy view + taffy lift on a .tui-indexed universal MAF, using a
hal2maf-rooted MAF (rooted at simHuman_chr6) as ground truth.

Fixtures (committed under tests/tui/, ~6 MB total):
  evolverMammals.uni.maf.gz       cactus-hal2maf --universal --refGenome Anc0
  evolverMammals.uni.maf.gz.tui   cactus-hal2maf --index (passes -u to taffy)
  evolverMammals.simHuman.maf.gz  hal2maf --refGenome simHuman_chr6
  evolverMammals.simHuman.maf.gz.tai  taffy index -i ...

Cases (each picked by surveying the universal MAF for the structural shape):
  forward             simHuman:0-4         all '+' strand, no paralogs
  reverse_nonref      simHuman:7659-7720   leaves on '-' strand
  forward_paralog     simHuman:187350-...  simCow & Anc2 appear twice on '+'
  reverse_paralog     simHuman:477815-...  simMouse & simRat each in +/- copies

Run from the repo root with:
    pytest tests/tui
"""
import os
import subprocess

import pytest


HERE         = os.path.dirname(os.path.abspath(__file__))
ROOT         = os.path.abspath(os.path.join(HERE, '..', '..'))
TAFFY        = os.path.join(ROOT, 'bin', 'taffy')
UNI_MAF      = os.path.join(HERE, 'evolverMammals.uni.maf.gz')
SIMH_MAF     = os.path.join(HERE, 'evolverMammals.simHuman.maf.gz')

REF_GENOME   = 'simHuman_chr6'
REF_SEQ      = 'simHuman.chr6'
REF_SEQ_FULL = f'{REF_GENOME}.{REF_SEQ}'

# (id, start, end) -- one block of the universal MAF picked per structural shape
CASES = [
    pytest.param('forward',          0,       4,       id='forward'),
    pytest.param('reverse_nonref',   7659,    7720,    id='reverse_nonref'),
    pytest.param('forward_paralog',  187350,  187354,  id='forward_paralog'),
    pytest.param('reverse_paralog',  477815,  477842,  id='reverse_paralog'),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(*args):
    """Run a subprocess, return (stdout, stderr).  Fails the test on non-zero."""
    p = subprocess.run(list(args), capture_output=True, text=True)
    assert p.returncode == 0, (
        f"\nCommand failed: {' '.join(args)}\nstderr: {p.stderr}")
    return p.stdout, p.stderr


def parse_maf_columns(maf_text):
    """
    Parse MAF text into a list of blocks; each block is a list of columns;
    each column is a set of (full_seq_name, forward_pos, strand_char,
    base_lowercased) tuples.  forward_pos is column-resolved (start +
    nongap_index on '+', seqlen-1-(start+nongap_index) on '-' — i.e. always
    the forward-strand coordinate of the base at that column).

    Robust to row ordering, so two backends can be compared per-column even
    when they emit rows in different orders.
    """
    blocks = []
    cur = None
    for line in maf_text.splitlines():
        if not line or line.startswith('#'):
            continue
        if line.startswith('a'):
            if cur:
                blocks.append(cur)
            cur = []
            continue
        if not line.startswith('s'):
            continue
        parts = line.split('\t')
        if len(parts) < 7:
            parts = line.split()
        seq, start, length, strand, seqlen, bases = (
            parts[1], int(parts[2]), int(parts[3]), parts[4],
            int(parts[5]), parts[6])
        if strand == '+':
            forward = lambda i, s=start: s + i
        else:
            forward = lambda i, s=start, L=seqlen: L - s - 1 - i
        while len(cur) < len(bases):
            cur.append(set())
        ng = 0
        for j, b in enumerate(bases):
            if b == '-':
                continue
            cur[j].add((seq, forward(ng), strand, b.lower()))
            ng += 1
    if cur:
        blocks.append(cur)
    return blocks


def _is_ancestor_genome(seq_name):
    """Return True if the genome part of 'genome.seq' is an evolverMammals
    ancestor (anything not a simXXX_chr6 leaf)."""
    g = seq_name.split('.', 1)[0]
    return g in ('Anc0', 'Anc1', 'Anc2', 'mr')


def _leaf_cols(blocks):
    """Flatten blocks to a set of (seq, forward_pos, strand) for leaf rows
    only.  Used for backend-agnostic comparisons -- the universal MAF carries
    one ancestor row per block (--noAncestors), hal2maf carries every
    intermediate, so an ancestors-inclusive comparison would be apples-to-
    oranges."""
    keys = set()
    for b in blocks:
        for col in b:
            for (s, p, st, _b) in col:
                if _is_ancestor_genome(s):
                    continue
                keys.add((s, p, st))
    return keys


def _read_wig(path):
    """Parse a wig (mixed variableStep + fixedStep) into a list of
    (chrom, 0-based pos).  taffy lift emits either format depending on
    whether consecutive same-seq target positions group into a long run."""
    pts = []
    mode = None       # 'var' | 'fixed'
    chrom = None
    fs_pos = 0        # next 0-based position for the current fixedStep block
    fs_step = 1
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('variableStep'):
                mode = 'var'
                chrom = line.split('chrom=', 1)[1].split()[0]
                continue
            if line.startswith('fixedStep'):
                mode = 'fixed'
                fs_step = 1
                for tok in line.split()[1:]:
                    if '=' not in tok: continue
                    k, v = tok.split('=', 1)
                    if   k == 'chrom': chrom = v
                    elif k == 'start': fs_pos = int(v) - 1  # 1-based -> 0-based
                    elif k == 'step':  fs_step = int(v)
                continue
            if not line or chrom is None: continue
            if mode == 'var' and line[0].isdigit():
                pts.append((chrom, int(line.split()[0]) - 1))
            elif mode == 'fixed':
                pts.append((chrom, fs_pos))
                fs_pos += fs_step
    return pts


def _write_wig(path, chrom, positions):
    """Emit a 1-based variableStep wig at the given positions (value=1.0)."""
    with open(path, 'w') as f:
        f.write(f'variableStep chrom={chrom}\n')
        for p in positions:
            f.write(f'{p + 1} 1.0\n')


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", autouse=True)
def _check_fixtures():
    """Sanity-check the binary + data files before any test runs."""
    assert os.path.isfile(TAFFY),    f"taffy not at {TAFFY}; run `make all` first"
    assert os.path.isfile(UNI_MAF),  f"missing fixture: {UNI_MAF}"
    assert os.path.isfile(SIMH_MAF), f"missing fixture: {SIMH_MAF}"


@pytest.fixture(scope="session")
def truth_blocks():
    """For each region, the simHuman .tai MAF parsed once into per-column sets.
    Used as ground truth by tests that compare against the universal MAF or
    against taffy lift output."""
    out = {}
    for p in CASES:
        name, start, end = p.values
        region = f'{REF_SEQ_FULL}:{start}-{end}'
        tai_out, _ = run(TAFFY, 'view', '-r', region, '-m', '-i', SIMH_MAF)
        out[name] = parse_maf_columns(tai_out)
    return out


# ---------------------------------------------------------------------------
# Tests -- parametrized across CASES where structure allows
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,start,end", CASES)
def test_view_returns_nonempty(name, start, end):
    """taffy view -U query at each region must return at least one block with
    at least two rows (row-0 + something).  Catches header-only regressions."""
    region = f'{REF_SEQ_FULL}:{start}-{end}'
    uni_out, _ = run(TAFFY, 'view', '-U', 'query', '-r', region, '-m', '-i', UNI_MAF)
    lines = uni_out.splitlines()
    n_blocks = sum(1 for L in lines if L.startswith('a'))
    n_rows   = sum(1 for L in lines if L.startswith('s'))
    assert n_blocks >= 1, f"region {region} returned no blocks"
    assert n_rows >= n_blocks * 2, (
        f"region {region} returned {n_blocks} blocks with only {n_rows} rows total")


@pytest.mark.parametrize("name,start,end", CASES)
def test_view_consistency_tui_vs_tai(name, start, end, truth_blocks):
    """taffy view -U query on the universal MAF (.tui-indexed) must produce
    the same per-column leaf (genome.seq, forward_pos, strand) set as
    taffy view on the simHuman-rooted MAF (.tai-indexed)."""
    region = f'{REF_SEQ_FULL}:{start}-{end}'
    uni_out, _ = run(TAFFY, 'view', '-U', 'query', '-r', region, '-m', '-i', UNI_MAF)
    uni_leaves = _leaf_cols(parse_maf_columns(uni_out))
    tai_leaves = _leaf_cols(truth_blocks[name])
    assert uni_leaves == tai_leaves


@pytest.mark.parametrize("name,start,end", CASES)
def test_leaf_to_leaf_lift(name, start, end, truth_blocks, tmp_path):
    """`taffy lift` between two leaves (simHuman -> simMouse) must agree, on
    the union of lifted simMouse forward positions, with the simMouse rows
    present in the simHuman-rooted truth MAF at the same columns."""
    region = f'{REF_SEQ_FULL}:{start}-{end}'

    # Expected: for each column in the simHuman MAF block(s) covering the
    # region, collect every simMouse forward position aligned at that column.
    expected = set()
    for b in truth_blocks[name]:
        for col in b:
            sh_in_col = any(s == REF_SEQ_FULL for (s, _p, _st, _bs) in col)
            if not sh_in_col:
                continue
            for (s, p, _st, _bs) in col:
                if s.startswith('simMouse_chr6.'):
                    expected.add((s.split('.', 1)[1], p))

    # Lift every simHuman base in the region to simMouse_chr6.
    in_wig  = str(tmp_path / 'in.wig')
    out_wig = str(tmp_path / 'out.wig')
    _write_wig(in_wig, REF_SEQ_FULL, range(start, end))
    run(TAFFY, 'lift', '-i', UNI_MAF, '-w', in_wig, '-g', 'simMouse_chr6', '-o', out_wig)
    got = set(_read_wig(out_wig))

    assert got == expected


def test_forward_paralog_record_count(tmp_path):
    """forward_paralog block: simCow appears twice on '+' for all 4 simHuman
    bases.  `taffy lift simHuman -> simCow_chr6` must emit exactly 4*2 = 8
    records (no merging, no de-dup) -- a strict count check that
    complements test_leaf_to_leaf_lift's union-based comparison."""
    start, end = 187350, 187354
    in_wig  = str(tmp_path / 'in.wig')
    out_wig = str(tmp_path / 'out.wig')
    _write_wig(in_wig, REF_SEQ_FULL, range(start, end))
    run(TAFFY, 'lift', '-i', UNI_MAF, '-w', in_wig, '-g', 'simCow_chr6', '-o', out_wig)
    got = _read_wig(out_wig)
    assert len(got) == 8, f"expected 8 records (4 src * 2 simCow paralogs), got {len(got)}"


def test_reverse_paralog_strand_flip(tmp_path):
    """reverse_paralog block: simMouse appears twice -- once '+' at
    simMouse:495907 (forward = 495907) and once '-' at maf_start 194591
    length 27 (forward = 636262 - 194591 - 1 = 441670).  Lifting simHuman:
    477815 -> simMouse_chr6 must emit BOTH forward positions.  Exercises
    the strand-flip math in the paralog scan."""
    in_wig  = str(tmp_path / 'in.wig')
    out_wig = str(tmp_path / 'out.wig')
    _write_wig(in_wig, REF_SEQ_FULL, [477815])
    run(TAFFY, 'lift', '-i', UNI_MAF, '-w', in_wig, '-g', 'simMouse_chr6', '-o', out_wig)
    got = set(_read_wig(out_wig))
    expected = {('simMouse.chr6', 495907), ('simMouse.chr6', 441670)}
    assert got == expected


# ---------------------------------------------------------------------------
# Wig-parser edge cases
# ---------------------------------------------------------------------------

def test_lift_input_position_0(tmp_path):
    """Regression: an input wig record at ancestor coord 0 must lift, not be
    silently dropped.  Earlier wig_parse stored coordinates as (void*)pos
    in an stHash and (void*)0 collided with the iterator's end sentinel,
    so any wig record at pos 0 vanished.  The streaming parser doesn't
    use stHash, but pin the behaviour with a test."""
    in_wig  = str(tmp_path / 'pos0.wig')
    out_wig = str(tmp_path / 'pos0.out.wig')
    # variableStep at pos 1 (1-based) == ancestor coord 0
    with open(in_wig, 'w') as f:
        f.write('variableStep chrom=Anc0.Anc0refChr0\n')
        f.write('1 1.0\n')
    run(TAFFY, 'lift', '-i', UNI_MAF, '-w', in_wig, '-g', 'simHuman_chr6', '-o', out_wig)
    got = _read_wig(out_wig)
    # From the simple block (chr6:0-4): Anc0.Anc0refChr0:0 column lifts to
    # exactly one simHuman position: simHuman.chr6:0.
    assert got == [('simHuman.chr6', 0)], (
        f"expected [(simHuman.chr6, 0)], got {got}")


def test_lift_fixedstep_wig(tmp_path):
    """A fixedStep wig input must lift the same way as the equivalent
    variableStep wig.  Exercises a separate parser code path."""
    out_var   = str(tmp_path / 'var.out.wig')
    out_fixed = str(tmp_path / 'fixed.out.wig')

    var_wig = str(tmp_path / 'var.wig')
    _write_wig(var_wig, REF_SEQ_FULL, [0, 1, 2, 3])
    run(TAFFY, 'lift', '-i', UNI_MAF, '-w', var_wig, '-g', 'simMouse_chr6', '-o', out_var)

    fixed_wig = str(tmp_path / 'fixed.wig')
    with open(fixed_wig, 'w') as f:
        f.write(f'fixedStep chrom={REF_SEQ_FULL} start=1 step=1\n')
        for _ in range(4):
            f.write('1.0\n')
    run(TAFFY, 'lift', '-i', UNI_MAF, '-w', fixed_wig, '-g', 'simMouse_chr6', '-o', out_fixed)

    assert _read_wig(out_var) == _read_wig(out_fixed)
    # Sanity: the four simHuman positions 0..3 map to simMouse 0..3 (the
    # simple-block alignment is 1:1 forward on this prefix).
    assert _read_wig(out_fixed) == [('simMouse.chr6', i) for i in range(4)]


def test_view_tcol_round_trip(tmp_path):
    """`taffy view -U tcol` prepends a `s tcol <col> ...` sentinel row to every
    emitted block.  Querying via that column with `-r tcol:N-N+1` must return
    the same block (modulo headers) as querying via the leaf coord that
    points at the same column."""
    # 1) Get the universal column for simHuman:0 via -U tcol.
    out1, _ = run(TAFFY, 'view', '-U', 'tcol', '-r', f'{REF_SEQ_FULL}:0-1', '-m', '-i', UNI_MAF)
    tcol_row = next((L for L in out1.splitlines() if L.startswith('s\ttcol\t')), None)
    assert tcol_row is not None, "no `tcol` sentinel row in -U tcol output"
    parts = tcol_row.split('\t')
    col, length = int(parts[2]), int(parts[3])
    assert length == 1, f"expected single-column query, got tcol length {length}"

    # 2) Re-query directly by column with `-r tcol:N-N+1`.  Should return the
    #    same block contents (rows, positions, bases) as the leaf-coord query
    #    above.
    out2, _ = run(TAFFY, 'view', '-U', 'tcol', '-r', f'tcol:{col}-{col+1}', '-m', '-i', UNI_MAF)

    # Compare per-column LEAF sets (column-keyed, row-order-agnostic, ignores
    # the prepended `tcol` sentinel which is identical in both outputs).
    assert _leaf_cols(parse_maf_columns(out1)) == _leaf_cols(parse_maf_columns(out2))


def test_lift_unknown_target_genome(tmp_path):
    """`taffy lift -g BogusGenome` must exit non-zero with a useful stderr
    message -- not segfault, hang, or quietly emit empty output that the
    caller would silently consume."""
    in_wig = str(tmp_path / 'in.wig')
    _write_wig(in_wig, REF_SEQ_FULL, [0])
    p = subprocess.run(
        [TAFFY, 'lift', '-i', UNI_MAF, '-w', in_wig, '-g', 'NoSuchGenome',
         '-o', str(tmp_path / 'out.wig')],
        capture_output=True, text=True)
    assert p.returncode != 0, (
        f"expected non-zero exit for unknown target genome, "
        f"got rc={p.returncode}, stderr={p.stderr!r}")
    assert 'NoSuchGenome' in p.stderr or 'not found' in p.stderr.lower(), (
        f"stderr should mention the offending genome or 'not found'; "
        f"got: {p.stderr!r}")


def test_lift_identity_simhuman_to_simhuman(tmp_path):
    """For every simHuman input position that's aligned in the universal
    MAF, lifting back to simHuman must include that same position in the
    output (possibly with extra paralog records at the same column).
    Positions that don't appear in the alignment (lineage-specific
    insertions) are excluded -- they legitimately produce no output."""
    # Drawn from CASES: each starts a block whose simHuman row appears.
    positions = [c.values[1] for c in CASES]   # 0, 7659, 187350, 477815
    in_wig  = str(tmp_path / 'in.wig')
    out_wig = str(tmp_path / 'out.wig')
    _write_wig(in_wig, REF_SEQ_FULL, positions)
    run(TAFFY, 'lift', '-i', UNI_MAF, '-w', in_wig, '-g', REF_GENOME, '-o', out_wig)
    got = {p for (_seq, p) in _read_wig(out_wig)}
    missing = set(positions) - got
    assert not missing, (
        f"identity lift dropped {len(missing)} simHuman position(s): "
        f"{sorted(missing)}")


def test_view_empty_and_out_of_range(tmp_path):
    """Region-parsing edge cases on the .tui code path:
      * a zero-length region (start == end) is intentionally rejected
        (exit 1, "Invalid region" on stderr) -- documented as half-open
        [start, end);
      * a syntactically valid range past the end of the sequence returns
        a clean header-only MAF with exit 0;
      * an unknown contig name behaves the same way as out-of-range.
    Mirrors the tai Bug 1 regression on the .tui side."""
    # zero-length: rejected with rc != 0
    p = subprocess.run(
        [TAFFY, 'view', '-U', 'query', '-r', f'{REF_SEQ_FULL}:100-100',
         '-m', '-i', UNI_MAF],
        capture_output=True, text=True)
    assert p.returncode != 0, "zero-length region should be rejected"
    assert 'Invalid region' in p.stderr, (
        f"expected 'Invalid region' on stderr, got: {p.stderr!r}")

    # past-end of seq: valid syntax, no overlap, header-only output, rc=0
    out, _ = run(TAFFY, 'view', '-U', 'query', '-r', f'{REF_SEQ_FULL}:601863-601864',
                 '-m', '-i', UNI_MAF)
    assert sum(1 for L in out.splitlines() if L.startswith('a')) == 0
    assert out.startswith('##maf')

    # unknown contig: same shape as out-of-range -- header-only, rc=0
    out, _ = run(TAFFY, 'view', '-U', 'query', '-r', 'NoSuchSeq:0-10',
                 '-m', '-i', UNI_MAF)
    assert sum(1 for L in out.splitlines() if L.startswith('a')) == 0
    assert out.startswith('##maf')


def test_lift_multichrom_wig(tmp_path):
    """A wig that switches between two ancestor chroms exercises the
    source-seq cache refresh in taf_lift.c: when a new `variableStep
    chrom=...` line is seen, the previous cached run table must be
    flushed and the new one loaded.  Both sources align to the same
    simMouse column (simMouse.chr6:0) so we expect TWO output records
    pointing there (one per input record, not merged)."""
    in_wig  = str(tmp_path / 'multi.wig')
    out_wig = str(tmp_path / 'multi.out.wig')
    # Two source chroms, one record each.  Both lift to the same simMouse
    # column because Anc0.Anc0refChr0:0, simHuman.chr6:0, and simMouse.chr6:0
    # are all at the same universal column (the simple block at chr6:0-4).
    with open(in_wig, 'w') as f:
        f.write('variableStep chrom=Anc0.Anc0refChr0\n')
        f.write('1 1.0\n')
        f.write(f'variableStep chrom={REF_SEQ_FULL}\n')
        f.write('1 1.0\n')
    run(TAFFY, 'lift', '-i', UNI_MAF, '-w', in_wig, '-g', 'simMouse_chr6', '-o', out_wig)
    got = _read_wig(out_wig)
    assert len(got) == 2, f"expected 2 records (one per source), got {len(got)}: {got}"
    assert all(p == ('simMouse.chr6', 0) for p in got), (
        f"both sources should lift to simMouse.chr6:0; got {got}")


def _read_bed(path):
    """Parse a BED file into a list of dicts."""
    out = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n').rstrip('\r')
            if not line or line.startswith('#'): continue
            f = line.split('\t')
            d = {'chrom': f[0], 'start': int(f[1]), 'end': int(f[2])}
            if len(f) >= 4: d['name']   = f[3]
            if len(f) >= 5: d['score']  = f[4]
            if len(f) >= 6: d['strand'] = f[5]
            out.append(d)
    return out


def test_lift_bed_equivalent_to_wig(tmp_path):
    """A BED interval covering [s,e) should produce target intervals whose
    base set matches what the same range gives via per-position wig lift.
    Exercises the -b path against the -w path (which has its own tests)."""
    region_start, region_end = 10000, 11000
    in_bed  = str(tmp_path / 'in.bed')
    out_bed = str(tmp_path / 'out.bed')
    in_wig  = str(tmp_path / 'in.wig')
    out_wig = str(tmp_path / 'out.wig')
    with open(in_bed, 'w') as f:
        f.write(f'{REF_SEQ_FULL}\t{region_start}\t{region_end}\tg1\t0\t+\n')
    _write_wig(in_wig, REF_SEQ_FULL, range(region_start, region_end))
    run(TAFFY, 'lift', '-i', UNI_MAF, '-b', in_bed, '-g', 'simMouse_chr6', '-o', out_bed)
    run(TAFFY, 'lift', '-i', UNI_MAF, '-w', in_wig, '-g', 'simMouse_chr6', '-o', out_wig)

    # Expand BED intervals into a (chrom, pos) set; compare to wig output's set.
    bed_set = set()
    for r in _read_bed(out_bed):
        for p in range(r['start'], r['end']):
            bed_set.add((r['chrom'], p))
    wig_set = set(_read_wig(out_wig))
    assert bed_set == wig_set, (
        f"bed vs wig base set differs: only-in-bed={len(bed_set - wig_set)} "
        f"only-in-wig={len(wig_set - bed_set)}")


def test_lift_bed_strand_xor(tmp_path):
    """An input BED on '-' lifted via a '-' alignment should output '+'.
    For [10000,11000) on REF_SEQ which lifts to simMouse via reverse-strand
    alignment, '-' input should produce all-'+' output records."""
    in_bed  = str(tmp_path / 'in.bed')
    out_bed = str(tmp_path / 'out.bed')
    with open(in_bed, 'w') as f:
        f.write(f'{REF_SEQ_FULL}\t10000\t11000\tg1\t99\t-\n')
    run(TAFFY, 'lift', '-i', UNI_MAF, '-b', in_bed, '-g', 'simMouse_chr6', '-o', out_bed)
    recs = _read_bed(out_bed)
    assert len(recs) > 0, "expected some output for [10000,11000)"
    strands = {r['strand'] for r in recs}
    assert strands == {'+'}, f"input '-' through reverse alignment should give '+'; got {strands}"
    # Score + name passthrough.
    assert all(r['name'] == 'g1' and r['score'] == '99' for r in recs)


def test_lift_bed_bed3_in_bed3_out(tmp_path):
    """BED3 in -> BED3 out (no name, score, strand columns added)."""
    in_bed  = str(tmp_path / 'in.bed')
    out_bed = str(tmp_path / 'out.bed')
    with open(in_bed, 'w') as f:
        f.write(f'{REF_SEQ_FULL}\t10000\t11000\n')
    run(TAFFY, 'lift', '-i', UNI_MAF, '-b', in_bed, '-g', 'simMouse_chr6', '-o', out_bed)
    with open(out_bed) as f:
        for line in f:
            line = line.rstrip()
            if not line: continue
            cols = line.split('\t')
            assert len(cols) == 3, f"BED3 in should give BED3 out, got {len(cols)} cols: {line}"


def test_lift_bed_unknown_source_seq(tmp_path):
    """Unknown source seq produces empty output (warning logged, no error)."""
    in_bed  = str(tmp_path / 'in.bed')
    out_bed = str(tmp_path / 'out.bed')
    with open(in_bed, 'w') as f:
        f.write('nonexistent_seq\t0\t100\tx\n')
    run(TAFFY, 'lift', '-i', UNI_MAF, '-b', in_bed, '-g', 'simMouse_chr6', '-o', out_bed)
    assert os.path.getsize(out_bed) == 0


def test_lift_bed_empty_input(tmp_path):
    """Empty BED input -> empty output, exit success."""
    in_bed  = str(tmp_path / 'in.bed')
    out_bed = str(tmp_path / 'out.bed')
    open(in_bed, 'w').close()
    run(TAFFY, 'lift', '-i', UNI_MAF, '-b', in_bed, '-g', 'simMouse_chr6', '-o', out_bed)
    assert os.path.getsize(out_bed) == 0


@pytest.mark.parametrize("n_cols", [4, 5])
def test_lift_bed_arity_preserved(tmp_path, n_cols):
    """BED4 / BED5 in -> BED4 / BED5 out: extras propagate."""
    in_bed  = str(tmp_path / 'in.bed')
    out_bed = str(tmp_path / 'out.bed')
    fields = [REF_SEQ_FULL, '10000', '11000', 'g1', '42']
    with open(in_bed, 'w') as f:
        f.write('\t'.join(fields[:n_cols]) + '\n')
    run(TAFFY, 'lift', '-i', UNI_MAF, '-b', in_bed, '-g', 'simMouse_chr6', '-o', out_bed)
    with open(out_bed) as fh:
        for line in fh:
            line = line.rstrip()
            if not line: continue
            cols = line.split('\t')
            assert len(cols) == n_cols, (
                f"BED{n_cols} in should give BED{n_cols} out, got {len(cols)}: {line}")
            assert cols[3] == 'g1'
            if n_cols >= 5: assert cols[4] == '42'


def test_lift_bed_strand_dot_omits_strand_column(tmp_path):
    """BED6 input with strand '.' (unstranded) should produce BED5-shaped
    output -- no strand col, since '.' isn't a definite orientation."""
    in_bed  = str(tmp_path / 'in.bed')
    out_bed = str(tmp_path / 'out.bed')
    with open(in_bed, 'w') as f:
        f.write(f'{REF_SEQ_FULL}\t10000\t11000\tg1\t0\t.\n')
    run(TAFFY, 'lift', '-i', UNI_MAF, '-b', in_bed, '-g', 'simMouse_chr6', '-o', out_bed)
    recs = _read_bed(out_bed)
    assert len(recs) > 0
    for r in recs:
        assert 'strand' not in r, f"input '.' should omit strand col, got {r}"


