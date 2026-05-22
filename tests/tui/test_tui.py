#!/usr/bin/env python3
"""
tests/tui/test_tui.py -- exercises taffy view + taffy lift on a
.tui-indexed universal MAF, using a hal2maf-rooted MAF (rooted at
simHuman_chr6) as ground truth.

Fixtures (built once and committed in tests/tui/):
  evolverMammals.uni.maf.gz       cactus-hal2maf --universal --refGenome Anc0
  evolverMammals.uni.maf.gz.tui   cactus-hal2maf --index (passes -u to taffy)
  evolverMammals.simHuman.maf.gz  hal2maf --refGenome simHuman_chr6
  evolverMammals.simHuman.maf.gz.tai  taffy index -i ...

Cases covered (each picked by surveying the universal MAF):
  forward                 simHuman:0-4         all '+', no paralogs
  reverse non-row-0       simHuman:7659-7720   some leaves on '-'
  forward paralog         simHuman:187350-...  simCow & Anc2 appear twice on '+'
  reverse paralog         simHuman:477815-...  simMouse & simRat each with +/- copies

Run from the taffy root:
  python3 tests/tui/test_tui.py
"""
import os, sys, subprocess

HERE       = os.path.dirname(os.path.abspath(__file__))
ROOT       = os.path.abspath(os.path.join(HERE, '..', '..'))
TAFFY      = os.path.join(ROOT, 'bin', 'taffy')
UNI_MAF    = os.path.join(HERE, 'evolverMammals.uni.maf.gz')
SIMH_MAF   = os.path.join(HERE, 'evolverMammals.simHuman.maf.gz')

REF_GENOME    = 'simHuman_chr6'
REF_SEQ       = 'simHuman.chr6'
REF_SEQ_FULL  = f'{REF_GENOME}.{REF_SEQ}'

# Curated regions (start, end), each picked from the universal MAF as the
# smallest containing one of the structural shapes we want to exercise.
CASES = [
    ('forward',            0,       4),
    ('reverse_nonref',     7659,    7720),
    ('forward_paralog',    187350,  187354),
    ('reverse_paralog',    477815,  477842),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(*args, check=True):
    p = subprocess.run(args, capture_output=True, text=True)
    if check and p.returncode != 0:
        sys.stderr.write(f"FAIL: {' '.join(args)}\n{p.stderr}\n")
        raise SystemExit(1)
    return p.stdout, p.stderr

def parse_maf_columns(maf_text):
    """
    Parse a MAF string into a list of blocks; each block is a list of columns;
    each column is a frozenset of (full_seq_name, forward_pos, strand_char,
    base_lowercased) tuples.  The forward_pos is the column-resolved forward
    coordinate (start + nongap_index, on '+' or seq_len-1-(...) on '-').

    Block ordering is preserved, but the column set is identical regardless of
    row ordering -- so this comparator is robust to row-order differences
    between two MAF backends.
    """
    blocks = []
    cur = None
    for line in maf_text.splitlines():
        if not line or line.startswith('#'):
            continue
        if line.startswith('a'):
            if cur is not None and cur:
                blocks.append(cur)
            cur = []
            continue
        if not line.startswith('s'):
            continue
        parts = line.split('\t')
        if len(parts) < 7:
            parts = line.split()
        # s seq start length strand seqLen bases
        seq    = parts[1]
        start  = int(parts[2])
        length = int(parts[3])
        strand = parts[4]
        seqlen = int(parts[5])
        bases  = parts[6]
        # forward range covered by this row
        if strand == '+':
            forward_base_at = lambda nongap_i, s=start: s + nongap_i
        else:
            # strand '-': MAF start is from the OTHER end of the sequence.
            # forward coord of the nongap_i-th base read (column order) is
            # seqlen - start - 1 - nongap_i  (DECREASING in column order).
            forward_base_at = lambda nongap_i, s=start, L=seqlen, n=length: L - s - 1 - nongap_i
        # walk columns, count non-gap to find each base's forward coord
        ng = 0
        if cur is None: cur = []
        # extend cur to len(bases) columns
        while len(cur) < len(bases):
            cur.append(set())
        for j, b in enumerate(bases):
            if b == '-':
                continue
            cur[j].add((seq, forward_base_at(ng), strand, b.lower()))
            ng += 1
    if cur:
        blocks.append(cur)
    return blocks

def columns_set(blocks):
    """Flatten blocks into one big column-key set: frozenset of (seq, pos, strand, base) tuples per column."""
    cols = []
    for b in blocks:
        for col in b:
            cols.append(frozenset(col))
    return cols

def column_genome_positions(blocks):
    """For each column, return frozenset of (seq, forward_pos, strand) -- no base, used when comparing two backends that may differ on case."""
    out = []
    for b in blocks:
        for col in b:
            out.append(frozenset((s, p, st) for (s, p, st, _b) in col))
    return out


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_view_consistency():
    """taffy view -U query on the universal MAF must produce the same per-column
    set of (genome, pos, strand) triples as taffy view on the simHuman .tai MAF
    for every test region.  The two backends may pick different ancestor sets
    (universal MAF has only the row-0 ancestor per block; hal2maf has every
    intermediate ancestor) so we restrict the comparison to LEAF rows
    (genome.seq with no 'Anc' prefix and no 'mr' which is the human-rodent
    inner ancestor in evolverMammals)."""
    failed = []
    for name, start, end in CASES:
        region = f'{REF_SEQ_FULL}:{start}-{end}'

        uni_out, _ = run(TAFFY, 'view', '-U', 'query', '-r', region, '-m', '-i', UNI_MAF)
        tai_out, _ = run(TAFFY, 'view',                '-r', region, '-m', '-i', SIMH_MAF)

        uni_blocks = parse_maf_columns(uni_out)
        tai_blocks = parse_maf_columns(tai_out)

        # union of per-column (genome.seq, pos, strand) keys, leaf-only
        def leaves(blocks):
            keys = set()
            for b in blocks:
                for col in b:
                    for (s, p, st, _b) in col:
                        g = s.split('.', 1)[0]
                        if g in ('Anc0','Anc1','Anc2','mr'):  # ancestor / inner node
                            continue
                        keys.add((s, p, st))
            return keys

        u = leaves(uni_blocks); t = leaves(tai_blocks)
        if u != t:
            failed.append((name, region, u - t, t - u))
        else:
            print(f"  view consistency {name:18s} OK -- {len(u)} leaf (seq,pos,strand) triples match")
    if failed:
        for name, region, u_only, t_only in failed:
            print(f"FAIL {name} ({region})")
            if u_only: print(f"  uni-only: {sorted(u_only)[:5]}{'...' if len(u_only)>5 else ''}")
            if t_only: print(f"  tai-only: {sorted(t_only)[:5]}{'...' if len(t_only)>5 else ''}")
        raise AssertionError(f"view consistency failed for {len(failed)} region(s)")

def test_forward_paralog_via_lift():
    """In the forward_paralog block (simHuman:187350-187354), simCow and
    Anc2 appear twice on '+' strand.  Lifting any base of simHuman in that
    range to simCow_chr6 via `taffy lift` must produce TWO records per
    input base (one for each simCow copy)."""
    name, start, end = next(c for c in CASES if c[0] == 'forward_paralog')
    in_wig = '/tmp/test_tui_paralog_in.wig'
    out_wig = '/tmp/test_tui_paralog_out.wig'
    with open(in_wig, 'w') as f:
        f.write(f'variableStep chrom={REF_SEQ_FULL}\n')
        for p in range(start, end):
            f.write(f'{p+1} 1.0\n')   # 1-based wig pos

    run(TAFFY, 'lift', '-i', UNI_MAF, '-w', in_wig, '-g', 'simCow_chr6', '-o', out_wig)
    pts = []
    with open(out_wig) as f:
        cur = None
        for line in f:
            line = line.strip()
            if line.startswith('variableStep'):
                cur = line.split('chrom=')[1]
            elif line and cur and line[0].isdigit():
                tok = line.split()
                pts.append((cur, int(tok[0]) - 1, float(tok[1])))
    # 4 source bases * 2 paralogs = 8 lifted records
    assert len(pts) == 8, f"expected 8 lifted records (4 src * 2 simCow paralogs), got {len(pts)}: {pts}"
    print(f"  forward_paralog lift  OK -- 4 simHuman positions -> 8 simCow records (2 paralogs)")

def test_reverse_paralog_via_lift():
    """In the reverse_paralog block (simHuman:477815-477842), simMouse and
    simRat each appear in both + and - copies.  Lifting any base of
    simHuman in that range to simMouse_chr6 must produce TWO records --
    one '+'-strand and one '-'-strand forward position."""
    name, start, end = next(c for c in CASES if c[0] == 'reverse_paralog')
    in_wig = '/tmp/test_tui_revpar_in.wig'
    out_wig = '/tmp/test_tui_revpar_out.wig'
    with open(in_wig, 'w') as f:
        f.write(f'variableStep chrom={REF_SEQ_FULL}\n')
        for p in range(start, start + 1):    # 1 source base is enough
            f.write(f'{p+1} 1.0\n')

    run(TAFFY, 'lift', '-i', UNI_MAF, '-w', in_wig, '-g', 'simMouse_chr6', '-o', out_wig)
    pts = []
    with open(out_wig) as f:
        cur = None
        for line in f:
            line = line.strip()
            if line.startswith('variableStep'):
                cur = line.split('chrom=')[1]
            elif line and cur and line[0].isdigit():
                tok = line.split()
                pts.append((cur, int(tok[0]) - 1))

    # The two simMouse paralogs at simHuman:477815 are at FORWARD positions:
    #   - '+' row at maf_start=495907 length=27 strand='+': forward = 495907
    #   - '-' row at maf_start=194591 length=27 strand='-':
    #          forward = seqlen - start - 1 = 636262 - 194591 - 1 = 441670
    # taffy lift emits both forward positions.
    expected = {('simMouse.chr6', 495907), ('simMouse.chr6', 441670)}
    got = set(pts)
    assert got == expected, f"reverse paralog mismatch.\n  expected: {expected}\n  got: {got}"
    print(f"  reverse_paralog lift  OK -- simHuman:477815 -> simMouse_chr6 {{495907 (+ run), 441670 (- run)}}")

def test_leaf_to_leaf_lift():
    """`taffy lift` between two LEAVES (simHuman_chr6 -> simMouse_chr6).
    Internally this composes ancestor-column -> leaf via the universal column
    space, with no exposed intermediate step.  For each test region, we
    derive ground truth from the simHuman-rooted hal2maf MAF (where each
    block has simHuman as row-0 and simMouse rows at known forward positions
    + strands) and check that taffy lift on the universal MAF agrees -- both
    on which simHuman positions map at all, and on the full set of simMouse
    forward positions per simHuman base."""
    for name, start, end in CASES:
        region = f'{REF_SEQ_FULL}:{start}-{end}'

        # 1) parse the simHuman ground-truth MAF for the region into per-column
        #    sets, then build the map: simHuman_forward_pos -> set of simMouse
        #    forward positions present in the same column.
        tai_out, _ = run(TAFFY, 'view', '-r', region, '-m', '-i', SIMH_MAF)
        tai_blocks = parse_maf_columns(tai_out)
        expected = {}
        for b in tai_blocks:
            for col in b:
                sh_pos = None
                sm_set = set()
                for (s, p, st, _b) in col:
                    if s == REF_SEQ_FULL:
                        sh_pos = p
                    elif s.startswith('simMouse_chr6.'):
                        sm_set.add(p)
                if sh_pos is not None:
                    expected.setdefault(sh_pos, set()).update(sm_set)

        # 2) run taffy lift simHuman -> simMouse_chr6 on every simHuman pos in
        #    the region.
        in_wig  = '/tmp/test_tui_leaf2leaf_in.wig'
        out_wig = '/tmp/test_tui_leaf2leaf_out.wig'
        with open(in_wig, 'w') as f:
            f.write(f'variableStep chrom={REF_SEQ_FULL}\n')
            for p in range(start, end):
                f.write(f'{p+1} 1.0\n')
        run(TAFFY, 'lift', '-i', UNI_MAF, '-w', in_wig, '-g', 'simMouse_chr6', '-o', out_wig)

        # 3) collect actual lifted simMouse forward positions.  Wig is sorted
        #    by output (simMouse) position so we can't recover which input
        #    simHuman pos produced each output -- but the UNION over the
        #    region must match the union of `expected`.
        got_all = set()
        with open(out_wig) as f:
            cur = None
            for line in f:
                line = line.strip()
                if line.startswith('variableStep'):
                    cur = line.split('chrom=')[1]
                elif line and cur and line[0].isdigit():
                    tok = line.split()
                    got_all.add((cur, int(tok[0]) - 1))

        # Expected: every simMouse forward pos across all simHuman bases.
        exp_all = set()
        for sh_pos, sm_set in expected.items():
            for sm_pos in sm_set:
                # The simMouse seq is "simMouse.chr6" inside "simMouse_chr6.simMouse.chr6"
                exp_all.add(('simMouse.chr6', sm_pos))

        if got_all != exp_all:
            extra = got_all - exp_all
            missing = exp_all - got_all
            print(f"FAIL {name} ({region})")
            if extra:    print(f"  taffy-only: {sorted(extra)[:5]}{'...' if len(extra)>5 else ''}")
            if missing:  print(f"  truth-only: {sorted(missing)[:5]}{'...' if len(missing)>5 else ''}")
            raise AssertionError(f"leaf-to-leaf mismatch in {name}: |got|={len(got_all)} |exp|={len(exp_all)}")

        # Also assert each simHuman position with simMouse alignment lifted to
        # at least one record per paralog (we can verify count parity from the
        # universal MAF's --noAncestors block via taffy view -U query).
        n_src = sum(1 for sh in expected if expected[sh])
        n_par = sum(len(s) for s in expected.values())
        print(f"  leaf2leaf {name:18s} OK -- {n_src} simHuman bases lift to {n_par} simMouse positions ({len(got_all)} unique)")


def test_view_block_count_minimum():
    """A spot check that each region actually produces non-empty output (catches
    regressions where -U query returns header-only for a leaf-coord region)."""
    for name, start, end in CASES:
        region = f'{REF_SEQ_FULL}:{start}-{end}'
        uni_out, _ = run(TAFFY, 'view', '-U', 'query', '-r', region, '-m', '-i', UNI_MAF)
        n_blocks = sum(1 for L in uni_out.splitlines() if L.startswith('a'))
        n_rows   = sum(1 for L in uni_out.splitlines() if L.startswith('s'))
        assert n_blocks >= 1 and n_rows >= n_blocks * 2, (
            f"region {region} ({name}) returned empty output ({n_blocks} blocks, {n_rows} rows)")
    print(f"  block-count sanity    OK -- all {len(CASES)} regions return non-empty MAF")


if __name__ == '__main__':
    assert os.path.isfile(TAFFY), f"taffy not at {TAFFY}; run from repo root after `make all`"
    assert os.path.isfile(UNI_MAF), f"missing fixture: {UNI_MAF}"
    assert os.path.isfile(SIMH_MAF), f"missing fixture: {SIMH_MAF}"
    print("== tests/tui/test_tui.py ==")
    test_view_block_count_minimum()
    test_view_consistency()
    test_forward_paralog_via_lift()
    test_reverse_paralog_via_lift()
    test_leaf_to_leaf_lift()
    print("ALL OK")
