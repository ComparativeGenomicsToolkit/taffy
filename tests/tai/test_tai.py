"""
Tests for taffy view + taffy index (.tai) on the evolverMammals MAF.  This is
the legacy non-universal indexing path; tests cover format permutations (MAF,
TAF, RLE-TAF), bgzip on/off, rename-via-name-map, plus three specific
regression cases.

Truth fixtures (committed in tests/tai/):
  evolverMammals_<start>_<stop>.maf   region MAF (mafExtractor output)
  evolverMammals_<start>_<stop>.taf   region TAF (taffy view output)
  evolverMammals_subregions.bed       region list
  evolverMammals_ancestor_name_map.tsv      Anc{N} -> Ancestor{N}
  evolverMammals_ancestor_rev_name_map.tsv  reverse map

Run from the repo root with:
    pytest tests/tai
"""
import os
import shutil
import subprocess

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
TAFFY = os.path.join(ROOT, 'bin', 'taffy')
SRC_MAF = os.path.join(ROOT, 'tests', 'evolverMammals.maf')
REGIONS_BED = os.path.join(HERE, 'evolverMammals_subregions.bed')
NAME_MAP     = os.path.join(HERE, 'evolverMammals_ancestor_name_map.tsv')
REV_NAME_MAP = os.path.join(HERE, 'evolverMammals_ancestor_rev_name_map.tsv')


# ---------------------------------------------------------------------------
# Static config
# ---------------------------------------------------------------------------

def _load_regions():
    out = []
    with open(REGIONS_BED) as f:
        for line in f:
            tok = line.split()
            if len(tok) >= 3:
                out.append((tok[0], int(tok[1]), int(tok[2])))
    return out

REGIONS = _load_regions()

# (variant_id, format, bgzip, blocksize, rename)
VARIANTS = [
    ('taf_b111',                  'taf',     False, 111, False),
    ('taf_bgzip_b200',            'taf',     True,  200, False),
    ('rle_b111',                  'rle_taf', False, 111, False),
    ('rle_bgzip_b200',            'rle_taf', True,  200, False),
    ('maf_b111',                  'maf',     False, 111, False),
    ('maf_bgzip_b200',            'maf',     True,  200, False),
    ('maf_bgzip_b200_rename',     'maf',     True,  200, True),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(*args):
    p = subprocess.run([str(a) for a in args], capture_output=True, text=True)
    assert p.returncode == 0, (
        f"\nCommand failed (rc={p.returncode}): {' '.join(str(a) for a in args)}\n"
        f"stderr: {p.stderr}")
    return p.stdout, p.stderr


def maf_strip(text, filter_empty=True):
    """Normalise a MAF for comparison: strip '#' and 'a' lines, collapse runs of
    whitespace to single spaces, drop blank lines, and optionally drop rows
    whose length field (column 4 of the s-line) is 0.  Matches the legacy
    maf_compare helper but in pure Python."""
    out = []
    for line in text.splitlines():
        if not line or line.startswith('#') or line.startswith('a'):
            continue
        # tab -> space, collapse multi-space
        line = line.replace('\t', ' ')
        while '  ' in line:
            line = line.replace('  ', ' ')
        if filter_empty:
            parts = line.split(' ')
            # s-line: ['s', seq, start, length, strand, seqlen, bases]
            if len(parts) >= 4 and parts[3] == '0':
                continue
        out.append(line)
    return '\n'.join(out)


# ---------------------------------------------------------------------------
# Session fixtures: source MAF + TAF + RLE-TAF prepared once
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", autouse=True)
def _check_binary():
    assert os.path.isfile(TAFFY),   f"taffy not at {TAFFY}; run `make all` first"
    assert os.path.isfile(SRC_MAF), f"missing source MAF: {SRC_MAF}"


@pytest.fixture(scope="session")
def src_files(tmp_path_factory):
    """Copy the source MAF + generate TAF and RLE-TAF in a session tmp dir."""
    base = tmp_path_factory.mktemp("tai_src")
    maf = base / "evolverMammals.maf"
    shutil.copy(SRC_MAF, maf)
    taf = base / "evolverMammals.taf"
    rle = base / "evolverMammals.rle.taf"
    run(TAFFY, 'view', '-i', maf, '-o', taf)
    run(TAFFY, 'view', '-i', maf, '-o', rle, '-u')
    return {'maf': str(maf), 'taf': str(taf), 'rle_taf': str(rle)}


@pytest.fixture(scope="session")
def indexed_variants(src_files, tmp_path_factory):
    """For each VARIANTS entry, prepare an indexed file (bgzip + rename as
    needed) and return a dict variant_id -> path."""
    base = tmp_path_factory.mktemp("tai_variants")
    out = {}
    for vid, fmt, bgz, bsize, rename in VARIANTS:
        src = src_files[fmt]
        ext = '.taf' if fmt in ('taf', 'rle_taf') else '.maf'
        path = str(base / f"{vid}{ext}")
        shutil.copy(src, path)
        if bgz:
            run('bash', '-c', f'bgzip -f {path}')
            path = path + '.gz'
        if rename:
            renamed = path + '.renamed'
            cmd = [TAFFY, 'view', '-i', path, '-n', NAME_MAP, '-o', renamed]
            if bgz: cmd.append('-c')
            run(*cmd)
            path = renamed
        run(TAFFY, 'index', '-i', path, '-b', str(bsize))
        out[vid] = path
    return out


# ---------------------------------------------------------------------------
# Region-extraction tests -- parametrized over 7 variants * 8 regions
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("variant", [v[0] for v in VARIANTS])
@pytest.mark.parametrize("region", REGIONS, ids=lambda r: f"{r[1]}-{r[2]}")
def test_region_extraction_maf(variant, region, indexed_variants, tmp_path):
    """Extract `region` as MAF and compare (clean-strip + diff) against the
    committed truth MAF.  Works for all 7 variants: the renamed variant
    needs the reverse name map at query time."""
    rename = next(v for v in VARIANTS if v[0] == variant)[4]
    contig, start, end = region
    in_path = indexed_variants[variant]
    out = str(tmp_path / f"{variant}_{start}_{end}.maf")
    cmd = [TAFFY, 'view', '-i', in_path, '-r', f'{contig}:{start}-{end}', '-m', '-o', out]
    if rename: cmd += ['-n', REV_NAME_MAP]
    run(*cmd)

    truth = os.path.join(HERE, f'evolverMammals_{start}_{end}.maf')
    assert os.path.getsize(truth) > 100, f"truth too small: {truth}"
    with open(out)   as f: a = maf_strip(f.read())
    with open(truth) as f: b = maf_strip(f.read())
    assert a == b


@pytest.mark.parametrize("variant", [v[0] for v in VARIANTS])
@pytest.mark.parametrize("region", REGIONS, ids=lambda r: f"{r[1]}-{r[2]}")
def test_region_extraction_taf(variant, region, indexed_variants, tmp_path):
    """Extract `region` as TAF and bytewise-diff against the committed truth TAF."""
    rename = next(v for v in VARIANTS if v[0] == variant)[4]
    contig, start, end = region
    in_path = indexed_variants[variant]
    out = str(tmp_path / f"{variant}_{start}_{end}.taf")
    cmd = [TAFFY, 'view', '-i', in_path, '-r', f'{contig}:{start}-{end}', '-o', out]
    if rename: cmd += ['-n', REV_NAME_MAP]
    run(*cmd)

    truth = os.path.join(HERE, f'evolverMammals_{start}_{end}.taf')
    assert os.path.getsize(truth) > 100, f"truth too small: {truth}"
    with open(out)   as f: a = f.read()
    with open(truth) as f: b = f.read()
    assert a == b


# ---------------------------------------------------------------------------
# `taffy stats -s` sanity check (TAF variants only)
# ---------------------------------------------------------------------------

ANC0_TRUE_STATS = '\n'.join(sorted([
    'Anc0.Anc0refChr0\t4151', 'Anc0.Anc0refChr10\t14504',
    'Anc0.Anc0refChr11\t38002', 'Anc0.Anc0refChr1\t3407',
    'Anc0.Anc0refChr2\t269145', 'Anc0.Anc0refChr3\t165',
    'Anc0.Anc0refChr4\t13557', 'Anc0.Anc0refChr5\t50896',
    'Anc0.Anc0refChr6\t22717', 'Anc0.Anc0refChr7\t1851',
    'Anc0.Anc0refChr8\t111467', 'Anc0.Anc0refChr9\t4824',
]))

@pytest.mark.parametrize("variant", [v[0] for v in VARIANTS if v[1] in ('taf', 'rle_taf')])
def test_taffy_stats_seq_anc0(variant, indexed_variants):
    """taffy stats -s on the TAF variants must report the expected Anc0
    sequence lengths (renamed variants would need different IDs, so we
    only test the un-renamed TAF/RLE variants)."""
    seq_stats, _ = run(TAFFY, 'stats', '-s', '-i', indexed_variants[variant])
    got = '\n'.join(sorted(seq_stats.strip().split('\n')))
    assert got == ANC0_TRUE_STATS


# ---------------------------------------------------------------------------
# 1-bp consistency: MAF vs TAF at every step-th position of Anc0refChr0
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("pos", list(range(0, 4151, 5)), ids=lambda p: str(p))
def test_1bp_taf_eq_maf(pos, src_files, tmp_path):
    """1-bp queries against the MAF and the TAF must produce byte-identical
    output -- catches divergences in the .tai-driven random-access path
    between the two formats."""
    # Re-index with blocksize 10000 (matches legacy test) once per session;
    # we don't have a session-scoped fixture for that so build per-call --
    # at this scale it's cheap (~ms) and avoids cross-test ordering deps.
    maf = src_files['maf']; taf = src_files['taf']
    # The tai is built once per session via fixture below.
    chr_name = 'Anc0.Anc0refChr0'
    maf_out = str(tmp_path / f"{pos}.maf.taf")
    taf_out = str(tmp_path / f"{pos}.taf.taf")
    run(TAFFY, 'view', '-i', maf, '-o', maf_out, '-r', f'{chr_name}:{pos}')
    run(TAFFY, 'view', '-i', taf, '-o', taf_out, '-r', f'{chr_name}:{pos}')
    with open(maf_out) as f: a = f.read()
    with open(taf_out) as f: b = f.read()
    assert a == b


@pytest.fixture(scope="session", autouse=True)
def _index_1bp_sources(src_files):
    """Build the .tai indexes the 1bp test depends on (b=10000)."""
    run(TAFFY, 'index', '-i', src_files['maf'], '-b', '10000')
    run(TAFFY, 'index', '-i', src_files['taf'], '-b', '10000')


# ---------------------------------------------------------------------------
# Regression tests for three documented `taffy view -r` bugs
# (see cactus-phast/taffy-bug.md in the original test_tai.py history)
# ---------------------------------------------------------------------------

def _count_alignment_blocks(maf_text):
    return sum(1 for line in maf_text.splitlines() if line.startswith('a'))


def test_view_missing_region_does_not_error(tmp_path):
    """Bug 1: querying a region not in the index must return header-only
    MAF with exit 0 (not exit 1 with a 'Region not found' stderr)."""
    maf_path = str(tmp_path / 'bug1.maf')
    with open(maf_path, 'w') as f:
        f.write('##maf version=1 scoring=N/A\n\n')
        f.write('a\ns\trefA.chrA\t0\t5\t+\t100\tACGTA\ns\tother.chrZ\t0\t5\t+\t50\tACGTA\n')
    run(TAFFY, 'index', '-i', maf_path)

    # contig not in index
    out, _ = run(TAFFY, 'view', '-i', maf_path, '-m', '-r', 'NoSuchContig:0-10')
    assert _count_alignment_blocks(out) == 0
    assert out.startswith('##maf')

    # contig in index but range doesn't overlap any block
    out, _ = run(TAFFY, 'view', '-i', maf_path, '-m', '-r', 'refA.chrA:50-60')
    assert _count_alignment_blocks(out) == 0

    # bgzipped empty-region output must be a valid header-only .maf.gz (not 0-byte).
    # The early-return path used to skip LW_destruct -> bgzip writer not flushed.
    gz_out = str(tmp_path / 'bug1.out.gz')
    run(TAFFY, 'view', '-i', maf_path, '-m', '-c', '-r', 'NoSuchContig:0-10', '-o', gz_out)
    assert os.path.getsize(gz_out) > 0, "bgzipped empty output should not be 0 bytes"
    p = subprocess.run(['bgzip', '-dc', gz_out], capture_output=True, text=True, check=True)
    assert p.stdout.startswith('##maf')
    assert _count_alignment_blocks(p.stdout) == 0


def test_view_clip_zero_length_row(tmp_path):
    """Bug 2: clip_alignment must not assert on length=0 (all-gap) rows."""
    maf_path = str(tmp_path / 'bug2.maf')
    with open(maf_path, 'w') as f:
        f.write('##maf version=1 scoring=N/A\n\n')
        # block 1 has an all-gap row (length 0); range fully covers the block
        # so neither left nor right trim runs in clip_alignment
        f.write('a\ns\trefA.chrA\t0\t5\t+\t100\tACGTA\ns\tother.chrZ\t0\t0\t+\t50\t-----\n\n')
        f.write('a\ns\trefA.chrA\t5\t5\t+\t100\tTTTTT\ns\tother.chrZ\t0\t5\t+\t50\tAAAAA\n')
    run(TAFFY, 'index', '-i', maf_path)

    out, _ = run(TAFFY, 'view', '-i', maf_path, '-m', '-r', 'refA.chrA:0-100')
    assert _count_alignment_blocks(out) == 2

    # exercise the partial-clip path on a sub-range to test gap-refill
    out, _ = run(TAFFY, 'view', '-i', maf_path, '-m', '-r', 'refA.chrA:2-7')
    assert _count_alignment_blocks(out) == 2


def test_view_first_block_past_zero_with_prev_contig(tmp_path):
    """Bug 3 coverage: queries on a contig whose first block starts past 0,
    with a lexicographically-smaller contig preceding it in the index.
    Both with and without the Bug 3 fix this passes; the test exists to
    keep the code path covered."""
    maf_path = str(tmp_path / 'bug3.maf')
    with open(maf_path, 'w') as f:
        f.write('##maf version=1 scoring=N/A\n\n')
        f.write('a\ns\ta.chrA\t0\t5\t+\t1000\tACGTA\ns\tother.chrZ\t0\t5\t+\t100\tACGTA\n\n')
        f.write('a\ns\tz.chrZ\t500\t10\t+\t2000\tACGTACGTAC\ns\tother.chrZ\t5\t10\t+\t100\tACGTACGTAC\n\n')
        f.write('a\ns\tz.chrZ\t1500\t10\t+\t2000\tACGTACGTAC\ns\tother.chrZ\t15\t10\t+\t100\tACGTACGTAC\n')
    run(TAFFY, 'index', '-i', maf_path, '-b', '1')

    out, _ = run(TAFFY, 'view', '-i', maf_path, '-m', '-r', 'z.chrZ:0-700')
    assert _count_alignment_blocks(out) == 1
    assert 'z.chrZ\t500' in out

    out, _ = run(TAFFY, 'view', '-i', maf_path, '-m', '-r', 'z.chrZ:0-2000')
    assert _count_alignment_blocks(out) == 2

    out, _ = run(TAFFY, 'view', '-i', maf_path, '-m', '-r', 'z.chrZ:0-100')
    assert _count_alignment_blocks(out) == 0
