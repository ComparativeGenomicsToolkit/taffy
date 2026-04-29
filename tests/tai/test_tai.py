#!/usr/bin/env python3

# test taf index by pulling out some regions frome evolver mammals.
#
# todo: should be added to some kind of proper framework, but will get the job done
#       in the meantime.
#
# MUST BE RUN FROM taf/, ie python tests/tai/tai_test.py

import sys, os, subprocess

def maf_compare(maf_path_1, maf_path_2, filter_empty=True):        
    assert os.path.isfile(maf_path_1)
    assert os.path.isfile(maf_path_2)

    # strip out whitespace
    clean_cmd = 'sed -e \'s/\\t/ /g\' | sed -e \'s/  */ /g\' | grep -v ^# | sed -r \'/^\\s*$/d\' | grep -v ^a'
    if filter_empty:
        clean_cmd += ' | awk \'$4 != 0 {print $0}\''
    subprocess.check_call('cat {} | {} > {}.clean'.format(maf_path_1, clean_cmd, maf_path_1), shell=True)
    subprocess.check_call('cat {} | {} > {}.clean'.format(maf_path_2, clean_cmd, maf_path_2), shell=True)

    # then run diff
    subprocess.check_call(['diff', maf_path_1 + '.clean', maf_path_2 + '.clean'])

    subprocess.check_call(['rm', '-f', '{}.clean'.format(maf_path_1)])
    subprocess.check_call(['rm', '-f', '{}.clean'.format(maf_path_2)])

        
def test_region(taf_path, contig, start, end, rev_name_map_path=None):
    assert os.path.isfile(maf_path)

    # extract the region right into maf
    out_path = './tests/tai/{}_{}_{}.taf.maf'.format(os.path.basename(os.path.splitext(maf_path)[0]), start, end)
    cmd = './bin/taffy view -i {} -r {}:{}-{} -m'.format(taf_path, contig, start, end)
    if rev_name_map_path:
        cmd += ' -n {}'.format(rev_name_map_path)
    cmd += '  > {}'.format(out_path)
    
    subprocess.check_call(cmd, shell=True)
    assert os.path.isfile(out_path)

    # run the comparison.
    # this maf was made with ./evolverMammalsMakeSubregions.py
    # note that we can run this query even on the renamed file because the rev_name_map contains both directions for Anc0!
    truth_path = './tests/tai/{}_{}_{}.maf'.format(os.path.basename(os.path.splitext(maf_path)[0]), start, end)
    assert os.path.isfile(truth_path)
    # sanity check to make sure we're not comparing empty files
    assert os.path.getsize(truth_path) > 100
    maf_compare(out_path, truth_path)

    subprocess.check_call(['rm', '-f', out_path])

    # repeat for taf
    out_path = './tests/tai/{}_{}_{}.out.taf'.format(os.path.basename(os.path.splitext(maf_path)[0]), start, end)
    cmd = './bin/taffy view -i {} -r {}:{}-{} -o {}'.format(taf_path, contig, start, end, out_path)
    if rev_name_map_path:
        cmd += ' -n {}'.format(rev_name_map_path)
    subprocess.check_call(cmd, shell=True)
    assert os.path.isfile(out_path)

    truth_path = './tests/tai/{}_{}_{}.taf'.format(os.path.basename(os.path.splitext(maf_path)[0]), start, end)
    assert os.path.getsize(truth_path) > 100
    subprocess.check_call(['diff', out_path, truth_path])

    subprocess.check_call(['rm', '-f', out_path])
    
    
def create_index(taf_path, block_size):
    assert os.path.isfile(taf_path)
    subprocess.check_call(['rm', '-f', taf_path + '.tai'])

    subprocess.check_call(['./bin/taffy', 'index', '-i', taf_path, '-b', str(block_size)])
    assert os.path.isfile(taf_path + '.tai')

def check_anc0_stats(stats_string, renamed=False):
    true_string = '''Anc0.Anc0refChr0	4151
Anc0.Anc0refChr10	14504
Anc0.Anc0refChr11	38002
Anc0.Anc0refChr1	3407
Anc0.Anc0refChr2	269145
Anc0.Anc0refChr3	165
Anc0.Anc0refChr4	13557
Anc0.Anc0refChr5	50896
Anc0.Anc0refChr6	22717
Anc0.Anc0refChr7	1851
Anc0.Anc0refChr8	111467
Anc0.Anc0refChr9	4824'''
    if renamed:
        true_string = true_string.replace('Anc0.Anc0', 'Ancestor0.Anc0')
    stats_string = '\n'.join(sorted(stats_string.strip().split('\n')))
    true_string = '\n'.join(sorted(true_string.strip().split('\n')))
    if stats_string != true_string:
        sys.stderr.write('\n    found stats:\n{}\n    different from true stats:\n{}\n'.format(stats_string, true_string))
    assert(stats_string == true_string)
    
def test_tai(regions_path, taf_path, bgzip, block_size, name_map_path=None, rev_name_map_path=None):
    sys.stderr.write(" * running indexing/extraction tests on {} with bzgip={} and blocksize={} and rename={}".format(taf_path, bgzip, block_size, name_map_path != None))
    sys.stderr.flush()
    assert (name_map_path == None) == (rev_name_map_path == None)
    
    if bgzip:
        gz_path = taf_path + '.gz'
        subprocess.check_call('bgzip -c {} > {}'.format(taf_path, gz_path), shell=True)
        taf_path = gz_path

    if name_map_path:
        # run the name conversion
        renamed_taf_path = taf_path + '.renamed'
        if bgzip:
            renamed_taf_path += '.gz'
        cmd = ['./bin/taffy', 'view', '-i', taf_path, '-n', name_map_path, '-o', renamed_taf_path]
        if bgzip:
            cmd += ['-c']
        subprocess.check_call(cmd)
        pre_renamed_path = taf_path
        taf_path = renamed_taf_path
    
    create_index(taf_path, block_size)

    with open(regions_path, 'r') as regions_file:
        for line in regions_file:
            contig, start, end = line.split()[:3]
            test_region(taf_path, contig, start, end, rev_name_map_path=rev_name_map_path)

    if taf_path.endswith('.taf') or taf_path.endswith('.taf.gz'):
        seq_stats = subprocess.check_output(['./bin/taffy', 'stats', '-s', '-i', taf_path]).decode('utf-8')
        check_anc0_stats(seq_stats, renamed=name_map_path != None)

    if bgzip:
        subprocess.check_call(['rm', '-f', taf_path])
    subprocess.check_call(['rm', '-f', taf_path + '.tai'])
    if name_map_path:
        subprocess.check_call(['rm', '-f', pre_renamed_path])
    sys.stderr.write("\t\t\tOK\n")

def test_tai_naming(regions_path, taf_path):
    sys.stderr.write(" * running name mapping tests on {}".format(taf_path))
    sys.stderr.flush()

    # rename the ancesstors
    mapping_path = './tests/name-mapping.tsv'
    with open(mapping_path, 'w') as mapping_file:
        mapping_file.write('\nAnc0\tAncestor0\n')
        mapping_file.write('Anc1\tAncestor1\n')
        mapping_file.write('Anc2\tAncestor2\n')
    renamed_taf_path = taf_path + '.renamed'        
    subprocess.check_call(['./bin/taffy', 'view', '-i', taf_path, '-o', renamed_taf_path, '-n', mapping_path])


def test_tai_taf_1bp_extractions(taf_path, maf_path, step):
    """ do a bunch of 1bp queries to make sure taf is same as maf """
    chr_name = 'Anc0.Anc0refChr0'
    chr_length = 4151
    sys.stderr.write(" * running TAF/MAF indexing/extraction comparison tests on {} with step {}".format(chr_name, step))
    sys.stderr.flush()

    create_index(taf_path, 10000)
    create_index(maf_path, 10000)
    for pos in range(0, chr_length, step):
        out_maf_path = './tests/tai/test_{}_{}.maf.taf'.format(chr_name, chr_length)
        out_taf_path = './tests/tai/test_{}_{}.taf.taf'.format(chr_name, chr_length)        
        subprocess.check_call(['./bin/taffy', 'view', '-i', maf_path, '-o', out_maf_path, '-r', '{}:{}'.format(chr_name, pos)])
        subprocess.check_call(['./bin/taffy', 'view', '-i', taf_path, '-o', out_taf_path, '-r', '{}:{}'.format(chr_name, pos)])
        subprocess.check_call(['diff', out_maf_path, out_taf_path])
        subprocess.check_call(['rm', '-f', out_maf_path, out_taf_path])
                     
    sys.stderr.write("\t\t\tOK\n")


def _count_alignment_blocks(maf_text):
    return sum(1 for line in maf_text.splitlines() if line.startswith('a'))


def test_view_missing_region_does_not_error():
    """Bug 1: querying a region not present in the index emits a header-only MAF and exits 0.

    Previously taffy view would exit 1 with "Region X not found in taffy index", which forced
    batch callers (e.g. cactus-phast's chunker) to wrap stderr-pattern matching around every call.
    """
    sys.stderr.write(" * Bug 1 regression: missing region returns header-only MAF, exit 0")
    sys.stderr.flush()
    maf_path = './tests/tai/test_bug1_missing.maf'
    with open(maf_path, 'w') as f:
        f.write('##maf version=1 scoring=N/A\n\n')
        f.write('a\ns\trefA.chrA\t0\t5\t+\t100\tACGTA\ns\tother.chrZ\t0\t5\t+\t50\tACGTA\n')
    subprocess.check_call(['./bin/taffy', 'index', '-i', maf_path])
    out_path = maf_path + '.out'
    # contig that is not in the index at all
    rc = subprocess.call('./bin/taffy view -i {} -m -r NoSuchContig:0-10 > {}'.format(maf_path, out_path), shell=True)
    assert rc == 0, "expected exit 0 for missing contig, got {}".format(rc)
    with open(out_path) as f:
        out = f.read()
    assert _count_alignment_blocks(out) == 0, "expected no alignment blocks, got: {!r}".format(out)
    assert out.startswith('##maf'), "expected MAF header, got: {!r}".format(out)
    # contig that IS in the index but the queried sub-range doesn't overlap any block
    rc = subprocess.call('./bin/taffy view -i {} -m -r refA.chrA:50-60 > {}'.format(maf_path, out_path), shell=True)
    assert rc == 0, "expected exit 0 for non-overlapping sub-range, got {}".format(rc)
    with open(out_path) as f:
        out = f.read()
    assert _count_alignment_blocks(out) == 0, "expected no alignment blocks, got: {!r}".format(out)
    # bgzipped empty-region output must be a valid header-only .maf.gz, NOT a 0-byte file.
    # The early-return path used to skip LW_destruct, leaving the bgzip writer unflushed.
    gz_out = maf_path + '.out.gz'
    rc = subprocess.call(['./bin/taffy', 'view', '-i', maf_path, '-m', '-c',
                          '-r', 'NoSuchContig:0-10', '-o', gz_out])
    assert rc == 0
    assert os.path.getsize(gz_out) > 0, "bgzipped empty output should not be 0 bytes"
    decompressed = subprocess.check_output(['bgzip', '-dc', gz_out]).decode('utf-8')
    assert decompressed.startswith('##maf'), "bgzipped output should decompress to a MAF header, got: {!r}".format(decompressed)
    assert _count_alignment_blocks(decompressed) == 0
    subprocess.check_call(['rm', '-f', maf_path, maf_path + '.tai', out_path, gz_out])
    sys.stderr.write("\t\t\tOK\n")


def test_view_clip_zero_length_row():
    """Bug 2: clip_alignment must not assert when a block contains a length=0 (all-gap) row.

    Reproduction: a MAF block with an `s` line whose length field is 0 and bases are all dashes.
    Previously this aborted with `Assertion strlen(row->bases) == 0 failed` whenever the queried
    range covered the block without trimming (so neither left nor right clip ran).
    """
    sys.stderr.write(" * Bug 2 regression: clip_alignment doesn't crash on length=0 rows")
    sys.stderr.flush()
    maf_path = './tests/tai/test_bug2_zero_len_row.maf'
    with open(maf_path, 'w') as f:
        f.write('##maf version=1 scoring=N/A\n\n')
        # block 1 contains an all-gap row (length 0); query range fully covers the block so
        # neither left nor right trim runs in clip_alignment
        f.write('a\ns\trefA.chrA\t0\t5\t+\t100\tACGTA\ns\tother.chrZ\t0\t0\t+\t50\t-----\n\n')
        f.write('a\ns\trefA.chrA\t5\t5\t+\t100\tTTTTT\ns\tother.chrZ\t0\t5\t+\t50\tAAAAA\n')
    subprocess.check_call(['./bin/taffy', 'index', '-i', maf_path])
    out_path = maf_path + '.out'
    rc = subprocess.call('./bin/taffy view -i {} -m -r refA.chrA:0-100 > {}'.format(maf_path, out_path), shell=True)
    assert rc == 0, "expected exit 0 (was SIGABRT), got {}".format(rc)
    with open(out_path) as f:
        out = f.read()
    assert _count_alignment_blocks(out) == 2, "expected 2 alignment blocks, got: {!r}".format(out)
    # also exercise the partial-clip path on a sub-range to make sure the gap-refill still works
    rc = subprocess.call('./bin/taffy view -i {} -m -r refA.chrA:2-7 > {}'.format(maf_path, out_path), shell=True)
    assert rc == 0
    with open(out_path) as f:
        out = f.read()
    assert _count_alignment_blocks(out) == 2, "expected 2 clipped blocks, got: {!r}".format(out)
    subprocess.check_call(['rm', '-f', maf_path, maf_path + '.tai', out_path])
    sys.stderr.write("\t\t\tOK\n")


def test_view_first_block_past_zero_with_prev_contig():
    """Smoke test for the Bug 3 code path: contig whose first block starts past 0, with a
    lexicographically-smaller contig preceding it in the index.

    NOTE: with Bug 1's fix in place this scenario is no longer user-visibly broken in either
    direction -- the OLD `tair_1 == NULL` lookup happens to scan forward through the prev
    contig and still finds the target blocks for any query that genuinely overlaps them, and
    a query that misses everything turns into a clean header-only exit either way. So this
    test is a *coverage* test for the new fallback-search code path rather than a regression
    test that fails without the Bug 3 fix specifically. The Bug 3 fix is still worth keeping:
    it skips wasted scanning of the prev contig and avoids relying on the scan loop to walk
    past unrelated blocks before reaching the target.
    """
    sys.stderr.write(" * Bug 3 coverage: queries on contig with first block past 0 + prev contig")
    sys.stderr.flush()
    maf_path = './tests/tai/test_bug3_first_past_zero.maf'
    with open(maf_path, 'w') as f:
        f.write('##maf version=1 scoring=N/A\n\n')
        f.write('a\ns\ta.chrA\t0\t5\t+\t1000\tACGTA\ns\tother.chrZ\t0\t5\t+\t100\tACGTA\n\n')
        f.write('a\ns\tz.chrZ\t500\t10\t+\t2000\tACGTACGTAC\ns\tother.chrZ\t5\t10\t+\t100\tACGTACGTAC\n\n')
        f.write('a\ns\tz.chrZ\t1500\t10\t+\t2000\tACGTACGTAC\ns\tother.chrZ\t15\t10\t+\t100\tACGTACGTAC\n')
    # use a small index block size so each block gets its own .tai entry; this is what makes
    # the "tair_2 lookup lands on the very block we want" pattern reachable in principle
    subprocess.check_call(['./bin/taffy', 'index', '-i', maf_path, '-b', '1'])
    out_path = maf_path + '.out'
    rc = subprocess.call('./bin/taffy view -i {} -m -r z.chrZ:0-700 > {}'.format(maf_path, out_path), shell=True)
    assert rc == 0
    with open(out_path) as f:
        out = f.read()
    assert _count_alignment_blocks(out) == 1, "expected 1 alignment block, got: {!r}".format(out)
    assert 'z.chrZ\t500' in out, "expected first z.chrZ block in output, got: {!r}".format(out)
    rc = subprocess.call('./bin/taffy view -i {} -m -r z.chrZ:0-2000 > {}'.format(maf_path, out_path), shell=True)
    assert rc == 0
    with open(out_path) as f:
        out = f.read()
    assert _count_alignment_blocks(out) == 2, "expected 2 alignment blocks, got: {!r}".format(out)
    rc = subprocess.call('./bin/taffy view -i {} -m -r z.chrZ:0-100 > {}'.format(maf_path, out_path), shell=True)
    assert rc == 0
    with open(out_path) as f:
        out = f.read()
    assert _count_alignment_blocks(out) == 0, "expected no alignment blocks, got: {!r}".format(out)
    subprocess.check_call(['rm', '-f', maf_path, maf_path + '.tai', out_path])
    sys.stderr.write("\tOK\n")


sys.stderr.write("Running tai tests...\n")
maf_path_in = './tests/evolverMammals.maf'
maf_path = './tests/tai/evolverMammals.maf'
subprocess.check_call(['cp', maf_path_in, maf_path])
taf_path = './tests/tai/evolverMammals.taf'
sys.stderr.write(" * creating evolver taf {}".format(taf_path))
subprocess.check_call(['./bin/taffy', 'view', '-i', maf_path, '-o', taf_path])
sys.stderr.write("\t\t\tOK\n")
taf_rle_path = './tests/tai/evolverMammals.rle.taf'
sys.stderr.write(" * creating run length encoded evolver taf {}".format(taf_rle_path))
subprocess.check_call(['./bin/taffy', 'view', '-i', maf_path, '-o', taf_rle_path, '-u'])
sys.stderr.write("\t\t\tOK\n")    
regions_path = './tests/tai/evolverMammals_subregions.bed'
name_map_path = './tests/tai/evolverMammals_ancestor_name_map.tsv'
rev_name_map_path = './tests/tai/evolverMammals_ancestor_rev_name_map.tsv'

test_tai(regions_path, taf_path, False, 111)
test_tai(regions_path, taf_path, True, 200)
test_tai(regions_path, taf_rle_path, False, 111)
test_tai(regions_path, taf_rle_path, True, 200)
test_tai(regions_path, maf_path, False, 111)
test_tai(regions_path, maf_path, True, 200)                         
test_tai(regions_path, maf_path, True, 200, name_map_path=name_map_path, rev_name_map_path=rev_name_map_path)
test_tai_taf_1bp_extractions(taf_path, maf_path, 5)

# regression tests for the three taffy view -r bugs documented in cactus-phast/taffy-bug.md
test_view_missing_region_does_not_error()
test_view_clip_zero_length_row()
test_view_first_block_past_zero_with_prev_contig()

subprocess.check_call(['rm', '-f', taf_path, taf_rle_path, maf_path])
