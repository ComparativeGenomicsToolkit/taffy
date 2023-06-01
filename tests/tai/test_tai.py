#!/usr/bin/env python3

# test taf index by pulling out some regions frome evolver mammals.
#
# todo: should be added to some kind of proper framework, but will get the job done
#       in the meantime.
#
# MUST BE RUN FROM taf/, ie python tests/tai/tai_test.py

import sys, os, subprocess

def maf_compare(maf_path_1, maf_path_2):        
    assert os.path.isfile(maf_path_1)
    assert os.path.isfile(maf_path_2)

    # strip out whitespace
    clean_cmd = 'sed -e \'s/\\t/ /g\' | sed -e \'s/  */ /g\' | grep -v ^# | sed -r \'/^\\s*$/d\' | grep -v ^a'
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
        cmd = ['taffy', 'view', '-i', taf_path, '-n', name_map_path, '-o', renamed_taf_path]
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

    # rename the ancesstors
    mapping_path = './tests/name-mapping.tsv'
    with open(mapping_path, 'w') as mapping_file:
        mapping_file.write('\nAnc0\tAncestor0\n')
        mapping_file.write('Anc1\tAncestor1\n')
        mapping_file.write('Anc2\tAncestor2\n')
    renamed_taf_path = taf_path + '.renamed'        
    subprocess.check_call(['./bin/taffy', 'view', '-i', taf_path, '-o', renamed_taf_path, '-n', mapping_path])

                     
    
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


subprocess.check_call(['rm', '-f', taf_path, taf_rle_path, maf_path])
