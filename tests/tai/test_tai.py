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
        
        
def test_region(taf_path, contig, start, end):
    assert os.path.isfile(maf_path)

    # extract the region right into maf
    out_path = './tests/tai/{}_{}_{}.taf.maf'.format(os.path.basename(os.path.splitext(maf_path)[0]), start, end)
    subprocess.check_call('./bin/taffy view -i {} -r {}:{}-{} -m  > {}'.format(taf_path, contig, start, end, out_path), shell=True)
    assert os.path.isfile(out_path)

    # run the comparison.
    # this maf was made with ./evolverMammalsMakeSubregions.py
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

def check_anc0_stats(stats_string):
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
    assert(stats_string.strip() == true_string)
    
def test_tai(regions_path, taf_path, bgzip, block_size):
    sys.stderr.write(" * running indexing/extraction tests on {} with bzgip={} and blocksize={}".format(taf_path, bgzip, block_size))
    if bgzip:
        gz_path = taf_path + '.gz'
        subprocess.check_call('bgzip -c {} > {}'.format(taf_path, gz_path), shell=True)
        taf_path = gz_path

    create_index(taf_path, block_size)

    with open(regions_path, 'r') as regions_file:
        for line in regions_file:
            contig, start, end = line.split()[:3]
            test_region(taf_path, contig, start, end)

    seq_stats = subprocess.check_output('taffy stats -s -i {} | sort -k1'.format(taf_path), shell=True).decode('utf-8')
    check_anc0_stats(seq_stats)

    if bgzip:
        subprocess.check_call(['rm', '-f', taf_path])
    subprocess.check_call(['rm', '-f', taf_path + '.tai'])
    sys.stderr.write("\t\t\tOK\n")    
    
sys.stderr.write("Running tai tests...\n")
maf_path = './tests/evolverMammals.maf'
taf_path = './tests/tai/evolverMammals.taf'
sys.stderr.write(" * creating evolver taf {}".format(taf_path))
subprocess.check_call(['./bin/taffy', 'view', '-i', maf_path, '-o', taf_path])
sys.stderr.write("\t\t\tOK\n")
taf_rle_path = './tests/tai/evolverMammals.rle.taf'
sys.stderr.write(" * creating run length encoded evolver taf {}".format(taf_rle_path))
subprocess.check_call(['./bin/taffy', 'view', '-i', maf_path, '-o', taf_rle_path, '-u'])
sys.stderr.write("\t\t\tOK\n")    
regions_path = './tests/tai/evolverMammals_subregions.bed'

test_tai(regions_path, taf_path, False, 111)
test_tai(regions_path, taf_path, True, 200)
test_tai(regions_path, taf_rle_path, False, 111)
test_tai(regions_path, taf_rle_path, True, 200)

subprocess.check_call(['rm', '-f', taf_path, taf_rle_path])
