#!/usr/bin/env python3

# mafTools needs to be installed to run this
# note that mafExtractor's stop position is inclusive, so we cut 1 off each region
# must be run from ./taf

# note: the mafExtractor output is used directly in tests, but it removes empty rows.  So .taf files are created as well -- make sure
#       there are no bugs in taffy when running this, then!
import subprocess

subprocess.check_call(['bin/taffy', 'index', '-i', './tests/evolverMammals.maf'])

with open('./tests/tai/evolverMammals_subregions.bed', 'r') as reg_file:
    for line in reg_file:
        contig, start, stop = line.split()[:3]
        maf_path = './tests/tai/evolverMammals_{}_{}.maf'.format(start, stop)
        taf_path = './tests/tai/evolverMammals_{}_{}.taf'.format(start, stop)        
        with open(maf_path, 'w') as maf_file, open(taf_path, 'w') as taf_file:
            stop = str(int(stop) - 1)
            subprocess.check_call(['mafExtractor', '-m', './tests/evolverMammals.maf', '-s', contig, '--start', start, '--stop', stop],
                                  stdout=maf_file)
            subprocess.check_call(['bin/taffy', 'view', '-i', './tests/evolverMammals.maf', '-r',
                                   '{}:{}-{}'.format(contig, start, int(stop)+1)], stdout=taf_file)
            
