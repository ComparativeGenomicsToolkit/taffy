#!/usr/bin/env python3

# mafTools needs to be installed to run this
# note that mafExtractor's stop position is inclusive, so we cut 1 off each region
# must be run from ./taf

import subprocess

with open('./tests/tai/evolverMammals_subregions.bed', 'r') as reg_file:
    for line in reg_file:
        contig, start, stop = line.split()[:3]
        with open('./tests/tai/evolverMammals_{}_{}.maf'.format(start, stop), 'w') as maf_file:
            stop = str(int(stop) - 1)
            subprocess.check_call(['mafExtractor', '-m', './tests/evolverMammals.maf', '-s', contig, '--start', start, '--stop', stop],
                                  stdout=maf_file)
            
