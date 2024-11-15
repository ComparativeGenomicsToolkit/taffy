# Taffy Scripts

All Taffy scripts are in the ./scripts directory.

## Alignment-Plot

For example, display an alignment:

    ./scripts/alignment_plot.py --max_sequences 6 ./tests/evolverMammals.maf --out_file ./out.pdf --bin_number 500 --sample_every_nth_block=5 --show_dot_plot --show_secondary_alignments --show_identity --show_sequence_boundaries --show_copy_number

Runs in (M1 Mac):

    6.47s user 2.36s system 175% cpu 5.042 total

And produces:

![image](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/examples/alignment_plot_example.png)

These options select six sequences stratified by alignment coverage to display and show coverage, identity, copy-number and dot-plots for each. Note the sample_every_nth_block param can be used to speed up displaying larger alignments by only sampling a subset of blocks. Surprisingly, this has little to no effect on visualizations at megabase scales.

TODO - documenting tree manipulation scripts
