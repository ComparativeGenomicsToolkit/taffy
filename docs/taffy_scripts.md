# Taffy Scripts

All Taffy scripts are in the ./scripts directory.

## Alignment-Plot

For example, display an alignment:

    ./scripts/alignment_plot.py --max_sequences 6 ./tests/evolverMammals.maf --out_file ./out.pdf --bin_number 500 --sample_every_nth_block=5 --show_dot_plot --show_secondary_alignments --show_identity --show_sequence_boundaries --show_copy_number

Runs in:

    6.47s user 2.36s system 175% cpu 5.042 total

And produces:

![image](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/examples/alignment_plot_example.png)

Showing the TAF file, selecting six sequences to display and showing covrage, identity, copy-number and dot-plots for each.

TODO - documenting tree manipulation scripts
