# Taffy Scripts

All Taffy scripts are in the ./scripts directory.

## Taffy-Alignment-Plot

For example, display an alignment:

    ./scripts/taffy_alignment_plot.py --max_sequences 6 ./tests/evolverMammals.maf --out_file ./out.pdf --bin_number 500 --sample_every_nth_block=5 --show_dot_plot --show_secondary_alignments --show_identity --show_sequence_boundaries --show_copy_number

Runs in (M1 Mac):

    8.04s user 2.35s system 199% cpu 5.207 total

And produces:

![image](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/examples/alignment_plot_example.png)

These options select six sequences stratified by alignment coverage to display and show coverage, identity, copy-number and dot-plots for each. Note the sample_every_nth_block param can be used to speed up displaying larger alignments by only sampling a subset of blocks. Surprisingly, this has little to no effect on visualizations at megabase scales.

You can also opt for a synteny plot, for example:

    time ./scripts/taffy_alignment_plot.py --max_sequences 3 --sample_highest_coverage ./tests/evolverMammals.maf --out_file ./out.png --bin_number 500 --sample_every_nth_block=5 --hide_coverage --show_synteny_plot --show_secondary_alignments --show_sequence_boundaries --output_format png

Runs in (M1 Mac):

    8.10s user 2.25s system 182% cpu 5.684 total

And produces:

![image](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/examples/synteny_example.png)

This works best for relatively closely related sequences.

TODO - documenting tree manipulation scripts
