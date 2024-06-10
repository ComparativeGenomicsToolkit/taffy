#!python3

""" Simple script to show pairwise alignment coverage of a taffy file
"""

import matplotlib.pyplot as plt
import numpy as np
import argparse
import taffy.lib


def main():
    # Set script arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "alignment_file",
        type=str,
        help="Path to the taf or maf file"
    )
    parser.add_argument(
        "--bin_number",
        type=int,
        default=500,
        help="Number of coverage bins to show."
    )
    parser.add_argument(
        "--all_sequences",
        dest="all_sequences",
        action='store_true',
        default=False,
        help="Instead of showing one row per species, show one row per sequence"
    )
    parser.add_argument(
        "--species_delimiter",
        default=".",
        help="The species name is determined as the prefix of the sequence name up to the given delimiter character, "
             "by default the . character."
    )
    parser.add_argument(
        "--out_file",
        type=str,
        default=None,
        help="Write to the given file, otherwise will render in a window."
    )
    parser.add_argument(
        "--max_sequences",
        type=int,
        default=20,
        help="Maximum number of rows to show."
    )
    parser.add_argument(
        "--sample_highest_identity",
        dest="sample_highest_identity",
        action='store_true',
        default=False,
        help="Show the highest identity rows instead of uniformly sampling across the identity distribution."
    )
    parser.add_argument(
        "--region",
        type=str,
        default=None,
        help="Print only SEQ:START-END, where SEQ is a row-0 sequence name, and START-END are 0-based "
             "open-ended like BED. Requires a .tai index file to exist"
    )
    parser.add_argument(
        "--print_coverage_distribution",
        action='store_true',
        default=False,
        help="Print the coverages and identities of each sample to standard out"
    )

    # Parse arguments
    args = parser.parse_args()

    if args.region:
        t = args.region.split(":")
        seq_name = t[0]
        t = t[1].split("-")
        start, end = int(t[0]), int(t[1])
        assert start <= end
        sequence_intervals = ((seq_name, start, end-start),)
        taf_index = taffy.lib.TafIndex(args.alignment_file + ".tai",
                                       is_maf=taffy.lib.TafIndex.is_maf(args.alignment_file))
    else:
        sequence_intervals, taf_index = None, None

    # Get the reference sequence intervals
    with taffy.lib.AlignmentReader(args.alignment_file, sequence_intervals=sequence_intervals,
                                   taf_index=taf_index) as ar:
        # Now get the intervals
        reference_sequence_intervals = list(taffy.lib.get_reference_sequence_intervals(ar))

    # The total length of the sequence intervals
    total_ref_length = sum([length for sequence_name, start, length in reference_sequence_intervals])

    # The size of each sequence bin
    bin_size = total_ref_length / args.bin_number

    # Build the coverage vectors for each species
    coverage_vectors = {}
    ref_offset = 0
    with taffy.lib.AlignmentReader(args.alignment_file, sequence_intervals=sequence_intervals,
                                   taf_index=taf_index) as ar:
        for block in ar:  # For each block in alignment
            ref_row = block.first_row()  # Get reference row
            ref_bases = ref_row.bases()
            for non_ref_row in list(block)[1:]:  # For each non-reference row
                seq_name = non_ref_row.sequence_name()
                if not args.all_sequences:  # Break the sequence name to reflect just the species name
                    seq_name = seq_name.split(args.species_delimiter)[0]

                # Get the coverage vector
                if seq_name in coverage_vectors:
                    coverage_vector = coverage_vectors[seq_name]
                else:
                    coverage_vector = np.zeros((2, args.bin_number))
                    coverage_vectors[seq_name] = coverage_vector

                # Iterate over alignment adding to coverage
                i = ref_offset
                for ref_base, non_ref_base in zip(ref_bases, non_ref_row.bases()):
                    if ref_base != '-':
                        if non_ref_base != '-':  # Coverage
                            j = int(i / bin_size)
                            coverage_vector[0, j] += 1
                            if non_ref_base.upper() == ref_base.upper() or non_ref_base == '*':  # Identity
                                coverage_vector[1, j] += 1
                        i += 1

            ref_offset += ref_row.length()

    # Normalize the coverage vectors by the bin length
    for seq_name in coverage_vectors:
        coverage_vectors[seq_name] /= bin_size

    # Sort the coverage plots by median identity (highest identity first, then decreasing)
    coverage_vectors = sorted([(sequence_name, coverage_vectors[sequence_name]) for sequence_name in coverage_vectors],
                              key=lambda x: -np.nanmedian(x[1][1, :]))

    # Subsample coverage vectors if fewer than total number of rows requested
    if args.sample_highest_identity:  # Pick the highest identity args.max_sequences row
        if len(coverage_vectors) > args.max_sequences:
            sampled_coverage_vectors = coverage_vectors[:args.max_sequences]
    else:  # Pick rows sampling rows from high to low identity
        step = max(1, len(coverage_vectors)//(args.max_sequences-1) if args.max_sequences > 1 else len(coverage_vectors))
        sampled_coverage_vectors = coverage_vectors[::step]
        #
        if len(sampled_coverage_vectors) < args.max_sequences <= len(coverage_vectors):  # If the step did not neatly
            # get the correct number of sequences
            sampled_coverage_vectors.append(coverage_vectors[-1])
        else:
            while len(sampled_coverage_vectors) > args.max_sequences:
                sampled_coverage_vectors.pop()

    # Now make the matplot lib plot
    fig, sub_plots = plt.subplots(len(sampled_coverage_vectors), 1, layout='constrained', figsize=(12, 8))
    if len(sampled_coverage_vectors) == 1:  # Stupid hack for when you only want one plot
        sub_plots = [sub_plots]

    bin_coordinates = [i*bin_size for i in range(args.bin_number)]
    for sub_plot, (species_name, coverage_vector) in zip(sub_plots, sampled_coverage_vectors):
        sub_plot.plot(bin_coordinates, coverage_vector[0, :], color='#0C7BDC', marker='o', linestyle="", markersize=1)
        sub_plot.axhline(y=np.nanmedian(coverage_vector[0, :]), linewidth=1, linestyle='dashed', color='grey')
        sub_plot.plot(bin_coordinates, coverage_vector[1, :], color='#FFC20A', marker='o', linestyle="",
                      markersize=1, alpha=0.6)  # Make the identity squiggle slightly opaque
        sub_plot.axhline(y=np.nanmedian(coverage_vector[1, :]), linewidth=1, linestyle=(0, (1, 1)), color='grey')
        sub_plot.set_xlim(0, total_ref_length)
        sub_plot.set_ylim(-0.05, min(1, max(coverage_vector[0, :]))+0.05)
        sub_plot.set_title(species_name, x=0.2, y=0.6)

        # Draw the sequence ends on the plot
        ref_offset = 0
        for sequence_name, start, length in reference_sequence_intervals:
            if ref_offset != 0:
                sub_plot.axvline(x=ref_offset, linewidth=0.5, linestyle='dashed', color='grey')
            ref_offset += length

    # Turn off the x-axis tick labels for everything except the last plot
    for sub_plot in sub_plots[:-1]:
        sub_plot.set_xticklabels([])

    # Label the bottom most axis
    sub_plots[-1].set_xlabel('Reference Position')
    sub_plots[-1].legend(("coverage", "median coverage", "identity", "median identity"))

    if args.out_file:
        plt.savefig(args.out_file, format="pdf", bbox_inches="tight")
    else:
        plt.show()

    if args.print_coverage_distribution:
        print("Species_name,Median_coverage,Average_coverage,Median_identity,Average_identity")
        for (species_name, coverage_vector) in coverage_vectors:
            print(f"{species_name}, {np.nanmedian(coverage_vector[0, :])}, {np.average(coverage_vector[0, :])}, "
                  f"{np.nanmedian(coverage_vector[1, :])}, {np.average(coverage_vector[1, :])}")


if __name__ == '__main__':
    main()
