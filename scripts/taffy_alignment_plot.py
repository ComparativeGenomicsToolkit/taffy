#!python3

""" Simple script to visualise taffy file
"""

import matplotlib.pyplot as plt
import numpy as np
import argparse
import taffy.lib

# Todos:
# Add support to show alignment annotations
# Make it easy to integrate into a Jupyter notebook


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
        "--sample_highest_coverage",
        dest="sample_highest_coverage",
        action='store_true',
        default=False,
        help="Show the highest coverage rows instead of uniformly sampling across the coverage distribution."
    )
    parser.add_argument(
        "--region",
        type=str,
        default=None,
        help="Show a limited region defined by SEQ:START-END, where SEQ is a row-0 sequence name, "
             "and START-END are 0-based open-ended like BED. Requires a .tai index file to exist for the alignment file"
    )
    parser.add_argument(
        "--print_coverage_distribution",
        action='store_true',
        default=False,
        help="Print the coverages and identities of each sample to standard out"
    )
    parser.add_argument(
        "--sample_every_nth_block",
        type=int,
        default=1,
        help="Speed up processing large alignments by only looking at every nth block"
    )
    parser.add_argument(
        "--hide_coverage",
        action='store_true',
        default=False,
        help="Don't show coverage"
    )
    parser.add_argument(
        "--show_copy_number",
        action='store_true',
        default=False,
        help="Show the copy number (number of alignments) per bin"
    )
    parser.add_argument(
        "--show_identity",
        action='store_true',
        default=False,
        help="Show the identity of the major alignment per bin"
    )
    parser.add_argument(
        "--show_dot_plot",
        action='store_true',
        default=False,
        help="Show a dot plot of the alignment"
    )
    parser.add_argument(
        "--show_synteny_plot",
        action='store_true',
        default=False,
        help="Show a synteny of the alignment"
    )
    parser.add_argument(
        "--show_secondary_alignments",
        action='store_true',
        default=False,
        help="Show points for minor alignments per bin (shown as triangles instead of circles)"
    )
    parser.add_argument(
        "--show_sequence_boundaries",
        action='store_true',
        default=False,
        help="Show dotted lines to illustrate the boundaries between contigs. Orange lines are breaks in the reference,"
             "red lines are breaks in the non-reference sequence"
    )
    parser.add_argument(
        "--output_format",
        type=str,
        default="pdf",
        help="Output format for plot, works with pdf/png"
    )

    # Parse arguments
    args = parser.parse_args()

    # If we have a specific subregion of the alignment we wish to visualize
    if args.region:
        t = args.region.split(":")
        seq_name = t[0]
        t = t[1].split("-")
        start, end = int(t[0]), int(t[1])
        assert start <= end
        sequence_intervals = ((seq_name, start, end-start),)
        taf_index = taffy.lib.TafIndex(args.alignment_file + ".tai",
                                       is_maf=taffy.lib.TafIndex.is_maf(args.alignment_file))
    else:  # Otherwise, we do not subset
        sequence_intervals, taf_index = None, None

    def get_seq_name_for_row(row):
        """ Quick function to get the desired name format for non-reference alignment rows"""
        return row.sequence_name().split(args.species_delimiter)[0] if not args.all_sequences else row.sequence_name()

    # Do two things by reading through the alignment:
    # (1) Get the reference sequence intervals
    # (2) Build a map of non-reference sequences to alignment coverage,  allowing us to filter to just those
    # considered sequences
    p_ref_seq, p_ref_start, p_ref_length = None, 0, 0  # Variables to track the reference sequence intervals
    total_ref_gap_length = 0  # Calculate bases in blocks that are skipped by any sampling
    reference_sequence_intervals = []  # List of reference sequence intervals that we'll build
    seq_name_to_coverage = {}  # Map of total aligned bases for each non-ref row, used to choose which rows to display
    with taffy.lib.AlignmentReader(args.alignment_file, sequence_intervals=sequence_intervals,
                                   taf_index=taf_index) as ar:
        for i, block in enumerate(ar):  # For each block in alignment
            if i % args.sample_every_nth_block:  # Only considered every nth block
                total_ref_gap_length += block.first_row().length()  # Add the length of the block to the total ref
                # gap length, we correct this at the end - although this is hacky because it assumes gaps are uniformly
                # distributed length wise (which should be approximately true, but is still heuristic)
                continue  # Now skip the block!
            # First (1), get the reference sequence intervals
            ref_row = block.first_row()  # Get reference row
            seq_name = ref_row.sequence_name()  # Get the reference row sequence name
            if seq_name != p_ref_seq:  # If we have a new reference sequence
                if p_ref_seq is not None:  # If there is a previous reference sequence, add it to the list
                    reference_sequence_intervals.append((p_ref_seq, p_ref_start, p_ref_length))
                # Set the new reference sequence interval
                p_ref_seq, p_ref_start, p_ref_length = seq_name, ref_row.start(), ref_row.length()
            else:
                p_ref_length += ref_row.length()  # If not a new ref sequence, just add the number of ref bases
                # in the block

            # (2) Now build a map of non-reference sequences to alignment coverage
            for non_ref_row in list(block)[1:]:  # For each non-reference row
                seq_name = get_seq_name_for_row(non_ref_row)
                if seq_name not in seq_name_to_coverage:
                    seq_name_to_coverage[seq_name] = 0
                # Add to coverage map
                seq_name_to_coverage[seq_name] += non_ref_row.length()

        if p_ref_seq is not None:  # Don't forget the last interval
            reference_sequence_intervals.append((p_ref_seq, p_ref_start, p_ref_length))

    # Sort the non-ref aligned sequences by coverage (with the highest coverage first, then decreasing)
    ordered_seq_names = sorted([seq_name for seq_name in seq_name_to_coverage],
                               key=lambda x: -seq_name_to_coverage[seq_name])

    # Subsample which sequences to show if fewer than total number of rows requested
    if args.sample_highest_coverage:  # Pick the highest coverage args.max_sequences row
        if len(ordered_seq_names) > args.max_sequences:
            sampled_seq_names = ordered_seq_names[:args.max_sequences]
        else:
            sampled_seq_names = ordered_seq_names
    else:  # Pick rows sampling rows from high to low coverage
        step = max(1, len(ordered_seq_names)//(args.max_sequences-1) if args.max_sequences > 1 else len(ordered_seq_names))
        sampled_seq_names = ordered_seq_names[::step]
        if len(sampled_seq_names) < args.max_sequences <= len(ordered_seq_names):  # If the step did not neatly
            # get the correct number of sequences
            sampled_seq_names.append(ordered_seq_names[-1])
        else:
            while len(sampled_seq_names) > args.max_sequences:  # If we got too many sequences/species remove
                # until we have the correct number
                sampled_seq_names.pop()

    # The total length of the sequence intervals
    total_ref_length = sum([length for sequence_name, start, length in reference_sequence_intervals])

    # The size of each sequence bin
    bin_size = total_ref_length / args.bin_number

    # For each species, build a partition of the alignment rows into bins
    seq_name_to_bins = {i: [[] for j in range(args.bin_number)] for i in sampled_seq_names}  # Map from seq
    # names to the bins for a row
    sampled_seq_names_set = set(sampled_seq_names)  # Make a set to make search quicker
    ref_offset = 0  # Coordinate of the first base on the reference, normalized to zero
    with taffy.lib.AlignmentReader(args.alignment_file, sequence_intervals=sequence_intervals,
                                   taf_index=taf_index) as ar:
        for i, block in enumerate(ar):  # For each block in alignment
            if i % args.sample_every_nth_block:  # Only considered every nth block
                continue  # Skip the block!
            ref_row = block.first_row()  # Get reference row
            for non_ref_row in list(block)[1:]:  # For each non-reference row
                seq_name = get_seq_name_for_row(non_ref_row)
                if seq_name in sampled_seq_names_set:
                    # Add the alignment row to the right bin, where we add both the ref and non-ref row together for
                    # comparison
                    seq_name_to_bins[seq_name][int(ref_offset / bin_size)].append((ref_row, non_ref_row))

            ref_offset += ref_row.length()

    # Partition alignment rows in each bin into colinear intervals, each termed a run, such that two consecutive
    # rows with the same sequence, strand and in order are in the same interval
    for seq_name, bins in seq_name_to_bins.items():
        # Each bin becomes a list of runs, each run being a sequence of co-linear rows.
        for i, rows in enumerate(bins):  # For each bin
            runs = []  # The set of runs is a list, representing a partition of the rows
            bins[i] = runs  # Replace the list of rows with list of runs
            if len(rows):  # If there is at least one alignment row
                # Sort alignment rows in bin so that consecutive rows are adjacent
                rows.sort(key=lambda r: (r[1].sequence_name(), r[1].strand(), r[1].start()))
                p_r = rows[0]  # The very first row
                run = [p_r]  # Add the first row to the run
                runs.append(run)  # Add the run to the list of runs
                for c_r in rows[1:]:  # For the remaining rows, if same sequence, strand and right order
                    if p_r[1].sequence_name() == c_r[1].sequence_name() and \
                            p_r[1].strand() == c_r[1].strand() and \
                            p_r[1].start() + p_r[1].length() <= c_r[1].start():
                        run.append(c_r)  # Append to current run
                    else:  # Otherwise, start a new run
                        run = [c_r]
                        runs.append(run)
                    p_r = c_r
                # Sort in place to order runs from largest to smallest
                runs.sort(key=lambda rs: sum([r[1].length() for r in rs]), reverse=True)

    seq_name_to_matches = {}  # Map from seq names to array storing aligned bases and matches for each run

    # Calculate the number of aligned and matched bases for each run
    for seq_name, bins in seq_name_to_bins.items():
        # First get max # of runs
        max_runs = max(len(runs) for runs in bins)
        # Now make arrays for matches / mismatches
        matches = np.zeros((args.bin_number, max_runs, 2))  # Indexed bin #, run #, [aligned bases, matched bases]
        seq_name_to_matches[seq_name] = matches
        for i, runs in enumerate(bins):  # For each bin
            for j, run in enumerate(runs):
                aligned_bases, matched_bases = 0, 0
                for ref_row, non_ref_row in run:
                    for ref_base, non_ref_base in zip(ref_row.bases(), non_ref_row.bases()):
                        if ref_base != '-':
                            if non_ref_base != '-':  # Coverage
                                aligned_bases += 1
                                if non_ref_base.upper() == ref_base.upper() or non_ref_base == '*':  # Identity
                                    matched_bases += 1
                matches[i, j, 0], matches[i, j, 1] = aligned_bases, matched_bases

    # Now make the matplot lib plot
    fig, sub_plots = plt.subplots(len(sampled_seq_names), 1, layout='constrained', figsize=(12, 8))
    if len(sampled_seq_names) == 1:  # Stupid hack for when you only want one plot
        sub_plots = [sub_plots]

    # Bin coordinates gives the x-axis scale, where we account for the length of the skipped blocks, assuming
    # equally distributed among bins
    bin_coordinates = [i*((total_ref_length + total_ref_gap_length) / args.bin_number) for i in range(args.bin_number)]
    for plot_no, (sub_plot, seq_name) in enumerate(zip(sub_plots, sampled_seq_names)):
        bins = seq_name_to_bins[seq_name]
        matches = seq_name_to_matches[seq_name]
        max_alignments = max(1, max([len(runs) for runs in bins]))  # Maximum number of runs in any bin
        sub_plot.set_title(seq_name, loc='left')  # Show the sequence name on the left of the plot

        if plot_no == 0:  # If first plot, add a secondary x-axis and label at the top to indicate
            # the coordinates of the reference
            secax = sub_plot.secondary_xaxis('top')
            secax.set_xlabel('Reference Position')
            secax.set_xlim(0, total_ref_length + total_ref_gap_length)  # Set the x-axis scale for this secondary axis

        sub_plot.set_xlim(0, total_ref_length + total_ref_gap_length)  # Set the x-axis scale
        sub_plot.xaxis.tick_top()  # Try to put it at the top
        sub_plot.xaxis.set_ticklabels([])  # Should hide the x-axis tick marks, but doesn't!

        # Plot the coverage as a blue points
        if not args.hide_coverage:
            color = '#0C7BDC'
            sub_plot.plot(bin_coordinates, matches[:, 0, 0]/bin_size, color=color, marker='o',
                          linestyle="", markersize=1)
            if args.show_secondary_alignments:
                for j in range(1, max_alignments):
                    secondary_coverage = matches[:, j, 0]/bin_size
                    # Replace all occurrences of 0 with -1
                    secondary_coverage = np.where(secondary_coverage == 0, -1, secondary_coverage)
                    sub_plot.plot(bin_coordinates, secondary_coverage, color=color, marker='v',
                                  linestyle="", markersize=1)  # Show the secondary points as triangles

            # Mess with the limits of the axes
            sub_plot.set_ylim(-0.05, 1+0.05)
            sub_plot.set_ylabel('coverage', color=color)

        # Similarly, plot the identity but use a yellow-ish line, as well as a dashed line showing the median
        y_axis_offset = 0
        if args.show_identity:
            color = '#FFC20A'
            sub_plot_2 = sub_plot.twinx()
            a = matches[:, 0, 1]/(matches[:, 0, 0] + 0.0001)
            sub_plot_2.plot(bin_coordinates, a, color=color, marker='o', linestyle="",
                            markersize=1)  # Make the identity squiggle slightly opaque

            if args.show_secondary_alignments:
                for j in range(1, max_alignments):
                    secondary_identity = matches[:, j, 1]/(matches[:, j, 0] + 0.0001)
                    # Replace all coverage 0 spots with a negative value to hide
                    secondary_identity = np.where(matches[:, j, 0] == 0, -1, secondary_identity)
                    sub_plot_2.plot(bin_coordinates, secondary_identity, color=color, marker='v',
                                    linestyle="", markersize=1)  # Show the secondary points as triangles

            sub_plot_2.set_ylim(-0.05, 1+0.05)
            sub_plot_2.set_ylabel('identity', color=color)
            y_axis_offset = 40

        # Plot the number of alignments per bin on the other axis
        if args.show_copy_number:
            color = 'xkcd:dark pink'
            sub_plot_3 = sub_plot.twinx()
            sub_plot_3.plot(bin_coordinates, [len(i) for i in bins], color=color, marker='o',
                            linestyle="", markersize=1)
            sub_plot_3.set_ylim(-0.05, max_alignments+0.05)
            sub_plot_3.set_ylabel('copy number', color=color)
            sub_plot_3.spines['right'].set_position(('outward', y_axis_offset))
            y_axis_offset += 40

        if args.show_dot_plot or args.show_synteny_plot:
            # First calculate the sequence offset for each non-reference sequence, laying out the non reference
            # sequences in the order they are first traversed in the alignment
            subsequence_offsets = {}
            non_ref_offset = 0
            for i, runs in enumerate(bins):  # For each bin
                for j, run in enumerate(runs):  # For each run
                    first_non_ref_row = run[0][1]  # Get the first row
                    if first_non_ref_row.sequence_name() not in subsequence_offsets:
                        subsequence_offsets[first_non_ref_row.sequence_name()] = non_ref_offset
                        non_ref_offset += first_non_ref_row.sequence_length()

        if args.show_dot_plot:
            # Make a squashed dot plot to show the approx coordinate of the bin on the non-ref sequence(s)
            sub_plot_4 = sub_plot.twinx()
            forward_color, reverse_color = 'xkcd:blue green', 'xkcd:periwinkle'
            max_coordinate = 1
            for j in range(max_alignments if args.show_secondary_alignments else 1):
                dot_plot = np.zeros((args.bin_number, 2))
                for i, runs in enumerate(bins):
                    if j < len(runs):
                        first_non_ref_row = runs[j][0][1]
                        if first_non_ref_row.strand():
                            k = subsequence_offsets[first_non_ref_row.sequence_name()] + first_non_ref_row.start()
                            dot_plot[i, :] = k, -1
                        else:
                            k = subsequence_offsets[first_non_ref_row.sequence_name()] + \
                                first_non_ref_row.sequence_length() - first_non_ref_row.start()
                            dot_plot[i, :] = -1, k
                    else:
                        dot_plot[i, :] = -1, -1
                sub_plot_4.plot(bin_coordinates, dot_plot[:, 0], color=forward_color, marker='o' if j == 0 else 'v',
                                linestyle="", markersize=1)
                sub_plot_4.plot(bin_coordinates, dot_plot[:, 1], color=reverse_color, marker='o' if j == 0 else 'v',
                                linestyle="", markersize=1)
                max_coordinate = max(max_coordinate, np.amax(dot_plot))
            sub_plot_4.set_ylim(-0.05, max_coordinate+0.05)
            sub_plot_4.set_ylabel('non-ref position', color=forward_color)
            sub_plot_4.spines['right'].set_position(('outward', y_axis_offset))

            # Optionally, draw the sequence ends on the plot
            if args.show_sequence_boundaries:
                for non_ref_seq_name, non_ref_offset in subsequence_offsets.items():
                    if non_ref_offset != 0:
                        sub_plot_4.axhline(y=non_ref_offset, linewidth=0.25, linestyle='dashed', color='red')

        if args.show_synteny_plot:
            # Make a squashed dot plot to show the approx coordinate of the bin on the non-ref sequence(s)
            sub_plot_5 = sub_plot.twiny()
            max_non_ref_coordinate = 1  # Track the largest non-reference coordinate so that we can lay out the x-axis
            colors = plt.cm.viridis(np.linspace(0, 1, args.bin_number))  # The color map for the lines
            # (so they have a gradient)
            lines = []  # A list of lines to plot
            for j in range(max_alignments if args.show_secondary_alignments else 1):
                for i, runs in enumerate(bins):
                    if j < len(runs):
                        first_non_ref_row = runs[j][0][1]
                        if first_non_ref_row.strand():
                            k = subsequence_offsets[first_non_ref_row.sequence_name()] + first_non_ref_row.start()
                            assert k >= 0 and bin_coordinates[i] >= 0
                            lines.append((bin_coordinates[i], k, i))
                        else:
                            k = subsequence_offsets[first_non_ref_row.sequence_name()] + \
                                first_non_ref_row.sequence_length() - first_non_ref_row.start()
                            assert k >= 0 and bin_coordinates[i] >= 0
                            lines.append((bin_coordinates[i], k, i))
                        max_non_ref_coordinate = max_non_ref_coordinate if max_non_ref_coordinate > k else k

            # Get the total span needed for the x-axis
            max_coordinate = total_ref_length + total_ref_gap_length
            if max_coordinate < max_non_ref_coordinate:
                max_coordinate = max_non_ref_coordinate

            # Set the axis limits
            sub_plot_5.set_xlim(-0.01, max_coordinate + 0.01)
            sub_plot_5.set_ylim(0, 1)

            # Mess with the labels
            sub_plot_5.xaxis.tick_bottom()  # Make it so the ticks are at the bottom
            sub_plot.xaxis.tick_top()  # Hide the default x-axis, again, because otherwise the call to twiny() (I think)
            # makes it reappear
            sub_plot.xaxis.set_ticklabels([])  # Again, should hide the x-axis tick marks (but doesn't)
            # If this is the last plot, add the label
            if plot_no+1 == len(sub_plots):
                sub_plot_5.set_xlabel('Non Reference Position')
                sub_plot_5.xaxis.set_label_position('bottom')

            # Now add the lines
            for (ref_pos, non_ref_pos, i) in lines:
                sub_plot_5.axline((max_coordinate * ref_pos/(total_ref_length + total_ref_gap_length), 1),
                                  (max_coordinate * non_ref_pos/max_non_ref_coordinate, 0), color=colors[i], alpha=0.6)

            # Optionally, draw the non-reference sequence ends on the plot, in this case each is a vertical line
            if args.show_sequence_boundaries:
                for non_ref_seq_name, non_ref_offset in subsequence_offsets.items():
                    if non_ref_offset != 0:
                        sub_plot_5.axvline(x=max_coordinate * non_ref_offset/max_non_ref_coordinate,
                                           linewidth=0.5, linestyle='dashed', color='red')

        # Optionally, draw the reference sequence ends on the plot (these should be accurate with respect to bins)
        if args.show_sequence_boundaries:
            sub_plot_6 = sub_plot.twinx()
            sub_plot_6.set_xlim(0, total_ref_length + total_ref_gap_length)
            sub_plot_6.yaxis.set_visible(False)  # Hide the y-axis
            ref_offset = 0
            for ref_seq_name, start, length in reference_sequence_intervals:
                if ref_offset != 0:
                    sub_plot_6.axvline(x=(ref_offset/total_ref_length) * (total_ref_length + total_ref_gap_length),
                                       linewidth=0.5, linestyle='dashed', color='orange')
                ref_offset += length

        if args.hide_coverage:
            sub_plot.yaxis.set_visible(False)

    # If an output filename is given then save it as a PDF in that file
    if args.out_file:
        plt.savefig(args.out_file, format=args.output_format, bbox_inches="tight")
    else:  # Otherwise, show it as a canvas
        plt.show()

    if args.print_coverage_distribution:
        print("Species_name,Median_coverage,Average_coverage,Median_identity,Average_identity")
        for seq_name in sampled_seq_names:
            matches = seq_name_to_matches[seq_name]
            print(f"{seq_name}, {np.nanmedian(matches[:, 0, 0]/bin_size)}, {np.average(matches[:, 0, 0]/bin_size)}, "
                  f"{np.nanmedian(matches[:, 0, 1]/(matches[:, 0, 0] + 0.0001))}, "
                  f"{np.average(matches[:, 0, 1]/(matches[:, 0, 0] + 0.0001))}")


if __name__ == '__main__':
    main()
