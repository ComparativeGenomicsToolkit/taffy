from taffy.lib import AlignmentReader, get_column_iterator, get_window_iterator, TafIndex
import torch
import math


def get_subsequence_intervals(sequence_intervals, number_of_partitions, partition_index):
    """ Partition a set of sequence intervals into non-overlapping sets.
    """
    total_length = sum([length for (name, start, length) in sequence_intervals])
    partition_length = int(math.ceil(total_length / float(number_of_partitions)))
    assert partition_length > 0

    # Get set of sub-intervals for the worker
    begin = partition_index * partition_length
    end = begin + partition_length
    offset = 0
    subsequence_intervals = []
    for seq, start, length in sequence_intervals:
        if begin < offset + length:  # Interval ends at or after desired beginning
            mod_length = length  # Variable representing the size of the interval after any shrinking (see
            # below)

            # If the offset (left end of the interval) starts before begin, update the start coordinate
            # and trim the length (we store this in mod_length so as not to alter the original length, which
            # we need to track the total interval length with respect to the end
            if offset < begin:
                start += begin - offset
                mod_length -= begin - offset
                assert mod_length >= 0  # This should be true because begin - offset < length

            if end <= offset + length:  # The end of the interval is beyond the end of the desired interval
                # Trim the end of the interval
                mod_length -= offset + length - end
                assert mod_length >= 0  # This similarly must be true assuming begin <= end
                subsequence_intervals.append((seq, start, mod_length))  # Add modified interval to set
                break  # Break at the point because no further intervals can overlap this one
            else:
                subsequence_intervals.append((seq, start, mod_length))  # End is beyond this interval, so add
                # to list and keep going
        offset += length
    # We're done
    return subsequence_intervals


class TorchDatasetAlignmentIterator(torch.utils.data.IterableDataset):
    """ PyTorch column or window iterator for indexed taf/maf file.

    Returns columns (or windows of columns) and labels.
    """
    def __init__(self, alignment_file, label_conversion_function,
                 taf_index_file=None, is_maf=False, sequence_intervals=None,
                 window_length=1, step=1,
                 include_sequence_names=False,
                 include_non_ref_columns=True,
                 include_column_tags=False,
                 column_one_hot=True):
        """
        param alignment_file: Taf or maf file
        :param label_conversion_function Function to convert tuple label to appropriate value for Pytorch
        :param taf_index_file: Index file for alignment file (optional if sequence intervals is None)
        :param is_maf: Used by Taf Index
        :param sequence_intervals: A sequence of reference sequence intervals, each a tuple of (sequence_name (string),
        start (int), length (int))
        :param window_length: Length of window of columns, must be > 0
        :param step: Number of bases between the start of each successive window, step <= window_length
        :param include_sequence_names:  Bool, include sequence names in the label
        :param include_non_ref_columns: Include non-reference columns in the iteration
        :param include_column_tags: Include any column tags as a final value in the label
        """
        super(TorchDatasetAlignmentIterator).__init__()
        self.alignment_file = alignment_file
        self.label_conversion_function = label_conversion_function
        self.taf_index_file = taf_index_file  # We pass in the file, not a TafIndex object to avoid trying to serialize
        # a TafIndex object, because pickling doesn't play nice with the C parts of the object through CFI
        self.is_maf = is_maf
        self.sequence_intervals = sequence_intervals
        self.window_length = window_length
        self.step = step
        self.include_sequence_names = include_sequence_names
        self.include_non_ref_columns = include_non_ref_columns
        self.include_column_tags = include_column_tags
        self.column_one_hot = column_one_hot

        # Checks
        if self.sequence_intervals is not None:
            assert self.taf_index_file  # Must be specified if there are sequence intervals to iterate on

    def __iter__(self):
        worker_info = torch.utils.data.get_worker_info()
        if worker_info is None or self.taf_index_file is None:  # If not doing anything in parallel or no taf index
            subsequence_intervals = self.sequence_intervals
        else:
            subsequence_intervals = get_subsequence_intervals(self.sequence_intervals,
                                                              number_of_partitions=worker_info.num_workers,
                                                              partition_index=worker_info.id)

        # Make the alignment reader
        alignment_reader = AlignmentReader(file=self.alignment_file,
                                           taf_index=TafIndex(file=self.taf_index_file, is_maf=self.is_maf) if
                                           self.taf_index_file is not None else None,
                                           sequence_intervals=subsequence_intervals)
        #  Create either a column or window iterator
        if self.window_length > 1:
            it = get_window_iterator(alignment_reader, window_length=self.window_length, step=self.step,
                                     include_sequence_names=self.include_sequence_names,
                                     include_non_ref_columns=self.include_non_ref_columns,
                                     include_column_tags=self.include_column_tags,
                                     column_as_int_array=not self.column_one_hot,
                                     column_as_int_array_one_hot=self.column_one_hot)
        else:
            it = get_column_iterator(alignment_reader, include_sequence_names=self.include_sequence_names,
                                     include_non_ref_columns=self.include_non_ref_columns,
                                     include_column_tags=self.include_column_tags,
                                     column_as_int_array=not self.column_one_hot,
                                     column_as_int_array_one_hot=self.column_one_hot)
        return self._convert_it(it)

    def _convert_it(self, it):
        for column, label in it:  # Convert the column and label to appropriate values
            yield torch.from_numpy(column), self.label_conversion_function(label)


def get_phyloP_label(label, default_value=0.0):
    """ Gets the phyloP tag from the label and converts it to a float. Otherwise, if not present
    returns default_value"""
    tags = label[1]
    if "phyloP" in tags:
        return float(tags["phyloP"])
    return default_value


def get_phyloP_labels(label, default_value=0.0):
    """ Gets the phyloP tags from the label of a window of columns and converts it to a torch tensor of floats.
    Any value not present returns default_value"""
    phyloPs = torch.empty(len(label))
    for i, column_label in enumerate(label):
        tags = column_label[1]
        phyloPs[i] = float(tags["phyloP"]) if "phyloP" in tags else default_value
    return phyloPs
