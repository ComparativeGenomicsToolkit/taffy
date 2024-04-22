from taffy.lib import AlignmentReader, get_column_iterator, get_window_iterator
import torch
import math


class TorchDatasetColumnIterator(torch.utils.data.IterableDataset):
    """ PyTorch column or window iterator for indexed taf/maf file.
    """
    def __init__(self, alignment_file, taf_index=None, sequence_intervals=None,
                 window_length=1, step=1, include_sequence_names=True,
                 include_non_ref_columns=True):
        """

        :param alignment_file: Taf or maf file
        :param taf_index: Index for alignment file
        :param sequence_intervals: A sequence of reference sequence intervals, each a tuple of (sequence_name (string),
        start (int), length (int))
        :param window_length: Length of window of columns, must be > 0
        :param step: Number of bases between the start of each successive window, step <= window_length
        :param include_sequence_names:  Bool, include sequence names in the columns
        :param include_non_ref_columns: Include non-reference columns in the iteration
        """
        super(TorchDatasetColumnIterator).__init__()
        self.alignment_file = alignment_file
        self.taf_index = taf_index
        self.sequence_intervals = sequence_intervals
        self.window_length = window_length
        self.step = step
        self.include_sequence_names = include_sequence_names
        self.include_non_ref_columns = include_non_ref_columns

    def __iter__(self):
        worker_info = torch.utils.data.get_worker_info()
        if worker_info is None or self.taf_index is None:  # If not doing anything in parallel or no taf index
            # Make the alignment reader
            alignment_reader = AlignmentReader(file=self.alignment_file,
                                               taf_index=self.taf_index,
                                               sequence_intervals=self.sequence_intervals)
        else:
            total_length = sum([length for (name, start, length) in self.sequence_intervals])
            per_worker_length = int(math.ceil(total_length / float(worker_info.num_workers)))
            assert per_worker_length > 0
            worker_id = worker_info.id

            # Get set of intervals for the worker
            begin = worker_id * per_worker_length
            end = begin + per_worker_length
            offset = 0
            sub_sequence_intervals = []
            for seq, start, length in self.sequence_intervals:
                if begin < offset + length:  # Interval ends at or after desired beginning
                    # If the offset (left end of the interval) starts before begin, update the start coordinate
                    # and trim the length (we store this in mod_length so as not to alter the original length, which
                    # we need to track the total interval length with respect to the end
                    if offset < begin:
                        start += begin - offset
                        mod_length -= begin - offset
                        assert mod_length >= 0  # This should be true because begin - offset < length
                    else:
                        mod_length = length

                    if end <= offset + length:  # The end of the interval is beyond the end of the desired interval
                        # Trim the end of the interval
                        mod_length -= offset + length - end
                        assert mod_length >= 0  # This similarly must be true assuming begin <= end
                        sub_sequence_intervals.append((seq, start, mod_length))  # Add modified interval to set
                        break  # Break at the point because no further intervals can overlap this one
                    else:
                        sub_sequence_intervals.append((seq, start, mod_length))  # End is beyond this interval, so add
                        # to list and keep going
                offset += length

            alignment_reader = AlignmentReader(file=self.alignment_file,
                                               taf_index=self.taf_index,
                                               sequence_intervals=sub_sequence_intervals)

        #  Return the column or window iterator
        if self.window_length > 1:
            return get_window_iterator(alignment_reader,
                                       window_length=self.window_length,
                                       step=self.step,
                                       include_sequence_names=self.include_sequence_names,
                                       include_non_ref_columns=self.include_non_ref_columns)
        else:
            return get_column_iterator(alignment_reader,
                                       include_sequence_names=self.include_sequence_names,
                                       include_non_ref_columns=self.include_non_ref_columns)
