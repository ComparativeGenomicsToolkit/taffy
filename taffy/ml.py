from lib import AlignmentReader, get_column_iterator
import torch
import math


class TorchDatasetColumnIterator(torch.utils.data.IterableDataset):
    """ PyTorch column iterator for taf file """
    def __init__(self, taf_file, taf_index=None, sequence_name=None, start=-1, length=-1):
        self.taf_file = taf_file
        self.taf_index = taf_index
        self.sequence_name = sequence_name
        self.start = start
        self.length = length

    def __iter__(self):
        # Make the alignment reader
        alignment_reader = AlignmentReader(file=self.taf_file,
                                           use_run_length_encoding=False,
                                           taf_index=self.taf_index,
                                           sequence_name=self.sequence_name,
                                           start=self.start,
                                           length=self.length)

        #  Return the column iterator
        return get_column_iterator(alignment_reader)

    # Define a `worker_init_fn` that configures each dataset copy differently
    @staticmethod
    def worker_init_fn(worker_id):
        worker_info = torch.utils.data.get_worker_info()
        dataset = worker_info.dataset  # the dataset copy in this worker process
        overall_start = dataset.start
        overall_length = dataset.length
        # configure the dataset to only process the split workload
        per_worker_length = int(math.ceil(overall_length / float(worker_info.num_workers)))
        worker_id = worker_info.id
        dataset.start = overall_start + worker_id * per_worker_length
        dataset.end = min(dataset.start + per_worker_length, overall_start + overall_length)
