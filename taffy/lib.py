from taffy._taffy_cffi import ffi, lib
import numpy as np
import torch
from collections import deque

def _to_py_string(s):
    """ Convert a cffi string into a Python string in Python3 """
    return ffi.string(s).decode("utf-8")


def _to_c_string(s):
    """ Convert a Python string into a C string """
    return bytes(s, "utf-8")


def _c_tags_to_dictionary(c_tag):
    """ Convert tags to a Python dictionary """
    tags = {}
    while c_tag != ffi.NULL:
        tags[_to_py_string(c_tag.key)] = _to_py_string(c_tag.value)
        c_tag = c_tag.n_tag
    return tags


def _dictionary_to_c_tags(tags):
    """ Convert a Python dictionary to c tags """
    first_c_tag, p_c_tag = ffi.NULL, ffi.NULL
    for key in tags:
        c_tag = lib.tag_construct(_to_c_string(str(key)), _to_c_string(str(tags[key])), ffi.NULL)
        if first_c_tag == ffi.NULL:
            first_c_tag = c_tag
        else:
            p_c_tag.n_tag = c_tag
        p_c_tag = c_tag
    return first_c_tag


def _get_c_file_handle(file_string_or_handle, modifier_string="r"):
    """ Gets the c file handle for file, which can be either a file string
     or a file handle. If file handle you can set the modifier string.
     Note the Python file handle is *way* slower
     """
    return lib.fopen(_to_c_string(file_string_or_handle), _to_c_string(modifier_string)) if \
        isinstance(file_string_or_handle, str) else ffi.cast("FILE *", file_string_or_handle)


class Alignment:
    """ Represents an alignment block. See taf.h """

    def __init__(self, c_alignment=None, py_row=None):
        self._c_alignment = c_alignment
        self._py_row = py_row

    def row_number(self):
        """ Number of rows in the alignment block """
        return self._c_alignment.row_number

    def column_number(self):
        """ Number of columns in the alignment block """
        return self._c_alignment.column_number

    def first_row(self):
        """ The first row in the alignment block """
        return self._py_row

    def column_tags(self, column_index):
        """ The tags for the given column index, represented as a dictionary of strings """
        assert 0 <= column_index < self.column_number()  # Check column index is valid
        return _c_tags_to_dictionary(self._c_alignment.column_tags[column_index]) \
            if self._c_alignment.column_tags != ffi.NULL else {}

    def set_column_tags(self, column_index, tags):
        """ Set the tags for a given column index using a dictionary of tags """
        assert 0 <= column_index < self.column_number()  # Check column index is valid
        lib.tag_destruct(self._c_alignment.column_tags[column_index])  # Clean up the old tags
        self._c_alignment.column_tags[column_index] = _dictionary_to_c_tags(tags)

    def get_column(self, column_index):
        """ Get a column of the alignment as a string. Use negative coordinates to get columns
         from the end of the block """
        column_index = column_index if column_index >= 0 else (self.column_number() + column_index)  # Correct if
        # requesting a column from the end of the alignment
        assert 0 <= column_index < self.column_number()
        column = lib.alignment_get_column(self._c_alignment, column_index)  # Get the column
        column_string = _to_py_string(column)  # Convert to Python string
        lib.free(column)  # Free C string
        return column_string

    def get_column_as_np_array(self, column_index):
        """ Get a column of the alignment as a numpy int32 array, where we map the bases to
        consecutive integers, e.g. A/a=0, C/c=1, G/g=2, T/t=3, -=4, everything else=5.

        Use negative coordinates to get columns from the end of the block """
        column_index = column_index if column_index >= 0 else (self.column_number() + column_index)  # Correct if
        # requesting a column from the end of the alignment
        assert 0 <= column_index < self.column_number()
        column = lib.alignment_get_column_as_int_array(self._c_alignment, column_index)  # Get the column
        column_length = self.row_number()
        # Convert the C array to a NumPy array, copying it in process
        column_np = np.copy(np.frombuffer(ffi.buffer(column, ffi.sizeof("int32_t") * column_length), dtype=np.int32))
        lib.free(column)  # Free C string
        return column_np

    def get_column_as_np_array_one_hot(self, column_index):
        """ As get_column_as_np_array, but encoded one hot, using float32 """
        column_np = self.get_column_as_np_array(column_index)
        column_np_one_hot = np.zeros(shape=(len(column_np),6), dtype=np.float32)
        column_np_one_hot[np.arange(len(column_np)), column_np] = 1.0
        return column_np_one_hot

    def get_column_sequences(self):
        """ Get the names of the sequences in the alignment in order as a list """
        row = self.first_row()
        sequence_names = [None]*self.row_number()
        for i in range(self.row_number()):
            sequence_names[i] = row.sequence_name()
            row = row.next_row()
        return sequence_names

    def __del__(self):
        lib.alignment_destruct(self._c_alignment, 0)  # Cleans up the underlying C alignment structure

    def __str__(self):
        return _to_py_string(lib.alignment_to_string(self._c_alignment))

    def __iter__(self):
        # Make a custom row iterator object to iterate over the rows in the
        class RowIter:
            def __init__(self, row):
                self.row = row

            def __next__(self):
                if self.row is not None:
                    x = self.row
                    self.row = self.row.next_row()
                    return x
                else:
                    raise StopIteration

        return RowIter(self.first_row())


class Row:
    """ Represents a row of an alignment block. See taf.h """

    def __init__(self, c_row=None, n_row=None):
        self._c_row = c_row  # The underlying C row
        self._n_row = n_row  # The next row in the sequence of rows
        # Note by choice we do not store pointers to the left and right rows

    def sequence_name(self):
        """ The name of the sequence for the row """
        return _to_py_string(self._c_row.sequence_name)

    def start(self):
        """ The start coordinate of the base in the row (if strand is False then will be with respect to
        reverse complement sequence """
        return self._c_row.start

    def length(self):
        """ The number of bases in the row of the alignment """
        return self._c_row.length

    def sequence_length(self):
        """ The length of the underlying sequence """
        return self._c_row.sequence_length

    def strand(self):
        """ The strand (boolean) of the sequence """
        return self._c_row.strand

    def bases(self):
        """ The alignment of the row, consisting of the sequence and gap characters """
        return _to_py_string(self._c_row.bases)

    def left_gap_sequence(self):
        """ The sequence of any unaligned bases in this row between this block and the previous
        (left) alignment block """
        return "" if self._c_row.left_gap_sequence == ffi.NULL else _to_py_string(self._c_row.self.left_gap_sequence)

    def next_row(self):
        """ Get the next row in the alignment block or None if last row """
        return self._n_row

    def __del__(self):
        lib.alignment_row_destruct(self._c_row)  # Cleans up the underlying C alignment structure

    def __str__(self):
        return _to_py_string(lib.alignment_row_to_string(self._c_row))


class TafIndex:
    """ Taf Index (.tai)
    """

    def __init__(self, file, is_maf):
        """ Load from a file. Can be a file name or a Python file handle """
        c_file_handle = _get_c_file_handle(file)
        self._c_taf_index = lib.tai_load(c_file_handle, is_maf)
        if isinstance(file, str):  # Close the underlying file handle if opened
            lib.fclose(c_file_handle)

    def __del__(self):
        lib.tai_destruct(self._c_taf_index)  # Clean up the underlying C


class AlignmentReader:
    """ Taf or maf alignment parser.
    """

    def __init__(self, file, taf_index=None,  sequence_intervals=None):
        """
        :param file: File can be either a Python file handle or a string giving a path to the file.
        The underlying file can be either maf or taf. Handing in the file name is much faster as it avoids using a
        Python file object. If using with a file string remember to close the file, either the with keyword or with
        the close method. If file is compressed with zip or bgzip will automatically detect that the file is
        compressed and read it okay.
        :param taf_index: A taf index object, which is specified allows the retrieval of subranges of the alignment.
        :param sequence_intervals: A sequence of one or more tuples, each of the form (sequence_name, start, length)
        that specify the range to retrieve. Will be retrieved in order.
        """
        self.p_c_alignment = ffi.NULL  # The previous C alignment returned
        self.p_c_rows_to_py_rows = {}  # Hash from C rows to Python rows of the previous
        # alignment block, allowing linking of rows between blocks
        self.file = file
        self.c_file_handle = _get_c_file_handle(file)
        self.file_string_not_handle = isinstance(file, str)  # Will be true if the file is a string, not a file handle
        self.c_li_handle = lib.LI_construct(self.c_file_handle)
        i = lib.check_input_format(lib.LI_peek_at_next_line(self.c_li_handle))
        if i not in (0, 1):
            raise RuntimeError("Input file is not a TAF or MAF file")
        self.taf_not_maf = i == 0  # Sniff the file header to determine if a taf file
        self.taf_index = taf_index  # Store the taf index (if there is one)
        self.header_tags = self._read_header()  # Read the header tags

        # Parse run length encoding
        self.use_run_length_encoding = False
        if "run_length_encode_bases" in self.header_tags:
            assert self.header_tags["run_length_encode_bases"] in ("0", "1")
            self.use_run_length_encoding = self.header_tags["run_length_encode_bases"] == "1"

        self.sequence_intervals = sequence_intervals
        self.sequence_interval_index = 0

        if taf_index:
            assert self.sequence_intervals is not None  # Must not be None
            if len(self.sequence_intervals) > 0:
                sequence_name, start, length = sequence_intervals[0]
                assert sequence_name  # The contig name can not be none if using a taf index
                assert start >= 0  # The start coordinate must be valid
                assert length >= 0  # The length must be valid
                assert len(sequence_intervals) > 0  # Must contain at least one interval
                self._c_taf_index_it = lib.tai_iterator(taf_index._c_taf_index,
                                                        self.c_li_handle,
                                                        self.use_run_length_encoding,
                                                        _to_c_string(sequence_name), start, length)
            else:
                self._c_taf_index_it = None  # Case we have an empty iteration
        else:
            assert not sequence_intervals  # Sequence intervals should not be specified if taf index not provided

    def get_header(self):
        """ Get tags from the header line as a dictionary of key:value pairs.
        Must be called if a header is present in the file before blocks are retrieved """
        return dict(self.header_tags)  # Make a copy

    def _read_header(self):
        # Internal method to read header tags
        c_tag = lib.taf_read_header(self.c_li_handle) if self.taf_not_maf else lib.maf_read_header(self.c_li_handle)
        p_tags = _c_tags_to_dictionary(c_tag)
        lib.tag_destruct(c_tag)  # Clean up tag
        return p_tags

    def __next__(self):
        """ Get the next alignment block """
        if self.taf_index and not self.sequence_intervals:  # Case we have a taf index but no sequence intervals, we're
            # immediately done
            raise StopIteration  # We're done

        # Read a taf/maf block
        if self.taf_not_maf:  # Is a taf block
            # Use the taf index if present
            c_alignment = lib.tai_next(self._c_taf_index_it, self.c_li_handle) if self.taf_index else \
                lib.taf_read_block(self.p_c_alignment, self.use_run_length_encoding, self.c_li_handle)
        else:  # Is a maf block
            c_alignment = lib.tai_next(self._c_taf_index_it, self.c_li_handle) if self.taf_index else \
                lib.maf_read_block(self.c_li_handle)

        if c_alignment == ffi.NULL:  # If the c_alignment is null
            self.sequence_interval_index += 1
            if self.sequence_intervals is None or self.sequence_interval_index >= len(self.sequence_intervals):
                raise StopIteration  # We're done
            # Otherwise, get next interval
            sequence_name, start, length = self.sequence_intervals[self.sequence_interval_index]
            # Make new iterator for next interval
            self.p_c_alignment = ffi.NULL  # Remove reference to prior alignment
            lib.tai_iterator_destruct(self._c_taf_index_it)  # Cleanup old iterator
            self._c_taf_index_it = lib.tai_iterator(self.taf_index._c_taf_index,
                                                    self.c_li_handle,
                                                    self.use_run_length_encoding,
                                                    _to_c_string(sequence_name), start, length)
            return self.__next__()

        # If maf use O(ND) algorithm to link to any prior alignment block
        if (not self.taf_not_maf) and self.p_c_alignment != ffi.NULL:
            lib.alignment_link_adjacent(self.p_c_alignment, c_alignment, 1)

        # Now add in the rows
        c_row, p_py_row, first_py_row = c_alignment.row, None, None
        while c_row != ffi.NULL:
            # Make the Pythob row
            py_row = Row(c_row=c_row)
            if first_py_row is None:
                first_py_row = py_row

            # Connect the row object to the chain of row objects for the block
            if p_py_row is not None:
                p_py_row._n_row = py_row

            p_py_row = py_row  # Update the previous python row object
            c_row = c_row.n_row  # Move to the next row

        # Now convert the new alignment into Python
        py_alignment = Alignment(c_alignment=c_alignment, py_row=first_py_row)

        # Set the new prior alignment / rows
        self.p_c_alignment = c_alignment

        return py_alignment

    def __iter__(self):
        return self  # Making this an iterable

    def close(self):
        """ Close any associated underlying file """
        if self.taf_index:
            lib.tai_iterator_destruct(self._c_taf_index_it)
        lib.LI_destruct(self.c_li_handle)  # Cleanup the allocated line iterator
        if self.file_string_not_handle:  # Close the underlying file handle
            lib.fclose(self.c_file_handle)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.close()


def get_reference_sequence_intervals(alignment_reader):
    """ Generates a sequence of reference sequence intervals by scanning through the alignment_reader.
    Each generated value is of the form (seq_name, start, length). Length is the number of
    reference bases in blocks within the interval (e.g. ignores unaligned ref gaps).
    """
    seq_intervals = []
    p_ref_seq, p_ref_start, p_ref_length = None, 0, 0
    for alignment in alignment_reader:  # For each alignment block
        # Get the reference row so we can keep track of the reference coordinates
        ref_row = alignment.first_row()
        if ref_row is None:  # Weird case we are on an empty block - technically possible in taf
            continue
        seq_name = ref_row.sequence_name()  # Get the sequence name
        assert seq_name is not None  # Seq name can not be unspecified
        if seq_name != p_ref_seq:  # If we have a new reference sequence
            if p_ref_seq is not None:  # If there is a previous reference sequence, add it to the list
                yield p_ref_seq, p_ref_start, p_ref_length
            p_ref_seq = seq_name  # Set the new reference sequence interval
            p_ref_start = ref_row.start()
            p_ref_length = ref_row.length()
        else:
            p_ref_length += ref_row.length()  # If not a new ref sequence, just add the number of ref bases
            # in the block
    if p_ref_seq is not None:
        yield p_ref_seq, p_ref_start, p_ref_length


def get_column_iterator(alignment_reader,
                        include_sequence_names=True,
                        include_non_ref_columns=True,
                        include_column_tags=False,
                        column_as_int_array=False,
                        column_as_int_array_one_hot=False):
    """ Create an alignment column iterator which returns successive
    columns from the alignment from an AlignmentReader object.

    If alignment_reader was created given a sequence interval only the columns from that given reference
    interval will be returned.

    Each column returned is a column string and a tuple representing the "label" of the column.
    The label tuple is a composed of the reference index of the column and:
        If include_sequence_names is True then the second value is an array of sequence names for the columns.
        If include_column_tags is True then the last value in the label tuple will be a dictionary of any tags.

    If include_non_ref_columns is True then columns not including the reference in the interval will also be returned.

    If column_as_int_array is True then columns will be returned as numpy arrays, see Alignment.get_column_as_np_array()
    """
    if column_as_int_array:
        assert not column_as_int_array_one_hot
    for alignment in alignment_reader:  # For each alignment block
        # Get the reference row so we can keep track of the reference coordinates
        ref_row = alignment.first_row()
        if ref_row is None:  # Weird case we are on an empty block - technically possible in taf
            continue
        ref_index = ref_row.start()
        ref_bases = ref_row.bases()

        # Determine the kind of column to return
        if column_as_int_array:
            get_column = alignment.get_column_as_np_array
        elif column_as_int_array_one_hot:
            get_column = alignment.get_column_as_np_array_one_hot
        else:
            get_column = alignment.get_column

        # If the output wants to match the column entries to the sequences
        if include_sequence_names:
            sequence_names = alignment.get_column_sequences()
            if include_non_ref_columns:
                if include_column_tags:
                    for i in range(alignment.column_number()):
                        yield get_column(i), (ref_index, sequence_names, alignment.column_tags(i))
                        if ref_bases[i] != '-':
                            ref_index += 1
                else:
                    for i in range(alignment.column_number()):
                        yield get_column(i), (ref_index, sequence_names)
                        if ref_bases[i] != '-':
                            ref_index += 1
            else:
                if include_column_tags:
                    for i in range(alignment.column_number()):
                        if ref_bases[i] != '-':
                            yield get_column(i), (ref_index, sequence_names, alignment.column_tags(i))
                            ref_index += 1
                else:
                    for i in range(alignment.column_number()):
                        if ref_bases[i] != '-':
                            yield get_column(i), (ref_index, sequence_names)
                            ref_index += 1
        else:
            if include_non_ref_columns:
                if include_column_tags:
                    for i in range(alignment.column_number()):
                        yield get_column(i), (ref_index, alignment.column_tags(i))
                        if ref_bases[i] != '-':
                            ref_index += 1
                else:
                    for i in range(alignment.column_number()):
                        yield get_column(i), (ref_index,)
                        if ref_bases[i] != '-':
                            ref_index += 1
            else:
                if include_column_tags:
                    for i in range(alignment.column_number()):
                        if ref_bases[i] != '-':
                            yield get_column(i), (ref_index, alignment.column_tags(i))
                            ref_index += 1
                else:
                    for i in range(alignment.column_number()):
                        if ref_bases[i] != '-':
                            yield get_column(i), (ref_index,)
                            ref_index += 1


def get_window_iterator(alignment_reader,
                        window_length=10, step=1,
                        include_sequence_names=True,
                        include_non_ref_columns=True,
                        include_column_tags=False,
                        column_as_int_array=False,
                        column_as_int_array_one_hot=False):
    """ Iterate over (overlapping) windows of the alignment. If columns are numpy arrays the return value
    will be a multidimensional matrix with the first index being the columns. Otherwise columns are returned
    as an array.

    :param alignment_reader: An alignment reader to iterate from
    :param window_length: The number of successive columns to include in a window, must be > 0
    :param step: The number of columns between successive windows, must be 0 < step <= window_length.
    If equal to the window length will mean windows are non-overlapping
    :param include_sequence_names: See get_column_iterator()
    :param include_non_ref_columns: See get_column_iterator()
    :param include_column_tags: See get_column_iterator()
    :param column_as_int_array: See get_column_iterator()
    :param column_as_int_array_one_hot: See get_column_iterator()
    :return: A numpy array of columns, each column from get_column_iterator(), and another numpy array of labels
    """
    assert window_length > 0  # Window length must be positive integer
    assert step <= window_length  # Step can not exceed the window length
    assert step >= 1  # Step can not be negative or 0
    q = deque()
    # If we have numpy arrays then we can concatenate them to create a single tensor
    column_it = get_column_iterator(alignment_reader,
                                    include_sequence_names=include_sequence_names,
                                    include_non_ref_columns=include_non_ref_columns,
                                    include_column_tags=include_column_tags,
                                    column_as_int_array=column_as_int_array,
                                    column_as_int_array_one_hot=column_as_int_array_one_hot)
    if column_as_int_array or column_as_int_array_one_hot:
        for column, labels in column_it:
            column.shape = (1,) + column.shape  # As we will be joining the column we add an extra "prefix" dimension
            q.append((column, labels))  # Add to the right end of the window
            assert len(q) <= window_length
            if len(q) == window_length:
                labels = np.empty(window_length, dtype=object)
                for i in range(window_length):  # Fill out the column and label arrays
                    labels[i] = q[i][1]
                columns = np.concatenate(tuple(i[0] for i in q))  # Concatenate together the column arrays
                yield columns, labels
                for i in range(step):
                    q.popleft()  # Remove from the left end of the window
    else:  # Otherwise, the columns are strings, and we return them as a 1-d array of strings for each window
        for column, labels in column_it:
            q.append((column, labels))  # Add to the right end of the window
            assert len(q) <= window_length
            if len(q) == window_length:
                columns, labels = np.empty(window_length, dtype=object), np.empty(window_length, dtype=object)
                for i in range(window_length):  # Fill out the column and label arrays
                    columns[i] = q[i][0]
                    labels[i] = q[i][1]
                yield columns, labels
                for i in range(step):
                    q.popleft()  # Remove from the left end of the window


def write_taf_index_file(taf_file, index_file, index_block_size=10000):
    """ Create a taf index file """
    c_taf_file_handle = _get_c_file_handle(taf_file)
    c_li_handle = lib.LI_construct(c_taf_file_handle)
    c_index_file_handle = _get_c_file_handle(index_file, "w")
    lib.tai_create(c_li_handle, c_index_file_handle, index_block_size)
    lib.LI_destruct(c_li_handle)  # Cleanup the allocated line iterator
    if isinstance(taf_file, str):  # Close the underlying file handle if opened
        lib.fclose(c_taf_file_handle)
    if isinstance(index_file, str):  # Close the underlying file handle
        lib.fclose(c_index_file_handle)


class AlignmentWriter:
    """ Taf or maf alignment writer.
    """

    def __init__(self, file, taf_not_maf=True, header_tags=None, repeat_coordinates_every_n_columns=-1,
                 use_compression=False):
        """ Use taf_not_maf to switch between MAF or TAF writing.

        File can be either a Python file handle or a file string.
        Handing in the file name is much faster as it avoids using a Python file object. If using
        with a file name remember to close the file, either the with keyword or with the close method.

        If using taf, to run length encode the bases include a tag in the header tags:
         "run_length_encode_bases"=1

        If use_compression is True then will use bgzf compression on output."""
        self.taf_not_maf = taf_not_maf
        self.header_tags = {} if header_tags is None else header_tags
        self.p_py_alignment = None  # The previous alignment
        self.c_lw_handle = lib.LW_construct(_get_c_file_handle(file, "w"), use_compression)
        self.file_string_not_handle = isinstance(file, str)
        self.repeat_coordinates_every_n_columns = repeat_coordinates_every_n_columns

    def write_header(self):
        """ Write the header line """
        c_tag = _dictionary_to_c_tags(self.header_tags)
        lib.taf_write_header(c_tag, self.c_lw_handle) if self.taf_not_maf else \
            lib.maf_write_header(c_tag, self.c_lw_handle)
        lib.tag_destruct(c_tag)

    def write_alignment(self, alignment):
        """ Writes the next alignment block """
        if self.taf_not_maf:
            lib.taf_write_block(self.p_py_alignment._c_alignment if self.p_py_alignment else ffi.NULL,
                                alignment._c_alignment,
                                "run_length_encode_bases" in self.header_tags and
                                int(self.header_tags["run_length_encode_bases"]),
                                self.repeat_coordinates_every_n_columns,
                                self.c_lw_handle)
        else:
            lib.maf_write_block(alignment._c_alignment,
                                self.c_lw_handle)
        self.p_py_alignment = alignment

    def close(self):
        """ Close any associated file """
        lib.LW_destruct(self.c_lw_handle, self.file_string_not_handle)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.close()

