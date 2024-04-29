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
        self.use_run_length_encoding = "run_length_encode_bases" in self.header_tags  # Use run length encoding
        self.sequence_intervals = sequence_intervals
        self.sequence_interval_index = 0

        if taf_index:
            assert sequence_intervals  # Must not be None
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
        # Read a taf/maf block
        if self.taf_not_maf:  # Is a taf block
            # Use the taf index if present
            c_alignment = lib.tai_next(self._c_taf_index_it, self.c_li_handle) if self.taf_index else \
                lib.taf_read_block(self.p_c_alignment, 0, self.c_li_handle)
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


def get_column_iterator(alignment_reader,
                        include_sequence_names=True,
                        include_non_ref_columns=True,
                        include_column_tags=False):
    """ Create an alignment column iterator which returns successive
    columns from the alignment from an AlignmentReader object.

    If alignment_reader was created given a sequence interval only the columns from that given reference
    interval will be returned.

    Each column returned is a reference index and string representing the column.
    If include_sequence_names is True then a third value, an array of sequence names for the columns, is included.

    If include_non_ref_columns is True then columns not including the reference in the interval will also be returned.

    If include_column_tags is True then the last returned value per column will be a dictionary of any tags.
    """
    for alignment in alignment_reader:  # For each alignment block
        # Get the reference row so we can keep track of the reference coordinates
        ref_row = alignment.first_row()
        if ref_row is None:  # Weird case we are on an empty block - technically possible in taf
            continue
        ref_index = ref_row.start()
        ref_bases = ref_row.bases()

        # If the output wants to match the column entries to the sequences
        if include_sequence_names:
            sequence_names = alignment.get_column_sequences()

        for i in range(alignment.column_number()):  # For each column
            if ref_bases[i] != '-':  # If a ref column
                if include_sequence_names:
                    if include_column_tags:
                        yield ref_index, alignment.get_column(i), sequence_names, alignment.column_tags(i)
                    else:
                        yield ref_index, alignment.get_column(i), sequence_names
                else:
                    if include_column_tags:
                        yield ref_index, alignment.get_column(i), alignment.column_tags(i)
                    else:
                        yield ref_index, alignment.get_column(i)
                ref_index += 1
            elif include_non_ref_columns:  # If we also want non-reference columns
                if include_sequence_names:
                    if include_column_tags:
                        yield ref_index, alignment.get_column(i), sequence_names, alignment.column_tags(i)
                    else:
                        yield ref_index, alignment.get_column(i), sequence_names
                else:
                    if include_column_tags:
                        yield ref_index, alignment.get_column(i), alignment.column_tags(i)
                    else:
                        yield ref_index, alignment.get_column(i)


def get_window_iterator(alignment_reader,
                        window_length=10, step=1,
                        include_sequence_names=True,
                        include_non_ref_columns=True,
                        include_column_tags=False):
    """ Iterate over (overlapping) windows of the alignment.

    :param alignment_reader: An alignment reader to iterate from
    :param window_length: The number of successive columns to include in a window, must be > 0
    :param step: The number of columns between successive windows, must be 0 < step <= window_length.
    If equal to the window length will mean windows are non-overlapping
    :param include_sequence_names: See get_column_iterator()
    :param include_non_ref_columns: See get_column_iterator()
    :param include_column_tags: See get_column_iterator()
    :return: A numpy array of columns, each column from get_column_iterator()
    """
    assert window_length > 0  # Window length must be positive integer
    assert step <= window_length  # Step can not exceed the window length
    assert step >= 1  # Step can not be negative or 0
    q = deque()
    for column in get_column_iterator(alignment_reader,
                                      include_sequence_names=include_sequence_names,
                                      include_non_ref_columns=include_non_ref_columns,
                                      include_column_tags=include_column_tags):
        q.append(column)  # Add to the right end of the window
        assert len(q) <= window_length
        if len(q) == window_length:
            columns = np.empty(window_length, dtype=object)
            for i in range(len(q)): # Fill out the columns
                columns[i] = q[i]
            yield columns
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

