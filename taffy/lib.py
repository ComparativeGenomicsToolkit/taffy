from taffy._taffy_cffi import ffi, lib


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

    def __del__(self):
        lib.alignment_destruct(self._c_alignment, 0)  # Cleans up the underlying C alignment structure


class Row:
    """ Represents a row of an alignment block. See taf.h """
    def __init__(self, c_row=None, l_row=None, r_row=None, n_row=None):
        self._c_row = c_row  # The underlying C row
        self._l_row = l_row  # The prior (left) row in the previous alignment block
        self._r_row = r_row  # The next (right) row in the next alignemnt clock
        self._n_row = n_row  # The next row in the sequence of rows

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

    def left_row(self):
        """ Get any left row in the prior alignment block """
        return self._l_row

    def right_row(self):
        """ Get any right row in the next alignment block """
        return self._r_row

    def __del__(self):
        lib.alignment_row_destruct(self._c_row)  # Cleans up the underlying C alignment structure


class AlignmentParser:
    """ Taf or maf alignment parser.
    """
    def __init__(self, file, taf_not_maf=True, use_run_length_encoding=False, file_string_not_handle=True):
        """ Use taf_not_maf to switch between MAF or TAF parsing.

        Set use_run_length_encoding to determine how
        to decode taf columns. If unknown can be set after construction by interrogating the header line after
        construction.

        Use the file_string_not_handle to determine if file is a file_handle (if file_string_not_handle=False) or
        a file name. Handing in the file name is much faster as it avoids using a Python file object. If using
        with a file name remember to close the file, either the with keyword or with the close method.
        """
        self.taf_not_maf = taf_not_maf
        self.use_run_length_encoding = use_run_length_encoding
        self.p_c_alignment = ffi.NULL  # The previous C alignment returned
        self.p_c_rows_to_py_rows = {}  # Hash from C rows to Python rows of the previous
        # alignment block, allowing linking of rows between blocks
        self.c_file_handle = lib.fopen(_to_c_string(file), _to_c_string("r")) if file_string_not_handle \
            else ffi.cast("FILE *", file)  # Either open the file or cast the python file handle to a C file handle
        # note the Python file handle is *way* slower
        self.file_string_not_handle = file_string_not_handle
        if taf_not_maf:  # If taf then we wrap the file handle in a C line iterator
            self.c_file_handle = lib.LI_construct(self.c_file_handle)

    def get_header(self):
        """ Get tags from the header line as a dictionary of key:value pairs. Must be called if a header is
         present in the file before blocks are retrieved """
        c_tag = lib.taf_read_header(self.c_file_handle) if self.taf_not_maf else lib.maf_read_header(self.c_file_handle)
        p_tags = _c_tags_to_dictionary(c_tag)
        lib.tag_destruct(c_tag)  # Clean up tag
        return p_tags

    def __next__(self):
        # Read a taf/maf block
        c_alignment = lib.taf_read_block(self.p_c_alignment, 0, self.c_file_handle) if self.taf_not_maf \
                      else lib.maf_read_block(self.c_file_handle)

        if c_alignment == ffi.NULL:  # If the c_alignment is null
            raise StopIteration  # We're done

        # If maf use O(ND) algorithm to link to any prior alignment block
        if (not self.taf_not_maf) and self.p_c_alignment != ffi.NULL:
            lib.alignment_link_adjacent(self.p_c_alignment, c_alignment, 1)

        # Now add in the rows
        c_row, p_py_row, c_rows_to_py_rows = c_alignment.row, None, {}
        while c_row != ffi.NULL:
            # The previous py row

            # Make the Python row object
            if c_row.l_row == ffi.NULL:  # If there is no prior left row to connect to
                py_row = Row(c_row=c_row)
            else:  # Otherwise, there is a prior left row to connect to
                l_py_row = self.p_c_rows_to_py_rows[c_row.l_row]
                py_row = Row(c_row=c_row, l_row=l_py_row)
                l_py_row._r_row = py_row

            # Add to the map of c rows to python rows
            c_rows_to_py_rows[c_row] = py_row

            # Connect the row object to the chain of row objects for the block
            if p_py_row is not None:
                p_py_row._n_row = py_row

            p_py_row = py_row  # Update the previous python row object
            c_row = c_row.n_row  # Move to the next row

        # Now convert the new alignment into Python
        py_alignment = Alignment(c_alignment=c_alignment, py_row=c_rows_to_py_rows[c_alignment.row])

        # Set the new prior alignment / rows
        self.p_c_alignment = c_alignment
        self.p_c_rows_to_py_rows = c_rows_to_py_rows

        return py_alignment

    def __iter__(self):
        return self  # Making this an iterable

    def close(self):
        """ Close any associated underlying file """
        if self.taf_not_maf:  # If taf, cleanup the allocated line iterator
            i = self.c_file_handle
            self.c_file_handle = self.c_file_handle.fh
            lib.LI_destruct(i)
        if self.file_string_not_handle:  # Close the underlying file handle
            lib.fclose(self.c_file_handle)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.close()


class AlignmentWriter:
    """ Taf or maf alignment writer.
    """

    def __init__(self, file, taf_not_maf=True, header_tags=None,
                 repeat_coordinates_every_n_columns=-1,
                 file_string_not_handle=True):
        """ Use taf_not_maf to switch between MAF or TAF writing. Set use_run_length_encoding to determine how
        to encode taf columns.

        Use the file_string_not_handle to determine if file is a file_handle (if file_string_not_handle=False) or
        a file name. Handing in the file name is much faster as it avoids using a Python file object. If using
        with a file name remember to close the file, either the with keyword or with the close method. """
        self.taf_not_maf = taf_not_maf
        self.header_tags = {} if header_tags is None else header_tags
        self.p_py_alignment = None  # The previous alignment
        self.c_file_handle = lib.fopen(_to_c_string(file), _to_c_string("w")) if file_string_not_handle \
            else ffi.cast("FILE *", file)  # Either open the file or cast the python file handle to a C file handle
        # note the Python file handle is *way* slower
        self.file_string_not_handle = file_string_not_handle
        self.repeat_coordinates_every_n_columns = repeat_coordinates_every_n_columns

    def write_header(self):
        """ Write the header line """
        c_tag = _dictionary_to_c_tags(self.header_tags)
        lib.taf_write_header(c_tag, self.c_file_handle) if self.taf_not_maf else \
            lib.maf_write_header(c_tag, self.c_file_handle)
        lib.tag_destruct(c_tag)

    def write_alignment(self, alignment):
        """ Writes the next alignment block """
        if self.taf_not_maf:
            lib.taf_write_block(self.p_py_alignment._c_alignment if self.p_py_alignment else ffi.NULL,
                                alignment._c_alignment,
                                "run_length_encode_bases" in self.header_tags and
                                int(self.header_tags["run_length_encode_bases"]),
                                self.repeat_coordinates_every_n_columns,
                                self.c_file_handle)
        else:
            lib.maf_write_block(alignment._c_alignment,
                                self.c_file_handle)
        self.p_py_alignment = alignment

    def close(self):
        """ Close any associated file """
        if self.file_string_not_handle:  # Close the underlying file handle
            lib.fclose(self.c_file_handle)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.close()
