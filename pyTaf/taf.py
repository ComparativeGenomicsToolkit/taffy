from _pyTaf_cffi import ffi, lib


class Alignment:
    def __init__(self, row_number=0, column_number=0, row=None, column_tags=None):
        self.row_number = row_number
        self.column_number = column_number
        self.row = row
        self.column_tags = column_tags  # Tags are optional, but specified as list of dictionaries, one for
        # each column


class Row:
    def __init__(self, sequence_name, start, length, sequence_length,
                 strand, bases, left_gap_sequence="", l_row=None, r_row=None, n_row=None):
        self.sequence_name = sequence_name
        self.start = start
        self.length = length
        self.sequence_length = sequence_length
        self.strand = strand
        self.bases = bases
        self.left_gap_sequence = left_gap_sequence
        self.l_row = l_row
        self.r_row = r_row
        self.n_row = n_row


def _to_string(s):
    """ Convert a cffi string into a Python string in Python3 """
    return ffi.string(s).decode("utf-8")


class AlignmentParser:
    """ Taf or maf alignment parser.

    Use with a context manager to ensure memory from C bindings is cleaned up """

    def __init__(self, file_handle, taf_not_maf=True, use_run_length_encoding=False):
        self.taf_not_maf = taf_not_maf
        self.use_run_length_encoding = use_run_length_encoding
        self.p_alignment = ffi.NULL  # The previous python alignment returned
        self.p_c_rows_to_python_rows = {}  # Hash from C rows to Python rows of the previous
        # alignment block, allowing linking of rows between blocks
        self.c_file_handle = ffi.cast("FILE *", file_handle)
        if taf_not_maf:
            self.c_file_handle = lib.LI_construct(self.c_file_handle)

    def get_header(self):
        tag = lib.taf_read_header(self.c_file_handle) if self.taf_not_maf else lib.maf_read_header(self.c_file_handle)
        # Convert tags to a Python dictionary
        tags = {}
        t = tag
        while t != ffi.NULL:
            tags[_to_string(t.key)] = _to_string(t.value)
            t = t.n_tag
        lib.tag_destruct(tag) # Clean up tag
        return tags

    def __next__(self):
        # Read a taf/maf block
        ca = lib.taf_read_block(self.p_alignment, 0, self.c_file_handle) if self.taf_not_maf \
            else lib.maf_read_block(self.c_file_handle)

        if ca == ffi.NULL:  # If it is null
            if self.p_alignment != ffi.NULL:  # Clean up the prior alignment, if it exists
                lib.alignment_destruct(self.p_alignment)
                self.p_alignment = ffi.NULL
            raise StopIteration

        # If maf use O(ND) algorithm to link to any prior alignment block
        if (not self.taf_not_maf) and self.p_alignment != ffi.NULL:
            lib.alignment_link_adjacent(self.p_alignment, ca, 1)

        # Convert the new alignment into Python
        pa = Alignment(row_number=ca.row_number,
                       column_number=ca.column_number,
                       row=None, column_tags=[{} for i in range(ca.column_number)])

        # Now add in the rows
        cr, ppr, c_rows_to_python_rows = ca.row, None, {}
        while cr != ffi.NULL:
            # Make the row object, linking it to any prior row
            pr = Row(sequence_name=_to_string(cr.sequence_name),
                     start=cr.start,
                     length=cr.length,
                     sequence_length=cr.sequence_length,
                     strand=cr.strand,
                     bases=_to_string(cr.bases),
                     left_gap_sequence="" if cr.left_gap_sequence == ffi.NULL else _to_string(cr.left_gap_sequence),
                     l_row=None if cr.l_row == ffi.NULL else self.p_c_rows_to_python_rows[cr.l_row])
            c_rows_to_python_rows[cr] = pr

            # Connect the row object to the chain of row objects
            if ppr is None:
                pa.row = pr
            else:
                ppr.n_row = pr

            ppr = pr  # Update the previous python row object
            cr = cr.n_row  # Move to the next row

        # Clean up any prior alignment and set the new prior alignment
        if self.p_alignment != ffi.NULL:
            lib.alignment_destruct(self.p_alignment)
        self.p_alignment = ca
        self.p_c_rows_to_python_rows = c_rows_to_python_rows

        return pa

    def __iter__(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        if self.taf_not_maf:
            lib.LI_destruct(self.c_file_handle)

