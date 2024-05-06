import unittest
import pathlib
import subprocess
from random import randint

import taffy.lib
from taffy.lib import AlignmentReader, AlignmentWriter, TafIndex, write_taf_index_file, \
    get_column_iterator, get_window_iterator
from taffy.newick import PhyloTree
from taffy.ml import TorchDatasetAlignmentIterator
from torch.utils.data import DataLoader


class TafTest(unittest.TestCase):
    def setUp(self):
        self.test_maf_file = (pathlib.Path().absolute() / "../tests/evolverMammals.maf").as_posix()
        self.test_taf_file = (pathlib.Path().absolute() / "../tests/evolverMammals.taf").as_posix()
        self.test_index_file = (pathlib.Path().absolute() / "../tests/evolverMammals.tai").as_posix()

    def test_maf_reader(self):
        """ Manually test the first couple blocks from the maf file """
        with AlignmentReader(self.test_maf_file) as mp:
            # Check the header is as expected
            header = mp.get_header()
            self.assertEqual(header, {"version": "1", "scoring": "N/A"})

            # Now check a couple of alignment blocks
            a = next(mp)

            # The first alignment block
            self.assertEqual(a.row_number(), 9)  # Should be nine rows in block
            self.assertEqual(a.column_number(), 50)
            for i in range(a.column_number()):
                self.assertEqual(a.column_tags(i), {})

            # Check the first row
            row = a.first_row()
            self.assertEqual(row.sequence_name(), "Anc0.Anc0refChr0")
            self.assertEqual(row.start(), 0)
            self.assertEqual(row.length(), 50)
            self.assertEqual(row.sequence_length(), 4151)
            self.assertEqual(row.strand(), True)
            self.assertEqual(row.bases(), "GTCAAGCTCAGTAGATACTGGATTAGGAATTCATGAGTTAAGCTGTAGCC")
            self.assertEqual(row.left_gap_sequence(), "")

            # Check the second row
            row = row.next_row()
            self.assertEqual(row.sequence_name(), "Anc1.Anc1refChr1")
            self.assertEqual(row.start(), 292714)
            self.assertEqual(row.length(), 50)
            self.assertEqual(row.sequence_length(), 296994)
            self.assertEqual(row.strand(), True)
            self.assertEqual(row.bases(), "GTCAAGCTCAGTAGATACTGGATTAGGAATTCATGAGTTAAGCTGTAGCC")
            self.assertEqual(row.left_gap_sequence(), "")

            # Check column retrieval
            self.assertEqual(a.get_column(0), "GGGGGGGGG")
            self.assertEqual(a.get_column(1), "TTTTTTTTT")
            self.assertEqual(a.get_column(2), "CCCCGCCCC")
            self.assertEqual(a.get_column(a.column_number()-1), "CCCTCCCTT")
            self.assertEqual(a.get_column(-1), "CCCTCCCTT")
            self.assertEqual(a.get_column(a.column_number()-2), "CCTCGTCCC")
            self.assertEqual(a.get_column(-2), "CCTCGTCCC")

            # The second alignment block
            b = next(mp)
            self.assertEqual(b.row_number(), 9)  # Should be nine rows in block
            self.assertEqual(b.column_number(), 1)  # Is a single base long
            for i in range(b.column_number()):
                self.assertEqual(b.column_tags(i), {})

            # Check the connections between adjacent blocks/rows
            row = b.first_row()
            self.assertEqual(row.sequence_name(), "Anc0.Anc0refChr0")
            self.assertEqual(row.start(), 50)
            self.assertEqual(row.length(), 1)
            self.assertEqual(row.sequence_length(), 4151)
            self.assertEqual(row.strand(), True)
            self.assertEqual(row.bases(), "A")
            self.assertEqual(row.left_gap_sequence(), "")

            # Check column retrieval
            self.assertEqual(b.get_column(0), "AAAATAAAA")
            self.assertEqual(b.get_column(-1), "AAAATAAAA")

    def test_column_iterator(self):
        """ Manually test the column iterator """
        with AlignmentReader(self.test_maf_file) as mp:
            column_it = get_column_iterator(mp)
            s = ["Anc0.Anc0refChr0", "Anc1.Anc1refChr1", "Anc2.Anc2refChr1", "mr.mrrefChr1",
                 "simCow_chr6.simCow.chr6", "simDog_chr6.simDog.chr6", "simHuman_chr6.simHuman.chr6",
                 "simMouse_chr6.simMouse.chr6", "simRat_chr6.simRat.chr6"]
            # First column
            column, (ref_index, sequence_names) = next(column_it)
            self.assertEqual(s, list(sequence_names))
            self.assertEqual(column, "GGGGGGGGG")

            # Second column
            column, (ref_index, sequence_names) = next(column_it)
            self.assertEqual(s, list(sequence_names))
            self.assertEqual(column, "TTTTTTTTT")

            # Third column
            column, (ref_index, sequence_names) = next(column_it)
            self.assertEqual(s, list(sequence_names))
            self.assertEqual(column, "CCCCGCCCC")

            # Check it works
            for column, label in column_it:
                pass

    def test_maf_to_taf(self, compress_file=False):
        """ Read a maf file, write a taf file, compress it with gzip and then read it back and check
        they are equal. Tests round trip read and write. Writes in random tags to the taf to test tag writing """
        def make_random_tags():  # Fn to make random tags
            return {str(randint(0, 1000)): str(randint(0, 1000)) for i in range(randint(0, 5))}
        column_tags = []  # List of tags per column

        # First read from the maf file and write the taf file
        with AlignmentReader(self.test_maf_file) as mp:
            maf_header_tags = mp.get_header()  # Get the maf header tags

            with AlignmentWriter(self.test_taf_file, header_tags=maf_header_tags, use_compression=compress_file) as tw:
                tw.write_header()  # Write the header

                for a in mp:  # For each alignment block in input
                    # Add random tags
                    for i in range(a.column_number()):
                        column_tags.append(make_random_tags())
                        a.set_column_tags(i, column_tags[-1])

                    tw.write_alignment(a)  # Write a corresponding output block

        # Now read back the taf file and check it is equivalent to the maf
        column_index = 0  # Used to track where in the list of columns we are
        with AlignmentReader(self.test_taf_file) as tp:
            with AlignmentReader(self.test_maf_file) as mp:
                self.assertEqual(mp.get_header(), tp.get_header())  # Check headers are equivalent

                for ma, ta in zip(mp, tp):  # For each of the two alignment blocks
                    # Check alignment blocks have same stats
                    self.assertEqual(ma.row_number(), ta.row_number())  # Should be nine rows in block
                    self.assertEqual(ma.column_number(), ta.column_number())
                    for i in range(ma.column_number()):
                        self.assertEqual(column_tags[column_index], ta.column_tags(i))
                        column_index += 1

                    # Check the columns are equal
                    for i in range(ma.column_number()):
                        self.assertEqual(ma.get_column(i), ta.get_column(i))

                    # Check the rows are equal
                    mr, tr = ma.first_row(), ta.first_row()
                    while mr:
                        self.assertEqual(mr.sequence_name(), tr.sequence_name())
                        self.assertEqual(mr.start(), tr.start())
                        self.assertEqual(mr.length(), tr.length())
                        self.assertEqual(mr.sequence_length(), tr.sequence_length())
                        self.assertEqual(mr.strand(), tr.strand())
                        self.assertEqual(mr.bases(), tr.bases())
                        self.assertEqual(mr.left_gap_sequence(), tr.left_gap_sequence())
                        mr = mr.next_row()
                        tr = tr.next_row()
                    self.assertTrue(not tr)

    def test_maf_to_taf_compressed(self):
        self.test_maf_to_taf(compress_file=True)

    def make_taf_and_taf_index(self, compress_file=False, taf_not_maf=True):
        # Convert MAF to TAF (or just write MAF if not taf_not_maf)
        with AlignmentReader(self.test_maf_file) as mp:
            maf_header_tags = mp.get_header()  # Get the maf header tags
            with AlignmentWriter(self.test_taf_file, header_tags=maf_header_tags,
                                 use_compression=compress_file, taf_not_maf=taf_not_maf) as tw:
                tw.write_header()  # Write the header
                for a in mp:  # For each alignment block in input
                    tw.write_alignment(a)  # Write a corresponding output block

        # Write the index file
        write_taf_index_file(taf_file=self.test_taf_file, index_file=self.test_index_file)

        # Make the Taf Index object
        taf_index = TafIndex(self.test_index_file, not taf_not_maf)

        return taf_index

    @staticmethod
    def get_random_sequence_intervals():
        # Gets a bunch of random sequence intervals for testing the iterators with
        sequence_intervals = []
        total_length = 0
        for i in range(randint(1, 10)):
            start = randint(0, 1000)
            length = randint(1, 50)  # Don't include 0 length intervals here, because makes accounting tricky
            total_length += length
            sequence_intervals.append(("Anc0.Anc0refChr0", start, length))
        return sequence_intervals, total_length

    def test_taf_index(self, compress_file=False, taf_not_maf=True):
        """ Index a taf file then load a portion of the file with the
         AlignmentReader """
        # Make the Taf Index object and the taf file
        taf_index = self.make_taf_and_taf_index(compress_file=False, taf_not_maf=True)

        # Create a taf/maf reader
        with AlignmentReader(self.test_taf_file, taf_index=taf_index, sequence_intervals=(("Anc0.Anc0refChr0",
                                                                                           100, 500),)) as tp:

            for a in tp:  # For each alignment block in input
                ref_row = a.first_row()  # Check the reference row is in the bounds of the search
                self.assertEqual(ref_row.sequence_name(), "Anc0.Anc0refChr0")
                self.assertTrue(ref_row.start() >= 100)
                self.assertTrue(ref_row.start() + ref_row.length() <= 600)

        # Test random intervals using the column iterator
        for test in range(100):
            # Make a bunch of random sequence intervals
            sequence_intervals, total_length = self.get_random_sequence_intervals()

            with AlignmentReader(self.test_taf_file, taf_index=taf_index, sequence_intervals=sequence_intervals) as tp:
                i, j, k = 0, 0, 0  # Index of total bases, sequence interval, and offset on sequence interval
                for column, (ref_index, seq_names, column_tags) in get_column_iterator(tp, include_non_ref_columns=False,
                                                                                       include_column_tags=True):
                    i += 1  # Increment total bases
                    seq_name, start, length = sequence_intervals[j]
                    self.assertEqual(seq_names[0], seq_name)
                    self.assertEqual(ref_index, k + start)
                    self.assertEqual(column_tags, {})
                    k += 1  # Increment index along sequence
                    if k >= length:  # If we have walked of the end of a sequence interval, move to the next
                        j, k = j+1, 0
                self.assertEqual(j, len(sequence_intervals))
                self.assertEqual(i, total_length)

        # Test random intervals using the window iterator
        for test in range(100):
            start = randint(0, 1000)
            length = randint(0, 50)
            window_length = randint(1, 10)
            step = randint(1, window_length)

            with AlignmentReader(self.test_taf_file, taf_index=taf_index, sequence_intervals=(("Anc0.Anc0refChr0",
                                                                                               start, length),)) as tp:
                j = start
                for columns, labels in get_window_iterator(tp, include_non_ref_columns=False,
                                                           window_length=window_length, step=step):
                    # Check that we have the expected number of columns
                    self.assertEqual(len(columns), window_length)
                    self.assertEqual(len(labels), window_length)

                    # Check coordinates of columns
                    k = 0
                    for ref_index, seq_names in labels:
                        self.assertEqual(seq_names[0], "Anc0.Anc0refChr0")
                        self.assertEqual(ref_index, j+k)
                        self.assertTrue(ref_index < start + length)
                        k += 1

                    j += step

    def test_taf_index_compression(self):
        self.test_taf_index(compress_file=True)

    def test_maf_index(self):
        self.test_taf_index(compress_file=False, taf_not_maf=False)

    def test_maf_index_compression(self):
        self.test_taf_index(compress_file=True, taf_not_maf=False)

    def test_newick_parser(self):
        """ Manually test newick tree parser """
        a = "((A:0.1,B:0.2)C:0.3,(D:0.4)E:0.5)F:0.6;"
        t = PhyloTree.newick_tree_parser(a)
        self.assertEqual(str(t), a)

        a_no_bl = "((A,B)C,(D)E)F;"
        self.assertEqual(t.tree_string(include_branch_lengths=False), a_no_bl)
        self.assertEqual(PhyloTree.newick_tree_parser(a_no_bl).tree_string(include_branch_lengths=False), a_no_bl)

        a_no_in = "((A:0.1,B:0.2):0.3,(D:0.4):0.5):0.6;"
        self.assertEqual(t.tree_string(include_internal_labels=False), a_no_in)
        self.assertEqual(PhyloTree.newick_tree_parser(a_no_in).tree_string(include_internal_labels=False), a_no_in)

        a_no_nl = "((:0.1,:0.2):0.3,(:0.4):0.5):0.6;"
        self.assertEqual(t.tree_string(include_internal_labels=False, include_leaf_labels=False), a_no_nl)
        self.assertEqual(PhyloTree.newick_tree_parser(a_no_nl).tree_string(include_internal_labels=False,
                                                                           include_leaf_labels=False), a_no_nl)

        a_no_to = "((,),());"
        self.assertEqual(t.tree_string(include_internal_labels=False, include_leaf_labels=False,
                                       include_branch_lengths=False), a_no_to)

    @staticmethod
    def identity_fn(i):
        return i

    def test_torchDatasetAlignmentIterator(self, compress_file=False, taf_not_maf=True):
        """ Tests the PyTorch dataset alignment iterator """
        # Make the Taf Index object and the taf file
        taf_index = self.make_taf_and_taf_index(compress_file=compress_file, taf_not_maf=taf_not_maf)

        # Make random examples
        for test in range(3):
            # Make a bunch of random sequence intervals
            sequence_intervals, total_length = self.get_random_sequence_intervals()

            ai = DataLoader(TorchDatasetAlignmentIterator(self.test_taf_file,
                                                          label_conversion_function=self.identity_fn,
                                                          taf_index_file=self.test_index_file,
                                                          is_maf=not taf_not_maf,
                                                          sequence_intervals=sequence_intervals,
                                                          include_non_ref_columns=False,
                                                          include_sequence_names=True,
                                                          include_column_tags=True),
                            num_workers=1)

            i, j, k = 0, 0, 0  # Index of total bases, sequence interval, and offset on sequence interval
            for column, (ref_index, seq_names, column_tags) in ai:
                i += 1  # Increment total bases
                seq_name, start, length = sequence_intervals[j]
                self.assertEqual(seq_names[0][0], seq_name)
                self.assertEqual(ref_index, k + start)
                self.assertEqual(column_tags, {})  # Columns tags are empty
                k += 1  # Increment index along sequence
                if k >= length:  # If we have walked of the end of a sequence interval, move to the next
                    j, k = j+1, 0
            self.assertEqual(j, len(sequence_intervals))
            self.assertEqual(i, total_length)

            # Now test a multi-processing version of the data loader - here we don't expect the columns to come out
            # in order, necessarily
            ai = DataLoader(TorchDatasetAlignmentIterator(self.test_taf_file,
                                                          label_conversion_function=self.identity_fn,
                                                          taf_index_file=self.test_index_file,
                                                          is_maf=not taf_not_maf,
                                                          sequence_intervals=sequence_intervals,
                                                          include_non_ref_columns=False,
                                                          include_sequence_names=True,
                                                          include_column_tags=True),
                            num_workers=5)
            i = 0
            for column, labels in ai:
                i += 1
            self.assertEqual(i, total_length)

    def test_get_reference_sequence_intervals(self):
        """ Tests the get_reference_sequence_intervals method.
        """
        # Make the Taf Index object and the taf file
        taf_index = self.make_taf_and_taf_index(compress_file=False, taf_not_maf=True)

        # Create a taf/maf reader for a series of random sequence intervals
        for test in range(100):
            start = randint(0, 1000)
            length = randint(1, 50)
            with AlignmentReader(self.test_taf_file, taf_index=taf_index, sequence_intervals=(("Anc0.Anc0refChr0",
                                                                                                start, length),)) as tp:
                seq_intervals = list(taffy.lib.get_reference_sequence_intervals(tp))
                self.assertEqual(seq_intervals, [("Anc0.Anc0refChr0", start, length)])


if __name__ == '__main__':
    unittest.main()
