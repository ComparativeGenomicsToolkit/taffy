import unittest
import pathlib
import subprocess
from random import randint
from taffy.lib import AlignmentReader, AlignmentWriter, TafIndex, write_taf_index_file, get_column_iterator


class TafTest(unittest.TestCase):
    def setUp(self):
        self.test_maf_file = (pathlib.Path().absolute() / "../tests/evolverMammals.maf").as_posix()
        self.test_taf_file = (pathlib.Path().absolute() / "../tests/evolverMammals.taf").as_posix()
        self.test_index_file = (pathlib.Path().absolute() / "../tests/evolverMammals.tai").as_posix()

    def test_maf_reader(self):
        """ Manually test the first couple blocks from the maf file """
        with AlignmentReader(self.test_maf_file, taf_not_maf=False) as mp:
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

            # Check left and right connections
            self.assertEqual(row.left_row(), a.first_row())
            self.assertEqual(a.first_row().right_row(), row)

            # Check column retrieval
            self.assertEqual(b.get_column(0), "AAAATAAAA")
            self.assertEqual(b.get_column(-1), "AAAATAAAA")

    def test_column_iterator(self):
        """ Manually test the column iterator """
        with AlignmentReader(self.test_maf_file, taf_not_maf=False) as mp:
            column_it = get_column_iterator(mp)
            s = ["Anc0.Anc0refChr0", "Anc1.Anc1refChr1", "Anc2.Anc2refChr1", "mr.mrrefChr1",
                 "simCow_chr6.simCow.chr6", "simDog_chr6.simDog.chr6", "simHuman_chr6.simHuman.chr6",
                 "simMouse_chr6.simMouse.chr6", "simRat_chr6.simRat.chr6"]
            # First column
            sequence_names, column = next(column_it)
            self.assertEqual(s, list(sequence_names))
            self.assertEqual(column, "GGGGGGGGG")

            # Second column
            sequence_names, column = next(column_it)
            self.assertEqual(s, list(sequence_names))
            self.assertEqual(column, "TTTTTTTTT")

            # Third column
            sequence_names, column = next(column_it)
            self.assertEqual(s, list(sequence_names))
            self.assertEqual(column, "CCCCGCCCC")

            # Check it works
            for sequence_names, column in column_it:
                pass

    def test_maf_to_taf(self, compress_file=False):
        """ Read a maf file, write a taf file, compress it with gzip and then read it back and check
        they are equal. Tests round trip read and write. Writes in random tags to the taf to test tag writing """
        def make_random_tags():  # Fn to make random tags
            return {str(randint(0, 1000)): str(randint(0, 1000)) for i in range(randint(0, 5))}
        column_tags = []  # List of tags per column

        # First read from the maf file and write the taf file
        with AlignmentReader(self.test_maf_file, taf_not_maf=False) as mp:
            maf_header_tags = mp.get_header()  # Get the maf header tags

            with AlignmentWriter(self.test_taf_file, header_tags=maf_header_tags) as tw:
                tw.write_header()  # Write the header

                for a in mp:  # For each alignment block in input
                    # Add random tags
                    for i in range(a.column_number()):
                        column_tags.append(make_random_tags())
                        a.set_column_tags(i, column_tags[-1])

                    tw.write_alignment(a)  # Write a corresponding output block

        # Compress the file
        if compress_file:
            subprocess.run(["gzip", "-f", self.test_taf_file])
            self.test_taf_file = f"{self.test_taf_file}.gz"

        # Now read back the taf file and check it is equivalent to the maf
        column_index = 0  # Used to track where in the list of columns we are
        with AlignmentReader(self.test_taf_file) as tp:
            with AlignmentReader(self.test_maf_file, taf_not_maf=False) as mp:
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

    def test_taf_index(self, compress_file=False):
        """ Index a taf file then load a portion of the file with the
         AlignmentReader """
        # Convert MAF to TAF
        with AlignmentReader(self.test_maf_file, taf_not_maf=False) as mp:
            maf_header_tags = mp.get_header()  # Get the maf header tags
            with AlignmentWriter(self.test_taf_file, header_tags=maf_header_tags) as tw:
                tw.write_header()  # Write the header
                for a in mp:  # For each alignment block in input
                    tw.write_alignment(a)  # Write a corresponding output block

        # Compress the file
        if compress_file:
            subprocess.run(["gzip", "-f", self.test_taf_file])
            self.test_taf_file = f"{self.test_taf_file}.gz"

        # Write the index file
        write_taf_index_file(taf_file=self.test_taf_file, index_file=self.test_index_file)

        # Make the Taf Index object
        taf_index = TafIndex(self.test_index_file)

        # Create a taf reader
        with AlignmentReader(self.test_taf_file,
                             taf_index=taf_index,
                             sequence_name="Anc0.Anc0refChr0",
                             start=100,
                             length=500) as tp:

            for a in tp:  # For each alignment block in input
                ref_row = a.first_row()  # Check the reference row is in the bounds of the search
                self.assertEqual(ref_row.sequence_name(), "Anc0.Anc0refChr0")
                self.assertTrue(ref_row.start() >= 100)
                self.assertTrue(ref_row.start() + ref_row.length() <= 600)

    def test_taf_index_compression(self):
        self.test_taf_index(compress_file=True)


if __name__ == '__main__':
    unittest.main()
