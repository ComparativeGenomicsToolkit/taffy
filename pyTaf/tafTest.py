import unittest
import pathlib
from random import randint
from taf import AlignmentParser, AlignmentWriter


class TafTest(unittest.TestCase):
    def setUp(self):
        self.test_maf_file = (pathlib.Path().absolute() / "../tests/evolverMammals.maf").as_posix()
        self.test_taf_file = (pathlib.Path().absolute() / "../tests/evolverMammals.taf").as_posix()

    def test_maf_parser(self):
        """ Manually test the first couple blocks from the maf file """
        with AlignmentParser(self.test_maf_file, taf_not_maf=False) as mp:
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

    def test_maf_to_taf(self):
        """ Read a maf file, write a taf file and then read it back and check
        they are equal. Tests round trip read and write. Writes in random tags to the taf to test tag writing """
        def make_random_tags():  # Fn to make random tags
            return {str(randint(0, 1000)): str(randint(0, 1000)) for i in range(randint(0, 5))}
        column_tags = []  # List of tags per column

        # First read from the maf file and write the taf file
        with AlignmentParser(self.test_maf_file, taf_not_maf=False) as mp:
            maf_header_tags = mp.get_header()  # Get the maf header tags

            with AlignmentWriter(self.test_taf_file, header_tags=maf_header_tags) as tw:
                tw.write_header()  # Write the header

                for a in mp:  # For each alignment block in input
                    # Add random tags
                    for i in range(a.column_number()):
                        column_tags.append(make_random_tags())
                        a.set_column_tags(i, column_tags[-1])

                    tw.write_alignment(a)  # Write a corresponding output block

        # Now read back the taf file and check it is equivalent to the maf
        column_index = 0  # Used to track where in the list of columns we are
        with AlignmentParser(self.test_taf_file) as tp:

            with AlignmentParser(self.test_maf_file, taf_not_maf=False) as mp:
                self.assertEqual(mp.get_header(), tp.get_header())  # Check headers are equivalent

                for ma, ta in zip(mp, tp):  # For each of the two alignment blocks
                    # Check alignment blocks have same stats
                    self.assertEqual(ma.row_number(), ta.row_number())  # Should be nine rows in block
                    self.assertEqual(ma.column_number(), ta.column_number())
                    for i in range(ma.column_number()):
                        self.assertEqual(column_tags[column_index], ta.column_tags(i))
                        column_index += 1

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


if __name__ == '__main__':
    unittest.main()
