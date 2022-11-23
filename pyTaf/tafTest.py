import unittest
import pathlib
from taf import AlignmentParser

class TafTest(unittest.TestCase):
    def setUp(self):
        self.test_maf_file = (pathlib.Path().absolute() / "../tests/evolverMammals.maf").as_posix()

    def test_maf_parser(self):
        """ Manually test the first couple blocks from the maf file """
        with open(self.test_maf_file) as f:
            mp = AlignmentParser(f, taf_not_maf=False)

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
            self.assertEqual(b.column_number(), 1) # Is a single base long
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
        they are equal. Tests round trip read and write. """
        with open(self.test_maf_file) as f:
            pass

if __name__ == '__main__':
    unittest.main()