import unittest
import pathlib
from taf import AlignmentParser

class TafTest(unittest.TestCase):
    def setUp(self):
        self.test_maf_file = (pathlib.Path().absolute() / "../tests/evolverMammals.maf").as_posix()

    def test_maf_parser(self):
        with open(self.test_maf_file) as f:
            with AlignmentParser(f, taf_not_maf=False) as mp:

                # Check the header is as expected
                header = mp.get_header()
                self.assertEqual(header, {"version": "1", "scoring": "N/A"})

                # Now check a couple of alignment blocks
                a = next(mp)

                self.assertEqual(a.row_number, 9)  # Should be nine rows in block
                self.assertEqual(a.column_number, 50)
                self.assertEqual(a.column_tags, [{} for i in range(a.column_number)])

                b = next(mp)
                self.assertEqual(b.row_number, 9)  # Should be nine rows in block
                self.assertEqual(b.column_number, 1) # Is a single base long
                self.assertEqual(b.column_tags, [{} for i in range(b.column_number)])

                #for c in mp:
                #    pass
                #print(c.column_number)


if __name__ == '__main__':
    unittest.main()