# Using Python Library

The following is a brief, interactive tutorial.
For the complete API see the taffy/lib.py module.
Having installed the Taffy Python library, you should be able to run the following import command without
error. It assumes you are running it from the tests directory of the package. If not, you will need to provide your own demo maf file:

```
# Import the AlignmentReader
from taffy.lib import AlignmentReader
```

Next, let's try opening a maf file:

```
import pathlib
test_maf_file = (pathlib.Path().absolute() / "./evolverMammals.maf").as_posix()
with AlignmentReader(test_maf_file) as mp:
    header = mp.get_header()
    print(header) 
...
{'version': '1', 'scoring': 'N/A'}
```

The dictionary you see printed on the last line represents the header
line of the file.

Okay, suppose instead of a MAF we want to load a TAF file. First, make a TAF file from the test maf file (in tests/) (from the command line). This file will also be created by the tests, so it may already exist.

```
taffy view --inputFile ./evolverMammals.maf --useCompression > ./evolverMammals.taf.gz
```

Now, let's load the file (it doesn't matter if you used compression or not):

```
test_taf_file = (pathlib.Path().absolute() / "./evolverMammals.taf.gz").as_posix()
with AlignmentReader(test_taf_file) as mp:
    print(mp.get_header())
```

Note the Python is identical and the output will be the same as with the maf file.
Next let's try iterating through
some alignment blocks:

```
with AlignmentReader(test_taf_file) as mp:
    for i, block in zip(range(5), mp): # Get the first five blocks in the file
        print(block, "\n") # Print the blocks, adding a newline between each block
...       
Anc0.Anc0refChr0        0       50      +       4151    GTCAAGCTCAGTAGATACTGGATTAGGAATTCATGAGTTAAGCTGTAGCC
Anc1.Anc1refChr1        292714  50      +       296994  GTCAAGCTCAGTAGATACTGGATTAGGAATTCATGAGTTAAGCTGTAGCC
Anc2.Anc2refChr1        5       50      +       4655    GTCAAGCTCAGTTGATGCTGGATTAGGAATTCATGAGTTAAGCTGTAGTC
mr.mrrefChr1    178277  50      +       182340  GTCAAGCTCTGTACATACTAGATTGGACATTCATGGATGAAACTGTGACT
simCow_chr6.simCow.chr6 5045    50      -       602619  GTGAAGCTCAGTTGATGCTGGATTGGGAACTCATGAGTTAAGCTGTAAGC
simDog_chr6.simDog.chr6 589129  50      +       593897  GTCAAGCTCAGTTGGTGCTGGATTAAGAATTCATGAGTTAGGCTGCAGTC
simHuman_chr6.simHuman.chr6     597375  50      +       601863  GTCAAGCTCAGTAGATATTGGATTAGGAATTCATAAGTTAACCTGTAGCC
simMouse_chr6.simMouse.chr6     630640  50      +       636262  GTCAAGCATTGTACATACTAGATTGGACATTCATGGATGACAATGTGACT
simRat_chr6.simRat.chr6 642153  50      +       647215  GTCAAGCTCTGTAAATAGTAGATTGGACATTCATGGATGAAACTGTGCCT
etc.
```

The alignment you see printed represents the first 'block' in the file. Next, let's iterate on the rows in a block:

```
with AlignmentReader(test_taf_file) as mp:
    block = next(mp) # Get the first block in the file
    print(f"Row number: {block.row_number()}") # Number of sequences aligned
    print(f"Column number: {block.column_number()}") # Number of alignment columns
    for row in block: # For each row in the block, print it and show some functionality of rows:
        print(f"Sequence name: {row.sequence_name()}")
        print(f"Start: {row.start()}") # The first position in the sequence
        print(f"Length: {row.length()}") # The number of bases in the alignment
        print(f"Sequence length: {row.sequence_length()}") # The number of bases in the actual sequence (ignoring gaps, etc. in the alignment)
        print(f"Strand: {row.strand()}") # The strand of the sequence
        print(f"Aligned bases {row.bases()}")
        print(row) # Or just conveniently print the information about the row as a string
...
Row number: 9
Column number: 50
Sequence name: Anc0.Anc0refChr0
Start: 0
Length: 50
Sequence length: 4151
Strand: True
Aligned bases GTCAAGCTCAGTAGATACTGGATTAGGAATTCATGAGTTAAGCTGTAGCC
Anc0.Anc0refChr0        0       50      +       4151    GTCAAGCTCAGTAGATACTGGATTAGGAATTCATGAGTTAAGCTGTAGCC
Sequence name: Anc1.Anc1refChr1
etc.
```

Suppose you want to iterate on the columns of the first five blocks in the
alignment (this time using the MAF file to show it doesn't matter which you start from):

```
with AlignmentReader(test_maf_file) as mp:
    for i, block in zip(range(5), mp): # Get the first five blocks in the file in order
        # Print the sequence names
        print(f'Sequence names: {" ".join(block.get_column_sequences())}')
        for j in range(block.column_number()): # For each column
            print(block.get_column(j)) # Get the column as a string
... 
Sequence names: Anc0.Anc0refChr0 Anc1.Anc1refChr1 Anc2.Anc2refChr1 mr.mrrefChr1 simCow_chr6.simCow.chr6 simDog_chr6.simDog.chr6 simHuman_chr6.simHuman.chr6 simMouse_chr6.simMouse.chr6 simRat_chr6.simRat.chr6
GGGGGGGGG
TTTTTTTTT
CCCCGCCCC
AAAAAAAAA
AAAAAAAAA
GGGGGGGGG
CCCCCCCCC
etc..
```

Now suppose you want to access a specific subalignment. For this you will need an index file, which you can build with taffy index, e.g.:

```
taffy index -i ./evolverMammals.taf.gz
```

Which creates the file ./tests/evolverMammals.taf.tai. Note, this will work with either a MAF or TAF file and with or without compression.
Given this index file, you can open it as follows:

```
from taffy.lib import TafIndex  # Import the TafIndex
taf_index = TafIndex(test_taf_file + ".tai", is_maf=False)
with AlignmentReader(test_taf_file, taf_index=taf_index, sequence_intervals=(("Anc0.Anc0refChr0", 1000, 50),)) as mp:
    print(mp.get_header())
    for block in mp:
        print(block, "\n")
...
{'version': '1', 'scoring': 'N/A'}
Anc0.Anc0refChr0        1000    6       +       4151    GCGCTT
Anc1.Anc1refChr1        293719  6       +       296994  GCGCTT
Anc2.Anc2refChr1        1058    6       +       4655    GCGCTT
mr.mrrefChr1    179071  6       +       182340  GCACTT
simCow_chr6.simCow.chr6 598470  6       +       602619  GCGCTT
simDog_chr6.simDog.chr6 590144  6       +       593897  GAGCTG
simHuman_chr6.simHuman.chr6     598418  6       +       601863  GCGCTT
simMouse_chr6.simMouse.chr6     631494  6       +       636262  GCACTT
simRat_chr6.simRat.chr6 642932  3       +       647215  GCA--- 
...
etc.
```

Which gets a particular subrange of blocks within the given reference sequence interval. To specify multiple intervals, just provide a sequence of multiple such intervals. The interval format is (seq_name, start, length).

If we want to iterate on the columns of the alignment without worrying about
blocks we can use the column iterator:

```
from taffy.lib import TafIndex, AlignmentReader, get_column_iterator, get_window_iterator
import pathlib
test_taf_file = (pathlib.Path().absolute() / "./evolverMammals.taf.gz").as_posix()
taf_index = TafIndex(test_taf_file + ".tai", is_maf=False)
with AlignmentReader(test_taf_file, taf_index=taf_index, sequence_intervals=(("Anc0.Anc0refChr0", 1000, 10),))  as mp:
    for column_string, label in get_column_iterator(mp, include_sequence_names=False):
        print(label, column_string)

(1000,) GGGGGGGGG
(1001,) CCCCCACCC
(1002,) GGGAGGGAA
(1003,) CCCCCCCC-
(1004,) TTTTTTTT-
(1005,) TTTTTGTT-
(1006,) AAAAAAAA
(1007,) CCCCTCTCC
(1008,) TTTTTTTTT
(1009,) AATACTAAC
```

Note, the label is a tuple containing the reference index and (optionally, here not shown) the names of the sequences.

Or if we wish to get windows of successive columns of the alignment. Here we get windows of five columns, stepping three columns between
successive windows:

```asm
from taffy.lib import TafIndex, AlignmentReader, get_column_iterator, get_window_iterator
import pathlib
test_taf_file = (pathlib.Path().absolute() / "./evolverMammals.taf.gz").as_posix()
taf_index = TafIndex(test_taf_file + ".tai", is_maf=False)
with AlignmentReader(test_taf_file, taf_index=taf_index, sequence_intervals=(("Anc0.Anc0refChr0", 1000, 10),)) as mp:
    for columns in get_window_iterator(mp, 
    window_length=5, step=3,
    include_sequence_names=False):
        print(columns)
...
(array(['GGGGGGGGG', 'CCCCCACCC', 'GGGAGGGAA', 'CCCCCCCC-', 'TTTTTTTT-'],
      dtype=object), array([(1000,), (1001,), (1002,), (1003,), (1004,)], dtype=object))
(array(['CCCCCCCC-', 'TTTTTTTT-', 'TTTTTGTT-', 'AAAAAAAA', 'CCCCTCTCC'],
      dtype=object), array([(1003,), (1004,), (1005,), (1006,), (1007,)], dtype=object))
```

The returned value of this iterator is a numpy array containing the sequence of the columns, and a numpy array of the labels for each column. Use the column_as_int_array option to get the columns using
a integer encoding or the column_as_int_array_one_hot to get each
column encoded using a one hot encoding.