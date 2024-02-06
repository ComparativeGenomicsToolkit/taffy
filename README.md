# Taffy

This is a C, Python and CLI library for manipulating/reading/writing TAF (described below) and 
[MAF](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) format multiple
sequence alignments. It allows conversion between the formats. The Python library is built
on top of the C library and is therefore quite fast.

# The TAF File Format

This is the specification for "transposed alignment format", or maybe "terrific alignment format" (.taf). The idea
is to describe a multiple sequence alignment as a series of columns, one per line, with bases in the
columns optionally run-length encoded, and row coordinates given as needed after the bases in each column.
Where coordinates are not given it is assumed the coordinates continue from the previous column.
User defined tags, as key:value pairs, are also included to allow the alignment
columns to be annotated.

The format supports line based indexing for rapid retrieval of any column or contiguous sequence of columns from
the file. The format is intentionally simple but should prove quite space efficient for large alignments 
(see stats below).

Its key potential benefits over the MAF format are that:
* it does not suffer the same issue with fragmentation as the number of sequences grows (there are no blocks!),
* is often less verbose (particularly for large alignments), 
* it supports extensible column annotations.
* apis support working with gzip compressed files (taf.gz, see below)
* we have integrated indexing giving efficient random access to both compressed and gzipped files (see below)

The format is composed of a header and then a sequence of columns.
Tokens are separated by white-space. The syntax is defined as follows:

    .taf -> header '\n' columns
    
    header -> '#taf' tags
    
    tags -> tag tags
         -> ''   
    
    tag -> key':'value
    
    key -> alphanumerical string
    value -> alphanumerical string

(The header is therefore just a possibly empty sequence of key:value pairs)

    columns -> column '\n' columns
               column 
    
    column -> bases coordinates_column tag_string 
           -> bases coordinates_column
           -> bases tag_string
           -> bases
           -> #string

(Columns are encoded one per line. Columns are either comment lines, starting with a "#" symbol and then containing an arbitrary comment string, which is ignored, or an encoding of an alignment of the bases in the column (see below) and then some optional coordinates and tags.)

    tag_string -> '@' tags

    bases -> run_length_encoded_bases
          -> bases

    run_length_encoded_bases -> base count run_length_encoded_bases
                             -> base count
    
    count -> integer > 0
    
    base -> alphabet character or '-' or '*' ([A-Z,a-z,-,*])

    bases -> string or alphabet characters or '-' or '*' ([A-Z,a-z,-,*]+)

(Aligned bases are either a run-length encoded
representation or simply as a sequence of aligned characters. To specify
which format to use we use the tag "run_length_encode_bases:1" in the header.
If "run_length_encode_bases:0" or the tag is not specified the format is
to NOT use run length encoding.)

    coordinates_column -> ';' coordinates
    
    coordinates -> coordinate_operation coordinates
                -> ''
    
    coordinate_operation -> 'i' row coordinate
                         -> 'd' row
                         -> 's' row coordinate
                         -> 'g' row gap_length
                         -> 'G' row gap_string

(The 'i' stands for insertion, the 'd' for deletion, the 's' for substitution, 'g' and 'G' for gap. These operations
allow us to update the coordinates of the sequences as we go, and work as their name suggests. Rows are indexed from zero.
The operations are affected in order, so, for example, inserting a row at position i will shift all remaining column indices by one place
for any operations that are specified after an insertion.
The substitute operation can be used to change the coordinates for a sequence in a row, or it can be used to periodically repeat
coordinates to prevent needing to scan back more than N rows to find the coordinates of a row.
Using the 'g' gap operation can be used to increment the coordinate of a sequence by a specified amount and is useful for inserting
long substrings. Using the 'G' instead of 'g' allows one to instead specify the gap substring.)

    row -> integer >= 0
    
    gap_length -> integer >= 0
    
    gap_string -> alphabet string ([A-Z,a-z]*)
    
    coordinate -> sequence_name offset strand sequence_length
    
    sequence_name -> string without white space
    
    offset -> integer >= 0
    
    strand -> '+'
           -> '-'

    sequence_length -> integer >= 0

Coordinates in TAF use the same conventions as MAF format. That is,
zero-based, half-open coordinates, with a negative strand meaning that 
coordinates are with respect to the reverse complement sequence. The 
inclusion of the sequence_length in the coordinate, while bloating,
does make it easy to convert back and forth between maf and taf.

# Taf Example

The following shows by example the translation between the MAF format and TAF format.

The MAF file (602 bytes):

    ##maf version=1 scoring=N/A
    
    a
    s       simDog.chr6     437451  11      +       593897  CCCGTCAGTGT
    s       simHuman.chr6   446327  11      +       601863  TCCGCCAAGGT
    s       simMouse.chr6   460751  11      +       636262  TTCATCAGAGT
    s       simRat.chr6     470339  11      +       647215  TTCATTAGGGT
    
    a
    s       simCow.chr6     445326  8       +       602619  TTTTCCCA
    s       simDog.chr6     437462  8       +       593897  TT-TTCCG
    s       simHuman.chr6   446338  8       +       601863  TTCTTCCG
    s       simMouse.chr6   460762  8       +       636262  TTTTACCG
    s       simRat.chr6     470355  8       +       647215  TTTTACCG

The corresponding TAF file (265 bytes):

    #taf version:1 scoring:N/A
    CTTT ; i 0 simDog.chr6 437451 + i 1 simHuman.chr6 446327 + i 2 simMouse.chr6 460751 + i 3 simRat.chr6 470339 11
    CCTT
    CCCC
    GGAA
    CCCT
    AAAA
    GAGG
    TGAG
    GGGG
    TTTT
    TTTTT ; i 0 simCow.chr6 445326 + g 4 5
    TTTTT
    T-CTT
    TTTTT
    CTTAA
    CCCCC
    CCCCC
    AGGGG

# Installing Taffy CLI/C Library

Do build this repo clone the repo as follows and then make:

    git clone https://github.com/ComparativeGenomicsToolkit/taffy.git --recursive
    cd taffy && make

To test the installation do:

    make test

This will run the unitests. You should see that all tests pass okay. You will 
then want to add the taf/bin directory to your path. 

# Taffy Utilities

All Taffy utilities are run using `taffy <command>`, where the available commands are:

```
    view           MAF / TAF conversion and region extraction
    norm           normalize TAF blocks
    add-gap-bases  add sequences from HAL or FASTA files into TAF gaps
    index          create a .tai index (required for region extraction)
    sort           sort the rows of a TAF file to a desired order           
    stats          print statistics of a TAF file
```

Taffy supports both uncompressed and [bgzipped](http://www.htslib.org/doc/bgzip.html) input (though Taffy must be built
with [htslib](http://www.htslib.org/) for bgzip support).

## Taffy View

For example, to convert a maf file to a taf use:

    taffy view -i MAF_FILE

To go back to maf use:

    taffy view -i TAF_FILE -m

Taffy, view also has a number of nice options, for example to display only mismatches to the chosen reference sequence (first row of the file), or to only show inferred lineage mutations. For example, if this is an input alignment block (in MAF format):

```asm
a
s	Anc0.Anc0refChr0	0	50	+	4151	GTCAAGCTCAGTAGATACTGGATTAGGAATTCATGAGTTAAGCTGTAGCC
s	Anc1.Anc1refChr1	292714	50	+	296994	GTCAAGCTCAGTAGATACTGGATTAGGAATTCATGAGTTAAGCTGTAGCC
s	Anc2.Anc2refChr1	5	50	+	4655	GTCAAGCTCAGTTGATGCTGGATTAGGAATTCATGAGTTAAGCTGTAGTC
s	mr.mrrefChr1	178277	50	+	182340	GTCAAGCTCTGTACATACTAGATTGGACATTCATGGATGAAACTGTGACT
s	simCow_chr6.simCow.chr6	5045	50	-	602619	GTGAAGCTCAGTTGATGCTGGATTGGGAACTCATGAGTTAAGCTGTAAGC
s	simDog_chr6.simDog.chr6	589129	50	+	593897	GTCAAGCTCAGTTGGTGCTGGATTAAGAATTCATGAGTTAGGCTGCAGTC
s	simHuman_chr6.simHuman.chr6	597375	50	+	601863	GTCAAGCTCAGTAGATATTGGATTAGGAATTCATAAGTTAACCTGTAGCC
s	simMouse_chr6.simMouse.chr6	630640	50	+	636262	GTCAAGCATTGTACATACTAGATTGGACATTCATGGATGACAATGTGACT
s	simRat_chr6.simRat.chr6	642153	50	+	647215	GTCAAGCTCTGTAAATAGTAGATTGGACATTCATGGATGAAACTGTGCCT
```

Then running 

    taffy view -i TAF_FILE -m -a

Using -a option to only show differences to the reference (first alignemnt row), will result in:

```asm
a
s	Anc0.Anc0refChr0	0	50	+	4151	GTCAAGCTCAGTAGATACTGGATTAGGAATTCATGAGTTAAGCTGTAGCC
s	Anc1.Anc1refChr1	292714	50	+	296994	**************************************************
s	Anc2.Anc2refChr1	5	50	+	4655	************T***G*******************************T*
s	mr.mrrefChr1	178277	50	+	182340	*********T***C*****A****G*AC*******GA*G**A****GA*T
s	simCow_chr6.simCow.chr6	5045	50	-	602619	**G*********T***G*******G****C*****************AG*
s	simDog_chr6.simDog.chr6	589129	50	+	593897	************T*G*G********A**************G****C**T*
s	simHuman_chr6.simHuman.chr6	597375	50	+	601863	*****************T****************A******C********
s	simMouse_chr6.simMouse.chr6	630640	50	+	636262	*******ATT***C*****A****G*AC*******GA*G*CAA***GA*T
s	simRat_chr6.simRat.chr6	642153	50	+	647215	*********T***A***G*A****G*AC*******GA*G**A****GC*T
```

Similarly running

    taffy view -i TAF_FILE -m -b -t NEWICK_TREE_STRING

Will result in showing only substitutions inferred along the lineages of the tree:

```asm
a
s	Anc0.Anc0refChr0	0	50	+	4151	GTCAAGCTCAGTAGATACTGGATTAGGAATTCATGAGTTAAGCTGTAGCC
s	Anc1.Anc1refChr1	292714	50	+	296994	**************************************************
s	Anc2.Anc2refChr1	5	50	+	4655	************T***G*******************************T*
s	mr.mrrefChr1	178277	50	+	182340	*********T***C*****A****G*AC*******GA*G**A****GA*T
s	simCow_chr6.simCow.chr6	5045	50	-	602619	**G*********************G****C*****************AG*
s	simDog_chr6.simDog.chr6	589129	50	+	593897	**************G**********A**************G****C****
s	simHuman_chr6.simHuman.chr6	597375	50	+	601863	*****************T****************A******C********
s	simMouse_chr6.simMouse.chr6	630640	50	+	636262	*******AT*******************************C*A*******
s	simRat_chr6.simRat.chr6	642153	50	+	647215	*************A***G*****************************C**
```

Note this requires specifying ancestor sequences in the tree and having them correspond to the names of the sequences in the input tree.

## Taffy Add-Gap-Bases

There is also a utility for adding sequences between blocks to a taf file

    taffy add-gap-bases SEQ_FILES -i TAF_FILE

## Taffy Norm

There is also a utility to merge together short alignment blocks to create a more
"normalized" maf/taf file:
 
    taffy norm

For example, to normalize a maf file do the following:

    taffy view -i MAF_FILE | taffy norm -k -b SEQUENCE_FILES out.maf

`taffy view` converts MAF_FILE into taf, `taffy norm` then merges together blocks. The 
`-k` option causes the output to be in maf format. The `-b`
 option reads in underlying sequence files and is used to 
retrieve any sequences that are unaligned between two blocks that is necessary to include in stitching together adjacent blocks. This uses the same method as `taffy add-gap-bases` to add these unaligned sequences.

## Taffy Sort

It can be useful to sort the rows of an alignment. For this we have `taffy sort`. For example:

    taffy view -i MAF_FILE | taffy sort -n SORT_FILE | taffy view -m

Taffy view first converts the input maf to TAF, taffy sort then sorts the rows of the taf according to a given file and finally the last taffy view pipes the output back to maf. The sort file is a sequence of sequence name prefixes, so that each sequence name in the alignment has, uniquely, one of the strings in the sort file as a prefix. The order of these prefixes is then used to sort the rows. Taffy sort also supports filtering, to remove selected rows of an alignment. For example:

    taffy view -i MAF_FILE | taffy sort -f FILTER_FILE -n SORT_FILE | taffy view -m

Where here we additionally remove any rows with sequence names with a prefix contained in the given FILTER_FILE.

# Referenced-based MAF/TAF and Indexing

Neither format specification requires it, but *in practice* TAF, like MAF, is used to specify alignments
along a single reference genome.  For instance, MAF files from both MultiZ and Cactus (via hal2maf)
are **referenced-based**, and have the following properties:

* The first row of each alignment block consists of a sequence from the reference genome.
* The first row is on the forward (+) strand
* The alignment blocks are ordered (in increasing order) by the coordinates of the sequences on the first row.
* The intervals on the first row of each block do not overlap

An **anchor line** in TAF is a column from which all sequence coordinates can be deduced without scanning
backwards to previous lines in the file. In other words, each coordinate must be specified with either
an `s` or `i` operation. If a TAF is reference-based and the lowest coordinate of each first-row
sequence is specified on an anchor line, then the TAF file is **indexable**.

A TAF file produced from a reference-based MAF file using `taf view` will be reference-based and indexable.

Indexable TAF files and reference-based MAF files can be indexed for random-access using `taffy index`:

    taffy index -i TAF_FILE (or MAF_FILE)

This command will create `TAF_FILE.tai`, which is a list mapping sequence names (first column) and
start positions (second column) to offsets in the TAF file (third column). If two consecutive
lines index the same sequence, the second one will be relative to the first (and the first column
will be `*`). The `-b` option specifies the interval length in the index.  But only anchor
lines (whose frequency is controlled by `taffy view -s` can ever be indexed.  Smaller index
intervals will result in faster lookup times at the cost of the index itself being slower
to load.

An indexed TAF or MAF file can be accessed using `taffy view -r` to quickly pull out a subregion. For
example, `taffy view -r hg38.chr10:550000-600000` will extract the 50000bp (0-based, open-ended)
interval on `hg38.chr10` in either TAF (default) or MAF (add `-m`) format. This works only if
the TAF/MAF is referenced on hg38.

Notes:

* The sequence names do not need to be ordered for TAF/MAF indexing (ie chr2 could come before chr1
in the file). Just the positions within each sequence must be in order

* The index could further be generalized to support out-of-order (but still non-overlapping!)
intervals without (I think) any changes to the interface or format. Instead of the
simple seek + extend approach used when querying now, it would need to recheck and
potentially re-seek after extending though any subsequent intervals.  This could potentially
allow indexing non-reference intervals, but the effiency will degrade with the number of
out-of-order intervals. 

# Using C Library

There is also a simple C library for working with taf/maf files. See taf.h in the
inc directory.

# Installing Python Library

Tl;dr: 

```
pip install taffy
```

Longer version: to avoid problems with conflicting versions of dependencies on your system, we strongly recommend installing
the package inside a Python 3 [virtual environment](https://virtualenv.pypa.io/en/stable/). 
To install the `virtualenv` command, if you don't have it already, run:

```
python3 -m pip install virtualenv
```

To set up a virtual environment in the directory `taffy_env`, run:

```
python3 -m virtualenv -p python3.XX taffy_env
```

Where XX is the specific version of Python3 that you wish to use (you can omit the .XX if you want to use the default). Also, note that I have tested this with 3.9.

Then, to enter the virtualenv, run:

```
source taffy_env/bin/activate
```

You can always exit out of the virtualenv by running `deactivate`.
Finally, install taffy with pip:

To install these notebooks in Python, clone the repo:
```
pip install taffy
```

To check that it worked try launching a python interpretor and importing
taffy.lib.

# Installing Python Library From Source

To build the Python library you must first install [htslib](http://www.htslib.org/) for bgzip support.

If you are building the C library from source you can also build the Python library as follows:

```
git clone https://github.com/benedictpaten/taf.git --recursive
cd taf && make test
```

This will build and the test the C installation. Then create a virtualenv
if you haven't already:

```
python3 -m pip install virtualenv
python3 -m virtualenv -p python3.XX taffy_env
source taffy_env/bin/activate
```

You can then build from source and test the distribution by running:

```
python3 -m pip install build
python3 -m build
pip install .
cd tests && python3 taffyTest.py
```

All these commands should succeed. Note, because of the current package structure
trying to import the library from the root directory will fail.

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
with AlignmentReader(test_taf_file, taf_index=taf_index, sequence_name="Anc0.Anc0refChr0",start=1000,length=50) as mp:
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

Which gets a particular subrange of blocks within the given reference sequence interval.

# Comparison Stats

Using the file:  https://hgwdev.gi.ucsc.edu/~markd/cactus/cactus241way/ucscNames/chr3_KI270777v1_alt.maf
(A randomly chosen alignment part of the Cactus 241-way alignment).
Running the command:

    taffy view -i ./chr3_KI270777v1_alt.maf > ./chr3_KI270777v1_alt.taf

Took: 
    
    22.29s user 0.34s system 99% cpu 22.685 total

And results in:

    -rw-r--r--   1 benedictpaten  staff  753714843 Dec 15  2020 chr3_KI270777v1_alt.maf
    -rw-r--r--   1 benedictpaten  staff   50840928 Sep 15 11:28 chr3_KI270777v1_alt.taf

That is a 14.8x compression. Gzipped:

    -rw-r--r--   1 benedictpaten  staff  89574699 Dec 15  2020 chr3_KI270777v1_alt.maf.gz
    -rw-r--r--   1 benedictpaten  staff  14038446 Sep 15 11:28 chr3_KI270777v1_alt.taf.gz

The .taf.gzip is 6.39x smaller than the .maf.gz and 53x smaller than the .maf
The .taf is 1.9x smaller than the .maf.gz.

Normalizing the maf file using the command:

    taffy view -i ./chr3_KI270777v1_alt.maf | taffy norm -k > ./chr3_KI270777v1_alt.norm.maf

Took: 

    23.60s user 0.23s system 97% cpu 24.460 total

Resulting in:

    -rw-r--r--   1 benedictpaten  staff  753714843 Dec 15  2020 chr3_KI270777v1_alt.maf
    -rw-r--r--   1 benedictpaten  staff  156736026 Sep 15 11:35 chr3_KI270777v1_alt.norm.maf
    -rw-r--r--   1 benedictpaten  staff   50840928 Sep 15 11:28 chr3_KI270777v1_alt.taf

So the normalized maf is 4.8x smaller than the maf and 3x larger than the taf.
The number of blocks in the maf file is also reduced as blocks are merged together by the 
normalization process:

    taf % grep '^a' ./chr3_KI270777v1_alt.maf | wc -l                                                                    
        45924
    taf % grep '^a' ./chr3_KI270777v1_alt.norm.maf | wc -l
        7157

Which is a 6.97x reduction in block number.

# TODOs

Things that are ongoing:

* Make `taffy add-gap-bases` use indexed fastas to avoid loading everything into memory

* Add a sorting option to taffy to get the rows of an alignment in a sorted order


