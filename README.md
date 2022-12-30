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
* is very easy to index, as each column of the alignment is a single line, 
* it supports extensible column annotations.

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

(Columns are encoded one per line. Each column encodes an alignment of the
bases in the column (see below) and then some optional coordinates and tags.)

    tag_string -> '#' tags

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

    git clone https://github.com/benedictpaten/taf.git --recursive
    cd taf && make

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
```

Taffy supports both uncompressed and [bgzipped](http://www.htslib.org/doc/bgzip.html) input (though Taffy must be built
with [htslib](http://www.htslib.org/) for bgzip support).

For example, to convert a maf file to a taf use:

    taffy view -i MAF_FILE

To go back to maf use:

    taffy view -i TAF_FILE -m

There is also a utility for adding sequences between blocks to a taf file

    taffy add-gap-bases SEQ_FILES -i TAF_FILE

And finally, a utility to merge together short alignment blocks to create a more
"normalized" maf/taf file:
 
    taffy norm

For example, to normalize a maf file do the following:

    taffy view -i MAF_FILE | taffy add-gap-bases SEQUENCE_FILES | taffy norm -k > out.maf

`taffy view` converts MAF_FILE into taf, `taffy add-gap-bases` adds in missing
unaligned sequences between maf blocks and `taffy norm` then merges together the blocks. The 
`-k` option causes the output to be in maf format.

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
python3 -m virtualenv -p python3.9 taffy_env
```

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
python3 -m virtualenv -p python3.9 taffy_env
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

See: [TODO]

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
