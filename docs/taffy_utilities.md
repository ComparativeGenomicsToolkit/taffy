# Taffy Utilities

All Taffy utilities are run using `taffy <command>`, where the available commands are:

```
    view           MAF / TAF conversion and region extraction
    norm           normalize TAF blocks
    add-gap-bases  add sequences from HAL or FASTA files into TAF gaps
    index          create a .tai index (required for region extraction)
    sort           sort the rows of a TAF file to a desired order           
    stats          print statistics of a TAF file
    coverage       print coverage statistics of a given genome in a TAF file
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

Note this requires specifying ancestor sequences in the tree and having them correspond uniquely to the prefixes of names of the sequences in the input tree.

## Taffy Add-Gap-Bases

There is also a utility for adding sequences between blocks to a taf file

    taffy add-gap-bases SEQ_FILES -i TAF_FILE

## Taffy Norm

There is also a utility to merge together short alignment blocks to create a more
"normalized" maf/taf file:

    taffy norm

For example, to normalize a maf file do the following:

    taffy view -i MAF_FILE | taffy norm -k -b SEQUENCE_FILES -o out.maf

`taffy view` converts MAF_FILE into taf, `taffy norm` then merges together blocks. The
`-k` option causes the output to be in maf format. The `-b`
option reads in underlying sequence files and is used to
retrieve any sequences that are unaligned between two blocks that is necessary to include in stitching together adjacent blocks. This uses the same method as `taffy add-gap-bases` to add these unaligned sequences.

Note(!), taffy norm will resort the rows alpha-numerically according to sequence name,  as is necessary to successfully merge all mergeable rows. Is the resorting is undesired, pipe the result to taffy sort (see below) to resort, e.g.

    taffy view -i MAF_FILE | taffy norm -b SEQUENCE_FILES | taffy sort -n SORT_FILE | taffy view -m

## Taffy Sort

It can be useful to sort the rows of an alignment. For this we have `taffy sort`. For example:

    taffy view -i MAF_FILE | taffy sort -n SORT_FILE | taffy view -m

Taffy view first converts the input maf to TAF, taffy sort then sorts the rows of the taf according to a given file and finally the last taffy view pipes the output back to maf. The sort file is a sequence of sequence name prefixes, so that each sequence name in the alignment has, uniquely, one of the strings in the sort file as a prefix. The order of these prefixes is then used to sort the rows. Taffy sort also supports filtering, to remove selected rows of an alignment. For example:

    taffy view -i MAF_FILE | taffy sort -f FILTER_FILE -n SORT_FILE | taffy view -m

Where here we additionally remove any rows with sequence names with a prefix contained in the given FILTER_FILE.

A common requirement is that a MAF/TAF has exactly one row for every species in the alignment. As Cactus, and other tools, output MAFs that don't necessarily follow this convention it is useful to be able to force this. Using:

    taffy view -i MAF_FILE | taffy sort -n SORT_FILE -p PAD_FILE -r DUP_FILE | taffy view -m

Where the -p specifies the prefixes to "pad", that is any block not containing a row matching a prefix
in the PAD_FILE will have that row added, using gaps and dummy coordinates to fill in the row. Similarly, the -r specifies that any set of two or more rows whose
names match a given prefix in the DUP_FILE will be pruned so that only one such row is kept in the block. The heuristic used for dropping dupes currently is intentionally very simple: all rows after the first occurrence of a row matching the given sequence prefix are dropped. Using these options (and optionally the filter option) allows you to construct a MAF ordered and with exactly the set of rows expected for every block.

In the taffy/scripts directory are some useful utilities for creating the sort/pad/dup-filter files given a guide tree. For example:

    ./scripts/tree_to_sort_file.py --traversal pre --reroot REF_NODE --out_file OUT_FILE NEWICK_TREE_FILE

Will create a sort order based upon a pre-order traversal of the  tree after rerooting the tree so that the given REF_NODE is the reference - this will place REF_NODE first in the sort order and then order remaining nodes from that node in a pre-order traversal.
Using options to exclude internal nodes or leaf nodes makes it easy
to use this to only include leaves, or to create a filter to exclude
internal nodes, say. For an example of usage see ./tests/447-way/example_norm.sh

## Taffy Coverage

This tool, `taffy coverage`, calculates basic coverage and percent identity statistics of a selected reference (defaults to first row) vs all other genomes in the alignment. Whole-genome and reference-contig-level statistics are provided. The values presented are:

* `ref-contig`: full name of reference contig (`_Total_` for whole-genome numbers)
* `max-gap`: reference bases in alignment gaps greater than this value are not counted in `len`.
* `len`: length of `ref-contig`
* `query`: name of query genome
* `aln-bp` : number of (non-`N`) bases in `ref-contig` that align to a (non-`N`) base in `query`
* `ident-bp` : number (non-`N`) bases in `ref-contig` that align to *THE SAME* (non-`N`) base in `query` at least once
* `1:1-aln-bp` : number of (non-`N`) bases in `ref-contig` that align to no other positions in `reference` genome and that align to exactly one (non-`N`) base in `query`
* `1:1-ident-bp` : number of (non-`N`) bases in `ref-contig` that align to no other positions in `reference` genome and that align to exactly one *OF THE SAME* (non-`N`) bases in `query`
* `aln` : `aln-bp / len`
* `1:1-aln` : `1:1-aln-bp / len`
* `ident` : `ident-bp / aln-bp`
* `1:1-ident` : `1:1-ident-bp / 1:1-aln-bp`

The coverage is broken down into overall statistics and just those corresponding to 1:1 alignments (ie where both genomes appear only once in the block). In the event of a column with multiple copies of the same genome, it will count as identical if any one of those copies matches the reference base. Total statistics as well as reference contig breakdowns are output. The percent identity

By default, the first `.` character is used to parse out the genome name from a sequence name.  So `hs1.chr1` would belong to genome `hs1`.  This does not work if the genome name itself contains a `.` character. In this case, it is best to supply the genome names with the `-g` flag, and they will be used to help the parser.  For example

    taffy view -i MAF_FILE | taffy coverage -g "$(halStats --genomes HAL_FILE)" > COV.tsv

The `-a` option can be used to add rows that ignore gaps greater than the specified size when computing coverage.  So `-a 10 -a 100` would report coverage statistics for the whole genome, as well as ignoring gaps `>10bp` and `>100bp`. There will be `3X` the number of output rows.

You can also use the `-s` option to add a breakdown of sex chromosomes and autosomes to the output table, ex `-s chrX -s chrY`.

## Taffy Stats

`taffy stats` can print some basic statistics about the reference contigs (first row) in the alignment. Options are
* List the reference contigs and their lengths (maximum end coordinates as found in the alignment -- could be smaller than true contig lengths). This is done quickly using the `.tai` index created with `taffy index`.
* List all the reference intervals in BED format found in the alignment.  This is done by scanning the whole alignment and is therefore much slower than the above option (but will produce more fine-grained output in the case where the alignment only covers subregions of the contigs).
  In addition, it is useful for getting course data about a MAF/TAF, e.g.:

```
taffy stats -i ./447-way/447-mammalian-2022v1_hg38_chr22_22000000_22100000.anc.norm.taf.gz -a
Total blocks:   3166
Total columns:  215664
Avg. columns/block:     68.118759
Total bases:    67026589
Total gaps:     226206389
Avg. column depth:      1359.675171
Avg. bases/column:      310.791718
Avg. gaps/column:       1048.883423
```

Note the -a option is required to print these aggregate stats.

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