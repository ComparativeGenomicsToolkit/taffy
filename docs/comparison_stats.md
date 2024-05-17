# Taf/Maf Comparison Stats

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