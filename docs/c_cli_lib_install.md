# Installing Taffy CLI/C Library

Do build this repo clone the repo as follows and then make:

    git clone https://github.com/ComparativeGenomicsToolkit/taffy.git --recursive
    cd taffy && make

To test the installation do:

    make test

This will run the unit tests. You should see that all tests pass okay. You will
then want to add the taf/bin directory to your path.

## htslib dependency

Taffy uses htslib (via `pkg-config --exists htslib` at build time) for bgzip
support and HTTP/S3 input. If `pkg-config` doesn't find an htslib install,
the build still succeeds but bgzipped input and URL input are disabled.

You can also point at a non-pkg-config htslib install by setting
`HTSLIB_CFLAGS` and `HTSLIB_LIBS` in the environment before `make`.

### URL/HTTPS/S3 input requires htslib with libcurl

`taffy view -r SEQ:s-e -i https://...` (or `s3://...`, etc.) routes through
htslib's libcurl-backed `hFILE` layer. Most distro/conda htslib packages
have this enabled. If you built htslib from source, you may need to rerun
configure with libcurl support:

    cd htslib
    autoreconf -i           # only if there's no ./configure yet
    ./configure --enable-libcurl
    make

If your htslib lacks libcurl, taffy will print a clear error suggesting the
above when a URL input is given.
