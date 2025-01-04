# Taffy

This is a MIT licensed C and Python library with a CLI for manipulating/reading/writing [TAF](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taf_format.md) (described below) and 
[MAF](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) format multiple
sequence alignments. It allows conversion between the formats and manipulation of the alignments with a number of useful utilities for preparing them for different use cases. The Python library is built
on top of the C library and is therefore quite fast.

## Taf Format Specification

See the [Taf format page](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taf_format.md) for a specification of the taf format and example.

## Installation

See [C/CLI Install](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/c_cli_lib_install.md) for how to build and install this source for using the C library and CLI utilities.

See [Python install](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/py_install.md) for how to install the Python library.

### CLI Utilities

See [taffy utilities](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taffy_utilities.md) for a description of the many useful taffy utilities, including:

 * [view](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taffy_utilities.md#taffy-view)     -      MAF / TAF conversion and region extraction
 * [norm](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taffy_utilities.md#taffy-norm)     -      normalize TAF blocks 
 * [add-gap-bases](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taffy_utilities.md#taffy-add-gap-bases) - add sequences from HAL or FASTA files into TAF gaps
 * [index](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taffy_utilities.md#referenced-based-maftaf-and-indexing)   -       create a .tai index (required for region extraction)
 * [sort](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taffy_utilities.md#taffy-sort)    -       sort the rows of a TAF file to a desired order           
 * [stats](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taffy_utilities.md#taffy-stats)   -       print statistics of a TAF file
 * [coverage](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taffy_utilities.md#taffy-coverage) -      print coverage statistics of a given genome in a TAF file

### Python scripts

See [taffy scripts](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taffy_scripts.md) for a description of useful Python scripts, including:

* [taffy alignment plot](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/taffy_scripts.md#taffy-alignment-plot)     -      A (relatively) fast MSA visualization, with coverage, copy-number, identity, and dotplot options.

### Using the Python API

See [using the Python API](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/py_usage.md) for how to work with MAF/TAF alignments using a convenient Python API designed to complement the CLI.

See [the example notebook](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/examples/learning_phlyoP.ipynb) for a quick worked example of using the Python API for machine learning with PyTorch.

### Third-party tools

- JBrowse 2 with mafviewer plugin has TAF support https://github.com/cmdcolin/jbrowse-plugin-mafviewer

## Using the C Library

There is also a simple C library for working with taf/maf files. See taf.h in the
inc directory.

## Comparing MAF and TAF file sizes

See [quick file size comparison](https://github.com/ComparativeGenomicsToolkit/taffy/blob/main/docs/comparison_stats.md).
