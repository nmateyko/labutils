# Labutils
This repository contains miscellaneous scripts that I find myself repeatedly
using in my lab work.

Install:

    pip install git+https://github.com/nmateyko/labutils.git


## sequence_analysis
This package contains scripts and functions related to DNA sequence analysis.

### utils.py
Contains useful functions like sequence distance metrics, fastq parser, and reverse complement.

### get_plasmid_inserts.py

This script is meant to extract inserts from full plasmid nanopore sequencing reads.
It's very useful for QCing plasmid libraries; send the
library to Plasmidsaurus, then run this script on the raw reads fastq.

If labutils is installed as described above, this can be run from the command line
like so:

    get_plasmid_inserts <args>

Simply supply the fastq file path, the upstream and downstream flanking regions,
and maximum expected insert size. Run

    get_plasmid_inserts -h

for help and to see optional arguments. A good flanking sequence length is
around 30 bp for the default distance threshold of 8.

The function that extracts the inserts, `get_insert`, can be imported like so:

    from sequence_analysis.get_plasmid_inserts import get_insert

