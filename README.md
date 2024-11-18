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

Simply supply the fastq file path, the upstream and downstream flanking regions, and other optional parameters. Run

    get_plasmid_inserts -h

for help and to see optional arguments. A good flanking sequence length is
around 30-50 bp. The default distance threshold for matching is 5% of the
length of the shortest flanking sequence; depending on the quality of your data,
this may need to be increased.

The program will output two files:  
1. A trimmed fastq of all inserts found. Inserts are oriented to match the direction
of the flanking sequences used, and the indices of the insert in the original full read
(after reverse complementing if this was done) and whether the read was reverse
complemented are added to the header.
2. A text file with three columns: upstream matched sequence, insert sequence,
and downstream matched sequence. This can be useful if you want to get an idea of whether
the regions being matched are what you expect, and if you want to visually compare insert
sequences easily as they are all aligned in a column. Optionally, quality can be included with
the `-q` argument, the reverse complement of the inserts can be added with the `-r` argument,
and empty inserts can be excluded with the `-n` argument.

## Running tests

Run the following from the top level of the directory to run all tests:
```
python -m unittest discover tests
```