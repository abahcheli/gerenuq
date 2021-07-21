# GERENUQ

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4281015.svg)](https://doi.org/10.5281/zenodo.4281015)

A simple commandline tool and python functions for filtering long reads from bam, sam and paf red alignment files according to various user-defined parameters.
# Installation
## Using Conda
```bash
  $ conda install -c conda-forge -c bioconda -c abahcheli gerenuq
```
## Using Pip
```bash
  $ pip install gerenuq
```
## Using Docker
```bash
  $ docker pull abahcheli/gerenuq
```
## Manual
```bash
  $ git clone https://github.com/abahcheli/gerenuq
  $ cd gerenuq
  $ python setup.py install
```
# Usage
```sh
gerenuq

$ Required inputs:
-i / --input <input raw samfile>
-o / --output <output filtered samfile>

Optional inputs:
-l / --length <minimum read length for cutoff (default 1000)>
-m / --matchlength <sequence identity, also known as minimum ratio of matches to read length (default 0.5)>
-s / --score <minimum score for the whole alignment (default 1)>
-q / --lengthscore <minimum ratio of length to score, may be considered as the fraction of bases that have a positive score (default 2)>
-t / --threads <number of processes to run (default 1)>
```

```python
gerenuq_filter_file(input_file, output_file, min_score = 1, min_len_to_score = 2, min_length = 1000, min_match_to_length = 0.5)
'''
Filters minimap2-mapped reads by mapping score, length, match-to-length and length-to-score ratios. Paf format files only filter by query cutoff.

Requires input_file in bam, sam or paf format and output_file (output in the same format as input).
'''

gerenuq_filter_read_list(read_list, format='sam', min_score = 1, min_len_to_score = 2, min_length = 1000, min_match_to_length = 0.5)
'''
Filters minimap2-mapped reads by mapping score, length, match-to-length and length-to-score ratios. Paf format files only filter by query cutoff.

Requires read_list as list of mapped read lines from sam or paf file (in tsv format). Returns a list of reads in sam or paf (tsv) format that passed filtering parameters. Headers will be ignored and not returned.
'''

```

# Getting Started
## Background and Theory
gerenuq.py is based off of a series of commands used to filter reads, originating from the filtering process used in the cmags paper. Instead of requiring a number of inputs and outputs, this script is a single line requiring a samfile input and returning a filtered samfile list.

The script filters reads mapped against a reference from a minimap2 results samfile. Required input parameters is a samfile (-i or --samfile) (see __Getting Started__) and an output file (-o or --output). 

The script will parse the samfile, filtering reads that are primary alignments, at least 1,000 bases long, meet a minimum ratio of 0.5 for the number of matches to the read length (sequence identity), and be less than a maximum ratio of 2 for the length divided by the score (inverse of the average score per base).

Optional parameters can change the filters in a number of ways (refer to the help command when running the script). It is highly recommended that you multi-thread to speed up the filtering process.

## Quick Start
For appropriate inputs, type ```python3 gerenuq.py --help```.

The samfile should be the output from a minimap2 alignment that may be filtered by samtools.

Required input parameters is a samfile (-i or --input) and output file (-o or --output). Output will a samfile filtered according to input or default parameters. For example, a simple input would be:

```
gerenuq -i raw_samfile.sam -o filtered_samfile.sam
```

Processing time increases exponentially for each additionally mapped read. Multi-processing is recommended, which the number of processes to run can be describe with the *-t* flag (or *--threads*).

```
gerenuq -i raw_samfile.sam -o filtered_samfile.sam --length 50000 --threads 20
```

