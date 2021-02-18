# GERENUQ

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4281015.svg)](https://doi.org/10.5281/zenodo.4281015)

A simple commandline tool and python functions for filtering long reads from bam, sam and paf red alignment files according to various user-defined parameters.
# Installation
## Using Conda
```bash
  $ conda install -c abahcheli gerenuq
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
