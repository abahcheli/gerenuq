# GERENUQ

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4281015.svg)](https://doi.org/10.5281/zenodo.4281015)

A simple commandline tool for filtering long reads from samfiles according to various user-defined parameters.
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
```bash
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
