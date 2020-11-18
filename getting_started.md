# gerenuq
### Background and Theory

gerenuq.py is based off of a series of commands used to filter reads, originating from the filtering process used in the cmags paper. Instead of requiring a number of inputs and outputs, this script is a single line requiring a samfile input and returning a filtered samfile list.

The script filters reads mapped against a reference from a minimap2 results samfile. Required input parameters is a samfile (-i or --samfile) (see __Getting Started__) and an output file (-o or --output). 

The script will parse the samfile, filtering reads that are primary alignments, at least 1,000 bases long, meet a minimum ratio of 0.5 for the number of matches to the read length (sequence identity), and be less than a maximum ratio of 2 for the length divided by the score (inverse of the average score per base).

Optional parameters can change the filters in a number of ways (refer to the help command when running the script). It is highly recommended that you multi-thread to speed up the filtering process.

### Getting Started

For appropriate inputs, type ```python3 gerenuq.py --help```.

The samfile should be the output from a minimap2 alignment that may be filtered by samtools.

Required input parameters is a samfile (-i or --input) and output file (-o or --output). Output will a samfile filtered according to input or default parameters. For example, a simple input would be:

```
python3 cigar-parse_phased.py -i raw_samfile.sam -o filtered_samfile.sam
```

Processing time increases exponentially for each additionally mapped read. Multi-processing is recommended, which the number of processes to run can be describe with the *-t* flag (or *--threads*).

```
python3 cigar-parse_phased.py -i raw_samfile.sam -o filtered_samfile.sam --length 50000 --threads 20
```
