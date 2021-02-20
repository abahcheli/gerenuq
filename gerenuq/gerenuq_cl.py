#! /usr/bin/env python
# Alec Bahcheli, Daniel Giguire from Gloor Lab, Western University, Canada

import sys, getopt, re, time, math, pysam
import concurrent.futures
import pandas as pd

# number of processes
worker_process_count = 1

# minimum ratio of length to score
min_len_to_score = 2

# minimum ratio of the number of matches to the length
min_match_to_length = 0.5

# minimum length for a read to be considered
min_length = 1000

# minimum score for an alignment to be considered
min_score = 1

# def filter_bamfile(read, min_score = 1, min_len_to_score = 2, min_length = 1000, min_match_to_length = 0.5):
def filter_bamfile(read):
    # (read[0], read[1], read[5], ((read[-1].split("), ("))[2]).split(", ")[1]))
    # split the read into the expected fields
    # read = read.split("\t")
    bitflag = read[1]
    cigar_string = read[2]
    alignment_score = read[3]
    if int(bitflag) == 16 or int(bitflag) == 0:
        # read mapping score
        score = int(alignment_score)
        # is the score high enough
        if score >= min_score:
            # read length 
            length = 0
            # number of matches
            num_of_matches = 0
            # get the cigar string
            cigar = re.findall("([0-9]*[MISH])", cigar_string)
            for element in cigar:
                # if the read doesn't match with the chromosome, the length of alignment increase but not the matches
                if re.search("[ISH]", element):
                    length += int(element.strip("[ISH]")) 
                # if the read at a position does match the chromosome sequence, the length and the number of matches increase
                elif re.search("M", element):
                    length += int(element.strip("M")) 
                    num_of_matches += int(element.strip("M"))
            # if it meets the filter cutoffs, return the whole read
            if bool(length > min_length and (num_of_matches / length) > min_match_to_length) and bool((int(length) / int(score)) < min_len_to_score):
                
                # return read name and cigar string
                return (read[0], read[2])

def filter_samfile(read):
    # split the read into the expected fields
    read = read.split("\t")
    bitflag = read[1]
    cigar_string = read[5]
    alignment_score = read[13]
    if int(bitflag) == 16 or int(bitflag) == 0:
        # read mapping score
        score = int(alignment_score[5:])
        # is the score high enough
        if score >= min_score:
            # read length 
            length = 0
            # number of matches
            num_of_matches = 0
            # get the cigar string
            cigar = re.findall("([0-9]*[MISH])", cigar_string)
            for element in cigar:
                # if the read doesn't match with the chromosome, the length of alignment increase but not the matches
                if re.search("[ISH]", element):
                    length += int(element.strip("[ISH]")) 
                # if the read at a position does match the chromosome sequence, the length and the number of matches increase
                elif re.search("M", element):
                    length += int(element.strip("M")) 
                    num_of_matches += int(element.strip("M"))
            # if it meets the filter cutoffs, return the whole read
            if bool(length > min_length and (num_of_matches / length) > min_match_to_length) and bool((int(length) / int(score)) < min_len_to_score):
                
                return ("\t".join(read))

def main():
    # version
    version = "version 0.2.0"

    t1 = time.time()

    # paf mode = False by default
    file_format = 'sam'

    # error code to return without necessary input
    error_code = '''
    gerenuq
    Filters minimap2-mapped reads by mapping score, length, match-to-length and length-to-score ratios. Paf format files only filter by query cutoff.

    Required inputs:
    -i / --input <input unfiltered sam, bam or paf file (auto-detected by file name)>
    -o / --output <output filtered sam, bam or paf file (auto-detected by file name)>

    Optional inputs:
    -l / --length <minimum read length for cutoff (default 1000)>
    -m / --matchlength <sequence identity, also known as minimum ratio of matches to read length (default 0.5)>
    -s / --score <minimum score for the whole alignment (default 1)>
    -q / --lengthscore <minimum ratio of length to score, may be considered as the fraction of bases that have a positive score (default 2)>
    -t / --threads <number of processes to run (default 1)>

    {vers}'''.format(vers=version)

    # get the options and files required for the input
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:l:m:s:q:t:v:",["input=","output=", "length=", "matchlength=", "score=", "lengthscore=", "threads=", "version="])
    except getopt.GetoptError:
        print (error_code)
        sys.exit()
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print (error_code)
            sys.exit()
        elif opt in ("-i", "--input"):
            input = str(arg)
            if re.search(r'\.[pP][aA][fF]', input):
                file_format = 'paf'
            elif re.search(r'\.[bB][aA][mM]', input):
                file_format = 'bam'
        elif opt in ("-o", "--output"):
            results_file = str(arg)
        elif opt in ("-l", "--length"):
            min_length = int(arg)
        elif opt in ("-m", "--matchlength"):
            min_match_to_length = float(arg)
        elif opt in ("-s", "--score"):
            min_score = int(arg)
        elif opt in ("-q", "--lengthscore"):
            min_len_to_score = float(arg)
        elif opt in ("-t", "--threads"):
            worker_process_count = int(arg)
        elif opt in ("-v", "--version"):
            print(version)

    # test if the minimum input parameters are defined
    try:
        input
        results_file
    except NameError:
        return (print(error_code))

    if file_format == 'paf': 
        df=pd.read_csv(input, sep="\t", header=None)
        
        print("File read into memory for parallelization")
        print(round(time.time() - t1, 2))
        # filter reads by query cutoff 
        filtered_reads = df[ (df[3]-df[2] ) / df[1] > min_match_to_length]
        
        print("Done filtering reads")
        print(round(time.time() - t1, 2))
        # output file
        filtered_reads.to_csv(results_file, sep = "\t", index=False, header=False)

        print("Finished writing filtered paf file")
        print(round(time.time() - t1, 2))

    elif file_format == 'bam':
        # make a list of just reads
        samfile = []
        infile = pysam.AlignmentFile(input, 'rb', threads=worker_process_count)
        for read in infile:
            # add in order: read name, alignment type (index 1), CIGAR string, and alignment score (final index after 'AS:i:') in that order
            # samfile.append(str(read))
            read = str(read).split("\t")
            samfile.append(list((read[0], read[1], read[5], ((read[-1].split("), ("))[2]).split(", ")[1])))
    
        print("File read into memory for parallelization")
        print(round(time.time() - t1, 2))
        # chunking the reads improves processing time by avoiding compilation congestion at the end
        if worker_process_count > 1:
            chunks = math.floor(0.2 * len(samfile) / worker_process_count)
        else:
            chunks = 1
        # parallelize read evaluations
        with concurrent.futures.ProcessPoolExecutor(max_workers = worker_process_count) as executor:
            # list of reads that satisfy the cutoff requirements
            # reads_of_interest = list(filter(None, executor.map(filter_bamfile, samfile, chunksize = chunks)))
            reads_of_interest = dict(list(filter(None, executor.map(filter_bamfile, samfile, chunksize = chunks))))
            
        print("Done filtering reads")
        print(round(time.time() - t1, 2))
        print(len(reads_of_interest))
        # compare to bam file by read name (dictionary key) and CIGAR score (dictionary value)
        infile = pysam.AlignmentFile(input, 'rb', threads=worker_process_count)
        bam_output = pysam.AlignmentFile(results_file, "wb", template=infile)
        # for read in reads_of_interest:
        #     bam_output.write(read)
        for read in infile:
            if read.query_name in reads_of_interest:
                if read.cigarstring == reads_of_interest.get(read.query_name):
                    bam_output.write(read)
        
        print("Finished writing filtered bam file")
        print(round(time.time() - t1, 2))

    else:
        # open the results file
        results = open(results_file, "w")
        # make a list of just reads
        samfile = []
        # get the headers
        with open(input) as input_raw:
            for read in input_raw:
                if read.startswith("@"):
                    results.write(read)
                else:
                    samfile.append(read) 

        print("File read into memory for parallelization")
        print(round(time.time() - t1, 2))
        # chunking the reads improves processing time by avoiding compilation congestion at the end
        if worker_process_count > 1:
            chunks = math.floor(0.2 * len(samfile) / worker_process_count)
        else:
            chunks = 1
        # parallelize read evaluations
        with concurrent.futures.ProcessPoolExecutor(max_workers = worker_process_count) as executor:
            # list of reads that satisfy the cutoff requirements
            reads_of_interest = list(filter(None, (executor.map(filter_samfile, samfile, chunksize = chunks))))
        
        print("Done filtering reads")
        print(round(time.time() - t1, 2))
        for read in reads_of_interest:
            # write the good reads to a file
            results.write(read)
        results.close()

        print("Finished writing filtered sam file")
        print(round(time.time() - t1, 2))

