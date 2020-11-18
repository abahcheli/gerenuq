# sam file filtering script
# Alec Bahcheli, Daniel Giguire from Gloor Lab, Western University, Canada

import sys, getopt, re, time, math
import concurrent.futures

# error code to return without necessary input
error_code = '''
gerenuq.py

Required inputs:
sam = "NULL" <input raw samfile>
results_file = "NULL" <output filtered samfile>

Optional inputs:
mls = 2 <sequence identity, also known as minimum ratio of matches to read length (default 0.5)>
mml = 0.5 <minimum ratio of length to score, may be considered as the fraction of bases that have a positive score (default 0.5)>
ml = 1000 <minimum read length for cutoff (default 1000)>
ms = 1 <minimum score for the whole alignment (default 1)>
wpc = 1 <number of processes to run (default 1)>

version 0.0.1'''

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

def it_meets_filters(length, num_of_matches):
    if int(length) > min_length and (int(num_of_matches) / int(length)) > min_match_to_length:
        return True
    else:
        return False

def it_is_good_score(length, score):
    if int(score) > min_score and (int(length) / int(score)) < min_len_to_score:
        return True
    else:
        return False

def filter_reads(read):
    # split the read into the expected fields
    read = read.split("\t")
    if int(read[1]) == 16 or int(read[1]) == 0:
        # read mapping score
        score = int(read[13][5:])
        # read length 
        length = 0
        # number of matches
        num_of_matches = 0
        # get the cigar string
        cigar = re.findall("([0-9]*[MISH])", read[5])
        for element in cigar:
            # if the read doesn't match with the chromosome, the length of alignment increase but not the matches
            if re.search("[ISH]", element):
                length += int(element.strip("[ISH]")) 
            # if the read at a position does match the chromosome sequence, the length and the number of matches increase
            elif re.search("M", element):
                length += int(element.strip("M")) 
                num_of_matches += int(element.strip("M"))
        # if it meets the filter cutoffs, return the whole read
        if it_is_good_score(length, score):
            if it_meets_filters(length, num_of_matches):
                return "\t".join(read)


def main(sam = "NULL", results_file = "NULL", mls = 0.5, mml = 0.5, ml = 1000, ms = 1, wpc = 1):
    t1 = time.time()
    
    # number of processes
    worker_process_count = wpc

    # minimum ratio of length to score
    min_len_to_score = mls

    # minimum ratio of the number of matches to the length
    min_match_to_length = mml

    # minimum length for a read to be considered
    min_length = ml

    # minimum score for an alignment to be considered
    min_score = ms

    print(error_code)

    # test if the minimmum input parameters are defined
    if sam == "NULL" or results_file == "NULL":
        return print(error_code)
    
    # open the samfile
    samfile_raw = open(sam)

    # open the results file
    results = open(results_file, "w")

    # make a list of just reads
    samfile = []

    # get the headers
    for read in samfile_raw:
        if read.startswith("@"):
            results.write(read)
        else:
            samfile.append(read)

    print("Samfile read into memory for parallelization")
    print(time.time() - t1)

    # chunking the reads improves processing time by avoiding compilation congestion at the end
    if worker_process_count > 1:
        chunks = math.floor(0.2 * len(samfile) / worker_process_count)
    else:
        chunks = 1

    # list of reads that satisfy the cutoff requirements
    reads_of_interest = []

    # parallelize read evaluations
    with concurrent.futures.ProcessPoolExecutor(max_workers = worker_process_count) as executor:
        for reads in executor.map(filter_reads, samfile, chunksize = chunks):
            if reads != None:
                reads_of_interest.append(reads)

    print("Done filtering reads")
    print(time.time() - t1)

    for read in reads_of_interest:
        # write the good reads to a file
        results.write(read)
    
    results.close()

    print("Finished writing filtered samfile")
    print(time.time() - t1)
