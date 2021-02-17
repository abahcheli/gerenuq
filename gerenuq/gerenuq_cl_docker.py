#! /usr/bin/env python
# sam file filtering script
# Alec Bahcheli, Daniel Giguire from Gloor Lab, Western University, Canada

import sys, getopt, re, time, math, os
import concurrent.futures

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

def main():
    # test if the minimmum input parameters are defined
    try:
        sam
        results_file
    except NameError:
        print(error_code)
        return False

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
    return True

if __name__ == '__main__':
    # version
    version = "version 0.0.1"

    t1 = time.time()

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

    # error code to return without necessary input
    error_code = '''
    gerenuq

    Required inputs:
    -i / --input <input raw samfile>
    -o / --output <output filtered samfile>

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
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print (error_code)
            sys.exit()
        elif opt in ("-i", "--input"):
            print(arg)
            sam = str(arg)
        elif opt in ("-o", "--output"):
            print(arg)
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

    restart_message = '''
    Gerenuq command executed. 
    To run another gerenuq command, input the following parameters as desired:

    Required inputs:
    -i / --input <input raw samfile>
    -o / --output <output filtered samfile>

    Optional inputs:
    -l / --length <minimum read length for cutoff (default 1000)>
    -m / --matchlength <sequence identity, also known as minimum ratio of matches to read length (default 0.5)>
    -s / --score <minimum score for the whole alignment (default 1)>
    -q / --lengthscore <minimum ratio of length to score, may be considered as the fraction of bases that have a positive score (default 2)>
    -t / --threads <number of processes to run (default 1)>
    '''

    # if the script did not execute with the defined inputs, print error message, otherwise print restart message and await inputs
    return_message = main()
    if return_message: 
        print(restart_message)

    # give the option to restart the script as required
    arguments = list(map(lambda x: ("-" + str(x)), input('Arguments:... ').split("-")))
    print(arguments)
    os.execv(sys.argv[0], arguments)

