# Alec Bahcheli, Daniel Giguire from Gloor Lab, Western University, Canada

import re, time, pysam
import pandas as pd

def filter_bamfile(read, min_score, min_len_to_score, min_length, min_match_to_length):
    # split the read into the expected fields
    bitflag = read[1]
    cigar_string = read[2]
    alignment_score = read[3]

    # # make sure it's a primary alignment
    # if int(bitflag) == 16 or int(bitflag) == 0:
    # read mapping score
    score = int(alignment_score[5:])
    # is the score high enough
    if score >= min_score:
        # read length 
        length = 0
        # number of matches
        num_of_matches = 0
        # length mapped on genome
        genome_mapped = 0

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
                genome_mapped += int(element.strip("M"))
            
            # if there is a match, length or mismatch 
            elif re.search("D", element):
                genome_mapped += int(element.strip("D"))


        # are we dealing with penalizing the ends or not
        if pen_ends:
            # evaluate whether the read is close to the end of the genome
            contig_mapped = read[2]

            # end position of the genome where the read maps to the genome
            mapping_end = int(read[3]) + genome_mapped

            # corrected length is length minus the last clipping
            corrected_length = length

            # double using the parameter of minimum length of read to evaluate whether the read qualifies for end of genome mapping consideration
            if (sequence_length_dict[contig_mapped] - mapping_end) <= min_length:
                # corrected length is length minus the last clipping
                if re.search("[SH]", cigar[-1]):
                    corrected_length = corrected_length - int(cigar[-1].strip("[SH]"))

            # beginning of genome also has different consideration
            elif int(read[3]) <= min_length:
                if re.search("[SH]", cigar[0]):
                    corrected_length = corrected_length - int(cigar[0].strip("[SH]"))

            # evaluate if the matches and minimum length to score with the corrected length
            if (length > min_length) and ((num_of_matches / corrected_length) > min_match_to_length) and ((int(corrected_length) / int(score)) < min_len_to_score):
                return ("\t".join(read))
        
        # if it meets the filter cutoffs, return the whole read
        elif (length > min_length) and ((num_of_matches / length) > min_match_to_length) and ((int(length) / int(score)) < min_len_to_score):
            
            return ("\t".join(read))


def filter_samfile(read, min_score, min_len_to_score, min_length, min_match_to_length):
    # split the read into the expected fields
    read = read.split("\t")
    bitflag = read[1]
    cigar_string = read[5]
    alignment_score = read[13]
    # if int(bitflag) == 16 or int(bitflag) == 0:
    # read mapping score
    score = int(alignment_score[5:])
    # is the score high enough
    if score >= min_score:
        # read length 
        length = 0
        # number of matches
        num_of_matches = 0
        # length mapped on genome
        genome_mapped = 0

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
                genome_mapped += int(element.strip("M"))
            
            # if there is a match, length or mismatch 
            elif re.search("D", element):
                genome_mapped += int(element.strip("D"))


        # are we dealing with penalizing the ends or not
        if pen_ends:
            # evaluate whether the read is close to the end of the genome
            contig_mapped = read[2]

            # end position of the genome where the read maps to the genome
            mapping_end = int(read[3]) + genome_mapped

            # corrected length is length minus the last clipping
            corrected_length = length

            # double using the parameter of minimum length of read to evaluate whether the read qualifies for end of genome mapping consideration
            if (sequence_length_dict[contig_mapped] - mapping_end) <= min_length:
                # corrected length is length minus the last clipping
                if re.search("[SH]", cigar[-1]):
                    corrected_length = corrected_length - int(cigar[-1].strip("[SH]"))

            # beginning of genome also has different consideration
            elif int(read[3]) <= min_length:
                if re.search("[SH]", cigar[0]):
                    corrected_length = corrected_length - int(cigar[0].strip("[SH]"))

            # evaluate if the matches and minimum length to score with the corrected length
            if (length > min_length) and ((num_of_matches / corrected_length) > min_match_to_length) and ((int(corrected_length) / int(score)) < min_len_to_score):
                return ("\t".join(read))
        
        # if it meets the filter cutoffs, return the whole read
        elif (length > min_length) and ((num_of_matches / length) > min_match_to_length) and ((int(length) / int(score)) < min_len_to_score):
            
            return ("\t".join(read))

def gerenuq_filter_file(input_file, output_file, min_score = 1, min_len_to_score = 2, min_length = 1000, min_match_to_length = 0.5, penalize_ends=True):
    '''
    gerenuq_filter_file(input_file, output_file, min_score = 1, min_len_to_score = 2, min_length = 1000, min_match_to_length = 0.5)

    Filters minimap2-mapped reads by mapping score, length, match-to-length and length-to-score ratios. Paf format files only filter by query cutoff.

    Requires input_file in bam, sam or paf format and output_file (output in the same format as input).

    version 0.2.6
    '''
    t1 = time.time()

    # globalize penalizing variable
    global pen_ends
    pen_ends = penalize_ends

    global sequence_length_dict
    sequence_length_dict = {}

    # default format is sam
    file_format = 'sam'
    if re.search(r'\.[pP][aA][fF]', input):
        file_format = 'paf'
    elif re.search(r'\.[bB][aA][mM]', input):
        file_format = 'bam'
    # run process
    if file_format == 'paf': 
        df=pd.read_csv(input_file, sep="\t", header=None)
        
        print("File read into memory for parallelization")
        print(round(time.time() - t1, 2))
        # filter reads by query cutoff 
        filtered_reads = df[ (df[3]-df[2] ) / df[1] > min_match_to_length]
        
        print("Done filtering reads")
        print(round(time.time() - t1, 2))
        # output file
        filtered_reads.to_csv(output_file, sep = "\t", index=False, header=False)

        print("Finished writing filtered paf file")
        print(round(time.time() - t1, 2))

    elif file_format == 'bam':
        # make a list of just reads
        samfile = []
        infile = pysam.AlignmentFile(input_file, 'rb')

        # iterate over header
        for read in str(infile.header).split('\n'):
            # append sequence length to dictionary
            if read.startswith("@SQ"):
                contig = read.split("\t")
                sequence_length_dict[contig[1].split(":")[1]] = int(contig[2].split(":")[1])

        for read in infile:
            # add in order: read name, alignment type (index 1), CIGAR string, and alignment score (final index after 'AS:i:') in that order
            read = str(read).split("\t")
            samfile.append(list((read[0], read[1], read[5], ((read[-1].split("), ("))[2]).split(", ")[1])))
    
        print("File read into memory for parallelization")
        print(round(time.time() - t1, 2))
        # filter reads
        reads_of_interest = dict(list(filter(None, map(lambda x: filter_bamfile(x, min_score=min_score, min_len_to_score=min_len_to_score, min_length=min_length, min_match_to_length=min_match_to_length), samfile))))

        print("Done filtering reads")
        print(round(time.time() - t1, 2))
        print(len(reads_of_interest.keys()))
        # compare to bam file by read name (dictionary key) and CIGAR score (dictionary value)
        infile = pysam.AlignmentFile(input_file, 'rb')
        bam_output = pysam.AlignmentFile(output_file, "wb", template=infile)
        for read in infile:
            if read.query_name in reads_of_interest:
                if read.cigarstring == reads_of_interest.get(read.query_name):
                    bam_output.write(read)
        
        print("Finished writing filtered bam file")
        print(round(time.time() - t1, 2))

    else:
        # open the results file
        results = open(output_file, "w")
        # make a list of just reads
        samfile = []
        with open(input_file) as input_raw:
            # get the headers
            for read in input_raw:
                if read.startswith("@"):
                    results.write(read)

                    # append sequence length to dictionary
                    if read.startswith("@SQ"):
                        contig = read.split("\t")
                        sequence_length_dict[contig[1].split(":")[1]] = int(contig[2].split(":")[1])

                else:
                    samfile.append(read) 

        print("File read into memory for parallelization")
        print(round(time.time() - t1, 2))
        # evaluate reads
        reads_of_interest = list(filter(None, map(lambda x: filter_bamfile(x, min_score=min_score, min_len_to_score=min_len_to_score, min_length=min_length, min_match_to_length=min_match_to_length), samfile)))
        
        print("Done filtering reads")
        print(round(time.time() - t1, 2))
        for read in reads_of_interest:
            # write the good reads to a file
            results.write(read)
        results.close()

        print("Finished writing filtered sam file")
        print(round(time.time() - t1, 2))
        print("done")

def gerenuq_filter_read_list(read_list, format='sam', min_score = 1, min_len_to_score = 2, min_length = 1000, min_match_to_length = 0.5):
    '''
    gerenuq_filter_read_list(read_list, format='sam', min_score = 1, min_len_to_score = 2, min_length = 1000, min_match_to_length = 0.5)

    Filters minimap2-mapped reads by mapping score, length, match-to-length and length-to-score ratios. Paf format files only filter by query cutoff.

    Requires read_list as list of mapped read lines from sam or paf file (in tsv format). Returns a list of reads in sam or paf (tsv) format that passed filtering parameters. Headers will be ignored and not returned.

    version 0.2.6
    '''
    # run process
    if format == 'paf': 
        read_list = list(filter(lambda x: not x.startswith("@"), read_list))
        # filter out headers
        df=pd.DataFrame(read_list)
        # filter reads by query cutoff 
        filtered_reads = df[ (df[3]-df[2] ) / df[1] > min_match_to_length]
        # return results
        return_list = []
        total_list = 0
        for i, row in filtered_reads.iterrows():
            total_list += 1
            return_list.append("\t".join(row.tolist()))
        return(return_list)
    else:
        # make a list of just reads, no headers
        read_list = list(filter(lambda x: not x.startswith("@"), read_list))
        # evaluate reads
        reads_of_interest = list(filter(None, map(lambda x: filter_bamfile(x, min_score=min_score, min_len_to_score=min_len_to_score, min_length=min_length, min_match_to_length=min_match_to_length), read_list)))
        # return results
        return(reads_of_interest)



