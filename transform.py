#!/usr/bin/env python3
import pysam
import sys

def reverse_read(read, reversepostion):
    if reversepostion == 0:
        softseq = read.query_sequence[:read.cigartuples[0][1]]
        softquality = read.query_qualities[:read.cigartuples[0][1]]
        matchseq = read.query_sequence[read.cigartuples[0][1]:]
        matchquality = read.query_qualities[read.cigartuples[0][1]:]
        read.query_sequence = matchseq + softseq
        read.query_qualities = matchquality + softquality
    if reversepostion == -1:
        softseq = read.query_sequence[-read.cigartuples[0][1]:]
        softquality = read.query_qualities[-read.cigartuples[0][1]:]
        matchseq = read.query_sequence[:-read.cigartuples[0][1]]
        matchquality = read.query_qualities[:-read.cigartuples[0][1]]
        read.query_sequence = softseq + matchseq
        read.query_qualities = softquality + matchquality
    return read


def checkcigar(read):
    if not read.cigartuples:
        return None
    elif read.cigartuples[0][0] == 4 and read.cigartuples[-1][0] == 4:
        if max(read.cigartuples[0][1], read.cigartuples[-1][1]) < 3:
            return None
        if read.cigartuples[0][1] > read.cigartuples[-1][1]:  # 44S22M6S
            return 0
        else:  # 6S22M44S
            return -1
    elif read.cigartuples[0][0] == 4 and read.cigartuples[0][1] >= 3:  # 44S31M
        return 0
    elif read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] >= 3:  # 31M44S
        return -1



def transform_sam(samname):
    '''
    Transfrom a sam file of candicate circRNA reads
    :param samname: a string of sam file path
    :return: tranformed sam file path
    '''
    samfile = pysam.AlignmentFile(samname, 'r')
    transformed_set = set()

    for read in samfile:
        softpos = checkcigar(read)
        if softpos != None:
            transformed_set.add(reverse_read(read, softpos))

    tranformed_sam_path = samname.replace(".sam", "_tranformed.sam")
    trans_file = pysam.AlignmentFile(tranformed_sam_path, 'w', template=samfile)
    for one in transformed_set:
        trans_file.write(one)

    trans_file.close()
    samfile.close()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('''
Usage:Transform the candidate circRNA reads' sequence(as well as quality) which may be mapped to fusion gene, drop linear mapping reads, A tranformed fastq file from the sam file in the same dictionary
Input:A string of a sam file
Output:Fastq path 
        ''')
        sys.exit()
    path = transform_sam(sys.argv[1])
