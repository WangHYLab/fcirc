#!/usr/bin/env python3
import sys
import os
import pysam
from get_graph import getgraph


def bool_paired(reada, readb):
    if reada.query_name != readb.query_name:    # exclude different reads
        return 0
    if not (reada.cigartuples and readb.cigartuples):  # unmapped end in paired-end reads
        return 0
    if len(reada.cigartuples) == 1 or len(readb.cigartuples) == 1:  # exclude only include M
        return 0
    if reada.flag & 16 == readb.flag & 16:  # mapped to the same strand
        if reada.cigartuples[0][0] == 0 and readb.cigartuples[-1][0] == 0 or \
                (reada.cigartuples[-1][0] == 0 and readb.cigartuples[0][0] == 0):   # 32S43M & 32M43S
            return 1
        if len(reada.cigartuples) == 3 and reada.cigartuples[-1][0] != 0 and reada.cigartuples[-1][1] <= 3 and readb.cigartuples[0][0] == 0 or \
                (reada.cigartuples[0][0] == 0 and len(readb.cigartuples) == 3 and readb.cigartuples[-1][0] != 0 and readb.cigartuples[-1][1] <= 3) or \
                (len(reada.cigartuples) == 3 and reada.cigartuples[0][0] != 0 and reada.cigartuples[0][1] <= 3 and readb.cigartuples[-1][0] == 0) or \
                (reada.cigartuples[-1][0] == 0 and len(readb.cigartuples) == 3 and readb.cigartuples[0][0] != 0 and readb.cigartuples[0][1] <= 3):  # 41S30M2S & 41MXXX | 2S30M41S & XXX41M
            return 1
    else:
        if reada.cigartuples[0][0] == 0 and readb.cigartuples[0][0] == 0 or \
                (reada.cigartuples[-1][0] == 0 and readb.cigartuples[-1][0] == 0):  # 32S43M & 43S32M
            return 1
        if len(reada.cigartuples) == 3 and reada.cigartuples[-1][0] != 0 and reada.cigartuples[-1][1] <= 3 and readb.cigartuples[-1][0] == 0 or \
                (reada.cigartuples[-1][0] == 0 and len(readb.cigartuples) == 3 and readb.cigartuples[-1][0] != 0 and readb.cigartuples[-1][1] <= 3) or \
                (len(reada.cigartuples) == 3 and reada.cigartuples[0][0] != 0 and reada.cigartuples[0][1] <= 3 and readb.cigartuples[0][0] == 0) or \
                (reada.cigartuples[0][0] == 0 and len(readb.cigartuples) == 3 and readb.cigartuples[0][0] != 0 and readb.cigartuples[0][1] <= 3):  # 41S30M2S & XXX41M | 2S30M41S & 41MXXX
            return 1


def bool_bridge(pairreada,pairreadb):
    if pairreada[0].cigartuples and not pairreada[1].cigartuples and not pairreadb[0].cigartuples and pairreadb[1].cigartuples:
        return 1
    elif not pairreada[0].cigartuples and pairreada[1].cigartuples and pairreadb[0].cigartuples and not pairreadb[1].cigartuples:
        return 1
    else:
        return 0



def getpairedreads(samplenameU, samplenameV, pattern='single-end', fusion_idx_dir=''):
    samfile_U = pysam.AlignmentFile(samplenameU, 'r')

    filtered_from_U = {}
    filtered_from_V = {}

    fusion_dic = getgraph(os.path.join(fusion_idx_dir, "fusion_table.tsv"))
    samfile_V = pysam.AlignmentFile(samplenameV, 'r')

    # single-end
    if pattern == 'single-end':
        reads_V = {}
        for one in samfile_V:
            if one.reference_name not in reads_V:
                reads_V[one.reference_name] = {}

            # store the primary alignment
            if one.query_name not in reads_V[one.reference_name] or not one.flag & 256:
                reads_V[one.reference_name][one.query_name] = one


        for read in samfile_U:
            partners = fusion_dic[read.reference_name]
            for pa in partners:
                if pa not in reads_V:
                    continue
                if read.query_name in reads_V[pa]:
                    anotherread = reads_V[pa][read.query_name]
                    if bool_paired(read, anotherread):
                        filtered_from_U[read.query_name] = read
                        filtered_from_V[anotherread.query_name] = anotherread
                        break


        print("Find " + str(len(filtered_from_U)) + " Reads in U! " + str(len(filtered_from_V)) + " Reads in V!")

        # write
        sam_filter_U_path = samplenameU.replace(".sam","_filtered.sam")
        Ureads = pysam.AlignmentFile(sam_filter_U_path, 'w', template=samfile_U)
        for oneU in filtered_from_U.keys():
            Ureads.write(filtered_from_U[oneU])
        Ureads.close()

        sam_filter_V_path = samplenameV.replace(".sam", "_filtered.sam")
        Vreads = pysam.AlignmentFile(sam_filter_V_path, 'w', template=samfile_V)
        for oneV in filtered_from_V.keys():
            Vreads.write(filtered_from_V[oneV])
        Vreads.close()



    # paired-end
    if pattern == 'paired-end':
        reads_V = {}
        temp_paired = []

        for one in samfile_V: # reads sorted by query_name works
            temp_paired.append(one)
            if len(temp_paired) == 2:
                if one.query_name == temp_paired[0].query_name:
                    if temp_paired[0].flag & 128 and temp_paired[1].flag & 64:
                        temp_paired = [temp_paired[1],temp_paired[0]]
                    if temp_paired[0].flag & 64 and temp_paired[1].flag & 128:
                        # add into reads_V
                        if temp_paired[0].reference_name != "*":
                            if temp_paired[0].reference_name not in reads_V:
                                reads_V[temp_paired[0].reference_name] = {}
                            reads_V[temp_paired[0].reference_name][temp_paired[0].query_name] = tuple(temp_paired)

                        if temp_paired[1].reference_name != "*":
                            if temp_paired[1].reference_name not in reads_V:
                                reads_V[temp_paired[1].reference_name] = {}
                            reads_V[temp_paired[1].reference_name][temp_paired[1].query_name] = tuple(temp_paired)
                        temp_paired = []
                    else:
                        temp_paired = [temp_paired[0] if temp_paired[1].flag & 256 else temp_paired[1]]
                else:
                    temp_paired = [one]



        temp_paired = []
        for read in samfile_U:
            temp_paired.append(read)
            if len(temp_paired) == 2:
                if read.query_name == temp_paired[0].query_name:
                    if temp_paired[0].flag & 128 and temp_paired[1].flag & 64:
                        temp_paired = [temp_paired[1],temp_paired[0]]
                    if temp_paired[0].flag & 64 and temp_paired[1].flag & 128:
                        if temp_paired[0].reference_name != "*" and temp_paired[1].reference_name != "*" and \
                                temp_paired[0].reference_name != temp_paired[1].reference_name:
                            temp_paired = []
                            continue
                        if temp_paired[0].reference_name != "*":
                            partners1 = fusion_dic[temp_paired[0].reference_name]
                            for pa1 in partners1:
                                if pa1 not in reads_V:
                                    continue
                                if read.query_name in reads_V[pa1]:
                                    anotherpairedread = reads_V[pa1][read.query_name]
                                    if bool_paired(anotherpairedread[0], temp_paired[0]) or bool_bridge(temp_paired, anotherpairedread):
                                        filtered_from_U[read.query_name] = tuple(temp_paired)
                                        filtered_from_V[read.query_name] = anotherpairedread
                                        break

                        if temp_paired[1].reference_name != "*":
                            partners2 = fusion_dic[temp_paired[1].reference_name]
                            for pa2 in partners2:
                                if pa2 not in reads_V:
                                    continue
                                if read.query_name in reads_V[pa2]:
                                    anotherpairedread = reads_V[pa2][read.query_name]
                                    if bool_paired(anotherpairedread[1],temp_paired[1]) or bool_bridge(temp_paired, anotherpairedread):
                                        filtered_from_U[read.query_name] = tuple(temp_paired)
                                        filtered_from_V[read.query_name] = anotherpairedread
                                        break
                        temp_paired = []
                    else:
                        temp_paired = [temp_paired[0] if temp_paired[1].flag & 256 else temp_paired[1]]
                else:
                    temp_paired = [read]

        print('Find ' + str(len(filtered_from_U)) + ' Reads in U! ' + str(len(filtered_from_V)) + ' Reads in V!')

        # write
        sam_filter_U_path = samplenameU.replace(".sam", "_filtered.sam")
        Ureads = pysam.AlignmentFile(sam_filter_U_path, 'w', template=samfile_U)
        for paironeU in filtered_from_U.keys():
            for oneU in filtered_from_U[paironeU]:
                Ureads.write(oneU)
        Ureads.close()

        sam_filter_V_path = samplenameV.replace(".sam", "_filtered.sam")
        Vreads = pysam.AlignmentFile(sam_filter_V_path, 'w', template=samfile_V)
        for paironeV in filtered_from_V.keys():
            for oneV in filtered_from_V[paironeV]:
                Vreads.write(oneV)
        Vreads.close()



    samfile_U.close()
    samfile_V.close()

    return sam_filter_U_path, sam_filter_V_path




if __name__ == '__main__':
    getpairedreads(sys.argv[1], sys.argv[2], sys.argv[3])
