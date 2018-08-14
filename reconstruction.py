#!/usr/bin/env python3
import pysam
import numpy as np
import os
import sys


class ReconGene:
    def __init__(self,tupleOfReadsA,tupleOfReadsB,fusion_idx_dir):
        self.readsA = tupleOfReadsA
        self.readsB = tupleOfReadsB
        self.GeneAName = tupleOfReadsA[0].reference_name
        self.GeneBName = tupleOfReadsB[0].reference_name
        self.GeneAPointPos = 0 # from 0 ,[start,stop), start(0-based inclusive), stop (0-based exclusive)
        self.GeneBPointPos = 0
        self.Aseq = ""
        self.Bseq = ""
        self.fusion_idx_dir = fusion_idx_dir
        self.reconstructed_seqs = {}
        '''
        self.GeneALength = 2 * max([one.infer_read_length() for one in tupleOfReadsA])
        self.GeneBLength = self.GeneALength
        self.GeneARefSeq = ['' for i in range(self.GeneALength)]
        self.GeneBRefSeq = ['' for i in range(self.GeneBLength)]
        self.GeneAInferSeq = ['' for i in range(self.GeneALength)]
        self.GeneBInferSeq = ['' for i in range(self.GeneBLength)]      
        self.fusionPointPos = [0,0] #left gene pos,right gene pos(start from 0)
        '''

    def find_fusionpoint(self,reads):
        candidate = []
        for one in reads:
            #print('...',one.reference_start , one.reference_end)
            if one.reference_start is None or one.reference_end is None:
                continue
            if one.query_alignment_start > one.infer_read_length() - one.query_alignment_end: # 11S64M0S, 11>0
                if len(one.cigartuples)>=4 and one.cigartuples[1][0] == 0 and one.cigartuples[2][0] == 3 and one.cigartuples[1][1] < 5:
                    candidate.append(one.reference_start+one.cigartuples[1][1]+one.cigartuples[2][1])
                else:
                    candidate.append(one.reference_start)
            else:   # 62M13S
                if len(one.cigartuples)>=4 and one.cigartuples[-2][0] == 0 and one.cigartuples[-3][0] == 3 and one.cigartuples[-2][1] < 5:
                    candidate.append(one.reference_end-one.cigartuples[-2][1]+one.cigartuples[-3][1])
                else:
                    candidate.append(one.reference_end)
            candidateGenePointPos = dict(zip(*np.unique(candidate, return_counts=True)))
        Pos = 0
        PosCount = 0
        for key in candidateGenePointPos.keys():
            if candidateGenePointPos[key] > PosCount:
                Pos = key
                PosCount = candidateGenePointPos[key]
        if PosCount == 1:
            return -1
        return Pos

    def get_paired_point(self):
        self.GeneAPointPos = self.find_fusionpoint(self.readsA)
        self.GeneBPointPos = self.find_fusionpoint(self.readsB)

    def readseq(self):
        Afasta = ''
        Bfasta = ''
        RECORD = False
        with open(os.path.join(self.fusion_idx_dir,"fusiongenes_ref_U.fa")) as Uf:
            for line in Uf:
                if line.startswith('>'):
                    if self.GeneAName in line:
                        RECORD = True
                        continue
                    else:
                        RECORD = False
                if RECORD == True:
                    Afasta += line
            RECORD = False

        with open(os.path.join(self.fusion_idx_dir,"fusiongenes_ref_V.fa")) as Vf:
            for line in Vf:
                if line.startswith('>'):
                    if self.GeneBName in line:
                        RECORD = True
                        continue
                    else:
                        RECORD = False
                if RECORD == True:
                    Bfasta += line
            RECORD = False
        return ''.join(Afasta.split('\n')),''.join(Bfasta.split('\n'))



    def write_inferred_fusion_gene(self):
        self.Aseq, self.Bseq = (seq.upper() for seq in self.readseq())
        self.reconstructed_seqs[self.GeneAName + '-' + self.GeneBName] = connect_sequence(self.Aseq[:self.GeneAPointPos], self.Bseq[self.GeneBPointPos:])
        self.reconstructed_seqs[self.GeneBName + '-' + self.GeneAName] = connect_sequence(self.Bseq[:self.GeneBPointPos], self.Aseq[self.GeneAPointPos:])
        with open("temp/inferred_fusion.fa",'a') as f:
            temp = ''
            for fusiongene in self.reconstructed_seqs.keys():
                if self.reconstructed_seqs[fusiongene] != "":
                    temp += '\n'.join(['>'+fusiongene,self.reconstructed_seqs[fusiongene],''])
            f.write(temp)


def connect_sequence(seq_head,seq_tail):
    if seq_tail == "" or seq_head == "":
        return ""
    seq_head = seq_head.upper()
    seq_tail = seq_tail.lower()
    if seq_head[-1] != seq_tail[0].upper():
        return seq_head + seq_tail
    else:
        for common_len in range(1, int(0.5*min(len(seq_head),len(seq_tail)))):
            if seq_head[-1-common_len] != seq_tail[0+common_len].upper():
                break
        return seq_head + seq_tail[common_len:]


def common_seq(seq_head,seq_tail):
    if seq_tail == "" or seq_head == "":
        return "."
    seq_head = seq_head.upper()
    seq_tail = seq_tail.lower()
    if seq_head[-1] != seq_tail[0].upper():
        return "."
    else:
        for common_len in range(1, int(0.5*min(len(seq_head),len(seq_tail)))):
            if seq_head[-1-common_len] != seq_tail[0+common_len].upper():
                break
        return seq_head[-common_len:]


def selectReads(queryNameList,Reads):
    selectedReads = []
    for one in Reads:
        if one.query_name in queryNameList:
            selectedReads.append(one)
    return tuple(selectedReads)


def build_fusion(fasta_path):
    os.system("hisat2-build -q {fasta} {idx}".format(fasta=fasta_path, idx=fasta_path))


def reconstruct(filtered_U_sam, filtered_V_sam, fusion_idx_dir = ''):
    Usam = tuple(pysam.AlignmentFile(filtered_U_sam, 'r'))
    Vsam = tuple(pysam.AlignmentFile(filtered_V_sam, 'r'))
    refNameDic = {}
    for Uread in Usam:
        if Uread.reference_name not in refNameDic and Uread.reference_name != '*':
            refNameDic[Uread.reference_name] = set([Uread.query_name])
        else:
            refNameDic[Uread.reference_name].add(Uread.query_name)
    for Vread in Vsam:
        if Vread.reference_name not in refNameDic and Vread.reference_name != '*':
            refNameDic[Vread.reference_name] = set([Vread.query_name])
        else:
            refNameDic[Vread.reference_name].add(Vread.query_name)
    #print(refNameDic)
    fusionnameDic = {}
    #print(os.path.join(fusion_idx_dir,"fusion_table.tsv"))
    with open(os.path.join(fusion_idx_dir,"fusion_table.tsv"), 'r') as f:
        for line in f:
            if not line.startswith('#'):
                temp = line.rstrip().split('\t')[0].split('-')
                start = temp[0]
                end = temp[1]
                if len(temp) == 3:
                    # HLA-A
                    start = temp[0] + '-' + temp[1]
                    end = temp[2]
                if start in refNameDic and end in refNameDic:
                    fusionnameDic[line.split('\t')[0]] = refNameDic[start] & refNameDic[end]

    if os.path.isfile("temp/inferred_fusion.fa"):
        os.system("rm -rf temp/inferred_fusion.fa")
    if os.path.isfile("temp/inferred_fusion_results.tsv"):
        os.system("rm -rf temp/inferred_fusion_results.tsv")

    for key in fusionnameDic.keys():
        if len(fusionnameDic[key]) >= 2:

            recgene = ReconGene(selectReads(fusionnameDic[key],Usam),selectReads(fusionnameDic[key],Vsam),fusion_idx_dir)
            recgene.get_paired_point()

            if recgene.GeneAPointPos == -1 or recgene.GeneBPointPos == -1:
                continue

            print(key, len(fusionnameDic[key]))
            print("Pos:",recgene.GeneAPointPos, recgene.GeneBPointPos)



            recgene.write_inferred_fusion_gene()
            with open("temp/inferred_fusion_results.tsv",'a') as f:
                f.write("{fusion5}\t{fusion3}\t{f5_pos}\t{f3_pos}\t{f5_seq}\t{f3_seq}\t{f53_com}\t{seqlen}\n".format(
                    fusion5=recgene.GeneAName,
                    fusion3=recgene.GeneBName,
                    f5_pos=recgene.GeneAPointPos,
                    f3_pos=recgene.GeneBPointPos,
                    f5_seq=recgene.Aseq[recgene.GeneAPointPos - 10:recgene.GeneAPointPos],
                    f3_seq=recgene.Bseq[recgene.GeneBPointPos:recgene.GeneBPointPos + 10],
                    f53_com=common_seq(recgene.Aseq[recgene.GeneAPointPos-10:recgene.GeneAPointPos],
                                       recgene.Bseq[recgene.GeneBPointPos:recgene.GeneBPointPos + 10]),
                    seqlen=str(len(recgene.reconstructed_seqs['-'.join((recgene.GeneAName, recgene.GeneBName))]))))
                f.write("{fusion5}\t{fusion3}\t{f5_pos}\t{f3_pos}\t{f5_seq}\t{f3_seq}\t{f53_com}\t{seqlen}\n".format(
                    fusion5=recgene.GeneBName,
                    fusion3=recgene.GeneAName,
                    f5_pos=recgene.GeneBPointPos,
                    f3_pos=recgene.GeneAPointPos,
                    f5_seq=recgene.Bseq[recgene.GeneBPointPos - 10:recgene.GeneBPointPos],
                    f3_seq=recgene.Aseq[recgene.GeneAPointPos:recgene.GeneAPointPos + 10],
                    f53_com=common_seq(recgene.Bseq[recgene.GeneBPointPos - 10:recgene.GeneBPointPos],
                                       recgene.Aseq[recgene.GeneAPointPos:recgene.GeneAPointPos + 10]),
                    seqlen=str(len(recgene.reconstructed_seqs['-'.join((recgene.GeneBName, recgene.GeneAName))]))))
    if fusionnameDic != {}:
        build_fusion("temp/inferred_fusion.fa")
        return 0
    return None



if __name__ == "__main__":
        reconstruct(sys.argv[1], sys.argv[2], sys.argv[3])











