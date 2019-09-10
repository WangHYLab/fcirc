#!/usr/bin/env python3
import pysam
import numpy as np
import os
import sys
import logging


class ReconGene:
    def __init__(self,tupleOfReadsA,tupleOfReadsB,fusion_idx_dir):
        self.readsA = tupleOfReadsA
        self.readsB = tupleOfReadsB
        self.GeneAName = tupleOfReadsA[0].reference_name
        self.GeneBName = tupleOfReadsB[0].reference_name
        self.fusionOrient=None
        self.GeneAPointPos = 0 # from 0 ,[start,stop), start(0-based inclusive), stop (0-based exclusive)
        self.GeneBPointPos = 0
        self.Aseq = ""
        self.Bseq = ""
        self.fusion_idx_dir = fusion_idx_dir
        self.reconstructed_seqs = {}
        '''
        self.GeneALength = 2 * max([one.infer_r  ead_length() for one in tupleOfReadsA])
        self.GeneBLength = self.GeneALength
        self.GeneARefSeq = ['' for i in range(self.GeneALength)]
        self.GeneBRefSeq = ['' for i in range(self.GeneBLength)]
        self.GeneAInferSeq = ['' for i in range(self.GeneALength)]
        self.GeneBInferSeq = ['' for i in range(self.GeneBLength)]
        self.fusionPointPos = [0,0] #left gene pos,right gene pos(start from 0)
        '''

    def find_fusionpoint(self,readsA,readsB):
        candidateA = []
        for one in readsA:
            if one.reference_start is None or one.reference_end is None:
                continue
            if one.query_alignment_start > one.infer_read_length()-one.query_alignment_end: # 11S64M0S, 11>0   #one.query_alignment_end-one.infer_read_length()
                if len(one.cigartuples)>=4 and one.cigartuples[1][0] == 0 and one.cigartuples[2][0] == 3 and one.cigartuples[1][1] < 5:
                    candidateA.append(-int(one.reference_start+one.cigartuples[1][1]+one.cigartuples[2][1]))
                else:
                    candidateA.append(-int(one.reference_start))
            else:   # 62M13S
                if len(one.cigartuples)>=4 and one.cigartuples[-2][0] == 0 and one.cigartuples[-3][0] == 3 and one.cigartuples[-2][1] < 5:
                    candidateA.append(one.reference_end-one.cigartuples[-2][1]+one.cigartuples[-3][1])
                else:
                    candidateA.append(one.reference_end)
        candidateGenePointPosA=Merge_breakpointer(dict(zip(*np.unique(candidateA, return_counts=True))))
        PlusOritetion_numA=sum([candidateGenePointPosA[key] for key in candidateGenePointPosA.keys() if key>0])
        MinusOritetion_numA=sum([candidateGenePointPosA[key] for key in candidateGenePointPosA.keys() if key<0])

        candidateB=[]
        for one in readsB:
            if one.reference_start is None or one.reference_end is None:
                continue
            if one.query_alignment_start > one.infer_read_length()-one.query_alignment_end: # 11S64M0S, 11>0   #one.query_alignment_end-one.infer_read_length()
                if len(one.cigartuples)>=4 and one.cigartuples[1][0] == 0 and one.cigartuples[2][0] == 3 and one.cigartuples[1][1] < 5:
                    candidateB.append(-int(one.reference_start+one.cigartuples[1][1]+one.cigartuples[2][1]))
                else:
                    candidateB.append(-int(one.reference_start))
            else:   # 62M13S
                if len(one.cigartuples)>=4 and one.cigartuples[-2][0] == 0 and one.cigartuples[-3][0] == 3 and one.cigartuples[-2][1] < 5:
                    candidateB.append(one.reference_end-one.cigartuples[-2][1]+one.cigartuples[-3][1])
                else:
                    candidateB.append(one.reference_end)
        candidateGenePointPosB=Merge_breakpointer(dict(zip(*np.unique(candidateB, return_counts=True))))
        PlusOritetion_numB=sum([candidateGenePointPosB[key] for key in candidateGenePointPosB.keys() if key>0])
        MinusOritetion_numB=sum([candidateGenePointPosB[key] for key in candidateGenePointPosB.keys() if key<0])

        if len(candidateGenePointPosA)>0 and len(candidateGenePointPosA)>0:
            if (PlusOritetion_numA+MinusOritetion_numB)>(PlusOritetion_numB+MinusOritetion_numA):
                #logging.info("Initial fusion oritentation:%s--->%s"%((readsA[0]).reference_name,(readsB[0].reference_name)))
                orient_Flag=1
            else:
                #logging.info("Initial fusion oritentation:%s--->%s"%(readsB[0].reference_name,readsA[0].reference_name))
                orient_Flag=-1

            if orient_Flag==1:
                candidateGenePointPosA_5=dict([ item for item in candidateGenePointPosA.items() if item[0]>0])
                candidateGenePointPosB_3=dict([ item for item in candidateGenePointPosB.items() if item[0]<0])
                if len(candidateGenePointPosA_5)>0 and len(candidateGenePointPosB_3)>0:
                    GeneAPointPos=(sorted(candidateGenePointPosA_5.items(),key=lambda item:item[1])[-1])[0]
                    GeneBPointPos=-int((sorted(candidateGenePointPosB_3.items(),key=lambda item:item[1])[-1])[0])
                    return orient_Flag,GeneAPointPos,GeneBPointPos
            else:
                candidateGenePointPosA_3=dict([item for item in candidateGenePointPosA.items() if item[0]<0])
                candidateGenePointPosB_5=dict([item for item in candidateGenePointPosB.items() if item[0]>0])
                if len(candidateGenePointPosA_3)>0 and len(candidateGenePointPosB_5)>0:
                    GeneAPointPos=-int(sorted(candidateGenePointPosA_3.items(),key=lambda item:item[1])[-1][0])
                    GeneBPointPos=(sorted(candidateGenePointPosB_5.items(),key=lambda item:item[1])[-1])[0]
                    return orient_Flag,GeneAPointPos,GeneBPointPos
        return None,0,0

    def get_paired_point(self):
        self.fusionOrient,self.GeneAPointPos,self.GeneBPointPos = self.find_fusionpoint(self.readsA,self.readsB)

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

def Merge_breakpointer(dict_var):
    New_dict={}
    Positions=list(dict_var.keys())
    for m in range(len(Positions)):
        for n in range(m+1,len(Positions)):
            if abs(int(Positions[m])-int(Positions[n]))>20:
                continue
            else:
                if dict_var[Positions[m]]>=dict_var[Positions[n]]:
                    dict_var[Positions[m]]=dict_var[Positions[m]]+dict_var[Positions[n]]
                    dict_var[Positions[n]]=0
                else:
                    dict_var[Positions[n]]=dict_var[Positions[n]]+dict_var[Positions[m]]
                    dict_var[Positions[m]]=0
    for i in dict_var.keys():
        if dict_var[i]!=0:
            New_dict[i]=dict_var[i]
    return New_dict


def build_fusion(fasta_path):
    os.system("hisat2-build -q {fasta} {idx}".format(fasta=fasta_path, idx=fasta_path))


def reconstruct(filtered_U_sam, filtered_V_sam,Gene_coordinate,fusion_idx_dir = ''):
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
                    
    with open(Gene_coordinate,'r') as f:
        Gene_coordinate={}
        for line in f:
            Line=((line.strip('\n')).strip()).split('\t')
            Chr,Strand,Start,End=Line[0],Line[1],Line[2],Line[3]
            Gene=Line[-1]
            Gene_coordinate[Gene]=(Chr,Strand,Start,End)                
    

    if os.path.isfile("temp/inferred_fusion.fa"):
        os.system("rm -rf temp/inferred_fusion.fa")
    if os.path.isfile("temp/inferred_fusion_results.tsv"):
        os.system("rm -rf temp/inferred_fusion_results.tsv")

    for key in fusionnameDic.keys():
        if len(fusionnameDic[key]) >= 2:
            recgene = ReconGene(selectReads(fusionnameDic[key],Usam),selectReads(fusionnameDic[key],Vsam),fusion_idx_dir)
            recgene.get_paired_point()
            if recgene.fusionOrient==None or recgene.GeneAPointPos == -1 or recgene.GeneBPointPos == -1: # Judge the breakpoints
                continue
            if Gene_coordinate[recgene.GeneAName][0]!=Gene_coordinate[recgene.GeneBName][0]:   #  A , B in  different Chromosomes
                Translocation_flag=True
            else:
                if int(Gene_coordinate[recgene.GeneAName][2])<int(Gene_coordinate[recgene.GeneBName][2]):
                    Dist=abs(int(Gene_coordinate[recgene.GeneAName][3])-int(Gene_coordinate[recgene.GeneBName][2]))
                else:
                    Dist=abs(int(Gene_coordinate[recgene.GeneAName][2])-int(Gene_coordinate[recgene.GeneBName][3]))
                if Dist>100000:  # A,B in Same Chromosome and separetes distant enough
                    Translocation_flag=True
                else:
                    Translocation_flag=False
            if Translocation_flag==True:
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
    logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(name)s %(levelname)s %(message)s",
                    datefmt = '%Y-%m-%d  %H:%M:%S %a')
    reconstruct(sys.argv[1], sys.argv[2],sys.argv[3],sys.argv[4])







