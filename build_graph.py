#!/usr/bin/env python3

import re
import os,getopt
import time
import sys
from pysam import FastaFile



def getgraph(location):
    '''
    Load the fusion gene table from the local file
    :param location: a string of fusion gene table location
    :return: a dict of fusion genes
    '''
    if not (os.path.exists(os.curdir+'/'+location) or os.path.exists(location)):
        print('Error: There is no %s' %location)
        exit(1)
    with open(location, 'r') as f:
        edge = {}
        nodes = []
        for line in f:
            if not line.startswith('#'):
                try:
                    temp = line.strip('\n').split('\t')[0].split('--')
                    start = temp[0]
                    end = temp[1]
                except IndexError as e:
                    print("%s: gene pairs placed in first columns should be joined with '--' not '-' "%e)
                    exit(1)

                if start in edge:
                    edge[start].append(end)
                else:
                    edge[start] = [end]

                if end in edge:
                    edge[end].append(start)
                else:
                    edge[end] = [start]
    return edge

def find_bool(graph):
    thisside = set()
    thatside = set()
    nodes = list(graph.keys())
    thisside.add(nodes[0])
    while(nodes!= []):
        edges_count = sum([len(graph[node]) for node in nodes])
        if edges_count == 0:
            break
        for start in graph.keys():
            if graph[start] == []:
                continue
            if start in thisside:
                for side in graph[start]:
                    thatside.add(side)
                graph[start] = []
                nodes.remove(start)
            elif start in thatside:
                for side in graph[start]:
                    thisside.add(side)
                graph[start] = []
                nodes.remove(start)
            if thisside & thatside != set():
                return None
        if sum([len(graph[node]) for node in nodes]) == edges_count:
            thisside.add(nodes[0])
    return list(thisside), list(thatside)






def Get_geneCoordinate(gtf):
    if not (os.path.exists(gtf) or (os.path.exists(os.curdir+'/'+gtf))):
        print("there is no gtf file,please check the path")
        exit(1)
    Out=open('fusion_genes_coordinate.txt','w')
    GeneCoordinate_dict={}
    with open(gtf,'r') as f:
        for line in f:
            if not line.startswith('#'):
                Line=line.strip('\n').split('\t')
                chr,Type,start,end,strand,gene=Line[0],Line[2],Line[3],Line[4],Line[6],Line[-1].split(';')[2].replace('gene_name "','').replace('"','').strip()
                if Type=='gene':
                    if strand=='+':
                        strand='1'
                    else:
                        strand='-1'
                    gene_coordinates='\t'.join([chr,strand,start,end,gene,'\n'])
                    if gene not in GeneCoordinate_dict:
                        GeneCoordinate_dict[gene]=[chr,strand,start,end,gene]
                    Out.write(gene_coordinates)
    Out.close()
    return GeneCoordinate_dict



def Get_fusionGene_seq(GenesU,GenesV,geneCoordniates_dict,reference_fa=''):
    if not (os.path.exists(reference_fa) or os.path.exists(os.curdir+'/'+reference_fa)):
        print('Error: There is no reference fasta files')
        exit(1)
    genome_name=os.path.basename(reference_fa)
    if not os.path.exists(genome_name):
        genome=os.system('ln -s %s %s'%(reference_fa,genome_name))
    genome_index=os.system('samtools faidx %s'%genome_name)
    genome=FastaFile(genome_name)


    fusiongenes_ref_U=open('fusion_total_index/fusiongenes_ref_U.fa','w')
    for gene in GenesU:
        chr, strand, start, end=None,None,None,None
        try:
            chr, strand, start, end, gene = geneCoordniates_dict[gene]
        except KeyError as e:
            print('%s: Input gene name wasnot found in Gtf,Check gene names'%e)
            exit(1)
        if strand!=None:
            if strand == '1':
                seq = genome.fetch(reference=chr, start=int(start), end=int(end))
            else:
                seq_plus = genome.fetch(reference=chr, start=int(start), end=int(end))
                trantab = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
                seq = seq_plus.translate(trantab)
                seq = seq[::-1]
            fusiongenes_ref_U.write('>%s \n' % gene)
            for line in re.findall(r'.{60}', seq):
                fusiongenes_ref_U.write('%s\n' % line)
    fusiongenes_ref_U.close()


    fusiongenes_ref_V=open('fusion_total_index/fusiongenes_ref_V.fa','w')
    for gene in GenesV:
        chr, strand, start, end = None, None, None, None
        try:
            chr, strand, start, end, gene = geneCoordniates_dict[gene]
        except KeyError as e:
            print('%s: Input gene name wasnot found in Gtf,Check gene names'%e)
            exit(1)
        if strand != None:
            if strand == '1':
                seq = genome.fetch(reference=chr, start=int(start), end=int(end))
            else:
                seq_plus = genome.fetch(reference=chr, start=int(start), end=int(end))
                trantab = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
                seq = seq_plus.translate(trantab)
                seq = seq[::-1]
            fusiongenes_ref_V.write('>%s \n' % gene)
            for line in re.findall(r'.{60}', seq):
                fusiongenes_ref_V.write('%s\n' % line)
    fusiongenes_ref_V.close()
    return 0



def main(fusion_tab,gtf,ref_genome):
    if not os.path.exists('fusion_total_index'):
        os.mkdir('fusion_total_index')
    graph = getgraph(fusion_tab)
    bipgraph=find_bool(graph)
    Flag=False
    if bipgraph is not None:
        if (len(bipgraph[0])+len(bipgraph[1])==len(graph.keys())):
            Flag=True
            genesU, genesV = bipgraph[0], bipgraph[1]
            geneCoordinate = Get_geneCoordinate(gtf)
            fusion_seqs = Get_fusionGene_seq(genesU, genesV, geneCoordinate, ref_genome)
            if fusion_seqs == 0:
                print("Sucessfully get the fusion genes sequence ")
            else:
                print("Failed to get the fusion genes sequence")
                exit(1)
    if not Flag:
        print("Failed to build the bipartite graph pairs")
    else:
        print("Sucessfully building the bipartite graph pairs")
        exit(1)






if __name__ == '__main__':
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
    usage_string = '''
       Usage:  python Preprocess.py -genome <genome-reference-path> -gtf <gtf-reference-path> -tab <fusion-table>
       Arguments:
         Required:
            -genome <genome-reference-path>
              genome reference fasta file, default: Homo_sapiens.GRCh38.dna.primary_assembly.fa
            -gtf <gtf-reference-path>
              gene transfer format file   default: Homo_sapiens.GRCh38.95.gtf
            -tab <fusion-table>
              Inputed fusion gene pairs table,file format should be tsv, fusion pairs Gene pairs should be connected with -- not -, default: fusion_total_index/fusion_table.tsv
         Other:
            -h --help
              help information
      '''
    if len(sys.argv[1:])==0:
        print(usage_string)
        sys.exit(1)
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h",["help","genome=","gtf=","tab="])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(1)

    genome='Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    gtf='Homo_sapiens.GRCh38.95.gtf'
    fusion_tab='fusion_total_index/fusion_table.tsv'

    for opt, arg in opts:
        if opt in ("-h", '--help'):
            print(usage_string)
            sys.exit()
        elif (opt =="--genome"):
            genome=arg
        elif (opt =="--gtf"):
            gtf=arg
        elif (opt =="--tab"):
            fusion_tab=arg
        else:
            assert False, "unhandled option"
    main(fusion_tab,gtf,genome)







