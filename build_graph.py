#!/usr/bin/env python3

import re
import os
import getopt
import sys
from pysam import FastaFile


def getgraph(location):
    """
    Load the fusion gene table from the local file
    :param location: a string of fusion gene table location
    :return: a dict of fusion genes
    """
    if not (os.path.exists(os.curdir + '/' + location) or os.path.exists(location)):
        print('Error: There is no %s' % location)
        exit(1)
    with open(location, 'r') as f:
        edge = {}
        for line in f:
            if not line.startswith('#'):
                try:
                    temp = line.strip('\n').split('\t')[0].split('--')
                    start = temp[0]
                    end = temp[1]
                except IndexError as e:
                    print("%s: gene pairs placed in first columns should be joined with '--' not '-' " % e)
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
    while nodes:
        edges_count = sum([len(graph[node]) for node in nodes])
        if edges_count == 0:
            break
        for start in graph.keys():
            if not graph[start]:
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


def Get_geneCoordinate(gtf, od):
    synonyms_dict = {"SEPTIN5": "SEPT5",
                     "SIKE1": "SIKE",
                     "AC011462.1": "TGFB1",
                     "SEPTIN9": "SEPT9",
                     "SEPTIN6": "SEPT6",
                     "H4C9": "HIST1H4I",
                     "GET1": "WRB",
                     "H2BC4": "HIST1H2BC",
                     "SEPTIN11": "SEPT11",
                     "ZFTA": "C11orf95",
                     "H2AC17": "HIST1H2AM",
                     "HLA-DMA": "HLA_DMA",
                     "HLA-DMB": "HLA_DMB",
                     "RNF213": "AC124319.1",
                     "PEDS1": "TMEM189",
                     "DYNLT2B": "TCTEX1D2",
                     "TGFB1": "AC011462.1",
                     "SEPTIN2": "SEPT2",
                     "SEPTIN8": "SEPT8",
                     "ERVK3-1": "ERVK3_1",
                     "HI-6": "HIST1H1T",
                     "CEP43": "FGFR1OP",
                     "H3C12": "HIST1H3J",
                     "CFAP251": "WDR66",
                     "MAGEA5P": "MAGEA5",
                     "BMERB1": "C16orf45",
                     "EZHIP": "CXorf67",
                     "H1-6": "HIST1H1T"}
    if not (os.path.exists(gtf) or (os.path.exists(os.curdir + '/' + gtf))):
        print("there is no gtf file, please check the path")
        exit(1)
    Out = open('{od}/fusion_genes_coordinate.txt'.format(od=od), 'w')
    GeneCoordinate_dict = {}
    with open(gtf, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                Line = line.strip('\n').split('\t')
                chrn, Type, start, end, strand, gene = Line[0], Line[2], Line[3], Line[4], Line[6], Line[-1].split(';')[
                    2].replace('gene_name "', '').replace('"', '').strip()
                if Type == 'gene':
                    if strand == '+':
                        strand = '1'
                    else:
                        strand = '-1'
                    if gene in synonyms_dict:
                        gene = synonyms_dict[gene]
                    gene_coordinates = '\t'.join([chrn, strand, start, end, gene, '\n'])
                    if gene not in GeneCoordinate_dict:
                        GeneCoordinate_dict[gene] = [chrn, strand, start, end, gene]
                    Out.write(gene_coordinates)
    Out.close()
    return GeneCoordinate_dict


def Get_fusionGene_seq(GenesU, GenesV, geneCoordniates_dict, reference_fa, od):
    if not (os.path.exists(reference_fa) or os.path.exists(os.curdir + '/' + reference_fa)):
        print('Error: There is no reference fasta files')
        exit(1)
    genome_name = os.path.basename(reference_fa)
    if not os.path.exists(genome_name):
        os.system('ln -s %s %s' % (os.path.abspath(os.curdir + '/' + reference_fa), od + '/' + genome_name))
    # fasta index
    os.system('samtools faidx %s' % od + '/' + genome_name)
    genome = FastaFile(od + '/' + genome_name)

    fusiongenes_ref_U = open('{od}/fusiongenes_ref_U.fa'.format(od=od), 'w')
    for gene in GenesU:
        chrn, strand, start, end = None, None, None, None
        try:
            chrn, strand, start, end, gene = geneCoordniates_dict[gene]
        except KeyError as e:
            print('%s: Input gene name wasnot found in Gtf, Check gene names' % e)
            exit(1)
        if strand:
            if strand == '1':
                seq = genome.fetch(reference=chrn, start=int(start), end=int(end))
            else:
                seq_plus = genome.fetch(reference=chrn, start=int(start), end=int(end))
                trantab = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
                seq = seq_plus.translate(trantab)
                seq = seq[::-1]
            fusiongenes_ref_U.write('>%s\n' % gene)
            for line in re.findall(r'.{60}', seq):
                fusiongenes_ref_U.write('%s\n' % line)
    fusiongenes_ref_U.close()

    fusiongenes_ref_V = open('{od}/fusiongenes_ref_V.fa'.format(od=od), 'w')
    for gene in GenesV:
        chrn, strand, start, end = None, None, None, None
        try:
            chrn, strand, start, end, gene = geneCoordniates_dict[gene]
        except KeyError as e:
            print('%s: Input gene name wasnot found in Gtf, Check gene names' % e)
            exit(1)
        if strand:
            if strand == '1':
                seq = genome.fetch(reference=chrn, start=int(start), end=int(end))
            else:
                seq_plus = genome.fetch(reference=chrn, start=int(start), end=int(end))
                trantab = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
                seq = seq_plus.translate(trantab)
                seq = seq[::-1]
            fusiongenes_ref_V.write('>%s\n' % gene)
            for line in re.findall(r'.{60}', seq):
                fusiongenes_ref_V.write('%s\n' % line)
    fusiongenes_ref_V.close()
    return 0


def main(fusion_tab, gtf, ref_genome, od):
    if not os.path.exists(od):
        os.mkdir(od)
    graph = getgraph(fusion_tab)
    bipgraph = find_bool(graph)
    Flag = False
    if bipgraph is not None:
        if len(bipgraph[0]) + len(bipgraph[1]) == len(graph.keys()):
            Flag = True
            genesU, genesV = bipgraph[0], bipgraph[1]
            geneCoordinate = Get_geneCoordinate(gtf, od)
            fusion_seqs = Get_fusionGene_seq(genesU, genesV, geneCoordinate, ref_genome, od)
            if fusion_seqs == 0:
                print("Sucessfully get the fusion genes sequence ")
            else:
                print("Failed to get the fusion genes sequence")
                exit(1)
    os.system('hisat2-build --quiet {od}/fusiongenes_ref_U.fa {od}/fusiongenes_ref_U'.format(od=od))
    os.system('hisat2-build --quiet {od}/fusiongenes_ref_V.fa {od}/fusiongenes_ref_V'.format(od=od))
    if not Flag:
        print("Failed to build the bipartite graph pairs")
    else:
        print("Sucessfully building the bipartite graph pairs")


if __name__ == '__main__':
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
    usage_string = '''
       Usage:  python build_graph.py --genome <genome-reference-path> --gtf <gtf-reference-path> --tab <fusion-table>
       Arguments:
         Required:
            -g/--genome <genome-reference-path>
              genome reference fasta file 
            -a/--gtf <gtf-reference-path>
              gene annotation GTF format file   
         Optional:
            -t/--tab <fusion-table>
              Inputed fusion gene pairs table,file format should be tsv, fusion pairs Gene pairs should be connected with --, default: reference_fusion_info/fusion_table.tsv
            -o/--outdir <outdir>
              output directory fusion reference information and index 
         Other:
            -h --help
              help information
      '''
    if len(sys.argv[1:]) == 0:
        print(usage_string)
        sys.exit(1)
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hg:a:t:o:", ["help", "genome=", "gtf=", "tab=", "outdir="])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(1)

    genome = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    gtf = 'Homo_sapiens.GRCh38.95.gtf'
    fusion_tab = 'reference_fusion_info/fusion_table.tsv'
    outdir = 'fusion_total_index'

    for opt, arg in opts:
        if opt in ("-h", '--help'):
            print(usage_string)
            sys.exit()

        if opt in ("-g", "--genome"):
            genome = arg
        if opt in ("-a", "--gtf"):
            gtf = arg
        if opt in ("-t", "--tab"):
            fusion_tab = arg
        if opt in ("-o", "--outdir"):
            if arg != '':
                outdir = arg

    print("Outdir is {od}!".format(od=outdir))
    main(fusion_tab, gtf, genome, outdir)
    movement = os.system("cp {ft} {od}/".format(ft=fusion_tab, od=outdir))
    if movement == 0:
        print("All processes done!")
    else:
        print("Fail to move fusion_tab into {od}".format(od=outdir))
        exit(1)
