#!/usr/bin/env python3
import urllib.request
import re
import os
import time
import sys


class Gene():
    def __init__(self, genename):
        '''
        Init a object of gene,including genename and its fasta
        :param genename: a string of genename
        '''
        self.genename = genename
        self.fastaformat = ''

    def __str__(self):
        '''
        Show this gene object
        :return: a string
        '''
        print(self.genename + ' Object')

    def getseqfasta(self):
        '''
        Download the fastq data from Cosmic
        :return: a string of fasta format or None if not exist
        '''
        header = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.96 Safari/537.36'
        }
        try:
            rq = urllib.request.Request(
                "http://cancer.sanger.ac.uk/cosmic/sequence?height=500&ln={gene}&type=cdna&width=700".format(
                    gene=self.genename), headers=header)
            rp = urllib.request.urlopen(rq)
        except:
            # it 's not good to be considered as a robot. If it fail, Wait 5 seconds to download.
            print('Try to download fasta again in 5 seconds later!')
            time.sleep(5)
            rq = urllib.request.Request(
                "http://cancer.sanger.ac.uk/cosmic/sequence?height=500&ln={gene}&type=cdna&width=700".format(
                    gene=self.genename), headers=header)
            rp = urllib.request.urlopen(rq)

        re_div = r'">(.*?)</div>'
        self.fastaformat = re.findall(re_div, rp.read().decode('utf8'), re.S)[0]
        if self.fastaformat != '':
            return self.fastaformat
        else:
            return None

def abs_path(path):
    return os.path.abspath(path.replace('~', os.getenv('HOME')))

def checkexist(pathfile):
    '''
    Check the file exist or not.If exists,remove it.
    :param pathfile: a string of file's path
    :return:None
    '''
    if os.path.exists(pathfile):
        print('The file "{pathfile}" has existed,it will be updated.'.format(pathfile=pathfile))
        os.remove(pathfile)


def getfusionlist():
    '''
    Download the fusion gene information table from Cosmic
    :return:None
    '''
    print('Download the fusion genes from cosmic.')
    rp = urllib.request.urlopen('http://cancer.sanger.ac.uk/cosmic/stats/fusion')
    re_tr = r'<tr class="c">(.*?)</tr>'
    fusion_list = re.findall(re_tr, rp.read().decode('utf8'), re.S)
    checkexist('fusion_table.tsv')
    with open('fusion_table.tsv', 'a') as f:
        f.write('#Update time:{now}, Source:http://cancer.sanger.ac.uk/cosmic/fusion, Total {num} fusion genes.\n'.format(now=str(time.asctime()), num=str(len(fusion_list))))
        f.write('#Genes\tSamples\tMutations\tPapers\tLink\n')
        re_link = r'href="(.*?)"'
        re_other = r'>(\b.*?)<'
        for one in fusion_list:
            line = re.findall(re_other, one, re.S)
            if line != []:
                line.append(re.findall(re_link, one, re.S)[0])
                f.write('\t'.join(line) + '\n')
            else:
                print("Error!")
    print('Finished!The fusion genes are in "{pth}/fusion_table.tsv"!'.format(pth=os.getcwd()))


def createbipgraph(location='fusion_table.tsv', ranked=True):
    '''
    Create a bipartite graph of two fusion groups
    :param location: a string of path
    :return: a list of two fusion gene groups
    '''
    print('Start creating fusion genes bipartite graph!')
    edge = getgraph(location)
    gene_col = getgenecol(location)
    bipgraph = find_bool(edge)
    if bipgraph is not None:
        print('Fusion genes bipartite graph created successfully!')
    if ranked is False:
        return bipgraph
    else:
        print('Start ranking the fusion genes by frequence(sample/mutation/paper)!')
        return rankbykey(*bipgraph, key=gene_col)


def getgraph(location):
    '''
    Load the fusion gene table from the local file
    :param location: a string of fusion gene table location
    :return: a dict of fusion genes
    '''
    with open(location, 'r') as f:
        edge = {}
        nodes = []
        for line in f:
            if not line.startswith('#'):
                temp = line.rstrip().split('\t')[0].split('-')
                start = temp[0]
                end = temp[1]
                if len(temp) == 3:
                    # HLA-A
                    start = temp[0] + '-' + temp[1]
                    end = temp[2]
                if start in edge:
                    edge[start].append(end)
                else:
                    edge[start] = [end]

                if end in edge:
                    edge[end].append(start)
                else:
                    edge[end] = [start]
    return edge


def getgenecol(location, key='Papers'):
    '''

    :param location: a string of fusion gene table location
    :param key:
    :return:
    '''
    Name_to_Col = {'Papers': -2, 'Mutations': -3, 'Samples': -4}
    keyloc = Name_to_Col[key]

    with open(location, 'r') as f:
        sum_gene_key = {}
        for line in f:
            if not line.startswith('#'):
                linelist = line.rstrip().split('\t')
                keynum = int(linelist[keyloc])
                temp = linelist[0].split('-')
                start = temp[0]
                end = temp[1]
                if len(temp) == 3:
                    # HLA-A
                    start = temp[0] + '-' + temp[1]
                    end = temp[2]
                # update sum paper number
                if start in sum_gene_key:
                    sum_gene_key[start] += keynum
                else:
                    sum_gene_key[start] = keynum
                if end in sum_gene_key:
                    sum_gene_key[end] += keynum
                else:
                    sum_gene_key[end] = keynum
    return sum_gene_key


def find_bool(graph):
    '''

    :param graph:
    :return:
    '''
    thisside = set()
    thatside = set()
    nodes = list(graph.keys())
    thisside.add(nodes[0])

    while(nodes != []):
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
    return [list(thisside), list(thatside)]


def rankbykey(list1, list2, key):
    '''
    Rank the fusion gene according the key col
    :param list1: first group of fusion genes
    :param list2: second group of fusion genes
    :param key: a string of col name
    :return: a list of two ranked fusion gene groups
    '''
    # sorted by key, reported papers is better
    list1 = sorted([tuple((key[gene], gene)) for gene in list1], reverse=True)
    list2 = sorted([tuple((key[gene], gene)) for gene in list2], reverse=True)
    list1 = [one[1] for one in list1]
    list2 = [one[1] for one in list2]
    return list1, list2


def getreferences(list1, list2, breaktime=0):
    '''
    Download the fasta file in bulk
    :param list1: first group of fusion genes
    :param list2: second group of fusion genes
    :return: None
    '''
    print('Start getting the reference fasta of fusion genes!')
    checkexist('fusiongenes_ref_U.fa')
    checkexist('fusiongenes_ref_V.fa')
    with open('fusiongenes_ref_U.fa', 'a') as f1:
        for k in range(len(list1)):
            tempgene = Gene(list1[k])
            tempgene.getseqfasta()
            if k != len(list1):
                f1.write(tempgene.fastaformat)
            else:
                f1.write(tempgene.fastaformat[:-1])
            time.sleep(breaktime)
    with open('fusiongenes_ref_V.fa', 'a') as f2:
        for l in range(len(list2)):
            tempgene = Gene(list2[l])
            tempgene.getseqfasta()
            if l != len(list2):
                f2.write(tempgene.fastaformat)
            else:
                f2.write(tempgene.fastaformat[:-1])
            time.sleep(breaktime)
    print('The reference fasta of fusion genes downloaded successfully!')
    print('Finished!The references are in "{pth}/fusiongenes_ref_U(or V).fa"!'.format(pth=os.getcwd()))


if __name__ == '__main__':
    # download fusion table
    argvs = sys.argv[1:]
    if len(argvs)==0 or len(argvs)>2:
        print('''Usage:     python downloadfusion.py <fusion_dir> [<breaktime>]
Arguments:
        <fusion_dir>
            create a directory of fusion genes and fasta data
    Optional:
        <breaktime> 
            seconds of break time when download fasta file, if the program has fault, try a big number(default:1)
    ''')
        exit()
    abspath = abs_path(argvs[0])
    if not os.path.isdir(abspath):
        try:
            os.mkdir(abspath)
        except:
            print("{path} exists, it will be rewrited.".format(path=abspath))
    os.chdir(abspath)
    getfusionlist()
    # get bip graph(sorted by col,default paper number,descending order)
    ranked_genelist1, ranked_genelist2 = createbipgraph()
    #print(ranked_genelist1)
    #print(ranked_genelist2)
    # download fasta reference
    breaktime = int(argvs[1]) if len(argvs) == 2 else 1
    getreferences(ranked_genelist1, ranked_genelist2, breaktime=breaktime)
