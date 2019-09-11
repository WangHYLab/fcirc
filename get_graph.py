#!/usr/bin/env python3
import re
import os
import time
import sys

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


if __name__ == '__main__':
    getgraph(sys.argv[1])
    find_bool(sys.argv[2])
