#!/usr/bin/env python3
import pysam
import sys
import scipy.stats

def test_uo(lengths,u):
    return scipy.stats.wilcoxon([x-u for x in lengths],zero_method="wilcox").pvalue


def verity_distribution(reads,inferred_inform):
    breakpointreads = []
    scanreads = []
    split_lengths = []
    read_lengths = []
    pvalue = -1

    breakpoint = int(inferred_inform[0])
    for read in reads:
        if read.flag & 4:
            continue
        if read.reference_start < breakpoint and read.reference_end > breakpoint:
            if read.reference_length > read.query_length + 5:
                scanreads.append(read)
            else:
                breakpointreads.append(read)
                split_lengths += [breakpoint - read.reference_start, read.reference_end - breakpoint]
                read_lengths.append(read.query_length)
    if breakpointreads != []:
        u_hat = sum(read_lengths)/len(read_lengths)/2
        pvalue = test_uo(split_lengths,u_hat)
    return breakpointreads, scanreads, pvalue


def write_fusion(inferred_sam, ifr="temp/inferred_fusion_results.tsv", out="fusion_results.tsv"):
    file = tuple(pysam.AlignmentFile(inferred_sam, 'r'))
    if len(file) == 0:
        return None

    fusion_reads = {}
    for read in file:
        if not read.flag & 4:
            temp = read.reference_name.split('-')
            fusion5 = temp[0]
            fusion3 = temp[1]
            if len(temp) == 3:
                # HLA-A
                fusion5 = temp[0] + '-' + temp[1]
                fusion3 = temp[2]
            if (fusion5, fusion3) not in fusion_reads:
                fusion_reads[(fusion5, fusion3)] = set()
            fusion_reads[(fusion5, fusion3)].add(read)

    fusion_result = {}
    string = '\t'.join((
         "#Fusion Name",
         "5'Gene",
         "3'Gene",
         "5'Gene BreakPoint Pos",
         "3'Gene BreakPoint Pos",
         "5'Gene Breakpoint Seq",
         "3'Gene Breakpoint Seq",
         "5'and 3'Common Breakpoint Seq",
         "BreakpointReads Count",
         "BreakpointReads",
         "BreakpointStrand Count(+,-)",
         "ScanningReads Count",
         "ScanningReads",
         "ScanningStrand Count(+,-)",
         "Fusion Seq Length"
         "P-Value"
    )) + '\n'
    result = ''
    with open(ifr,'r') as f:
        for line in f:
            # temp: fusion5,fusion3,fusion5pos,fusion3pos,fusion5seq,fusion3seq,commonseq
            temp = line.rstrip().split('\t')
            fusionpair = (temp[0],temp[1])
            fusion_result[fusionpair] = temp[2:]
            if fusionpair not in fusion_reads:
                continue
            bp, sc, p = verity_distribution(fusion_reads[fusionpair],fusion_result[fusionpair])
            if p != -1:
                bpminu = sum([1 for x in bp if x.flag & 16 == 16])
                bpplus = len(bp) - bpminu
                scminu = sum([1 for x in sc if x.flag & 16 == 16])
                scplus = len(sc) - scminu
                result = '\t'.join((
                    '-'.join(fusionpair),
                    fusionpair[0],
                    fusionpair[1],
                    fusion_result[fusionpair][0],
                    fusion_result[fusionpair][1],
                    fusion_result[fusionpair][2],
                    fusion_result[fusionpair][3],
                    fusion_result[fusionpair][4],
                    str(len(bp)),
                    ','.join([x.query_name for x in bp]),
                    ','.join((str(bpplus), str(bpminu))),
                    str(len(sc)),
                    ','.join([x.query_name for x in sc]),
                    ','.join((str(scplus), str(scminu))),
                    fusion_result[fusionpair][5],
                    str(p)
                )) + '\n'
                string += result
    with open(out,'w') as f:
        f.write(string)

    if result != '':
        return 0

    return None

def write_fcirc(circ_sam, fragmentcount, fr="fusion_results.tsv", out="fcircRNA_results.tsv"):
    file = pysam.AlignmentFile(circ_sam, 'r')
    circ = {}
    for one in file:
        if one.flag & 4:
            continue
        fusionname = one.reference_name
        start = one.reference_start
        end = one.reference_end
        if (fusionname, start, end) not in circ:
            circ[(fusionname, start, end)] = [one]
        else:
            circ[(fusionname, start, end)].append(one)
    if circ == {}:
        return None

    fusion = {}
    with open(fr,'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            temp = line.rstrip().split('\t')
            fusion[temp[0]] = temp

    string = '\t'.join((
        "#FcircRNA_NO",
        "Fusion Name",
        "Backsplice start",
        "Backsplice end",
        "FusionBreakPoint Pos",
        "FusionSeq Length",
        "Support FcircRNA Reads Count",
        "Support FcircRNA Reads",
        "FcircRNA Strand Count(+,-)",
        "FCI"
    )) + '\n'

    fcirccount = 0
    for key in circ.keys():
        fcirccount += 1
        circminu = sum([1 for one in circ[key] if one.flag & 16 == 16])
        circplus = len(circ[key]) - circminu

        string += '\t'.join((
            "No_" + str(fcirccount),
            key[0],
            str(key[1]),
            str(key[2]),
            fusion[key[0]][3],
            fusion[key[0]][-2],
            str(circminu + circplus),
            ','.join([x.query_name for x in circ[key]]),
            ','.join((str(circplus), str(circminu))),
            str((circminu + circplus)/int(fusion[key[0]][-2])/int(fragmentcount))
        )) + '\n'

    with open(out,'w') as f:
        f.write(string)
    print("Find {num} kind(s) of fcircRNAs!".format(num=fcirccount))
    return 0

if __name__ == "__main__":
    #write_fusion(sys.argv[1])
    write_fcirc(sys.argv[1], sys.argv[2])