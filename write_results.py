#!/usr/bin/env python3
import pysam
import sys
import scipy.stats


def test_uo(lengths, u):
    return scipy.stats.wilcoxon([x - u for x in lengths], zero_method="wilcox").pvalue


def verity_distribution(reads, inferred_inform):
    breakpointreads = []
    scanreads = []  # scanreads
    split_lengths = []
    read_lengths = []
    pvalue = -1

    breakpoint = int(inferred_inform[3])  # inferred_inform[0]
    for read in reads:
        if read.flag & 4:
            continue
        if read.reference_start < breakpoint and read.reference_end > breakpoint:
            if 4 in read.cigartuples[0][:] or 4 in read.cigartuples[-1][:]:
                if read.cigartuples[0][1] > 3 or read.cigartuples[-1][1] > 3:
                    scanreads.append(read)
            else:
                if abs(read.reference_start - breakpoint) > 10 and abs(read.reference_end - breakpoint) > 10:
                    breakpointreads.append(read)
                    split_lengths += [breakpoint - read.reference_start, read.reference_end - breakpoint]
                    read_lengths.append(read.query_length)
    if breakpointreads != []:
        u_hat = sum(read_lengths) / len(read_lengths) / 2
        pvalue = test_uo(split_lengths, u_hat)
    return breakpointreads, scanreads, pvalue


def filter_similar(dict_sam):
    Sam_dict = dict_sam
    keys = list(Sam_dict.keys())
    fusions = [i[0] for i in keys]
    for fusion in set(fusions):
        if fusions.count(fusion) >= 2:
            Target_key = [key for key in keys if fusion in key]
            for m in range(len(Target_key)):
                for n in range(m + 1, len(Target_key)):
                    if abs(int(Target_key[m][1]) - int(Target_key[n][1])) > 50 or abs(
                            int(Target_key[m][2]) - int(Target_key[n][2])) > 50:
                        continue
                    else:
                        if len(Sam_dict[Target_key[m]]) > len(Sam_dict[Target_key[n]]):
                            for i in Sam_dict[Target_key[n]]:
                                (Sam_dict[Target_key[m]]).append(i)
                            Sam_dict[Target_key[n]] = ''
                        else:
                            for i in Sam_dict[Target_key[m]]:
                                (Sam_dict[Target_key[n]]).append(i)
                            Sam_dict[Target_key[m]] = ''
    New_dict = {}
    for fusion in keys:
        if Sam_dict[fusion] is not '':
            New_dict[fusion] = Sam_dict[fusion]
    return New_dict


def write_fusion(inferred_sam, ifr="temp/inferred_fusion_results.tsv", out="fusion_results.tsv"):
    file = tuple(pysam.AlignmentFile(inferred_sam, 'r'))
    if len(file) == 0:
        return -1

    fusion_reads = {}
    for read in file:
        if not read.flag & 4:
            temp = read.reference_name.split('-')
            fusion5 = temp[0]
            fusion3 = temp[1]
            if len(temp) == 3:
                fusion5 = temp[0] + '-' + temp[1]
                fusion3 = temp[2]
            if (fusion5, fusion3) not in fusion_reads:
                fusion_reads[(fusion5, fusion3)] = set()
            fusion_reads[(fusion5, fusion3)].add(read)

    fusion_result = {}
    string = '\t'.join((
        "#Fusion_Name",
        "5'Gene",
        "3'Gene",
        "5'Gene_chr",
        "5'Gene_strand",
        "3'Gene_chr",
        "3'Gene_strand",
        "5'Gene_BreakPoint_Pos",
        "3'Gene_BreakPoint_Pos",
        "5'and3'_Common_Breakpoint_Seq",
        "BreakpointReads_Count",
        "BreakpointReads",
        "BreakpointStrand_Count(+,-)",
        "ScanningReads_Count",
        "ScanningReads",
        "ScanningStrand_Count(+,-)",
        "P-Value"
    )) + '\n'
    # "5'Gene Breakpoint Seq","3'Gene Breakpoint Seq" "Fusion Seq Length"
    result = ''
    with open(ifr, 'r') as f:
        for line in f:
            # temp: f5_chr f5_strand  f5_coords f5_pos f3_chr f3_strand f3_coords f3_pos f5_seq f3_seq f53_com seqlen
            temp = line.rstrip().split('\t')
            fusionpair = (temp[0], temp[1])
            fusion_result[fusionpair] = temp[2:]
            if fusionpair not in fusion_reads:
                continue
            bp, sc, p = verity_distribution(fusion_reads[fusionpair], fusion_result[fusionpair])
            #
            if p != -1:
                bpminu = sum([1 for x in bp if x.flag & 16 == 16])
                bpplus = len(bp) - bpminu
                scminu = sum([1 for x in sc if x.flag & 16 == 16])
                scplus = len(sc) - scminu
                strand5=[ '+' if i=='1' else '-' for i in fusion_result[fusionpair][1]][0]
                strand3=['+' if j=='1' else '-' for j in  fusion_result[fusionpair][5]][0]
                f5_start,f5_end=fusion_result[fusionpair][2].split(':')
                f3_start,f3_end=fusion_result[fusionpair][6].split(':')

                if strand5=='+':
                    breakpos5=str(int(f5_start) + int(fusion_result[fusionpair][3]))
                else:
                    breakpos5 = str(int(f5_end) - int(fusion_result[fusionpair][3]))
                if strand3=='+':
                    breakpos3=str(int(f3_start) + int(fusion_result[fusionpair][7]))
                else:
                    breakpos3=str(int(f3_end) - int(fusion_result[fusionpair][7]))

                result = '\t'.join((
                    '-'.join(fusionpair),
                    fusionpair[0],
                    fusionpair[1],
                    fusion_result[fusionpair][0],
                    strand5,
                    fusion_result[fusionpair][4],
                    strand3,
                    breakpos5,
                    breakpos3,
                    fusion_result[fusionpair][10],
                    str(len(bp)),
                    ','.join([x.query_name for x in bp]),
                    ','.join((str(bpplus), str(bpminu))),
                    str(len(sc)),
                    ','.join([x.query_name for x in sc]),
                    ','.join((str(scplus), str(scminu))),
                    str(p)
                )) + '\n'
                string += result
    with open(out, 'w') as f:
        f.write(string)

    if result == '':
        return -1

    return fusion_result


def write_fcirc(circ_sam, fragmentcount, fusion_result, fr="fusion_results.tsv", out="fcircRNA_results.tsv"):
    file = pysam.AlignmentFile(circ_sam, 'r')
    circ = {}
    for one in file:
        if one.flag & 4:
            continue
        fusionname = one.reference_name
        start = one.reference_start
        end = one.reference_end
        if (fusionname, start, end) not in circ:
            circ[fusionname, start, end] = [one]
        else:
            circ[fusionname, start, end].append(one)
    circ = filter_similar(circ)
    if circ == {}:
        return -1
    fusion = {}
    with open(fr, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            temp = line.rstrip().split('\t')
            fusion[temp[0]] = temp

    string = '\t'.join((
        "#FcircRNA_NO",
        "Fusion Name",
        "Backsplice_start",
        "Backsplice_end",
        "Fusion5'_BreakPoint_Pos",
        "Fusion3'_BreakPoint_Pos",
        "Support_FcircRNA_Reads_Count",
        "FcircRNA_Strand_Count(+, -)",
        "Support_FcircRNA_Reads"
    )) + '\n'
    # "FusionSeq Length",  "FCI"
    fcirccount = 0
    for key in circ.keys():
        if key[0] in fusion.keys():
            fcirccount += 1
            circminu = sum([1 for one in circ[key] if one.flag & 16 == 16])
            circplus = len(circ[key]) - circminu
            strand_bkstart=fusion[key[0]][4]
            strand_bkend=fusion[key[0]][6]
            f5_start,f5_end=fusion_result[tuple(key[0].split('-'))][2].split(':')
            f5_total=int(f5_end)-int(f5_start)
            f3_start,f3_end=fusion_result[tuple(key[0].split('-'))][6].split(':')
            f3_total=int(f3_end)-int(f3_start)
            if strand_bkstart=='+':
                bkstart_pos=str(int(f5_start)+int(key[1]))
            else:
                bkstart_pos=str(int(f5_start)+f5_total-int(key[1]))
            if strand_bkend=='+':
                bkend_pos=str(int(fusion[key[0]][8])+(int(key[2])-int(fusion_result[tuple(key[0].split('-'))][3])))
            else:
                bkend_pos=str(int(fusion[key[0]][8])-(int(key[2])-int(fusion_result[tuple(key[0].split('-'))][3])))
            string += '\t'.join((
            "No_" + str(fcirccount),
            key[0],
            fusion[key[0]][3] + ":" + bkstart_pos+":"+strand_bkstart,
            fusion[key[0]][5] + ":" + bkend_pos+':'+strand_bkend,
            fusion[key[0]][3] + ":" + fusion[key[0]][7] +':'+fusion[key[0]][4],
            fusion[key[0]][5] + ":" + fusion[key[0]][8] +':'+fusion[key[0]][6],
            str(circminu + circplus),
            ','.join((str(circplus), str(circminu))),
            ','.join([x.query_name for x in circ[key]]),
            )) + '\n'  #str((circminu + circplus) / int(fusion[key[0]][-2]) / int(fragmentcount))
            with open(out, 'w') as f:
                f.write(string)
        else:
            continue  # There is a circRNA But not a Fusion_circRNA!
    print("Find {num} kind(s) of fcircRNAs!".format(num=fcirccount))
    return 0



if __name__ == "__main__":
    # write_fusion(sys.argv[1])
    write_fcirc(sys.argv[1], sys.argv[2])

