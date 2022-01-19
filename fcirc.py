#!/usr/bin/env python3
# -*- coding: utf8 -*-
import sys
import getopt
import os
import time
# import logging
from filter_reads import drop_mapped_single_tofastq, drop_unmapped_single_tosam, drop_mapped_paired_tofastq, \
    drop_unmapped_paired_tosam, sam_single_tofastq, sam_paired_tofastq, count_reads
from quality_control import quality_filter_SE, quality_filter_PE
from getpairedreads import getpairedreads
from reconstruction import reconstruct
from transform import transform_sam
from write_results import write_fusion, write_fcirc


def version():
    print("fcirc (fusion circRNA identifier)    Version: 1.0.2")


def usage():
    usage_string = '''
Program:    fcirc (fusion circRNA identifier)
Version:    1.0.2(written by python3)
Contact:    Zhaoqing Cai <caizhaoqingeawt@163.com> HongZhang Xue <xuehzh95@foxmail.com>

Usage:      python fcirc.py [options] -x <ht2-trans-idx> -f <ht2-fusion-idx-dir> {-1 <fastq1> | -1 <fastq1> -2 <fastq2>}

Arguments:
    Required:
        -x <ht2-trans-idx>, --trans_idx <ht2-trans-idx>
            transcription index filename prefix (minus trailing .X.ht2)
        -f <ht2-fusion-idx-dir>, --fusion_idx_dir <ht2-fusion-idx-dir>
            fusion index directory (contains fusiongenes_ref_U and fusiongenes_ref_V)
        -1 <fastq1>, --file1 <fastq1>
            fastq file 1 (single-end pattern:only -1)
        -2 <fastq2>, --file2 <fastq2>
            fastq file 2 (paired-end pattern:-1 and -2, files should be like -1 xxx_1.fastq -2 xxx_2.fastq)

    Optional:
        -q <quality_val>
            the minimum phred qulaity of read(default:0)
        -c <fusion-genes-coordinates> --fusion_genes_coord
            fusion genes coordinates file(defalut: fusion_genes_coordinate.txt in fusion index directory)
        -o <output_dir>, --output <output_dir>
            output file directory (default: .)
        -t <int>, --thread <int>
            number of hisat2 alignment and pysam filter threads to launch (default:1)

    Others:
        -h, --help
            help information
        -v, --version
            version information
    '''
    print(usage_string)


def abs_path(path):
    return os.path.abspath(os.path.join(os.curdir, path.replace('~', os.getenv('HOME'))))


def str_current_time():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


def main(argvs):
    # logging.basicConfig(level=logging.INFO,
    #                     format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s')
    if len(argvs) == 0:
        usage()
        sys.exit()

    try:
        opts, args = getopt.getopt(argvs, "hv1:2:x:f:c:q:o:t:",
                                   ["help", "version", "file1=", "file2=", "trans_idx=", "fusion_idx_dir=",
                                    "fusion_genes_coord=", "output=", "thread="])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    # default arguments:
    quality = 0
    thread = '1'
    output_path = os.getcwd()

    # required arguments:
    fastq1_path = ''
    fastq2_path = ''
    trans_idx_path = ''
    fusion_idx_dir_path = ''
    pattern = "single-end"


    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-v', '--version'):
            version()
            sys.exit()

        if opt in ('-1', '--file1'):
            pattern = "single-end"
            abs_arg = abs_path(arg)
            if os.path.isfile(abs_arg):
                fastq1_path = abs_arg

        if opt in ('-2', '--file2'):
            pattern = "paired-end"
            abs_arg = abs_path(arg)
            if os.path.isfile(abs_arg):
                fastq2_path = abs_arg

        if opt in ('-x', '--trans_idx'):
            abs_arg = abs_path(arg)
            trans_idx_path = abs_arg
            if not os.path.isfile(trans_idx_path + '.1.ht2'):
                print("Error:Fail to locate hisat2 transcription index prefix, please check {path}".format(
                    path=trans_idx_path))
                sys.exit(1)

        if opt in ('-f', '--fusion_idx_dir'):
            abs_arg = abs_path(arg)
            if os.path.isdir(abs_arg):
                fusion_idx_dir_path = abs_arg.rstrip('/')
            elif fusion_idx_dir_path == '':
                print("Error:Fail to locate hisat2 fusion index prefix, please check {path}".format(
                    path=fusion_idx_dir_path))
                sys.exit(1)

        fusion_genes_coord = os.path.join(fusion_idx_dir_path, 'fusion_genes_coordinate.txt')

        if opt in ('-c', '--fusion_genes_coord'):
            abs_arg = abs_path(arg)
            if os.path.isfile(abs_arg):
                fusion_genes_coord = abs_arg
            else:
                print('Error to locate fusion genes coordinate file, please check{path}'.format(path=fusion_genes_coord))
                sys.exit(1)

        if opt == '-q':
            qulaity = int(arg)
            if qulaity < 0:
                print("Error: The quality value can not be less than 0")
                exit(1)

        if opt in ('-o', '--output'):
            abs_arg = abs_path(arg)
            if not os.path.isdir(abs_arg):
                try:
                    os.mkdir(abs_arg)
                except:
                    print("Error:Fail to make dir at {path}".format(path=arg))
                    sys.exit(1)
            output_path = arg.rstrip('/')
        if opt in ('-t', '--thread'):
            thread = arg


    try:
        os.mkdir("{od}/temp".format(od=output_path))
    except:
        print("temp dir exists, it will be rewrited.")

    # start
    os.chdir(output_path)
    # single-end pattern
    if pattern == "single-end":
        if fastq1_path == '':
            print("Error:Fail to locate fastq file,please check the path!")
            sys.exit(1)

        projectname = fastq1_path.split('/')[-1].rstrip('.fastq').rstrip('.fq')

        hisat2_trans_command = '''hisat2 -p {thread} \
            --mp 6,2 \
            --rdg 5,3 \
            --rfg 5,3 \
            --no-softclip \
            -x {hisat2_trans_idx} \
            -U {singlefastq} \
            -S {output_transsam} \
            2>{output_translog}'''.format(thread=str(thread),
                                          hisat2_trans_idx=trans_idx_path,
                                          singlefastq=fastq1_path,
                                          output_transsam="temp/trans_" + projectname + ".sam",
                                          output_translog="temp/trans_" + projectname + ".log")
        return_state = os.system(hisat2_trans_command)
        if return_state != 0:
            print("[{now}] Error mapping reads to transcription!".format(now=str_current_time()))
            exit(return_state)

        print("[{now}] Finish mapping reads to transcription!".format(now=str_current_time()))
        unmapped_fastq_path = drop_mapped_single_tofastq("temp/trans_" + projectname + ".sam")

        if quality > 0:
            quality_cmd = quality_filter_SE(quality, unmapped_fastq_path, 'single-end')
            if quality_cmd != 0:
                print("Error: Fail to filter quality reads")
                exit(1)

        hisat2_fusionU_command = '''hisat2 -p {thread} \
            --mp 6,2 \
            --rdg 5,2 \
            --rfg 5,2 \
            --sp 1,1 \
            --score-min L,0,-0.8 \
            -x {hisat2_fusion_idx}/fusiongenes_ref_U \
            -U {unmappedfastq} \
            -S {output_fusionsam} \
            2>{output_fusionlog}'''.format(thread=str(thread),
                                           hisat2_fusion_idx=fusion_idx_dir_path,
                                           unmappedfastq=unmapped_fastq_path,
                                           output_fusionsam="temp/fusionU_" + projectname + ".sam",
                                           output_fusionlog="temp/fusionU_" + projectname + ".log")
        return_state = os.system(hisat2_fusionU_command)
        if return_state != 0:
            exit(return_state)
        print("[{now}] Finish mapping reads to fusion references U!".format(now=str_current_time()))

        hisat2_fusionV_command = '''hisat2 -p {thread} \
            --mp 6,2 \
            --rdg 5,2 \
            --rfg 5,2 \
            --sp 1,1 \
            --score-min L,0,-0.8 \
            -x {hisat2_fusion_idx}/fusiongenes_ref_V \
            -U {unmappedfastq} \
            -S {output_fusionsam} \
            2>{output_fusionlog}'''.format(thread=str(thread),
                                           hisat2_fusion_idx=fusion_idx_dir_path,
                                           unmappedfastq=unmapped_fastq_path,
                                           output_fusionsam="temp/fusionV_" + projectname + ".sam",
                                           output_fusionlog="temp/fusionV_" + projectname + ".log")
        return_state = os.system(hisat2_fusionV_command)
        if return_state != 0:
            exit(return_state)
        print("[{now}] Finish mapping reads to fusion references V!".format(now=str_current_time()))

        mapped_samU_path = drop_unmapped_single_tosam("temp/fusionU_" + projectname + ".sam", thread=str(thread))
        mapped_samV_path = drop_unmapped_single_tosam("temp/fusionV_" + projectname + ".sam", thread=str(thread))
        print("[{now}] Finish dropping unmapped read in fusion references U and V!".format(now=str_current_time()))

        mapped_filtered_samU_path, mapped_filtered_samV_path = getpairedreads(mapped_samU_path, mapped_samV_path,
                                                                              pattern, fusion_idx_dir_path)
        print("[{now}] Finish filtering fusion-related reads in fusion references U and V!".format(
            now=str_current_time()))

        return_state = reconstruct(mapped_filtered_samU_path, mapped_filtered_samV_path, fusion_genes_coord,
                                   fusion_idx_dir_path)
        if return_state != 0:
            print("[{now}] No fusion-related read is detected!".format(now=str_current_time()))
            exit(return_state)

        # reads both exists in U and V
        sam_single_tofastq(mapped_filtered_samV_path, "temp/inferred_fusion.fastq")

        hisat2_verify_fusion_command = '''hisat2 -p {thread} \
            --mp 6,2 \
            --rdg 5,2 \
            --rfg 5,2 \
            --sp 1,1 \
            --score-min L,0,-0.8 \
            -x {inferred_fusion_fasta} \
            -U {inferred_fusion_fastq} \
            -S {output_inferred_fusionsam} \
            2>{output_inferred_fusionlog}'''.format(thread=str(thread),
                                                    inferred_fusion_fasta="temp/inferred_fusion.fa",
                                                    inferred_fusion_fastq="temp/inferred_fusion.fastq",
                                                    output_inferred_fusionsam="temp/inferred_fusion.sam",
                                                    output_inferred_fusionlog="temp/inferred_fusion.log")
        return_state = os.system(hisat2_verify_fusion_command)
        if return_state != 0:
            exit(return_state)
        print("[{now}] Finish mapping reads to inferred fusion references!".format(now=str_current_time()))

        # check the distribution near the fusion break point
        fusion_result = write_fusion("temp/inferred_fusion.sam")
        if fusion_result == -1:
            print("[{now}] No fusion gene is detected!".format(now=str_current_time()))
            exit(return_state)

        # filter and get the circ-related reads
        transform_sam("temp/inferred_fusion.sam")
        sam_single_tofastq("temp/inferred_fusion_tranformed.sam", "temp/inferred_fusion_tranformed.fastq")

        hisat2_fusion_circ_command = '''hisat2 -p {thread} \
            --no-softclip \
            --pen-noncansplice 0 \
            -x {inferred_fusion_fasta} \
            -U {inferred_fusion_circ_fastq} \
            -S {output_inferred_fusion_circsam} \
            2>{output_inferred_fusion_circlog}'''.format(thread=str(thread),
                                                         inferred_fusion_fasta="temp/inferred_fusion.fa",
                                                         inferred_fusion_circ_fastq="temp/inferred_fusion_tranformed.fastq",
                                                         output_inferred_fusion_circsam="temp/inferred_fusion_circ.sam",
                                                         output_inferred_fusion_circlog="temp/inferred_fusion_circ.log")
        return_state = os.system(hisat2_fusion_circ_command)
        if return_state != 0:
            exit(return_state)
        # count input reads
        readcount = count_reads("temp/trans_" + projectname + ".sam")
        return_state = write_fcirc("temp/inferred_fusion_circ.sam", readcount, fusion_result)
        if return_state == -1:
            print("[{now}] No fusion circRNA is detected!".format(now=str_current_time()))
            exit(return_state)
        print("[{now}] Finish all!See the result in 'fcircRNA_results.tsv'!".format(now=str_current_time()))







    # paired-end pattern
    elif pattern == "paired-end":
        if fastq1_path == '' or fastq2_path == '':
            print("Error:Fail to locate fastq files,please check the path")
            sys.exit(1)

        projectname = fastq1_path.split('/')[-1].split('_1')[0]

        hisat2_trans_command = '''hisat2 -p {thread} \
            --mp 6,2 \
            --rdg 5,3 \
            --rfg 5,3 \
            --no-softclip \
            -x {hisat2_trans_idx} \
            -1 {mate1fastq} \
            -2 {mate2fastq} \
            -S {output_transsam} \
            2>{output_translog}'''.format(thread=str(thread),
                                          hisat2_trans_idx=trans_idx_path,
                                          mate1fastq=fastq1_path,
                                          mate2fastq=fastq2_path,
                                          output_transsam="temp/trans_" + projectname + ".sam",
                                          output_translog="temp/trans_" + projectname + ".log")
        return_state = os.system(hisat2_trans_command)
        if return_state != 0:
            print("[{now}] Error: Fail mapping reads to transcription!".format(now=str_current_time()))
            exit(return_state)
        print("[{now}] Finish mapping reads to transcription!".format(now=str_current_time()))

        unmapped_fastq1_path, unmapped_fastq2_path = drop_mapped_paired_tofastq("temp/trans_" + projectname + ".sam")
        if quality > 0:
            quality_cmd = quality_filter_PE(unmapped_fastq1_path, unmapped_fastq2_path, quality, 'pair-end')
            if quality_cmd != 0:
                print("Error: Fail to filter quality reads")
                exit(1)

        hisat2_fusionU_command = '''hisat2 -p {thread} \
            --mp 12,4 \
            --rdg 10,4 \
            --rfg 10,4 \
            --sp 1,1 \
            --score-min L,0,-0.8 \
            -x {hisat2_fusion_idx}/fusiongenes_ref_U \
            -1 {mate1fastq} \
            -2 {mate2fastq} \
            -S {output_fusionsam} \
            2>{output_fusionlog}'''.format(thread=str(thread),
                                           hisat2_fusion_idx=fusion_idx_dir_path,
                                           mate1fastq=unmapped_fastq1_path,
                                           mate2fastq=unmapped_fastq2_path,
                                           output_fusionsam="temp/fusionU_" + projectname + ".sam",
                                           output_fusionlog="temp/fusionU_" + projectname + ".log")
        return_state = os.system(hisat2_fusionU_command)
        if return_state != 0:
            print("[{now}] Error: Fail mapping reads to fusion references U!".format(now=str_current_time()))
            exit(return_state)
        print("[{now}] Finish mapping reads to fusion references U!".format(now=str_current_time()))

        hisat2_fusionV_command = '''hisat2 -p {thread} \
            --mp 12,4 \
            --rdg 10,4 \
            --rfg 10,4 \
            --sp 1,1 \
            --score-min L,0,-0.8 \
            -x {hisat2_fusion_idx}/fusiongenes_ref_V \
            -1 {mate1fastq} \
            -2 {mate2fastq} \
            -S {output_fusionsam} \
            2>{output_fusionlog}'''.format(thread=str(thread),
                                           hisat2_fusion_idx=fusion_idx_dir_path,
                                           mate1fastq=unmapped_fastq1_path,
                                           mate2fastq=unmapped_fastq2_path,
                                           output_fusionsam="temp/fusionV_" + projectname + ".sam",
                                           output_fusionlog="temp/fusionV_" + projectname + ".log")
        return_state = os.system(hisat2_fusionV_command)
        if return_state != 0:
            print("[{now}] Error: Fail mapping reads to fusion references V!".format(now=str_current_time()))
            exit(return_state)
        print("[{now}] Finish mapping reads to fusion references V!".format(now=str_current_time()))

        mapped_samU_path = drop_unmapped_paired_tosam("temp/fusionU_" + projectname + ".sam", thread=str(thread))
        mapped_samV_path = drop_unmapped_paired_tosam("temp/fusionV_" + projectname + ".sam", thread=str(thread))

        mapped_filtered_samU_path, mapped_filtered_samV_path = getpairedreads(mapped_samU_path, mapped_samV_path,
                                                                              pattern, fusion_idx_dir_path)
        print("[{now}] Finish filtering fusion-related reads in fusion references U and V!".format(
            now=str_current_time()))

        return_state = reconstruct(mapped_filtered_samU_path, mapped_filtered_samV_path, fusion_genes_coord,
                                   fusion_idx_dir_path)
        if return_state != 0:
            print("[{now}] No fusion-related read is detected!".format(now=str_current_time()))
            exit(return_state)

        # reads both exists in U and V
        sam_paired_tofastq(mapped_filtered_samV_path, "temp/inferred_fusion_1.fastq", "temp/inferred_fusion_2.fastq")

        hisat2_verify_fusion_command = '''hisat2 -p {thread} \
                    --mp 6,2 \
                    --rdg 5,2 \
                    --rfg 5,2 \
                    --sp 1,1 \
                    --score-min L,0,-0.8 \
                    -x temp/inferred_fusion.fa \
                    -1 {inferred_fusion1_fastq} \
                    -2 {inferred_fusion2_fastq} \
                    -S {output_inferred_fusionsam} \
                    2>{output_inferred_fusionlog}'''.format(thread=str(thread),
                                                            inferred_fusion1_fastq="temp/inferred_fusion_1.fastq",
                                                            inferred_fusion2_fastq="temp/inferred_fusion_2.fastq",
                                                            output_inferred_fusionsam="temp/inferred_fusion.sam",
                                                            output_inferred_fusionlog="temp/inferred_fusion.log")
        return_state = os.system(hisat2_verify_fusion_command)
        if return_state == -1:
            print("[{now}] No fusion gene is detected!".format(now=str_current_time()))
            exit(return_state)
        print("[{now}] Finish mapping reads to inferred fusion references!".format(now=str_current_time()))

        # check the distribution near the fusion break point
        fusion_result = write_fusion("temp/inferred_fusion.sam")
        if fusion_result == -1:
            exit(return_state)

        # filter and get the circ-related reads
        transform_sam("temp/inferred_fusion.sam")
        sam_single_tofastq("temp/inferred_fusion_tranformed.sam", "temp/inferred_fusion_tranformed.fastq")
        sam_single_tofastq("temp/inferred_fusion_tranformed.sam", "temp/inferred_fusion_tranformed.fastq")

        hisat2_fusion_circ_command = '''hisat2 -p {thread} \
                    --no-softclip \
                    --pen-noncansplice 0 \
                    -x temp/inferred_fusion.fa \
                    -U {inferred_fusion_circ_fastq} \
                    -S {output_inferred_fusion_circsam} \
                    2>{output_inferred_fusion_circlog}'''.format(thread=str(thread),
                                                                 inferred_fusion_circ_fastq="temp/inferred_fusion_tranformed.fastq",
                                                                 output_inferred_fusion_circsam="temp/inferred_fusion_circ.sam",
                                                                 output_inferred_fusion_circlog="temp/inferred_fusion_circ.log")
        return_state = os.system(hisat2_fusion_circ_command)
        if return_state != 0:
            exit(return_state)
        # count input reads
        readcount = count_reads("temp/trans_" + projectname + ".sam")
        return_state = write_fcirc("temp/inferred_fusion_circ.sam", readcount, fusion_result)
        if return_state == -1:
            print("[{now}] No fusion circRNA is detected!".format(now=str_current_time()))
            exit(return_state)
        print("[{now}] Finish all!See the result in 'fcircRNA_results.tsv'!".format(now=str_current_time()))


if __name__ == "__main__":
    #logging.basicConfig(level=logging.INFO,
    #                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s')
    print("[{now}] Start running # {command}".format(now=str_current_time(), command=' '.join(sys.argv)))
    main(sys.argv[1:])
