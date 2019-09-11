#!/usr/bin/env python3
# -*- coding: utf8 -*-
import sys, getopt
import os
import time
import logging
from filter_reads import drop_mapped_single_tofastq, drop_unmapped_single_tosam, drop_mapped_paired_tofastq, drop_unmapped_paired_tosam, sam_single_tofastq, sam_paired_tofastq, count_reads
from getpairedreads import getpairedreads
from reconstruction import reconstruct
from transform import transform_sam
from write_results import write_fusion, write_fcirc

def version():
    print("fcirc (fusion circRNA identifier)    Version: 1.0.1")


def usage():
    usage_string = '''
Program:    fcirc (fusion circRNA identifier)
Version:    1.0.1(written by python3)
Contact:    HongZhang Xue <xuehzh95@foxmail.com>

Usage:      python fcirc.py [options] -x <ht2-trans-idx> -f <ht2-fusion-idx-dir> {-1 <fastq1> | -1 <fastq1> -2 <fastq2>}

Arguments:
    Required:
        -x <ht2-trans-idx>, --trans_idx <ht2-trans-idx>
            transcription index filename prefix (minus trailing .X.ht2)
        -f <ht2-fusion-idx-dir>, --fusion_idx_dir <ht2-fusion-idx-dir>
            fusion index directory (contains fusiongenes_ref_U and fusiongenes_ref_V)
        -c <fusion-genes-coordinates> --fusion_genes_coord
            fusion genes coordinates file(defalut: fusion_genes_coordinate.txt)
        -1 <fastq1>, --file1 <fastq1>
            fastq file 1 (single-end pattern:only -1)
        -2 <fastq2>, --file2 <fastq2>
            fastq file 2 (paired-end pattern:-1 and -2, files should be like -1 xxx_1.fastq -2 xxx_2.fastq)

    Optional:
        -o <output_dir>, --output <outout_dir>
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
    return os.path.abspath(path.replace('~', os.getenv('HOME')))


def str_current_time():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


def main(argvs):
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s')
    if len(argvs) == 0:
        usage()
        sys.exit()

    try:
        opts, args = getopt.getopt(argvs, "hv1:2:x:f:c:o:t:", ["help", "version", "file1", "file2" ,"trans_idx", "fusion_idx_dir","fusion_genes_coord", "output", "thread"])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    # default arguments:
    thread = '1'
    output_path = os.getcwd()

    # required arguments:
    fastq1_path = ''
    fastq2_path = ''
    trans_idx_path = ''
    fusion_idx_dir_path = ''
    software_dir=os.path.split(os.path.abspath(sys.argv[0]))[0]
    fusion_genes_coord=software_dir+'/'+'fusion_genes_coordinate.txt'
    #if not os.path.exists(Genes_coord):
        #print('ERROR in finding Genes_Coordinate.txt,Please check Genes_Coordinate.txt exit in fcirc.py directory')
        #sys.exit(1)


    for opt, arg in opts:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ('-v', '--version'):
            version()
            sys.exit()

        if opt in ('-1','--file1'):
            pattern = "single-end"
            abs_arg = abs_path(arg)
            if os.path.isfile(abs_arg):
                fastq1_path = abs_arg

        if opt in ('-2','--file2'):
            pattern = "paired-end"
            abs_arg = abs_path(arg)
            if os.path.isfile(abs_arg):
                fastq2_path = abs_arg

        if opt in ('-x','--trans_idx'):
            abs_arg = abs_path(arg)
            if os.path.isfile(abs_arg+'.1.ht2'):
                trans_idx_path = abs_arg
            elif trans_idx_path == '' :
                print("Error:Fail to locate hisat2 transcription index prefix,please check {path}".format(path=abs_arg))
                sys.exit(1)

        if opt in ('-f','--fusion_idx_dir'):
            abs_arg = abs_path(arg)
            if os.path.isdir(abs_arg):
                fusion_idx_dir_path = abs_arg.rstrip('/')
            elif fusion_idx_dir_path == '' :
                print("Error:Fail to locate hisat2 fusion index prefix,please check {path}".format(path=abs_arg))
                sys.exit(1)

        if opt in ('-c','--fusion_genes_coord'):
            abs_arg=abs_path(arg)
            if os.path.isfile(abs_arg):
                fusion_genes_coord=abs_arg
            else:
                print('Error to locate fusion genes coordinate file,Please check{path}'.format(path=abs_arg))
                sys.exit(1)


        if opt in ('-o','--output'):
            abs_arg = abs_path(arg)
            if not os.path.isdir(abs_arg):
                try:
                    os.mkdir(abs_arg)
                except:
                    print("Error:Fail to make dir at {path}".format(path=abs_arg))
                    sys.exit(1)
            output_path = abs_arg.rstrip('/')
            os.chdir(output_path)
        if opt in ('-t','--thread'):
            thread = arg

    try:
        os.mkdir("temp")
    except:
        print("temp dir exists, it will be rewrited.")

    # single-end pattern
    if pattern == "single-end":
        if fastq1_path == '':
            print("Error:Fail to locate fastq file,please check the path!")
            sys.exit(1)

        projectname = fastq1_path.split('/')[-1].rstrip('.fastq').rstrip('.fq')

        hisat2_trans_command = '''hisat2 -p {thread} \
            --mp 12,4 \
            --rdg 10,6 \
            --rfg 10,6 \
            --no-softclip \
            -x {hisat2_trans_idx} \
            -U {singlefastq} \
            -S {output_transsam} \
            2>{output_translog}'''.format(thread=str(thread),
                                          hisat2_trans_idx=trans_idx_path,
                                          singlefastq=fastq1_path,
                                          output_transsam=output_path + "/temp/trans_" + projectname + ".sam",
                                          output_translog=output_path + "/temp/trans_" + projectname + ".log")
        return_state = os.system(hisat2_trans_command)
        if return_state != 0:
            exit(return_state)

        print("[{now}] Finish mapping reads to transcription!".format(now=str_current_time()))
        unmapped_fastq_path = drop_mapped_single_tofastq("temp/trans_" + projectname + ".sam")

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
                                           output_fusionsam=output_path + "/temp/fusionU_" + projectname + ".sam",
                                           output_fusionlog=output_path + "/temp/fusionU_" + projectname + ".log")
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
                                           output_fusionsam=output_path + "/temp/fusionV_" + projectname + ".sam",
                                           output_fusionlog=output_path + "/temp/fusionV_" + projectname + ".log")
        return_state = os.system(hisat2_fusionV_command)
        if return_state != 0:
            exit(return_state)
        print("[{now}] Finish mapping reads to fusion references V!".format(now=str_current_time()))

        mapped_samU_path = drop_unmapped_single_tosam("temp/fusionU_" + projectname + ".sam", thread=str(thread))
        mapped_samV_path = drop_unmapped_single_tosam("temp/fusionV_" + projectname + ".sam", thread=str(thread))
        print("[{now}] Finish dropping unmapped read in fusion references U and V!".format(now=str_current_time()))

        mapped_filtered_samU_path, mapped_filtered_samV_path = getpairedreads(mapped_samU_path, mapped_samV_path, pattern, fusion_idx_dir_path)
        print("[{now}] Finish filtering fusion-related reads in fusion references U and V!".format(now=str_current_time()))

        return_state = reconstruct(mapped_filtered_samU_path, mapped_filtered_samV_path,fusion_genes_coord, fusion_idx_dir_path)
        if return_state != 0:
            print("[{now}] No fusion-related read is detected!".format(now=str_current_time()))
            exit(return_state)

        # reads both exists in U and V
        sam_single_tofastq(mapped_filtered_samV_path,"temp/inferred_fusion.fastq")

        hisat2_verify_fusion_command = '''hisat2 -p {thread} \
            --mp 6,2 \
            --rdg 5,2 \
            --rfg 5,2 \
            --sp 1,1 \
            --score-min L,0,-0.8 \
            -x temp/inferred_fusion.fa \
            -U {inferred_fusion_fastq} \
            -S {output_inferred_fusionsam} \
            2>{output_inferred_fusionlog}'''.format(thread=str(thread),
                                                    inferred_fusion_fastq="temp/inferred_fusion.fastq",
                                                    output_inferred_fusionsam="temp/inferred_fusion.sam",
                                                    output_inferred_fusionlog="temp/inferred_fusion.log")
        return_state = os.system(hisat2_verify_fusion_command)
        if return_state != 0:
            exit(return_state)
        print("[{now}] Finish mapping reads to inferred fusion references!".format(now=str_current_time()))

        # check the distribution near the fusion break point
        return_state = write_fusion("temp/inferred_fusion.sam")
        if return_state != 0:
            print("[{now}] No fusion gene is detected!".format(now=str_current_time()))
            exit(return_state)

        # filter and get the circ-related reads
        transform_sam("temp/inferred_fusion.sam")
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
        return_state  = write_fcirc("temp/inferred_fusion_circ.sam", readcount)
        if return_state != 0:
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
            --mp 12,4 \
            --rdg 10,6 \
            --rfg 10,6 \
            --no-softclip \
            -x {hisat2_trans_idx} \
            -1 {mate1fastq} \
            -2 {mate2fastq} \
            -S {output_transsam} \
            2>{output_translog}'''.format(thread=str(thread),
                                            hisat2_trans_idx=trans_idx_path,
                                            mate1fastq=fastq1_path,
                                            mate2fastq=fastq2_path,
                                            output_transsam=output_path + "/temp/trans_"  + projectname + ".sam",
                                            output_translog=output_path + "/temp/trans_"  + projectname + ".log")
        return_state = os.system(hisat2_trans_command)
        if return_state != 0:
            exit(return_state)
        print("[{now}] Finish mapping reads to transcription!".format(now=str_current_time()))

        unmapped_fastq1_path, unmapped_fastq2_path = drop_mapped_paired_tofastq("temp/trans_" + projectname + ".sam")

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
                                           output_fusionsam=output_path + "/temp/fusionU_" + projectname + ".sam",
                                           output_fusionlog=output_path + "/temp/fusionU_" + projectname + ".log")
        return_state = os.system(hisat2_fusionU_command)
        if return_state != 0:
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
                                            output_fusionsam=output_path + "/temp/fusionV_" + projectname + ".sam",
                                            output_fusionlog=output_path + "/temp/fusionV_" + projectname + ".log")
        return_state = os.system(hisat2_fusionV_command)
        if return_state != 0:
            exit(return_state)
        print("[{now}] Finish mapping reads to fusion references V!".format(now=str_current_time()))

        mapped_samU_path = drop_unmapped_paired_tosam("temp/fusionU_" + projectname + ".sam", thread=str(thread))
        mapped_samV_path = drop_unmapped_paired_tosam("temp/fusionV_" + projectname + ".sam", thread=str(thread))

        mapped_filtered_samU_path, mapped_filtered_samV_path = getpairedreads(mapped_samU_path, mapped_samV_path, pattern, fusion_idx_dir_path)
        print("[{now}] Finish filtering fusion-related reads in fusion references U and V!".format(now=str_current_time()))

        return_state= reconstruct(mapped_filtered_samU_path, mapped_filtered_samV_path,fusion_genes_coord,fusion_idx_dir_path)
        if return_state != 0:
            print("[{now}] No fusion-related read is detected!".format(now=str_current_time()))
            exit(return_state)

        # reads both exists in U and V
        sam_paired_tofastq(mapped_filtered_samV_path,"temp/inferred_fusion_1.fastq","temp/inferred_fusion_2.fastq")

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
        if return_state != 0:
            print("[{now}] No fusion gene is detected!".format(now=str_current_time()))
            exit(return_state)
        print("[{now}] Finish mapping reads to inferred fusion references!".format(now=str_current_time()))

        # check the distribution near the fusion break point
        return_state = write_fusion("temp/inferred_fusion.sam")
        if return_state != 0:
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
        return_state  = write_fcirc("temp/inferred_fusion_circ.sam", readcount)
        if return_state != 0:
            print("[{now}] No fusion circRNA is detected!".format(now=str_current_time()))
            exit(return_state)
        print("[{now}] Finish all!See the result in 'fcircRNA_results.tsv'!".format(now=str_current_time()))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s')
    logging.info("[{now}] Start running # {command}".format(now=str_current_time(), command=' '.join(sys.argv)))
    main(sys.argv[1:])
