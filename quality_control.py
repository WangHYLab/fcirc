#!/usr/bin/env python3
# -*- coding: utf8 -*-

import os
import sys


def quality_filter_SE(quality_var,input_file,pattern="single-end"):
    dir,file=os.path.dirname(input_file),os.path.basename(input_file)
    out=file.split('.')[0]+"_temp"+".fastq"
    cutadapt='cutadapt -q %i,%i --trim-n  -o %s '%(quality_var,quality_var,out)
    qulity_cmd=os.system(cutadapt)
    if qulity_cmd!=0:
        print("Error: Fail to filter low quality reads")
        sys.exit(1)
    os.remove(input_file)
    os.rename(out,input_file)
    return 0


def quality_filter_PE(quality_var,R1,R2,pattern="pair-end"):
    dir, file1,file2 = os.path.dirname(R1), os.path.basename(R1),os.path.basename(R1)
    outR1,outR2=file1.split('.')[0]+'_temp1'+'.fastq',file2.split('.')[0]+'_temp1'+'.fastq'
    cutadapt = 'cutadapt -q %i,%i --trim-n  -o %s -p %s %s %s' % (quality_var, quality_var,outR1,outR2,R1,R2)
    qulity_cmd = os.system(cutadapt)
    if qulity_cmd != 0:
        print("Error: Fail to filter low quality reads")
        sys.exit(1)
    os.remove(R1)
    os.remove(R2)
    os.rename(outR1,R1)
    os.rename(outR2,R2)
    return 0




if __name__=='__main__':
    quality_filter_SE(sys.argv[1:])


