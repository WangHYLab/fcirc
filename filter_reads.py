#!/usr/bin/env python3
#import pysam
import os


def drop_mapped_single_tofastq(samfile_path):
    outfastq_path = "temp/unmapped_" + samfile_path.split('/')[-1].replace('sam', 'fastq')
    #status = pysam.fastq('-f', '4', samfile_path,'-0',outfastq_path, catch_stdout=False)
    os.system("samtools fastq -f 4 {samfile}>{outfastq} 2>/dev/null".format(samfile=samfile_path,
                                                                  outfastq=outfastq_path))
    return outfastq_path


def drop_unmapped_single_tosam(samfile_path,thread='1'):
    outsam_path = "temp/mapped_" + samfile_path.split('/')[-1]
    #status = pysam.view('-@', thread, '-h', '-F', '4','-o', outsam_path, samfile_path, catch_stdout=False)
    os.system("samtools view -@ {thread} -h -F 4 -o {outsam} {samfile} 2>/dev/null".format(thread=thread,
                                                                              outsam=outsam_path,
                                                                              samfile=samfile_path))
    return outsam_path


def drop_mapped_paired_tofastq(samfile_path):
    outfastq_paths = "temp/unmapped_" + samfile_path.split('/')[-1].replace('.sam', '_1.fastq'), "temp/unmapped_" + samfile_path.split('/')[-1].replace('.sam', '_2.fastq')
    #status = pysam.fastq('-F', '2', samfile_path ,'-1', outfastq_paths[0], '-2', outfastq_paths[1])
    os.system("samtools fastq -F 2 {samfile} -1 {outfastq_1} -2 {outfastq_2} 2>/dev/null".format(samfile=samfile_path,
                                                                                    outfastq_1=outfastq_paths[0],
                                                                                    outfastq_2=outfastq_paths[1]))
    return outfastq_paths


def drop_unmapped_paired_tosam(samfile_path,thread='1'):
    '''
    only for paired-end, select the reads mapped at least once
    '''
    #samfile_dir = '/'.join(samfile_path.split('/')[:-1])
    #print(samfile_path,samfile_dir)
    # one end mapped
    #status = pysam.view('-@', thread, '-h', '-F', '4', '-f', '8', '-o', "temp/temp1.sam", samfile_path, catch_stdout=False)
    #status = pysam.view('-@', thread, '-h', '-F', '8', '-f', '4', '-o', "temp/temp2.sam", samfile_path, catch_stdout=False)
    os.system("samtools view -@ {thread} -h -F 4 -f 8 -o temp/temp1.sam {samfile} 2>/dev/null".format(thread=thread,
                                                                                         samfile=samfile_path))
    os.system("samtools view -@ {thread} -h -F 8 -f 4 -o temp/temp2.sam {samfile} 2>/dev/null".format(thread=thread,
                                                                                         samfile=samfile_path))
    # two ends mapped
    #status = pysam.view('-@', thread, '-h', '-F', '12','-o', "temp/temp3.sam", samfile_path, catch_stdout=False)
    os.system("samtools view -@ {thread} -h -F 12 -o temp/temp3.sam {samfile} 2>/dev/null".format(thread=thread,
                                                                                         samfile=samfile_path))

    outsam_path = "temp/mapped_" + samfile_path.split('/')[-1]
    # merge mapped files
    #os.system("cat temp/temp1.sam temp/temp2.sam temp/temp3.sam>{outsam}".format(outsam=outsam_path))
    #status = pysam.merge('-@', thread, '-n', '-f', '-u', outsam_path, 'temp/temp1.sam', 'temp/temp2.sam', 'temp/temp3.sam')
    os.system("samtools merge -@ {thread} -n -f -u {outsam} temp/temp1.sam temp/temp2.sam temp/temp3.sam 2>/dev/null".format(thread=thread,outsam="temp/mapped_all.sam"))
    os.system("samtools sort -@ %s -n -O sam -o %s temp/mapped_all.sam"%(thread,outsam_path))
    # delete temp files
    os.system("rm -rf temp/mapped_all.sam")
    os.system("rm -rf temp/temp1.sam")
    os.system("rm -rf temp/temp2.sam")
    os.system("rm -rf temp/temp3.sam")

    return outsam_path


def sam_single_tofastq(samfile_path, outfastq_path):
    #with open(outfastq_path,'w') as f:
    #    f.write(pysam.fastq(samfile_path))
    os.system("samtools fastq {samfile}>{outfastq} 2>/dev/null".format(samfile=samfile_path,
                                                             outfastq=outfastq_path))

def sam_paired_tofastq(samfile_path, outfastq_path1, outfastq_path2):
    #pysam.fastq(samfile_path, '-1', outfastq_path1, '-2', outfastq_path2, catch_stdout=False)
    os.system("samtools fastq {samfile} -1 {outfastq_1} -2 {outfastq_2} 2>/dev/null".format(samfile=samfile_path,
                                                                                    outfastq_1=outfastq_path1,
                                                                                    outfastq_2=outfastq_path2))
def count_reads(samfile_path):
    logfile = samfile_path.replace(".sam",".log")
    with open(logfile, 'r') as f:
        count = f.readline().split(' ')[0]
    return int(count)

if __name__ == '__main__':
    pass
