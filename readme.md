# Fcirc
**Fcirc** is a pipeline for exploring transcripts and circRNAs of known fusions. The sources of known fusion genes are from the multiple databases (COSMIC,ChimerDB, TicDB, FARE-CAFE and FusionCancer) or user-added gene-pairs. It costs less time to find fusion-related (fusion forward splicing and back-splicing transcripts) reads with higher sensitivity than novel detecting fusion methods. The steps of **Fcirc** are as follows:

![Fcirc pipeline](https://github.com/WangHYLab/supplementary_files/blob/master/Images/Figure_1.png "fcirc pipeline")

## Installation
*Fcirc* is written in **python3**, requiring [HISAT2](http://ccb.jhu.edu/software/hisat2/index.shtml) for aligning reads, [samtools](http://www.htslib.org/download/) for selecting reads and python packages numpy,scipy,pysam.
#### Hardware requirements
For running fcirc a computer with the following is needed:
* minimum 8 GB of RAM(aligning to hunman genome reference)
* 1 CPU (minimum)

#### Installing fcirc from Github Clone
```
    git clone https://github.com/WangHYLab/fcirc
```
#### Installing required dependencies
* hisat2 [http://ccb.jhu.edu/software/hisat2/index.shtml](http://ccb.jhu.edu/software/hisat2/index.shtml)
* samtools [http://www.htslib.org/download/](http://www.htslib.org/download/)
* some python packages that can be installed by pip are as follows:
```
pip install -r requirements.txt
```
or
```
pip install numpy
pip install scipy
pip install pysam
```
Make sure that **hisat2** and **samtools** are added to environment variables so that **Fcirc** can invoke them.

#### Preparing genome resource and known fusion-pairs
* The genome resource is hisat2 index, which can be downloaded from hisat2 websites[http://ccb.jhu.edu/software/hisat2/index.shtml]. For human fusion transcript detection, it's recommended to use genome_tran of GRCh38 or GRCh37. It can also be finished with FASTA sequence file and annotation GTF file by hisat2 script.

* Known fusion-pairs can be downloaded from [Github page](https://github.com/WangHYLab/fcirc) and bipartite fusions index can be built by hisat2-build in reference_fusion_info directory as follows:

```
cd reference_fusion_info
hisat2-build fusiongenes_ref_U.fa fusiongenes_ref_U
hisat2-build fusiongenes_ref_V.fa fusiongenes_ref_V
```

* (optional) add gene pairs of fusions you interested

## Usage
#### Input data
The input data shall be single-end or paired-end RNA-Seq in FASTQ format, which can be raw data or trimmed data.

#### Command line options
**Fcirc** can be run with a simple command line.
```
python fcric.py [options] -x <ht2-trans-idx> -f <ht2-fusion-idx-dir> -c <fusion-genes-coordinates> {-1 <fastq1> | -1 <fastq1> -2 <fastq2>} 
```
Arguments can be used as following:
```
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
```

#### Output data
The output includes 
* fusion information
* f-circRNA information

**1. Fusion information** is stored in file **'fusion_results.tsv'** as the format:
```
#Fusion Name    5'Gene  3'Gene  5'Gene BreakPoint Pos   3'Gene BreakPoint Pos   5'Gene Breakpoint Seq   3'Gene Breakpoint Seq   5'and 3'Common Breakpoint Seq   BreakpointReads Count   BreakpointReads         BreakpointStrand Count(+,-)     ScanningReads Count     ScanningReads           ScanningStrand Count(+,-)       Fusion Seq Length   P-Value
PML-RARA        PML     RARA    28736                  39134                    CAGGGGAAAG              AGCCATTGAG             .                               171                    SRR3239817.18109433...  2,169                           25                       SRR3239817.15285364...     5,2                             38066     0.00438668413546
...
...
```
Description of each column's values:
**#Fusion Name** - -The name of the fusion transcipt,from which fusion genes can be infered.
**5'Gene** - -The gene located at the 5' end of the fusion transcript.
**3'Gene** - - The gene located at the 3' end of the fusion transcript.
**5'Gene BreakPoint Pos** - - The position of the break point of the 5' end of the fusion transcript
**3'Gene BreakPoint Pos** - - The position of the break point of the 3' end of the fusion transcript
**5'Gene Breakpoint Seq** - - Sequence of the 5' end of the gene in the fusion breakpoint attachment.
**3'Gene Breakpoint Seq** - - Sequence of the 3' end of the gene in the fusion breakpoint attachment
**5'and 3'Common Breakpoint Seq** - - The same sequence at the breakpoint of the 3' end of the transcript and the 5' end of the gene
**BreakpointReads Count** - -The number of read spanning the fusion breakpoint.
**BreakpointReads** - -The reads spanning the fusion breakpoint
**BreakpointStrand Count(+,-)** - - The number of reads located in forward Strand and reverse strand respectively
**ScanningReads Count(+,-)** - - The number of pair of reads are located on both sides of the breakpoint.
**Fusion Seq Length** - - The length of the fusion transcripts sequence.
**P-Value** - - Under the assumption of fusion inferred by **Fcirc**, the p value of length cut by breakpoint distribution. If it is smaller than 0.05, the inferred fusion is not correct.


**2. FcircRNA information** is stored in file **'fcircRNA_results.tsv'** as the format:
```
#FcircRNA_NO	Fusion Name	Backsplice start	Backsplice end	FusionBreakPoint Pos	FusionSeq Length	Support FcircRNA Reads Count	Support FcircRNA Reads	                FcircRNA Strand Count(+,-)	FCI
No_1        	PML-RARA	28736            	39134	        1183	                2394	            	1	                            SRR3239817.17199266	                    0,1	                        8.407767446938172e-12
No_2	        PML-RARA	155	                1301        	1183	                2394	            	2	                            SRR3239817.2305597,SRR3239817.17076681  0,2	                        1.6815534893876344e-11
No_3        	PML-RARA	198	                1725	        1183	                2394	            	2	                            SRR3239817.4206524,SRR3239817.46897792	0,2	                        1.6815534893876344e-11
No_4        	PML-RARA	493	                1298	        1183	                2394            	1	                            SRR3239817.1361558                  	0,1	                        8.407767446938172e-12
...
...
```
Description of each column's values:
**#FcircRNA_NO** - The id of fusion circRNA.
**Fusion Name** - The fusion gene name constituting the fusion circular RNA
**Backsplice start** - The starting position of back-spliced.
**Backsplice end** - The end position of back-spliced..
**FusionBreakPoint Pos** - The position of fusion breakpoint. 
**FusionSeq Length** - The length of fusion sequence.
**Support FcircRNA Reads Count** - The number of reads Supporting the f-circRNA.
**Support FcircRNA Reads** - The reads Supporting the f-circRNA.
**FcircRNA Strand Count(+,-)** - The number of reads supporting f-circRNA on positive and negative strand.
**FCI** - The normalized expression of f-circRNA(reads count/sequencing depth/fusion gene length). 

## Quick start
SRR3239817 from NCBI SRA database, a NB4 sample（test.fastq has been reduced in test_data） sequenced with single-end, 75 bp can be tested here.
```
python fcirc.py -t 4 -o fcirc_out -x grch37_tran/genome_tran -f known_fusion_dir -1 test_data/test.fastq
```
It costs about half an hour. If it runs successfully, some log information will be printed as following:
```
[2018-06-02 16:06:59] Start running # python fcirc.py -t 4 -o fcirc_out -x grch37_tran/genome_tran -f known_fusion_dir -1 SRR3239817.fastq
[2018-06-02 16:13:07] Finish mapping reads to transcription!
[2018-06-02 16:16:19] Finish mapping reads to fusion references U!
[2018-06-02 16:17:43] Finish mapping reads to fusion references V!
[2018-06-02 16:18:45] Finish dropping unmapped read in fusion references U and V!
Find 274 Reads in U! 274 Reads in V!
[2018-06-02 16:18:47] Finish filtering fusion-related reads in fusion references U and V!
[2018-06-02 16:18:49] Finish mapping reads to inferred fusion references!
Find 22 kind(s) of fcircRNAs!
[2018-06-02 16:18:49] Finish all!See the result in 'fcircRNA_results.tsv'!
```

