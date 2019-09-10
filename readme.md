# Fcirc
*Fcirc* is a pipeline for transcripts and circRNAs of known fusions exploration. The sourcse of known fusion genes is from the multiple databases [COSMIC,ChimerDB,TicDB,FARE-CAFE and FusionCancer]and gene-pairs user added. It cost less time to find fusion-related(fusion forward splicing and back-splicing transripts) reads with higher sensitity than novel detecting fusion methods. The steps of *fcirc* are as folowing:

![Fcirc pipeline](https://github.com/WangHYLab/fcirc/raw/master/Fig 1. Fcirc pipeline.png  "fcirc pipeline")

## Installation
*Fcirc* is written by **python3**, requiring [hisat2](http://ccb.jhu.edu/software/hisat2/index.shtml) for aligning reads, [samtools](http://www.htslib.org/download/) for selecting reads and python packages numpy,scipy,pysam.
#### Hardware requirements
For running fcirc it is needed a computer with:
* minimum 8 GB of RAM(aligning to hunman genome reference)
* 1 CPU (minimum)

#### Installing fcirc from Github Clone
```
    git clone https://github.com/WangHYLab/fcirc
```
#### Installing required dependencies
* hisat2 [http://ccb.jhu.edu/software/hisat2/index.shtml](http://ccb.jhu.edu/software/hisat2/index.shtml)
* samtools [http://www.htslib.org/download/](http://www.htslib.org/download/)
* some python packages which be insatalled by pip as following:
```
pip install -r requirements.txt
```
or
```
pip install numpy
pip install scipy
pip install pysam
```
Make sure that **hisat2** and **samtools** are add to envionment variables so that *fcirc* can invoke them.

#### Preparing genome resource and known fusion-pairs
* Genome resource is hisat2 index, which can be download from [hisat2 websites](http://ccb.jhu.edu/software/hisat2/index.shtml). For human fusion transcript detection, it's recommanded to use *genome_tran* of GRCh38 or GRCh37. It also can be finished with FASTA sequence file and annotation GTF file by hisat2's script.

* Known fusion-pairs can be downloaded from Github page(https://cancer.sanger.ac.uk/cosmic/fusion) and bipartite fusions index can be build by hisat2-build as following:
```
python3 downloadfusion.py known_fusion_dir
cd known_fusion_dir
hisat2-build fusiongenes_ref_U.fa fusiongenes_ref_U
hiast2-build fusiongenes_ref_V.fa fusiongenes_ref_V
```

* (optional) add gene_pairs you interested
```
python3 addfusion.py known_fusion_dir fusion_list.txt fusion_seq.fa 
```
## Usage
#### Input data
The input data shall be single-end or paired -end FASTQ files which can be raw data as well as trimmed data.

#### Command line options
*Fcirc* can be run with a simple command line.
```
python fcric.py [options] -x <ht2-trans-idx> -f <ht2-fusion-idx-dir> {-1 <fastq1> | -1 <fastq1> -2 <fastq2>} 
```
Arguments can be used as following:
```
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
PML-RARA        PML     RARA    1183                    178                      CAGGGGAAAG              CCATTGAGAC             .                               199                     SRR3239817.48859220...  6,193                           1                       SRR3239817.49196217     1,0                             2394                0.41116390271
...
...
```
The columns' meaning is as following:  
**P-Value** - Under the assumption of fusion inferred by *fcirc*, the p value of length cut by breakpoint distribution. If it is smaller than 0.05, the inferred fusion is not correct.


**2. FcircRNA information** is stored in file **'fcircRNA_results.tsv'** as the format:
```
#FcircRNA_NO	Fusion Name	Backsplice start	Backsplice end	FusionBreakPoint Pos	FusionSeq Length	Support FcircRNA Reads Count	Support FcircRNA Reads	                FcircRNA Strand Count(+,-)	FCI
No_1        	PML-RARA	1014            	1635	        1183	                2394	            	1	                            SRR3239817.17199266	                    0,1	                        8.407767446938172e-12
No_2	        PML-RARA	155	                1301        	1183	                2394	            	2	                            SRR3239817.2305597,SRR3239817.17076681  0,2	                        1.6815534893876344e-11
No_3        	PML-RARA	198	                1725	        1183	                2394	            	2	                            SRR3239817.4206524,SRR3239817.46897792	0,2	                        1.6815534893876344e-11
No_4        	PML-RARA	493	                1298	        1183	                2394            	1	                            SRR3239817.1361558                  	0,1	                        8.407767446938172e-12
...
...
```
The columns' meaning is as following:  
**FCI** - the normalized expression of f-circRNA(reads count/sequencing depth/fusion gene length). 

## Quick start
SRR3239817 from NCBI, a NB4 sample（test.fastq has been reduced in test_data） sequenced with single-end, 75bp can be test here.
```
python fcirc.py -t 4 -o fcirc_out -x grch37_tran/genome_tran -f known_fusion_dir -1 test_data/test.fastq
```
It may cost about half an hour,if it runs succcessfully, some log information will be printed as following:
```
[2018-06-02 16:06:59] Start running # python fcirc.py -t 4 -o fcirc_out -x grch37_tran/genome_tran -f known_fusion_dir -1 SRR3239817.fastq
[2018-06-02 16:13:07] Finish mapping reads to transcription!
[2018-06-02 16:16:19] Finish mapping reads to fusion references U!
[2018-06-02 16:17:43] Finish mapping reads to fusion references V!
[2018-06-02 16:18:45] Finish dropping unmapped read in fusion references U and V!
Find 293 Reads in U! 293 Reads in V!
[2018-06-02 16:18:47] Finish filtering fusion-related reads in fusion references U and V!
PML-RARA 293
Pos: 1183 178
[2018-06-02 16:18:49] Finish mapping reads to inferred fusion references!
Find 35 kind(s) of fcircRNAs!
[2018-06-02 16:18:49] Finish all!See the result in 'fcircRNA_results.tsv'!
```

