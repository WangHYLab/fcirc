# Fcirc
**Fcirc** is a pipeline for exploring linear transcripts and circRNAs of known fusions from RNA-Seq data. Known fusion genes are from the multiple databases (COSMIC, ChimerDB, TicDB, FARE-CAFE and FusionCancer) or user-added gene-pairs. It costs less time to find fusions with higher sensitivity than existing methods for detecting fusion. The steps of Fcirc are as follows:

![Fcirc pipeline](https://github.com/WangHYLab/supplementary_files/blob/master/Images/Figure_1.png "fcirc pipeline")

## Installation
**Fcirc** is written in **python3**, requiring [HISAT2](http://ccb.jhu.edu/software/hisat2/index.shtml) for aligning reads, [samtools](http://www.htslib.org/download/) for selecting reads and python packages numpy, scipy, pysam.
#### Hardware requirements
For running Fcirc a computer with the following configuration is needed:
* minimum 8 GB of RAM (aligning to human genome reference)
* 1 CPU (minimum)

#### Installing Fcirc from Github Clone
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
Make sure that **hisat2** and **samtools** are added to environment variables so that Fcirc can invoke them.

#### Preparing genome resource and known fusion-pairs
* The genome resource is hisat2 index, which can be downloaded from hisat2 websites[http://ccb.jhu.edu/software/hisat2/index.shtml]. For human fusion transcript detection, it's recommended to use genome_tran of GRCh38 or GRCh37. It can also be finished with FASTA sequence file and annotation GTF file by hisat2 script.

* Known fusion-pairs can be downloaded from [Github page](https://github.com/WangHYLab/fcirc) and bipartite fusions index can be built by hisat2-build in reference_fusion_info directory as follows:

```
cd reference_fusion_info
hisat2-build fusiongenes_ref_U.fa fusiongenes_ref_U
hisat2-build fusiongenes_ref_V.fa fusiongenes_ref_V
```

* (optional) Add gene pairs of fusions you are interested in.

## Usage
#### Input data
The input data shall be single-end or paired-end RNA-Seq in FASTQ format, which can be raw data or trimmed data.

#### Command line options
Fcirc can be run with a simple command line.
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
            fusion genes coordinates file (defalut: fusion_genes_coordinate.txt)    
        -1 <fastq1>, --file1 <fastq1>
            fastq file 1 (single-end pattern:only -1)
        -2 <fastq2>, --file2 <fastq2>
            fastq file 2 (paired-end pattern: -1 and -2, files should be like -1 xxx_1.fastq -2 xxx_2.fastq)

    Optional:     
        -o <output_dir>, --output <outout_dir>
            output file directory (default: .)
        -t <int>, --thread <int>
            number of hisat2 alignment and pysam filter threads to launch (default: 1)    

    Others:
        -h, --help
            help information  
        -v, --version
            version information 
```

#### Output data
The output includes: 
* fusion information
* f-circRNA information

**1. Fusion information** is stored in a file **'fusion_results.tsv'** as the following format:
```
#Fusion Name    5'Gene  3'Gene  5'Gene BreakPoint Pos   3'Gene BreakPoint Pos   5'Gene Breakpoint Seq   3'Gene Breakpoint Seq   5'and 3'Common Breakpoint Seq   BreakpointReads Count   BreakpointReads         BreakpointStrand Count(+,-)     ScanningReads Count     ScanningReads           ScanningStrand Count(+,-)       Fusion Seq Length   P-Value
PML-RARA	PML	RARA	28736	39134	CAGGGGAAAG	AGCCATTGAG	.	171	SRR3239817.18109433...	2,169	25	SRR3239817.15285364,SRR3239817.48232615...	5,20	38066	0.00439
...
...
```
The description of each column:
</br>**#Fusion Name** - - The name of the fusion
</br>**5'Gene** - - The gene encoding the 5' end of the fusion transcript
</br>**3'Gene** - - The gene encoding the 3' end of the fusion transcript
</br>**5'Gene BreakPoint Pos** - - The position of the breakpoint for the 5' end of the fusion transcript
</br>**3'Gene BreakPoint Pos** - - The position of the breakpoint for the 3' end of the fusion transcript
</br>**5'Gene Breakpoint Seq** - - Sequence of the 5'Gene at the fusion breakpoint 
</br>**3'Gene Breakpoint Seq** - - Sequence of the 3'Gene at the fusion breakpoint 
</br>**5'and 3'Common Breakpoint Seq** - - The same sequence at the breakpoint of the 3' end of the gene and the 5' end of the gene
</br>**BreakpointReads Count** - -The number of reads spanning the fusion breakpoint
</br>**BreakpointReads** - -The reads spanning the fusion breakpoint
</br>**BreakpointStrand Count(+,-)** - - The number of reads located in forward strand and reverse strand respectively
</br>**ScanningReads Count(+,-)** - - The number of pair of reads are located on both sides of the breakpoint
</br>**Fusion Seq Length** - - The length of the fusion transcript
</br>**P-Value** - - A p value indicating if reads around the breakpoint are evenly distributed


**2. FcircRNA information** is stored in a file **'fcircRNA_results.tsv'** as the following format:
```
#FcircRNA_NO	Fusion Name	Backsplice start	Backsplice end	FusionBreakPoint Pos	FusionSeq Length	Support FcircRNA Reads Count	Support FcircRNA Reads	                FcircRNA Strand Count(+,-)	FCI
No_1    PML-RARA	3400	32843	28736	38066	2	SRR3239817.46897792,SRR3239817.4206524	0,2	0.000
No_2	PML-RARA	28596	32347	28736	38066	1	SRR3239817.23906640	0,1	0.000
No_3	PML-RARA	3357	28856	28736	38066	2	SRR3239817.17076681,SRR3239817.2305597	0,2	0.000
...
...
```
The description of each column:
</br>**#FcircRNA_NO** - - The id of fusion circRNA
</br>**Fusion Name** - - The name of fusion gene
</br>**Backsplice start** - - The starting position of back-spliced end
</br>**Backsplice end** - - The end position of back-spliced end
</br>**FusionBreakPoint Pos** - - The position of fusion breakpoint
</br>**FusionSeq Length** - - The length of fusion sequence
</br>**Support FcircRNA Reads Count** - - The number of reads supporting the f-circRNA
</br>**Support FcircRNA Reads** - - The reads supporting the f-circRNA
</br>**FcircRNA Strand Count(+,-)** - - The number of reads supporting f-circRNA on positive and negative strand
</br>**FCI** - - The normalized expression of f-circRNA (FCI=(reads count * 10^6)/(sequencing depth * fusion gene length))

## Quick start
You can start this pipeline using a testing RNA-Seq data, whose reads are partially from a RNA-Seq dataset SRR3239817 (NCBI SRA database), an acute leukaemia cell line NB4.
```
python fcirc.py -t 4 -o fcirc_out -x transcriptome_HISAT2_index_path -f known_fusion_directory_path -1  test_fastq_path
```
It costs few minutes. If it runs successfully, some log information will be printed as following:
```
2019-11-11 22:37:46,453 fcirc.py[line:426][2019-11-11 22:37:46] Start running # fcirc.py -t 4 -o fcirc_out0 -x transcriptome_HISAT2_index_path -f known_fusion_directory_path -1 test_fastq_path
[2019-11-11 22:37:51] Finish mapping reads to transcription!
[2019-11-11 22:37:51] Finish mapping reads to fusion references U!
[2019-11-11 22:37:52] Finish mapping reads to fusion references V!
[2019-11-11 22:37:52] Finish dropping unmapped read in fusion references U and V!
Find 274 Reads in U! 274 Reads in V!
[2019-11-11 22:37:52] Finish filtering fusion-related reads in fusion references U and V!
[2019-11-11 22:37:53] Finish mapping reads to inferred fusion references!
Find 22 kind(s) of fcircRNAs!
[2019-11-11 22:37:53] Finish all!See the result in 'fcircRNA_results.tsv'!
```

