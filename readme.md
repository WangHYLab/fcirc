# Fcirc
**Fcirc** is a pipeline for exploring linear transcripts and circRNAs of known fusions based on RNA-Seq data. Known fusion genes are from the multiple databases (COSMIC, ChimerDB, TicDB, FARE-CAFE and FusionCancer) or user-added gene-pairs. It costs less time to find fusions with higher sensitivity than existing methods for detecting fusions. The steps of Fcirc are as follows:

![Fcirc pipeline](https://github.com/WangHYLab/supplementary_files/blob/master/Images/Figure_1small.png "fcirc pipeline")

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
pip install cutadapt
```
Make sure that **hisat2** and **samtools** are added to environment variables so that Fcirc can invoke them.

#### Preparing genome resource and known fusion-pairs
* The genome resource is hisat2 index, which can be downloaded from [hisat2 website](http://ccb.jhu.edu/software/hisat2/index.shtml). For human fusion transcript detection, it's recommended to use genome_tran of GRCh38 or GRCh37. It can also be finished with FASTA sequence file and annotation GTF file by hisat2 script.

* Known fusion-pairs can be downloaded from [Github page](https://github.com/WangHYLab/fcirc) and bipartite fusions index can be built by hisat2-build in reference_fusion_info directory as follows:

```
unzip fusion_total_index.zip
cd fusion_total_index
hisat2-build fusiongenes_ref_U.fa fusiongenes_ref_U
hisat2-build fusiongenes_ref_V.fa fusiongenes_ref_V
```

* (optional) Add gene pairs of fusions in **reference_fusion_info/fusion_table.tsv** with 3 columns (gene-pairs, 5'gene, 3'gene) you are interested in. For human, [genome FASTA file and annotation GTF file](http://asia.ensembl.org/info/data/ftp/index.html) can be downloaded here.
```
python3 build_graph.py --genome absolute__to_genome --gtf absolute__to_gtf --tab absolute_path_to_fusionpairs_table
# e.g. python3 build_graph.py --genome ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf ../ref/Homo_sapiens.GRCh38.105.gtf --tab reference_fusion_info/fusion_table.tsv  --outdir GRCh38_fusion_index  
```

##### <font color='red'>  Warning: gene pairs joined with'--' not '-' should be placed in the first column of fusionpairs_table </font>

## Usage
#### Input data
The input data shall be single-end or paired-end RNA-Seq in FASTQ format, which can be raw data or trimmed data.

#### Command line options
Fcirc can be run with a simple command line.
```
python fcirc.py [options] -x <ht2-trans-idx> -f <ht2-fusion-idx-dir> -c <fusion-genes-coordinates> {-1 <fastq1> | -1 <fastq1> -2 <fastq2>} 
```
Arguments can be used as following:
```
    Required:
        -x <ht2-trans-idx>, --trans_idx <ht2-trans-idx>
            transcription index filename prefix (minus trailing .X.ht2)
        -f <ht2-fusion-idx-dir>, --fusion_idx_dir <ht2-fusion-idx-dir>
            fusion index directory (contains fusiongenes_ref_U and fusiongenes_ref_V) 
        -1 <fastq1>, --file1 <fastq1>
            fastq file 1 (single-end pattern: only -1)
        -2 <fastq2>, --file2 <fastq2>
            fastq file 2 (paired-end pattern: -1 and -2, files should be like -1 xxx_1.fastq -2 xxx_2.fastq)

    Optional:
        -q <quality_val>
           the minimum phred qulaity of read(default:0)
        -c <fusion-genes-coordinates> --fusion_genes_coord
           fusion genes coordinates file (defalut: fusion_genes_coordinate.txt in fusion index directory)   
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
#Fusion_Name	5'Gene	3'Gene	5'Gene_chr	5'Gene_strand	3'Gene_chr	3'Gene_strand	5'Gene_BreakPoint_Pos	3'Gene_BreakPoint_Pos	5'and3'_Common_Breakpoint_Seq	BreakpointReads_Count	BreakpointReads	BreakpointStrand_Count(+,-)	ScanningReads_Count	ScanningReads	ScanningStrand_Count(+,-)	P-Value
PML--RARA	PML	RARA	15	+	17	+	74023408	40348313	.	117	SRR3239817.48728782,SRR3239817.46047306,SRR3239817.46553524,SRR3239817.16929141,SRR3239817.19547854,SRR3239817.24567755,SRR3239817,......
```
The description of each column:
</br>**#Fusion Name** - - The name of the fusion
</br>**5'Gene** - - The gene encoding the 5' end of the fusion transcript
</br>**3'Gene** - - The gene encoding the 3' end of the fusion transcript
</br>**5'Gene_chr**- - The chromosome of 5'end gene
</br>**5'Gene_strand**- - The strand of 5'end gene
</br>**3'Gene_chr** -- The chromosome of 3'end gene
</br>**3'Gene_strand**- - The strand of 3'end gene
</br>**5'Gene BreakPoint Pos** - - The position of the breakpoint for the 5' end of the fusion transcript
</br>**3'Gene BreakPoint Pos** - - The position of the breakpoint for the 3' end of the fusion transcript
</br>**5'and 3'Common Breakpoint Seq** - - The same sequence at the breakpoint of the 3' end of the gene and the 5' end of the gene
</br>**BreakpointReads Count** - -The number of reads spanning the fusion breakpoint
</br>**BreakpointReads** - -The reads spanning the fusion breakpoint
</br>**BreakpointStrand Count(+,-)** - - The number of reads located in forward strand and reverse strand respectively
</br>**ScanningReads Count(+,-)** - - The number of pair of reads are located on both sides of the breakpoint
</br>**ScanningReads**- - The reads located on both sides of the breakpoint
</br>**ScanningStrand_Count(+,-)** -- The number of Scanning reads located in forward strand and reverse strand respectively
</br>**P-Value** - - A p value indicating if reads around the breakpoint are evenly distributed


**2. FcircRNA information** is stored in a file **'fcircRNA_results.tsv'** as the following format:
```
#FcircRNA_NO	Fusion Name	Backsplice_start	Backsplice_end	Fusion5'_BreakPoint_Pos	Fusion3'_BreakPoint_Pos	Support_FcircRNA_Reads_Count	FcircRNA_Strand_Count(+, -)	Support_FcircRNA_Reads
No_1	PML--RARA	15:74023268:+	17:40351924:+	15:74023408:+	17:40348313:+	1	0,1	SRR3239817.23906640
No_2	PML--RARA	15:73998438:+	17:40352058:+	15:74023408:+	17:40348313:+	3	0,3	SRR3239817.6429653,SRR3239817.3386413,SRR3239817.31829112
No_3	PML--RARA	15:73998328:+	17:40354455:+	15:74023408:+	17:40348313:+	1	0,1	SRR3239817.3123010
No_4	PML--RARA	15:73998193:+	17:40352044:+	15:74023408:+	17:40348313:+	5	0,5	SRR3239817.29876711,SRR3239817.36732283,SRR3239817.47058005,SRR3239817.32495621,SRR3239817.13611951
No_5	PML--RARA	15:73998454:+	17:40352406:+	15:74023408:+	17:40348313:+	1	0,1	SRR3239817.28808693
No_6	PML--RARA	15:74022909:+	17:40355327:+	15:74023408:+	17:40348313:+	3	0,3	SRR3239817.42010789,SRR3239817.11495312,SRR3239817.33451057
......
......
...
```
The description of each column:
</br>**#FcircRNA_NO** - - The id of fusion circRNA
</br>**Fusion Name** - - The name of fusion gene
</br>**Backsplice start** - - The starting position of back-spliced end
</br>**Backsplice end** - - The end position of back-spliced end
</br>**Fusion5'_BreakPoint_Pos** - - The position of fusion breakpoint on 5'end
</br>**Fusion3'_BreakPoint_Pos** - - The position of fusion breakpoint on 3'end
</br>**Support_FcircRNA_Reads_Count** - - The number of reads supporting the f-circRNA
</br>**FcircRNA_Strand_Count(+, -)** - - The number of reads supporting f-circRNA on positive and negative strand
</br>**Support FcircRNA Reads** - - The reads supporting the f-circRNA

## Quick start
You can start this pipeline using a testing RNA-Seq data, whose reads are partially from a RNA-Seq dataset SRR3239817 (NCBI SRA database), for an acute leukaemia cell line NB4.
```
python fcirc.py -t 4 -o fcirc_out -x transcriptome_HISAT2_index_path -f known_fusion_directory_path -c fusion_genes_coordinate.txt -1 test_fastq_path
```
It costs few minutes. If it runs successfully, some log information will be printed as following:
```
[2022-01-13 21:48:19] Start running # fcirc/fcirc.py -t 1 -o fcirc_test -x ref/grch38_tran/genome_tran -f fcirc/fusion_total_index/ -1 fcirc/test_data/test.fastq.gz
[2022-01-13 21:48:27] Finish mapping reads to transcription!
[2022-01-13 21:48:27] Finish mapping reads to fusion references U!
[2022-01-13 21:48:27] Finish mapping reads to fusion references V!
[2022-01-13 21:48:27] Finish dropping unmapped read in fusion references U and V!
Find 274 Reads in U! 274 Reads in V!
[2022-01-13 21:48:27] Finish filtering fusion-related reads in fusion references U and V!
[2022-01-13 21:48:29] Finish mapping reads to inferred fusion references!
Find 22 kind(s) of fcircRNAs!
[2022-01-13 21:48:29] Finish all!See the result in 'fcircRNA_results.tsv'!
```

## Citation
Cai Z, Xue H, Xu Y, et al. Fcirc: A comprehensive pipeline for the exploration of fusion linear and circular RNAs. Gigascience. 2020;9(6):giaa054. doi: [10.1093/gigascience/giaa054](https://doi.org/10.1093/gigascience/giaa054)
