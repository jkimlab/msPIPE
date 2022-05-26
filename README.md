# msPIPE
- Methylation analysis pipeline for WGBS data


## Requirements

- [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [Samtools](http://www.htslib.org/)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml), [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Bismark](https://github.com/FelixKrueger/Bismark)
- [BS-Seeker2](https://github.com/BSSeeker/BSseeker2)
- Other required R packages can be installed via msPIPE/bin/script/Package_install.R.
<br />

*Or you can use msPIPE on docker without having to prepare the environment.* ***\< Recommended\>***  
:point_right: [HOW TO USE msPIPE on docker](#using-docker) ![Docker](https://img.shields.io/badge/Docker-%230db7ed.svg?&logo=Docker&logoColor=white)
<br />
<br />


## Download

```
git clone https://github.com/jkimlab/msPIPE.git
```



## Running
#### Running command
```
/PATH/TO/msPIPE/msPIPE.py -p params.conf -o OUTDIR 
```


#### Preparing an input parameter file
The parameter file must contain information necessary for pipeline execution.
* params_format.conf   

	> ###INPUT PARAMETER FILE FORMAT###
	>
	> [DMR]   
	> ANALYSIS1 = sample1, sample2 (Two sample names for DMR analysis)   
	> ANALYSIS2 = sample1, sample3 (Two sample names for DMR analysis)   
	>
	> [REFERENCE]  
	> UCSC_NAME = UCSC reference version name   
	> FASTA = [path to reference fasta file (not required)]  
	> GTF= [path to reference gtf file (not required)]  
	>
	> [LIB1]  
	> SAMPLE_NAME = sample name  
	> LIB_NAME = library name  
	> LIB_TYPE = P or S (Paired-End or Single-Read)  
	> FILE_1 = [path to sequencing read file]  
	> FILE_2 = [path to sequencing read file]  



### Additional options  
#### msPIPE options  
 ```
 usage: msPIPE.py [-h] --param params.conf --out PATH [--core int] [--qvalue float] [--skip_trimming] [--program bismark or bs2]
                            [--bsmooth] [--skip_mapping] [--skip_calling] [--calling_data PATH] [--skip_GMA] [--skip_DMR]

optional arguments:
  -h, --help            show this help message and exit
  --param params.conf, -p params.conf
                        config format parameter file
  --out PATH, -o PATH   output directory
  --core int, -c int    core (default:5)
  --qvalue float, -q float
                        q-value cutoff (default:0.5)
  --skip_trimming       skip the trimgalore trimming
  --program bismark or bs2
                        program option for mapping & calling
  --bsmooth             use bsmooth for DMR analysis
  --skip_mapping        skip the bismark mapping
  --skip_calling        skip the methylation calling
  --calling_data PATH, -m PATH
                        methylCALL directory
  --skip_GMA            skip the Gene-Methyl analysis
  --skip_DMR            skip the DMR analysis

 ```
 

#### Software options
You can choose the software to be used for analysis by selecting the module.
1. **<span style="color:DarkslateBlue">WGBS reads mapping & calling (--program bismark or bs2)<span>**.  
	Use this option to select the tools to be used for mapping and calling of wgbs reads.  
    --program bismark : use Bismark with bowtie2 (default)  
    --program bs2 : use BS-Seeker2 with bowtie   
    
2. **<span style="color:DarkslateBlue">Differential methylation analysis (--bsmooth)<span>**  
	With this option, BSmooth is used to analyze instead of methylKit and inhouse script.


#### Skip options
You can leave out some pipeline steps with the *--skip_\<STEP\>* option.  
The main steps of the entire pipeline and the steps that can be omitted are as follows.
1. **<span style="color:DarkslateGray">check all input<span>**. 
2. **<span style="color:DarkslateGray">Prepare bisulfite-converted reference genome (Bismark or BS-Seeker2)<span>**  
	- It will be skipped if the same assembly name of the bisulfite genome has already been created under msPIPE/reference/ directory.
3. **<span style="color:DarkslateGray">WGBS reads trimming  (TrimGalore)<span>**
	- Can drop with *--skip_trimming* option.
	- Trimmed reads to be used in mapping can be delivered through the TRIMMED_FILE_* parameters. ([LIB1] on below format)
	- Without TRIMMED_FILE_* parameters, the pipeline searches the files on the output directory.
4. **<span style="color:DarkslateGray">WGBS reads mapping  (Bismark or BS-Seeker2)<span>**
	- Can drop with *--skip_mapping* option.
	- Mapping file to be used in the next step can be delivered through the BAM_FILE parameter. ([LIB2] on below format)
	- Without BAM_FILE parameter, the pipeline searches the file on the output directory.
5. **<span style="color:DarkslateGray">Methylation calling (Bismark or BS-Seeker2)<span>**
	- Can drop with *--skip_calling* option.
	- Pipeline use calling output on the output directory.
	- Other msPIPE calling output can be given with the *--calling_data* option.
6. **<span style="color:DarkslateGray">Gene-Methylation analysis (Methylation profiling and Hypomethylated region analysis)<span>**
	- Can drop with *--skip_GMA* option.
7. **<span style="color:DarkslateGray">Differential methylation analysis<span>**
	- Can drop with *--skip_DMR* option.  

* params_format.conf

  > ###INPUT WITH SKIP OPTIONS###  
  > 
  > **...**
  >   
  > 
  > [LIB1]  
  > SAMPLE_NAME = sample name  
  > LIB_NAME = library name  
  > LIB_TYPE = P or S (Paired-End or Single-Read)  
  > TRIMMED_FILE_1 = [path to preprocessed read file (with --skip_trimming option)]  
  > TRIMMED_FILE_2 = [path to preprocessed read file (with --skip_trimming option)]  
  >
  > [LIB2]  
  > SAMPLE_NAME = sample name  
  > LIB_NAME = library name  
  > LIB_TYPE = P or S (Paired-End or Single-Read)  
  > BAM_FILE = [path to bismark mapping file (with --skip_mapping option)]  


<br /> 
  
___

## Running example

- Running example using mouse rod WGBS data from [Corso-Díaz, Ximena et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7228806/)  

    | GEO accession | sample-24M |  sample-3M |
    | ------------| ------------- | ------------- |
    | [GSE134873](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/bioproject/556668) | [rep1](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361857), [rep2](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361856), [rep3](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361855) | [rep1](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361836), [rep2](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361837), [rep3](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361839)|

- params_mouse.conf  
    *Replace the '/PATH/TO/DATA' with a data path on your local server.*
    
    ```
    [DMR]
    ANALYSIS1 = 24M, 3M

    [REFERENCE]
    UCSC_NAME = mm10

    [LIB1]
    SAMPLE_NAME = 24M
    LIB_NAME = 24M_rep1
    LIB_TYPE = P
    FILE_1 = /PATH/TO/DATA/SRX6589858_1.fastq.gz
    FILE_2 = /PATH/TO/DATA/SRX6589858_2.fastq.gz

    [LIB2]
    SAMPLE_NAME = 24M
    LIB_NAME = 24M_rep2
    LIB_TYPE = P
    FILE_1 = /PATH/TO/DATA/SRX6589859_1.fastq.gz
    FILE_2 = /PATH/TO/DATA/SRX6589859_2.fastq.gz

    [LIB3]
    SAMPLE_NAME = 24M
    LIB_NAME = 24M_rep3
    LIB_TYPE = P
    FILE_1 = /PATH/TO/DATA/SRX6589860_1.fastq.gz
    FILE_2 = /PATH/TO/DATA/SRX6589860_2.fastq.gz

    [LIB4]
    SAMPLE_NAME = 3M
    LIB_NAME = 3M_rep1
    LIB_TYPE = P
    FILE_1 = /PATH/TO/DATA/SRX6589850_1.fastq.gz
    FILE_2 = /PATH/TO/DATA/SRX6589850_2.fastq.gz
    
    [LIB5]
    SAMPLE_NAME = 3M
    LIB_NAME = 3M_rep2
    LIB_TYPE = P
    FILE_1 = /PATH/TO/DATA/SRX6589851_1.fastq.gz
    FILE_2 = /PATH/TO/DATA/SRX6589851_2.fastq.gz

    [LIB6]
    SAMPLE_NAME = 3M
    LIB_NAME = 3M_rep3
    LIB_TYPE = P
    FILE_1 = /PATH/TO/DATA/SRX6589852_1.fastq.gz
    FILE_2 = /PATH/TO/DATA/SRX6589852_2.fastq.gz
    ```
- Running command

    ```
     ./msPIPE/msPIPE.py -p params_mouse.conf -o mouse_result -c 5 -q 0.5
    ```

___

## Analysis Output

All output created by msPIPE will be written to the `methylCALL` and `Analysis` directories in the given `output` directory.
The output of pre-processing (read files processed by trimming and quality control), alignment, and methylation calling for each input library (named with LIB_NAME in config file) will be in `methylCALL` directory.
The output of methylation analysis will be in `Analysis` directory. 

* An example output structure of `Analysis` directory

	```
	[Analysis]
	|- avg_methlevel.pdf
	|- [annotations]
	|- [sample1]
	   |- ALL_TEXTFILES_AND_PLOTS_FOR_SAMPLE1
	   |- [AroundTSS]
	   |- [MethylSeekR]
	|- [sample2]
	   |- ...
	|- [DMR]
	   |- [sample1.sample2]
	```

* Output files and directories in `Analysis`

	* avg_methlevel.pdf : a bar plot of average methylation level for CpG, CHG, and CHH context
	* `annotations` : a directory with information of genes, exons, introns, promoters, and intergenic regions in BED format files
	* `sample1` : a directory with all results of methylation analysis for *sample1*
		* Average_methyl_lv.txt : average methylation level for each gene and its promoter
		* Avg_Genomic_Context_CpG.txt : average methylation level for each genomic context (gene, exon, intron, promoter, and intergenic)
		* CXX_methylCalls.bed : all methylation calls for each CX context (CXX is one of CpG, CHG, and CHH)
		* `AroundTSS`/meth_lv_3M.txt : for each gene, average methylation levels in bins around TSS (+/- 1500 bp)
		* `MethylSeekR` : a directory with all results for running MethylSeekR
		* UMR-Promoter.cnt.bed : the number of UMRs in each promoter region
		* UMR-Promoter.pos.bed : the genomic coordinates of UMRs in each promoter region
		* Circos.CpG_UMRs_LMRs.pdf : a circos plot for methylation level in whole-genome scale
		* Genomic_Context_CpG.pdf : a bar plot for average methylation level of each genomic context (gene, exon, intron, promoter, and intergenic)
		* hist_sample1_CXX.pdf : the distribution of methylation in CX context (CXX is one of CpG, CHG, and CHH)

If DMC/DMR analysis is performed, `DMR` directory will be created in `Analysis` directory. When the DMC/DMR analysis is performed, GO enrichment test will be carried out for the gene set with DMCs or DMRs in their promoter. 

* Examples of output files and directories in `DMR` for comparison pair sample1 and sample2

	* `sample1.sample2` : a directory with all results of DMC/DMR analysis, in this case sample1 will be treated as *control* and sample2 will be treated as *case*
		* DMR_q0.5.bed : information of differentially methylated regions
		* `methylkit` : output of running methylKit
		* DMC_q0.5.bed : filtered DMCs with q-value 0.5
		* hypoDMC_detailed_count_methyl.txt : the number of hypomethylated DMCs in each promoter (methylation level *case* < *control*)
		* hyperDMC_detailed_count_methyl.txt : the number of hypermethylated DMCs in each promoter (methylation level *case* > *control*)
		* intersection.DMC2Promoter.txt : a list of intersection between genes and DMCs
		* DMC_genelist.txt : a list of genes with DMCs overlapped their promoter region
		* DMC_gene.GOresult.txt : a text output of GO enrichment test for genes with DMCs from methylKit using g:Profiler
		* DMC_gene.GOresult.pdf : Plots of GO enrichment test for genes with DMCs from methylKit using g:Profiler 
		* DMR_gene.GOresult.txt : a text output of GO enrichment test for genes with DMRs from BSmooth using g:Profiler
		* DMR_gene.GOresult.pdf : Plots of GO enrichment test for genes with DMRs from BSmooth using g:Profiler 
___

## Using Docker

### build msPIPE docker image

```
git clone https://github.com/jkimlab/msPIPE.git
cd msPIPE
docker build -t jkimlab/mspipe:latest .
```
    
- or you can pull docker image from the docker hub
    ```
    docker pull jkimlab/mspipe:latest
    ```

### Preparing an input parameter file for Docker
 - The parameter file must be written based on the internal path of the docker container and placed within the output dir.
 - Mount the volumes with '-v' options to deliver input data and receive output results.
 - params_docker.conf  
    ```
    [DMR]
    ANALYSIS1 = 24M, 3M

    [REFERENCE]
    UCSC_NAME = mm10

    [LIB1]
    SAMPLE_NAME = 24M
    LIB_NAME = 24M_rep1
    LIB_TYPE = P
    FILE_1 = /msPIPE/data/SRX6589858_1.fastq.gz
    FILE_2 = /msPIPE/data/SRX6589858_2.fastq.gz

    [LIB2]
    SAMPLE_NAME = 24M
    LIB_NAME = 24M_rep2
    LIB_TYPE = P
    FILE_1 = /msPIPE/data/SRX6589859_1.fastq.gz
    FILE_2 = /msPIPE/data/SRX6589859_2.fastq.gz
    
    ...
    
    
    ```

### Running pipeline on Docker
    
 ```
 #docker run -v [local path]:[docker path] [docker image name] [msPIPE command]

 docker run -v /PATH/TO/INPUT/DATA:/msPIPE/data:ro -v /PATH/TO/REUSABLE/REFERENCE:/msPIPE/reference -v /PATH/TO/OUTDIR:/work_dir/ jkimlab/mspipe:latest msPIPE.py -p params_docker.conf -o result
 ```
 
 - Mount the volumes with '-v' options to deliver input data and receive output results.
    - input data dir → /msPIPE/data
    - reusable references dir → /msPIPE/reference
    - output dir → /work_dir
 - All local paths to mount volumes are must be expressed as absolute paths.
 - Replace the '/PATH/TO/*' with a directory path on your local server.
 
 
 
 <br />

    
 ## CONTACT

[bioinfolabkr@gmail.com](mailto:bioinfolabkr@gmail.com)
