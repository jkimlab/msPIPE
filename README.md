# msPIPE
- Methylation analysis pipeline for WGBS data


## Requirements

- Trim Galore ([https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
- Samtools ([http://www.htslib.org/](http://www.htslib.org/))
- Bismark ([https://github.com/FelixKrueger/Bismark](https://github.com/FelixKrueger/Bismark))
- cutadapt ([https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/))
<br />

*Or you can use msPIPE on docker without having to prepare the environment.* 

:point_right: [HOW TO USE msPIPE on docker](#using-docker) ![Docker](https://img.shields.io/badge/Docker-%230db7ed.svg?&logo=Docker&logoColor=white)
<br />
<br />


## Download

```
git clone https://github.com/jkimlab/msPIPE.git
```


## Running

### Preparing an input parameter file

The parameter file contains the information necessary for pipeline execution.

```
### INPUT PARAMETER FILE FORMAT ###

## [DMR]
## ANALYSIS1 = Two sample names for DMR analysis
## ANALYSIS2 = Two sample names for DMR analysis

## [REFERENCE]
## UCSC_NAME = UCSC reference version name
## FASTA = [path to reference fasta file(not required)]
## GTF= [path to reference gtf file(not required)]

## [LIB1]
## SAMPLE_NAME = sample name
## LIB_NAME = library name
## LIB_TYPE = P or S (Paired-End or Single-Read)
## FILE_1 = path to sequencing read file
## FILE_2 = path to sequencing read file
```


### Running pipeline

- Running command 
    ```
    /PATH/TO/msPIPE/msPIPE.py -p params.conf -o OUTDIR &> logs
    ```
    
- msPIPE options

    ```
    ./msPIPE/msPIPE.py -h

    usage: msPIPE.py [-h] --param params.conf --out PATH [--core int]
                     [--qvalue float] [--skip_trimming] [--skip_mapping]
                     [--skip_calling] [--skip_HMR] [--skip_DMR]

    optional arguments:
      -h, --help            show this help message and exit
      --param params.conf, -p params.conf
                            config format parameter file
      --out PATH, -o PATH   output directory
      --core int, -c int    core (default:5)
      --qvalue float, -q float
                            q-value cutoff (default:0.5)
      --skip_trimming       skip the trimgalore trimming
      --skip_mapping        skip the bismark mapping
      --skip_calling        skip the methylation calling
      --skip_HMR            skip the HMR analysis
      --skip_DMR            skip the DMR analysis
    ```

<br />


## Running example

- Running example using mouse rod WGBS data from [Corso-Díaz, Ximena et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7228806/)
- GEO accession : [GSE134873](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/bioproject/556668)
- Used samples
    - 24M [rep1](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361857), [rep2](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361856), [rep3](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361855)  
    - 3M [rep1](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361836), [rep2](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361837), [rep3](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/biosample/12361839)

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





