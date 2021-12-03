# msPIPE
- Methylation analysis pipeline for WGBS data

## Download

    ```
    git clone https://github.com/jkimlab/msPIPE.git
    ```

## Download and installation with docker
- Build msPIPE docker image

    ```
    git clone https://github.com/jkimlab/msPIPE.git
    cd msPIPE
    docker build -t [image_name] .
    ```

- or you can pull docker image from the docker hub

    ```
    docker pull jkimlab/mspipe:latest
    ```


## Prepare input parameter file

```basic
### INPUT PARAMETER FORMAT ###

## [DMR]
## ANAYSIS1 = Two sample names for DMR analysis

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

- parameter_example.conf

    ```
    [DMR]
    ANALYSIS1 = SAMPLE1,SAMPLE2
    ANALYSIS2 = SAMPLE1,SAMPLE3

    [REFERENCE]
    UCSC_NAME = hg38

    [LIB1]
    SAMPLE_NAME = SAMPLE1
    LIB_NAME = SAMPLE1_lib1
    LIB_TYPE = P
    FILE_1 = /PATH/TO/DATA/SAMPLE1_lib1-1.fq.gz
    FILE_2 = /PATH/TO/DATA/SAMPLE1_lib1-2.fq.gz

    [LIB2]
    SAMPLE_NAME = SAMPLE1
    LIB_NAME = SAMPLE1_lib1
    LIB_TYPE = P
    FILE_1 = /PATH/TO/DATA/SAMPLE1_lib2-1.fq.gz
    FILE_2 = /PATH/TO/DATA/SAMPLE1_lib2-2.fq.gz

    [LIB3]
    SAMPLE_NAME = SAMPLE2
    LIB_NAME = SAMPLE2
        LIB_TYPE = P
    FILE_1 = /PATH/TO/DATA/SAMPLE2-1.fq.gz
    FILE_2 = /PATH/TO/DATA/SAMPLE2-2.fq.gz

    [LIB4]
    SAMPLE_NAME = SAMPLE3
    LIB_NAME = SAMPLE3
    LIB_TYPE = S
    FILE_1 = /PATH/TO/DATA/SAMPLE3.fq.gz

    ```


## Running
- Running command

    ```
    /PATH/TO/msPIPE/msPIPE.py -p parameter.conf -o OUTDIR &> logs
    ```


- msPIPE running options

    ```
    ./msPIPE.py

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
 

## **Third party tools**


- Trim Galore ([https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
- Samtools ([http://www.htslib.org/](http://www.htslib.org/))
- Bismark ([https://github.com/FelixKrueger/Bismark](https://github.com/FelixKrueger/Bismark))
- cutadapt ([https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/))



## Contact

[bioinfolabkr@gmail.com](mailto:bioinfolabkr@gmail.com)





