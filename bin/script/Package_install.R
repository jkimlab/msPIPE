#!/usr/bin/env Rscript
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager");

BiocManager::install("BiocGenerics");
BiocManager::install("S4Vectors");
BiocManager::install("Rhtslib");
BiocManager::install("IRanges");
BiocManager::install("GenomeInfoDb");
BiocManager::install("GenomicRanges");
BiocManager::install("Biostrings");
BiocManager::install("rtracklayer");
BiocManager::install("BSgenome");
BiocManager::install("methylKit");
BiocManager::install("MethylSeekR");
BiocManager::install("bsseq")
BiocManager::install("BiocParallel")

install.packages("circlize");
install.packages("MASS");
install.packages("stringi");
install.packages("ggplot2");
install.packages("gprofiler2");
install.packages("doParallel");
install.packages("foreach");
