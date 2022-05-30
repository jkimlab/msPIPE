#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager",quietly = TRUE) )
    install.packages("BiocManager", repos='http://cran.us.r-project.org', quiet=TRUE);

BiocManager::install("BiocGenerics", quiet=TRUE);
BiocManager::install("S4Vectors", quiet=TRUE);
BiocManager::install("Rhtslib", quiet=TRUE);
BiocManager::install("IRanges", quiet=TRUE);
BiocManager::install("GenomeInfoDb", quiet=TRUE);
BiocManager::install("GenomicRanges", quiet=TRUE);
BiocManager::install("Biostrings", quiet=TRUE);
BiocManager::install("rtracklayer", quiet=TRUE);
BiocManager::install("BSgenome", quiet=TRUE);
BiocManager::install("methylKit", quiet=TRUE);
BiocManager::install("MethylSeekR", quiet=TRUE);
BiocManager::install("bsseq", quiet=TRUE)
BiocManager::install("BiocParallel", quiet=TRUE)

install.packages("circlize", repos='http://cran.us.r-project.org',quiet=TRUE);
install.packages("MASS", repos='http://cran.us.r-project.org',quiet=TRUE);
install.packages("stringi", repos='http://cran.us.r-project.org',quiet=TRUE);
install.packages("ggplot2", repos='http://cran.us.r-project.org',quiet=TRUE);
install.packages("gprofiler2", repos='http://cran.us.r-project.org',quiet=TRUE);
install.packages("doParallel", repos='http://cran.us.r-project.org',quiet=TRUE);
install.packages("foreach", repos='http://cran.us.r-project.org',quiet=TRUE);
