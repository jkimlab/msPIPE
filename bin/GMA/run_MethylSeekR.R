#!/usr/bin/env Rscript

library(BSgenome);
library(MethylSeekR);
library(parallel);
library("stringi");
set.seed(123);

args <- commandArgs(TRUE);
assembly_version <- args[1];
cpu <- args[2];
chr <- args[3];
file <- args[4];
out_dir <- args[5];

av_gen <- available.genomes(splitNameParts=TRUE);
index<-grep(assembly_version,av_gen$pkgname);
library(av_gen$pkgname[index[1]],character.only = TRUE);
organism_name <- get(stri_split_lines1(av_gen$organism[index[1]]));
av_gen$pkgname[index[1]];
sLengths=seqlengths(organism_name);
methFname <- file;
message(paste0("File: ",methFname));
message(paste0("chromosome: ",chr));
message(paste0("output directory: ",out_dir));
meth.gr <-readMethylome(FileName=methFname, seqLengths=sLengths);
PMDsegments.gr <- NA;

library(rtracklayer);
session <- browserSession();
genome(session) <- assembly_version;
query <- ucscTableQuery(session, "cpgIslandExt");
CpGislands.gr <- track(query);
genome(CpGislands.gr) <- NA;
CpGislands.gr <- suppressWarnings(resize(CpGislands.gr, 5000, fix="center"));
stats <- calculateFDRs(m=meth.gr, CGIs=CpGislands.gr, PMDs=PMDsegments.gr, num.cores=cpu);
FDR.cutoff <- 5;
m.sel <- 0.5;
n.sel=as.integer(names(stats$FDRs[as.character(m.sel), ] [stats$FDRs[as.character(m.sel), ]<FDR.cutoff])[1]);
UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, PMDs=PMDsegments.gr, num.cores=cpu, myGenomeSeq=organism_name, seqLengths=sLengths);
final_pdf <- paste0(out_dir,"/",chr,".UMRsLMRs.gr.pdf");
umrslmrs_gr <- paste0(out_dir,"/",chr,".UMRsLMRs.gr.rds");
umrslmrs_tab <- paste0(out_dir,"/",chr,".UMRsLMRs.tsv");
message("Print UMRs and LMRs\n");
message("\tumrslmrs_gr\n\tumrslmrs_tab\n");
saveUMRLMRSegments(segs=UMRLMRsegments.gr, GRangesFilename=umrslmrs_gr,  TableFilename=umrslmrs_tab);
dev.off();
