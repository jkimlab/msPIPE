#!/usr/bin/env Rscript

library(BSgenome);

args <- commandArgs(TRUE);
assembly_version <- args[1];

av_gen <- available.genomes(splitNameParts=TRUE);
index<-grep(assembly_version,av_gen$pkgname);
avail <- grep(assembly_version, av_gen$pkgname[index]);
if(length(avail)!=0){
		BiocManager::install(av_gen$pkgname[index]);
		print(1);
}else{
		print(-1);
		message(paste0(av_gen$pkgname[index]," is not available."));
}
