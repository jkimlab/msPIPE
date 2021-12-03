#!/usr/bin/env Rscript
#install.packages("circlize")
library("circlize");
library("MASS");
##### arguments
args <- commandArgs(TRUE);
cytoband <- args[1];
path1 <- args[2];
out_pdf <- args[3];


umrs <- paste(path1, "sort.UMRs.bed", sep="")
lmrs <- paste(path1, "sort.LMRs.bed", sep = "")


CpGF <- paste(path1, "CpG_methylLev_window.bed", sep='') 
message("loading data...\n");
message("\tCpG context...");
CpG <- read.table(CpGF, sep="\t",quote = "\"", as.is=TRUE, header=FALSE, skipNul=TRUE)
names(CpG) <- c("chr", "start", "end", "value")
CpG$value <- CpG$value/100
#CpG <- CpG[!is.na(CpG$value),]


message("\tUMRs...");
UMRs <- read.table(umrs, header=F);
message("done\n\tLMRs...");
LMRs <- read.table(lmrs, header=F);




message("done\ndrawing pdf..\n");
pdf(out_pdf, width=20, height=20);
circos.clear();

message("\tcytoband..\n");
cytoband.file = cytoband;
cytoband.df = read.table(cytoband.file, colClasses = c("character", "numeric", "numeric", "character", "character"), sep = "\t");
circos.initializeWithIdeogram(cytoband.df);
#circos.initializeWithIdeogram(cytoband.df, chromosome.index = paste0("chr", c(18)));

message("\tmethylation levels ...\n"); 
#col_fun = colorRamp2(breaks = c(0, 1), colors = c("green", "red"))
col_fun = colorRamp2(c( -0.5, 0 ), c("grey","red"))
CpG$value[is.na(CpG$value)] <- -1
circos.genomicTrack(CpG, ylim=c(0,1),  
					panel.fun = function(region, value, ...) {
						i = getI(...)
						circos.genomicRect(region, abs(value), ytop.column = i, ybottom = i-0.98, 
										   col = col_fun(value[[1]]), border= NA, ...)
					
					cell.xlim = get.cell.meta.data("cell.xlim")
					for(h in c(0, 0.2,0.4,0.6,0.8,1)) {
						    circos.lines(cell.xlim, c(h, h), col = "#00000040")
					  }
					})

#---------------------------------------------------------------------------------------------------

message("\tUMRs ...\n");
names(UMRs) <- c("chr", "start", "end", "value")
circos.genomicTrackPlotRegion(UMRs, ylim = c(0, 100), panel.fun = function(region, value, ...) {
  col = ifelse(value[[1]] > 0, "#99EDC3", "#E41A1C")
  circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
  cell.xlim = get.cell.meta.data("cell.xlim")
  for(h in c(0, 20, 40, 60, 80, 100)) {
    circos.lines(cell.xlim, c(h, h), col = "#00000040")
  }
}, track.height = 0.2);

message("\tLMRs ...\n"); 
names(LMRs) <- c("chr", "start", "end", "value")
circos.genomicTrackPlotRegion(LMRs, ylim = c(0, 100), panel.fun = function(region, value, ...) {
  col = ifelse(value[[1]] > 0, "skyblue", "#E41A1C")
  circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
  cell.xlim = get.cell.meta.data("cell.xlim")
  for(h in c(0, 20, 40, 60, 80, 100)) {
    circos.lines(cell.xlim, c(h, h), col = "#00000040")
  }
}, track.height = 0.2);


#-----------------------------------------------------------------------------------------------------
dev.off()
