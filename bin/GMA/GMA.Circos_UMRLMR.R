#!/usr/bin/env Rscript
#install.packages("circlize")
library("circlize");
library("MASS");
##### arguments
args <- commandArgs(TRUE);
cytoband <- args[1];
umrs <- args[2];
lmrs <- args[3];
out_pdf <- args[4];

message("loading data...\n");
message("\tUMRs...");
bed1 <- read.table(umrs, header=F);
message("done\n\tLMRs...");
bed2 <- read.table(lmrs, header=F);
message("done\ndrawing pdf..\n");
pdf(out_pdf, width=20, height=20);
circos.clear();

message("\tcytoband..\n");
cytoband.file = cytoband;
cytoband.df = read.table(cytoband.file, colClasses = c("character", "numeric", "numeric", "character", "character"), sep = "\t");
circos.initializeWithIdeogram(cytoband.df);
#circos.initializeWithIdeogram(cytoband.df, chromosome.index = paste0("chr", c(18)));

message("\tUMRs ...\n");
names(bed1) <- c("chr", "start", "end", "value")
circos.genomicTrackPlotRegion(bed1, ylim = c(0, 100), panel.fun = function(region, value, ...) {
  col = ifelse(value[[1]] > 0, "#99EDC3", "#E41A1C")
  circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
  cell.xlim = get.cell.meta.data("cell.xlim")
  for(h in c(0, 20, 40, 60, 80, 100)) {
    circos.lines(cell.xlim, c(h, h), col = "#00000040")
  }
}, track.height = 0.1);

message("\tLMRs ...\n"); 
names(bed2) <- c("chr", "start", "end", "value")
circos.genomicTrackPlotRegion(bed2, ylim = c(0, 100), panel.fun = function(region, value, ...) {
  col = ifelse(value[[1]] > 0, "skyblue", "#E41A1C")
  circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
  cell.xlim = get.cell.meta.data("cell.xlim")
  for(h in c(0, 20, 40, 60, 80, 100)) {
    circos.lines(cell.xlim, c(h, h), col = "#00000040")
  }
}, track.height = 0.1);
dev.off();
