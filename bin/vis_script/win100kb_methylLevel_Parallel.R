library("MASS");

if(!require(doParallel)) {
	      install.packages("doParallel")
}
if(!require(foreach)) {
	      install.packages("foreach")
}

library(doParallel)
library(foreach)

options(scipen=999)

##### arguments
args <- commandArgs(TRUE);
cytoband.file <- args[1];
path <- args[2]
core_num <- args[3]

myCluster <- parallel::makeCluster(core_num,type = "FORK")
doParallel::registerDoParallel(myCluster)

bin_size <- 100000

input_path <- paste(path, "union_CpG/", sep = '')
cytoband.df = read.table(cytoband.file, colClasses = c("character", "numeric", "numeric", "character", "character"), sep = "\t");
colnames(cytoband.df) <- c("chr", "start", "end", "ex1", "ex2")

out_CpGF <- paste(path, "window.CpG_methylcalls.bed", sep ="")


foreach::foreach( ch =1:nrow(cytoband.df), .combine=c, packages=c("MASS") )%dopar% {
	message("here\n")
	input_file <- paste(input_path, cytoband.df$chr[ch], ".txt" ,sep = "" )
	chr_data <- read.table(input_file, header=F, colClasses=c("character", "numeric", "numeric", "numeric"))
	colnames(chr_data) <- c("chr", "pos", "all", "meth")
	chr_data$level <- (data$methyl/data$all)*100

	CpG_levels <- c()
	win <- seq(0,cytoband.df$end[ch]+bin_size, by=bin_size)
	start <- win[1:length(win)-1]
	end <- win[2:length(win)]
	end[length(end)] <- cytoband.df$end[ch]
	
	for(i in 1:length(win)-1)
	{ 
		CpG_levels[i] <- mean( chr_data$level[ (chr_data$pos-1) < end[i] & (chr_data$pos -1)>= start[i]] )
	}
	
	chr <- as.factor( rep(ch, length(start)) )
	CpG_windows <- data.frame(chr,start, end, CpG_levels)
	write.table( CpG_windows, file=out_CpGF, append=TRUE,row.names=FALSE, col.names = FALSE,sep="\t" )

}

