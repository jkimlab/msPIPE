library("MASS");
##### arguments
args <- commandArgs(TRUE);
cytoband.file <- args[1];
path <- args[2];

outpath <- args[3]
#outpath <- path
bin_size <- 100000


options(scipen=999)
CpGF <- paste(path, "CpG_methylCalls.bed", sep='')

message("intput data:\n")
message(CpGF)


cytoband.df = read.table(cytoband.file, colClasses = c("character", "numeric", "numeric", "character", "character"), sep = "\t");

message("\nloading data...\n")
CpGdata <- read.table(CpGF, header=F, colClasses = c("character", "numeric", "numeric", "character", "character", "character"));

message('done\n')

getLevels <- function(fraction){
	unmethyl<-sapply(strsplit(as.character(fraction),';'), "[", 1)
	unmethyl<-as.numeric( sapply(strsplit(as.character(unmethyl),':'), "[", 2) )
	

	methyl<-sapply(strsplit(as.character(fraction),';'), "[", 2)
	methyl<-as.numeric( sapply(strsplit(as.character(methyl),':'), "[", 2) )
	
	values <- methyl/(unmethyl+methyl)
	return(values)
}

out_CpGF <- paste(outpath,"window_CpGMethylCalls.bed",sep="")

message("output file:\n")
message(out_CpGF)

for(i in 1:length(cytoband.df$V1)){
	ch = cytoband.df$V1[i]
	CpG <- CpGdata[CpGdata$V1 == ch,]
	
	CpG$level <- getLevels(CpG$V4)

	win <- seq(0,cytoband.df$V3[i]+bin_size, by=bin_size)
	start <- win[1:length(win)-1]
	end <- win[2:length(win)]
	end[length(end)] <- cytoband.df$V3[i]
	
	CpG_levels <- c()
	
	for(i in 1:length(win)-1)
	{ 
		if(i==1){
			cat( CpG$level[CpG$V2 < end[i]& CpG$V2>= start[i]] )
			stop()
		}
		CpG_levels[i] <- mean( CpG$level[CpG$V2 < end[i]& CpG$V2>= start[i]] )
	}
	
	chr <- as.factor( rep(ch, length(start)) )
	CpG_windows <- data.frame(chr,start, end, CpG_levels)


	write.table( CpG_windows, file=out_CpGF, append=TRUE,row.names=FALSE, col.names = FALSE,sep="\t" )

}


warnings()

