#args[1] = parameter file, args[2] = <control,case> , args[3] = outpath, args[4] = core(defalt=2)

suppressMessages(library(bsseq))
suppressMessages(library(BiocParallel))


#------------------------------------------------------
args=commandArgs(TRUE)
allsamples <- read.table(args[1])
#out_name = paste0(args[2], '-', args[3], sep='')

compare_set <- args[2]
outpath <- args[3]

q_cutoff <- as.numeric( args[4] )

core_num <- 2
if (length(args)>=5){
	core_num <- as.numeric(args[5])
}


# -------------------------------------------------------------
colnames(allsamples) = c('group','name', 'file')


chr_list <- readLines( paste(outpath, '/chr_list.txt',sep='') )


control_name <- strsplit(compare_set, ',')[[1]][1]
case_name <- strsplit(compare_set, ',')[[1]][2]

control_set <- allsamples[allsamples[,1]==control_name,2]
case_set <- allsamples[allsamples[,1]==case_name,2]

samples <-  allsamples[allsamples[,1]==control_name,]
samples <- rbind(samples,  allsamples[allsamples[,1]==case_name,])
f_name_list <- samples[,2]

# -----------------------------------------------------
print('< Read bismark >')
bs_list <- c()

bs.tmp0 <- read.bismark(files = samples[,3], colData = DataFrame(row.names = f_name_list))
bs.all <- chrSelectBSseq(bs.tmp0, seqnames = chr_list, order = TRUE)


print('< Smoothing >')
bs.all.fit <- BSmooth(BSseq=bs.all, BPPARAM=MulticoreParam(workers=core_num, progressbar=TRUE))

print('< t-Stat >')
bs.all.tstat <- BSmooth.tstat(bs.all.fit, group1 = control_set, group2 = case_set, estimate.var = "group2", verbose=TRUE,  local.correct = FALSE)


print('< Find DMR >')
print(q_cutoff)
bs.dmrs0 <- dmrFinder(bs.all.tstat, qcutoff=c(q_cutoff,1-q_cutoff) , stat='tstat')
bs.dmrs <- subset(bs.dmrs0, n>=3 & abs(meanDiff)>=0.1)
#bs.dmrs <- subset(bs.dmrs0)
rm(bs.all.tstat, bs.dmrs0)
write.table(bs.dmrs, file=paste0(outpath,'/DMR/',control_name,".",case_name,"/","DMR_q",q_cutoff,".bed"), quote=FALSE, row.names=FALSE, sep="\t" )


print('< Finish >')


