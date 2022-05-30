  
suppressMessages(library('gprofiler2'))
suppressMessages(library('BSgenome'))

args <- commandArgs(TRUE);

f_genelist <- args[1];
assembly_ver <- args[2];
out <- args[3];
threshold <- 0.05;


## convert ucsc name -> organism name
av_gen <- available.genomes(splitNameParts=TRUE);
spc <- av_gen[av_gen[,4]==assembly_ver,2][1]
spc <- tolower(as.character(spc))

## read gene list
genelist <- read.table(f_genelist,header=F);

## gost
GO_result <- gost(query = as.vector(genelist$V1), organism=spc, user_threshold=threshold,  evcodes=TRUE,significant = FALSE );
df.GO_result <- apply(GO_result$result, 2, as.character);
write.table(df.GO_result,file=paste0(out, ".GOresult.txt"), quote=F, sep="\t");

res_pdf <- gostplot(GO_result, capped = FALSE, interactive = FALSE);
pdf(paste0(out, ".GOresult.pdf"));
res_pdf;
dev.off();
