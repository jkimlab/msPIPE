

if(!require(ggplot2)) {
      install.packages("ggplot2")
}
library(ggplot2)

args <- commandArgs(TRUE)

file <- args[1]
sample_name <- args[2]
out_name <- args[3]

message("reading file...\t")
df <- read.table(file)



message("calculate average meth levels")
meth <- data.frame(do.call('rbind', strsplit(as.character(df$V10), split=';', fixed=T)))

unmeth <- data.frame(as.numeric(data.frame(do.call('rbind', strsplit(as.character(meth$X1), split=':', fixed=T)))$X2))
meth <- data.frame(as.numeric(data.frame(do.call('rbind', strsplit(as.character(meth$X2), split=':', fixed=T)))$X2))

all_methlevel<- cbind(as.factor(df$V4), 100*(meth/(meth+unmeth)))
colnames(all_methlevel) <- c("context", "level")

genomic_context <- data.frame(levels(factor(all_methlevel$context)))
colnames(genomic_context)<- c("context")

for(i in seq(1,5)){
  genomic_context$avg_level[i] <- mean(all_methlevel$level[
    all_methlevel$context==genomic_context$context[i]])
}

print(genomic_context)

p <- ggplot(genomic_context, aes(context,avg_level))+geom_bar(stat="identity", width=0.7) + theme_minimal()+
  lims(y=c(0,100)) +
  ggtitle( paste("\n",sample_name,"\n", sep=" ") )+
  theme(plot.title=element_text(hjust=0.5, size=55,face="bold"),
        axis.title=element_text(size=40,face="bold"),
        axis.text =element_text(size=40)
  ) +
  labs(y="avg levels (%)") +
  scale_x_discrete(limits=c("promoter","gene","exon", "intron", "intergenic"))

out_pdf <- paste(out_name ,sep="")
ggsave(out_pdf, width=20, height = 20)
