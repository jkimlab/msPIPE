library(scales)
library(ggplot2)

if(!require(doParallel)) {
	  install.packages("doParallel")
}
if(!require(foreach)) {
	  install.packages("foreach")
}

library(doParallel)
library(foreach)


args <- commandArgs(TRUE)

vis_param <- args[1]
Ideo_file <- args[2]
outpath <- args[3]
core_num <- as.numeric( args[4] )

#win_size <-as.numeric( args[5] )
win_size <- 100000

outpath <- gsub( "/$","",  outpath)

myCluster <- parallel::makeCluster(core_num,type = "FORK")
doParallel::registerDoParallel(myCluster)

options(digits=3)

draw_hist <-  function(df, cx, name ){
  if(cx=="CpG"){
    cx_col = "#E41A1C"
  } else if (cx == "CHG"){
    cx_col = "#377EB8"
  }else if (cx=="CHH"){
    cx_col = "#4DAF4A"
  }
	
  ggplot(df, aes(x=level)) + theme_minimal()+
	  geom_histogram(aes(y=..density..), alpha=0.9, col="black", fill=cx_col,binwidth = 11, boundary = 10)+
	  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
	  ggtitle( paste("\n",name, cx, "histogram\n", sep=" ") )+ 
	  theme(plot.title=element_text(hjust=0.5, size=55,face="bold"),
			axis.title=element_text(size=40,face="bold"),
			axis.text =element_text(size=30),
			panel.grid.minor.x = element_blank()
			) + 
    labs(x="levels (%)", y="density")
  
  out_pdf <- paste(outpath, "/", name, "/hist_",name,"_",cx,".pdf" ,sep="")
  ggsave(out_pdf, width=20, height = 20)

}

window_levels <- function( cyto.df, data.df, cx,sample_name){
	out_CpGF <- paste(outpath, "/",sample_name , "/",cx, '_methylLev_window.bed', sep ='')
	file.create(out_CpGF)
	
	for(i in 1:length(cyto.df$chr)){
    	ch = cyto.df$chr[i]
    	chr_data <- data.df[data.df$chr == ch,]


    	win <- seq( 1, cyto.df$end[i]+win_size, by=win_size)
    	start <- win[1:length(win)-1]
    	end <- win[2:length(win)]
    	end[length(end)] <- cyto.df$end[i]

		mean_levels <- c()
	
		for(i in 1:length(win)-1){
			mean_levels[i] <- round( mean( chr_data$level[chr_data$pos-1 < end[i]& chr_data$pos-1 >= start[i] ]),3 )
		}
		chr <- as.factor( rep(ch, length(start)) )
		start <- start-1
		end <- end -1
		windows_data <- format( data.frame(chr,start, end, mean_levels), scientific = FALSE)
		write.table( windows_data, file=out_CpGF, append=TRUE,row.names=FALSE, col.names = FALSE,sep="\t", quote=FALSE)
	}

}

data_param <- read.table(vis_param, header=TRUE)
Ideo.df <- read.table(Ideo_file, colClasses = c("character", "numeric", "numeric", "character", "character"), sep = "\t");
colnames(Ideo.df) <- c('chr','start','end','ex1', 'ex2')


## READ FILES for PARALLEL
mean_levels <- c()

mean_levels <- foreach::foreach( i =1:nrow(data_param), .combine=c, .packages=c("ggplot2" ) ) %dopar% {
  file <- toString(data_param$file[i])
  data <- read.table( file , header=FALSE, col.names = c("chr", "pos", "all", "methyl"))
  data$level <- (data$methyl/data$all)*100
  
  if (data_param$context[i] == "CpG"){
	  
	  #make file
	  name <- data_param$sample[i]
	  wig_file_name <- paste(outpath, "/", name, "/CpG_methylLev.wig", sep='')
	  file.create(wig_file_name)
  	  
	  #window_levels( Ideo.df, data, data_param$context[i] ,data_param$sample[i] )
	  data$level <- round(data$level, 3)
	  for(chr_i in 1:length(Ideo.df$chr)){
		  chr_data <- data[data$chr == Ideo.df$chr[chr_i],]
		  write(paste("variableStep chrom=", Ideo.df$chr[chr_i], sep=''), file= wig_file_name ,append=TRUE)
		  write.table(data.frame(chr_data$pos, chr_data$level), file= wig_file_name, append=TRUE,row.names=FALSE, col.names = FALSE,sep="\t", quote=FALSE )
	  }
  }

  draw_hist(data, data_param$context[i],data_param$sample[i])
  mean(data$level)

}


data_param$avg_level <- mean_levels

p <- ggplot(data_param,aes(x=context,y=avg_level,fill=sample) )+theme_minimal()+
	geom_bar(stat="identity",position="dodge", width=0.65 ,colour="grey70")+
	scale_x_discrete(limits=c("CpG", "CHG", "CHH")) +
	lims(y=c(0,100)) +
	scale_fill_brewer(palette = "Greens",direction = -1)

p1 <- p + theme(axis.title=element_text(size=40,face="bold"),
				axis.text = element_text(size=30),
				legend.title = element_text(size=40,face = "bold"),
				legend.text = element_text(size=40),
				legend.key.size = unit(1, 'cm'),
				legend.position = c(0.8, 0.8),
				legend.margin = margin(10, 10, 10, 10),
				legend.box.background = element_rect(fill = "white")
				) +
		labs(y="avg level(%)")

outpdf <- paste(outpath ,"/","avg_methlevel.pdf", sep = "")
ggsave(outpdf, plot=p1, width=30, height = 20)


