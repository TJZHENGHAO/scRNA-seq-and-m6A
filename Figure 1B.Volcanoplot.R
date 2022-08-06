library('getopt')

command=matrix(c('inputfile', 'i', 1,'character', 'edgeR_diff.xls file',
		 'pq','q', 1, 'character', 'available options are : pvalue or qvalue',
		 'sig','s', 1, 'double', 'significant threshold',
		 'FC','f', 1, 'double', 'fold change threshold',
                 'ncolFC', 'd', 1, 'double', 'FC column',
                 'ncolpq', 'c', 1, 'double', 'pvalue or qvalue column'),
	       byrow=T,ncol=5)

args=getopt(command)

if (is.null(args$inputfile) || is.null(args$pq) ||  is.null(args$sig) ||  is.null(args$FC) || is.null(args$ncolFC) || is.null(args$ncolpq) ) {
	cat(paste(getopt(command, usage = T), "\n"))
	q(status=1)
}

print(args$inputfile)
print(args$pq)
print(args$sig)
print(args$FC)
print(args$ncolFC)
print(args$ncolpq)


##2
library("ggplot2")
library("pheatmap")

data <- read.delim(args$inputfile,header=T,sep="\t",as.is=T,check.names = F)

Pvalue_or_Qvalue = args$pq
significant_threshold=args$sig
FC_threshold = args$FC
FC_column = args$ncolFC
pvalue_or_qvalue_column = args$ncolpq

#Pvalue_or_Qvalue_threshold = 0.05
#log2FC_threshold = 1

#data$log2fc <- log2(data[,9])
data$log2fc <- data[,FC_column]

if(Pvalue_or_Qvalue == "pvalue")
{
	data$minuslog10pvalue <- log10(data[,pvalue_or_qvalue_column])*-1
#	data$threshold <- as.factor(ifelse(data[,10] < 0.05 & abs(log2(data[,9])) >=1 ,ifelse((log2(data[,9])) >1 ,'Up','Down'),'Not'))
	data$threshold <- as.factor(ifelse(data[,pvalue_or_qvalue_column] < significant_threshold & (data[,FC_column] > FC_threshold | data[,FC_column] < -1*FC_threshold) ,ifelse(data[,FC_column] > FC_threshold ,'Up','Down'),'Not'))
}else{
	data$minuslog10pvalue <- log10(data[,pvalue_or_qvalue_column])*-1
	data$threshold <- as.factor(ifelse(data[,pvalue_or_qvalue_column] < significant_threshold & (data[,FC_column] > FC_threshold | data[,FC_column] < -1*FC_threshold) ,ifelse(data[,FC_column] > FC_threshold ,'Up','Down'),'Not'))
}

# as.factor(ifelse(data[,10] < 0.05 & abs(log2(data[,9])) >=1 ,ifelse((log2(data[,9])) >1 ,'Up','Down'),'Not'),levels=c('Up','Down','Not'))
data$threshold1 <- ""
data$threshold1 <- factor(ifelse(data[,ncol(data)-1] == "Up" | data[,ncol(data)-1] == "Down" ,
                                    ifelse( data[,ncol(data)-1] == "Up" ,paste("Up: ",table(data[,ncol(data)-1])["Up"]),paste("Down: ",table(data[,ncol(data)-1])["Down"])),
                                    paste("no-Peaks: ",table(data[,ncol(data)-1])["Not"])),
                          levels = c(paste("Up: ",table(data[,ncol(data)-1])["Up"]),paste("Down: ",table(data[,ncol(data)-1])["Down"]),
                                     paste("no-Peaks: ",table(data[,ncol(data)-1])["Not"])))

data <- data[order(data[,ncol(data)],decreasing = TRUE),]

#mark <- na.omit(data[data[,21] == "TXNIP",])
#print(mark)
## volcano plot
p <- ggplot(data=data,mapping=aes(x=log2fc, y=minuslog10pvalue,
           colour=threshold1,fill=threshold1))+geom_point(size=0.01)+
  scale_color_manual(values=c("red", "blue","gray"))+
  #geom_point() +
  #xlim(c(-4, 4)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-FC_threshold,FC_threshold),lty=4,col="black",lwd=0.6)+
  geom_hline(yintercept = -log10(significant_threshold),lty=4,col="black",lwd=0.6)+
  labs(x=expression(log["2"]*"(fold change)"),y=bquote(paste('-log'['10']*'('*'', .(Pvalue_or_Qvalue),")")),title="",family = "Times")+

  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=10),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold" , color="black", size=12),
        axis.text.y = element_text(face="bold" ,  color="black", size=12),
        axis.title.x = element_text(face="bold" ,color="black", size=12),
        axis.title.y = element_text(face="bold" ,color="black", size=12)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
#p + geom_point(data=mark,mapping=aes(x=log2fc, y=minuslog10pvalue, colour=threshold1,fill=threshold1 ),  size=1, color = "orange" )  
ggsave("volcano-plot-3.pdf", width = 6, height = 5)


