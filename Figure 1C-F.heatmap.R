library("Cairo")
library("ggplot2")


Args <- commandArgs(trailingOnly=T)
DF <- Args[1]
out <- Args[2]

kegg2 <- read.table(DF, sep="\t",header = T,as.is = T)
kegg2 <- kegg2[order(kegg2[,5]),]

for (i in 1:nrow(kegg2)) {kegg2[i,3] <- as.numeric(strsplit(kegg2[i,3],"/")[[1]][1])/as.numeric(strsplit(kegg2[i,3],"/")[[1]][2])}
kegg2[,3] <- as.numeric(kegg2[,3])

Description_sort <- factor(kegg2$Description,levels = rev(unique(kegg2$Description)))
 p <- ggplot(data = kegg2,mapping = aes(x = Description_sort,y= GeneRatio,size = Count,colour=-log10(pvalue)))+
  geom_point()  + 
  scale_size(range=c(5,14)) + 
  scale_colour_gradient(low="green",high = "red") +coord_flip()+
  theme(legend.text = element_text( size = 14),legend.title=element_text( size = 18), legend.box.background=element_blank() , legend.background=element_blank(), legend.key = element_rect(fill = "transparent")) + theme(legend.key.height = unit(1,'cm'))+
  theme(axis.text.x=element_text(size=18,color="black"),axis.text.y=element_text(size=18,color="black"),axis.title = element_text(size=22,color="black")) + 
  theme(panel.background=element_rect(color="black",fill="white")) +
  theme(panel.grid.major  = element_line(colour = "gray")) +
  labs(x="Description", colour = "-log10(Pvalue)")
 


pdf(paste0(out, ".dotplot.pdf"),width = 11, height = 8)
p
dev.off()

tiff(paste0(out, ".dotplot.tiff"), height = 8 ,width = 11, compression="lzw", units="in", res=300,pointsize=8)
print(p)
dev.off()

