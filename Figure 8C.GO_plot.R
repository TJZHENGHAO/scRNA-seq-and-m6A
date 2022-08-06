library("ggplot2")
library("Cairo")


Args <- commandArgs(trailingOnly=T)
DF <- Args[1]
out <- Args[2]

GO <- read.delim(DF, as.is=T)
GO1 <- GO[order(GO[,1]),]
print(GO1[,1:2])
Description_sort <- factor(GO1$Description,levels = rev(unique(GO1$Description)))

p <- ggplot(data = GO1,mapping = aes(x = Description_sort, y = -log10(pvalue),fill=Ontology) )+ geom_bar(stat = "identity") + coord_flip() +
	   geom_text(aes(label = Count), size = 6, hjust=-0.2) +
	   theme(panel.background=element_rect(color="black",fill="white")) +
	#facet_grid(. ~ Ontology,scales = "free", space = "free") + theme(axis.text.x=element_text(angle=60,colour="black",size=18,hjust=1)) +
	theme(axis.text=element_text(colour="black",size=18)) +
	theme(axis.title=element_text(colour="black",size=20)) + labs(x="Description", y= "-log10(Pvalue)", fill = "") +
	theme(legend.text = element_text( size = 14),legend.title=element_text( size = 18))+
	scale_fill_manual(values = c("biological process" = "red", "cellular compoment" = "green", "molecular function" ="blue")) +
	ylim(0, max(-1*log10(GO1$pvalue))+0.5)

pdf(paste0(out, ".barplot.pdf"), width = 11, height = 8)
p
dev.off()

tiff(paste0(out, ".barplot.tiff"), height = 8 ,width = 11, compression="lzw", units="in", res=300,pointsize=8)
print(p)
dev.off()

