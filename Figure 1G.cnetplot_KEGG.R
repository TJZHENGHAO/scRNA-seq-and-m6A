library("ggplot2")
library("multienrichjam")
library(clusterProfiler)
Args <- commandArgs(trailingOnly=T)
DF <- Args[1]
gene_m <- Args[2]
out <- Args[3]

enrichDF <- read.table(DF, sep="\t", header = T, as.is = T)
gene <- read.delim(gene_m, header = T,as.is = T)

colnames(enrichDF)[c(8,11)] <- c("SYMBOL", "geneID")
print(head(enrichDF,2))
enrichResult <- enrichDF2enrichResult(enrichDF)
head(enrichResult)

#tmp <- log2(gene[,2])
tmp <- gene[,2]
names(tmp) <- gene[,1]
FClist <- tmp
head(FClist)

CategoryNum <- nrow(enrichDF)

p <- cnetplot(enrichResult,showCategory=CategoryNum,foldChange=FClist,circular=TRUE,colorEdge=TRUE) +scale_color_gradient2(name = "Log2 fold change", low = "green", mid = "gray", high = "red",midpoint = 0)  + guides(colour = guide_colourbar(order = 1),  size = guide_legend(order = 2) ) +  theme(plot.margin =  margin(0.1, 0.5, 0.1, 0.5, "cm"))  #scale_colour_continuous(name = "log2fc") #guides(color = guide_legend(order = 1)) + guides(color = guide_colourbar(title = "Log2 fold change"))

pdf(paste0(out, ".cnetplot.pdf"), width=10,height=8)
print(p)
dev.off()


tiff(paste0(out, ".cnetplot.tiff"), height = 8 ,width = 10, compression="lzw", units="in", res=300,pointsize=8)
print(p)
dev.off()

