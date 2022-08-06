library("pheatmap")

Args <- commandArgs(trailingOnly=TRUE)
file=Args[1]
phe=Args[2]

data <- read.delim(file, header=T, sep="\t", row.names=1,check.names=F)
colData <- read.delim(phe,sep=",",header = T,row.names=1)

data <- data[,rownames(colData)]
heatmap_data <- as.matrix(data)
col <- colorRampPalette(c("blue", "white", "red"))(256)

annotation_col =data.frame(CellType = factor(colData[,1]))  
rownames(annotation_col) = rownames(colData)

pdf("heatmap.pdf",width = 10, height =8 , onefile = FALSE)
pheatmap(heatmap_data,cellwidth = NA, cellheight =NA,  treeheight_col = 30 ,treeheight_row = 30 ,
         color = col, scale="row", legend=TRUE,border_color=NA ,
         fontsize_row=10, fontsize_col=15,main="Heatmap",annotation_col = annotation_col,
         #annotation_colors = ann_colors,
         show_rownames = F,cluster_cols = F, cluster_rows=T, fontfamily="serif",show_colnames = FALSE)
dev.off()
