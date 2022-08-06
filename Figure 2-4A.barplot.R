library("reshape2")
library("ggpubr")
library("ggplot2")
library("ggsci")
library("dplyr")
library("plyr")


Args <- commandArgs(trailingOnly=TRUE)
phe <- Args[1]
dat <- Args[2]
out <- Args[3]

#../SuHx_vs_Control/pheno_data.csv ../m6a_expr.txt
pheno <- read.csv(phe, row.names=1, header=T, as.is=T, check.names=F)
Dat <- read.delim(dat, header=T, as.is=T, check.names=F)
Dat <- Dat[,c("gene",row.names(pheno))]

print(head(pheno))
print(head(Dat[1:2,1:2]))

Dat_m <- melt(Dat)
Dat_m$group <- pheno[Dat_m[,2],]
colnames(Dat_m)[3] <- "Expression"
head(Dat_m)
Stat <- ddply(Dat_m, .(gene ,group), summarise, mean = mean(Expression), sd = sd(Expression))
Stat$mean_sd <- Stat$mean + Stat$sd
#ceiling(max(Stat[,ncol(Stat)]))
Ymax <- max(Stat[,ncol(Stat)]) +0.3

p <- ggbarplot(Dat_m, x = "gene", y = "Expression", color = "group", palette = "nejm", add = "mean_sd", position = position_dodge(0.78),  add.params = list(width = 0.5))  + #add = "jitter" outlier.shape = NA
    theme(axis.text.x=element_text(angle=45,hjust=1, size=10))

# Use significance symbol as label
p1 <- p + stat_compare_means(aes(group = group), label = "p.signif", label.y=Ymax, symnum.args  = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")))

ggsave(paste0(out,".barplot3.pdf"), height = 5, width=8)


tiff(paste0(out, ".barplot3.tiff"), height = 5 ,width = 8, compression="lzw", units="in", res=300,pointsize=8)
print(p1)
dev.off()

