library(monocle)
library(gplots)
library(ggplot2)
library(dplyr)
library(ggsci)
library(patchwork)

#getOption('timeout')
#options(timeout=10000)
#InstallData("pbmc3k") 

pdmat <- read.csv("pheno_data.csv", check.names=F, as.is=T, row.names=1)
#fdmat <- read.delim("fd.txt", sep="\t", check.names=F, as.is=T)
ctmat <- read.delim("filter_value.txt", sep="\t", check.names=F, as.is=T, row.names=1)
fdmat <- ctmat[,1,drop=F]
fdmat[,1] <- row.names(fdmat)
colnames(fdmat) <- "gene_short_name"

print(head(pdmat))
print(head(ctmat[1:2,1:3]))
print(head(fdmat))

pd <- new("AnnotatedDataFrame", data=pdmat)
fd <- new("AnnotatedDataFrame", data=fdmat)
#ct=as.data.frame(sce@assays$RNA@counts)
ct <- ctmat
ct[1:4,1:4]


sc_cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData =fd,
expressionFamily = tobit(Lower=0.1), ## diff in RPKM count TPM format
lowerDetectionLimit=0.1)


cds <- sc_cds
save(cds,file = 'input_cds.Rdata')


##Step 2
rm(list=ls())
options(stringsAsFactors = F)

load(file = 'input_cds.Rdata')
ordering_genes <- read.delim("id2", header=F) 
cds <- setOrderingFilter(cds, ordering_genes)
#plot_ordering_genes(cds)
cg <- head(ordering_genes)

# 第二步: 降维。降维的目的是为了更好的展示数据。函数里提供了很多种方法,
# 不同方法的最后展示的图都不太一样, 其中“DDRTree”是Monocle2使用的默认方法
cds <- reduceDimension(cds, max_components = 2,
					    method = 'DDRTree')
# 第三步: 对细胞进行排序
cds <- orderCells(cds)
# 最后两个可视化函数，对结果进行可视化
plot_cell_trajectory(cds, color_by = "Cluster") 
ggsave('monocle_cell_trajectory_for_seurat.pdf')

length(cg)
plot_genes_in_pseudotime(cds[cg,],
						  color_by = "Cluster") 
ggsave('monocle_plot_genes_in_pseudotime_for_seurat.pdf')

# https://davetang.org/muse/2017/10/01/getting-started-monocle/

my_cds_subset=cds
# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(my_cds_subset))
# 这个differentialGeneTest会比较耗费时间
my_pseudotime_de <- differentialGeneTest(my_cds_subset,
										  fullModelFormulaStr = "~sm.ns(Pseudotime)",
										   cores = 8)
# 不知道为什么无法开启并行计算了

head(my_pseudotime_de)
save( my_cds_subset,my_pseudotime_de,
	  file = 'output_of_phe2_monocle.Rdata')

##Step3
## 后面是对前面的结果进行精雕细琢 
rm(list=ls())
options(stringsAsFactors = F)
load(file = 'output_of_phe2_monocle.Rdata')
cds=my_cds_subset
colnames(pData(cds))
table(pData(cds)$State,pData(cds)$Cluster)
p1=plot_cell_trajectory(cds, color_by = "Cluster") + scale_color_nejm() 
ggsave('trajectory_by_cluster.pdf')
#plot_cell_trajectory(cds, color_by = "celltype") 
p2=plot_cell_trajectory(cds, color_by = "Pseudotime") 
ggsave('trajectory_by_Pseudotime.pdf')
p3=plot_cell_trajectory(cds, color_by = "State") + scale_color_npg()

ggsave('trajectory_by_State.pdf')
#p1+p2/p3
#ggsave('trajectory_by_all.pdf')

