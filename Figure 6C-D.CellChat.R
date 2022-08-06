library("CellChat")

x <- read.delim("data.unique.xls",row.names=1,header=T,sep="\t",check.names=F)
colData <- read.delim("Control_meta.txt", head=T, as.is=T,row.names=1,check.names=F)

x <- x[, rownames(colData)]
#keep <- rowSums(x) > 0
cellchat <- createCellChat(object = as.matrix(x), meta = colData , group.by = "cell_type")

#cellchat <- cco.pbmc
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat, "Control.rds")


Control <- readRDS("Control.rds")
SuHx <- readRDS("SuHx.rds")
MCT <- readRDS("MCT.rds")

cco.list <- list(Control=Control, SuHx=SuHx)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)

#p <- netVisual_bubble(cellchat,  comparison = c(1, 2), angle.x = 60)
#ggsave("Compare_LR_bubble.pdf")

gg1 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in SuHx", angle.x = 90, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in SuHx", angle.x = 90, remove.isolate = T)
#> Comparing communications on a merged object

pdf("SuHx_vs_Control.Compare_LR_bubble.pdf", width =15, height=7)
print(gg1 + gg2)
dev.off()
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = "SuHx", features.name = "SuHx", only.pos = FALSE, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = "SuHx")
write.table(net, "SuHx_vs_Control.xls", quote=F, sep= "\t")

cco.list <- list(Control=Control, MCT=MCT)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)
gg1 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in MCT", angle.x = 90, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in MCT", angle.x = 90, remove.isolate = T)
#> Comparing communications on a merged object

pdf("MCT_vs_Control.Compare_LR_bubble.pdf", width =15, height=7)
print(gg1 + gg2)
dev.off()

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = "MCT", features.name = "MCT", only.pos = FALSE, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = "MCT")
write.table(net, "MCT_vs_Control.xls", quote=F, sep= "\t", row.names=F)

