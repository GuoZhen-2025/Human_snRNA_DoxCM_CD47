#### 1. Load packages ####
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(patchwork)
library(tidydr)
library(ggplot2)
library(cowplot)
library(magrittr)
library(CellChat)

#### 2. Set pathway ####
setwd("//storage2.ris.wustl.edu/guozhen/Active/DoxCM_CD47/") #local
#### 3. Read files ####

CM <- readRDS("./Final RDS/Final_CM_clustered.rds")
DimPlot(CM, dims=c(1,2), reduction = "umap", label = F, pt.size = 0.005, repel = T)
DotPlot(CM, features = "CD47", group.by = "celltype1")
FeaturePlot(CM, features = "CD47")
FB <- readRDS("./Final RDS/Final_FB_clustered.rds")
DimPlot(CM, dims=c(1,2), reduction = "umap", label = F, pt.size = 0.005, repel = T)


# DoxCM_Final <- readRDS("./DoxCM Output/Final_RDS/DoxCM_Final_cleanup_used_for_paper.rds")
# table(DoxCM_Final$celltype)

DoxCM_Final <- merge(CM, y = FB, project = "combined")

Idents(DoxCM_Final) <- "group"
DoxCM.object <- subset(DoxCM_Final, idents = "DoxCM")
saveRDS(DoxCM.object, "./Final RDS/Final_DoxCM_CM&FB_communication.rds")

Donor.object <- subset(DoxCM_Final, idents = "Donor")
saveRDS(Donor.object, "./Final RDS/Final_Donor_CM&FB_communication.rds")

NICM.object <- subset(DoxCM_Final, idents = "NICM")
saveRDS(NICM.object, "./Final RDS/Final_NICM_CM&FB_communication.rds")

#DoxCM.object <- subset(DoxCM.object, downsample = 20000)


DoxCM.object <- readRDS("./Final RDS/Final_DoxCM_CM&FB_communication.rds")
DoxCM.object <- readRDS("./Final RDS/Final_Donor_CM&FB_communication.rds")
DoxCM.data.input <- GetAssayData(DoxCM.object, assay = "SCT", slot = "data") # assay = "RNA", slot = "scale.data" or "counts"
DoxCM.meta <- DoxCM.object@meta.data[, c("celltype1","group")]

DoxCM.cellchat <- createCellChat(DoxCM.data.input)
DoxCM.cellchat <- addMeta(DoxCM.cellchat, meta = DoxCM.meta)
DoxCM.cellchat <- setIdent(DoxCM.cellchat, ident.use = "celltype1")

levels(DoxCM.cellchat@idents)
groupSize <- as.numeric(table(DoxCM.cellchat@idents))
groupSize

#. Cellchat database choose
DoxCM.cellchat@DB <- CellChatDB.human
#showDatabaseCategory(CellChatDB.human) # can not show in server
dplyr::glimpse(CellChatDB.human$interaction)

#. Standard
DoxCM.cellchat <- subsetData(DoxCM.cellchat, features = NULL)
future::plan("multisession", workers = 4) # modify according to your PC core
DoxCM.cellchat <- identifyOverExpressedGenes(DoxCM.cellchat)
options(future.globals.maxSize = 10000 * 1024^2) # options(future.globals.maxSize = 1000 * 1024^2)
DoxCM.cellchat <- identifyOverExpressedInteractions(DoxCM.cellchat)
# DoxCM.cellchat <- projectData(DoxCM.cellchat, PPI.human) optinal


DoxCM.cellchat <- computeCommunProb(DoxCM.cellchat, raw.use = T)
DoxCM.cellchat <- filterCommunication(DoxCM.cellchat, min.cells = 10)
DoxCM.cellchat <- computeCommunProbPathway(DoxCM.cellchat)
DoxCM.cellchat <- aggregateNet(DoxCM.cellchat)
DoxCM.cellchat <- netAnalysis_computeCentrality(DoxCM.cellchat, slot.name = "netP")

group.net <- subsetCommunication(DoxCM.cellchat)
write.csv(group.net, file = "./Output/Cell&Cell/Donor_Cellchat(CM&FB).csv")
saveRDS(DoxCM.cellchat, "./Output/Cell&Cell/Donor_Cellchat(CM&FB).rds")

DoxCM.object <- readRDS("./Final RDS/Final_NICM_CM&FB_communication.rds")
DoxCM.data.input <- GetAssayData(DoxCM.object, assay = "SCT", slot = "data") # assay = "RNA", slot = "scale.data" or "counts"
DoxCM.meta <- DoxCM.object@meta.data[, c("celltype1","group")]

DoxCM.cellchat <- createCellChat(DoxCM.data.input)
DoxCM.cellchat <- addMeta(DoxCM.cellchat, meta = DoxCM.meta)
DoxCM.cellchat <- setIdent(DoxCM.cellchat, ident.use = "celltype1")

levels(DoxCM.cellchat@idents)
groupSize <- as.numeric(table(DoxCM.cellchat@idents))
groupSize

#. Cellchat database choose
DoxCM.cellchat@DB <- CellChatDB.human
#showDatabaseCategory(CellChatDB.human) # can not show in server
dplyr::glimpse(CellChatDB.human$interaction)

#. Standard
DoxCM.cellchat <- subsetData(DoxCM.cellchat, features = NULL)
future::plan("multisession", workers = 7) # modify according to your PC core
DoxCM.cellchat <- identifyOverExpressedGenes(DoxCM.cellchat)
options(future.globals.maxSize = 10000 * 1024^2) # options(future.globals.maxSize = 1000 * 1024^2)
DoxCM.cellchat <- identifyOverExpressedInteractions(DoxCM.cellchat)
# DoxCM.cellchat <- projectData(DoxCM.cellchat, PPI.human) optinal


DoxCM.cellchat <- computeCommunProb(DoxCM.cellchat, raw.use = T)
DoxCM.cellchat <- filterCommunication(DoxCM.cellchat, min.cells = 10)
DoxCM.cellchat <- computeCommunProbPathway(DoxCM.cellchat)
DoxCM.cellchat <- aggregateNet(DoxCM.cellchat)
DoxCM.cellchat <- netAnalysis_computeCentrality(DoxCM.cellchat, slot.name = "netP")

group.net <- subsetCommunication(DoxCM.cellchat)
write.csv(group.net, file = "./Output/Cell&Cell/NICM_Cellchat(CM&FB).csv")
saveRDS(DoxCM.cellchat, "./Output/Cell&Cell/NICM_Cellchat(CM&FB).rds")



library(svglite)
#. Visulation
DoxCM.cellchat <- readRDS("./Output/Cell&Cell/Donor_Cellchat(CM&FB).rds")
pdf("./Output/Cell&Cell/Overall_interactions_number.pdf", width = 5, height = 4)
par(mfrow = c(1, 1), xpd = T)
netVisual_circle(DoxCM.cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = "Number of interactions")
dev.off()

pdf("./Output/Cell&Cell/Overall_interactions_weight.pdf", width = 5, height = 4)
netVisual_circle(DoxCM.cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = "Interaction weights/strength")
dev.off()



pdf("Overall_Interaction1.pdf", width = 15, height = 8)
dev.off()


pdf("./Output/Cell&Cell/Individual_interactions.pdf", width = 5, height = 4)
mat <- DoxCM.cellchat@net$weight
par(mfrow = c(3,3), mar = c(1,1,1,1))
for (i in 1:nrow(mat)){
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,
                   vertex.weight = groupSize,
                   arrow.width = 0.2, arrow.size = 0.1,
                   weight.scale = T,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}
dev.off()

dev.off()
par(mfrow = c(2,2), mar = c(1,1,1,1))
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[15, ] <- mat[15, ]
netVisual_circle(mat2,
                 vertex.weight = groupSize,
                 arrow.width = 0.2, arrow.size = 0.1,
                 weight.scale = T,
                 edge.weight.max = max(mat),
                 title.name = rownames(mat)[15])
# ggsave("Interested_cellchat (15x8).pdf", width = 15, height = 8, path = "./DoxCM Output/FBs/Revision/")

#. Pathway
DoxCM.cellchat@netP$pathways
pathway.show <- group.net$pathway_name
pathway.show <- c("IGF")
levels(DoxCM.cellchat@idents)
vertex.receiver <- seq(1,4)
#vertex.receiver <- c(1,7,8)
pdf("./Output/Cell&Cell/Overall_pathway_hierarchy.pdf", width = 7, height = 5)
netVisual_aggregate(DoxCM.cellchat,
                    signaling = pathway.show,
                    vertex.receiver = vertex.receiver,
                    layout = "hierarchy")
dev.off()

pdf("./Output/Cell&Cell/Overall_pathway_circle.pdf", width = 7, height = 5)
par(mfrow = c(1,1))
netVisual_aggregate(DoxCM.cellchat,
                    signaling = pathway.show,
                    layout = "circle")
dev.off()



pdf("./Output/Cell&Cell/FGF_pathway_chord.pdf", width = 4, height = 4)
par(mfrow = c(1,1))
netVisual_aggregate(DoxCM.cellchat,
                    signaling = "FGF",
                    layout = "chord",
                    title.space = 0)
dev.off()

netVisual_heatmap(DoxCM.cellchat,
                  signaling = pathway.show,
                  color.heatmap = "Reds")

pdf("./Output/Cell&Cell/IGF_pathway_contribution.pdf", width = 3, height = 2)
netAnalysis_contribution(DoxCM.cellchat, signaling = pathway.show)
dev.off()

pair <- extractEnrichedLR(DoxCM.cellchat, signaling = pathway.show, geneLR.return = F)
LR.show <- pair[1,]
pdf("./Output/Cell&Cell/IGF_pathway_top1_pair.pdf", width = 4, height = 4)
netVisual_individual(DoxCM.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "circle")
dev.off()
#netVisual_individual(DoxCM.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "chord")


levels(DoxCM.cellchat@idents)
pathway.show <- group.net$pathway_name
pdf("./Output/Cell&Cell/IGF_pathway_top1_pair.pdf", width = 4, height = 4)
netVisual_bubble(DoxCM.cellchat,
                 #signaling = pathway.show,
                 sources.use = "Fib_4",
                 remove.isolate = F,
                 font.size = 14, 
                 title.name = "Fib_4 vs others",
                 font.size.title = 20) + coord_flip()# + RotatedAxis()
dev.off()

#ggsave("IGF (20x20).png", width = 20, height = 40, path = "./DoxCM Output/FBs/Revision/")

pdf("./Output/Cell&Cell/CM_to_Fib_LR_pairs.pdf", width = 7, height = 9)
netVisual_bubble(DoxCM.cellchat,angle.x = 45,
                 sources.use = c("CM_1","CM_2","CM_3","CM_4"),
                 targets.use = c("Fib_4","Fib_3","Fib_2","Fib_1"),
                 #signaling = c("THBS","FN1","IGF","TGFb","PDGF","PERIOSTIN"),
                 remove.isolate = F,
                 font.size = 14) #+ coord_flip() + RotatedAxis()
dev.off()

pdf("./Output/Cell&Cell/IGF_pathway_gene_expression.pdf", width = 3, height = 3)
plotGeneExpression(DoxCM.cellchat, signaling = "IGF")
dev.off()

colors <- brewer.pal(14,"Oranges")
pdf("./Output/Cell&Cell/IGF_pathway_dotplot.pdf", width = 4, height = 2.5)
plotGeneExpression(DoxCM.cellchat, signaling = "IGF", type = "dot", col = colors)
dev.off()

DoxCM.cellchat <- netAnalysis_computeCentrality(DoxCM.cellchat, slot.name = "netP")

pdf("./Output/Cell&Cell/PERIOSTIN_pathway_influencer.pdf", width = 7, height = 5)
netAnalysis_signalingRole_network(DoxCM.cellchat, signaling = "PERIOSTIN",
                                  width = 15, height = 6, font.size = 10)
dev.off()
#ggsave("IGF influencer.png", width = 15, height = 6, path = "./DoxCM Output/FBs/Revision/")

pdf("./Output/Cell&Cell/Overall_scatter.pdf", width = 4, height = 3)
netAnalysis_signalingRole_scatter(DoxCM.cellchat, do.label = T)
dev.off()

netAnalysis_signalingRole_scatter(DoxCM.cellchat, signaling = c("COLLAGEN","NRXN","IGF"))
pdf("./Output/Cell&Cell/All_pathways_outgoing_heatmap.pdf", width = 6, height = 12)
netAnalysis_signalingRole_heatmap(DoxCM.cellchat, pattern = "outgoing", width = 6, height = 12)
dev.off()
pdf("./Output/Cell&Cell/All_pathways_incoming_heatmap.pdf", width = 6, height = 12)
netAnalysis_signalingRole_heatmap(DoxCM.cellchat, pattern = "incoming", width = 6, height = 12)
dev.off()

pdf("./Output/Cell&Cell/Selected_pathways_outgoing_heatmap.pdf", width = 7, height = 6)
netAnalysis_signalingRole_heatmap(DoxCM.cellchat, signaling = c("THBS","FN1","TENASCIN","NRXN","IGF"), pattern = "outgoing")
dev.off()
pdf("./Output/Cell&Cell/Selected_pathways_incoming_heatmap.pdf", width = 7, height = 6)
netAnalysis_signalingRole_heatmap(DoxCM.cellchat, signaling = c("THBS","FN1","TENASCIN","NRXN","IGF"), pattern = "incoming")
dev.off()

library(NMF)
library(ggalluvial)

selectK(DoxCM.cellchat, pattern = "outgoing")
nPatterns = 3
dev.off()
DoxCM.cellchat <- identifyCommunicationPatterns(DoxCM.cellchat, pattern = "outgoing",
                                                k = nPatterns, width = 7, height = 15, font.size = 6)
dev.off()
pdf("./Output/Cell&Cell/All_pathways_outgoing_patterns.pdf", width = 8, height = 6)
netAnalysis_river(DoxCM.cellchat, pattern = "outgoing")# + coord_flip() + RotatedAxis()
dev.off()


netAnalysis_dot(DoxCM.cellchat, pattern = "outgoing")

selectK(DoxCM.cellchat, pattern = "incoming")
nPatterns = 4
dev.off()

pdf("./Output/Cell&Cell/Overall_outgoing_CommunicationPatterns.pdf", width = 6, height = 8)
DoxCM.cellchat <- identifyCommunicationPatterns(DoxCM.cellchat, pattern = "outgoing",
                                                k = nPatterns, width = 7, height = 15, font.size = 6)
dev.off()

netAnalysis_river(DoxCM.cellchat, pattern = "incoming")
netAnalysis_dot(DoxCM.cellchat, pattern = "incoming")



Donor.cellchat <- readRDS("./Output/Cell&Cell/Donor_Cellchat(CM&FB).rds")
DoxCM.cellchat <- readRDS("./Output/Cell&Cell/DoxCM_Cellchat(CM&FB).rds")
NICM.cellchat <- readRDS("./Output/Cell&Cell/NICM_Cellchat(CM&FB).rds")


object.list <- list(Donor = Donor.cellchat, DoxCM = DoxCM.cellchat)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

pdf("./Output/Cell&Cell/Fib_vs_CM_comparisoin_numbers (2x3).pdf", width = 2, height = 3)
compareInteractions(cellchat, show.legend = F, group = c(1,2))
dev.off()
pdf("./Output/Cell&Cell/Fib_vs_CM_comparisoin_strengths (2x3).pdf", width = 2, height = 3)
compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weigh")
dev.off()

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weigh")
p <- gg1+gg2
p
ggsave(p,"Fib_vs_CM_comparisoin (4.5x5.5).pdf", width = 4.5, height = 5.5, path = "./Output/Cell&Cell/")


par(mfrow = c(1,2), xpd = T)
pdf("./Output/Cell&Cell/Fib_vs_CM_interactions_number.pdf", width = 5, height = 4)
netVisual_diffInteraction(cellchat, weight.scale = T)
dev.off()
pdf("./Output/Cell&Cell/Fib_vs_CM_interactions_weight.pdf", width = 5, height = 4)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

par(mfrow = c(1,1), xpd = T)
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2

cc1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T)
cc2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = T)
cc1+cc2


DoxCM.cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("./Output/Cell&Cell/Donor_overall_scatter.pdf", width = 4, height = 2.2)
netAnalysis_signalingRole_scatter(DoxCM.cellchat, do.label = T)
dev.off()



pdf("./Output/Cell&Cell/Fib_vs_CM_information_flow.pdf", width = 4, height = 6)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T, color.use = c("lightblue","red"))
dev.off()
pdf("./Output/Cell&Cell/Fib_vs_CM_information_flow_1.pdf", width = 4, height = 6)
rankNet(cellchat, mode = "comparison", stacked = F, do.stat = T, color.use = c("lightblue","red"))
dev.off()

diff.count <- cellchat@net$DoxCM$count - cellchat@net$Donor$count
write.csv(cellchat@net$DoxCM$count, "./DoxCM Output/FBs/Revision/output_DoxCMcounts.rds")
write.csv(cellchat@net$Donor$count, "./DoxCM Output/FBs/Revision/output_Donorcounts.rds")

library(pheatmap)
pheatmap(diff.count,
         treeheight_row = "0", treeheight_col = "0",
         cluster_rows = T, cluster_cols = T)

View(cellchat)
