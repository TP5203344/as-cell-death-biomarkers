# step6_cellchat.R - Cell-cell communication analysis (CellChat)
# Supplementary Figure S10
# Data: GSE115469 (seurat_processed.rds, 8439 cells, 7 cell types)
# Environment: R 4.4.x, CellChat 1.6.1, Seurat 5.5.0

library(Seurat)
library(CellChat)

# 1. Load Seurat object and prepare CellChat input
seurat <- readRDS("seurat_processed.rds")
data.input <- GetAssayData(seurat, assay = "RNA", layer = "data")
meta <- data.frame(group = seurat$celltype, row.names = colnames(seurat))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

# 2. Set ligand-receptor interaction database (human)
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# 3. Preprocessing
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 4. Communication probability inference
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 5. Save object
saveRDS(cellchat, "cellchat_object.rds")

# 6. Visualization
# 6.1 Network overview (count + strength)
png("S10_cellchat_count.png", width = 1400, height = 1200, res = 200)
netVisual_circle(net = cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
dev.off()
png("S10_cellchat_strength.png", width = 1400, height = 1200, res = 200)
netVisual_circle(net = cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength")
dev.off()
# 6.2 Pathway bubble (top 25)
netP <- cellchat@netP$prob
path.strength <- apply(netP, 3, sum)
top25 <- names(sort(path.strength, decreasing = TRUE))[1:25]
png("S10_cellchat_pathway_bubble.png", width = 3200, height = 2600, res = 260)
netVisual_bubble(cellchat, signaling = top25, remove.isolate = TRUE)
dev.off()

# 7. Export communication pairs
write.csv(cellchat@net$count > 0, "cellchat_communication_pairs.csv")
write.csv(data.frame(pathway = names(path.strength), strength = as.numeric(path.strength)),
          "cellchat_pathway_strength.csv", row.names = FALSE)
