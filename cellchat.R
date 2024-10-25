setwd("C:/Users/S236282/Desktop/Itaconate data/CellChat")

library(Seurat)
library(tidyverse)
library(ggplot2)
library(CellChat)
library(patchwork)
library(presto)
library(future)
library(NMF)
library(ggalluvial)
library(cowplot)
library(svglite)
library(stringr)
library(plyr)
options(stringsAsFactors = FALSE)

#reading the .rds file
neu = readRDS("C:/Users/S236282/Desktop/Itaconate data/CellChat/seurat_annotation.rds")
View(neu@meta.data)

#loading the data
data.input = GetAssayData(neu, assay = "SCT", layer = "data")
labels = neu@meta.data$type
meta <- data.frame(type = labels, row.names = colnames(data.input))

# Create a new column 'batch' in the meta data frame
meta$location = neu@meta.data$batch

#View(meta)
#unique(meta$type)

#print(dim(data.input))  
#print(dim(meta)) 

#create the cell chat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "type")

#setting up the mouse database
#CellChatDB <- CellChatDB.mouse
#showDatabaseCategory(CellChatDB)
#View(CellChatDB)
#unique(CellChatDB$interaction$annotation)

# Show the structure of the database
#dplyr::glimpse(CellChatDB$interaction)

################# Creating new database using CellTalkDB ###################

#txt to rds
db.lr <- read.table("mouse_lr_pair.txt", header = TRUE, sep = "\t")

genedata <- read.csv("mouse_gene_info.csv")

#modifying the colnames as they are incompatible with cell chat
colnames(db.lr) <- plyr::mapvalues(colnames(db.lr), from = c("ligand_gene_symbol","receptor_gene_symbol","lr_pair"), 
                                     to = c("ligand","receptor","interaction_name"), warn_missing = TRUE)

# Add pathway name
db.lr$pathway_name <- str_extract(db.lr$interaction_name, "^[^_]*")

# Create a new database by typing one of the following commands
# Use user-provided gene information
db.new <- updateCellChatDB(db = db.lr, gene_info = genedata)

#Using the manual database for specific cell-cell communication analysis
CellChatDB = db.new
cellchat@DB <- CellChatDB

######## Preprocessing the expression data for cell-cell communication analysis #########

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
print(cellchat@data.signaling)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Communication probability and infer cellular communication network

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 5)

#Extract the inferred cellular communication network as a data frame

df.net <- subsetCommunication(cellchat)
View(df.net)
write.csv(df.net,file="Ligand Receptor interactions.csv")

#Infer the cell-cell communication at a signaling pathway level

cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

####### Visualization according to the pathways ########################
pathways.show <- c("Tnf") 
vertex.receiver = seq(1,4) # Here we define `vertex.receive` so that the left portion of the hierarchy plot 
                           
#shows signaling to fibroblast and the right portion shows signaling to immune cells
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")


############### Visualization according to the receptor #################

tlr2_interactions <- df.net[df.net$receptor == "Tlr2" | df.net$ligand == "Tlr2", ]
print(tlr2_interactions)
