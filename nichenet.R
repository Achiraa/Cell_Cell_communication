#setting the path for the directory
setwd("C:/Users/S236282/Desktop/Itaconate data/NicheNet")

#Loading necessary libraries
library(nichenetr) 
library(Seurat)
library(SeuratObject)
library(tidyverse)

#reading the .rds file
seuratObj = readRDS("C:/Users/S236282/Desktop/Itaconate data/NicheNet/seurat_annotation.rds")
seuratObj@meta.data %>% head()

set.seed(1234)

#seuratObj@meta.data$type %>% table() 
#DimPlot(seuratObj, reduction = "umap")

#seuratObj@meta.data$batch %>% table()
#DimPlot(seuratObj, reduction = "umap", group.by = "batch")

#Loading the organism of interest
#ligand-receptor network contains information on potential ligand-receptor bindings, and it is used to identify potential ligands
#the weighted ligand-receptor network contains weights representing the potential that a ligand will bind to a receptor, and it is used for visualization.
organism <- "mouse"
if(organism == "human"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}

#matrix describing the potential that a ligand may regulate a target gene, and it is used to run the ligand activity analysis
lr_network <- lr_network %>% distinct(from, to)

######################## NicheNet Analysis ################################

#Keeping the receiver Neutrophils as they are of cell of interest
#The receiver cell population can only consist of one cell type,
receiver="Neu"

#getting the genes expressed in at least 5% Neutophils cells
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)

# getting the list of all the receptors available in the network
all_receptors <- unique(lr_network$to)

#define expressed receptors as genes that are in the ligand-receptor network and expressed in the receiver.
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

#define the potential ligands as all ligands whose cognate receptors are expressed.
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

#defining the sender population
sender_celltypes <- c("MPC", "SMC", "EC", "Schwann", "Epi", "Mac", "T cell 1", "T cell 2", "Mast")

#Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

#filter potential ligands to those that are expressed in sender cells
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

#Also checking the length
length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)

######################### Define the gene set of interest #########################

#defining the test and reference state
condition_oi <-  "HO"
condition_reference <- "Blood"

#subsetting Neutrophils
seurat_obj_receiver <- subset(seuratObj, idents = receiver)

#extracting the DE genes within neutrophils
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "batch",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Defining the background genes i.e, all genes in the neutrophils
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
length(background_expressed_genes)


##################### Perform NicheNet ligand activity analysis ########################

#ain step of NicheNet where the potential ligands are ranked based on the presence 
#of their target genes in the gene set of interest (compared to the background set of genes)
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

#Ligands are ranked based on the area under the precision-recall curve (AUPR) between a 
#ligand's target predictions and the observed transcriptional response.
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

best_upstream_ligands <- ligand_activities %>% top_n(50, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

#visualize the ligand activity measure (AUPR) of these top-ranked ligands:
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

#Active target genes are defined as genes in the gene set of interest that 
#have the highest regulatory potential for each top-ranked ligand.
#The function get_weighted_ligand_target_links will return genes that are in the gene set of interest and are the top `n` targets of a ligand
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>%
  bind_rows() %>% drop_na()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

#we identify which receptors have the highest interaction potential with the top-ranked ligands.
#Then, we create a heatmap for ligand-receptor interactions.

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 
nrow(active_ligand_target_links)
head(active_ligand_target_links)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")


ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig)

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))



############################# Sender focused approach ##################################

#simply subset the ligand activities to only contain expressed ligands from all populations
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()


#We can then perform target gene and receptor inference as above.

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr

# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target


# Receptor plot
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor

best_upstream_ligands_all %in% rownames(seuratObj) %>% table()

### Visualizing expression and log-fold change in sender cells

# Dotplot of sender-focused approach
p_dotplot <- DotPlot(subset(seuratObj, celltype %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot

celltype_order <- levels(Idents(seuratObj)) 
celltype_order

# Use this if cell type labels are the identities of your Seurat object
# if not: indicate the celltype_col properly
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = seuratObj,
  condition_colname = "batch",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  celltype_col = "type",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE]) 

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")

p_lfc

(make_line_plot(ligand_activities = ligand_activities_all,
                potential_ligands = potential_ligands_focused) +
    theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))


## 7. Summary visualizations of the NicheNet analysis

figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  p_lfc + theme(legend.position = "none",
                axis.title.y = element_blank()),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc)+7, ncol(vis_ligand_lfc)+8, ncol(vis_ligand_target)))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot