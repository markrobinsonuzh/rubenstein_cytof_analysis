# 
# knitr::spin("cytofWorkflow_Cancer_20clusts_merge3_PDFs.R")
## ----setup_knitr, include = FALSE, cache = FALSE-------------------------
# Decide whether to display parts for BioC (TRUE) or F1000 (FALSE)
on.bioc <- TRUE
# Use fig.width = 7 for html and fig.width = 6 for pdf
fig.width <- ifelse(on.bioc, 7, 6)

library(knitr)
# Use fig.width = 7 for html and fig.width = 6 for pdf
knitr::opts_chunk$set(cache = 2, warning = FALSE, message = FALSE, error = FALSE,
cache.path = "cache/", fig.path = "figure/", fig.width = 15, fig.height=12)


## ----libraries, echo = FALSE, results = "hide"---------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(flowCore)
  library(matrixStats)
  library(ggplot2)
  library(ggridges)
  library(reshape2)
  library(dplyr)
  library(limma)
  library(ggrepel)
  library(RColorBrewer)
  library(pheatmap)
  library(ComplexHeatmap)
  library(FlowSOM)
  library(ConsensusClusterPlus)
  library(Rtsne)
  library(cowplot)
  library(lme4)
  library(multcomp)
  library(cytofWorkflow)
})


## ----load-metadata-------------------------------------------------------
library(readxl)
metadata_filename <- "metadata_CK_12-18-2017.xlsx"
md <- read_excel(metadata_filename)

## Make sure condition variables are factors with the right levels
md$condition <- factor(md$condition, levels = c("healthy", "LC"))
data.frame(md)

dir("002_fcsclean")

file.exists(file.path("002_fcsclean",md$file_name))


## Define colors for conditions
color_conditions <- c("#6A3D9A", "#FF7F00")
names(color_conditions) <- levels(md$condition)


## ----load-fcs, include=FALSE---------------------------------------------
library(flowCore)
fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, 
  truncate_max_range = FALSE, path="002_fcsclean")
fcs_raw


## ----load-panel----------------------------------------------------------
panel_filename <- "panel_CK_12-18-2017.xlsx"
panel <- read_excel(panel_filename)
head(data.frame(panel))
# Replace problematic characters 
panel$Antigen <- gsub("-", "_", panel$Antigen)

panel_fcs <- pData(parameters(fcs_raw[[1]]))
head(panel_fcs)
# Replace problematic characters 
panel_fcs$desc <- unname(gsub("-", "_", panel_fcs$desc))
panel_fcs$desc <- sapply(strsplit(panel_fcs$desc,"_"), function(u) paste0(u[-1],collapse="_"))
# Lineage markers
(lineage_markers <- panel$Antigen[panel$Transform == 1])

lineage_markers <- setdiff(lineage_markers, "CD45")

# Functional markers
#(functional_markers <- panel$Antigen[panel$Functional == 1])

#cbind(lineage_markers, lineage_markers %in% panel_fcs$desc)

# Spot checks
all(lineage_markers %in% panel_fcs$desc)
#all(functional_markers %in% panel_fcs$desc)


## ----arcsinh-transformation----------------------------------------------
## arcsinh transformation and column subsetting
fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  #expr <- asinh(expr[, c(lineage_markers, functional_markers)] / cofactor)
  expr <- asinh(expr[, c(lineage_markers)] / cofactor)
  exprs(x) <- expr
  x
})
fcs

## ----extract-expression--------------------------------------------------
## Extract expression
expr <- fsApply(fcs, exprs)
dim(expr)

## ----01-transformation---------------------------------------------------
library(matrixStats)
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

## ----sample-ids----------------------------------------------------------
## Generate sample IDs corresponding to each cell in the `expr` matrix
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))

## ----plot-merker-expression-distribution, fig.cap = "Per-sample smoothed densities of marker expression (arcsinh-transformed) of 10 lineage markers and 14 functional markers in the PBMC dataset. Two conditions: unstimulated (Ref) and stimulated with BCR/FcR-XL (BCRXL) for each of the 8 healthy donors are presented and colored by experimental condition."----
library(ggplot2)
library(reshape2)

ggdf <- data.frame(sample_id = sample_ids, expr)
ggdf <- melt(ggdf, id.var = "sample_id", 
  value.name = "expression", variable.name = "antigen")
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = expression, color = condition, 
  group = sample_id)) +
  geom_density() +
  facet_wrap(~ antigen, nrow = 4, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
    strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  scale_color_manual(values = color_conditions)


## ----plot-number-of-cells, fig.height = 4, fig.cap = "Barplot showing the number of cells measured for each sample in the PBMC dataset. Bars are colored by experimental condition: unstimulated (Ref) and stimulated with BCR/FcR-XL (BCRXL). Numbers in the names on the x-axis indicate patient IDs. Numbers on top of the bars indicate the cell counts."----
cell_table <- table(sample_ids)

ggdf <- data.frame(sample_id = names(cell_table), 
  cell_counts = as.numeric(cell_table))
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = condition)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = cell_counts), hjust=0.5, vjust=-0.5, size = 2.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color_conditions, drop = FALSE) +
  scale_x_discrete(drop = FALSE)


## ----plot-mds, fig.cap = "MDS plot for the unstimulated (Ref) and stimulated with BCR/FcR-XL (BCRXL) samples obtained for each of the 8 healthy donors in the PBMC dataset. Calculations are based on the median (arcsinh-transformed) marker expression of 10 lineage markers and 14 functional markers across all cells measured for each sample.  Distances between samples in the plot approximate the typical change in medians. Numbers in the label names indicate patient IDs."----
library(dplyr)
# Get the median marker expression per sample
expr_median_sample_tbl <- data.frame(sample_id = sample_ids, expr) %>%
  group_by(sample_id) %>% summarize_all(funs(median))

expr_median_sample <- t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id

library(limma)
mds <- plotMDS(expr_median_sample, plot = FALSE)

library(ggrepel)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y, 
  sample_id = colnames(expr_median_sample))
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = MDS1, y = MDS2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw() +
  scale_color_manual(values = color_conditions) +
  coord_fixed()


## ----plot-dendogram, fig.cap = "Heatmap of the median (arcsinh-transformed) marker expression of 10 lineage markers and 14 functional markers across all cells measured for each sample in the PBMC dataset. Color-coded with yellow for lower expression and blue for higher expression. The numbers in the heatmap represent the actual expression values. Dendrograms present clustering of samples (columns) and markers (rows) which is based on hierarchical clustering with Euclidean distance metric and average linkage. The two conditions: unstimulated (Ref) and stimulated with BCR/FcR-XL (BCRXL) for each of the 8 healthy donors are presented with a bar colored by experimental condition on top of the heatmap. Numbers in the column label names indicate patient IDs."----
library(RColorBrewer)
library(pheatmap)

# Column annotation for the heatmap
mm <- match(colnames(expr_median_sample), md$sample_id)
annotation_col <- data.frame(condition = md$condition[mm],
  row.names = colnames(expr_median_sample))
annotation_colors <- list(condition = color_conditions)

# Colors for the heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)

pheatmap(expr_median_sample, color = color, display_numbers = TRUE, 
  number_color = "black", fontsize_number = 5, annotation_col = annotation_col, 
  annotation_colors = annotation_colors, clustering_method = "average")


## ----nrs, fig.height = 4, fig.cap = "Non-redundancy scores for each of the 10 lineage markers and all samples in the PBMC dataset. The full points represent the per-sample NR scores (colored by experimental conditions), while empty black circles indicate the mean NR scores from all the samples. Markers on the x-axis are sorted according to the decreasing average NRS."----
## Define a function that calculates the NRS per sample 
NRS <- function(x, ncomp = 3){
  pr <- prcomp(x, center = TRUE, scale. = FALSE) 
  score <- rowSums(outer(rep(1, ncol(x)), 
    pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
  return(score)
}

## Calculate the score
nrs_sample <- fsApply(fcs[, lineage_markers], NRS, use.exprs = TRUE)
rownames(nrs_sample) <- md$sample_id
nrs <- colMeans(nrs_sample, na.rm = TRUE)

## Plot the NRS for ordered markers
lineage_markers_ord <- names(sort(nrs, decreasing = TRUE))
nrs_sample <- data.frame(nrs_sample)
nrs_sample$sample_id <- rownames(nrs_sample)

ggdf <- melt(nrs_sample, id.var = "sample_id", 
  value.name = "nrs", variable.name = "antigen")

ggdf$antigen <- factor(ggdf$antigen, levels = lineage_markers_ord)
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = antigen, y = nrs)) +
  geom_point(aes(color = condition), alpha = 0.9, 
    position = position_jitter(width = 0.3, height = 0)) +
  geom_boxplot(outlier.color = NA, fill = NA) +
  stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = color_conditions)

lineage_markers <- lineage_markers_ord[1:12]


library(FlowSOM)
fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)
## Get the cell clustering into 100 SOM codes
cell_clustering_som <- som$map$mapping[,1]

## ----flowsom-meta-clustering, message = FALSE----------------------------
## Metaclustering into 20 clusters with ConsensusClusterPlus
library(ConsensusClusterPlus)

codes <- som$map$codes
plot_outdir <- "consensus_plots"
nmc <- 20

mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100, 
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png", 
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", 
                           distance = "euclidean", seed = 1234)

## Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[cell_clustering_som]


## ----color-clusters------------------------------------------------------
# Define cluster colors (here there are 30 colors)
color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", 
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

## ----plot-clustering-heatmap1, fig.cap = "Heatmap of the median marker intensities of the 10 lineage markers across the 20 cell populations obtained with FlowSOM after the metaclustering step with ConsensusClusterPlus (PBMC data). The color in the heatmap represents the median of the arcsinh, 0-1 transformed marker expression calculated over cells from all the samples, varying from blue for lower expression to red for higher expression. The numbers indicate the actual expression values. The dendrogram on the left represents the hierarchical similarity between the 20 metaclusters (metric: Euclidean distance; linkage: average). Each cluster has a unique color assigned (bar on the left) which is identical in other visualizations of these 20 clusters (e.g. Figure t-SNE map in \\@ref(fig:tsne-plot-one-clustering1)) facilitating the figure interpretation. Values in the brackets next to the cluster numbers indicate the relative size of clusters."----

plot_clustering_heatmap_wrapper <- function(expr, expr01, 
                                            cell_clustering, color_clusters, cluster_merging = NULL){
  
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  legend_breaks = seq(from = 0, to = max(expr01_median), length.out = 6)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop , 
                       "%)")
  
  # Annotation for the original clusters
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  # Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Cluster_merging <- cluster_merging$new_cluster 
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Cluster_merging <- color_clusters2
  }
  
  pheatmap(expr_heat, color = color_heat, cluster_cols = FALSE, 
           cluster_rows = cluster_rows, labels_row = labels_row, 
           display_numbers = TRUE, number_color = "black", 
           fontsize = 8, fontsize_number = 6, legend_breaks = legend_breaks,
           annotation_row = annotation_row, annotation_colors = annotation_colors)
  
}

plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers], 
                                expr01 = expr01[, lineage_markers], 
                                cell_clustering = cell_clustering1, color_clusters = color_clusters)

plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord], 
                                expr01 = expr01[, lineage_markers_ord], 
                                cell_clustering = cell_clustering1, color_clusters = color_clusters)

## ----plot-clustering-distribution1, fig.height = 8, fig.cap = "Distributions of marker intensities (arcsinh-transformed) of the 10 lineage markers in the 20 cell populations obtained with FlowSOM after the metaclustering step with ConsensusClusterPlus (PBMC data). Red densities represent marker expression for cells in a given cluster. Green densities are calculated over all the cells and serve as a reference."----

library(ggridges)

plot_clustering_distr_wrapper <- function(expr, cell_clustering){
  # Calculate the median expression
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  # Calculate cluster frequencies
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering, 
                            labels = paste0(levels(cell_clustering), " (", freq_clust, "%)"))
  ### Data organized per cluster
  ggd <- melt(data.frame(cluster = cell_clustering, expr), 
              id.vars = "cluster", value.name = "expression", 
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  ### The reference data
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster, 
                             levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
  
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster, 
                                             color = reference, fill = reference), alpha = 0.3) +
    facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
    theme_ridges() +
    theme(axis.text = element_text(size = 7), 
          strip.text = element_text(size = 7), legend.position = "none") 
  
}

plot_clustering_distr_wrapper(expr = expr[, lineage_markers], 
                              cell_clustering = cell_clustering1)



library(ComplexHeatmap)

plot_clustering_heatmap_wrapper2 <- function(expr, expr01,
                                             lineage_markers, functional_markers = NULL, sample_ids = NULL,
                                             cell_clustering, color_clusters, cluster_merging = NULL, 
                                             plot_cluster_annotation = TRUE){
  
  # Calculate the median expression of lineage markers
  expr_median <- data.frame(expr[, lineage_markers], 
                            cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%  summarize_all(funs(median))
  expr01_median <- data.frame(expr01[, lineage_markers], 
                              cell_clustering = cell_clustering) %>% 
    group_by(cell_clustering) %>%  summarize_all(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, lineage_markers], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, lineage_markers])
  
  # Median expression of functional markers in each sample per cluster
  expr_median_sample_cluster_tbl <- data.frame(expr01[, functional_markers, 
                                                      drop = FALSE], sample_id = sample_ids, cluster = cell_clustering) %>%
    group_by(sample_id, cluster) %>% summarize_all(funs(median))
  
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  legend_breaks = round(seq(from = 0, to = max(expr_heat), length.out = 6),2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop , "%)")
  
  ### Annotation for the original clusters
  annotation_row1 <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  color_clusters1 <- color_clusters[1:nlevels(annotation_row1$Cluster)]
  names(color_clusters1) <- levels(annotation_row1$Cluster)
  
  ### Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    mm <- match(annotation_row1$Cluster, cluster_merging$original_cluster)
    annotation_row2 <- data.frame(Cluster_merging = 
                                    factor(cluster_merging$new_cluster[mm]))
    color_clusters2 <- color_clusters[1:nlevels(annotation_row2$Cluster_merging)]
    names(color_clusters2) <- levels(annotation_row2$Cluster_merging)
  }
  
  ### Heatmap annotation for the original clusters
  ha1 <- Heatmap(annotation_row1, name = "Cluster", 
                 col = color_clusters1, cluster_columns = FALSE, 
                 cluster_rows = cluster_rows, row_dend_reorder = FALSE, 
                 show_row_names = FALSE, width = unit(0.5, "cm"), 
                 rect_gp = gpar(col = "grey"))
  ### Heatmap annotation for the merged clusters
  if(!is.null(cluster_merging)){
    ha2 <- Heatmap(annotation_row2, name = "Cluster \nmerging", 
                   col = color_clusters2, cluster_columns = FALSE, 
                   cluster_rows = cluster_rows, row_dend_reorder = FALSE, 
                   show_row_names = FALSE, width = unit(0.5, "cm"), 
                   rect_gp = gpar(col = "grey"))
  }
  ### Cluster names and sizes - text
  ha_text <- rowAnnotation(text = row_anno_text(labels_row, 
                                                gp = gpar(fontsize = 6)), width = max_text_width(labels_row))
  ### Cluster sizes - barplot
  ha_bar <- rowAnnotation("Frequency (%)" = row_anno_barplot(
    x = clustering_prop, border = FALSE, axis = TRUE, 
    axis_gp = gpar(fontsize = 5), gp = gpar(fill = "#696969", col = "#696969"), 
    bar_width = 0.9), width = unit(0.7, "cm"), show_annotation_name = TRUE, 
    annotation_name_rot = 0, annotation_name_offset = unit(5, "mm"), 
    annotation_name_gp = gpar(fontsize = 7))
  ### Heatmap for the lineage markers
  ht1 <- Heatmap(expr_heat, name = "Expr",  column_title = "Lineage markers", 
                 col = color_heat, cluster_columns = FALSE, cluster_rows = cluster_rows, 
                 row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, 
                                                                       labels = legend_breaks, color_bar = "continuous"), 
                 show_row_names = FALSE, row_dend_width = unit(2, "cm"), 
                 rect_gp = gpar(col = "grey"), column_names_gp = gpar(fontsize = 8))
  
  if(plot_cluster_annotation){
    draw_out <- ha1
  }else{
    draw_out <- NULL  
  }
  if(!is.null(cluster_merging)){
    draw_out <- draw_out + ha2 + ht1 + ha_bar + ha_text
  }else{
    draw_out <- draw_out + ht1 + ha_bar + ha_text
  }
  
  ### Heatmaps for the signaling markers
  if(!is.null(functional_markers)){
    for(i in 1:length(functional_markers)){
      ## Rearange so the rows represent clusters
      expr_heat_fun <- as.matrix(dcast(expr_median_sample_cluster_tbl[, 
                                                                      c("sample_id", "cluster", functional_markers[i])], 
                                       cluster ~ sample_id, value.var = functional_markers[i])[, -1])
      
      draw_out <- draw_out + Heatmap(expr_heat_fun, 
                                     column_title = functional_markers[i], col = color_heat, 
                                     cluster_columns = FALSE, cluster_rows = cluster_rows, 
                                     row_dend_reorder = FALSE, show_heatmap_legend = FALSE, 
                                     show_row_names = FALSE, rect_gp = gpar(col = "grey"), 
                                     column_names_gp = gpar(fontsize = 8))
    }
  }
  draw(draw_out, row_dend_side = "left")
}

#plot_clustering_heatmap_wrapper2(expr = expr, expr01 = expr01,
#                                 lineage_markers = lineage_markers, functional_markers = "pS6", 
#                                 sample_ids = sample_ids, cell_clustering = cell_clustering1, 
#                                 color_clusters = color_clusters, cluster_merging = NULL)


dups <- which(!duplicated(expr[, lineage_markers]))

## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids) 

## How many cells to downsample per-sample
tsne_ncells <- pmin(table(sample_ids), 5000)  

## Get subsampled indices
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})

tsne_inds <- unlist(tsne_inds)

tsne_expr <- expr[tsne_inds, lineage_markers]

## ----tsne-run------------------------------------------------------------
## Run t-SNE
library(Rtsne)

set.seed(1234)
tsne_out <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE, verbose = TRUE)


## ----tsne-plot-one-expr-CD4, fig.cap = "t-SNE plot based on the arcsinh-transformed expression of the 10 lineage markers in the cells from the PBMC dataset. t-SNE was run with no PCA step, perplexity equal to 30 and 1000 iterations. From each of the 16 samples, 2000 cells were randomly selected. Cells are colored according to the expression level of the CD4 marker."----
## Plot t-SNE colored by CD4 expression
dr <- data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2], 
                 expr[tsne_inds, lineage_markers_ord])

ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = CD4)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))
ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = CD8)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))
ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = HLA_DR)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))


dr$sample_id <- sample_ids[tsne_inds]
mm <- match(dr$sample_id, md$sample_id)
dr$condition <- md$condition[mm]
dr$condition_time <- md$time[mm]
dr$patient <- md$patient[mm]
dr$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)

## Plot t-SNE colored by clusters
ggp <- ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering1)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggp


## ----tsne-plot-facet-sample, fig.height = 12----
## Facet per sample
ggp + facet_wrap(~ sample_id) 

## ----tsne-plot-facet-condition, fig.height = 12----
## Facet per condition
ggp + facet_wrap(~ condition) 

## ----tsne-plot-facet-time, fig.height = 12----
## Facet per condition
ggp + facet_wrap(~ condition_time) 


## ----tsne-plot-facet-patient-time, fig.height = 12----
## Facet per patient/time
ggp + facet_wrap(patient ~ condition_time) 


## ----cluster-merging1----------------------------------------------------
cluster_merging1_filename <- "CK_2017-12-18merge3.xlsx"
cluster_merging1 <- read_excel(cluster_merging1_filename)
data.frame(cluster_merging1)
## Convert to factor with merged clusters in desired order
#levels_clusters_merged <- c("B-cells IgM+", "B-cells IgM-", "CD4 T-cells", 
#                            "CD8 T-cells", "DC", "NK cells", "monocytes", "surface-")
levels_clusters_merged <- unique(cluster_merging1$new_cluster)
cluster_merging1$new_cluster <- factor(cluster_merging1$new_cluster, 
                                       levels = levels_clusters_merged)
## New clustering1m
mm <- match(cell_clustering1, cluster_merging1$original_cluster)
cell_clustering1m <- cluster_merging1$new_cluster[mm]

## export Tregs into FCS
ncells <- fsApply(fcs_raw, nrow)[,1]
inds_subset <- lapply(ncells, function(u) seq_len(u))
ks <- split(cell_clustering1m=="Tregs", rep(md$file_name, ncells))
sapply(ks, sum)

# check order of everything
all(names(ks)==names(inds_subset))
all( md$file_name == names(ks) )

fcs_tregs <- fcs_raw
for(i in 1:length(md$file_name))
  fcs_tregs[[i]] <- fcs_tregs[[i]][ks[[i]]]

tregs_fn <- gsub(".fcs","_tregs.fcs", md$file_name)
write.flowSet(fcs_tregs, outdir="003_out_fcs", filename=tregs_fn)



## ----tsne-plot-one-clustering1m, fig.cap = "t-SNE plot for the PBMC dataset, where cells are colored according to the manual merging of the 20 cell populations, obtained with FlowSOM, into 8 PBMC populations. As in Figure \\@ref(fig:tsne-plot-one-clustering1), t-SNE analysis uses the arcsinh-transformed expression of the 10 lineage markers in 2000 randomly selected cells from each of the 16 samples."----
dr$cell_clustering1m <- cell_clustering1m[tsne_inds]
ggp <- ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering1m)) +
  geom_point(size = 0.3) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggp

ggsave(filename = "tSNE_panel.pdf", plot = ggp, width = 8, height = 6)

## ----tsne-plot-facet-merge-patient-time-1, fig.height = 12----
## Facet per patient/time
ggp + facet_wrap(patient ~ condition_time) 

## ----plot-clustering-heatmap1m-merging, fig.height = 12, fig.cap = "Heatmap as in Figure \\@ref(fig:plot-clustering-heatmap1), where the additional color bar on the left indicates how the 20 metaclusters, obtained with FlowSOM, are merged into the 8 PBMC populations."----
pchw <- plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
                                expr01 = expr01[, lineage_markers_ord], cell_clustering = cell_clustering1,
                                color_clusters = color_clusters, cluster_merging = cluster_merging1)

ggsave(filename = "heatmap_panel_01.pdf", pchw$gtable, width = 10, height = 6)


## ----plot-clustering-heatmap1m-merging-originalscale, fig.height = 12, fig.cap = "Heatmap as in Figure \\@ref(fig:plot-clustering-heatmap1), where the additional color bar on the left indicates how the 20 metaclusters, obtained with FlowSOM, are merged into the 8 PBMC populations."----
pchw <- plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
                                        expr01 = expr[, lineage_markers_ord], cell_clustering = cell_clustering1,
                                        color_clusters = color_clusters, cluster_merging = cluster_merging1)

ggsave(filename = "heatmap_panel_original.pdf", pchw$gtable, width = 10, height = 6)



o <- pchw$tree_row$order
cell_counts <- table(cell_clustering1, sample_ids)
props <- t(t(cell_counts) / colSums(cell_counts))

ph <- pheatmap(props[o,], color = color, display_numbers = TRUE, scale = "row",
         number_color = "black", fontsize_number = 5, annotation_col = annotation_col, 
         annotation_colors = annotation_colors, cluster_rows = FALSE, cluster_cols = FALSE)

ph <- pheatmap(props[o,], color = color, display_numbers = TRUE,
               number_color = "black", fontsize_number = 5, annotation_col = annotation_col, 
               annotation_colors = annotation_colors, cluster_rows = FALSE, cluster_cols = FALSE)

ph <- pheatmap(sqrt(props[o,]), color = color, display_numbers = TRUE,
               number_color = "black", fontsize_number = 5, annotation_col = annotation_col, 
               annotation_colors = annotation_colors, cluster_rows = FALSE, cluster_cols = FALSE)

ggsave(filename = "cluster_freq_sqrt.pdf", ph$gtable, width=6, height=8)



## ----plot-clustering-heatmap1m, fig.height = 12, fig.cap = "Heatmap of the median marker intensities of the 10 lineage markers in the 8 PBMC cell populations obtained by manual merging of the 20 metaclusters generated by FlowSOM. As in Figure \\@ref(fig:plot-clustering-heatmap1), the heat represents the median of arcsinh and 0-1 transformed marker expression calculated over cells from all the samples. The dendrogram on the left represents the hierarchical similarity between the 8 populations calculated using Euclidean distance and average linkage. Values in the brackets indicate the relative size of each of the cell populations across all the samples."----
# plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
#                                 expr01 = expr01[, lineage_markers_ord], cell_clustering = cell_clustering1m,
#                                 color_clusters = color_clusters)


lineage_markers_ord


sample_ids1 <- sample_ids
sample_ids1[grep("Healthy", sample_ids1)] <- "healthy"
## ----plot-clustering-heatmap2-1, fig.height = 12
plot_clustering_heatmap_wrapper2(expr = expr, expr01 = expr01,
                                 lineage_markers = lineage_markers, functional_markers = "HLA_DR", 
                                 sample_ids = sample_ids1, cell_clustering = cell_clustering1m, 
                                 color_clusters = color_clusters, cluster_merging = NULL)

## ----plot-clustering-heatmap2-1-unnormalized, fig.height = 12
pdf("heatmap_lineage+HLADR_original.pdf", width = 10, height = 5)
pchw2 <- plot_clustering_heatmap_wrapper2(expr = expr, expr01 = expr,
                                 lineage_markers = lineage_markers, functional_markers = "HLA_DR", 
                                 sample_ids = sample_ids1, cell_clustering = cell_clustering1m, 
                                 color_clusters = color_clusters, cluster_merging = NULL)
dev.off()
pchw2

pdf("heatmap_lineage+HLADR_original_clusters3,4,13.pdf", width = 10, height = 3)
k <- cell_clustering1 %in% c(3,4,13)
pchw2 <- plot_clustering_heatmap_wrapper2(expr = expr[k,], expr01 = expr[k,],
                                          lineage_markers = lineage_markers, functional_markers = "HLA_DR", 
                                          sample_ids = sample_ids1[k], cell_clustering = cell_clustering1[k], 
                                          color_clusters = color_clusters[c(3,4,13)], cluster_merging = NULL)
dev.off()
pchw2


## ----plot-clustering-heatmap2-2, fig.height = 12
plot_clustering_heatmap_wrapper2(expr = expr, expr01 = expr01,
                                 lineage_markers = lineage_markers, functional_markers = "Grz_B", 
                                 sample_ids = sample_ids1, cell_clustering = cell_clustering1m, 
                                 color_clusters = color_clusters, cluster_merging = NULL)

## ----plot-clustering-heatmap2-2-unnormalized, fig.height = 12
pdf("heatmap_lineage+GrzB_original.pdf", width=10, height=5)
pchw2 <- plot_clustering_heatmap_wrapper2(expr = expr, expr01 = expr,
                                 lineage_markers = lineage_markers, functional_markers = "Grz_B", 
                                 sample_ids = sample_ids1, cell_clustering = cell_clustering1m, 
                                 color_clusters = color_clusters, cluster_merging = NULL)
dev.off()
pchw2


## ----plot-clustering-heatmap2-3, fig.height = 12
plot_clustering_heatmap_wrapper2(expr = expr, expr01 = expr01,
                                 lineage_markers = lineage_markers, functional_markers = "ki_67", 
                                 sample_ids = sample_ids1, cell_clustering = cell_clustering1m, 
                                 color_clusters = color_clusters, cluster_merging = NULL)

## ----plot-clustering-heatmap2-3-unnormalized, fig.height = 12
pdf("heatmap_lineage+ki67_original.pdf", width=10, height=5)
pchw2 <- plot_clustering_heatmap_wrapper2(expr = expr, expr01 = expr,
                                 lineage_markers = lineage_markers, functional_markers = "ki_67", 
                                 sample_ids = sample_ids1, cell_clustering = cell_clustering1m, 
                                 color_clusters = color_clusters, cluster_merging = NULL)
dev.off()
pchw2



mm <- match(dr$cell_clustering1, cluster_merging1$original_cluster)
dr$cell_clustering1m <- dr$cell_clustering1[mm]

dr1 <- dr[dr$condition_time != "healthy",]

## Plot t-SNE colored by clusters
ggp <- ggplot(dr1,  aes(x = tSNE1, y = tSNE2, color = ki_67)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

## ----tsne-plot-facet-patient-time-1, fig.height = 12----
## Facet per patient/time
ggp + facet_wrap(patient ~ condition_time) 

ggsave(filename = "tSNE_split_ki67.pdf", 
       ggp + facet_wrap(patient ~ condition_time),
       width=12, height=8)



## Plot t-SNE colored by clusters
ggp <- ggplot(dr1,  aes(x = tSNE1, y = tSNE2, color = Grz_B)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

## ----tsne-plot-facet-patient-time-2, fig.height = 12----
## Facet per patient/time
ggp + facet_wrap(patient ~ condition_time) 

ggsave(filename = "tSNE_split_GrzB.pdf", 
       ggp + facet_wrap(patient ~ condition_time),
       width=12, height=8)




## Plot t-SNE colored by clusters
ggp <- ggplot(dr1,  aes(x = tSNE1, y = tSNE2, color = HLA_DR)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

## ----tsne-plot-facet-patient-time-3, fig.height = 12----
## Facet per patient/time
ggp + facet_wrap(patient ~ condition_time) 


ggsave(filename = "tSNE_split_HLADR.pdf", 
       ggp + facet_wrap(patient ~ condition_time),
       width=12, height=8)


## Plot t-SNE colored by clusters
ggp <- ggplot(dr1,  aes(x = tSNE1, y = tSNE2, color = CD38)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

## ----tsne-plot-facet-patient-time-3, fig.height = 12----
## Facet per patient/time
ggp + facet_wrap(patient ~ condition_time) 


ggsave(filename = "tSNE_split_CD38.pdf", 
       ggp + facet_wrap(patient ~ condition_time),
       width=12, height=8)

