library(Seurat)
library(dplyr)
library(future)
options(future.globals.maxSize = 150 * 1024^3)
plan(sequential)
library(tidyverse)

ROOT <- "/rds/general/user/ts3225/ephemeral/R"

atlas_ref_dir <- file.path(ROOT, "atlas")

PDAC_primary <- readRDS(file = "/rds/general/user/ts3225/ephemeral/R/QC4/PDAC_primary.rds")
PDAC_met <- readRDS(file = "/rds/general/user/ts3225/ephemeral/R/QC4/PDAC_met.rds")


atlas_counts <- Read10X(
  data.dir = atlas_ref_dir,
  gene.column = 1,
  unique.features = TRUE
)
atlas_meta <- read.table(
  file.path(atlas_ref_dir, "meta.tsv"),
  header = TRUE,
  sep = "\t",
  as.is = TRUE,
  row.names = 1)

atlas_meta <- atlas_meta[colnames(atlas_counts), , drop = FALSE]

stopifnot(identical(rownames(atlas_meta), colnames(atlas_counts)))


atlas_ref <- CreateSeuratObject(
  counts = atlas_counts,
  project = "atlas",
  meta.data = atlas_meta)

DefaultAssay(atlas_ref) <- "RNA"

atlas_ref <- NormalizeData(atlas_ref)
atlas_ref <- FindVariableFeatures(atlas_ref)
atlas_ref <- ScaleData(atlas_ref)
atlas_ref <- RunPCA(atlas_ref)

DefaultAssay(PDAC_primary) <- "decontXcounts"

atlas_anchors_primary <- FindTransferAnchors(
  reference = atlas_ref,
  query = PDAC_primary,
  reference.assay = "RNA",
  query.assay = "decontXcounts",
  reference.reduction = "pca",
  dims = 1:30
)


primary_pred_l1 <- TransferData(
  anchorset = atlas_anchors_primary,
  refdata = atlas_ref$Level_1,
  dims = 1:30
)
primary_pred_l2 <- TransferData(
  anchorset = atlas_anchors_primary,
  refdata = atlas_ref$Level_2,
  dims = 1:30
)

PDAC_primary$atlas_prediction1 <- primary_pred_l1$predicted.id
PDAC_primary$atlas_score1 <- primary_pred_l1$prediction.score.max
PDAC_primary$atlas_prediction2 <- primary_pred_l2$predicted.id
PDAC_primary$atlas_score2 <- primary_pred_l2$prediction.score.max


DefaultAssay(PDAC_met) <- "decontXcounts"

atlas_anchors_met <- FindTransferAnchors(
  reference = atlas_ref,
  query = PDAC_met,
  reference.assay = "RNA",
  query.assay = "decontXcounts",
  reference.reduction = "pca",
  dims = 1:30
)


met_pred_l1 <- TransferData(
  anchorset = atlas_anchors_met,
  refdata = atlas_ref$Level_1,
  dims = 1:30
)
met_pred_l2 <- TransferData(
  anchorset = atlas_anchors_met,
  refdata = atlas_ref$Level_2,
  dims = 1:30
)

PDAC_met$atlas_prediction1 <- met_pred_l1$predicted.id
PDAC_met$atlas_score1 <- met_pred_l1$prediction.score.max
PDAC_met$atlas_prediction2 <- met_pred_l2$predicted.id
PDAC_met$atlas_score2 <- met_pred_l2$prediction.score.max

cluster_annotation_primary <- PDAC_primary@meta.data %>%
  dplyr::count(seurat_clusters,atlas_prediction1 ,atlas_prediction2) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::arrange(seurat_clusters, desc(frac))

write.csv(cluster_annotation_primary, file = "/rds/general/user/ts3225/ephemeral/R/QC4/cluster_annotation_primary.csv")


top_cluster_annotation_primary <- cluster_annotation_primary %>%
  slice_max(order_by = frac, n = 1, with_ties = FALSE)

write.csv(top_cluster_annotation_primary, file = "/rds/general/user/ts3225/ephemeral/R/QC4/top_cluster_annotation_primary.csv")


cluster_annotation_met <- PDAC_met@meta.data %>%
  dplyr::count(seurat_clusters,atlas_prediction1 ,atlas_prediction2) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::arrange(seurat_clusters, desc(frac))

write.csv(cluster_annotation_met, file = "/rds/general/user/ts3225/ephemeral/R/QC4/cluster_annotation_met.csv")


top_cluster_annotation_met <- cluster_annotation_met %>%
  slice_max(order_by = frac, n = 1, with_ties = FALSE)

write.csv(top_cluster_annotation_met, file = "/rds/general/user/ts3225/ephemeral/R/QC4/top_cluster_annotation_met.csv")


rm(atlas_counts, atlas_meta, atlas_ref,
   atlas_anchors_primary, atlas_anchors_met,
   primary_pred_l1, primary_pred_l2,
   met_pred_l1, met_pred_l2)
gc()



atlas_rds <- file.path(ROOT, "scAtlas.rds.gz")
outer_con <- gzfile(atlas_rds, "rb")
inner_con <- gzcon(outer_con)

atlas_ref <- readRDS(inner_con)

close(inner_con)

print(class(atlas_ref))
print(Assays(atlas_ref))
print(Reductions(atlas_ref))
print(colnames(atlas_ref@meta.data))


print(grep(
  "level|cell|type|annot|cluster|subtype",
  colnames(atlas_ref@meta.data),
  value = TRUE,
  ignore.case = TRUE
))

DefaultAssay(atlas_ref) <- "RNA"

stopifnot("Clusters" %in% colnames(atlas_ref@meta.data))

if (length(VariableFeatures(atlas_ref)) == 0) {
  atlas_ref <- FindVariableFeatures(atlas_ref, nfeatures = 3000)
}

if (!"pca" %in% Reductions(atlas_ref)) {
  atlas_ref <- NormalizeData(atlas_ref)
  atlas_ref <- ScaleData(atlas_ref)
  atlas_ref <- RunPCA(atlas_ref, npcs = 50)
}


dims_use <- 1:30

DefaultAssay(PDAC_primary) <- "decontXcounts"

atlas2_anchors_primary <- FindTransferAnchors(
  reference = atlas_ref,
  query = PDAC_primary,
  reference.assay = "RNA",
  query.assay = "decontXcounts",
  reference.reduction = "pca",
  reduction = "pcaproject",
  dims = dims_use
)



primary2_pred_clusters <- TransferData(
  anchorset = atlas2_anchors_primary,
  refdata = atlas_ref$Clusters,
  dims = dims_use
)



PDAC_primary$atlas2_clusters <- primary2_pred_clusters$predicted.id
PDAC_primary$atlas2_clusters_score <- primary2_pred_clusters$prediction.score.max


DefaultAssay(PDAC_met) <- "decontXcounts"

atlas2_anchors_met <- FindTransferAnchors(
  reference = atlas_ref,
  query = PDAC_met,
  reference.assay = "RNA",
  query.assay = "decontXcounts",
  reference.reduction = "pca",
  reduction = "pcaproject",
  dims = dims_use
)

met2_pred_clusters <- TransferData(
  anchorset = atlas2_anchors_met,
  refdata = atlas_ref$Clusters,
  dims = dims_use
)


PDAC_met$atlas2_clusters <- met2_pred_clusters$predicted.id
PDAC_met$atlas2_clusters_score <- met2_pred_clusters$prediction.score.max


cluster_annotation_primary2 <- PDAC_primary@meta.data %>%
  count(seurat_clusters, atlas2_clusters) %>%
  group_by(seurat_clusters) %>%
  mutate(frac = n / sum(n)) %>%
  arrange(seurat_clusters, desc(frac))

write.csv(
  cluster_annotation_primary2,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/cluster_annotation_primary2.csv",
  row.names = FALSE
)

top_cluster_annotation_primary2 <- cluster_annotation_primary2 %>%
  slice_max(order_by = frac, n = 1, with_ties = FALSE)


write.csv(
  top_cluster_annotation_primary2,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/top_cluster_annotation_primary2.csv",
  row.names = FALSE
)

cluster_annotation_met2 <- PDAC_met@meta.data %>%
  count(seurat_clusters, atlas2_clusters) %>%
  group_by(seurat_clusters) %>%
  mutate(frac = n / sum(n)) %>%
  arrange(seurat_clusters, desc(frac))

write.csv(
  cluster_annotation_met2,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/cluster_annotation_met2.csv",
  row.names = FALSE
)

top_cluster_annotation_met2 <- cluster_annotation_met2 %>%
  slice_max(order_by = frac, n = 1, with_ties = FALSE)

write.csv(
  top_cluster_annotation_met2,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/top_cluster_annotation_met2.csv",
  row.names = FALSE
)

rm(atlas_ref, atlas2_anchors_primary, atlas2_anchors_met,
   primary2_pred_clusters, met2_pred_clusters)
gc()


library(SingleR)
library(celldex)

blueref <- BlueprintEncodeData()
blue_cell <- SingleR(
  test = GetAssayData(PDAC_primary, assay = "decontXcounts", layer = "data"),
  ref = blueref, de.method="wilcox",
  labels = blueref$label.fine)

PDAC_primary$SingleR_blue_cell <- blue_cell$labels

PDAC_primary$SingleR_blue_pruned <- blue_cell$pruned.labels

cluster_col <- "seurat_clusters"

meta_tmp <- PDAC_primary@meta.data

meta_tmp$cluster_id <- as.character(meta_tmp[[cluster_col]])

meta_tmp$SingleR_blue_pruned_clean <- ifelse(
  is.na(meta_tmp$SingleR_blue_pruned),
  "Unassigned",
  meta_tmp$SingleR_blue_pruned
)

cluster_annotation_SingleR_blue <- meta_tmp %>%
  dplyr::count(
    cluster_id,
    SingleR_blue_cell,
    SingleR_blue_pruned_clean,
    name = "n"
  ) %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::arrange(as.numeric(cluster_id), desc(frac))

top_cluster_annotation_SingleR_blue <- cluster_annotation_SingleR_blue %>%
  dplyr::slice_max(order_by = frac, n = 1, with_ties = FALSE)

write.csv(cluster_annotation_SingleR_blue, file = "/rds/general/user/ts3225/ephemeral/R/QC4/primary5cluster_annotation_SingleR_blue.csv")
write.csv(top_cluster_annotation_SingleR_blue, file = "/rds/general/user/ts3225/ephemeral/R/QC4/primary5top_cluster_annotation_SingleR_blue.csv")


blue_cell2 <- SingleR(
  test = GetAssayData(PDAC_met, assay = "decontXcounts", layer = "data"),
  ref = blueref, de.method="wilcox",
  labels = blueref$label.fine)

PDAC_met$SingleR_blue_cell <- blue_cell2$labels

PDAC_met$SingleR_blue_pruned <- blue_cell2$pruned.labels


cluster_col <- "seurat_clusters"

meta_tmp <- PDAC_met@meta.data

meta_tmp$cluster_id <- as.character(meta_tmp[[cluster_col]])

meta_tmp$SingleR_blue_pruned_clean <- ifelse(
  is.na(meta_tmp$SingleR_blue_pruned),
  "Unassigned",
  meta_tmp$SingleR_blue_pruned
)

cluster_annotation_SingleR_blue <- meta_tmp %>%
  dplyr::count(
    cluster_id,
    SingleR_blue_cell,
    SingleR_blue_pruned_clean,
    name = "n"
  ) %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::arrange(as.numeric(cluster_id), desc(frac))

top_cluster_annotation_SingleR_blue <- cluster_annotation_SingleR_blue %>%
  dplyr::slice_max(order_by = frac, n = 1, with_ties = FALSE)

write.csv(cluster_annotation_SingleR_blue, file = "/rds/general/user/ts3225/ephemeral/R/QC4/met5cluster_annotation_SingleR_blue.csv")
write.csv(top_cluster_annotation_SingleR_blue, file = "/rds/general/user/ts3225/ephemeral/R/QC4/met5top_cluster_annotation_SingleR_blue.csv")

rm(blueref, blue_cell, blue_cell2, meta_tmp)
gc()

marker_genes = list(
  "T cells" = c("CD3D", "CD3E", "IL7R", "GNLY"),
  "B cells" = c("MS4A1", "CD79A"),
  "Mast cells" = c("KIT"),
  "Macrophages" = c("CD68", "CD14"),
  "Fibroblasts" = c("COL1A1", "DCN"),
  "Hepatocytes" = c("ALB"),
  "Epithelial cells" = c("EPCAM", "KRT19", "KRT18"),
  "Classical subtype" = c("GATA6" , "CLDN18"),
  "Basal subtype" = c("KRT17", "S100A2"),
  "Endothelial cells" = c("PECAM1", "VWF"),
  "Acinar cells" = "REG1A",
  "Endocrine" = "INS",
  "Dendritic cells" = c("FCER1A", "TSPAN13", "GPR183"),
  "Proliferation" = c("MKI67"))

Validation_primary <- DotPlot(
  PDAC_primary,
  features = marker_genes,
  group.by = "seurat_clusters", dot.scale = 8) +
  theme(
    strip.text.x = element_text(size = 8),
    axis.text.x = element_text(
      size = 8,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)) +
  labs(
    x = "Marker genes",
    y = "Seurat Clusters",
    colour = "Average expression",
    size = "Percent expressed")

saveRDS(Validation_primary, file = "/rds/general/user/ts3225/ephemeral/R/QC4/Validation_primary.rds")
ggsave(file = "/rds/general/user/ts3225/ephemeral/R/QC4/Validation_primary.pdf", plot = Validation_primary, width = 20, height = 20)

Validation_met <- DotPlot(
  PDAC_met,
  features = marker_genes,
  group.by = "seurat_clusters", dot.scale = 8) +
  theme(
    strip.text.x = element_text(size = 8),
    axis.text.x = element_text(
      size = 8,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)) +
  labs(
    x = "Marker genes",
    y = "Seurat Clusters",
    colour = "Average expression",
    size = "Percent expressed")

saveRDS(Validation_met, file = "/rds/general/user/ts3225/ephemeral/R/QC4/Validation_met.rds")
ggsave(file = "/rds/general/user/ts3225/ephemeral/R/QC4/Validation_met.pdf", plot = Validation_met, width = 20, height = 20)


meta_met <- as.data.frame(PDAC_primary@meta.data)

sample_col <- "orig.ident"
cluster_col <- "seurat_clusters"
dataset_col <- "dataset"

meta_met[[cluster_col]] <- factor(as.character(meta_met[[cluster_col]]))

meta_met$sample_key <- paste(
  meta_met[[dataset_col]],
  meta_met[[sample_col]],
  sep = "___"
)

counts_mat_met <- table(
  sample_key = meta_met$sample_key,
  seurat_cluster = meta_met[[cluster_col]]
)

counts_df_met <- as.data.frame.matrix(counts_mat_met)

counts_df_met <- data.frame(
  sample_key = rownames(counts_df_met),
  counts_df_met,
  row.names = NULL,
  check.names = FALSE
)

sample_info_primary <- unique(
  meta_met[, c("sample_key", dataset_col, sample_col)]
)

sample_cluster_table_primary <- merge(
  sample_info_primary,
  counts_df_met,
  by = "sample_key",
  all.x = TRUE,
  sort = FALSE
)

sample_cluster_table_primary$sample_key <- NULL

sample_cluster_table_primary <- sample_cluster_table_primary[
  order(
    sample_cluster_table_primary[[dataset_col]],
    sample_cluster_table_primary[[sample_col]]
  ),
]

sample_cluster_table_primary

write.csv(
  sample_cluster_table_primary,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/PDAC_primary_sample_cluster_table.csv",
  row.names = FALSE
)

Seurat_plot_primary <- DimPlot(PDAC_primary, reduction = "umap_harmony", raster = FALSE, group.by = "seurat_clusters", label = T, repel = T)
saveRDS(Seurat_plot_primary, file = "/rds/general/user/ts3225/ephemeral/R/QC4/Seurat_plot_primary.rds")
ggsave(file = "/rds/general/user/ts3225/ephemeral/R/QC4/Seurat_plot_primary.pdf", plot = Seurat_plot_primary, width = 10, height = 10)


meta_met <- as.data.frame(PDAC_met@meta.data)

sample_col <- "orig.ident"
cluster_col <- "seurat_clusters"
dataset_col <- "dataset"

meta_met[[cluster_col]] <- factor(as.character(meta_met[[cluster_col]]))

meta_met$sample_key <- paste(
  meta_met[[dataset_col]],
  meta_met[[sample_col]],
  sep = "___"
)

counts_mat_met <- table(
  sample_key = meta_met$sample_key,
  seurat_cluster = meta_met[[cluster_col]]
)

counts_df_met <- as.data.frame.matrix(counts_mat_met)

counts_df_met <- data.frame(
  sample_key = rownames(counts_df_met),
  counts_df_met,
  row.names = NULL,
  check.names = FALSE
)

sample_info_met <- unique(
  meta_met[, c("sample_key", dataset_col, sample_col)]
)

sample_cluster_table_met <- merge(
  sample_info_met,
  counts_df_met,
  by = "sample_key",
  all.x = TRUE,
  sort = FALSE
)

sample_cluster_table_met$sample_key <- NULL

sample_cluster_table_met <- sample_cluster_table_met[
  order(
    sample_cluster_table_met[[dataset_col]],
    sample_cluster_table_met[[sample_col]]
  ),
]

sample_cluster_table_met

write.csv(
  sample_cluster_table_met,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/PDAC_met_sample_cluster_table.csv",
  row.names = FALSE
)

Seurat_plot_met <- DimPlot(PDAC_met, reduction = "umap_harmony", raster = FALSE, group.by = "seurat_clusters", label = T, repel = T)
saveRDS(Seurat_plot_met, file = "/rds/general/user/ts3225/ephemeral/R/QC4/Seurat_plot_met.rds")
ggsave(file = "/rds/general/user/ts3225/ephemeral/R/QC4/Seurat_plot_met.pdf", plot = Seurat_plot_met, width = 10, height = 10)


saveRDS(
  PDAC_primary,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/PDAC_primary_annotated.rds",
  compress = FALSE)


saveRDS(
  PDAC_met,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/PDAC_met_annotated.rds",
  compress = FALSE)
