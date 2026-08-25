#### azimuth
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Matrix)
library(ggplot2)
library(grid)
library(patchwork)
library(Seurat)
library(liana)
library(tidyverse)
library(magrittr)
library(S4Vectors)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)


PDAC_met <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/PDAC_met_annotated.rds")

library(AzimuthAPI)

DefaultAssay(PDAC_met) <- "decontXcounts"

PDAC_met <- CloudAzimuth(
  object = PDAC_met,
  assay = "decontXcounts"
)
cluster_col <- "seurat_clusters"

meta_tmp <- PDAC_met@meta.data
meta_tmp$cluster_id <- as.character(meta_tmp[[cluster_col]])

meta_tmp$azimuth_label_clean <- as.character(meta_tmp$azimuth_label)
meta_tmp$azimuth_label_clean[
  is.na(meta_tmp$azimuth_label_clean) | meta_tmp$azimuth_label_clean == ""
] <- "Unassigned"

cluster_annotation_azimuth <- meta_tmp %>%
  dplyr::count(cluster_id, azimuth_label_clean, name = "n") %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::arrange(as.numeric(cluster_id), desc(frac))

top_cluster_annotation_azimuth <- cluster_annotation_azimuth %>%
  dplyr::slice_max(order_by = frac, n = 1, with_ties = FALSE)

write.csv(
  cluster_annotation_azimuth,
  "/Users/thirisantracy/Desktop/thesis/met_azimuth_cluster_annotation.csv",
  row.names = FALSE
)

write.csv(
  top_cluster_annotation_azimuth,
  "/Users/thirisantracy/Desktop/thesis/met_azimuth_top_cluster_annotation.csv",
  row.names = FALSE
)


PDAC_primary <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/PDAC_primary_annotated.rds")

library(AzimuthAPI)

DefaultAssay(PDAC_primary) <- "decontXcounts"

PDAC_primary <- CloudAzimuth(
  object = PDAC_primary,
  assay = "decontXcounts"
)


cluster_col <- "seurat_clusters"

meta_tmp <- PDAC_primary@meta.data
meta_tmp$cluster_id <- as.character(meta_tmp[[cluster_col]])

meta_tmp$azimuth_label_clean <- as.character(meta_tmp$azimuth_label)
meta_tmp$azimuth_label_clean[
  is.na(meta_tmp$azimuth_label_clean) | meta_tmp$azimuth_label_clean == ""
] <- "Unassigned"

cluster_annotation_azimuth <- meta_tmp %>%
  dplyr::count(cluster_id, azimuth_label_clean, name = "n") %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::arrange(as.numeric(cluster_id), desc(frac))

top_cluster_annotation_azimuth <- cluster_annotation_azimuth %>%
  dplyr::slice_max(order_by = frac, n = 1, with_ties = FALSE)

write.csv(
  cluster_annotation_azimuth,
  "/Users/thirisantracy/Desktop/thesis/primary_azimuth_cluster_annotation.csv",
  row.names = FALSE
)

write.csv(
  top_cluster_annotation_azimuth,
  "/Users/thirisantracy/Desktop/thesis/primary_azimuth_top_cluster_annotation.csv",
  row.names = FALSE
)

metseuratmarkers0.5 <- read.csv(file = "/Users/thirisantracy/Desktop/thesis/metseuratmarkers0.5.csv")

metseuratsig <- metseuratmarkers0.5 |> 
  filter(p_val_adj < 0.05,
         avg_log2FC > 0.58)

metseuratsig2 <- metseuratmarkers0.5 |> 
  filter(p_val_adj < 0.01,
         avg_log2FC > 1.5)

metseuratsig2 <- metseuratsig2 |> 
  arrange(cluster, desc(avg_log2FC))

write.csv(metseuratsig, file = "/Users/thirisantracy/Desktop/thesis/metseuratsig.csv")
write.csv(metseuratsig2, file = "/Users/thirisantracy/Desktop/thesis/metseuratsig2.csv")






primaryseuratmarkers0.5 <- read.csv(file = "/Users/thirisantracy/Desktop/thesis/primaryseuratmarkers0.5.csv")
primaryseuratsig2 <- primaryseuratmarkers0.5 |> 
  filter(p_val_adj < 0.01,
         avg_log2FC > 1.5)

primaryseuratsig2 <- primaryseuratsig2 |> arrange(cluster, desc(avg_log2FC))
write.csv(primaryseuratsig2, file = "/Users/thirisantracy/Desktop/thesis/primaryseuratsig2.csv")



metseuratsig2 <- read.csv(file = "/Users/thirisantracy/Desktop/thesis/metseuratsig2.csv")
primaryseuratsig2 <- read.csv(file = "/Users/thirisantracy/Desktop/thesis/primaryseuratsig2.csv")


primaryseuratsig2_clean <- primaryseuratsig2 %>%
  select(
    cluster,
    gene,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  ) %>%
  arrange(
    cluster,
    desc(avg_log2FC)
  )

metseuratsig2_clean <- metseuratsig2 %>%
  select(
    cluster,
    gene,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  ) %>%
  arrange(
    cluster,
    desc(avg_log2FC)
  )


primary_top20 <- primaryseuratsig2_clean %>%
  group_by(cluster) %>%
  slice_head(n = 20) %>%
  ungroup()

met_top20 <- metseuratsig2_clean %>%
  group_by(cluster) %>%
  slice_head(n = 20) %>%
  ungroup()


primary_top20_summary <- primary_top20 %>%
  group_by(cluster) %>%
  summarise(
    top_markers = paste(gene, collapse = ", "),
    .groups = "drop"
  )

met_top20_summary <- met_top20 %>%
  group_by(cluster) %>%
  summarise(
    top_markers = paste(gene, collapse = ", "),
    .groups = "drop"
  )


write.csv(
  primary_top20_summary,
  file = "/Users/thirisantracy/Desktop/thesis/primary_top20_marker_summary.csv",
  row.names = FALSE
)

write.csv(
  met_top20_summary,
  file = "/Users/thirisantracy/Desktop/thesis/met_top20_marker_summary.csv",
  row.names = FALSE
)

write.csv(
  primary_top20,
  file = "/Users/thirisantracy/Desktop/thesis/primary_top20.csv",
  row.names = FALSE
)

write.csv(
  met_top20,
  file = "/Users/thirisantracy/Desktop/thesis/met_top20.csv",
  row.names = FALSE
)




primary_ <- primaryseuratsig2_clean %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()

met_ <- metseuratsig2_clean %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()


primary__summary <- primary_ %>%
  group_by(cluster) %>%
  summarise(
    top_markers = paste(gene, collapse = ", "),
    .groups = "drop"
  )

met__summary <- met_ %>%
  group_by(cluster) %>%
  summarise(
    top_markers = paste(gene, collapse = ", "),
    .groups = "drop"
  )


write.csv(
  primary__summary,
  file = "/Users/thirisantracy/Desktop/thesis/primary__marker_summary.csv",
  row.names = FALSE
)

write.csv(
  met__summary,
  file = "/Users/thirisantracy/Desktop/thesis/met__marker_summary.csv",
  row.names = FALSE
)

write.csv(
  primary_,
  file = "/Users/thirisantracy/Desktop/thesis/primary_.csv",
  row.names = FALSE
)

write.csv(
  met_,
  file = "/Users/thirisantracy/Desktop/thesis/met_.csv",
  row.names = FALSE
)



pdac_validation_markers <- list(
  "Ductal cells" = c("EPCAM", "KRT19", "KRT7"),
  "Normal-like ductal cells" = c("AMBP", "CFTR", "SLC4A4"),
  "Tumour epithelial cells" = c("CEACAM5", "CEACAM6", "S100P"),
  "T cells" = c("CD3D", "CD3E", "TRAC"),
  "NK cells" = c("GNLY", "NKG7", "KLRD1"),
  "Dendritic cells" = c("CD1C", "FCER1A", "CLEC10A"),
  "TAMs" = c("CD68", "C1QA", "APOE"),
  "Monocytes" = c("FCN1", "VCAN", "S100A8"),
  "CAFs" = c("COL1A1", "DCN", "LUM"),
  "Endothelial cells" = c("VWF", "PECAM1", "KDR"),
  "B/plasma cells" = c("MS4A1", "CD79A","MZB1", "JCHAIN"),
  "Mast cells" = c("CPA3", "TPSAB1", "KIT"),
  "Endocrine cells" = c("CHGA", "CHGB", "SYP"),
  "Hepatocytes" = c("ALB", "APOA1", "TTR"),
  "Acinar cells" = c("PRSS1", "CTRB1", "CTRB2"))

names(pdac_validation_markers) <- stringr::str_wrap(names(pdac_validation_markers), width = 18)

DotPlot(PDAC_primary,
        features = pdac_validation_markers,
        scale = TRUE,
        dot.scale = 10,
        group.by = "seurat_clusters") +
  labs(x = "Marker Genes", y = "Seurat Cluster") +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        strip.text.x = element_text(size = 12, lineheight = 0.9))

DotPlot(PDAC_met,
        features = pdac_validation_markers,
        scale = TRUE,
        dot.scale = 10,
        group.by = "seurat_clusters") +
  labs(x = "Marker Genes", y = "Seurat Cluster") +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        strip.text.x = element_text(size = 12, lineheight = 0.9))


PDAC_primary <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/PDAC_primary_annotated.rds")
PDAC_met <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/PDAC_met_annotated.rds")



primary_annotations <- c(
  "0" = "Tumour epithelial cells",
  "1" = "CD4 T cells",
  "2" = "CD8 T cells",
  "3" = "TAMs",
  "4" = "NK cells",
  "5" = "CAFs",
  "6" = "Monocytes", #MDSCs
  "7" = "Low-quality cells",
  "8" = "Tumour epithelial cells",
  "9" = "Cycling tumour epithelial cells", #cycling
  "10" = "B cells",
  "11" = "Monocytes",
  "12" = "Normal-like ductal cells",
  "13" = "Dendritic cells",
  "14" = "Cycling T/NK cells",
  "15" = "Endothelial cells",
  "16" = "CD4 T cells", #P21
  "17" = "Pericytes",
  "18" = "Mast cells",
  "19" = "Plasma cells",
  "20" = "Cycling TAMs", #cycling
  "21" = "Acinar cells",
  "22" = "Tumour epithelial cells", #P17_basal-like_EMT
  "23" = "Tumour epithelial cells",
  "24" = "Mesothelial cells",
  "25" = "Tuft-like epithelial cells") #P12_onefrom_P18_earlyADM


primary_DE_annotations <- c(
  "0" = "Tumour epithelial cells",
  "1" = "T/NK cells",
  "2" = "T/NK cells",
  "3" = "TAMs",
  "4" = "T/NK cells",
  "5" = "CAFs",
  "6" = "Monocytes", #MDSCs
  "7" = "Low-quality cells", #getrid
  "8" = "Tumour epithelial cells",
  "9" = "Tumour epithelial cells", #cycling
  "10" = "B/plasma cells",
  "11" = "Monocytes",
  "12" = "Normal-like ductal cells",
  "13" = "Dendritic cells",
  "14" = "T/NK cells",
  "15" = "Endothelial cells",
  "16" = "T/NK cells", #P21
  "17" = "Pericytes",
  "18" = "Mast cells",
  "19" = "B/plasma cells",
  "20" = "TAMs", #cycling
  "21" = "Acinar cells",
  "22" = "Tumour epithelial cells", #P17_basal-like_EMT
  "23" = "Tumour epithelial cells",
  "24" = "Mesothelial cells",
  "25" = "Tuft-like epithelial cells") #P12_onefrom_P18_earlyADM


met_annotations <- c(
  "0" = "Tumour epithelial cells",
  "1" = "T/NK cells",
  "2" = "CD4 T cells",
  "3" = "TAMs",
  "4" = "Cycling tumour epithelial cells", #cycling
  "5" = "Tumour epithelial cells",
  "6" = "Tumour epithelial cells",
  "7" = "Monocytes", #MDSC
  "8" = "CD4 T cells", #Treg
  "9" = "Tumour epithelial cells",
  "10" = "TAMs",
  "11" = "Tumour epithelial cells", #M14
  "12" = "Neutrophils",
  "13" = "CAFs", #pericytestoo
  "14" = "B/plasma cells",
  "15" = "Cycling immune cells", #cycling
  "16" = "Tumour epithelial cells", #M14
  "17" = "Endothelial cells",
  "18" = "Hepatocytes",
  "19" = "Tumour epithelial cells",
  "20" = "Tumour epithelial cells", #M13
  "21" = "Tumour epithelial cells", #M08
  "22" = "Tumour epithelial cells", #M24
  "23" = "Tumour epithelial cells", #M19
  "24" = "Mast cells", #M24
  "25" = "Cycling tumour epithelial cells", #M14_cycling
  "26" = "Mixed epithelial/immune cells", #M08&M09
  "27" = "Mixed epithelial-myeloid cells", #M13
  "28" = "Tumour epithelial cells") #M11

met_DE_annotations <- c(
  "0" = "Tumour epithelial cells",
  "1" = "T/NK cells",
  "2" = "T/NK cells",
  "3" = "TAMs",
  "4" = "Tumour epithelial cells", #cycling
  "5" = "Tumour epithelial cells",
  "6" = "Tumour epithelial cells",
  "7" = "Monocytes", #MDSC
  "8" = "T/NK cells", #Treg
  "9" = "Tumour epithelial cells",
  "10" = "TAMs",
  "11" = "Tumour epithelial cells", #M14
  "12" = "Neutrophils",
  "13" = "CAFs", #pericytestoo
  "14" = "B/plasma cells",
  "15" = "Cycling immune cells", 
  "16" = "Tumour epithelial cells", #M14
  "17" = "Endothelial cells",
  "18" = "Hepatocytes",
  "19" = "Tumour epithelial cells",
  "20" = "Tumour epithelial cells", #M13
  "21" = "Tumour epithelial cells", #M08
  "22" = "Tumour epithelial cells", #M24
  "23" = "Tumour epithelial cells", #M19
  "24" = "Mast cells", #M24
  "25" = "Tumour epithelial cells", #M14_cycling
  "26" = "Mixed epithelial/immune cells", #M08&M09
  "27" = "Mixed epithelial-myeloid cells", #M13
  "28" = "Tumour epithelial cells") #M11



PDAC_primary$celltype_annotation <- unname(primary_annotations[as.character(PDAC_primary$seurat_clusters)])
PDAC_met$celltype_annotation <- unname(met_annotations[as.character(PDAC_met$seurat_clusters)])

table(PDAC_primary$celltype_annotation, useNA = "ifany")
table(PDAC_met$celltype_annotation, useNA = "ifany")

PDAC_primary$DEcelltype_annotation <- unname(primary_DE_annotations[as.character(PDAC_primary$seurat_clusters)])
PDAC_met$DEcelltype_annotation <- unname(met_DE_annotations[as.character(PDAC_met$seurat_clusters)])

table(PDAC_primary$DEcelltype_annotation, useNA = "ifany")
table(PDAC_met$DEcelltype_annotation, useNA = "ifany")
dataset_names <- c("one" = "GSE154778", "two" = "GSE205013", "three" = "GSE263733", "four" = "GSE197177")
PDAC_primary$dataset_GSE <- dplyr::recode(as.character(PDAC_primary$dataset),
                                          one = "GSE154778", two = "GSE205013", three = "GSE263733", four = "GSE197177")

PDAC_met$dataset_GSE <- dplyr::recode(as.character(PDAC_met$dataset),
                                      one = "GSE154778", two = "GSE205013", three = "GSE263733", four = "GSE197177")

table(PDAC_primary$dataset_GSE)
table(PDAC_met$dataset_GSE)
primary_types <- sort(unique(na.omit(
  as.character(PDAC_primary_filtered$celltype_annotation)
)))

met_types <- sort(unique(na.omit(
  as.character(PDAC_met_filtered$celltype_annotation)
)))
primary_types
met_types
length(primary_types)
length(met_types)

shared_types <- sort(intersect(primary_types, met_types))
primary_only_types <- sort(setdiff(primary_types, met_types))
met_only_types <- sort(setdiff(met_types, primary_types))

all_celltypes <- c(
  shared_types,
  primary_only_types,
  met_only_types
)
celltype_colours <- setNames(
  Seurat::DiscretePalette(
    length(all_celltypes),
    palette = "glasbey"
  ),
  all_celltypes
)


celltype_colours <- colorspace::lighten(
  celltype_colours,
  amount = 0.22,
  space = "HCL"
)

manual_colours <- c(
  "Tumour epithelial cells"          = "#4C78A8",
  "Cycling tumour epithelial cells" = "#5BC0BE",
  
  "CD4 T cells"                      = "#E76FAD",
  "CD8 T cells"                      = "#59A14F",
  "T/NK cells"                       = "#9C89B8",
  "Cycling T/NK cells"               = "#7AC7E8",
  "Cycling immune cells"             = "#D6A5E8",
  "NK cells"                         = "#2FA84F",
  
  "TAMs"                             = "#D95F4C",
  "Cycling TAMs"                     = "#F4A582",
  "Monocytes"                        = "#A05A9C",
  "Neutrophils"                      = "#E5B82E",
  
  "CAFs"                             = "#667BC6",
  "Endothelial cells"                = "#2A9D8F",
  "Dendritic cells"                  = "#9B79B6",
  
  "Low-quality cells"                = "#C8C8C8",
  "Mixed epithelial/immune cells"    = "#999999",
  "Mixed epithelial-myeloid cells"   = "#6F6F6F"
)

types_to_override <- intersect(
  names(manual_colours),
  names(celltype_colours)
)

celltype_colours[types_to_override] <-
  manual_colours[types_to_override]

celltype_colours
PDAC_primary$celltype_annotation <- factor(
  PDAC_primary$celltype_annotation,
  levels = all_celltypes
)

PDAC_met$celltype_annotation <- factor(
  PDAC_met$celltype_annotation,
  levels = all_celltypes
)

primary_plot <- DimPlot(PDAC_primary, reduction = "umap_harmony", group.by = "celltype_annotation", pt.size = 0.1, raster = FALSE)
primary_plot <- LabelClusters(primary_plot,  id = "celltype_annotation",repel = TRUE, size = 3.5, fontface = "bold") +
  labs(title = NULL, colour = NULL) +
  scale_colour_manual(values = celltype_colours) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 9))
primary_plot


met_plot <- DimPlot(PDAC_met, reduction = "umap_harmony", group.by = "celltype_annotation", pt.size = 0.1, raster = FALSE)
met_plot <- LabelClusters(met_plot, id = "celltype_annotation", repel = TRUE, size = 3.5, fontface = "bold") +
  labs(title = NULL, colour = NULL) +
  scale_colour_manual(values = celltype_colours) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), legend.text = element_text(size = 9))
met_plot



primary_cells_keep <- colnames(PDAC_primary)[
  as.character(PDAC_primary$seurat_clusters) != "7"
]

met_cells_keep <- colnames(PDAC_met)[
  !as.character(PDAC_met$seurat_clusters) %in% c("15", "26", "27")
]

PDAC_primary_filtered <- subset(
  PDAC_primary,
  cells = primary_cells_keep
)

PDAC_met_filtered <- subset(
  PDAC_met,
  cells = met_cells_keep
)

table(PDAC_primary_filtered$celltype_annotation, useNA = "ifany")
table(PDAC_met_filtered$celltype_annotation, useNA = "ifany")


primary_filtered_plot <- DimPlot(PDAC_primary_filtered, reduction = "umap_harmony", group.by = "celltype_annotation", pt.size = 0.1, raster = FALSE)
primary_filtered_plot <- LabelClusters(primary_filtered_plot, id = "celltype_annotation", repel = TRUE, size = 3, fontface = "bold") +
  labs(title = NULL, colour = NULL) +
  scale_colour_manual(values = celltype_colours) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), legend.text = element_text(size = 9))

primary_filtered_plot

met_filtered_plot <- DimPlot(PDAC_met_filtered, reduction = "umap_harmony", group.by = "celltype_annotation", pt.size = 0.1, raster = FALSE)
met_filtered_plot <- LabelClusters(met_filtered_plot, id = "celltype_annotation", repel = TRUE, size = 3, fontface = "bold") +
  labs(title = NULL, colour = NULL) +
  scale_colour_manual(values = celltype_colours) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), legend.text = element_text(size = 9))
met_filtered_plot

ggsave(filename = "/Users/thirisantracy/Desktop/thesis/images/primary_filtered_plot.png",
       primary_filtered_plot, width = 10, height = 10, units = "in", dpi = 400, bg = "white")

ggsave(filename = "/Users/thirisantracy/Desktop/thesis/images/met_filtered_plot.png",
       met_filtered_plot, width = 10, height = 10, units = "in", dpi = 400, bg = "white")

# PDF copies used in the thesis figures supplied with this script
ggsave(
  filename = "/Users/thirisantracy/Desktop/thesis/images/filteredprimary.pdf",
  plot = primary_filtered_plot,
  width = 10,
  height = 10,
  units = "in",
  bg = "white"
)

ggsave(
  filename = "/Users/thirisantracy/Desktop/thesis/images/filteredmet.pdf",
  plot = met_filtered_plot,
  width = 10,
  height = 10,
  units = "in",
  bg = "white"
)

saveRDS(PDAC_primary, file = "/Users/thirisantracy/Desktop/thesis/PDAC_primary_annotated.rds")
saveRDS(PDAC_met, file = "/Users/thirisantracy/Desktop/thesis/PDAC_met_annotated.rds")



saveRDS(PDAC_primary_filtered, file = "/Users/thirisantracy/Desktop/thesis/PDAC_primary_filtered.rds")
saveRDS(PDAC_met_filtered, file = "/Users/thirisantracy/Desktop/thesis/PDAC_met_filtered.rds")
PDAC_primary_filtered <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/PDAC_primary_filtered.rds")
PDAC_met_filtered <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/PDAC_met_filtered.rds")

map_symbols <- function(obj) {
  obj <- JoinLayers(obj, assay = "decontXcounts")
  counts <- LayerData(obj, assay = "decontXcounts", layer = "counts")
  old <- rownames(counts); ensg <- sub("\\..*$", "", old); is_ensg <- grepl("^ENSG", ensg)
  
  map <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(ensg[is_ensg]),
                               "SYMBOL", "ENSEMBL", multiVals = "first")
  genes <- old
  genes[is_ensg] <- unname(map[ensg[is_ensg]])
  
  keep <- !is.na(genes) & genes != ""
  counts <- counts[keep, ]; genes <- genes[keep]; unique_genes <- unique(genes)
  
  collapse <- Matrix::sparseMatrix(i = match(genes, unique_genes), j = seq_along(genes), x = 1,
                                   dims = c(length(unique_genes), length(genes)),
                                   dimnames = list(unique_genes, NULL))
  
  obj[["decontXsymbols"]] <- CreateAssay5Object(counts = collapse %*% counts)
  DefaultAssay(obj) <- "decontXsymbols"
  obj
}

PDAC_primary_filtered <- map_symbols(PDAC_primary_filtered)
PDAC_met_filtered <- map_symbols(PDAC_met_filtered)

sum(grepl("^ENSG", rownames(PDAC_primary_filtered[["decontXsymbols"]])))
sum(duplicated(rownames(PDAC_primary_filtered[["decontXsymbols"]])))

sum(grepl("^ENSG", rownames(PDAC_met_filtered[["decontXsymbols"]])))
sum(duplicated(rownames(PDAC_met_filtered[["decontXsymbols"]])))

old_counts <- LayerData(PDAC_primary_filtered, assay = "decontXcounts", layer = "counts")
new_counts <- LayerData(PDAC_primary_filtered, assay = "decontXsymbols", layer = "counts")

c(same_cells = identical(colnames(old_counts), colnames(new_counts)),
  old_features = nrow(old_counts),
  new_features = nrow(new_counts),
  ENSG_remaining = sum(grepl("^ENSG", rownames(new_counts))),
  duplicates = sum(duplicated(rownames(new_counts))),
  counts_removed = sum(old_counts) - sum(new_counts))

saveRDS(PDAC_primary_filtered, file = "/Users/thirisantracy/Desktop/thesis/PDAC_primary_filtered.rds")
saveRDS(PDAC_met_filtered, file = "/Users/thirisantracy/Desktop/thesis/PDAC_met_filtered.rds")

PDAC_primary_filtered <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/PDAC_primary_filtered.rds")
PDAC_met_filtered <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/PDAC_met_filtered.rds")
subject_key <- c(
  setNames(sprintf("P%02d", 1:10), sprintf("P%02d", 1:10)),
  setNames(sprintf("P%02d", 16:26), sprintf("P%02d", 11:21)),
  setNames(sprintf("P%02d", 36:56), sprintf("P%02d", 22:42)),
  setNames(sprintf("P%02d", 57:59), sprintf("P%02d", 43:45)),
  setNames(sprintf("P%02d", c(11:15, 27:31)), sprintf("M%02d", 1:10)),
  setNames(sprintf("P%02d", 32:35), sprintf("M%02d", 11:14)),
  setNames(sprintf("P%02d", c(38, 39, 44, 51, 53, 55, 56)), sprintf("M%02d", 15:21)),
  setNames(sprintf("P%02d", 57:60), sprintf("M%02d", 22:25))
)

PDAC_primary$subject_ID <- unname(subject_key[as.character(PDAC_primary$patient)])
PDAC_met$subject_ID <- unname(subject_key[as.character(PDAC_met$patient)])

PDAC_primary_filtered$subject_ID <- unname(subject_key[as.character(PDAC_primary_filtered$patient)])
PDAC_met_filtered$subject_ID <- unname(subject_key[as.character(PDAC_met_filtered$patient)])

primary_mapping <- unique(PDAC_primary_filtered@meta.data[, c("patient", "subject_ID")])
met_mapping <- unique(PDAC_met_filtered@meta.data[, c("patient", "subject_ID")])

primary_mapping <- primary_mapping[order(primary_mapping$patient), ]
met_mapping <- met_mapping[order(met_mapping$patient), ]

View(primary_mapping)
View(met_mapping)

library(dplyr)
library(tibble)

add_n_cells <- function(obj) {
  
  meta <- as.data.frame(obj@meta.data)
  
  required_columns <- c(
    "dataset_GSE",
    "tumor",
    "patient",
    "DEcelltype_annotation"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(meta)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  # Convert list/factor metadata columns into ordinary character vectors
  convert_to_character <- function(x, column_name) {
    
    if (is.data.frame(x) || is.matrix(x)) {
      x <- x[, 1]
    }
    
    if (is.list(x)) {
      x <- vapply(
        x,
        function(value) {
          if (length(value) == 0 || all(is.na(value))) {
            return(NA_character_)
          }
          
          as.character(value[[1]])
        },
        character(1)
      )
    } else {
      x <- as.character(x)
    }
    
    if (length(x) != nrow(meta)) {
      stop(
        column_name,
        " does not have one value per cell."
      )
    }
    
    x
  }
  
  columns_to_convert <- intersect(
    c(
      "dataset_GSE",
      "tumor",
      "patient",
      "subject_ID",
      "DEcelltype_annotation"
    ),
    colnames(meta)
  )
  
  for (column_name in columns_to_convert) {
    meta[[column_name]] <- convert_to_character(
      meta[[column_name]],
      column_name
    )
  }
  
  # Check that global grouping information is present for every cell
  if (
    anyNA(meta$dataset_GSE) ||
    anyNA(meta$tumor) ||
    anyNA(meta$patient)
  ) {
    stop(
      "Some cells have missing dataset_GSE, tumor or patient values."
    )
  }
  
  # Number of cells in each:
  # dataset × tumour condition × patient
  global_group <- interaction(
    meta$dataset_GSE,
    meta$tumor,
    meta$patient,
    drop = TRUE,
    lex.order = TRUE
  )
  
  meta$n_cells_global <- as.integer(
    ave(
      rep.int(1L, nrow(meta)),
      global_group,
      FUN = sum
    )
  )
  
  # Number of cells in each:
  # dataset × tumour condition × patient × broad cell type
  meta$n_cells_celltype <- NA_integer_
  
  use_for_celltype <- (
    !is.na(meta$DEcelltype_annotation) &
      nzchar(meta$DEcelltype_annotation)
  )
  
  celltype_group <- interaction(
    meta$dataset_GSE[use_for_celltype],
    meta$tumor[use_for_celltype],
    meta$patient[use_for_celltype],
    meta$DEcelltype_annotation[use_for_celltype],
    drop = TRUE,
    lex.order = TRUE
  )
  
  meta$n_cells_celltype[use_for_celltype] <- as.integer(
    ave(
      rep.int(1L, sum(use_for_celltype)),
      celltype_group,
      FUN = sum
    )
  )
  
  obj@meta.data <- meta
  
  return(obj)
}

PDAC_primary_filtered <- add_n_cells(
  PDAC_primary_filtered
)

PDAC_met_filtered <- add_n_cells(
  PDAC_met_filtered
)


head(
  PDAC_primary_filtered@meta.data[
    ,
    c(
      "dataset_GSE",
      "tumor",
      "patient",
      "subject_ID",
      "DEcelltype_annotation",
      "n_cells_global",
      "n_cells_celltype"
    )
  ]
)



library(dplyr)

primary_global_ncells <- PDAC_primary_filtered@meta.data %>%
  distinct(
    dataset_GSE,
    tumor,
    patient,
    subject_ID,
    n_cells_global
  ) %>%
  arrange(dataset_GSE, patient)

met_global_ncells <- PDAC_met_filtered@meta.data %>%
  distinct(
    dataset_GSE,
    tumor,
    patient,
    subject_ID,
    n_cells_global
  ) %>%
  arrange(dataset_GSE, patient)

View(primary_global_ncells)
View(met_global_ncells)

primary_celltype_ncells <- PDAC_primary_filtered@meta.data %>%
  filter(!is.na(DEcelltype_annotation)) %>%
  distinct(
    dataset_GSE,
    tumor,
    patient,
    subject_ID,
    DEcelltype_annotation,
    n_cells_celltype
  ) %>%
  arrange(
    DEcelltype_annotation,
    dataset_GSE,
    patient
  )

met_celltype_ncells <- PDAC_met_filtered@meta.data %>%
  filter(!is.na(DEcelltype_annotation)) %>%
  distinct(
    dataset_GSE,
    tumor,
    patient,
    subject_ID,
    DEcelltype_annotation,
    n_cells_celltype
  ) %>%
  arrange(
    DEcelltype_annotation,
    dataset_GSE,
    patient
  )

View(primary_celltype_ncells)
View(met_celltype_ncells)

saveRDS(PDAC_primary_filtered, file = "/Users/thirisantracy/Desktop/thesis/PDAC_primary_filtered.rds")
saveRDS(PDAC_met_filtered, file = "/Users/thirisantracy/Desktop/thesis/PDAC_met_filtered.rds")
