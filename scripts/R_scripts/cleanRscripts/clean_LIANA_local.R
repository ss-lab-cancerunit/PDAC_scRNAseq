library(Seurat)
library(dplyr)
library(Matrix)
library(liana)
library(tidyverse)
library(magrittr)
library(S4Vectors)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(patchwork)
library(tidyr)

PDAC_primary_filtered <- readRDS(file = "/Users/thirisantracy/Desktop/PBS2/PDAC_primary_filtered.rds")
PDAC_met_filtered <- readRDS(file = "/Users/thirisantracy/Desktop/PBS2/PDAC_met_filtered.rds")


DefaultAssay(PDAC_primary_filtered) <- "decontXcounts"
DefaultAssay(PDAC_met_filtered) <- "decontXcounts"

sum(grepl("^ENSG", rownames(PDAC_primary_filtered)))
sum(grepl("^ENSG", rownames(PDAC_met_filtered)))

ncol(PDAC_primary_filtered)
ncol(PDAC_met_filtered)

head(colnames(PDAC_primary_filtered))

colnames(PDAC_primary_filtered@meta.data)

Layers(PDAC_primary_filtered[["decontXcounts"]])
Layers(PDAC_met_filtered[["decontXcounts"]])


PDAC_primary_filtered <- NormalizeData(PDAC_primary_filtered)
PDAC_met_filtered <- NormalizeData(PDAC_met_filtered)




primary_liana_annotations <- c(
  "0" = "Tumour epithelial cells",
  "1" = "CD4 T cells",
  "2" = "CD8 T cells",
  "3" = "TAMs",
  "4" = "NK cells",
  "5" = "CAFs",
  "6" = "Monocytes", 
  "8" = "Tumour epithelial cells",
  "9" = "Tumour epithelial cells",
  "10" = "B cells",
  "11" = "Monocytes",
  "13" = "Dendritic cells",
  "15" = "Endothelial cells",
  "16" = "CD4 T cells", 
  "17" = "Pericytes", 
  "18" = "Mast cells",
  "19" = "Plasma cells",
  "20" = "TAMs", 
  "22" = "Tumour epithelial cells", 
  "23" = "Tumour epithelial cells") 

met_liana_annotations <- c(
  "0" = "Tumour epithelial cells",
  "2" = "CD4 T cells",
  "3" = "TAMs",
  "4" = "Tumour epithelial cells", #cycling
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
  "16" = "Tumour epithelial cells", #M14
  "17" = "Endothelial cells",
  "19" = "Tumour epithelial cells",
  "20" = "Tumour epithelial cells", #M13
  "21" = "Tumour epithelial cells", #M08
  "22" = "Tumour epithelial cells", #M24
  "23" = "Tumour epithelial cells", #M19
  "24" = "Mast cells",
  "25" = "Tumour epithelial cells",
  "28" = "Tumour epithelial cells")

PDAC_primary_filtered$liana_celltype <- unname(
  primary_liana_annotations[
    as.character(PDAC_primary_filtered$seurat_clusters)
  ]
)

PDAC_met_filtered$liana_celltype <- unname(
  met_liana_annotations[
    as.character(PDAC_met_filtered$seurat_clusters)
  ]
)




primary_liana_obj <- subset(
  PDAC_primary_filtered,
  cells = colnames(PDAC_primary_filtered)[
    !is.na(PDAC_primary_filtered$liana_celltype)
  ]
)

met_liana_obj <- subset(
  PDAC_met_filtered,
  cells = colnames(PDAC_met_filtered)[
    !is.na(PDAC_met_filtered$liana_celltype)
  ]
)

table(primary_liana_obj$liana_celltype)
table(met_liana_obj$liana_celltype)


primary_compare_annotations <- c(
  "0" = "Tumour epithelial cells",
  "1" = "CD4 T cells",
  "2" = "T/NK cells",
  "3" = "TAMs",
  "4" = "T/NK cells",
  "5" = "CAFs",
  "6" = "Monocytes",
  "8" = "Tumour epithelial cells",
  "9" = "Tumour epithelial cells",
  "10" = "B/plasma cells",
  "11" = "Monocytes",
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
  "25" = "Tuft-like epithelial cells") #P12_onefrom_P18_ea

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





saveRDS(PDAC_met_filtered, file = "/Users/thirisantracy/Desktop/thesis/PDAC_met_filtered.rds")
saveRDS(PDAC_primary_filtered, file = "/Users/thirisantracy/Desktop/thesis/PDAC_primary_filtered.rds")

rm(PDAC_primary_filtered)

saveRDS(primary_liana_obj, file = "/Users/thirisantracy/Desktop/thesis/primary_liana_obj.rds")
saveRDS(met_liana_obj, file = "/Users/thirisantracy/Desktop/thesis/met_liana_obj.rds")

primary_liana_obj <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/primary_liana_obj.rds")
met_liana_obj <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/met_liana_obj.rds")
###### hpc script but running with PDAC_primary_liana and met, no ENSG in there

######### actual hpc script with liana.R
#### preparation

assay_name <- "decontXcounts"
celltype_col <- "liana_celltype"


DefaultAssay(primary_liana_obj) <- assay_name
DefaultAssay(met_liana_obj) <- assay_name



to_liana_sce <- function(obj, assay_name = "decontXcounts") {
  counts_mat <- LayerData(
    object = obj,
    assay = assay_name,
    layer = "counts"
  )
  
  data_mat <- LayerData(
    object = obj,
    assay = assay_name,
    layer = "data"
  )
  
  sce <- SingleCellExperiment(
    assays = list(
      counts = counts_mat,
      logcounts = data_mat
    ),
    colData = obj@meta.data
  )
  
  return(sce)
}

primary_liana_sce <- to_liana_sce(primary_liana_obj, assay_name)
met_liana_sce <- to_liana_sce(met_liana_obj, assay_name)

assayNames(primary_liana_sce)
assayNames(met_liana_sce)

sum(grepl("^ENSG", rownames(primary_liana_sce)))
sum(grepl("^ENSG", rownames(met_liana_sce)))

table(colData(primary_liana_sce)[[celltype_col]])
table(colData(met_liana_sce)[[celltype_col]])

saveRDS(primary_liana_sce, "/Users/thirisantracy/Desktop/thesis/primary_liana_sce.rds")

saveRDS(met_liana_sce, "/Users/thirisantracy/Desktop/thesis/met_liana_sce.rds")

######common cell types between the two and use the primary_common and met_common
ROOT<-"/Users/thirisantracy/Desktop/thesis/liana_results"

primary_res<-readRDS(file.path(ROOT,"primary_liana_aggregated.rds"))
met_res<-readRDS(file.path(ROOT,"met_liana_aggregated.rds"))

dim(primary_res)
dim(met_res)

colnames(primary_res)
head(primary_res)

colnames(met_res)
head(met_res)

primary_sig<-primary_res%>%
  dplyr::filter(aggregate_rank<=0.01)

met_sig<-met_res%>%
  dplyr::filter(aggregate_rank<=0.01)

nrow(primary_sig)
nrow(met_sig)


library(ComplexHeatmap)

primary_heat<-liana::heat_freq(primary_sig)
ComplexHeatmap::draw(primary_heat)


met_heat<-liana::heat_freq(met_sig)
ComplexHeatmap::draw(met_heat)

primary_heat<-liana::heat_freq(
  primary_sig,
  grid_text=TRUE,
  font_size=10,
  name="Frequency"
)

met_heat<-liana::heat_freq(
  met_sig,
  grid_text=TRUE,
  font_size=10,
  name="Frequency"
)
ComplexHeatmap::draw(primary_heat)
ComplexHeatmap::draw(met_heat)
############################

key <- c(
  "source",
  "target",
  "ligand.complex",
  "receptor.complex"
)

primary_compare <- primary_sig %>%
  dplyr::select(
    source,
    target,
    ligand.complex,
    receptor.complex,
    aggregate_rank,
    sca.LRscore,
    natmi.edge_specificity
  ) %>%
  dplyr::rename(
    rank_primary = aggregate_rank,
    magnitude_primary = sca.LRscore,
    specificity_primary = natmi.edge_specificity
  )

met_compare <- met_sig %>%
  dplyr::select(
    source,
    target,
    ligand.complex,
    receptor.complex,
    aggregate_rank,
    sca.LRscore,
    natmi.edge_specificity
  ) %>%
  dplyr::rename(
    rank_metastasis = aggregate_rank,
    magnitude_metastasis = sca.LRscore,
    specificity_metastasis = natmi.edge_specificity
  )

liana_comparison <- dplyr::full_join(
  primary_compare,
  met_compare,
  by = key
)
#########liana_analysis
############################################################
# LIANA sender-receiver circular network
############################################################

library(dplyr)
library(tidyr)
library(igraph)
library(scales)

############################################################
# Cell types present in LIANA results
############################################################

primary_cells <- sort(
  unique(
    c(
      as.character(primary_sig$source),
      as.character(primary_sig$target)
    )
  )
)

met_cells <- sort(
  unique(
    c(
      as.character(met_sig$source),
      as.character(met_sig$target)
    )
  )
)


primary_cells
met_cells

############################################################
# Cell types shared between the two conditions
############################################################

common_cells <- intersect(
  primary_cells,
  met_cells
)

common_cells

cell_order <- c(
  "Tumour epithelial cells",
  "CAFs",
  "Endothelial cells",
  "TAMs",
  "Monocytes",
  "T/NK cells",
  "B/plasma cells",
  "Mast cells"
)

cell_order
make_network_edges <- function(
    sig_results,
    cells_to_keep
) {
  
  sig_results %>%
    
    dplyr::filter(
      source %in% cells_to_keep,
      target %in% cells_to_keep
    ) %>%
    
    dplyr::distinct(
      source,
      target,
      ligand.complex,
      receptor.complex
    ) %>%
    
    dplyr::count(
      source,
      target,
      name = "n_interactions"
    ) %>%
    
    dplyr::filter(
      n_interactions > 0
    )
}

primary_edges <- make_network_edges(
  primary_sig,
  cell_order
)

met_edges <- make_network_edges(
  met_sig,
  cell_order
)

met_edges <- make_network_edges(
  met_sig,
  cell_order
)

primary_edges %>%
  arrange(desc(n_interactions)) %>%
  print(n = Inf)

met_edges %>%
  arrange(desc(n_interactions)) %>%
  print(n = Inf)

library(igraph)
library(scales)

############################################################
# Remove self-signalling for main figure
############################################################

primary_edges_no_self <- primary_edges %>%
  dplyr::filter(
    source != target
  )

met_edges_no_self <- met_edges %>%
  dplyr::filter(
    source != target
  )

primary_edges %>%
  dplyr::filter(source == target)

met_edges %>%
  dplyr::filter(source == target)



############################################################
# Graph objects
############################################################

primary_graph_full <- make_liana_graph(
  primary_edges,
  cell_order
)

met_graph_full <- make_liana_graph(
  met_edges,
  cell_order
)

primary_graph_no_self <- make_liana_graph(
  primary_edges_no_self,
  cell_order
)

met_graph_no_self <- make_liana_graph(
  met_edges_no_self,
  cell_order
)

############################################################
# Cleaner LIANA circular network plot
############################################################

############################################################
# Cleaner LIANA circular network plot
# with outward-facing labels
############################################################
############################################################
# LIANA circular network
# manual labels positioned outside network
############################################################

plot_liana_network_clean <- function(
    graph,
    title_text
) {
  
  ##########################################################
  # Fixed node positions
  ##########################################################
  
  plot_layout <- fixed_layout[
    igraph::V(graph)$name,
    ,
    drop = FALSE
  ]
  
  
  ##########################################################
  # Edge information
  ##########################################################
  
  edge_counts <- igraph::E(graph)$n_interactions
  
  edge_sources <- igraph::ends(
    graph,
    igraph::E(graph),
    names = TRUE
  )[, 1]
  
  
  ##########################################################
  # SAME edge-width scale for primary + metastatic
  ##########################################################
  
  edge_widths <- 0.3 +
    3.5 *
    sqrt(
      edge_counts / max_interactions
    )
  
  
  ##########################################################
  # Edge colour = sender cell type
  ##########################################################
  
  edge_colours <- node_colours[
    edge_sources
  ]
  
  edge_colours <- grDevices::adjustcolor(
    edge_colours,
    alpha.f = 0.45
  )
  
  
  ##########################################################
  # Draw network WITHOUT igraph labels
  ##########################################################
  
  plot(
    graph,
    
    layout = plot_layout,
    
    vertex.color =
      node_colours[
        igraph::V(graph)$name
      ],
    
    vertex.frame.color = "grey30",
    
    vertex.size = 18,
    
    vertex.label = NA,
    
    edge.width = edge_widths,
    
    edge.color = edge_colours,
    
    edge.arrow.size = 0.08,
    
    edge.arrow.width = 0.6,
    
    edge.curved = 0.20,
    
    asp = 1,
    
    margin = 0.85,
    
    main = title_text
  )
  
  
  ##########################################################
  # MANUAL outward label positions
  ##########################################################
  
  label_scale <- 1.28
  
  label_x <- plot_layout[, 1] * label_scale
  label_y <- plot_layout[, 2] * label_scale
  
  
  ##########################################################
  # Manual x/y nudges for each label
  ##########################################################
  
  x_nudge <- c(
    "Tumour epithelial cells" =  0.00,
    "CAFs"                    = -0.05,
    "Endothelial cells"       =  0.35,
    "TAMs"                    =  0.00,
    "Monocytes"               =  0.00,
    "T/NK cells"              =  0.00,
    "B/plasma cells"          = -0.30,
    "Mast cells"              = -0.05
  )
  
  y_nudge <- c(
    "Tumour epithelial cells" = -0.05,
    "CAFs"                    = -0.05,
    "Endothelial cells"       =  0.00,
    "TAMs"                    = -0.05,
    "Monocytes"               =  0.00,
    "T/NK cells"              = -0.05,
    "B/plasma cells"          =  0.00,
    "Mast cells"              = -0.05
  )
  
  label_x <- label_x + x_nudge[
    igraph::V(graph)$name
  ]
  
  label_y <- label_y + y_nudge[
    igraph::V(graph)$name
  ]
  
  
  ##########################################################
  # Horizontal alignment
  ##########################################################
  
  label_adj_x <- ifelse(
    label_x > 0.20,
    0,
    ifelse(
      label_x < -0.20,
      1,
      0.5
    )
  )
  
  
  ##########################################################
  # Vertical alignment
  ##########################################################
  
  label_adj_y <- ifelse(
    label_y > 0.70,
    0,
    ifelse(
      label_y < -0.70,
      1,
      0.5
    )
  )
  
  
  ##########################################################
  # Add labels manually
  ##########################################################
  
  graphics::text(
    x = label_x,
    y = label_y,
    
    labels = igraph::V(graph)$name,
    
    cex = 0.68,
    
    col = "black",
    
    adj = cbind(
      label_adj_x,
      label_adj_y
    ),
    
    xpd = TRUE
  )
}

par(
  mfrow = c(1, 2),
  xpd = TRUE
)

plot_liana_network_clean(
  primary_graph_no_self,
  "Primary PDAC"
)

plot_liana_network_clean(
  met_graph_no_self,
  "Metastatic PDAC"
)

par(
  mfrow = c(1, 1)
)
OUTDIR <- "/Users/thirisantracy/Desktop/thesis/images"

dir.create(
  OUTDIR,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# PRIMARY — no self signalling
############################################################

quartz(
  type = "pdf",
  file = file.path(
    OUTDIR,
    "LIANA_primary_network_no_self.pdf"
  ),
  width = 6,
  height = 6
)

plot_liana_network_clean(
  primary_graph_no_self,
  "Primary PDAC"
)

dev.off()


############################################################
# METASTATIC — no self signalling
############################################################

quartz(
  type = "pdf",
  file = file.path(
    OUTDIR,
    "LIANA_metastatic_network_no_self.pdf"
  ),
  width = 6,
  height = 6
)

plot_liana_network_clean(
  met_graph_no_self,
  "Metastatic PDAC"
)

dev.off()

############################################################
# Exact interaction key
############################################################

interaction_key <- c(
  "source",
  "target",
  "ligand.complex",
  "receptor.complex"
)


############################################################
# Interactions retained in metastasis but not primary
############################################################

met_specific_interactions <- met_sig %>%
  
  dplyr::anti_join(
    primary_sig,
    by = interaction_key
  ) %>%
  
  dplyr::distinct(
    source,
    target,
    ligand.complex,
    receptor.complex,
    .keep_all = TRUE
  )

met_specific_interactions_no_self <- met_specific_interactions %>%
  
  dplyr::filter(
    source != target
  )

nrow(met_specific_interactions)

nrow(met_specific_interactions_no_self)

############################################################
# Metastasis-specific interactions by sender cell type
############################################################

met_specific_by_sender <- met_specific_interactions_no_self %>%
  
  dplyr::group_by(
    source
  ) %>%
  
  dplyr::summarise(
    
    n_LR_interactions = dplyr::n(),
    
    n_distinct_ligands =
      dplyr::n_distinct(ligand.complex),
    
    n_distinct_receptors =
      dplyr::n_distinct(receptor.complex),
    
    n_receiver_celltypes =
      dplyr::n_distinct(target),
    
    receiver_celltypes =
      paste(
        sort(unique(target)),
        collapse = "; "
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_LR_interactions)
  )
met_specific_by_sender %>%
  print(
    n = Inf,
    width = Inf
  )


############################################################
# Metastasis-specific interactions by sender → receiver
############################################################

met_specific_sender_receiver <- met_specific_interactions_no_self %>%
  
  dplyr::group_by(
    source,
    target
  ) %>%
  
  dplyr::summarise(
    
    n_LR_interactions =
      dplyr::n(),
    
    n_distinct_ligands =
      dplyr::n_distinct(ligand.complex),
    
    n_distinct_receptors =
      dplyr::n_distinct(receptor.complex),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    source,
    dplyr::desc(n_LR_interactions)
  )

met_specific_sender_receiver %>%
  print(
    n = Inf,
    width = Inf
  )

caf_met_all <- met_sig %>%
  
  dplyr::filter(
    source == "CAFs",
    source != target
  ) %>%
  
  dplyr::distinct(
    source,
    target,
    ligand.complex,
    receptor.complex,
    .keep_all = TRUE
  )

sort(
  unique(caf_met_all$target)
)

caf_met_ligand_summary <- caf_met_all %>%
  
  dplyr::group_by(
    ligand.complex
  ) %>%
  
  dplyr::summarise(
    
    n_receiver_celltypes =
      dplyr::n_distinct(target),
    
    n_receptors =
      dplyr::n_distinct(receptor.complex),
    
    n_LR_contexts =
      dplyr::n_distinct(
        paste(
          target,
          receptor.complex,
          sep = " | "
        )
      ),
    
    receiver_celltypes =
      paste(
        sort(unique(target)),
        collapse = "; "
      ),
    
    receptors =
      paste(
        sort(unique(receptor.complex)),
        collapse = "; "
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_receiver_celltypes),
    dplyr::desc(n_receptors),
    dplyr::desc(n_LR_contexts)
  )

caf_met_ligand_summary %>%
  print(
    n = Inf,
    width = Inf
  )

caf_met_all %>%
  
  dplyr::select(
    ligand.complex,
    receptor.complex,
    target,
    aggregate_rank
  ) %>%
  
  dplyr::arrange(
    ligand.complex,
    target,
    receptor.complex
  ) %>%
  
  print(
    n = Inf,
    width = Inf
  )

caf_met_specific <- met_specific_interactions_no_self %>%
  dplyr::filter(
    source == "CAFs"
  )

nrow(caf_met_specific)

dplyr::n_distinct(
  caf_met_specific$ligand.complex
)

sort(
  unique(caf_met_specific$target)
)

caf_met_specific_ligands <- caf_met_specific %>%
  
  dplyr::group_by(
    ligand.complex
  ) %>%
  
  dplyr::summarise(
    
    n_receiver_celltypes =
      dplyr::n_distinct(target),
    
    n_receptors =
      dplyr::n_distinct(receptor.complex),
    
    n_LR_contexts =
      dplyr::n_distinct(
        paste(
          target,
          receptor.complex,
          sep = " | "
        )
      ),
    
    receiver_celltypes =
      paste(
        sort(unique(target)),
        collapse = "; "
      ),
    
    receptors =
      paste(
        sort(unique(receptor.complex)),
        collapse = "; "
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_receiver_celltypes),
    dplyr::desc(n_receptors),
    dplyr::desc(n_LR_contexts)
  )

caf_met_specific_ligands %>%
  print(
    n = Inf,
    width = Inf
  )

caf_recurrent <- caf_met_specific_ligands %>%
  dplyr::filter(
    n_receiver_celltypes >= 2
  )

caf_recurrent_names <- caf_recurrent$ligand.complex


caf_recurrent_targets <- caf_met_specific %>%
  
  dplyr::filter(
    ligand.complex %in% caf_recurrent_names
  ) %>%
  
  dplyr::distinct(
    ligand.complex,
    target
  ) %>%
  
  dplyr::arrange(
    ligand.complex,
    target
  )

caf_recurrent_targets %>%
  print(
    n = Inf,
    width = Inf
  )

caf_recurrent_LR <- caf_met_specific %>%
  
  dplyr::filter(
    ligand.complex %in% caf_recurrent_names
  ) %>%
  
  dplyr::select(
    ligand.complex,
    target,
    receptor.complex,
    aggregate_rank
  ) %>%
  
  dplyr::arrange(
    ligand.complex,
    target,
    receptor.complex
  )

caf_recurrent_LR %>%
  print(
    n = Inf,
    width = Inf
  )


############################################################
# CAF metastasis-specific ligand UpSet plot
#
# Input:
# caf_met_specific =
# exact CAF-derived LR interactions retained only in
# metastatic LIANA analysis
#
# All 32 ligands are included.
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)


############################################################
# 1. Receiver-cell order
############################################################

receiver_order <- c(
  "Tumour epithelial cells",
  "Endothelial cells",
  "TAMs",
  "Monocytes",
  "Mast cells",
  "T/NK cells",
  "B/plasma cells"
)


############################################################
# 2. Ligand × receiver membership
#
# One row = one ligand associated with one receiver cell type
############################################################

caf_upset_long <- caf_met_specific %>%
  
  dplyr::distinct(
    ligand.complex,
    target
  )


# Check number of ligands
dplyr::n_distinct(
  caf_upset_long$ligand.complex
)


############################################################
# 3. Create receiver combination for every ligand
############################################################

caf_upset_membership <- caf_upset_long %>%
  
  dplyr::group_by(
    ligand.complex
  ) %>%
  
  dplyr::summarise(
    
    receiver_combination =
      paste(
        receiver_order[
          receiver_order %in% unique(target)
        ],
        collapse = " + "
      ),
    
    n_receivers =
      dplyr::n_distinct(target),
    
    .groups = "drop"
  )


############################################################
# 4. Build unique UpSet intersections
#
# Ordered first by number of receiver cell types,
# then by intersection size
############################################################

caf_upset_intersections <- caf_upset_membership %>%
  
  dplyr::group_by(
    receiver_combination,
    n_receivers
  ) %>%
  
  dplyr::summarise(
    
    intersection_size =
      dplyr::n(),
    
    ligands =
      paste(
        sort(ligand.complex),
        collapse = "\n"
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_receivers),
    dplyr::desc(intersection_size),
    receiver_combination
  ) %>%
  
  dplyr::mutate(
    intersection_id =
      paste0(
        "I",
        dplyr::row_number()
      )
  )


caf_upset_intersections %>%
  print(
    n = Inf,
    width = Inf
  )


############################################################
# 5. Fix intersection order ONCE
############################################################

intersection_order <- as.character(
  caf_upset_intersections$intersection_id
)

caf_upset_intersections$intersection_id <-
  factor(
    caf_upset_intersections$intersection_id,
    levels = intersection_order
  )


############################################################
# 6. Build dot-matrix data
############################################################

caf_upset_matrix <- tidyr::crossing(
  
  intersection_id =
    intersection_order,
  
  receiver =
    receiver_order
) %>%
  
  dplyr::left_join(
    
    caf_upset_intersections %>%
      dplyr::mutate(
        intersection_id =
          as.character(intersection_id)
      ) %>%
      dplyr::select(
        intersection_id,
        receiver_combination
      ),
    
    by = "intersection_id"
  ) %>%
  
  dplyr::rowwise() %>%
  
  dplyr::mutate(
    
    active =
      receiver %in%
      strsplit(
        receiver_combination,
        " + ",
        fixed = TRUE
      )[[1]]
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(
    
    intersection_id =
      factor(
        intersection_id,
        levels = intersection_order
      ),
    
    receiver =
      factor(
        receiver,
        levels = rev(receiver_order)
      )
  )


############################################################
# 7. Set sizes
#
# Number of DISTINCT CAF ligands associated with each receiver
############################################################

caf_set_sizes <- caf_upset_long %>%
  
  dplyr::distinct(
    ligand.complex,
    target
  ) %>%
  
  dplyr::count(
    target,
    name = "set_size"
  ) %>%
  
  tidyr::complete(
    target = receiver_order,
    fill = list(
      set_size = 0
    )
  ) %>%
  
  dplyr::mutate(
    
    receiver =
      factor(
        target,
        levels = rev(receiver_order)
      )
  )


caf_set_sizes


############################################################
# 8. TOP PANEL
# Intersection-size bars + ligand names
############################################################

caf_upset_bar <- ggplot(
  
  caf_upset_intersections,
  
  aes(
    x = intersection_id,
    y = intersection_size
  )
  
) +
  
  geom_col(
    width = 0.65
  ) +
  
  # Ligand names above each intersection
  geom_text(
    aes(
      y = intersection_size + 0.30,
      label = ligands
    ),
    vjust = 0,
    size = 2.5,
    lineheight = 0.90
  ) +
  
  ##########################################################
# IMPORTANT:
# Explicit x scale = identical to matrix x scale
##########################################################

scale_x_discrete(
  limits = intersection_order,
  drop = FALSE,
  expand = expansion(
    add = 0.6
  )
) +
  
  scale_y_continuous(
    
    breaks =
      0:max(
        caf_upset_intersections$intersection_size
      ),
    
    expand = expansion(
      add = c(
        0.18,
        1.8
      )
    )
  ) +
  
  labs(
    x = NULL,
    y = "Intersection size"
  ) +
  
  theme_classic() +
  
  theme(
    
    axis.text.x =
      element_blank(),
    
    axis.ticks.x =
      element_blank(),
    
    axis.line.x =
      element_blank(),
    
    axis.title.y =
      element_text(
        size = 11
      ),
    
    axis.text.y =
      element_text(
        size = 9
      ),
    
    plot.margin =
      margin(
        t = 10,
        r = 10,
        b = 0,
        l = 10
      )
  )


############################################################
# 9. RIGHT-BOTTOM PANEL
# UpSet dot matrix
############################################################

caf_upset_matrix_plot <- ggplot(
  
  caf_upset_matrix,
  
  aes(
    x = intersection_id,
    y = receiver
  )
  
) +
  
  ##########################################################
# Inactive dots
##########################################################

geom_point(
  size = 2.6,
  alpha = 0.10
) +
  
  
  ##########################################################
# Connecting lines
##########################################################

geom_line(
  
  data =
    caf_upset_matrix %>%
    dplyr::filter(active),
  
  aes(
    group = intersection_id
  ),
  
  linewidth = 0.6
) +
  
  
  ##########################################################
# Active dots
##########################################################

geom_point(
  
  data =
    caf_upset_matrix %>%
    dplyr::filter(active),
  
  size = 3.2
) +
  
  
  ##########################################################
# IMPORTANT:
# EXACT same x scale as intersection bars
##########################################################

scale_x_discrete(
  limits = intersection_order,
  drop = FALSE,
  expand = expansion(
    add = 0.6
  )
) +
  
  
  ##########################################################
# IMPORTANT:
# EXACT same y scale as set-size plot
##########################################################

scale_y_discrete(
  limits = rev(receiver_order),
  drop = FALSE,
  expand = expansion(
    add = 0.6
  )
) +
  
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic() +
  
  theme(
    
    axis.text.x =
      element_blank(),
    
    axis.ticks.x =
      element_blank(),
    
    axis.line.x =
      element_blank(),
    
    # Cell-type labels are shown in set-size plot,
    # so don't repeat them here
    axis.text.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    axis.line.y =
      element_blank(),
    
    plot.margin =
      margin(
        t = 0,
        r = 10,
        b = 10,
        l = 0
      )
  )


############################################################
# 10. LEFT-BOTTOM PANEL
# Set-size bars
############################################################

caf_set_size_plot <- ggplot(
  
  caf_set_sizes,
  
  aes(
    x = set_size,
    y = receiver
  )
  
) +
  
  geom_col(
    width = 0.65
  ) +
  
  
  ##########################################################
# Integer set-size axis
##########################################################

scale_x_reverse(
  
  breaks =
    seq(
      0,
      max(caf_set_sizes$set_size),
      by = 2
    ),
  
  expand = expansion(
    mult = c(
      0.05,
      0
    )
  )
) +
  
  
  ##########################################################
# IMPORTANT:
# EXACT same y scale as dot matrix
##########################################################

scale_y_discrete(
  limits = rev(receiver_order),
  position = "right",
  drop = FALSE,
  expand = expansion(
    add = 0.6
  )
) +
  
  
  labs(
    x = "Set size",
    y = NULL
  ) +
  
  theme_classic() +
  
  theme(
    
    axis.line.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    axis.text.y =
      element_text(
        size = 10,
        margin = margin(
          l = 5,
          r = 5
        )
      ),
    
    axis.title.x =
      element_text(
        size = 10
      ),
    
    axis.text.x =
      element_text(
        size = 9
      ),
    
    plot.margin =
      margin(
        t = 0,
        r = 0,
        b = 10,
        l = 10
      )
  )


############################################################
# 11. ALIGN TOP BARS WITH DOT-MATRIX COLUMNS
############################################################

right_aligned <- cowplot::align_plots(
  
  caf_upset_bar,
  caf_upset_matrix_plot,
  
  align = "v",
  
  axis = "lr"
)

caf_upset_bar_aligned <-
  right_aligned[[1]]

caf_upset_matrix_xaligned <-
  right_aligned[[2]]


############################################################
# 12. ALIGN SET-SIZE BARS WITH DOT-MATRIX ROWS
############################################################

bottom_aligned <- cowplot::align_plots(
  
  caf_set_size_plot,
  caf_upset_matrix_xaligned,
  
  align = "h",
  
  axis = "tb"
)

caf_set_size_aligned <-
  bottom_aligned[[1]]

caf_upset_matrix_aligned <-
  bottom_aligned[[2]]


############################################################
# 13. Build top row
#
# Blank area occupies same width as left set-size panel
############################################################

caf_upset_top <- cowplot::plot_grid(
  
  NULL,
  caf_upset_bar_aligned,
  
  nrow = 1,
  
  rel_widths = c(
    0.28,
    0.72
  )
)


############################################################
# 14. Build bottom row
############################################################

caf_upset_bottom <- cowplot::plot_grid(
  
  caf_set_size_aligned,
  caf_upset_matrix_aligned,
  
  nrow = 1,
  
  rel_widths = c(
    0.28,
    0.72
  )
)


############################################################
# 15. FINAL UpSet figure
############################################################

caf_met_specific_upset <- cowplot::plot_grid(
  
  caf_upset_top,
  caf_upset_bottom,
  
  ncol = 1,
  
  rel_heights = c(
    1.45,
    1
  )
)


############################################################
# Preview
############################################################

caf_met_specific_upset

ggsave(
  filename = file.path(
    OUTDIR,
    "LIANA_CAF_met_specific_ligand_UpSet.pdf"
  ),
  plot = caf_met_specific_upset,
  width = 14,
  height = 7,
  units = "in",
  device = "pdf",
  useDingbats = FALSE,
  bg = "white"
)

library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)

############################################################
# 1. Output folder
############################################################

OUTDIR <- "/Users/thirisantracy/Desktop/thesis/images"

dir.create(
  OUTDIR,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# 2. CAF targets to plot
############################################################

caf_targets <- c(
  "Tumour epithelial cells",
  "Endothelial cells",
  "Mast cells",
  "Monocytes",
  "T/NK cells",
  "TAMs",
  "B/plasma cells"
)

############################################################
# 3. Prepare ALL CAF metastasis-specific exact LR pairs
#    that will be used across the 6 plots
############################################################

caf_all_plot_df <- caf_met_specific %>%
  
  dplyr::filter(
    source == "CAFs",
    target %in% caf_targets
  ) %>%
  
  dplyr::mutate(
    interaction = paste(
      ligand.complex,
      receptor.complex,
      sep = " \u2192 "
    )
  )

############################################################
############ cell interactions only in met

############################################################
# Selected sender-receiver pairs to inspect
############################################################

selected_pairs <- tibble::tribble(
  ~source,                  ~target,
  
  # Existing metastasis-only pairs already plotted
  "CAFs",                   "Endothelial cells",
  "CAFs",                   "Tumour epithelial cells",
  "TAMs",                   "Endothelial cells",
  "CAFs",                   "TAMs",
  "Endothelial cells",      "CAFs",
  "CAFs",                   "T/NK cells",
  "Endothelial cells",      "T/NK cells",
  "B/plasma cells",         "Monocytes",
  "Endothelial cells",      "TAMs",
  
  # Newly noticed metastasis-only pairs
  "Mast cells",             "CAFs",
  "Tumour epithelial cells","Monocytes",
  "Mast cells",             "TAMs"
)

############################################################
# Rebuild comparison object for ONLY the selected pairs
############################################################

liana_focus <- liana_comparison %>%
  dplyr::inner_join(
    selected_pairs,
    by = c("source", "target")
  )

############################################################
# Keep metastasis-only exact LR pairs
############################################################

liana_met_only <- liana_focus %>%
  dplyr::filter(
    !is.na(rank_metastasis),
    is.na(rank_primary)
  ) %>%
  dplyr::mutate(
    LR_pair = paste(
      ligand.complex,
      receptor.complex,
      sep = " → "
    )
  )

# Prepare the 4 additional metastasis-only plotting datasets
# using the SAME column names expected by make_lr_plot()
############################################################

mast_to_caf_plot_df <- liana_met_only %>%
  dplyr::filter(
    source == "Mast cells",
    target == "CAFs"
  ) %>%
  dplyr::transmute(
    interaction = LR_pair,
    natmi.edge_specificity = specificity_metastasis,
    sca.LRscore = magnitude_metastasis
  )


bplasma_to_monocytes_plot_df <- liana_met_only %>%
  dplyr::filter(
    source == "B/plasma cells",
    target == "Monocytes"
  ) %>%
  dplyr::transmute(
    interaction = LR_pair,
    natmi.edge_specificity = specificity_metastasis,
    sca.LRscore = magnitude_metastasis
  )


tumourepi_to_monocytes_plot_df <- liana_met_only %>%
  dplyr::filter(
    source == "Tumour epithelial cells",
    target == "Monocytes"
  ) %>%
  dplyr::transmute(
    interaction = LR_pair,
    natmi.edge_specificity = specificity_metastasis,
    sca.LRscore = magnitude_metastasis
  )


mast_to_tams_plot_df <- liana_met_only %>%
  dplyr::filter(
    source == "Mast cells",
    target == "TAMs"
  ) %>%
  dplyr::transmute(
    interaction = LR_pair,
    natmi.edge_specificity = specificity_metastasis,
    sca.LRscore = magnitude_metastasis
  )
combined_specificity_limits <- c(
  min(caf_all_plot_df$natmi.edge_specificity,
      mast_to_caf_plot_df$natmi.edge_specificity,
      bplasma_to_monocytes_plot_df$natmi.edge_specificity,
      tumourepi_to_monocytes_plot_df$natmi.edge_specificity,
      mast_to_tams_plot_df$natmi.edge_specificity),
  max(caf_all_plot_df$natmi.edge_specificity,
      mast_to_caf_plot_df$natmi.edge_specificity,
      bplasma_to_monocytes_plot_df$natmi.edge_specificity,
      tumourepi_to_monocytes_plot_df$natmi.edge_specificity,
      mast_to_tams_plot_df$natmi.edge_specificity)
)

combined_specificity_breaks <- c(0.2, 0.4, 0.6, 0.8)

combined_magnitude_limits <- c(
  min(caf_all_plot_df$sca.LRscore,
      mast_to_caf_plot_df$sca.LRscore,
      bplasma_to_monocytes_plot_df$sca.LRscore,
      tumourepi_to_monocytes_plot_df$sca.LRscore,
      mast_to_tams_plot_df$sca.LRscore),
  max(caf_all_plot_df$sca.LRscore,
      mast_to_caf_plot_df$sca.LRscore,
      bplasma_to_monocytes_plot_df$sca.LRscore,
      tumourepi_to_monocytes_plot_df$sca.LRscore,
      mast_to_tams_plot_df$sca.LRscore)
)

combined_magnitude_breaks <- c(0.90, 0.95)
library(ggplot2)
library(dplyr)
library(cowplot)
library(stringr)
library(grid)

############################################################
# Clean legend breaks
############################################################

combined_specificity_breaks <- c(0.2, 0.4, 0.6, 0.8)
combined_magnitude_breaks   <- c(0.90, 0.95)

############################################################
# Force all strip titles to have the same height
# so panel boxes stay visually uniform
############################################################

format_source_label <- function(x) {
  
  dplyr::case_when(
    x == "Tumour epithelial cells" ~ "Tumour epithelial\ncells",
    x == "B/plasma cells"          ~ "B/plasma\ncells",
    TRUE                           ~ paste0(x, "\n")
  )
}

############################################################
# Main plotting function
############################################################

make_lr_plot <- function(
    plot_df,
    source_name,
    target_name,
    specificity_limits,
    specificity_breaks,
    magnitude_limits,
    magnitude_breaks,
    show_legend = FALSE
) {
  
  source_label <- format_source_label(source_name)
  
  plot_df <- plot_df %>%
    
    dplyr::arrange(
      natmi.edge_specificity,
      dplyr::desc(sca.LRscore)
    ) %>%
    
    dplyr::mutate(
      
      interaction = factor(
        interaction,
        levels = unique(interaction)
      ),
      
      ######################################################
      # More vertical spacing between rows so dots do not clash
      ######################################################
      y_pos = rev(
        seq(
          from = 1,
          by = 2,
          length.out = dplyr::n()
        )
      ),
      
      source = factor(
        source_label,
        levels = source_label
      ),
      
      target = factor(
        target_name,
        levels = target_name
      )
    )
  
  ggplot(
    plot_df,
    aes(
      x = target,
      y = y_pos
    )
  ) +
    
    geom_point(
      aes(
        size   = natmi.edge_specificity,
        colour = sca.LRscore
      ),
      alpha = 0.98
    ) +
    
    facet_grid(
      . ~ source,
      scales = "free_x",
      space  = "free_x"
    ) +
    
    scale_x_discrete(
      expand = expansion(add = c(0.90, 0.90))
    ) +
    
    scale_y_continuous(
      breaks = plot_df$y_pos,
      labels = plot_df$interaction,
      expand = expansion(add = c(1.2, 1.2))
    ) +
    
    scale_size_continuous(
      name   = "Interaction\nspecificity",
      limits = specificity_limits,
      breaks = specificity_breaks,
      labels = sprintf("%.1f", specificity_breaks),
      range  = c(4.0, 11.0)
    ) +
    
    scale_colour_viridis_c(
      name      = "Interaction\nmagnitude",
      limits    = magnitude_limits,
      breaks    = magnitude_breaks,
      labels    = sprintf("%.2f", magnitude_breaks),
      option    = "D",
      direction = -1
    ) +
    
    labs(
      x = NULL,
      y = NULL
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      ######################################################
      # Facet strip
      ######################################################
      strip.background = element_rect(
        fill      = "grey90",
        colour    = "grey40",
        linewidth = 0.5
      ),
      
      strip.text = element_text(
        size      = 9,
        lineheight = 0.90,
        margin    = margin(t = 3, r = 2, b = 3, l = 2)
      ),
      
      ######################################################
      # Axes
      ######################################################
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      axis.text.x = element_text(
        size  = 9,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      
      axis.text.y = element_text(
        size = 8.3
      ),
      
      ######################################################
      # Gridlines
      ######################################################
      panel.grid.major.x = element_line(
        colour    = "grey85",
        linewidth = 0.35
      ),
      
      panel.grid.major.y = element_line(
        colour    = "grey85",
        linewidth = 0.35
      ),
      
      panel.grid.minor = element_blank(),
      
      ######################################################
      # Legends
      ######################################################
      legend.position = if (show_legend) "right" else "none",
      legend.title    = element_text(size = 9),
      legend.text     = element_text(size = 8),
      
      ######################################################
      # Margins
      ######################################################
      plot.margin = margin(
        t = 6, r = 6, b = 6, l = 6
      )
    )
}

############################################################
# Build CAF plots again
############################################################

caf_targets_for_panel <- c(
  "Tumour epithelial cells",
  "Endothelial cells",
  "Mast cells",
  "Monocytes",
  "T/NK cells",
  "TAMs",
  "B/plasma cells"
)

caf_plot_list <- lapply(
  caf_targets_for_panel,
  function(tg) {
    make_lr_plot(
      plot_df = caf_all_plot_df %>%
        dplyr::filter(target == tg),
      source_name = "CAFs",
      target_name = tg,
      specificity_limits = combined_specificity_limits,
      specificity_breaks = combined_specificity_breaks,
      magnitude_limits = combined_magnitude_limits,
      magnitude_breaks = combined_magnitude_breaks,
      show_legend = FALSE
    )
  }
)

names(caf_plot_list) <- caf_targets_for_panel

############################################################
# Build the 4 new metastasis-only plots again
############################################################

new_plot_list <- list(
  
  make_lr_plot(
    plot_df = mast_to_caf_plot_df,
    source_name = "Mast cells",
    target_name = "CAFs",
    specificity_limits = combined_specificity_limits,
    specificity_breaks = combined_specificity_breaks,
    magnitude_limits = combined_magnitude_limits,
    magnitude_breaks = combined_magnitude_breaks,
    show_legend = FALSE
  ),
  
  make_lr_plot(
    plot_df = bplasma_to_monocytes_plot_df,
    source_name = "B/plasma cells",
    target_name = "Monocytes",
    specificity_limits = combined_specificity_limits,
    specificity_breaks = combined_specificity_breaks,
    magnitude_limits = combined_magnitude_limits,
    magnitude_breaks = combined_magnitude_breaks,
    show_legend = FALSE
  ),
  
  make_lr_plot(
    plot_df = tumourepi_to_monocytes_plot_df,
    source_name = "Tumour epithelial cells",
    target_name = "Monocytes",
    specificity_limits = combined_specificity_limits,
    specificity_breaks = combined_specificity_breaks,
    magnitude_limits = combined_magnitude_limits,
    magnitude_breaks = combined_magnitude_breaks,
    show_legend = FALSE
  ),
  
  make_lr_plot(
    plot_df = mast_to_tams_plot_df,
    source_name = "Mast cells",
    target_name = "TAMs",
    specificity_limits = combined_specificity_limits,
    specificity_breaks = combined_specificity_breaks,
    magnitude_limits = combined_magnitude_limits,
    magnitude_breaks = combined_magnitude_breaks,
    show_legend = FALSE
  )
)

############################################################
# One shared legend
############################################################

legend_plot <- make_lr_plot(
  plot_df = caf_all_plot_df %>%
    dplyr::filter(target == "Tumour epithelial cells"),
  source_name = "CAFs",
  target_name = "Tumour epithelial cells",
  specificity_limits = combined_specificity_limits,
  specificity_breaks = combined_specificity_breaks,
  magnitude_limits = combined_magnitude_limits,
  magnitude_breaks = combined_magnitude_breaks,
  show_legend = TRUE
)

shared_legend <- cowplot::get_legend(
  legend_plot +
    theme(
      legend.position = "right"
    )
)

############################################################
# VERY IMPORTANT:
# Align ALL real plots first, not row-by-row
############################################################

all_real_plots <- c(
  caf_plot_list,
  new_plot_list
)

aligned_plots <- cowplot::align_plots(
  plotlist = all_real_plots,
  align = "hv",
  axis  = "tblr"
)

caf_aligned <- aligned_plots[1:6]
new_aligned <- aligned_plots[7:10]

############################################################
# Blank placeholders for row 2 right-hand side
############################################################

blank_plot <- ggplot() +
  theme_void()

############################################################
# Build one single 3 x 4 grid
# Row 1 = four CAF plots
# Row 2 = two CAF plots left-aligned + two blanks
# Row 3 = four met-only plots
############################################################

main_panel <- cowplot::plot_grid(
  
  caf_aligned[[1]], caf_aligned[[2]], caf_aligned[[3]], caf_aligned[[4]],
  caf_aligned[[5]], caf_aligned[[6]], blank_plot,        blank_plot,
  new_aligned[[1]], new_aligned[[2]], new_aligned[[3]],  new_aligned[[4]],
  
  ncol = 4,
  
  align = "hv",
  axis  = "tblr",
  
  rel_widths  = c(1, 1, 1, 1),
  rel_heights = c(1, 1, 1)
)

############################################################
# Add shared legend on the right
############################################################

final_combined_panel <- cowplot::plot_grid(
  main_panel,
  shared_legend,
  nrow = 1,
  rel_widths = c(0.90, 0.10)
)

final_combined_panel

############################################################
# Build one single 2 x 5 grid
#
# Row 1 = first 5 CAF plots
# Row 2 = last CAF plot + 4 met-only plots
############################################################

main_panel <- cowplot::plot_grid(
  
  caf_aligned[[1]], caf_aligned[[2]], caf_aligned[[3]], caf_aligned[[4]], caf_aligned[[5]],
  caf_aligned[[6]], new_aligned[[1]], new_aligned[[2]], new_aligned[[3]], new_aligned[[4]],
  
  ncol = 5,
  
  align = "hv",
  axis  = "tblr",
  
  rel_widths  = c(1, 1, 1, 1, 1),
  rel_heights = c(1, 1)
)

############################################################
# Add shared legend on the right
############################################################

final_combined_panel <- cowplot::plot_grid(
  main_panel,
  shared_legend,
  nrow = 1,
  rel_widths = c(0.90, 0.10)
)

final_combined_panel

ggsave(
  filename = file.path(
    OUTDIR,
    "LIANA_metastasis_specific_LR_panel_2rows5cols_aligned.pdf"
  ),
  plot = final_combined_panel,
  width = 15,
  height = 10,
  units = "in",
  device = "pdf",
  useDingbats = FALSE,
  bg = "white"
)



shared_legend <- cowplot::get_legend(
  
  make_lr_plot(
    plot_df = caf_all_plot_df %>%
      dplyr::filter(
        target == "Tumour epithelial cells"
      ),
    
    source_name = "CAFs",
    target_name = "Tumour epithelial cells",
    
    specificity_limits = combined_specificity_limits,
    specificity_breaks = combined_specificity_breaks,
    
    magnitude_limits = combined_magnitude_limits,
    magnitude_breaks = combined_magnitude_breaks,
    
    show_legend = TRUE
  ) +
    
    theme(
      legend.position = "right",
      legend.box = "vertical"
    )
)


############################################################
# All 11 real plots
############################################################

all_real_plots <- c(
  
  unname(
    caf_plot_list[
      c(
        "Tumour epithelial cells",
        "Endothelial cells",
        "Mast cells",
        "Monocytes",
        "T/NK cells",
        "TAMs",
        "B/plasma cells"
      )
    ]
  ),
  
  new_plot_list
)

############################################################
# Align them BEFORE making the grid
############################################################

aligned_plots <- cowplot::align_plots(
  plotlist = all_real_plots,
  align = "hv",
  axis = "tblr"
)

caf_aligned <- aligned_plots[1:7]
new_aligned <- aligned_plots[8:11]


############################################################
# FINAL 2 ROW × 6 COLUMN PANEL
############################################################

final_combined_panel <- cowplot::plot_grid(
  
  # ROW 1 — six CAF receiver populations
  caf_aligned[[1]],
  caf_aligned[[2]],
  caf_aligned[[3]],
  caf_aligned[[4]],
  caf_aligned[[5]],
  caf_aligned[[6]],
  
  # ROW 2 — final CAF target + four additional pairs + legend
  caf_aligned[[7]],
  new_aligned[[1]],
  new_aligned[[2]],
  new_aligned[[3]],
  new_aligned[[4]],
  shared_legend,
  
  ncol = 6,
  nrow = 2,
  
  rel_widths = c(
    1, 1, 1, 1, 1, 1
  ),
  
  rel_heights = c(
    1, 1
  ),
  
  align = "hv",
  axis = "tblr"
)

final_combined_panel

ggsave(
  filename = file.path(
    OUTDIR,
    "LIANA_metastasis_specific_LR_panel_2rows6cols.pdf"
  ),
  plot = final_combined_panel,
  width = 19,
  height = 10.5,
  units = "in",
  device = "pdf",
  useDingbats = FALSE,
  bg = "white"
)
