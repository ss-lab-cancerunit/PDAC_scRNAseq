#### cell to cell communication with lianna
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


key<-c("source","target","ligand.complex","receptor.complex")



shared<-dplyr::inner_join(
  primary_sig,
  met_sig,
  by=key,
  suffix=c("_primary","_met")
)

primary_only<-dplyr::anti_join(
  primary_sig,
  met_sig,
  by=key
)

met_only<-dplyr::anti_join(
  met_sig,
  primary_sig,
  by=key
)


nrow(shared)
nrow(primary_only)
nrow(met_only)

pair_comparison<-dplyr::bind_rows(
  shared%>%dplyr::count(source,target,name="n")%>%dplyr::mutate(group="Shared"),
  primary_only%>%dplyr::count(source,target,name="n")%>%dplyr::mutate(group="Primary"),
  met_only%>%dplyr::count(source,target,name="n")%>%dplyr::mutate(group="Metastasis")
)%>%
  tidyr::pivot_wider(names_from=group,values_from=n,values_fill=0)

pair_comparison%>%
  dplyr::arrange(dplyr::desc(Primary+Metastasis))


pair_long<-pair_comparison%>%
  tidyr::pivot_longer(c(Shared,Primary,Metastasis),names_to="group",values_to="n")%>%
  dplyr::mutate(pair=paste(source,"→",target))

library(ggplot2)

p_pairs<-ggplot(pair_long,aes(x=group,y=pair,fill=n))+
  geom_tile(color="white")+
  geom_text(aes(label=n),size=3)+
  labs(x=NULL,y=NULL,fill="Interactions")+
  theme_bw()+
  theme(
    axis.text.y=element_text(size=8),
    panel.grid=element_blank()
  )

p_pairs

all_long<-dplyr::bind_rows(
  primary_sig%>%dplyr::mutate(condition="Primary"),
  met_sig%>%dplyr::mutate(condition="Metastasis")
)%>%
  dplyr::mutate(
    LR_pair=paste(ligand.complex,receptor.complex,sep=" → "),
    condition=factor(condition,levels=c("Primary","Metastasis"))
  )


make_liana_pair_plot<-function(source_cell,target_cell){
  
  plot_df<-all_long%>%
    dplyr::filter(source==source_cell,target==target_cell)%>%
    dplyr::mutate(
      target_condition=paste(target,condition,sep=" | "),
      target=target_condition
    )
  
  target_order<-c(
    paste(target_cell,"Primary",sep=" | "),
    paste(target_cell,"Metastasis",sep=" | ")
  )
  
  liana::liana_dotplot(
    plot_df,
    source_groups=source_cell,
    target_groups=target_order,
    magnitude="sca.LRscore",
    specificity="natmi.edge_specificity",
    ntop=Inf
  )+
    labs(
      title=paste(source_cell,"→",target_cell),
      x=NULL,
      y="Ligand-receptor interaction"
    )+
    theme(
      axis.text.x=element_text(angle=45,hjust=1),
      axis.text.y=element_text(size=6)
    )
}

pairs<-all_long%>%
  dplyr::distinct(source,target)%>%
  dplyr::arrange(source,target)


dir.create(file.path(ROOT,"pair_dotplots"),showWarnings=FALSE)

for(i in seq_len(nrow(pairs))){
  
  source_cell<-pairs$source[i]
  target_cell<-pairs$target[i]
  
  p<-make_liana_pair_plot(source_cell,target_cell)
  
  n_lr<-all_long%>%
    dplyr::filter(source==source_cell,target==target_cell)%>%
    dplyr::distinct(ligand.complex,receptor.complex)%>%
    nrow()
  
  safe_source<-gsub("[^A-Za-z0-9]+","_",source_cell)
  safe_target<-gsub("[^A-Za-z0-9]+","_",target_cell)
  
  ggsave(
    file.path(ROOT,"pair_dotplots",paste0(safe_source,"_to_",safe_target,".pdf")),
    p,
    width=8,
    height=max(6,n_lr*0.22)
  )
}



all_long_spec<-all_long%>%
  dplyr::filter(natmi.edge_specificity>0.5)


nrow(all_long)
nrow(all_long_spec)

dplyr::n_distinct(all_long_spec$LR_pair)


pairs_spec<-all_long_spec%>%
  dplyr::distinct(source,target)%>%
  dplyr::arrange(source,target)
nrow(pairs_spec)


make_liana_spec_plot<-function(source_cell,target_cell){
  
  plot_df<-all_long_spec%>%
    dplyr::filter(source==source_cell,target==target_cell)%>%
    dplyr::mutate(
      target_condition=paste(target,condition,sep=" | "),
      target=target_condition
    )
  
  target_order<-c(
    paste(target_cell,"Primary",sep=" | "),
    paste(target_cell,"Metastasis",sep=" | ")
  )
  
  target_order<-target_order[target_order%in%plot_df$target]
  
  liana::liana_dotplot(
    plot_df,
    source_groups=source_cell,
    target_groups=target_order,
    magnitude="sca.LRscore",
    specificity="natmi.edge_specificity",
    ntop=Inf
  )+
    labs(
      title=paste(source_cell,"→",target_cell),
      x=NULL,
      y="Ligand-receptor interaction"
    )+
    theme(
      plot.title=element_text(size=9,face="bold",hjust=0.5),
      axis.title.y=element_text(size=8),
      axis.text.x=element_text(size=6,angle=90,hjust=1,vjust=0.5),
      axis.text.y=element_text(size=9),
      strip.text=element_text(size=7),
      legend.title=element_text(size=6),
      legend.text=element_text(size=6)
    )
}
dir.create(file.path(ROOT,"pair_dotplots_specificity05"),showWarnings=FALSE)


for(i in seq_len(nrow(pairs_spec))){
  
  source_cell<-pairs_spec$source[i]
  target_cell<-pairs_spec$target[i]
  
  p<-make_liana_spec_plot(source_cell,target_cell)
  
  n_lr<-all_long_spec%>%
    dplyr::filter(source==source_cell,target==target_cell)%>%
    dplyr::distinct(ligand.complex,receptor.complex)%>%
    nrow()
  
  safe_source<-gsub("[^A-Za-z0-9]+","_",source_cell)
  safe_target<-gsub("[^A-Za-z0-9]+","_",target_cell)
  
  ggsave(
    file.path(
      ROOT,
      "pair_dotplots_specificity05",
      paste0(safe_source,"_to_",safe_target,"_specificity_gt_0.5.pdf")
    ),
    p,
    width=7,
    height=max(6,n_lr*0.3)
  )
}

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

selected_pairs <- tibble::tribble(
  ~source,              ~target,
  
  # Existing
  "CAFs",               "Endothelial cells",
  "CAFs",               "Tumour epithelial cells",
  "TAMs",               "Endothelial cells",
  "CAFs",               "TAMs",
  "Endothelial cells",  "CAFs",
  
  # New
  "CAFs",               "T/NK cells",
  "Endothelial cells",  "T/NK cells",
  "B/plasma cells",     "Monocytes",
  "Endothelial cells",  "TAMs"
)

liana_focus <- liana_comparison %>%
  dplyr::inner_join(
    selected_pairs,
    by = c("source", "target")
  )


liana_plot_data <- liana_focus %>%
  dplyr::mutate(
    LR_pair = paste(
      ligand.complex,
      receptor.complex,
      sep = " → "
    )
  ) %>%
  
  dplyr::select(
    source,
    target,
    LR_pair,
    magnitude_primary,
    magnitude_metastasis,
    specificity_primary,
    specificity_metastasis
  ) %>%
  
  tidyr::pivot_longer(
    cols = c(
      magnitude_primary,
      magnitude_metastasis,
      specificity_primary,
      specificity_metastasis
    ),
    names_to = c(".value", "tumour"),
    names_pattern = "(magnitude|specificity)_(primary|metastasis)"
  ) %>%
  
  dplyr::mutate(
    tumour = factor(
      tumour,
      levels = c(
        "primary",
        "metastasis"
      ),
      labels = c(
        "Primary",
        "Metastasis"
      )
    )
  )


magnitude_limits <- range(
  liana_plot_data$magnitude,
  na.rm = TRUE
)

specificity_limits <- range(
  liana_plot_data$specificity,
  na.rm = TRUE
)

magnitude_limits
specificity_limits

make_and_save_liana_plot <- function(source_cell, target_cell, outdir) {
  
  plot_df <- liana_plot_data %>%
    dplyr::filter(
      source == source_cell,
      target == target_cell
    )
  
  LR_order <- plot_df %>%
    dplyr::group_by(LR_pair) %>%
    dplyr::summarise(
      max_specificity = max(specificity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(max_specificity) %>%
    dplyr::pull(LR_pair)
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      LR_pair = factor(LR_pair, levels = LR_order)
    )
  
  p <- ggplot(plot_df, aes(x = tumour, y = LR_pair)) +
    geom_point(aes(size = specificity, colour = magnitude), na.rm = TRUE) +
    scale_size_continuous(
      name = "Interaction\nspecificity",
      range = c(2, 8),
      limits = c(0, 1),
      breaks = c(0.25, 0.50, 0.75, 1.00)
    ) +
    scale_colour_viridis_c(
      name = "Interaction\nmagnitude",
      option = "viridis",
      direction = -1,
      limits = magnitude_limits
    ) +
    guides(
      size = guide_legend(order = 1),
      colour = guide_colourbar(order = 2)
    ) +
    labs(
      title = paste(source_cell, "→", target_cell),
      x = "Tumour",
      y = "Ligand–receptor interaction"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    )
  
  n_lr <- plot_df %>%
    dplyr::distinct(LR_pair) %>%
    nrow()
  
  safe_source <- gsub("[^A-Za-z0-9]+", "_", source_cell)
  safe_target <- gsub("[^A-Za-z0-9]+", "_", target_cell)
  
  ggsave(
    filename = file.path(outdir, paste0(safe_source, "_to_", safe_target, ".pdf")),
    plot = p,
    width = 9,
    height = max(6, n_lr * 0.25),
    units = "in",
    bg = "white"
  )
  
  return(p)
}

OUTDIR <- "/Users/thirisantracy/Desktop/thesis/images"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

p_CAF_Endo <- make_and_save_liana_plot("CAFs", "Endothelial cells", OUTDIR)
p_CAF_Epi  <- make_and_save_liana_plot("CAFs", "Tumour epithelial cells", OUTDIR)
p_TAM_Endo <- make_and_save_liana_plot("TAMs", "Endothelial cells", OUTDIR)
p_CAF_TAM  <- make_and_save_liana_plot("CAFs", "TAMs", OUTDIR)
p_Endo_CAF <- make_and_save_liana_plot("Endothelial cells", "CAFs", OUTDIR)
p_CAF_TNK <- make_and_save_liana_plot(
  "CAFs",
  "T/NK cells",
  OUTDIR
)

p_Endo_TNK <- make_and_save_liana_plot(
  "Endothelial cells",
  "T/NK cells",
  OUTDIR
)

p_BPlasma_Mono <- make_and_save_liana_plot(
  "B/plasma cells",
  "Monocytes",
  OUTDIR
)

p_Endo_TAM <- make_and_save_liana_plot(
  "Endothelial cells",
  "TAMs",
  OUTDIR
)


#
###########################################################
# Metastasis-only significant interactions
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

liana_met_only %>%
  dplyr::select(
    source,
    target,
    LR_pair,
    rank_metastasis,
    magnitude_metastasis,
    specificity_metastasis
  ) %>%
  dplyr::arrange(
    source,
    target,
    rank_metastasis
  ) %>%
  View()

met_magnitude_limits <- range(
  liana_met_only$magnitude_metastasis,
  na.rm = TRUE
)

met_magnitude_limits


############################################################
# Plot metastasis-only interactions for one cell-pair
############################################################

make_met_only_plot <- function(
    source_cell,
    target_cell
) {
  
  plot_df <- liana_met_only %>%
    dplyr::filter(
      source == source_cell,
      target == target_cell
    )
  
  if (nrow(plot_df) == 0) {
    stop(paste("No metastasis-only interactions found for", source_cell, "→", target_cell))
  }
  
  LR_order <- plot_df %>%
    dplyr::arrange(dplyr::desc(rank_metastasis)) %>%
    dplyr::pull(LR_pair)
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      LR_pair = factor(LR_pair, levels = unique(LR_order)),
      target_label = factor(target_cell, levels = target_cell),
      source_label = factor(source_cell, levels = source_cell)
    )
  
  p <- ggplot(
    plot_df,
    aes(
      x = target_label,
      y = LR_pair
    )
  ) +
    geom_point(
      aes(
        size = specificity_metastasis,
        colour = magnitude_metastasis
      )
    ) +
    facet_grid(. ~ source_label) +
    scale_size_continuous(
      name = "Interaction\nspecificity",
      range = c(2.5, 8),
      limits = c(0, 1),
      breaks = c(0.25, 0.50, 0.75, 1.00)
    ) +
    scale_colour_viridis_c(
      name = "Interaction\nmagnitude",
      option = "viridis",
      direction = -1,
      limits = met_magnitude_limits
    ) +
    guides(
      size = guide_legend(order = 1),
      colour = guide_colourbar(order = 2)
    ) +
    labs(
      x = NULL,
      y = "Ligand–receptor interaction"
    ) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(
        fill = "grey90",
        colour = "grey40",
        linewidth = 0.5
      ),
      strip.text = element_text(
        size = 11,
        face = "plain"
      ),
      
      axis.text.x = element_text(
        size = 10,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.text.y = element_text(size = 8.5),
      
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.35),
      panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      
      plot.margin = margin(t = 8, r = 10, b = 8, l = 8)
    )
  
  return(p)
}
p_met_CAF_Endo <- make_met_only_plot(
  "CAFs",
  "Endothelial cells"
)

p_met_CAF_Endo

p_met_CAF_Epi <- make_met_only_plot(
  "CAFs",
  "Tumour epithelial cells"
)

p_met_CAF_Epi

p_met_TAM_Endo <- make_met_only_plot(
  "TAMs",
  "Endothelial cells"
)

p_met_TAM_Endo

p_met_CAF_TAM <- make_met_only_plot(
  "CAFs",
  "TAMs"
)

p_met_CAF_TAM

p_met_Endo_CAF <- make_met_only_plot(
  "Endothelial cells",
  "CAFs"
)

p_met_Endo_CAF

p_met_Endo_CAF <- make_met_only_plot(
  "Endothelial cells",
  "CAFs"
)

p_met_Endo_CAF

p_met_CAF_TNK <- make_met_only_plot(
  "CAFs",
  "T/NK cells"
)

p_met_Endo_TNK <- make_met_only_plot(
  "Endothelial cells",
  "T/NK cells"
)

p_met_BPlasma_Mono <- make_met_only_plot(
  "B/plasma cells",
  "Monocytes"
)

p_met_Endo_TAM <- make_met_only_plot(
  "Endothelial cells",
  "TAMs"
)

p_met_Mast_CAF <- make_met_only_plot(
  "Mast cells",
  "CAFs"
)

p_met_BPlasma_Mono <- make_met_only_plot(
  "B/plasma cells",
  "Monocytes"
)

p_met_Epi_Mono <- make_met_only_plot(
  "Tumour epithelial cells",
  "Monocytes"
)

p_met_Mast_TAM <- make_met_only_plot(
  "Mast cells",
  "TAMs"
)

p_met_Mast_CAF
p_met_BPlasma_Mono
p_met_Epi_Mono
p_met_Mast_TAM

save_met_only_plot(
  p_met_Mast_CAF,
  "Mast cells",
  "CAFs",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_BPlasma_Mono,
  "B/plasma cells",
  "Monocytes",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_Epi_Mono,
  "Tumour epithelial cells",
  "Monocytes",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_Mast_TAM,
  "Mast cells",
  "TAMs",
  OUTDIR_MET
)

library(cowplot)

shared_legend_new <- cowplot::get_legend(
  p_met_Mast_CAF +
    theme(
      legend.position = "right",
      legend.box = "vertical"
    )
)

plot_list_new <- list(
  p_met_Mast_CAF + theme(legend.position = "none"),
  p_met_BPlasma_Mono + theme(legend.position = "none"),
  p_met_Epi_Mono + theme(legend.position = "none"),
  p_met_Mast_TAM + theme(legend.position = "none")
)

new_grid <- cowplot::plot_grid(
  plotlist = plot_list_new,
  ncol = 2,
  labels = c("a", "b", "c", "d"),
  label_size = 11,
  label_fontface = "bold",
  align = "hv"
)

new_panel <- cowplot::plot_grid(
  new_grid,
  shared_legend_new,
  nrow = 1,
  rel_widths = c(0.88, 0.12)
)

new_panel


save_met_only_plot <- function(
    p,
    source_cell,
    target_cell,
    outdir
) {
  
  ##########################################################
  # Count LR interactions for this source → target pair
  ##########################################################
  
  n_lr <- liana_met_only %>%
    dplyr::filter(
      source == source_cell,
      target == target_cell
    ) %>%
    dplyr::distinct(LR_pair) %>%
    nrow()
  
  
  ##########################################################
  # Safe file names
  ##########################################################
  
  safe_source <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    source_cell
  )
  
  safe_target <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    target_cell
  )
  
  
  ##########################################################
  # Save plot
  ##########################################################
  
  ggsave(
    filename = file.path(
      outdir,
      paste0(
        "Metastasis_only_",
        safe_source,
        "_to_",
        safe_target,
        ".pdf"
      )
    ),
    plot = p,
    width = 6.5,
    height = max(
      4.5,
      n_lr * 0.35
    ),
    units = "in",
    bg = "white"
  )
}
OUTDIR_MET <- "/Users/thirisantracy/Desktop/thesis/images/metastasis_specific"

dir.create(
  OUTDIR_MET,
  showWarnings = FALSE,
  recursive = TRUE
)



save_met_only_plot(
  p_met_CAF_Endo,
  "CAFs",
  "Endothelial cells",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_CAF_Epi,
  "CAFs",
  "Tumour epithelial cells",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_TAM_Endo,
  "TAMs",
  "Endothelial cells",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_CAF_TAM,
  "CAFs",
  "TAMs",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_Endo_CAF,
  "Endothelial cells",
  "CAFs",
  OUTDIR_MET
)
save_met_only_plot(
  p_met_CAF_TNK,
  "CAFs",
  "T/NK cells",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_Endo_TNK,
  "Endothelial cells",
  "T/NK cells",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_BPlasma_Mono,
  "B/plasma cells",
  "Monocytes",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_Endo_TAM,
  "Endothelial cells",
  "TAMs",
  OUTDIR_MET
)

#CAFs and T/NK cells
#endo and T/Nk cells
#B/plasma and monocytes
#endothelial and TAM


selected_pairs <- tibble::tribble(
  ~source,                ~target,
  
  # Existing
  "CAFs",                 "Endothelial cells",
  "CAFs",                 "Tumour epithelial cells",
  "TAMs",                 "Endothelial cells",
  "CAFs",                 "TAMs",
  "Endothelial cells",    "CAFs",
  "CAFs",                 "T/NK cells",
  "Endothelial cells",    "T/NK cells",
  "B/plasma cells",       "Monocytes",
  "Endothelial cells",    "TAMs",
  
  # Newly noticed metastasis-only pairs
  "Mast cells",           "CAFs",
  "Tumour epithelial cells", "Monocytes",
  "Mast cells",           "TAMs"
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

make_liana_graph <- function(edges, cell_order) {
  
  vertices <- data.frame(
    name = cell_order
  )
  
  igraph::graph_from_data_frame(
    d = edges,
    directed = TRUE,
    vertices = vertices
  )
}


primary_graph <- make_liana_graph(
  primary_edges,
  cell_order
)

met_graph <- make_liana_graph(
  met_edges,
  cell_order
)

n_cells <- length(cell_order)

angles <- seq(
  pi / 2,
  pi / 2 - 2 * pi,
  length.out = n_cells + 1
)[1:n_cells]

fixed_layout <- cbind(
  cos(angles),
  sin(angles)
)

rownames(fixed_layout) <- cell_order

max_interactions <- max(
  c(
    primary_edges$n_interactions,
    met_edges$n_interactions
  ),
  na.rm = TRUE
)

max_interactions

scale_edge_width <- function(x) {
  
  0.4 +
    6 *
    sqrt(x / max_interactions)
}

if (exists("celltype_colours")) {
  
  node_colours <- rep(
    "grey75",
    length(cell_order)
  )
  
  names(node_colours) <- cell_order
  
  matching_colours <- intersect(
    cell_order,
    names(celltype_colours)
  )
  
  node_colours[matching_colours] <-
    celltype_colours[matching_colours]
  
} else {
  
  node_colours <- scales::hue_pal()(
    length(cell_order)
  )
  
  names(node_colours) <- cell_order
}

node_colours

plot_liana_network <- function(
    graph,
    title_text
) {
  
  ##########################################################
  # Fixed layout
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
  # More restrained edge-width scaling
  ##########################################################
  
  edge_widths <- scales::rescale(
    sqrt(edge_counts),
    to = c(0.4, 4)
  )
  
  
  ##########################################################
  # Colour edge according to SENDER cell
  ##########################################################
  
  edge_colours <- node_colours[
    edge_sources
  ]
  
  edge_colours <- grDevices::adjustcolor(
    edge_colours,
    alpha.f = 0.45
  )
  
  
  ##########################################################
  # Plot
  ##########################################################
  
  plot(
    graph,
    
    layout = plot_layout,
    
    
    ########################################################
    # Nodes
    ########################################################
    
    vertex.color =
      node_colours[
        igraph::V(graph)$name
      ],
    
    vertex.frame.color = "grey30",
    
    vertex.size = 24,
    
    
    ########################################################
    # Labels
    ########################################################
    
    vertex.label =
      igraph::V(graph)$name,
    
    vertex.label.cex = 0.8,
    
    vertex.label.color = "black",
    
    vertex.label.dist = 1.35,
    
    
    ########################################################
    # Edges
    ########################################################
    
    edge.width = edge_widths,
    
    edge.color = edge_colours,
    
    edge.arrow.size = 0.12,
    
    edge.arrow.width = 0.7,
    
    edge.curved = 0.22,
    
    
    ########################################################
    # Self loops
    ########################################################
    
    loop.angle = pi / 3,
    
    
    ########################################################
    # Plot appearance
    ########################################################
    
    asp = 1,
    
    margin = 0.4,
    
    main = title_text
  )
}

plot_liana_network(
  primary_graph,
  "Primary PDAC"
)


plot_liana_network(
  met_graph,
  "Metastatic PDAC"
)

par(
  mfrow = c(1, 2),
  xpd = TRUE
)

plot_liana_network(
  primary_graph,
  "Primary PDAC"
)

plot_liana_network(
  met_graph,
  "Metastatic PDAC"
)

par(
  mfrow = c(1, 1)
)

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
# PRIMARY — full network
############################################################

quartz(
  type = "pdf",
  file = file.path(
    OUTDIR,
    "LIANA_primary_network_full.pdf"
  ),
  width = 6,
  height = 6
)

plot_liana_network_clean(
  primary_graph_full,
  "Primary PDAC"
)

dev.off()


############################################################
# METASTATIC — full network
############################################################

quartz(
  type = "pdf",
  file = file.path(
    OUTDIR,
    "LIANA_metastatic_network_full.pdf"
  ),
  width = 6,
  height = 6
)

plot_liana_network_clean(
  met_graph_full,
  "Metastatic PDAC"
)

dev.off()

##### recurrence
############################################################
# Ligand recurrence / promiscuity analysis
############################################################

library(dplyr)
library(tidyr)
library(stringr)

liana_all_sig <- dplyr::bind_rows(
  
  primary_sig %>%
    dplyr::mutate(
      condition = "Primary"
    ),
  
  met_sig %>%
    dplyr::mutate(
      condition = "Metastasis"
    )
)

liana_all_sig <- liana_all_sig %>%
  
  dplyr::distinct(
    condition,
    source,
    target,
    ligand.complex,
    receptor.complex,
    .keep_all = TRUE
  )

nrow(liana_all_sig)

dplyr::n_distinct(
  liana_all_sig$ligand.complex
)

dplyr::n_distinct(
  liana_all_sig$receptor.complex
)

############################################################
# Ligand recurrence within each condition
############################################################

ligand_recurrence_by_condition <- liana_all_sig %>%
  
  dplyr::mutate(
    sender_receiver = paste(
      source,
      target,
      sep = " → "
    ),
    
    ligand_receptor = paste(
      ligand.complex,
      receptor.complex,
      sep = " → "
    )
  ) %>%
  
  dplyr::group_by(
    condition,
    ligand.complex
  ) %>%
  
  dplyr::summarise(
    
    n_sender_celltypes =
      dplyr::n_distinct(source),
    
    n_receiver_celltypes =
      dplyr::n_distinct(target),
    
    n_sender_receiver_contexts =
      dplyr::n_distinct(sender_receiver),
    
    n_receptors =
      dplyr::n_distinct(receptor.complex),
    
    n_interaction_contexts =
      dplyr::n_distinct(
        paste(
          source,
          target,
          receptor.complex,
          sep = " | "
        )
      ),
    
    sender_celltypes =
      paste(
        sort(unique(source)),
        collapse = "; "
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
  )


ligand_recurrence_by_condition %>%
  
  dplyr::arrange(
    dplyr::desc(n_receiver_celltypes),
    dplyr::desc(n_sender_receiver_contexts),
    dplyr::desc(n_receptors)
  ) %>%
  
  print(
    n = Inf,
    width = Inf
  )

View(
  ligand_recurrence_by_condition %>%
    dplyr::arrange(
      dplyr::desc(n_receiver_celltypes),
      dplyr::desc(n_sender_receiver_contexts),
      dplyr::desc(n_receptors)
    )
)

primary_ligand_recurrence <- ligand_recurrence_by_condition %>%
  
  dplyr::filter(
    condition == "Primary"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_receiver_celltypes),
    dplyr::desc(n_sender_receiver_contexts),
    dplyr::desc(n_receptors)
  )

View(
  primary_ligand_recurrence
)

met_ligand_recurrence <- ligand_recurrence_by_condition %>%
  
  dplyr::filter(
    condition == "Metastasis"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_receiver_celltypes),
    dplyr::desc(n_sender_receiver_contexts),
    dplyr::desc(n_receptors)
  )

View(
  met_ligand_recurrence
)


############################################################
# Primary vs metastatic ligand recurrence comparison
############################################################

ligand_recurrence_compare <- ligand_recurrence_by_condition %>%
  
  dplyr::select(
    condition,
    ligand.complex,
    n_sender_celltypes,
    n_receiver_celltypes,
    n_sender_receiver_contexts,
    n_receptors,
    n_interaction_contexts
  ) %>%
  
  tidyr::pivot_wider(
    
    names_from = condition,
    
    values_from = c(
      n_sender_celltypes,
      n_receiver_celltypes,
      n_sender_receiver_contexts,
      n_receptors,
      n_interaction_contexts
    ),
    
    values_fill = 0
  )

ligand_recurrence_compare <- ligand_recurrence_compare %>%
  
  dplyr::mutate(
    
    delta_receiver_celltypes =
      n_receiver_celltypes_Metastasis -
      n_receiver_celltypes_Primary,
    
    delta_sender_receiver_contexts =
      n_sender_receiver_contexts_Metastasis -
      n_sender_receiver_contexts_Primary,
    
    delta_receptors =
      n_receptors_Metastasis -
      n_receptors_Primary,
    
    delta_interaction_contexts =
      n_interaction_contexts_Metastasis -
      n_interaction_contexts_Primary
  )

View(
  ligand_recurrence_compare
)

############################################################
# Ligand condition classification
############################################################

ligand_condition_status <- liana_all_sig %>%
  
  dplyr::distinct(
    ligand.complex,
    condition
  ) %>%
  
  dplyr::mutate(
    present = TRUE
  ) %>%
  
  tidyr::pivot_wider(
    names_from = condition,
    values_from = present,
    values_fill = FALSE
  ) %>%
  
  dplyr::mutate(
    
    condition_status =
      dplyr::case_when(
        
        Primary & Metastasis ~
          "Shared",
        
        Primary & !Metastasis ~
          "Primary only",
        
        !Primary & Metastasis ~
          "Metastasis only"
      )
  )

ligand_condition_status %>%
  
  dplyr::count(
    condition_status
  )

ligand_recurrence_compare <- ligand_recurrence_compare %>%
  
  dplyr::left_join(
    ligand_condition_status %>%
      dplyr::select(
        ligand.complex,
        condition_status
      ),
    by = "ligand.complex"
  )


ligand_recurrence_compare %>%
  
  dplyr::arrange(
    dplyr::desc(
      pmax(
        n_receiver_celltypes_Primary,
        n_receiver_celltypes_Metastasis
      )
    ),
    dplyr::desc(
      pmax(
        n_sender_receiver_contexts_Primary,
        n_sender_receiver_contexts_Metastasis
      )
    )
  ) %>%
  
  print(
    n = Inf,
    width = Inf
  )

############################################################
# Ligands with broader communication in metastasis
############################################################

met_more_recurrent <- ligand_recurrence_compare %>%
  
  dplyr::filter(
    delta_receiver_celltypes > 0 |
      delta_sender_receiver_contexts > 0 |
      delta_receptors > 0
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(delta_receiver_celltypes),
    dplyr::desc(delta_sender_receiver_contexts),
    dplyr::desc(delta_receptors)
  )

View(
  met_more_recurrent
)

primary_more_recurrent <- ligand_recurrence_compare %>%
  
  dplyr::filter(
    delta_receiver_celltypes < 0 |
      delta_sender_receiver_contexts < 0 |
      delta_receptors < 0
  ) %>%
  
  dplyr::arrange(
    delta_receiver_celltypes,
    delta_sender_receiver_contexts,
    delta_receptors
  )

View(
  primary_more_recurrent
)


met_only_ligands <- ligand_condition_status %>%
  
  dplyr::filter(
    condition_status == "Metastasis only"
  ) %>%
  
  dplyr::pull(
    ligand.complex
  )


met_only_interactions <- liana_all_sig %>%
  
  dplyr::filter(
    condition == "Metastasis",
    ligand.complex %in% met_only_ligands
  )

length(met_only_ligands)

dplyr::n_distinct(met_only_interactions$ligand.complex)

nrow(met_only_interactions)

met_only_ligand_summary <- met_only_interactions %>%
  
  dplyr::mutate(
    sender_receiver = paste(source, target, sep = " → ")
  ) %>%
  
  dplyr::group_by(
    ligand.complex
  ) %>%
  
  dplyr::summarise(
    
    n_sender_celltypes =
      dplyr::n_distinct(source),
    
    n_receiver_celltypes =
      dplyr::n_distinct(target),
    
    n_sender_receiver_contexts =
      dplyr::n_distinct(sender_receiver),
    
    n_receptors =
      dplyr::n_distinct(receptor.complex),
    
    sender_celltypes =
      paste(sort(unique(source)), collapse = "; "),
    
    receiver_celltypes =
      paste(sort(unique(target)), collapse = "; "),
    
    receptors =
      paste(sort(unique(receptor.complex)), collapse = "; "),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_receiver_celltypes),
    dplyr::desc(n_sender_receiver_contexts),
    dplyr::desc(n_receptors)
  )


View(met_only_ligand_summary)

met_only_ligand_summary %>%
  print(n = Inf, width = Inf)


met_only_interactions <- liana_all_sig %>%
  
  dplyr::filter(
    condition == "Metastasis",
    ligand.complex %in% met_only_ligands,
    source != target
  )

length(met_only_ligands)

dplyr::n_distinct(
  met_only_interactions$ligand.complex
)


met_only_ligand_summary <- met_only_interactions %>%
  
  dplyr::mutate(
    sender_receiver = paste(
      source,
      target,
      sep = " → "
    )
  ) %>%
  
  dplyr::group_by(
    ligand.complex
  ) %>%
  
  dplyr::summarise(
    
    n_sender_celltypes =
      dplyr::n_distinct(source),
    
    n_receiver_celltypes =
      dplyr::n_distinct(target),
    
    n_sender_receiver_contexts =
      dplyr::n_distinct(sender_receiver),
    
    n_receptors =
      dplyr::n_distinct(receptor.complex),
    
    sender_celltypes =
      paste(
        sort(unique(source)),
        collapse = "; "
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
    dplyr::desc(n_sender_receiver_contexts),
    dplyr::desc(n_receptors)
  )
met_only_ligand_summary %>%
  print(
    n = Inf,
    width = Inf
  )

met_only_interactions %>%
  dplyr::distinct(
    ligand.complex,
    source
  ) %>%
  dplyr::count(
    source,
    sort = TRUE
  )

met_only_interactions %>%
  dplyr::distinct(
    ligand.complex,
    target
  ) %>%
  dplyr::count(
    target,
    sort = TRUE
  )

met_only_interactions %>%
  
  dplyr::distinct(
    ligand.complex,
    source,
    target
  ) %>%
  
  dplyr::arrange(
    ligand.complex,
    source,
    target
  ) %>%
  
  print(
    n = Inf
  )

View(met_only_interactions)



############################################################
# Metastasis-only ligands
# remove self-signalling for this analysis
############################################################

met_only_ligands <- ligand_condition_status %>%
  dplyr::filter(
    condition_status == "Metastasis only"
  ) %>%
  dplyr::pull(
    ligand.complex
  )

met_only_interactions <- liana_all_sig %>%
  dplyr::filter(
    condition == "Metastasis",
    ligand.complex %in% met_only_ligands,
    source != target
  )


############################################################
# Summary of metastasis-only ligands
############################################################

met_only_ligand_summary <- met_only_interactions %>%
  
  dplyr::mutate(
    sender_receiver = paste(
      source,
      target,
      sep = " → "
    )
  ) %>%
  
  dplyr::group_by(
    ligand.complex
  ) %>%
  
  dplyr::summarise(
    n_sender_celltypes =
      dplyr::n_distinct(source),
    
    n_receiver_celltypes =
      dplyr::n_distinct(target),
    
    n_sender_receiver_contexts =
      dplyr::n_distinct(sender_receiver),
    
    n_receptors =
      dplyr::n_distinct(receptor.complex),
    
    sender_celltypes =
      paste(
        sort(unique(source)),
        collapse = "; "
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
    dplyr::desc(n_sender_receiver_contexts),
    dplyr::desc(n_receptors)
  )
met_only_ligand_summary %>%
  print(n = Inf, width = Inf)

View(met_only_ligand_summary)

############################################################
# Ligands targeting >=2 receiver cell types
############################################################

met_only_recurrent_ligands <- met_only_ligand_summary %>%
  dplyr::filter(
    n_receiver_celltypes >= 2
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_receiver_celltypes),
    dplyr::desc(n_sender_receiver_contexts),
    ligand.complex
  )

met_only_recurrent_ligands %>%
  print(n = Inf, width = Inf)


############################################################
# UpSet input:
# rows = ligands
# columns = receiver cell types
############################################################

met_only_upset_df <- met_only_interactions %>%
  
  dplyr::filter(
    ligand.complex %in%
      met_only_recurrent_ligands$ligand.complex
  ) %>%
  
  dplyr::distinct(
    ligand.complex,
    target
  ) %>%
  
  dplyr::mutate(
    present = 1L
  ) %>%
  
  tidyr::pivot_wider(
    names_from = target,
    values_from = present,
    values_fill = 0L
  )

met_only_upset_df

############################################################
# Convert to UpSetR input
############################################################

met_only_upset_input <- as.data.frame(
  met_only_upset_df
)

rownames(met_only_upset_input) <-
  met_only_upset_input$ligand.complex

met_only_upset_input$ligand.complex <- NULL

receiver_order <- c(
  "Endothelial cells",
  "Tumour epithelial cells",
  "TAMs",
  "Monocytes",
  "T/NK cells",
  "CAFs",
  "Mast cells",
  "B/plasma cells"
)

receiver_order <- intersect(
  receiver_order,
  colnames(met_only_upset_input)
)

receiver_order

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
library(UpSetR)
############################################################
# UpSet plot
############################################################

UpSetR::upset(
  met_only_upset_input,
  
  sets = receiver_order,
  
  keep.order = TRUE,
  
  order.by = "freq",
  
  nsets = length(receiver_order),
  
  nintersects = NA,
  
  mb.ratio = c(0.6, 0.4),
  
  point.size = 3,
  
  line.size = 1,
  
  text.scale = c(
    1.4,  # main bar title
    1.2,  # main bar ticks
    1.2,  # set size title
    1.0,  # set size ticks
    1.2,  # set names
    1.2   # numbers above bars
  ),
  
  mainbar.y.label =
    "Number of metastasis-only ligands\nshared across receiver cell types",
  
  sets.x.label =
    "Number of metastasis-only ligands"
)

############################################################
# Save UpSet plot
############################################################

pdf(
  file.path(
    OUTDIR,
    "LIANA_metastasis_only_ligands_receiver_upset.pdf"
  ),
  width = 8,
  height = 5.5
)

UpSetR::upset(
  met_only_upset_input,
  
  sets = receiver_order,
  
  keep.order = TRUE,
  
  order.by = "freq",
  
  nsets = length(receiver_order),
  
  nintersects = NA,
  
  mb.ratio = c(0.6, 0.4),
  
  point.size = 3,
  
  line.size = 1,
  
  text.scale = c(1.4, 1.2, 1.2, 1.0, 1.2, 1.2),
  
  mainbar.y.label =
    "Number of metastasis-only ligands\nshared across receiver cell types",
  
  sets.x.label =
    "Number of metastasis-only ligands"
)

dev.off()

write.csv(
  met_only_ligand_summary,
  file.path(
    OUTDIR,
    "LIANA_metastasis_only_ligand_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  met_only_recurrent_ligands,
  file.path(
    OUTDIR,
    "LIANA_metastasis_only_recurrent_ligands.csv"
  ),
  row.names = FALSE
)

################# doing last figure plot
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
# 4. Global scales so all 6 plots are comparable
############################################################

specificity_limits <- range(
  caf_all_plot_df$natmi.edge_specificity,
  na.rm = TRUE
)

lrscore_limits <- range(
  caf_all_plot_df$sca.LRscore,
  na.rm = TRUE
)

specificity_breaks <- pretty(
  specificity_limits,
  n = 4
)

lrscore_breaks <- pretty(
  lrscore_limits,
  n = 4
)

############################################################
# 5. Small helper for safe file names
############################################################

clean_filename <- function(x) {
  
  x %>%
    stringr::str_replace_all("/", "_") %>%
    stringr::str_replace_all(" ", "_")
}

############################################################
# 6. Function to make one CAF -> target plot
############################################################

make_caf_target_dotplot <- function(
    target_name,
    show_legend = TRUE
) {
  
  plot_df <- caf_all_plot_df %>%
    
    dplyr::filter(
      target == target_name
    ) %>%
    
    dplyr::arrange(
      natmi.edge_specificity,
      sca.LRscore
    ) %>%
    
    dplyr::mutate(
      
      ######################################################
      # Put strongest interactions higher up
      ######################################################
      
      interaction = factor(
        interaction,
        levels = unique(interaction)
      ),
      
      ######################################################
      # Keep source strip at the top
      ######################################################
      
      source = factor(
        source,
        levels = "CAFs"
      ),
      
      ######################################################
      # X-axis target label
      ######################################################
      
      target = factor(
        target,
        levels = target_name
      )
    )
  
  
  p <- ggplot(
    plot_df,
    aes(
      x = target,
      y = interaction
    )
  ) +
    
    geom_point(
      aes(
        size = natmi.edge_specificity,
        colour = sca.LRscore
      )
    ) +
    
    ########################################################
  # Top strip with source name, like LIANA-style figure
  ########################################################
  
  facet_grid(
    . ~ source,
    scales = "free_x",
    space = "free_x"
  ) +
    scale_x_discrete(
      expand = expansion(add = 1)
    ) +
    
    ########################################################
  # SAME scales across all 6 plots
  ########################################################
  
  scale_size_continuous(
    name = "Interaction\nspecificity",
    limits = specificity_limits,
    breaks = specificity_breaks,
    range = c(2.2, 8)
  ) +
    
    scale_colour_viridis_c(
      name = "Interaction\nmagnitude",
      limits = lrscore_limits,
      breaks = lrscore_breaks,
      option = "D",
      direction = -1
    ) +
    
    labs(
      x = "Target",
      y = "Ligand–receptor interactions"
    ) +
    ########################################################
  # Theme to resemble LIANA-style dot plots
  ########################################################
  
  theme_bw(base_size = 11) +
    
    theme(
      
      ######################################################
      # Facet strip
      ######################################################
      
      strip.background = element_rect(
        fill = "grey90",
        colour = "grey40",
        linewidth = 0.5
      ),
      
      strip.text = element_text(
        size = 11,
        face = "plain"
      ),
      
      ######################################################
      # Axes
      ######################################################
      
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      
      axis.text.x = element_text(
        size = 10,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      
      axis.text.y = element_text(size = 8),
      
      ######################################################
      # Gridlines
      ######################################################
      
      panel.grid.major.x = element_line(
        colour = "grey85"
      ),
      
      panel.grid.major.y = element_line(
        colour = "grey90"
      ),
      
      panel.grid.minor = element_blank(),
      
      ######################################################
      # Legends
      ######################################################
      
      legend.position = if (show_legend) "right" else "none",
      
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      
      ######################################################
      # Margins
      ######################################################
      
      plot.margin = margin(
        t = 8,
        r = 10,
        b = 8,
        l = 8
      )
    )
  
  return(p)
}

############################################################
# 7. Make the 6 plots
############################################################

caf_dotplots <- lapply(
  caf_targets,
  make_caf_target_dotplot
)

names(caf_dotplots) <- caf_targets

############################################################
# 8. View one if you want
############################################################

caf_dotplots[["Tumour epithelial cells"]]
caf_dotplots[["Endothelial cells"]]
caf_dotplots[["Mast cells"]]
caf_dotplots[["Monocytes"]]
caf_dotplots[["T/NK cells"]]
caf_dotplots[["TAMs"]]

############################################################
# 9. Save all 6 with identical dimensions
############################################################

############################################################
# Save all six CAF dot plots as PDF
############################################################

OUTDIR <- "/Users/thirisantracy/Desktop/thesis/images"

dir.create(
  OUTDIR,
  showWarnings = FALSE,
  recursive = TRUE
)

plot_width  <- 4.2
plot_height <- 5.6

for (tg in caf_targets) {
  
  outfile <- file.path(
    OUTDIR,
    paste0(
      "LIANA_CAF_to_",
      clean_filename(tg),
      "_dotplot.pdf"
    )
  )
  
  ggsave(
    filename = outfile,
    plot = caf_dotplots[[tg]],
    width = plot_width,
    height = plot_height,
    units = "in",
    device = "pdf",
    useDingbats = FALSE,
    bg = "white"
  )
  
  message(
    tg,
    " -> ",
    outfile,
    " | saved = ",
    file.exists(outfile)
  )
}


library(ggplot2)
library(cowplot)

############################################################
# 1. Order of the seven CAF target panels
############################################################

panel_order <- c(
  "Tumour epithelial cells",
  "Endothelial cells",
  "Mast cells",
  "Monocytes",
  "T/NK cells",
  "TAMs",
  "B/plasma cells"
)

############################################################
# 2. Rebuild CAF dotplots
############################################################

caf_dotplots <- lapply(
  panel_order,
  make_caf_target_dotplot
)

names(caf_dotplots) <- panel_order

############################################################
# 3. Shared legend only
############################################################

shared_legend <- cowplot::get_legend(
  caf_dotplots[["Tumour epithelial cells"]] +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
)

############################################################
# 4. Clean panel versions
# - no repeated legends
# - no repeated strip title
# - no axis titles
# - keep target labels vertical
############################################################

caf_panel_plots <- lapply(
  panel_order,
  function(tg) {
    caf_dotplots[[tg]] +
      scale_y_discrete(
        expand = expansion(add = c(0.6, 0.6))
      ) +
      theme(
        legend.position = "none",
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(
          size = 10,
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          margin = margin(t = 4)
        ),
        axis.text.y = element_text(size = 8.5),
        panel.grid.major.y = element_line(
          colour = "grey85",
          linewidth = 0.35
        ),
        panel.grid.major.x = element_line(
          colour = "grey88",
          linewidth = 0.35
        ),
        panel.grid.minor = element_blank(),
        plot.margin = margin(
          t = 4,
          r = 4,
          b = 4,
          l = 4
        )
      )
  }
)

names(caf_panel_plots) <- panel_order

############################################################
# TOP ROW — 4 equal-sized panels
############################################################

top_row <- cowplot::plot_grid(
  
  caf_panel_plots[["Tumour epithelial cells"]],
  caf_panel_plots[["Endothelial cells"]],
  caf_panel_plots[["Mast cells"]],
  caf_panel_plots[["Monocytes"]],
  
  nrow = 1,

  label_size = 11,
  label_fontface = "bold",
  
  hjust = -0.25,
  vjust = 1.1,
  
  align = "hv",
  axis = "tblr",
  
  rel_widths = c(
    1, 1, 1, 1
  )
)


############################################################
# BOTTOM ROW — same panel width as top row
#
# Half-width blank spacer on each side:
#
# 0.5 + 1 + 1 + 1 + 0.5 = 4
#
# So each actual plot still occupies 1/4 of total width
############################################################

bottom_row <- cowplot::plot_grid(
  
  NULL,
  
  caf_panel_plots[["T/NK cells"]],
  caf_panel_plots[["TAMs"]],
  caf_panel_plots[["B/plasma cells"]],
  
  NULL,
  
  nrow = 1,
  
  rel_widths = c(
    0.5,
    1,
    1,
    1,
    0.5
  ),
  
  label_size = 11,
  label_fontface = "bold",
  
  hjust = -0.25,
  vjust = 1.1,
  
  align = "hv",
  axis = "tblr"
)


############################################################
# STACK THE TWO ROWS
############################################################

caf_seven_grid <- cowplot::plot_grid(
  
  top_row,
  bottom_row,
  
  ncol = 1,
  
  rel_heights = c(
    1,
    1
  )
)


############################################################
# ADD ONE SHARED LEGEND
############################################################

caf_lr_panel <- cowplot::plot_grid(
  
  caf_seven_grid,
  shared_legend,
  
  nrow = 1,
  
  rel_widths = c(
    0.90,
    0.10
  )
)


############################################################
# PREVIEW
############################################################

caf_lr_panel


ggsave(
  filename = file.path(
    OUTDIR,
    "LIANA_CAF_met_specific_LR_panel_with_Bplasma.pdf"
  ),
  plot = caf_lr_panel,
  width = 13,
  height = 10,
  units = "in",
  device = "pdf",
  useDingbats = FALSE,
  bg = "white"
)



caf_receptor_recurrence <- caf_met_specific %>%
  
  dplyr::group_by(
    receptor.complex
  ) %>%
  
  dplyr::summarise(
    
    n_receiver_celltypes =
      dplyr::n_distinct(target),
    
    n_ligands =
      dplyr::n_distinct(ligand.complex),
    
    n_LR_contexts =
      dplyr::n_distinct(
        paste(
          target,
          ligand.complex,
          sep = " | "
        )
      ),
    
    receiver_celltypes =
      paste(
        sort(unique(target)),
        collapse = "; "
      ),
    
    ligands =
      paste(
        sort(unique(ligand.complex)),
        collapse = "; "
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_receiver_celltypes),
    dplyr::desc(n_ligands),
    dplyr::desc(n_LR_contexts)
  )

caf_receptor_recurrence %>%
  print(
    n = Inf,
    width = Inf
  )

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

############################################################
# Check what is present
############################################################

liana_met_only %>%
  dplyr::count(source, target, name = "n_LR_pairs") %>%
  dplyr::arrange(source, target) %>%
  print(n = Inf)

############################################################
# Shared colour scale across met-only plots
############################################################

met_magnitude_limits <- range(
  liana_met_only$magnitude_metastasis,
  na.rm = TRUE
)

############################################################
# Plot metastasis-only interactions in the SAME style
# as your CAF-specific met-only plots
############################################################

make_met_only_plot <- function(
    source_cell,
    target_cell
) {
  
  plot_df <- liana_met_only %>%
    dplyr::filter(
      source == source_cell,
      target == target_cell
    )
  
  if (nrow(plot_df) == 0) {
    stop(
      paste(
        "No metastasis-only interactions found for",
        source_cell, "→", target_cell
      )
    )
  }
  
  LR_order <- plot_df %>%
    dplyr::arrange(dplyr::desc(rank_metastasis)) %>%
    dplyr::pull(LR_pair)
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      LR_pair = factor(
        LR_pair,
        levels = unique(LR_order)
      ),
      target_label = factor(
        target_cell,
        levels = target_cell
      ),
      source_label = factor(
        source_cell,
        levels = source_cell
      )
    )
  
  p <- ggplot(
    plot_df,
    aes(
      x = target_label,
      y = LR_pair
    )
  ) +
    geom_point(
      aes(
        size = specificity_metastasis,
        colour = magnitude_metastasis
      )
    ) +
    facet_grid(. ~ source_label) +
    scale_size_continuous(
      name = "Interaction\nspecificity",
      range = c(2.5, 8),
      limits = c(0, 1),
      breaks = c(0.25, 0.50, 0.75, 1.00)
    ) +
    scale_colour_viridis_c(
      name = "Interaction\nmagnitude",
      option = "viridis",
      direction = -1,
      limits = met_magnitude_limits
    ) +
    guides(
      size = guide_legend(order = 1),
      colour = guide_colourbar(order = 2)
    ) +
    labs(
      x = NULL,
      y = "Ligand–receptor interaction"
    ) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(
        fill = "grey90",
        colour = "grey40",
        linewidth = 0.5
      ),
      strip.text = element_text(
        size = 11,
        face = "plain"
      ),
      
      axis.text.x = element_text(
        size = 10,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.text.y = element_text(size = 8.5),
      
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      
      panel.grid.major.x = element_line(
        colour = "grey88",
        linewidth = 0.35
      ),
      panel.grid.major.y = element_line(
        colour = "grey85",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      
      plot.margin = margin(
        t = 8, r = 10, b = 8, l = 8
      )
    )
  
  return(p)
}

############################################################
# Save function
############################################################

save_met_only_plot <- function(
    p,
    source_cell,
    target_cell,
    outdir
) {
  
  n_lr <- liana_met_only %>%
    dplyr::filter(
      source == source_cell,
      target == target_cell
    ) %>%
    dplyr::distinct(LR_pair) %>%
    nrow()
  
  safe_source <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    source_cell
  )
  
  safe_target <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    target_cell
  )
  
  ggsave(
    filename = file.path(
      outdir,
      paste0(
        "Metastasis_only_",
        safe_source,
        "_to_",
        safe_target,
        ".pdf"
      )
    ),
    plot = p,
    width = 6.5,
    height = max(4.5, n_lr * 0.35),
    units = "in",
    device = "pdf",
    bg = "white"
  )
}

############################################################
# Output folder
############################################################

OUTDIR_MET <- "/Users/thirisantracy/Desktop/thesis/images/metastasis_specific"

dir.create(
  OUTDIR_MET,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# Make the FOUR new plots
############################################################

p_met_Mast_CAF <- make_met_only_plot(
  "Mast cells",
  "CAFs"
)

p_met_BPlasma_Mono <- make_met_only_plot(
  "B/plasma cells",
  "Monocytes"
)

p_met_Epi_Mono <- make_met_only_plot(
  "Tumour epithelial cells",
  "Monocytes"
)

p_met_Mast_TAM <- make_met_only_plot(
  "Mast cells",
  "TAMs"
)

############################################################
# Show them
############################################################

p_met_Mast_CAF
p_met_BPlasma_Mono
p_met_Epi_Mono
p_met_Mast_TAM

############################################################
# Save them individually
############################################################

save_met_only_plot(
  p_met_Mast_CAF,
  "Mast cells",
  "CAFs",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_BPlasma_Mono,
  "B/plasma cells",
  "Monocytes",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_Epi_Mono,
  "Tumour epithelial cells",
  "Monocytes",
  OUTDIR_MET
)

save_met_only_plot(
  p_met_Mast_TAM,
  "Mast cells",
  "TAMs",
  OUTDIR_MET
)


library(cowplot)

shared_legend_new <- cowplot::get_legend(
  p_met_Mast_CAF +
    theme(
      legend.position = "right",
      legend.box = "vertical"
    )
)

plot_list_new <- list(
  p_met_Mast_CAF + theme(legend.position = "none"),
  p_met_BPlasma_Mono + theme(legend.position = "none"),
  p_met_Epi_Mono + theme(legend.position = "none"),
  p_met_Mast_TAM + theme(legend.position = "none")
)

new_grid <- cowplot::plot_grid(
  plotlist = plot_list_new,
  ncol = 2,
  labels = c("a", "b", "c", "d"),
  label_size = 11,
  label_fontface = "bold",
  align = "hv"
)

new_panel <- cowplot::plot_grid(
  new_grid,
  shared_legend_new,
  nrow = 1,
  rel_widths = c(0.88, 0.12)
)

new_panel

ggsave(
  filename = file.path(
    OUTDIR_MET,
    "Metastasis_only_new_pairs_panel.pdf"
  ),
  plot = new_panel,
  width = 10,
  height = 10,
  units = "in",
  device = "pdf",
  bg = "white"
)


library(ggplot2)
library(dplyr)
library(cowplot)
library(stringr)

############################################################
# Helper to wrap long strip labels
############################################################

wrap_strip_label <- function(x) {
  dplyr::case_when(
    x == "Tumour epithelial cells" ~ "Tumour epithelial\ncells",
    x == "B/plasma cells" ~ "B/plasma\ncells",
    TRUE ~ x
  )
}

############################################################
# Shared magnitude limits
############################################################

met_magnitude_limits <- range(
  liana_met_only$magnitude_metastasis,
  na.rm = TRUE
)

############################################################
# Rebuild metastasis-only plot function
############################################################

make_met_only_plot <- function(
    source_cell,
    target_cell
) {
  
  plot_df <- liana_met_only %>%
    dplyr::filter(
      source == source_cell,
      target == target_cell
    )
  
  if (nrow(plot_df) == 0) {
    stop(
      paste(
        "No metastasis-only interactions found for",
        source_cell, "→", target_cell
      )
    )
  }
  
  LR_order <- plot_df %>%
    dplyr::arrange(dplyr::desc(rank_metastasis)) %>%
    dplyr::pull(LR_pair)
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      LR_pair = factor(
        LR_pair,
        levels = unique(LR_order)
      ),
      target_label = factor(
        target_cell,
        levels = target_cell
      ),
      source_label = factor(
        wrap_strip_label(source_cell),
        levels = wrap_strip_label(source_cell)
      )
    )
  
  ggplot(
    plot_df,
    aes(
      x = target_label,
      y = LR_pair
    )
  ) +
    geom_point(
      aes(
        size = specificity_metastasis,
        colour = magnitude_metastasis
      )
    ) +
    facet_grid(. ~ source_label) +
    scale_size_continuous(
      name = "Interaction\nspecificity",
      range = c(2.5, 8),
      limits = c(0, 1),
      breaks = c(0.25, 0.50, 0.75)   # removed 1.00
    ) +
    scale_colour_viridis_c(
      name = "Interaction\nmagnitude",
      option = "viridis",
      direction = -1,
      limits = met_magnitude_limits
    ) +
    guides(
      size = guide_legend(order = 1),
      colour = guide_colourbar(order = 2)
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(
        fill = "grey90",
        colour = "grey40",
        linewidth = 0.5
      ),
      strip.text = element_text(
        size = 9,
        face = "plain",
        lineheight = 0.9
      ),
      
      axis.text.x = element_text(
        size = 9,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.text.y = element_text(size = 8.5),
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      panel.grid.major.x = element_line(
        colour = "grey88",
        linewidth = 0.35
      ),
      panel.grid.major.y = element_line(
        colour = "grey85",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      
      plot.margin = margin(
        t = 6, r = 6, b = 6, l = 6
      )
    )
}

p_met_Mast_CAF <- make_met_only_plot(
  "Mast cells",
  "CAFs"
)

p_met_BPlasma_Mono <- make_met_only_plot(
  "B/plasma cells",
  "Monocytes"
)

p_met_Epi_Mono <- make_met_only_plot(
  "Tumour epithelial cells",
  "Monocytes"
)

p_met_Mast_TAM <- make_met_only_plot(
  "Mast cells",
  "TAMs"
)

############################################################
# CAF plots to include
# Use the 6 you want
############################################################

caf_targets_for_panel <- c(
  "Tumour epithelial cells",
  "Endothelial cells",
  "Mast cells",
  "Monocytes",
  "T/NK cells",
  "TAMs"
)

############################################################
# Clean / standardise CAF plots
############################################################

clean_caf_plot <- function(p) {
  p +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme(
      strip.text = element_text(
        size = 9,
        lineheight = 0.9
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        size = 9,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.text.y = element_text(size = 8.5),
      legend.position = "none",
      plot.margin = margin(
        t = 6, r = 6, b = 6, l = 6
      )
    )
}

caf_plot_list <- lapply(
  caf_targets_for_panel,
  function(tg) clean_caf_plot(caf_dotplots[[tg]])
)

names(caf_plot_list) <- caf_targets_for_panel

############################################################
# Data from the four additional metastasis-specific pairs
############################################################

new_four_df <- liana_met_only %>%
  dplyr::filter(
    (source == "Mast cells" &
       target == "CAFs") |
      
      (source == "B/plasma cells" &
         target == "Monocytes") |
      
      (source == "Tumour epithelial cells" &
         target == "Monocytes") |
      
      (source == "Mast cells" &
         target == "TAMs")
  )

############################################################
# SAME specificity and magnitude scales for ALL 10 panels
############################################################

combined_specificity_limits <- range(
  c(
    caf_all_plot_df$natmi.edge_specificity,
    new_four_df$specificity_metastasis
  ),
  na.rm = TRUE
)

combined_magnitude_limits <- range(
  c(
    caf_all_plot_df$sca.LRscore,
    new_four_df$magnitude_metastasis
  ),
  na.rm = TRUE
)

############################################################
# Sensible legend breaks based on actual data range
############################################################

combined_specificity_breaks <- pretty(
  combined_specificity_limits,
  n = 4
)

combined_specificity_breaks <-
  combined_specificity_breaks[
    combined_specificity_breaks >= combined_specificity_limits[1] &
      combined_specificity_breaks <= combined_specificity_limits[2]
  ]

combined_magnitude_breaks <- pretty(
  combined_magnitude_limits,
  n = 4
)

combined_magnitude_breaks <-
  combined_magnitude_breaks[
    combined_magnitude_breaks >= combined_magnitude_limits[1] &
      combined_magnitude_breaks <= combined_magnitude_limits[2]
  ]

combined_specificity_limits
combined_specificity_breaks

combined_magnitude_limits
combined_magnitude_breaks

############################################################
# Apply SAME scales to every plot
############################################################

add_common_liana_scales <- function(p) {
  
  p +
    
    scale_size_continuous(
      name = "Interaction\nspecificity",
      limits = combined_specificity_limits,
      breaks = combined_specificity_breaks,
      range = c(2.2, 8)
    ) +
    
    scale_colour_viridis_c(
      name = "Interaction\nmagnitude",
      limits = combined_magnitude_limits,
      breaks = combined_magnitude_breaks,
      option = "D",
      direction = -1
    )
}

caf_plot_list <- lapply(
  caf_plot_list,
  add_common_liana_scales
)


############################################################
# Clean four additional metastasis-specific plots
############################################################

clean_new_plot <- function(p) {
  
  p +
    
    facet_grid(
      . ~ source_label,
      labeller = labeller(
        source_label = label_wrap_gen(width = 15)
      )
    ) +
    
    labs(
      x = NULL,
      y = NULL
    ) +
    
    theme(
      strip.text = element_text(
        size = 9,
        lineheight = 0.9
      ),
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      axis.text.x = element_text(
        size = 9,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      
      axis.text.y = element_text(
        size = 8.5
      ),
      
      legend.position = "none",
      
      plot.margin = margin(
        t = 6,
        r = 6,
        b = 6,
        l = 6
      )
    )
}

new_plot_list <- list(
  
  clean_new_plot(
    add_common_liana_scales(
      p_met_Mast_CAF
    )
  ),
  
  clean_new_plot(
    add_common_liana_scales(
      p_met_BPlasma_Mono
    )
  ),
  
  clean_new_plot(
    add_common_liana_scales(
      p_met_Epi_Mono
    )
  ),
  
  clean_new_plot(
    add_common_liana_scales(
      p_met_Mast_TAM
    )
  )
)

shared_legend_combined <- cowplot::get_legend(
  
  add_common_liana_scales(
    p_met_Mast_CAF
  ) +
    
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
)

all_liana_plots <- c(
  
  unname(
    caf_plot_list[
      c(
        "Tumour epithelial cells",
        "Endothelial cells",
        "Mast cells",
        "Monocytes",
        "T/NK cells",
        "TAMs"
      )
    ]
  ),
  
  new_plot_list
)

all_liana_plots <- lapply(
  all_liana_plots,
  function(p) {
    p +
      theme(
        legend.position = "none"
      )
  }
)


############################################################
# Align all 10 plots
############################################################

all_liana_plots_aligned <- cowplot::align_plots(
  plotlist = all_liana_plots,
  align = "hv",
  axis = "tblr"
)

############################################################
# 2 rows × 5 columns
# No panel letters
############################################################

combined_liana_grid <- cowplot::plot_grid(
  plotlist = all_liana_plots_aligned,
  ncol = 5,
  nrow = 2,
  align = "hv",
  axis = "tblr"
)

############################################################
# Add ONE shared legend on the right
############################################################

combined_liana_panel <- cowplot::plot_grid(
  combined_liana_grid,
  shared_legend_combined,
  nrow = 1,
  rel_widths = c(
    0.91,
    0.09
  )
)

############################################################
# Preview
############################################################

combined_liana_panel

ggsave(
  filename = file.path(
    OUTDIR,
    "LIANA_metastasis_specific_LR_combined_panel.pdf"
  ),
  plot = combined_liana_panel,
  width = 16,
  height = 8.5,
  units = "in",
  device = "pdf",
  useDingbats = FALSE,
  bg = "white"
)


format_source_label <- function(x) {
  dplyr::case_when(
    x == "Tumour epithelial cells" ~ "Tumour epithelial\ncells",
    x == "B/plasma cells"          ~ "B/plasma\ncells",
    x == "Mast cells"              ~ "Mast\ncells",
    x == "CAFs"                    ~ "CAFs\n ",
    TRUE                           ~ paste0(x, "\n ")
  )
}

combined_specificity_breaks <- round(
  seq(
    combined_specificity_limits[1],
    combined_specificity_limits[2],
    length.out = 4
  ),
  2
)

combined_magnitude_breaks <- round(
  seq(
    combined_magnitude_limits[1],
    combined_magnitude_limits[2],
    length.out = 3
  ),
  2
)


library(dplyr)
library(ggplot2)

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
  
  ##########################################################
  # Order interactions
  ##########################################################
  
  plot_df <- plot_df %>%
    dplyr::arrange(
      natmi.edge_specificity,
      dplyr::desc(sca.LRscore)
    ) %>%
    dplyr::mutate(
      interaction = rev(unique(interaction))[match(interaction, rev(unique(interaction)))]
    )
  
  interaction_levels <- unique(plot_df$interaction)
  
  ##########################################################
  # Give MORE vertical spacing between rows
  ##########################################################
  
  y_map <- data.frame(
    interaction = interaction_levels,
    y_pos = seq(
      from = 1,
      by = 1.35,   # increase this to 1.4 if still too tight
      length.out = length(interaction_levels)
    )
  )
  
  plot_df <- plot_df %>%
    dplyr::left_join(
      y_map,
      by = "interaction"
    ) %>%
    dplyr::mutate(
      source_display = format_source_label(source_name),
      target = factor(target_name, levels = target_name)
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
        size = natmi.edge_specificity,
        colour = sca.LRscore
      )
    ) +
    
    facet_grid(
      . ~ source_display
    ) +
    
    scale_y_continuous(
      breaks = y_map$y_pos,
      labels = y_map$interaction,
      expand = expansion(add = c(0.6, 0.6))
    ) +
    
    scale_x_discrete(
      expand = expansion(add = c(0.55, 0.55))
    ) +
    
    scale_size_continuous(
      name = "Interaction\nspecificity",
      limits = specificity_limits,
      breaks = specificity_breaks,
      range = c(2.0, 6.0)   # smaller max size to reduce overlap
    ) +
    
    scale_colour_viridis_c(
      name = "Interaction\nmagnitude",
      limits = magnitude_limits,
      breaks = magnitude_breaks,
      option = "D",
      direction = -1
    ) +
    
    labs(
      x = NULL,
      y = NULL
    ) +
    
    theme_bw(base_size = 10) +
    
    theme(
      strip.background = element_rect(
        fill = "grey90",
        colour = "grey40",
        linewidth = 0.5
      ),
      strip.text = element_text(
        size = 8.5,
        lineheight = 0.9,
        margin = margin(3, 2, 3, 2)
      ),
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      axis.text.x = element_text(
        size = 8.5,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      
      axis.text.y = element_text(size = 8),
      
      panel.grid.major.x = element_line(colour = "grey85"),
      panel.grid.major.y = element_line(colour = "grey85"),
      panel.grid.minor = element_blank(),
      
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      
      plot.margin = margin(
        t = 5, r = 5, b = 5, l = 5
      )
    )
}

caf_targets_for_panel <- c(
  "Tumour epithelial cells",
  "Endothelial cells",
  "Mast cells",
  "Monocytes",
  "T/NK cells",
  "TAMs"
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
# ONE shared legend
############################################################

shared_legend_combined <- cowplot::get_legend(
  
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
  )
)

############################################################
# ROW 1
############################################################

row1 <- cowplot::plot_grid(
  
  caf_plot_list[["Tumour epithelial cells"]],
  caf_plot_list[["Endothelial cells"]],
  caf_plot_list[["Mast cells"]],
  caf_plot_list[["Monocytes"]],
  
  nrow = 1,
  
  rel_widths = c(
    1, 1, 1, 1
  ),
  
  align = "hv",
  axis = "tblr"
)

############################################################
# ROW 2
# Two CAF panels centred in FOUR equal slots
############################################################

row2 <- cowplot::plot_grid(
  
  NULL,
  
  caf_plot_list[["T/NK cells"]],
  caf_plot_list[["TAMs"]],
  
  NULL,
  
  nrow = 1,
  
  rel_widths = c(
    1, 1, 1, 1
  ),
  
  align = "hv",
  axis = "tblr"
)


############################################################
# ROW 3
############################################################

row3 <- cowplot::plot_grid(
  
  new_plot_list[[1]],
  new_plot_list[[2]],
  new_plot_list[[3]],
  new_plot_list[[4]],
  
  nrow = 1,
  
  rel_widths = c(
    1, 1, 1, 1
  ),
  
  align = "hv",
  axis = "tblr"
)

############################################################
# STACK ALL 3 ROWS
############################################################

combined_rows <- cowplot::plot_grid(
  
  row1,
  row2,
  row3,
  
  ncol = 1,
  
  rel_heights = c(
    1,
    1,
    1
  ),
  
  align = "v"
)


############################################################
# FINAL PANEL + SHARED LEGEND
############################################################

final_combined_panel <- cowplot::plot_grid(
  
  combined_rows,
  shared_legend_combined,
  
  nrow = 1,
  
  rel_widths = c(
    0.91,
    0.09
  )
)

final_combined_panel

ggsave(
  filename = file.path(
    OUTDIR,
    "LIANA_metastasis_specific_LR_panel_3rows.pdf"
  ),
  plot = final_combined_panel,
  width = 14,
  height = 12,
  units = "in",
  device = "pdf",
  useDingbats = FALSE,
  bg = "white"
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
  
  source_label <- stringr::str_wrap(
    source_name,
    width = 18
  )
  
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
      # Give each interaction extra vertical space
      ######################################################
      y_pos = rev(
        seq(
          from = 1,
          by = 1.35,
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
        size = natmi.edge_specificity,
        colour = sca.LRscore
      )
    ) +
    
    facet_grid(
      . ~ source,
      scales = "free_x",
      space = "free_x"
    ) +
    
    scale_x_discrete(
      expand = expansion(add = c(0.55, 0.55))
    ) +
    
    scale_y_continuous(
      breaks = plot_df$y_pos,
      labels = plot_df$interaction,
      expand = expansion(add = c(0.7, 0.7))
    ) +
    
    scale_size_continuous(
      name = "Interaction\nspecificity",
      limits = specificity_limits,
      breaks = specificity_breaks,
      labels = sprintf("%.1f", specificity_breaks),
      range = c(2.5, 6.2)
    ) +
    
    scale_colour_viridis_c(
      name = "Interaction\nmagnitude",
      limits = magnitude_limits,
      breaks = magnitude_breaks,
      labels = sprintf("%.2f", magnitude_breaks),
      option = "D",
      direction = -1
    ) +
    
    labs(
      x = NULL,
      y = NULL
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      strip.background = element_rect(
        fill = "grey90",
        colour = "grey40",
        linewidth = 0.5
      ),
      
      strip.text = element_text(
        size = 9,
        lineheight = 0.9
      ),
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      axis.text.x = element_text(
        size = 9,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      
      axis.text.y = element_text(
        size = 8.3
      ),
      
      panel.grid.major.x = element_line(
        colour = "grey85",
        linewidth = 0.35
      ),
      
      panel.grid.major.y = element_line(
        colour = "grey85",
        linewidth = 0.35
      ),
      
      panel.grid.minor = element_blank(),
      
      legend.position = if (show_legend) "right" else "none",
      
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      
      plot.margin = margin(
        t = 6, r = 6, b = 6, l = 6
      )
    )
}

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

row2 <- cowplot::plot_grid(
  
  caf_plot_list[["T/NK cells"]],
  caf_plot_list[["TAMs"]],
  NULL,
  NULL,
  
  nrow = 1,
  
  rel_widths = c(
    1, 1, 1, 1
  ),
  
  align = "hv",
  axis = "tblr"
)

row1 <- cowplot::plot_grid(
  caf_plot_list[["Tumour epithelial cells"]],
  caf_plot_list[["Endothelial cells"]],
  caf_plot_list[["Mast cells"]],
  caf_plot_list[["Monocytes"]],
  nrow = 1,
  rel_widths = c(1, 1, 1, 1),
  align = "hv",
  axis = "tblr"
)

row3 <- cowplot::plot_grid(
  new_plot_list[[1]],
  new_plot_list[[2]],
  new_plot_list[[3]],
  new_plot_list[[4]],
  nrow = 1,
  rel_widths = c(1, 1, 1, 1),
  align = "hv",
  axis = "tblr"
)


shared_legend_combined <- cowplot::get_legend(
  make_lr_plot(
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
)

combined_rows <- cowplot::plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(1, 1, 1),
  align = "v"
)

final_combined_panel <- cowplot::plot_grid(
  combined_rows,
  shared_legend_combined,
  nrow = 1,
  rel_widths = c(0.91, 0.09)
)

final_combined_panel

ggsave(
  filename = file.path(
    OUTDIR,
    "LIANA_metastasis_specific_LR_panel_3rows_leftaligned.pdf"
  ),
  plot = final_combined_panel,
  width = 14,
  height = 12,
  units = "in",
  device = "pdf",
  useDingbats = FALSE,
  bg = "white"
)




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
