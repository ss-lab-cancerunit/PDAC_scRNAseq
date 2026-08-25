library(DESeq2)
library(Seurat)
library(tidyverse)
library(SeuratObject)
library(Matrix)
library(pheatmap)
library(fgsea)
library(msigdbr)
library(EnhancedVolcano)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(apeglm)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)
library(enrichplot)
library(ComplexHeatmap)
library(grid)
library(gprofiler2)
library(EnhancedVolcano)


pseudo_PDAC_all <- readRDS("/Users/thirisantracy/Desktop/thesis/pseudo_PDAC_all.rds")
pb_counts <- LayerData(pseudo_PDAC_all, assay = "decontXcounts", layer = "counts")
pb_meta <- pseudo_PDAC_all@meta.data[colnames(pb_counts), , drop = FALSE]

required_meta <- c("dataset_GSE", "tumor", "patient", "subject_ID", "n_cells_global")
stopifnot(all(required_meta %in% colnames(pb_meta)))
stopifnot(identical(colnames(pb_counts), rownames(pb_meta)))

pb_meta$dataset_GSE <- factor(pb_meta$dataset_GSE)
pb_meta$tumor <- factor(pb_meta$tumor, levels = c("primary", "metastasis"))
pb_meta$tumor <- relevel(factor(pb_meta$tumor), ref = "primary")


keep_samples <- !is.na(pb_meta$n_cells_global) & pb_meta$n_cells_global >= 10
cat("Samples before filtering:", ncol(pb_counts), "\n")
cat("Samples removed:", sum(!keep_samples), "\n")

pb_counts <- pb_counts[, keep_samples, drop = FALSE]
pb_meta <- droplevels(pb_meta[keep_samples, , drop = FALSE])

stopifnot(identical(colnames(pb_counts), rownames(pb_meta)))
print(table(pb_meta$dataset_GSE, pb_meta$tumor))

design_matrix <- model.matrix(~ dataset_GSE + tumor, data = pb_meta)
stopifnot(qr(design_matrix)$rank == ncol(design_matrix))


dds_global <- DESeqDataSetFromMatrix(
  countData = pb_counts,
  colData = pb_meta,
  design = ~ dataset_GSE + tumor
)

keep_genes <- rowSums(counts(dds_global) >= 10) >= 3
table(keep_genes)
cat("Genes before filtering:", nrow(dds_global), "\n")
cat("Genes retained:", sum(keep_genes), "\n")
dds_global <- dds_global[keep_genes, ]
dds_global <- DESeq(dds_global)

res_global <- results(
  dds_global,
  contrast = c("tumor", "metastasis", "primary"),
  alpha = 0.05
)

coef_name <- "tumor_metastasis_vs_primary"
stopifnot(coef_name %in% resultsNames(dds_global))

res_global_shrunk <- lfcShrink(dds_global, coef = coef_name, type = "apeglm")

test_tbl <- as.data.frame(res_global) %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, baseMean, stat, pvalue, padj)

lfc_tbl <- as.data.frame(res_global_shrunk) %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, log2FoldChange, lfcSE)

res_global_tbl <- left_join(test_tbl, lfc_tbl, by = "gene") %>%
  mutate(direction = case_when(
    !is.na(padj) & padj < 0.05 & log2FoldChange > 0.58 ~ "Higher in metastasis",
    !is.na(padj) & padj < 0.05 & log2FoldChange < -0.58 ~ "Higher in primary",
    TRUE ~ "Not significant"
  )) %>%
  arrange(padj)

print(table(res_global_tbl$direction))

metastasis_up_global <- res_global_tbl %>% filter(direction == "Higher in metastasis")
primary_up_global <- res_global_tbl %>% filter(direction == "Higher in primary")

write.csv(res_global_tbl, "/Users/thirisantracy/Desktop/thesis/global_DESeq2_results.csv", row.names = FALSE)
write.csv(metastasis_up_global, "/Users/thirisantracy/Desktop/thesis/global_metastasis_up.csv", row.names = FALSE)
write.csv(primary_up_global, "/Users/thirisantracy/Desktop/thesis/global_primary_up.csv", row.names = FALSE)
saveRDS(dds_global, "/Users/thirisantracy/Desktop/thesis/dds_global.rds")


vsd_global <- vst(dds_global, blind = FALSE)

pca_data <- plotPCA(vsd_global, intgroup = c("dataset_GSE", "tumor"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_global_pca <- ggplot(pca_data, aes(PC1, PC2, colour = dataset_GSE, shape = tumor)) +
  geom_point(size = 3, alpha = 0.85) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  labs(colour = "Dataset", shape = "Tumour") +
  theme_classic()

p_global_pca
volcano_labels <- dplyr::bind_rows(
  res_global_tbl %>%
    dplyr::filter(direction == "Higher in metastasis") %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 5),
  
  res_global_tbl %>%
    dplyr::filter(direction == "Higher in primary") %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 5)
) %>%
  dplyr::pull(gene)

volcano_labels



# Number of genes in each category
n_metastasis <- sum(
  res_global_tbl$direction == "Higher in metastasis"
)

n_primary <- sum(
  res_global_tbl$direction == "Higher in primary"
)

n_not_sig <- sum(
  res_global_tbl$direction == "Not significant"
)

# Legend labels
legend_met <- paste0(
  "Upregulated in metastatic tumours (n=",
  n_metastasis,
  ")"
)

legend_primary <- paste0(
  "Downregulated in metastatic tumours (n=",
  n_primary,
  ")"
)

legend_ns <- paste0(
  "Not significant (n=",
  n_not_sig,
  ")"
)

volcano_cols <- ifelse(
  res_global_tbl$direction == "Higher in metastasis",
  "#D73027",
  ifelse(
    res_global_tbl$direction == "Higher in primary",
    "#4575B4",
    "grey75"
  )
)

names(volcano_cols) <- ifelse(
  res_global_tbl$direction == "Higher in metastasis",
  legend_met,
  ifelse(
    res_global_tbl$direction == "Higher in primary",
    legend_primary,
    legend_ns
  )
)


p_global_volcano <- EnhancedVolcano(
  res_global_tbl,
  
  lab = res_global_tbl$gene,
  selectLab = volcano_labels,
  
  x = "log2FoldChange",
  y = "padj",
  
  pCutoff = 0.05,
  FCcutoff = 0.58,
  
  colCustom = volcano_cols,
  colAlpha = 0.7,
  
  pointSize = 1.5,
  labSize = 3.5,
  
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  colConnectors = "grey40",
  arrowheads = FALSE,
  
  xlab = expression(
    "Log"[2] * " fold change"
  ),
  
  ylab = expression(
    -log[10] * "(adjusted " * italic(P) * ")"
  ),
  
  axisLabSize = 12,
  
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  legendDropLevels = TRUE,
  legendPosition = "right",
  legendLabSize = 10,
  legendIconSize = 3,
  
  gridlines.major = TRUE,
  gridlines.minor = FALSE,
  
  cutoffLineType = "dashed",
  cutoffLineWidth = 0.5,
  cutoffLineCol = "grey40",
  
  border = "partial",
  borderWidth = 0.7,
  borderColour = "black"
)

p_global_volcano <- p_global_volcano +
  scale_colour_manual(
    values = setNames(
      c(
        "#D73027",
        "#4575B4",
        "grey75"
      ),
      c(
        legend_met,
        legend_primary,
        legend_ns
      )
    ),
    breaks = c(
      legend_met,
      legend_primary,
      legend_ns
    )
  )

graphics.off()
p_global_volcano
ggsave(
  "/Users/thirisantracy/Desktop/thesis/images/global_volcano_final.pdf",
  plot = p_global_volcano,
  width = 8,
  height = 7,
  units = "in",
  bg = "white"
)
global_background <- res_global_tbl %>%
  dplyr::filter(!is.na(pvalue)) %>%
  dplyr::pull(gene) %>%
  unique()

global_metastasis_up <- metastasis_up_global$gene
global_primary_up <- primary_up_global$gene
cat("ORA background:", length(global_background), "\n")
cat("Higher in metastasis:", length(global_metastasis_up), "\n")
cat("Higher in primary:", length(global_primary_up), "\n")


writeLines(global_background, "/Users/thirisantracy/Desktop/thesis/global/global_ORA_background.txt")
writeLines(global_metastasis_up, "/Users/thirisantracy/Desktop/thesis/global/global_metastasis_up_ORA.txt")
writeLines(global_primary_up, "/Users/thirisantracy/Desktop/thesis/global/global_primary_up_ORA.txt")
OUT <- "/Users/thirisantracy/Desktop/thesis/global"
run_gprofiler <- function(genes) {
  if (length(genes) == 0) return(NULL)
  gost(query = genes, organism = "hsapiens", ordered_query = FALSE, significant = TRUE,
       user_threshold = 0.05, correction_method = "fdr",
       domain_scope = "custom_annotated", custom_bg = global_background,
       sources = c("GO:BP", "REAC", "KEGG", "WP"))
}

gprofiler_metastasis <- run_gprofiler(global_metastasis_up)
gprofiler_primary <- run_gprofiler(global_primary_up)

gprofiler_met_tbl <- if (is.null(gprofiler_metastasis)) data.frame() else gprofiler_metastasis$result
gprofiler_pri_tbl <- if (is.null(gprofiler_primary)) data.frame() else gprofiler_primary$result

list_columns <- names(gprofiler_met_tbl)[vapply(gprofiler_met_tbl, is.list, logical(1))]
print(list_columns)

flatten_list_columns <- function(x) {
  x[] <- lapply(x, function(column) {
    if (is.list(column)) vapply(column, function(value) paste(value, collapse = ";"), character(1)) else column
  })
  x
}

gprofiler_met_csv <- flatten_list_columns(gprofiler_met_tbl)
gprofiler_pri_csv <- flatten_list_columns(gprofiler_pri_tbl)

write.csv(gprofiler_met_csv, "/Users/thirisantracy/Desktop/thesis/global/global_gprofiler_metastasis.csv", row.names = FALSE)
write.csv(gprofiler_pri_csv, "/Users/thirisantracy/Desktop/thesis/global/global_gprofiler_primary.csv", row.names = FALSE)

saveRDS(gprofiler_met_tbl, "/Users/thirisantracy/Desktop/thesis/global/global_gprofiler_metastasis.rds")
saveRDS(gprofiler_pri_tbl, "/Users/thirisantracy/Desktop/thesis/global/global_gprofiler_primary.rds")


go_metastasis <- gprofiler_met_tbl %>% filter(source == "GO:BP", !is.na(term_id), !is.na(p_value)) %>% distinct(term_id, .keep_all = TRUE)
go_primary <- gprofiler_pri_tbl %>% filter(source == "GO:BP", !is.na(term_id), !is.na(p_value)) %>% distinct(term_id, .keep_all = TRUE)

run_rrvgo <- function(go_tbl, threshold = 0.7) {
  if (nrow(go_tbl) < 2) return(NULL)
  sim <- rrvgo::calculateSimMatrix(go_tbl$term_id, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
  scores <- setNames(-log10(pmax(go_tbl$p_value, .Machine$double.xmin)), go_tbl$term_id)
  reduced <- rrvgo::reduceSimMatrix(sim, scores = scores[rownames(sim)], threshold = threshold, orgdb = "org.Hs.eg.db")
  list(sim = sim, reduced = reduced)
}

rrvgo_metastasis <- run_rrvgo(go_metastasis)
rrvgo_primary <- run_rrvgo(go_primary)


if (!is.null(rrvgo_metastasis)) write.csv(rrvgo_metastasis$reduced, file.path(OUT, "global_rrvgo_metastasis.csv"), row.names = FALSE)
if (!is.null(rrvgo_primary)) write.csv(rrvgo_primary$reduced, file.path(OUT, "global_rrvgo_primary.csv"), row.names = FALSE)


if (!is.null(rrvgo_metastasis)) {
  p_rrvgo_metastasis <- rrvgo::scatterPlot(rrvgo_metastasis$sim, rrvgo_metastasis$reduced, onlyParents = TRUE, size = "score", addLabel = TRUE, labelSize = 2)
  p_rrvgo_metastasis
}

if (!is.null(rrvgo_primary)) {
  p_rrvgo_primary <- rrvgo::scatterPlot(rrvgo_primary$sim, rrvgo_primary$reduced, onlyParents = TRUE, size = "score", addLabel = TRUE, labelSize = 2)
  p_rrvgo_primary
}
collections <- tribble(
  ~name, ~collection, ~subcollection,
  "Hallmark", "H", NA_character_,
  "Reactome", "C2", "CP:REACTOME",
  "WikiPathways", "C2", "CP:WIKIPATHWAYS",
  "KEGG", "C2", "CP:KEGG_MEDICUS",
  "TFT_GTRD", "C3", "TFT:GTRD",
  "TFT_Legacy", "C3", "TFT:TFT_LEGACY",
  "C4", "C4", NA_character_,
  "GOBP", "C5", "GO:BP",
  "C6", "C6", NA_character_,
  "C7", "C7", "IMMUNESIGDB"
)

msigdb_sets <- list()

for (i in seq_len(nrow(collections))) {
  nm <- collections$name[i]
  if (is.na(collections$subcollection[i])) {
    db <- msigdbr(db_species = "HS", species = "Homo sapiens", collection = collections$collection[i])
  } else {
    db <- msigdbr(db_species = "HS", species = "Homo sapiens", collection = collections$collection[i], subcollection = collections$subcollection[i])
  }
  msigdb_sets[[nm]] <- db %>% distinct(gs_name, gene_symbol)
}

sapply(msigdb_sets, \(x) n_distinct(x$gs_name))


fgsea_tbl <- res_global %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(!is.na(stat), is.finite(stat)) %>%
  arrange(desc(stat))


fgsea_ranks <- setNames(fgsea_tbl$stat, fgsea_tbl$gene)
fgsea_ranks <- sort(fgsea_ranks, decreasing = TRUE)
cat("Ranked genes:", length(fgsea_ranks), "\n")
cat("Duplicated gene names:", anyDuplicated(names(fgsea_ranks)), "\n")
head(fgsea_ranks)
tail(fgsea_ranks)


fgsea_results <- list()

for (nm in names(msigdb_sets)) {
  pathways <- split(msigdb_sets[[nm]]$gene_symbol, msigdb_sets[[nm]]$gs_name)
  set.seed(123)
  fgsea_results[[nm]] <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = fgsea_ranks,
    minSize = 15,
    maxSize = 500,
    eps = 0,
    nPermSimple = 1000,
    scoreType = "std"
  ) %>%
    as.data.frame() %>%
    arrange(padj) %>%
    mutate(
      collection = nm,
      direction = ifelse(NES > 0, "Enriched in metastasis", "Enriched in primary")
    )
}

fgsea_sig <- lapply(fgsea_results, function(x) {
  x %>% filter(!is.na(padj), padj < 0.05) %>% arrange(padj)
})

sapply(fgsea_sig, nrow)


OUT <- "/Users/thirisantracy/Desktop/thesis/global"

saveRDS(fgsea_results, file.path(OUT, "global_fgsea_all_results.rds"))
saveRDS(fgsea_sig, file.path(OUT, "global_fgsea_significant_results.rds"))

fgsea_all_csv <- bind_rows(fgsea_results, .id = "database")
fgsea_sig_csv <- bind_rows(fgsea_sig, .id = "database")

fgsea_all_csv$leadingEdge <- vapply(fgsea_all_csv$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
fgsea_sig_csv$leadingEdge <- vapply(fgsea_sig_csv$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))

write.csv(fgsea_all_csv, file.path(OUT, "global_fgsea_all_results.csv"), row.names = FALSE)
write.csv(fgsea_sig_csv, file.path(OUT, "global_fgsea_significant_results.csv"), row.names = FALSE)

for (nm in names(fgsea_results)) {
  all_tbl <- fgsea_results[[nm]]
  sig_tbl <- fgsea_sig[[nm]]
  all_tbl$leadingEdge <- vapply(all_tbl$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  sig_tbl$leadingEdge <- vapply(sig_tbl$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  write.csv(all_tbl, file.path(OUT, "fgsea_collections", paste0(nm, "_all.csv")), row.names = FALSE)
  write.csv(sig_tbl, file.path(OUT, "fgsea_collections", paste0(nm, "_significant.csv")), row.names = FALSE)
}
plot_fgsea_hallmark <- function(x) {
  
  plot_data <- x %>%
    mutate(
      pathway_label = pathway %>%
        str_remove("^HALLMARK_") %>%
        str_replace_all("_", " "),
      
      direction = factor(
        direction,
        levels = c(
          "Enriched in metastasis",
          "Enriched in primary"
        ),
        labels = c(
          "Upregulated in metastatic tumours",
          "Downregulated in metastatic tumours"
        )
      ),
      
      pathway_label = forcats::fct_reorder(
        pathway_label,
        NES
      )
    )
  
  ggplot(
    plot_data,
    aes(
      x = NES,
      y = pathway_label,
      fill = direction
    )
  ) +
    geom_col(
      width = 0.72
    ) +
    
    geom_vline(
      xintercept = 0,
      linewidth = 0.4
    ) +
    
    scale_fill_manual(
      values = c(
        "Upregulated in metastatic tumours" = "#B2182B",
        "Downregulated in metastatic tumours" = "#2166AC"
      ),
      breaks = c(
        "Upregulated in metastatic tumours",
        "Downregulated in metastatic tumours"
      )
    ) +
    
    labs(
      x = "Normalised enrichment score",
      y = NULL,
      fill = NULL
    ) +
    
    theme_classic() +
    
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    )
}

p_hallmark_all <- plot_fgsea_hallmark(
  fgsea_sig$Hallmark
)

p_hallmark_all
View(rrvgo_metastasis$reduced)

options(ggrepel.max.overlaps = 10)
p_rrvgo_metastasis <- rrvgo::scatterPlot(
  rrvgo_metastasis$sim,
  rrvgo_metastasis$reduced,
  onlyParents = FALSE,
  size = "score",
  addLabel = TRUE,
  labelSize = 3
)

p_rrvgo_metastasis

ggsave(
  file.path(
    OUT,
    "global_rrvgo_metastasisthesis.pdf"
  ),
  plot = p_rrvgo_metastasis,
  width = 10,
  height = 8,
  units = "in"
)
# Significant Hallmark pathways enriched in metastasis
hallmark_met_top10 <- fgsea_sig$Hallmark %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10) %>%
  mutate(
    Pathway = pathway) %>%
  dplyr::select(
    Pathway,
    NES,
    padj,
    size
  )

hallmark_met_top10


gprofiler_met_top10 <- rrvgo_metastasis$reduced %>%
  
  # Keep only the representative term from each rrvgo cluster
  filter(go == parent) %>%
  
  # Add the original g:Profiler statistics back in
  left_join(
    go_metastasis %>%
      dplyr::select(
        term_id,
        p_value,
        intersection_size,
        term_size
      ),
    by = c("go" = "term_id")
  ) %>%
  
  # Most significant representative terms first
  arrange(p_value) %>%
  
  slice_head(n = 10) %>%
  
  dplyr::select(
    `GO biological process` = term,
    `GO ID` = go,
    `Adjusted p-value` = p_value,
    `Intersection size` = intersection_size,
    `Term size` = term_size
  )

View(gprofiler_met_top10)
