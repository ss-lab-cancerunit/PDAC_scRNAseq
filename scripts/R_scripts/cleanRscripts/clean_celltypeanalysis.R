library(DESeq2); library(Seurat); library(tidyverse); library(apeglm)
library(ggplot2); library(ggrepel); library(pheatmap)
library(RColorBrewer); library(ComplexHeatmap)

ROOT <- "/Users/thirisantracy/Desktop/thesis"
OUT <- file.path(ROOT, "DEcelltypes")

pseudo_PDAC_DEcelltype <- readRDS(file.path(ROOT, "pseudo_PDAC_DEcelltype.rds"))
ct_counts <- LayerData(pseudo_PDAC_DEcelltype, assay = "decontXcounts", layer = "counts")
ct_meta <- pseudo_PDAC_DEcelltype@meta.data[colnames(ct_counts), , drop = FALSE]

required_meta <- c("dataset_GSE", "tumor", "patient", "subject_ID", "DEcelltype_annotation", "n_cells_celltype")
stopifnot(all(required_meta %in% colnames(ct_meta)))
stopifnot(identical(colnames(ct_counts), rownames(ct_meta)))

ct_meta$dataset_GSE <- as.character(ct_meta$dataset_GSE)
ct_meta$tumor <- tolower(as.character(ct_meta$tumor))
ct_meta$DEcelltype_annotation <- as.character(ct_meta$DEcelltype_annotation)

celltypes <- sort(unique(na.omit(ct_meta$DEcelltype_annotation)))
print(celltypes)
print(table(ct_meta$DEcelltype_annotation, ct_meta$tumor))

safe_name <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x))

run_celltype_DESeq2 <- function(celltype) {
  message("\nRunning: ", celltype)
  
  keep_samples <- ct_meta$DEcelltype_annotation == celltype &
    !is.na(ct_meta$DEcelltype_annotation) &
    !is.na(ct_meta$n_cells_celltype) &
    ct_meta$n_cells_celltype >= 10
  
  counts_ct <- ct_counts[, keep_samples, drop = FALSE]
  meta_ct <- ct_meta[keep_samples, , drop = FALSE]
  
  meta_ct$dataset_GSE <- droplevels(factor(meta_ct$dataset_GSE))
  meta_ct$tumor <- relevel(factor(meta_ct$tumor, levels = c("primary", "metastasis")), ref = "primary")
  meta_ct <- droplevels(meta_ct)
  
  stopifnot(identical(colnames(counts_ct), rownames(meta_ct)))
  
  sample_numbers <- table(meta_ct$tumor)
  print(sample_numbers)
  
  if (length(sample_numbers) < 2 || any(sample_numbers < 3)) {
    message("Skipped ", celltype, ": fewer than 3 pseudobulk samples in one condition.")
    return(NULL)
  }
  
  design_formula <- if (nlevels(meta_ct$dataset_GSE) > 1) ~ dataset_GSE + tumor else ~ tumor
  design_matrix <- model.matrix(design_formula, data = meta_ct)
  
  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    message("Skipped ", celltype, ": design matrix is not full rank.")
    return(NULL)
  }
  
  dds_ct <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts_ct)), colData = meta_ct, design = design_formula)
  
  keep_genes <- rowSums(counts(dds_ct) >= 10) >= 3
  cat("Samples retained:", ncol(dds_ct), "\n")
  cat("Genes before filtering:", nrow(dds_ct), "\n")
  cat("Genes retained:", sum(keep_genes), "\n")
  
  dds_ct <- dds_ct[keep_genes, ]
  
  if (nrow(dds_ct) == 0) {
    message("Skipped ", celltype, ": no genes passed filtering.")
    return(NULL)
  }
  
  dds_ct <- DESeq(dds_ct)
  res_ct <- results(dds_ct, contrast = c("tumor", "metastasis", "primary"), alpha = 0.05)
  
  coef_name <- "tumor_metastasis_vs_primary"
  if (!coef_name %in% resultsNames(dds_ct)) stop("Coefficient not found for ", celltype)
  
  res_ct_shrunk <- lfcShrink(dds_ct, coef = coef_name, type = "apeglm")
  
  test_tbl <- as.data.frame(res_ct) %>% tibble::rownames_to_column("gene") %>%
    dplyr::select(gene, baseMean, stat, pvalue, padj)
  
  lfc_tbl <- as.data.frame(res_ct_shrunk) %>% tibble::rownames_to_column("gene") %>%
    dplyr::select(gene, log2FoldChange, lfcSE)
  
  res_tbl <- left_join(test_tbl, lfc_tbl, by = "gene") %>%
    mutate(direction = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 0.58 ~ "Higher in metastasis",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -0.58 ~ "Higher in primary",
      TRUE ~ "Not significant"
    )) %>% arrange(padj)
  
  metastasis_up <- res_tbl %>% filter(direction == "Higher in metastasis")
  primary_up <- res_tbl %>% filter(direction == "Higher in primary")
  
  print(table(res_tbl$direction))
  
  ct_name <- safe_name(celltype)
  ct_out <- file.path(OUT, ct_name)
  dir.create(ct_out, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(res_tbl, file.path(ct_out, paste0(ct_name, "_DESeq2_results.csv")), row.names = FALSE)
  write.csv(metastasis_up, file.path(ct_out, paste0(ct_name, "_metastasis_up.csv")), row.names = FALSE)
  write.csv(primary_up, file.path(ct_out, paste0(ct_name, "_primary_up.csv")), row.names = FALSE)
  saveRDS(dds_ct, file.path(ct_out, paste0(ct_name, "_dds.rds")))
  saveRDS(res_ct, file.path(ct_out, paste0(ct_name, "_unshrunk_results.rds")))
  
  vsd_ct <- vst(dds_ct, blind = FALSE)
  pca_data <- plotPCA(vsd_ct, intgroup = c("dataset_GSE", "tumor"), returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  
  p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = dataset_GSE, shape = tumor)) +
    geom_point(size = 3, alpha = 0.85) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    labs(title = celltype, colour = "Dataset", shape = "Tumour") +
    theme_classic()
  
  ggsave(file.path(ct_out, paste0(ct_name, "_PCA.pdf")), p_pca, width = 7, height = 5)
  
  pdf(file.path(ct_out, paste0(ct_name, "_dispersion.pdf")), width = 7, height = 5)
  plotDispEsts(dds_ct)
  dev.off()
  
  pdf(file.path(ct_out, paste0(ct_name, "_MA.pdf")), width = 7, height = 5)
  plotMA(res_ct_shrunk, alpha = 0.05)
  dev.off()
  
  sample_cor <- cor(assay(vsd_ct))
  cor_anno <- as.data.frame(colData(vsd_ct)[, c("dataset_GSE", "tumor")])
  colnames(cor_anno) <- c("Dataset", "Tumour")
  sample_labels <- paste(colData(vsd_ct)$dataset_GSE, colData(vsd_ct)$patient, sep = "_")
  rownames(cor_anno) <- sample_labels
  rownames(sample_cor) <- sample_labels
  colnames(sample_cor) <- sample_labels
  
  pdf(file.path(ct_out, paste0(ct_name, "_sample_correlation.pdf")), width = 10, height = 9)
  pheatmap(sample_cor, annotation_col = cor_anno, annotation_row = cor_anno,
           show_colnames = FALSE, show_rownames = FALSE, border_color = NA)
  dev.off()
  
  volcano_labels <- bind_rows(
    metastasis_up %>% arrange(padj) %>% slice_head(n = 5),
    primary_up %>% arrange(padj) %>% slice_head(n = 5)
  ) %>% pull(gene)
  
  volcano_tbl <- res_tbl %>% mutate(minus_log10_padj = -log10(pmax(padj, .Machine$double.xmin)))
  
  p_volcano <- ggplot(volcano_tbl, aes(log2FoldChange, minus_log10_padj, colour = direction)) +
    geom_point(size = 1.3, alpha = 0.7) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    ggrepel::geom_text_repel(data = volcano_tbl %>% filter(gene %in% volcano_labels),
                             aes(label = gene), size = 3, max.overlaps = Inf) +
    scale_colour_manual(values = c("Higher in metastasis" = "#D73027",
                                   "Higher in primary" = "#4575B4", "Not significant" = "grey75")) +
    labs(title = celltype, x = expression(Log[2]~fold~change),
         y = expression(-Log[10]~adjusted~italic(P)), colour = NULL) +
    theme_classic()
  
  ggsave(file.path(ct_out, paste0(ct_name, "_volcano.pdf")), p_volcano, width = 8, height = 7)
  
  top_genes <- bind_rows(
    metastasis_up %>% arrange(padj) %>% slice_head(n = 10),
    primary_up %>% arrange(padj) %>% slice_head(n = 10)
  ) %>% pull(gene) %>% unique()
  
  if (length(top_genes) >= 2) {
    heatmap_mat <- assay(vsd_ct)[top_genes, , drop = FALSE]
    heat_anno <- as.data.frame(colData(dds_ct)[, c("dataset_GSE", "tumor")])
    colnames(heat_anno) <- c("Dataset", "Tumour")
    sample_labels <- paste(colData(dds_ct)$dataset_GSE, colData(dds_ct)$patient, sep = "_")
    colnames(heatmap_mat) <- sample_labels
    rownames(heat_anno) <- sample_labels
    heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
    
    p_heatmap <- ComplexHeatmap::pheatmap(heatmap_mat, name = "Z-score",
                                          color = heat_colors, scale = "row", cluster_rows = TRUE,
                                          cluster_cols = TRUE, annotation_col = heat_anno,
                                          show_rownames = TRUE, show_colnames = FALSE,
                                          border_color = NA, fontsize_row = 8)
    
    pdf(file.path(ct_out, paste0(ct_name, "_top20_heatmap.pdf")), width = 10, height = 8)
    ComplexHeatmap::draw(p_heatmap, heatmap_legend_side = "right",
                         annotation_legend_side = "right", newpage = TRUE)
    dev.off()
  }
  
  summary_tbl <- tibble(
    celltype = celltype,
    n_primary = sum(meta_ct$tumor == "primary"),
    n_metastasis = sum(meta_ct$tumor == "metastasis"),
    genes_tested = nrow(res_tbl),
    metastasis_up = nrow(metastasis_up),
    primary_up = nrow(primary_up)
  )
  
  list(celltype = celltype, res = res_ct, res_shrunk = res_ct_shrunk,
       results = res_tbl, metastasis_up = metastasis_up,
       primary_up = primary_up, summary = summary_tbl)
}

celltype_results <- setNames(lapply(celltypes, run_celltype_DESeq2), celltypes)
celltype_results <- celltype_results[!vapply(celltype_results, is.null, logical(1))]

celltype_summary <- bind_rows(lapply(celltype_results, function(x) x$summary))
print(celltype_summary)

write.csv(celltype_summary, file.path(OUT, "DEcelltype_DESeq2_summary.csv"), row.names = FALSE)
saveRDS(celltype_results, file.path(OUT, "DEcelltype_DESeq2_results.rds"))

cat("\nCompleted cell types:\n")
print(names(celltype_results))

celltype_summary
names(celltype_results)

library(dplyr)

OUT <- "/Users/thirisantracy/Desktop/thesis/DEcelltypes"

# Reload only if celltype_results is not already in your environment
if (!exists("celltype_results")) {
  celltype_results <- readRDS(file.path(OUT, "DEcelltype_DESeq2_results.rds"))
}

all_celltype_results <- dplyr::bind_rows(lapply(names(celltype_results), function(ct) {
  x <- as.data.frame(celltype_results[[ct]]$results)
  x$celltype <- ct
  x[, c("celltype", setdiff(colnames(x), "celltype")), drop = FALSE]
}))

all_celltype_results$celltype <- as.character(all_celltype_results$celltype)
all_celltype_results$direction <- as.character(all_celltype_results$direction)

all_celltype_DEGs <- all_celltype_results %>%
  dplyr::filter(direction != "Not significant")

celltype_DEG_summary <- all_celltype_DEGs %>%
  dplyr::group_by(celltype, direction) %>%
  dplyr::summarise(n_DEGs = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = direction, values_from = n_DEGs, values_fill = 0)

print(celltype_DEG_summary)

celltype_metastasis_up <- all_celltype_DEGs %>%
  dplyr::filter(direction == "Higher in metastasis") %>%
  dplyr::arrange(celltype, padj)

celltype_primary_up <- all_celltype_DEGs %>%
  dplyr::filter(direction == "Higher in primary") %>%
  dplyr::arrange(celltype, padj)

write.csv(all_celltype_results, file.path(OUT, "all_celltype_DESeq2_results.csv"), row.names = FALSE)
write.csv(all_celltype_DEGs, file.path(OUT, "all_celltype_significant_DEGs.csv"), row.names = FALSE)
write.csv(celltype_metastasis_up, file.path(OUT, "all_celltype_metastasis_up_DEGs.csv"), row.names = FALSE)
write.csv(celltype_primary_up, file.path(OUT, "all_celltype_primary_up_DEGs.csv"), row.names = FALSE)
write.csv(celltype_DEG_summary, file.path(OUT, "celltype_DEG_summary.csv"), row.names = FALSE)



library(msigdbr)
library(dplyr)
library(tibble)

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
  msigdb_sets[[nm]] <- db %>% dplyr::distinct(gs_name, gene_symbol)
}

saveRDS(msigdb_sets, "/Users/thirisantracy/Desktop/thesis/global/msigdb_sets.rds")

library(fgsea); library(tidyverse); library(ggplot2)
ROOT <- "/Users/thirisantracy/Desktop/thesis"
OUT <- file.path(ROOT, "DEcelltypes")



if (!exists("celltype_results")) celltype_results <- readRDS(file.path(OUT, "DEcelltype_DESeq2_results.rds"))
if (!exists("msigdb_sets")) msigdb_sets <- readRDS(file.path(ROOT, "global", "msigdb_sets.rds"))

safe_name <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x))

plot_celltype_fgsea <- function(x, celltype, collection) {
  if (nrow(x) == 0) return(NULL)
  plot_data <- x %>% mutate(
    pathway_label = stringr::str_remove(pathway, "^(HALLMARK_|REACTOME_|WP_|KEGG_MEDICUS_|GOBP_|GTRD_|TFT_|MODULE_|MORF_|GNF2_|IMMUNESIGDB_)"),
    pathway_label = stringr::str_replace_all(pathway_label, "_", " "),
    pathway_label = forcats::fct_reorder(pathway_label, NES)
  )
  ggplot(plot_data, aes(pathway_label, NES, fill = direction)) +
    geom_col(width = 0.75) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    coord_flip() +
    scale_fill_manual(values = c("Enriched in metastasis" = "#D73027", "Enriched in primary" = "#2166AC")) +
    labs(title = paste(celltype, "—", collection), x = NULL, y = "Normalised enrichment score", fill = NULL) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 11))
}

run_celltype_fgsea <- function(celltype) {
  message("\nRunning fgsea: ", celltype)
  ct_name <- safe_name(celltype)
  ct_out <- file.path(OUT, ct_name, "fgsea")
  dir.create(file.path(ct_out, "collections"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(ct_out, "plots"), recursive = TRUE, showWarnings = FALSE)
  
  res_ct <- celltype_results[[celltype]]$res
  rank_tbl <- as.data.frame(res_ct) %>% tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(stat), is.finite(stat)) %>% dplyr::arrange(desc(stat))
  
  ranks <- setNames(rank_tbl$stat, rank_tbl$gene)
  ranks <- sort(ranks, decreasing = TRUE)
  ranks <- ranks[!duplicated(names(ranks))]
  
  cat("Ranked genes:", length(ranks), "\n")
  
  results <- list()
  
  for (nm in names(msigdb_sets)) {
    pathways <- split(msigdb_sets[[nm]]$gene_symbol, msigdb_sets[[nm]]$gs_name)
    set.seed(123)
    results[[nm]] <- fgsea::fgseaMultilevel(
      pathways = pathways, stats = ranks, minSize = 15,
      maxSize = 500, eps = 0, nPermSimple = 1000,
      scoreType = "std"
    ) %>% as.data.frame() %>% dplyr::arrange(padj) %>%
      dplyr::mutate(collection = nm, direction = ifelse(NES > 0, "Enriched in metastasis", "Enriched in primary"))
  }
  
  significant <- lapply(results, function(x) {
    x %>% dplyr::filter(!is.na(padj), padj < 0.05) %>% dplyr::arrange(padj)
  })
  
  top <- lapply(significant, function(x) {
    top_met <- x %>% dplyr::filter(NES > 0) %>% dplyr::arrange(desc(NES)) %>% dplyr::slice_head(n = 5)
    top_pri <- x %>% dplyr::filter(NES < 0) %>% dplyr::arrange(NES) %>% dplyr::slice_head(n = 5)
    dplyr::bind_rows(top_pri, top_met) %>% dplyr::arrange(NES)
  })
  
  print(sapply(significant, nrow))
  
  saveRDS(results, file.path(ct_out, paste0(ct_name, "_fgsea_all.rds")))
  saveRDS(significant, file.path(ct_out, paste0(ct_name, "_fgsea_significant.rds")))
  saveRDS(top, file.path(ct_out, paste0(ct_name, "_fgsea_top5_each_direction.rds")))
  
  all_csv <- dplyr::bind_rows(results, .id = "database")
  sig_csv <- dplyr::bind_rows(significant, .id = "database")
  top_csv <- dplyr::bind_rows(top, .id = "database")
  
  all_csv$leadingEdge <- vapply(all_csv$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  sig_csv$leadingEdge <- vapply(sig_csv$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  top_csv$leadingEdge <- vapply(top_csv$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  
  write.csv(all_csv, file.path(ct_out, paste0(ct_name, "_fgsea_all.csv")), row.names = FALSE)
  write.csv(sig_csv, file.path(ct_out, paste0(ct_name, "_fgsea_significant.csv")), row.names = FALSE)
  write.csv(top_csv, file.path(ct_out, paste0(ct_name, "_fgsea_top5_each_direction.csv")), row.names = FALSE)
  
  for (nm in names(results)) {
    all_tbl <- results[[nm]]
    sig_tbl <- significant[[nm]]
    all_tbl$leadingEdge <- vapply(all_tbl$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    sig_tbl$leadingEdge <- vapply(sig_tbl$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    write.csv(all_tbl, file.path(ct_out, "collections", paste0(nm, "_all.csv")), row.names = FALSE)
    write.csv(sig_tbl, file.path(ct_out, "collections", paste0(nm, "_significant.csv")), row.names = FALSE)
    
    p <- plot_celltype_fgsea(top[[nm]], celltype, nm)
    if (!is.null(p)) ggsave(file.path(ct_out, "plots", paste0(nm, "_GSEA.pdf")), p, width = 9, height = 7)
  }
  
  list(all = results, significant = significant, top = top)
}

celltype_fgsea_results <- setNames(lapply(names(celltype_results), run_celltype_fgsea), names(celltype_results))
saveRDS(celltype_fgsea_results, file.path(OUT, "all_DEcelltype_fgsea_results.rds"))

celltype_fgsea_summary <- dplyr::bind_rows(lapply(names(celltype_fgsea_results), function(ct) {
  tibble::tibble(
    celltype = ct,
    collection = names(celltype_fgsea_results[[ct]]$significant),
    n_significant = sapply(celltype_fgsea_results[[ct]]$significant, nrow)
  )
}))

print(celltype_fgsea_summary)
write.csv(celltype_fgsea_summary, file.path(OUT, "DEcelltype_fgsea_summary.csv"), row.names = FALSE)

View(celltype_fgsea_results[["Tumour epithelial cells"]]$significant$Hallmark)
View(celltype_fgsea_results[["Tumour epithelial cells"]]$significant$TFT_GTRD)
View(celltype_fgsea_results[["CAFs"]]$significant$Reactome)
View(celltype_fgsea_results[["TAMs"]]$significant$GOBP)

library(tidyverse)

ROOT <- "/Users/thirisantracy/Desktop/thesis"
OUT <- file.path(ROOT, "DEcelltypes")

celltype_fgsea_results <- readRDS(file.path(OUT, "all_DEcelltype_fgsea_results.rds"))

all_fgsea_sig <- dplyr::bind_rows(lapply(names(celltype_fgsea_results), function(ct) {
  dplyr::bind_rows(lapply(names(celltype_fgsea_results[[ct]]$significant), function(db) {
    x <- celltype_fgsea_results[[ct]]$significant[[db]]
    if (nrow(x) == 0) return(NULL)
    x$celltype <- ct
    x$database <- db
    x
  }))
}))

all_fgsea_sig <- all_fgsea_sig %>%
  dplyr::select(celltype, database, pathway, NES, pval, padj, size, leadingEdge, dplyr::everything())

saveRDS(all_fgsea_sig, file.path(OUT, "all_celltype_significant_fgsea.rds"))
#################### results section in thesis

library(tidyverse)
library(ComplexHeatmap)

plot_pathway_heatmap <- function(
    database_name,
    title_name,
    remove_prefix,
    pathway_pattern = NULL
) {
  
  pathway_sig <- all_fgsea_sig %>%
    dplyr::filter(database == database_name)
  
  # Used for Gavish: retain only pathways containing "GAVISH"
  if (!is.null(pathway_pattern)) {
    pathway_sig <- pathway_sig %>%
      dplyr::filter(str_detect(
        pathway,
        regex(pathway_pattern, ignore_case = TRUE)
      ))
  }
  
  pathway_sig <- pathway_sig %>%
    dplyr::mutate(
      pathway = pathway %>%
        str_remove(remove_prefix) %>%
        str_replace_all("_", " ")
    ) %>%
    dplyr::select(pathway, celltype, NES)
  
  if (nrow(pathway_sig) == 0) {
    stop("No significant pathways found.")
  }
  
  pathway_matrix <- pathway_sig %>%
    pivot_wider(
      names_from = celltype,
      values_from = NES
    ) %>%
    column_to_rownames("pathway") %>%
    as.matrix()
  
  n_celltypes <- ncol(pathway_matrix)
  
  pathway_order <- data.frame(
    pathway = rownames(pathway_matrix),
    n_metastasis = rowSums(pathway_matrix > 0, na.rm = TRUE),
    n_primary = rowSums(pathway_matrix < 0, na.rm = TRUE),
    mean_NES = rowMeans(pathway_matrix, na.rm = TRUE)
  ) %>%
    dplyr::mutate(
      group = case_when(
        
        n_metastasis > 0 & n_primary == 0 ~
          "Consistently upregulated in metastatic tumours",
        
        n_primary > 0 & n_metastasis == 0 ~
          "Consistently downregulated in metastatic tumours",
        
        n_metastasis > 0 & n_primary > 0 ~
          "Mixed direction"
      ),
      
      group = factor(
        group,
        levels = c(
          "Consistently upregulated in metastatic tumours",
          "Mixed direction",
          "Consistently downregulated in metastatic tumours")
      )
    ) %>%
    arrange(
      group,
      desc(n_metastasis),
      desc(mean_NES)
    )
  
  pathway_matrix <- pathway_matrix[
    pathway_order$pathway,
    ,
    drop = FALSE
  ]
  
  max_NES <- max(
    abs(pathway_matrix),
    na.rm = TRUE
  )
  
  NES_limit <- ceiling(max_NES)
  
  pathway_colours <- circlize::colorRamp2(
    c(-NES_limit, 0, NES_limit),
    c("#2166AC", "white", "#B2182B")
  )
  
  pathway_plot <- Heatmap(
    pathway_matrix,
    
    name = "NES",
    col = pathway_colours,
    na_col = "grey90",
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    row_split = pathway_order$group,
    cluster_row_slices = FALSE,
    
    # Keep the heatmap block narrow
    width = grid::unit(7, "cm"),
    
    # Give each pathway enough vertical space
    height = grid::unit(
      nrow(pathway_matrix) * 0.28,
      "cm"
    ),
    
    row_title_gp = grid::gpar(
      fontsize = 7,
      fontface = "plain"
    ),
    row_title_rot = 0,
    
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 6),
    
    row_names_max_width = grid::unit(
      17,
      "cm"
    ),
    
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 8),
    
    column_title = paste0(
      "Significant ",
      title_name,
      " pathways across cell types"
    ),
    
    heatmap_legend_param = list(
      at = c(
        -NES_limit,
        -NES_limit / 2,
        0,
        NES_limit / 2,
        NES_limit
      ),
      title_gp = grid::gpar(fontsize = 9),
      labels_gp = grid::gpar(fontsize = 8),
      legend_height = grid::unit(3, "cm")
    ),
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      
      value <- pathway_matrix[i, j]
      
      if (!is.na(value)) {
        grid::grid.text(
          sprintf("%.1f", value),
          x,
          y,
          gp = grid::gpar(
            fontsize = 5.5,
            col = ifelse(
              abs(value) >= 3,
              "white",
              "black"
            )
          )
        )
      }
    }
  )
  draw(pathway_plot)
  
  invisible(
    list(
      plot = pathway_plot,
      matrix = pathway_matrix,
      order = pathway_order
    )
  )
}

kegg_heatmap <- plot_pathway_heatmap(
  database_name = "KEGG",
  title_name = "KEGG",
  remove_prefix = "^KEGG_MEDICUS_"
)


gtrd_heatmap <- plot_pathway_heatmap(
  database_name = "TFT_GTRD",
  title_name = "GTRD transcription-factor target",
  remove_prefix = "^(GTRD_|TFT_)"
)



gavish_heatmap <- plot_pathway_heatmap(
  database_name = "C4",
  title_name = "Gavish",
  remove_prefix = "^GAVISH_3CA_",
  pathway_pattern = "GAVISH"
)

hallmark_heatmap <- plot_pathway_heatmap(
  database_name = "Hallmark",
  title_name = "Hallmark",
  remove_prefix = "^HALLMARK_"
)


plot_clustered_heatmap <- function(
    heatmap_result,
    title_name
) {
  
  ########################################################
  # KEEP THE COMPLETE ORIGINAL MATRIX
  ########################################################
  
  pathway_matrix <- heatmap_result$matrix
  
  ########################################################
  # MATRIX USED ONLY TO CALCULATE CLUSTERING
  #
  # NA means not significantly enriched.
  # Replace with 0 only for distance calculation.
  ########################################################
  
  clustering_matrix <- pathway_matrix
  
  clustering_matrix[
    is.na(clustering_matrix)
  ] <- 0
  
  
  ########################################################
  # HIERARCHICAL CLUSTERING
  ########################################################
  
  row_clustering <- hclust(
    dist(clustering_matrix),
    method = "complete"
  )
  
  column_clustering <- hclust(
    dist(t(clustering_matrix)),
    method = "complete"
  )
  
  
  ########################################################
  # COLOUR SCALE
  ########################################################
  
  max_NES <- max(
    abs(pathway_matrix),
    na.rm = TRUE
  )
  
  NES_limit <- ceiling(max_NES)
  
  pathway_colours <- circlize::colorRamp2(
    c(
      -NES_limit,
      0,
      NES_limit
    ),
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )
  
  
  ########################################################
  # HEATMAP
  ########################################################
  
  pathway_plot <- ComplexHeatmap::Heatmap(
    
    pathway_matrix,
    
    name = "NES",
    
    col = pathway_colours,
    
    na_col = "grey90",
    
    
    ######################################################
    # CLUSTER BOTH DIMENSIONS
    ######################################################
    
    cluster_rows = row_clustering,
    
    cluster_columns = column_clustering,
    
    show_row_dend = TRUE,
    
    show_column_dend = TRUE,
    
    
    ######################################################
    # SIZE
    ######################################################
    
    width = grid::unit(
      7,
      "cm"
    ),
    
    height = grid::unit(
      nrow(pathway_matrix) * 0.28,
      "cm"
    ),
    
    
    ######################################################
    # ROW LABELS
    ######################################################
    
    row_names_side = "right",
    
    row_names_gp = grid::gpar(
      fontsize = 6
    ),
    
    row_names_max_width = grid::unit(
      17,
      "cm"
    ),
    
    
    ######################################################
    # COLUMN LABELS
    ######################################################
    
    column_names_rot = 90,
    
    column_names_gp = grid::gpar(
      fontsize = 8
    ),
    
    
    column_title = paste0(
      "Significant ",
      title_name,
      " pathways across cell types"
    ),
    
    
    ######################################################
    # LEGEND
    ######################################################
    
    heatmap_legend_param = list(
      
      at = c(
        -NES_limit,
        -NES_limit / 2,
        0,
        NES_limit / 2,
        NES_limit
      ),
      
      title_gp = grid::gpar(
        fontsize = 9
      ),
      
      labels_gp = grid::gpar(
        fontsize = 8
      ),
      
      legend_height = grid::unit(
        3,
        "cm"
      )
    ),
    
    
    ######################################################
    # PRINT NES VALUES
    ######################################################
    
    cell_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
    ) {
      
      value <- pathway_matrix[i, j]
      
      if (!is.na(value)) {
        
        grid::grid.text(
          
          sprintf(
            "%.1f",
            value
          ),
          
          x,
          y,
          
          gp = grid::gpar(
            
            fontsize = 5.5,
            
            col = ifelse(
              abs(value) >= 3,
              "white",
              "black"
            )
          )
        )
      }
    }
  )
  
  
  ComplexHeatmap::draw(
    pathway_plot,
    heatmap_legend_side = "right"
  )
  
  
  invisible(
    list(
      plot = pathway_plot,
      matrix = pathway_matrix,
      row_clustering = row_clustering,
      column_clustering = column_clustering
    )
  )
}
hallmark_clustered <- plot_clustered_heatmap(
  hallmark_heatmap,
  "Hallmark"
)

pdf(
  "/Users/thirisantracy/Desktop/thesis/images/Hallmark_clustered_heatmap.pdf",
  width = 11,
  height = 10
)

ComplexHeatmap::draw(
  hallmark_clustered$plot,
  heatmap_legend_side = "right"
)

dev.off()

gavish_clustered <- plot_clustered_heatmap(
  gavish_heatmap,
  "Gavish"
)

gavish_full_height <- max(
  12,
  nrow(gavish_clustered$matrix) * 0.14
)

pdf(
  "/Users/thirisantracy/Desktop/thesis/images/Gavish_full_clustered_heatmap_supplementary.pdf",
  width = 12,
  height = gavish_full_height
)

ComplexHeatmap::draw(
  gavish_clustered$plot,
  heatmap_legend_side = "right"
)

dev.off()



############################################################
# Direction-only clustered heatmap
#
# Keep pathways whose significant NES values point in only
# one direction across cell types, then cluster the retained
# pathway and cell-type profiles.
############################################################

plot_direction_only_heatmap <- function(
    heatmap_result,
    title_name
) {
  pathway_matrix <- heatmap_result$matrix
  
  direction_only <- apply(
    pathway_matrix,
    1,
    function(values) {
      values <- values[!is.na(values)]
      
      length(values) > 0 &&
        (all(values > 0) || all(values < 0))
    }
  )
  
  pathway_matrix <- pathway_matrix[
    direction_only,
    ,
    drop = FALSE
  ]
  
  if (nrow(pathway_matrix) < 2) {
    stop(
      "Fewer than two direction-consistent pathways were found for ",
      title_name,
      "."
    )
  }
  
  clustering_matrix <- pathway_matrix
  clustering_matrix[is.na(clustering_matrix)] <- 0
  
  row_clustering <- hclust(
    dist(clustering_matrix),
    method = "complete"
  )
  
  column_clustering <- hclust(
    dist(t(clustering_matrix)),
    method = "complete"
  )
  
  max_NES <- max(abs(pathway_matrix), na.rm = TRUE)
  NES_limit <- ceiling(max_NES)
  
  pathway_colours <- circlize::colorRamp2(
    c(-NES_limit, 0, NES_limit),
    c("#2166AC", "white", "#B2182B")
  )
  
  pathway_plot <- ComplexHeatmap::Heatmap(
    pathway_matrix,
    name = "NES",
    col = pathway_colours,
    na_col = "grey90",
    cluster_rows = row_clustering,
    cluster_columns = column_clustering,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    width = grid::unit(7, "cm"),
    height = grid::unit(nrow(pathway_matrix) * 0.28, "cm"),
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 6),
    row_names_max_width = grid::unit(17, "cm"),
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 8),
    column_title = paste0(
      "Direction-consistent ",
      title_name,
      " pathways across cell types"
    ),
    heatmap_legend_param = list(
      at = c(-NES_limit, -NES_limit / 2, 0, NES_limit / 2, NES_limit),
      title_gp = grid::gpar(fontsize = 9),
      labels_gp = grid::gpar(fontsize = 8),
      legend_height = grid::unit(3, "cm")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      value <- pathway_matrix[i, j]
      
      if (!is.na(value)) {
        grid::grid.text(
          sprintf("%.1f", value),
          x,
          y,
          gp = grid::gpar(
            fontsize = 5.5,
            col = ifelse(abs(value) >= 3, "white", "black")
          )
        )
      }
    }
  )
  
  ComplexHeatmap::draw(
    pathway_plot,
    heatmap_legend_side = "right"
  )
  
  invisible(
    list(
      plot = pathway_plot,
      matrix = pathway_matrix,
      row_clustering = row_clustering,
      column_clustering = column_clustering
    )
  )
}


gavish_direction_heatmap <- plot_direction_only_heatmap(
  gavish_heatmap,
  "Gavish"
)

pdf(
  "/Users/thirisantracy/Desktop/thesis/images/Gavish_direction_clustered_heatmap.pdf",
  width = 12,
  height = 14
)

ComplexHeatmap::draw(
  gavish_direction_heatmap$plot,
  heatmap_legend_side = "right"
)

dev.off()


hallmark_direction_heatmap <- plot_direction_only_heatmap(
  hallmark_heatmap,
  "Hallmark"
)

library(tidyverse)

deg_plot_data <- celltype_summary %>%
  mutate(
    total_DEGs = metastasis_up + primary_up,
    
    celltype_label = paste0(
      celltype,
      " (P=", n_primary,
      ", M=", n_metastasis, ")"
    )
  ) %>%
  arrange(total_DEGs) %>%
  mutate(
    celltype_label = factor(
      celltype_label,
      levels = celltype_label
    )
  ) %>%
  dplyr::select(
    celltype,
    celltype_label,
    n_primary,
    n_metastasis,
    total_DEGs,
    metastasis_up,
    primary_up
  ) %>%
  pivot_longer(
    cols = c(metastasis_up, primary_up),
    names_to = "direction",
    values_to = "n_DEGs"
  ) %>%
  mutate(
    plot_count = ifelse(
      direction == "primary_up",
      -n_DEGs,
      n_DEGs
    ),
    
    direction = recode(
      direction,
      "metastasis_up" = "Upregulated in metastatic tumours",
      "primary_up" = "Downregulated in metastatic tumours"
    )
  )

p_celltype_DEGs <- ggplot(
  deg_plot_data,
  aes(
    x = plot_count,
    y = celltype_label,
    fill = direction
  )
) +
  
  geom_col(width = 0.72) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_text(
    data = deg_plot_data %>%
      filter(plot_count > 0),
    aes(label = n_DEGs),
    hjust = -0.2,
    size = 3.5
  ) +
  
  geom_text(
    data = deg_plot_data %>%
      filter(plot_count < 0),
    aes(label = n_DEGs),
    hjust = 1.2,
    size = 3.5
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
  
  scale_x_continuous(
    labels = abs,
    expand = expansion(mult = c(0.10, 0.10))
  ) +
  
  labs(
    x = "Number of significant DEGs",
    y = NULL,
    fill = NULL
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

p_celltype_DEGs
####CAF gprofiler
library(dplyr)
library(gprofiler2)
library(rrvgo)
library(org.Hs.eg.db)
library(ggplot2)

ROOT <- "/Users/thirisantracy/Desktop/thesis"
OUT <- file.path(ROOT, "DEcelltypes")

# Load the cell-type DESeq2 results
celltype_results <- readRDS(
  file.path(
    OUT,
    "DEcelltype_DESeq2_results.rds"
  )
)

# Extract CAF results
caf_results <- celltype_results[["CAFs"]]$results

head(caf_results)
table(caf_results$direction)


caf_background <- caf_results %>%
  filter(!is.na(pvalue)) %>%
  pull(gene) %>%
  unique()

caf_metastasis_up <- caf_results %>%
  filter(
    direction == "Higher in metastasis"
  ) %>%
  pull(gene) %>%
  unique()

caf_primary_up <- caf_results %>%
  filter(
    direction == "Higher in primary"
  ) %>%
  pull(gene) %>%
  unique()



cat(
  "CAF ORA background:",
  length(caf_background),
  "\n"
)

cat(
  "Higher in metastasis:",
  length(caf_metastasis_up),
  "\n"
)

CAF_OUT <- file.path(
  OUT,
  "CAFs",
  "gprofiler_rrvgo"
)

dir.create(
  CAF_OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


writeLines(
  caf_background,
  file.path(
    CAF_OUT,
    "CAF_ORA_background.txt"
  )
)

writeLines(
  caf_metastasis_up,
  file.path(
    CAF_OUT,
    "CAF_metastasis_up_ORA.txt"
  )
)

writeLines(
  caf_primary_up,
  file.path(
    CAF_OUT,
    "CAF_primary_up_ORA.txt"
  )
)



run_caf_gprofiler <- function(genes) {
  
  if (length(genes) == 0) {
    return(NULL)
  }
  
  gprofiler2::gost(
    query = genes,
    organism = "hsapiens",
    ordered_query = FALSE,
    significant = TRUE,
    user_threshold = 0.05,
    correction_method = "fdr",
    domain_scope = "custom_annotated",
    custom_bg = caf_background,
    sources = c(
      "GO:BP",
      "REAC",
      "KEGG",
      "WP"
    )
  )
}
caf_gprofiler_metastasis <- run_caf_gprofiler(
  caf_metastasis_up
)

caf_gprofiler_primary <- run_caf_gprofiler(
  caf_primary_up
)


caf_gprofiler_met_tbl <-
  if (is.null(caf_gprofiler_metastasis)) {
    data.frame()
  } else {
    caf_gprofiler_metastasis$result
  }

caf_gprofiler_pri_tbl <-
  if (is.null(caf_gprofiler_primary)) {
    data.frame()
  } else {
    caf_gprofiler_primary$result
  }

dim(caf_gprofiler_met_tbl)
dim(caf_gprofiler_pri_tbl)

View(caf_gprofiler_met_tbl)
View(caf_gprofiler_pri_tbl)



flatten_list_columns <- function(x) {
  
  x[] <- lapply(
    x,
    function(column) {
      
      if (is.list(column)) {
        
        vapply(
          column,
          function(value) {
            paste(value, collapse = ";")
          },
          character(1)
        )
        
      } else {
        
        column
        
      }
    }
  )
  
  x
}

caf_gprofiler_met_csv <-
  flatten_list_columns(
    caf_gprofiler_met_tbl
  )

caf_gprofiler_pri_csv <-
  flatten_list_columns(
    caf_gprofiler_pri_tbl
  )

write.csv(
  caf_gprofiler_met_csv,
  file.path(
    CAF_OUT,
    "CAF_gprofiler_metastasis.csv"
  ),
  row.names = FALSE
)

write.csv(
  caf_gprofiler_pri_csv,
  file.path(
    CAF_OUT,
    "CAF_gprofiler_primary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  caf_gprofiler_met_tbl,
  file.path(
    CAF_OUT,
    "CAF_gprofiler_metastasis.rds"
  )
)

saveRDS(
  caf_gprofiler_pri_tbl,
  file.path(
    CAF_OUT,
    "CAF_gprofiler_primary.rds"
  )
)

caf_go_metastasis <- caf_gprofiler_met_tbl %>%
  filter(
    source == "GO:BP",
    !is.na(term_id),
    !is.na(p_value)
  ) %>%
  distinct(
    term_id,
    .keep_all = TRUE
  )

caf_go_primary <- caf_gprofiler_pri_tbl %>%
  filter(
    source == "GO:BP",
    !is.na(term_id),
    !is.na(p_value)
  ) %>%
  distinct(
    term_id,
    .keep_all = TRUE
  )


nrow(caf_go_metastasis)
nrow(caf_go_primary)

run_caf_rrvgo <- function(
    go_tbl,
    threshold = 0.7) {
  
  if (nrow(go_tbl) < 2) {
    return(NULL)
  }
  
  sim <- rrvgo::calculateSimMatrix(
    go_tbl$term_id,
    orgdb = "org.Hs.eg.db",
    ont = "BP",
    method = "Rel"
  )
  
  scores <- setNames(
    -log10(
      pmax(
        go_tbl$p_value,
        .Machine$double.xmin
      )
    ),
    go_tbl$term_id
  )
  
  reduced <- rrvgo::reduceSimMatrix(
    sim,
    scores = scores[rownames(sim)],
    threshold = threshold,
    orgdb = "org.Hs.eg.db"
  )
  
  list(
    sim = sim,
    reduced = reduced
  )
}

caf_rrvgo_metastasis <- run_caf_rrvgo(
  caf_go_metastasis
)

caf_rrvgo_primary <- run_caf_rrvgo(
  caf_go_primary
)



if (!is.null(caf_rrvgo_metastasis)) {
  
  View(
    caf_rrvgo_metastasis$reduced
  )
}

if (!is.null(caf_rrvgo_primary)) {
  
  View(
    caf_rrvgo_primary$reduced
  )
}

if (!is.null(caf_rrvgo_metastasis)) {
  
  write.csv(
    caf_rrvgo_metastasis$reduced,
    file.path(
      CAF_OUT,
      "CAF_rrvgo_metastasis.csv"
    ),
    row.names = FALSE
  )
}

if (!is.null(caf_rrvgo_primary)) {
  
  write.csv(
    caf_rrvgo_primary$reduced,
    file.path(
      CAF_OUT,
      "CAF_rrvgo_primary.csv"
    ),
    row.names = FALSE
  )
}

options(
  ggrepel.max.overlaps = 10
)

if (!is.null(caf_rrvgo_metastasis)) {
  
  p_caf_rrvgo_metastasis <-
    rrvgo::scatterPlot(
      caf_rrvgo_metastasis$sim,
      caf_rrvgo_metastasis$reduced,
      onlyParents = FALSE,
      size = "score",
      addLabel = TRUE,
      labelSize = 3
    )
  
  p_caf_rrvgo_metastasis
}


if (!is.null(caf_rrvgo_metastasis)) {
  
  ggsave(
    file.path(
      CAF_OUT,
      "CAF_rrvgo_metastasis.pdf"
    ),
    plot = p_caf_rrvgo_metastasis,
    width = 10,
    height = 8,
    units = "in"
  )
}

library(tidyverse)

# Hallmark
all_fgsea_sig %>%
  filter(
    database == "Hallmark",
    pathway %in% c(
      "HALLMARK_INTERFERON_ALPHA_RESPONSE",
      "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
      "HALLMARK_FATTY_ACID_METABOLISM"
    )
  ) %>%
  dplyr::select(
    pathway,
    celltype,
    NES,
    padj
  ) %>%
  arrange(
    pathway,
    desc(NES)
  )


all_fgsea_sig %>%
  filter(
    database == "C4",
    str_detect(
      pathway,
      regex(
        "GAVISH.*(UNFOLDED|PROTEASOM|HYPOXIA|PDAC)",
        ignore_case = TRUE
      )
    )
  ) %>%
  dplyr::select(
    pathway,
    celltype,
    NES,
    padj
  ) %>%
  arrange(
    pathway,
    desc(NES)
  )

library(tidyverse)

hallmark_selected <- c(
  "Interferon-alpha response" = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "Unfolded protein response" = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "Fatty-acid metabolism" = "HALLMARK_FATTY_ACID_METABOLISM"
)

get_hallmark_core <- function(pathway_name, pathway_label) {
  
  pathway_data <- dplyr::bind_rows(
    lapply(names(celltype_fgsea_results), function(ct) {
      
      celltype_fgsea_results[[ct]]$significant$Hallmark %>%
        dplyr::filter(pathway == pathway_name) %>%
        dplyr::select(leadingEdge) %>%
        tidyr::unnest(cols = leadingEdge) %>%
        dplyr::mutate(celltype = ct)
    })
  ) %>%
    dplyr::rename(gene = leadingEdge) %>%
    dplyr::distinct(celltype, gene)
  
  hit_celltypes <- unique(pathway_data$celltype)
  
  core_genes <- pathway_data %>%
    dplyr::count(gene, name = "n_celltypes") %>%
    dplyr::filter(
      n_celltypes == length(hit_celltypes)
    ) %>%
    dplyr::pull(gene)
  
  tidyr::expand_grid(
    pathway = pathway_label,
    gene = core_genes,
    celltype = hit_celltypes
  )
}

hallmark_core_plot_data <- dplyr::bind_rows(
  get_hallmark_core(
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "Interferon-alpha response"
  ),
  
  get_hallmark_core(
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "Unfolded protein response"
  ),
  
  get_hallmark_core(
    "HALLMARK_FATTY_ACID_METABOLISM",
    "Fatty-acid metabolism"
  )
)

pathway_order <- c(
  "Interferon-alpha response",
  "Unfolded protein response",
  "Fatty-acid metabolism"
)

hallmark_core_plot_data$pathway <- factor(
  hallmark_core_plot_data$pathway,
  levels = pathway_order
)

gene_order <- hallmark_core_plot_data %>%
  dplyr::distinct(pathway, gene) %>%
  dplyr::arrange(pathway, gene) %>%
  dplyr::pull(gene)

hallmark_core_plot_data$gene <- factor(
  hallmark_core_plot_data$gene,
  levels = rev(gene_order)
)


celltype_order <- names(celltype_fgsea_results)

hallmark_core_plot_data$celltype <- factor(
  hallmark_core_plot_data$celltype,
  levels = celltype_order
)

p_hallmark_core <- ggplot(
  hallmark_core_plot_data,
  aes(
    x = celltype,
    y = gene,
    colour = pathway
  )
) +
  geom_point(
    size = 4
  ) +
  scale_colour_manual(
    values = c(
      "Interferon-alpha response" = "#4472C4",
      "Unfolded protein response" = "#E15759",
      "Fatty-acid metabolism" = "#59A14F"
    )
  ) +
  labs(
    x = "Cell types",
    y = "Shared leading-edge gene",
    colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_line(
      colour = "grey90",
      linewidth = 0.4
    ),
    
    panel.grid.major.y = element_line(
      colour = "grey93",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank(),
    
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 9,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      margin = margin(t = 8)
    ),
    
    axis.title.y = element_text(
      size = 9,
      margin = margin(r = 8)
    ),
    
    legend.position = "right",
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.key.height = unit(
      0.45,
      "cm"
    ),
    
    legend.key.width = unit(
      0.45,
      "cm"
    ),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )

p_hallmark_core

ggsave(
  "/Users/thirisantracy/Desktop/thesis/images/hallmark_shared_leadingedge_dotplot.pdf",
  p_hallmark_core,
  width = 11,
  height = 8,
  units = "in",
  bg = "white"
)

library(tidyverse)

get_gavish_core <- function(pathway_name, pathway_label) {
  
  pathway_data <- dplyr::bind_rows(
    lapply(
      names(celltype_fgsea_results),
      function(ct) {
        
        celltype_fgsea_results[[ct]]$significant$C4 %>%
          dplyr::filter(
            pathway == pathway_name
          ) %>%
          dplyr::select(
            leadingEdge
          ) %>%
          tidyr::unnest(
            cols = leadingEdge
          ) %>%
          dplyr::mutate(
            celltype = ct
          )
      }
    )
  ) %>%
    dplyr::rename(
      gene = leadingEdge
    ) %>%
    dplyr::distinct(
      celltype,
      gene
    )
  
  hit_celltypes <- unique(
    pathway_data$celltype
  )
  
  core_genes <- pathway_data %>%
    dplyr::count(
      gene,
      name = "n_celltypes"
    ) %>%
    dplyr::filter(
      n_celltypes == length(hit_celltypes)
    ) %>%
    dplyr::pull(
      gene
    )
  
  tidyr::expand_grid(
    pathway = pathway_label,
    gene = core_genes,
    celltype = hit_celltypes
  )
}

gavish_core_plot_data <- dplyr::bind_rows(
  
  get_gavish_core(
    "GAVISH_3CA_MALIGNANT_METAPROGRAM_9_UNFOLDED_PROTEIN_RESPONSE",
    "Metaprogram 9: unfolded protein response"
  ),
  
  get_gavish_core(
    "GAVISH_3CA_MALIGNANT_METAPROGRAM_8_PROTEASOMAL_DEGRADATION",
    "Metaprogram 8: proteasomal degradation"
  )
)

gavish_pathway_order <- c(
  "Metaprogram 9: unfolded protein response",
  "Metaprogram 8: proteasomal degradation"
)

gavish_core_plot_data$pathway <- factor(
  gavish_core_plot_data$pathway,
  levels = gavish_pathway_order
)

celltype_order <- names(
  celltype_fgsea_results
)

gavish_core_plot_data$celltype <- factor(
  gavish_core_plot_data$celltype,
  levels = celltype_order
)

gavish_gene_order <- gavish_core_plot_data %>%
  dplyr::distinct(
    pathway,
    gene
  ) %>%
  dplyr::arrange(
    pathway,
    gene
  ) %>%
  dplyr::pull(
    gene
  )

gavish_core_plot_data$gene <- factor(
  gavish_core_plot_data$gene,
  levels = rev(gavish_gene_order)
)

p_gavish_core <- ggplot(
  gavish_core_plot_data,
  aes(
    x = celltype,
    y = gene,
    colour = pathway
  )
) +
  geom_point(
    size = 4
  ) +
  scale_colour_manual(
    values = c(
      "Metaprogram 9: unfolded protein response" = "#4472C4",
      "Metaprogram 8: proteasomal degradation" = "#E15759"
    )
  ) +
  labs(
    x = "Cell types",
    y = "Shared leading-edge gene",
    colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_line(
      colour = "grey90",
      linewidth = 0.4
    ),
    panel.grid.major.y = element_line(
      colour = "grey93",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 9,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      margin = margin(t = 8)
    ),
    
    axis.title.y = element_text(
      size = 9,
      margin = margin(r = 8)
    ),
    
    legend.position = "right",
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.key.height = unit(
      0.45,
      "cm"
    ),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )

p_gavish_core

ggsave(
  "/Users/thirisantracy/Desktop/thesis/images/gavish_shared_leadingedge_dotplot.pdf",
  p_gavish_core,
  width = 10,
  height = 6.5,
  units = "in",
  bg = "white"
)



gavish_tcell_pathways <- all_fgsea_sig %>%
  dplyr::filter(
    database == "C4",
    stringr::str_detect(
      pathway,
      stringr::regex(
        "CD4|CD8|CYTOTOX|DYSFUNCTION|REG",
        ignore_case = TRUE
      )
    )
  ) %>%
  dplyr::distinct(pathway) %>%
  dplyr::arrange(pathway)

gavish_tcell_pathways


all_fgsea_sig %>%
  dplyr::filter(
    database == "C4",
    pathway %in% gavish_tcell_pathways$pathway
  ) %>%
  dplyr::select(
    pathway,
    celltype,
    NES,
    padj
  ) %>%
  dplyr::arrange(pathway, celltype)


tcell_gavish_pathways <- c(
  "CD4 T-cell regulatory" =
    "GAVISH_3CA_METAPROGRAM_CD4_T_CELLS_T_REG",
  
  "CD8 T-cell cytotoxic" =
    "GAVISH_3CA_METAPROGRAM_CD8_T_CELLS_CYTOTOXIC",
  
  "CD4 T-cell cytotoxic" =
    "GAVISH_3CA_METAPROGRAM_CD4_T_CELLS_CYTOTOXIC",
  
  "CD8 T-cell dysfunction" =
    "GAVISH_3CA_METAPROGRAM_CD8_T_CELLS_DYSFUNCTION"
)


tcell_leadingedge <- dplyr::bind_rows(
  lapply(
    names(tcell_gavish_pathways),
    function(label) {
      
      pathway_id <- tcell_gavish_pathways[[label]]
      
      celltype_fgsea_results[["T/NK cells"]]$significant$C4 %>%
        dplyr::filter(
          pathway == pathway_id
        ) %>%
        dplyr::select(
          leadingEdge
        ) %>%
        tidyr::unnest(
          cols = leadingEdge
        ) %>%
        dplyr::mutate(
          programme = label
        )
    }
  )
) %>%
  dplyr::rename(
    gene = leadingEdge
  ) %>%
  dplyr::distinct(
    programme,
    gene
  )

View(tcell_leadingedge)

tcell_leadingedge %>%
  dplyr::filter(
    gene == "FOXP3"
  )

programme_order <- c(
  "CD4 T-cell regulatory",
  "CD8 T-cell cytotoxic",
  "CD4 T-cell cytotoxic",
  "CD8 T-cell dysfunction"
)

tcell_leadingedge$programme <- factor(
  tcell_leadingedge$programme,
  levels = programme_order
)

gene_order <- tcell_leadingedge %>%
  dplyr::count(
    gene,
    name = "n_programmes"
  ) %>%
  dplyr::arrange(
    n_programmes,
    gene
  ) %>%
  dplyr::pull(
    gene
  )

tcell_leadingedge$gene <- factor(
  tcell_leadingedge$gene,
  levels = gene_order
)

p_tcell_leadingedge <- ggplot(
  tcell_leadingedge,
  aes(
    x = programme,
    y = gene,
    colour = programme
  )
) +
  geom_point(
    size = 3
  ) +
  scale_colour_manual(
    values = c(
      "CD4 T-cell regulatory" = "#7B61A8",
      "CD8 T-cell cytotoxic" = "#D95F4C",
      "CD4 T-cell cytotoxic" = "#E6A23C",
      "CD8 T-cell dysfunction" = "#4C78A8"
    ))+ 
  labs(
    x = "Gavish metaprogram of T cells ",
    y = "Leading-edge genes",
    colour = NULL,
    
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(
      colour = "grey92",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(
      size = 10,
      colour = "black",
      margin = margin(t = 8)
    ),
    
    axis.title.y = element_text(
      size = 10,
      colour = "black",
      margin = margin(r = 8)
    ),
    
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    
    legend.position = "right",
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.key.height = unit(
      0.5,
      "cm"
    ),
    
    plot.margin = margin(
      8, 12, 8, 8
    )
  )

p_tcell_leadingedge

ggsave(
  "/Users/thirisantracy/Desktop/thesis/images/Tcell_Gavish_leadingedge.pdf",
  plot = p_tcell_leadingedge,
  width = 9,
  height = 9,
  units = "in",
  bg = "white"
)


caf_go_rrvgo_table <- caf_go_metastasis %>%
  dplyr::select(
    term_id,
    term_name,
    p_value,
    intersection_size,
    term_size
  ) %>%
  dplyr::left_join(
    caf_rrvgo_metastasis$reduced %>%
      dplyr::select(
        go,
        cluster,
        parent,
        parentTerm,
        score
      ),
    by = c("term_id" = "go")
  ) %>%
  dplyr::arrange(p_value)

View(caf_go_rrvgo_table)

caf_go_rrvgo_representatives <- caf_go_rrvgo_table %>%
  dplyr::filter(term_id == parent) %>%
  dplyr::arrange(p_value)

View(caf_go_rrvgo_representatives)
hallmark_up_all <- all_fgsea_sig %>%
  dplyr::filter(
    database == "Hallmark",
    NES > 0,
    padj < 0.05
  ) %>%
  dplyr::select(
    celltype,
    pathway,
    NES,
    padj,
    size
  ) %>%
  dplyr::arrange(
    celltype,
    dplyr::desc(NES)
  )

View(hallmark_up_all)

hallmark_up_by_celltype <- hallmark_up_all %>%
  split(.$celltype)

names(hallmark_up_by_celltype)

View(hallmark_up_by_celltype[["Tumour epithelial cells"]])
View(hallmark_up_by_celltype[["CAFs"]])
View(hallmark_up_by_celltype[["TAMs"]])
View(hallmark_up_by_celltype[["Monocytes"]])
View(hallmark_up_by_celltype[["Endothelial cells"]])
View(hallmark_up_by_celltype[["Mast cells"]])
View(hallmark_up_by_celltype[["T/NK cells"]])
View(hallmark_up_by_celltype[["B/plasma cells"]])


######upset plot for the pathways

UPSET_OUT <- "/Users/thirisantracy/Desktop/thesis/images/leading_edge_upsets"

dir.create(
  UPSET_OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


safe_file <- function(x) {
  gsub(
    "^_|_$",
    "",
    gsub(
      "[^A-Za-z0-9]+",
      "_",
      x
    )
  )
}


plot_pathway_leadingedge_upset <- function(
    label,
    database,
    pathway
) {
  
  ########################################################
  # Get significant metastasis-enriched result
  # for this pathway in each cell type
  ########################################################
  
  pathway_results <- all_fgsea_sig %>%
    dplyr::filter(
      database == .env$database,
      pathway == .env$pathway,
      !is.na(padj),
      padj < 0.05,
      NES > 0
    )
  
  if (nrow(pathway_results) == 0) {
    stop(
      "No significant metastasis-enriched results found for: ",
      label
    )
  }
  
  ########################################################
  # Extract leading-edge genes separately per cell type
  ########################################################
  
  leading_edges <- setNames(
    pathway_results$leadingEdge,
    pathway_results$celltype
  )
  
  leading_edges <- lapply(
    leading_edges,
    function(x) {
      unique(
        as.character(x)
      )
    }
  )
  
  leading_edges <- leading_edges[
    lengths(leading_edges) > 0
  ]
  
  if (length(leading_edges) < 2) {
    stop(
      label,
      " is significantly metastasis-enriched in fewer than 2 cell types."
    )
  }
  
  ########################################################
  # UpSet matrix
  ########################################################
  
  overlap_matrix <- ComplexHeatmap::make_comb_mat(
    leading_edges
  )
  
  ########################################################
  # Genes shared across ALL included cell types
  ########################################################
  
  shared_all <- Reduce(
    intersect,
    leading_edges
  )
  
  ########################################################
  # Table showing recurrence of each gene
  ########################################################
  
  gene_overlap <- dplyr::bind_rows(
    lapply(
      names(leading_edges),
      function(ct) {
        
        tibble::tibble(
          celltype = ct,
          gene = leading_edges[[ct]]
        )
      }
    )
  ) %>%
    dplyr::distinct(
      celltype,
      gene
    ) %>%
    dplyr::group_by(
      gene
    ) %>%
    dplyr::summarise(
      n_celltypes = dplyr::n_distinct(celltype),
      
      celltypes = paste(
        sort(
          unique(celltype)
        ),
        collapse = "; "
      ),
      
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      dplyr::desc(n_celltypes),
      gene
    )
  
  ########################################################
  # Plot
  ########################################################
  
  upset_plot <- ComplexHeatmap::UpSet(
    
    overlap_matrix,
    
    set_order = names(
      leading_edges
    ),
    
    comb_order = order(
      ComplexHeatmap::comb_degree(
        overlap_matrix
      ),
      ComplexHeatmap::comb_size(
        overlap_matrix
      ),
      decreasing = TRUE
    ),
    
    top_annotation =
      ComplexHeatmap::upset_top_annotation(
        overlap_matrix,
        add_numbers = TRUE,
        numbers_rot = 0
      ),
    
    right_annotation =
      ComplexHeatmap::upset_right_annotation(
        overlap_matrix,
        add_numbers = TRUE
      ),
    
    column_title = paste0(
      label,
      "\nLeading-edge genes shared across metastatic cell types"
    )
  )
  
  ########################################################
  # Save PDF
  ########################################################
  
  file_name <- safe_file(label)
  
  pdf(
    file.path(
      UPSET_OUT,
      paste0(
        file_name,
        "_UpSet.pdf"
      )
    ),
    width = 10,
    height = 7
  )
  
  ComplexHeatmap::draw(
    upset_plot
  )
  
  dev.off()
  
  ########################################################
  # Save gene-overlap table
  ########################################################
  
  write.csv(
    gene_overlap,
    file.path(
      UPSET_OUT,
      paste0(
        file_name,
        "_gene_overlap.csv"
      )
    ),
    row.names = FALSE
  )
  
  ########################################################
  # Return everything
  ########################################################
  
  invisible(
    list(
      label = label,
      pathway = pathway,
      database = database,
      leading_edges = leading_edges,
      shared_all = shared_all,
      gene_overlap = gene_overlap,
      overlap_matrix = overlap_matrix,
      plot = upset_plot
    )
  )
}

hallmark_interferon_upset <-
  plot_pathway_leadingedge_upset(
    
    label =
      "Hallmark - Interferon alpha response",
    
    database =
      "Hallmark",
    
    pathway =
      "HALLMARK_INTERFERON_ALPHA_RESPONSE"
  )

hallmark_interferon_upset$shared_all

View(
  hallmark_interferon_upset$gene_overlap
)

hallmark_upr_upset <-
  plot_pathway_leadingedge_upset(
    
    label =
      "Hallmark - Unfolded protein response",
    
    database =
      "Hallmark",
    
    pathway =
      "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
  )

hallmark_upr_upset$shared_all

hallmark_fattyacid_upset <-
  plot_pathway_leadingedge_upset(
    
    label =
      "Hallmark - Fatty acid metabolism",
    
    database =
      "Hallmark",
    
    pathway =
      "HALLMARK_FATTY_ACID_METABOLISM"
  )
hallmark_fattyacid_upset$shared_all
gavish_proteasome_upset <-
  plot_pathway_leadingedge_upset(
    
    label =
      "Gavish - Malignant proteasomal degradation",
    
    database =
      "C4",
    
    pathway =
      "GAVISH_3CA_MALIGNANT_METAPROGRAM_8_PROTEASOMAL_DEGRADATION"
  )

gavish_proteasome_upset$shared_all

gavish_upr_upset <-
  plot_pathway_leadingedge_upset(
    
    label =
      "Gavish - Malignant unfolded protein response",
    
    database =
      "C4",
    
    pathway =
      "GAVISH_3CA_MALIGNANT_METAPROGRAM_9_UNFOLDED_PROTEIN_RESPONSE"
  )

#########rerunning progeny and viper

########################################################
# PROGENy + VIPER
# FINAL ANALYSIS
#
# Multiple-testing correction:
# BH separately within each cell type
########################################################


########################################################
# 1. Packages and paths
########################################################

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

library(decoupleR)
library(progeny)

library(viper)
library(dorothea)

library(ComplexHeatmap)
library(circlize)
library(grid)


ROOT <- "/Users/thirisantracy/Desktop/thesis"

OUT <- file.path(
  ROOT,
  "DEcelltypes"
)

ACTIVITY_OUT <- file.path(
  OUT,
  "PROGENy_VIPER_per_celltype_BH"
)

FIG_OUT <- file.path(
  ROOT,
  "images"
)

dir.create(
  ACTIVITY_OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIG_OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


########################################################
# 2. Load cell-type DESeq2 results FRESH
########################################################

celltype_results <- readRDS(
  file.path(
    OUT,
    "DEcelltype_DESeq2_results.rds"
  )
)

names(celltype_results)

length(celltype_results)
########################################################
# 3. Function to make DESeq2-stat matrix
########################################################

make_stat_matrix <- function(celltype) {
  
  res <- as.data.frame(
    celltype_results[[celltype]]$res
  )
  
  
  stat <- res$stat
  
  names(stat) <- rownames(res)
  
  
  # Remove missing/infinite values
  stat <- stat[
    !is.na(stat) &
      is.finite(stat)
  ]
  
  
  # Remove duplicate gene names if present
  stat <- stat[
    !duplicated(names(stat))
  ]
  
  
  mat <- matrix(
    stat,
    ncol = 1,
    dimnames = list(
      names(stat),
      celltype
    )
  )
  
  
  return(mat)
}

test_mat <- make_stat_matrix(
  "Tumour epithelial cells"
)

dim(test_mat)

head(test_mat)

range(test_mat)

########################################################
# 4. PROGENy network
########################################################

progeny_net <- progeny::getModel(
  organism = "Human",
  top = 500,
  decoupleR = TRUE
)


progeny_net <- progeny_net %>%
  dplyr::transmute(
    source = as.character(variable),
    target = as.character(gene),
    weight = mor * likelihood
  )


head(progeny_net)

colnames(progeny_net)

dim(progeny_net)

unique(progeny_net$source)
########################################################
# 5. Run PROGENy separately for each cell type
#    + BH correction WITHIN that cell type
########################################################

run_progeny_celltype <- function(celltype) {
  
  message(
    "Running PROGENy: ",
    celltype
  )
  
  
  ######################################################
  # DESeq2 statistics for this cell type
  ######################################################
  
  mat <- make_stat_matrix(
    celltype
  )
  
  
  ######################################################
  # PROGENy activity inference
  ######################################################
  
  result <- decoupleR::run_mlm(
    mat = mat,
    network = progeny_net,
    .source = "source",
    .target = "target",
    .mor = "weight",
    minsize = 5
  )
  
  
  ######################################################
  # BH correction WITHIN this cell type
  ######################################################
  
  result <- result %>%
    
    dplyr::mutate(
      
      celltype = celltype,
      
      padj = p.adjust(
        p_value,
        method = "BH"
      ),
      
      direction = dplyr::case_when(
        
        score > 0 ~
          "Higher activity in metastatic tumours",
        
        score < 0 ~
          "Lower activity in metastatic tumours",
        
        TRUE ~
          "No direction"
      )
    )
  
  
  return(result)
}

########################################################
# 6. Run all cell types
########################################################

progeny_results <- dplyr::bind_rows(
  
  lapply(
    names(celltype_results),
    run_progeny_celltype
  )
)

dim(progeny_results)

colnames(progeny_results)

head(progeny_results)
########################################################
# 7. Verify per-cell-type BH correction
########################################################

progeny_BH_check <- progeny_results %>%
  
  dplyr::group_by(celltype) %>%
  
  dplyr::summarise(
    
    maximum_difference =
      max(
        abs(
          padj -
            p.adjust(
              p_value,
              method = "BH"
            )
        ),
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

progeny_BH_check


########################################################
# 8. Significant PROGENy activities
#    BH-adjusted P < 0.05 within each cell type
########################################################

progeny_sig <- progeny_results %>%
  
  dplyr::filter(
    !is.na(padj),
    padj < 0.05
  ) %>%
  
  dplyr::select(
    celltype,
    source,
    score,
    p_value,
    padj,
    direction
  ) %>%
  
  dplyr::arrange(
    celltype,
    padj
  )

progeny_sig %>%
  dplyr::count(
    celltype,
    direction
  ) %>%
  print(n = Inf)


View(
  progeny_sig
)

saveRDS(
  progeny_results,
  file.path(
    ACTIVITY_OUT,
    "PROGENy_all_per_celltype_BH.rds"
  )
)

saveRDS(
  progeny_sig,
  file.path(
    ACTIVITY_OUT,
    "PROGENy_significant_per_celltype_BH.rds"
  )
)


write.csv(
  progeny_results,
  file.path(
    ACTIVITY_OUT,
    "PROGENy_all_per_celltype_BH.csv"
  ),
  row.names = FALSE
)

write.csv(
  progeny_sig,
  file.path(
    ACTIVITY_OUT,
    "PROGENy_significant_per_celltype_BH.csv"
  ),
  row.names = FALSE
)

########################################################
# 9. Clustered PROGENy heatmap
########################################################

celltype_order <- names(
  celltype_results
)


plot_clustered_progeny_heatmap <- function(
    selected_pathways,
    title_text
) {
  
  selected_pathways <- unique(
    selected_pathways
  )
  
  
  ######################################################
  # DISPLAY MATRIX:
  # significant activity only
  ######################################################
  
  display_data <- tidyr::expand_grid(
    source = selected_pathways,
    celltype = celltype_order
  ) %>%
    
    dplyr::left_join(
      
      progeny_sig %>%
        dplyr::filter(
          source %in% selected_pathways
        ) %>%
        dplyr::select(
          source,
          celltype,
          score
        ),
      
      by = c(
        "source",
        "celltype"
      )
    )
  
  
  display_matrix <- display_data %>%
    
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = score
    ) %>%
    
    tibble::column_to_rownames(
      "source"
    ) %>%
    
    as.matrix()
  
  
  display_matrix <- display_matrix[
    selected_pathways,
    celltype_order,
    drop = FALSE
  ]
  
  
  ######################################################
  # CLUSTERING MATRIX:
  # all actual activity scores
  ######################################################
  
  clustering_data <- tidyr::expand_grid(
    source = selected_pathways,
    celltype = celltype_order
  ) %>%
    
    dplyr::left_join(
      
      progeny_results %>%
        dplyr::filter(
          source %in% selected_pathways
        ) %>%
        dplyr::select(
          source,
          celltype,
          score
        ),
      
      by = c(
        "source",
        "celltype"
      )
    )
  
  
  clustering_matrix <- clustering_data %>%
    
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = score
    ) %>%
    
    tibble::column_to_rownames(
      "source"
    ) %>%
    
    as.matrix()
  
  
  clustering_matrix <- clustering_matrix[
    selected_pathways,
    celltype_order,
    drop = FALSE
  ]
  
  
  # Only needed if an activity estimate was not returned
  clustering_matrix[
    is.na(clustering_matrix)
  ] <- 0
  
  
  ######################################################
  # Hierarchical clustering
  ######################################################
  
  row_clustering <- hclust(
    dist(clustering_matrix),
    method = "complete"
  )
  
  
  column_clustering <- hclust(
    dist(t(clustering_matrix)),
    method = "complete"
  )
  
  
  ######################################################
  # Colour scale
  ######################################################
  
  max_score <- max(
    abs(display_matrix),
    na.rm = TRUE
  )
  
  score_limit <- ceiling(
    max_score
  )
  
  
  activity_colours <- circlize::colorRamp2(
    c(
      -score_limit,
      0,
      score_limit
    ),
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )
  
  
  ######################################################
  # Heatmap
  ######################################################
  
  heatmap <- ComplexHeatmap::Heatmap(
    
    display_matrix,
    
    name = "PROGENy\nactivity score",
    
    col = activity_colours,
    
    na_col = "grey92",
    width = grid::unit(
      ncol(display_matrix) * 0.70,
      "cm"
    ),
    
    height = grid::unit(
      nrow(display_matrix) * 0.55,
      "cm"
    ),
    
    # hierarchical clustering
    cluster_rows = row_clustering,
    cluster_columns = column_clustering,
    
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    
    
    row_names_side = "left",
    
    row_names_gp = grid::gpar(
      fontsize = 9
    ),
    
    
    column_names_rot = 90,
    
    column_names_gp = grid::gpar(
      fontsize = 8
    ),
    
    
    column_title = title_text,
    
    column_title_gp = grid::gpar(
      fontsize = 11
    ),
    
    
    rect_gp = grid::gpar(
      col = "white",
      lwd = 0.5
    ),
    
    
    heatmap_legend_param = list(
      
      at = c(
        -score_limit,
        -score_limit / 2,
        0,
        score_limit / 2,
        score_limit
      )
    ),
    
    
    cell_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
    ) {
      
      value <- display_matrix[i, j]
      
      if (!is.na(value)) {
        
        grid::grid.text(
          
          sprintf(
            "%.1f",
            value
          ),
          
          x,
          y,
          
          gp = grid::gpar(
            
            fontsize = 7,
            
            col = ifelse(
              abs(value) >=
                0.6 * score_limit,
              "white",
              "black"
            )
          )
        )
      }
    }
  )
  
  
  invisible(
    list(
      plot = heatmap,
      matrix = display_matrix,
      clustering_matrix = clustering_matrix,
      row_clustering = row_clustering,
      column_clustering = column_clustering
    )
  )
}

progeny_main_pathways <- progeny_sig %>%
  dplyr::distinct(source) %>%
  dplyr::pull(source)

length(
  progeny_main_pathways
)

progeny_main_pathways

progeny_heatmap_main_clustered <-
  plot_clustered_progeny_heatmap(
    
    selected_pathways =
      progeny_main_pathways,
    
    title_text =
      "PROGENy pathway activity in metastatic versus primary PDAC"
  )
ComplexHeatmap::draw(
  progeny_heatmap_main_clustered$plot,
  heatmap_legend_side = "right"
)

pdf(
  file.path(
    FIG_OUT,
    "PROGENy_clustered_per_celltype_BH.pdf"
  ),
  width = 10,
  height = 7
)

ComplexHeatmap::draw(
  progeny_heatmap_main_clustered$plot,
  heatmap_legend_side = "right"
)

dev.off()


########################################################
# 10. DoRothEA network for VIPER
########################################################

data(
  dorothea_hs_pancancer,
  package = "dorothea"
)


table(
  dorothea_hs_pancancer$confidence
)


dorothea_net <- dorothea_hs_pancancer %>%
  
  dplyr::filter(
    confidence %in% c(
      "A",
      "B",
      "C"
    )
  ) %>%
  
  dplyr::transmute(
    source = tf,
    target = target,
    mor = mor
  )


dim(
  dorothea_net
)

head(
  dorothea_net
)

length(
  unique(
    dorothea_net$source
  )
)

########################################################
# 11. Run VIPER separately for each cell type
#     + BH correction WITHIN that cell type
########################################################

run_viper_celltype <- function(celltype) {
  
  message(
    "Running VIPER: ",
    celltype
  )
  
  
  ######################################################
  # DESeq2 statistics for this cell type
  ######################################################
  
  mat <- make_stat_matrix(
    celltype
  )
  
  
  ######################################################
  # VIPER activity inference
  ######################################################
  
  result <- decoupleR::run_viper(
    
    mat = mat,
    
    network = dorothea_net,
    
    .source = "source",
    
    .target = "target",
    
    .mor = "mor",
    
    minsize = 5,
    
    pleiotropy = TRUE,
    
    eset.filter = FALSE,
    
    verbose = FALSE
  )
  
  
  ######################################################
  # BH correction WITHIN this cell type
  ######################################################
  
  result <- result %>%
    
    dplyr::mutate(
      
      celltype = celltype,
      
      padj = p.adjust(
        p_value,
        method = "BH"
      ),
      
      direction = dplyr::case_when(
        
        score > 0 ~
          "Higher activity in metastatic tumours",
        
        score < 0 ~
          "Lower activity in metastatic tumours",
        
        TRUE ~
          "No direction"
      )
    )
  
  
  return(result)
}

########################################################
# 12. Run all cell types
########################################################

viper_results <- dplyr::bind_rows(
  
  lapply(
    names(celltype_results),
    run_viper_celltype
  )
)

dim(
  viper_results
)

colnames(
  viper_results
)

head(
  viper_results
)

########################################################
# 13. Verify BH within each cell type
########################################################

viper_BH_check <- viper_results %>%
  
  dplyr::group_by(celltype) %>%
  
  dplyr::summarise(
    
    maximum_difference =
      max(
        abs(
          padj -
            p.adjust(
              p_value,
              method = "BH"
            )
        ),
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

viper_BH_check

########################################################
# 14. Significant VIPER activities
########################################################

viper_sig <- viper_results %>%
  
  dplyr::filter(
    !is.na(padj),
    padj < 0.05
  ) %>%
  
  dplyr::select(
    celltype,
    source,
    score,
    p_value,
    padj,
    direction
  ) %>%
  
  dplyr::arrange(
    celltype,
    padj
  )

nrow(
  viper_sig
)


viper_sig %>%
  
  dplyr::count(
    celltype,
    direction
  ) %>%
  
  dplyr::arrange(
    celltype,
    direction
  ) %>%
  
  print(
    n = Inf
  )

saveRDS(
  viper_results,
  file.path(
    ACTIVITY_OUT,
    "VIPER_all_per_celltype_BH.rds"
  )
)

saveRDS(
  viper_sig,
  file.path(
    ACTIVITY_OUT,
    "VIPER_significant_per_celltype_BH.rds"
  )
)


write.csv(
  viper_results,
  file.path(
    ACTIVITY_OUT,
    "VIPER_all_per_celltype_BH.csv"
  ),
  row.names = FALSE
)

write.csv(
  viper_sig,
  file.path(
    ACTIVITY_OUT,
    "VIPER_significant_per_celltype_BH.csv"
  ),
  row.names = FALSE
)

########################################################
# 15. Clustered VIPER heatmap
########################################################

########################################################
# Clustered VIPER heatmap
# Significant values displayed
# Full scores used for clustering
# Values printed using vectorised layer_fun
########################################################

plot_clustered_viper_heatmap <- function(
    selected_TFs,
    title_text
) {
  
  selected_TFs <- unique(
    selected_TFs
  )
  
  
  if (length(selected_TFs) < 2) {
    stop(
      "Fewer than two TFs selected."
    )
  }
  
  
  ######################################################
  # DISPLAY MATRIX:
  # significant values only
  ######################################################
  
  display_data <- tidyr::expand_grid(
    source = selected_TFs,
    celltype = celltype_order
  ) %>%
    
    dplyr::left_join(
      
      viper_sig %>%
        dplyr::filter(
          source %in% selected_TFs
        ) %>%
        dplyr::select(
          source,
          celltype,
          score
        ),
      
      by = c(
        "source",
        "celltype"
      )
    )
  
  
  display_matrix <- display_data %>%
    
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = score
    ) %>%
    
    tibble::column_to_rownames(
      "source"
    ) %>%
    
    as.matrix()
  
  
  display_matrix <- display_matrix[
    selected_TFs,
    celltype_order,
    drop = FALSE
  ]
  
  
  ######################################################
  # CLUSTERING MATRIX:
  # use ALL VIPER activity scores
  ######################################################
  
  clustering_data <- tidyr::expand_grid(
    source = selected_TFs,
    celltype = celltype_order
  ) %>%
    
    dplyr::left_join(
      
      viper_results %>%
        dplyr::filter(
          source %in% selected_TFs
        ) %>%
        dplyr::select(
          source,
          celltype,
          score
        ),
      
      by = c(
        "source",
        "celltype"
      )
    )
  
  
  clustering_matrix <- clustering_data %>%
    
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = score
    ) %>%
    
    tibble::column_to_rownames(
      "source"
    ) %>%
    
    as.matrix()
  
  
  clustering_matrix <- clustering_matrix[
    selected_TFs,
    celltype_order,
    drop = FALSE
  ]
  
  
  ######################################################
  # Only if an activity estimate was not returned
  ######################################################
  
  clustering_matrix[
    is.na(clustering_matrix)
  ] <- 0
  
  
  ######################################################
  # Hierarchical clustering
  ######################################################
  
  row_clustering <- hclust(
    dist(clustering_matrix),
    method = "complete"
  )
  
  
  column_clustering <- hclust(
    dist(t(clustering_matrix)),
    method = "complete"
  )
  
  
  ######################################################
  # Symmetric score scale
  ######################################################
  
  max_score <- max(
    abs(display_matrix),
    na.rm = TRUE
  )
  
  
  score_limit <- ceiling(
    max_score
  )
  
  
  activity_colours <- circlize::colorRamp2(
    c(
      -score_limit,
      0,
      score_limit
    ),
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )
  
  
  ######################################################
  # Heatmap
  ######################################################
  
  heatmap <- ComplexHeatmap::Heatmap(
    
    display_matrix,
    
    name = "VIPER\nactivity score",
    
    col = activity_colours,
    
    # Non-significant combinations
    na_col = "grey92",
    
    
    ####################################################
    # Compact tile proportions
    ####################################################
    
    width = grid::unit(
      ncol(display_matrix) * 0.75,
      "cm"
    ),
    
    height = grid::unit(
      nrow(display_matrix) * 0.46,
      "cm"
    ),
    
    ####################################################
    # Clustering
    ####################################################
    
    cluster_rows = row_clustering,
    
    cluster_columns = column_clustering,
    
    show_row_dend = TRUE,
    
    show_column_dend = TRUE,
    
    
    ####################################################
    # Row labels
    ####################################################
    
    row_names_side = "right",
    
    row_names_gp = grid::gpar(
      
      fontsize = ifelse(
        length(selected_TFs) > 50,
        6,
        8
      )
    ),
    
    
    ####################################################
    # Column labels
    ####################################################
    
    column_names_rot = 90,
    
    column_names_gp = grid::gpar(
      fontsize = 8
    ),
    
    
    ####################################################
    # Title
    ####################################################
    
    column_title = title_text,
    
    column_title_gp = grid::gpar(
      fontsize = 11
    ),
    
    
    ####################################################
    # Cell borders
    ####################################################
    
    rect_gp = grid::gpar(
      col = "white",
      lwd = 0.5
    ),
    
    
    ####################################################
    # Legend
    ####################################################
    
    heatmap_legend_param = list(
      
      at = c(
        -score_limit,
        -score_limit / 2,
        0,
        score_limit / 2,
        score_limit
      )
    ),
    
    
    ####################################################
    # VECTORISED VALUE LABELS
    #
    # Same purpose as cell_fun, but much faster for
    # >100 TFs.
    ####################################################
    
    ####################################################
    # Activity-score labels
    ####################################################
    
    cell_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
    ) {
      
      value <- display_matrix[i, j]
      
      if (!is.na(value)) {
        
        grid::grid.text(
          
          sprintf(
            "%.1f",
            value
          ),
          
          x,
          y,
          
          gp = grid::gpar(
            
            fontsize = 6.5,
            
            col = ifelse(
              
              abs(value) >=
                score_limit * 0.5,
              
              "white",
              
              "black"
            )
          )
        )
      }
    }
  )
  
  ######################################################
  # Return
  ######################################################
  
  invisible(
    list(
      plot = heatmap,
      matrix = display_matrix,
      clustering_matrix = clustering_matrix,
      row_clustering = row_clustering,
      column_clustering = column_clustering
    )
  )
}
########################################################
# 16. Main VIPER TF selection
########################################################

viper_main_TFs <- viper_sig %>%
  
  dplyr::filter(
    score > 0
  ) %>%
  
  dplyr::distinct(
    source
  ) %>%
  
  dplyr::pull(
    source
  )


length(
  viper_main_TFs
)

viper_main_TFs

viper_heatmap_main_clustered <-
  plot_clustered_viper_heatmap(
    
    selected_TFs =
      viper_main_TFs,
    
    title_text =
      "Transcription-factor activity increased in metastatic PDAC"
  )


pdf(
  file.path(
    FIG_OUT,
    "VIPER_main_clustered_per_celltype_BH.pdf"
  ),
  width = 8,
  height = max(
    7,
    main_height_in + 3
  )
)

ComplexHeatmap::draw(
  viper_heatmap_main_clustered$plot,
  heatmap_legend_side = "right"
)

dev.off()

########################################################
# 17. Supplementary VIPER
########################################################

viper_supplementary_TFs <- viper_sig %>%
  
  dplyr::distinct(
    source
  ) %>%
  
  dplyr::pull(
    source
  )


length(
  viper_supplementary_TFs
)

viper_heatmap_supplementary_clustered <-
  plot_clustered_viper_heatmap(
    
    selected_TFs =
      viper_supplementary_TFs,
    
    title_text =
      "Significant transcription-factor activity across cell types"
  )


pdf(
  file.path(
    FIG_OUT,
    "VIPER_full_clustered_supplementary_per_celltype_BH.pdf"
  ),
  width = 10,
  height = max(
    8,
    nrow(
      viper_heatmap_supplementary_clustered$matrix
    ) * 0.28
  )
)

ComplexHeatmap::draw(
  viper_heatmap_supplementary_clustered$plot,
  heatmap_legend_side = "right"
)

dev.off()

########################################################
# FINAL CHECKS
########################################################

cat(
  "\nPROGENy significant results:",
  nrow(progeny_sig),
  "\n"
)

cat(
  "PROGENy pathways significant in >=1 cell type:",
  length(progeny_main_pathways),
  "\n\n"
)


cat(
  "VIPER significant results:",
  nrow(viper_sig),
  "\n"
)

cat(
  "VIPER TFs in main figure:",
  length(viper_main_TFs),
  "\n"
)

cat(
  "VIPER TFs in supplementary figure:",
  length(viper_supplementary_TFs),
  "\n"
)


progeny_sig %>%
  dplyr::count(
    celltype,
    direction
  ) %>%
  print(n = Inf)


viper_sig %>%
  dplyr::count(
    celltype,
    direction
  ) %>%
  print(n = Inf)
