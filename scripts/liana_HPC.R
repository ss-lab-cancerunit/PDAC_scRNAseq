library(Seurat)
library(SingleCellExperiment)
library(liana)
library(dplyr)
library(readr)

ROOT <- "/rds/general/user/ts3225/ephemeral/R/QC4"
OUTDIR <- file.path(ROOT, "liana_results")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

primary_liana_sce <- readRDS(file.path(ROOT, "primary_liana_sce.rds"))
met_liana_sce <- readRDS(file.path(ROOT, "met_liana_sce.rds"))


primary_liana <- liana_wrap(
  sce = primary_liana_sce,
  idents_col = "liana_celltype",
  resource = "Consensus",
  method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
  expr_prop = 0.1,
  min_cells = 10,
  assay.type = "logcounts",
  base = exp(1)
)

met_liana <- liana_wrap(
  sce = met_liana_sce,
  idents_col = "liana_celltype",
  resource = "Consensus",
  method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
  expr_prop = 0.1,
  min_cells = 10,
  assay.type = "logcounts",
  base = exp(1)
)

primary_res <- liana_aggregate(primary_liana)
met_res <- liana_aggregate(met_liana)

saveRDS(primary_liana, file.path(OUTDIR, "primary_liana_raw.rds"))
saveRDS(met_liana, file.path(OUTDIR, "met_liana_raw.rds"))
saveRDS(primary_res, file.path(OUTDIR, "primary_liana_aggregated.rds"))
saveRDS(met_res, file.path(OUTDIR, "met_liana_aggregated.rds"))

write_csv(primary_res, file.path(OUTDIR, "primary_liana_aggregated.csv"))
write_csv(met_res, file.path(OUTDIR, "met_liana_aggregated.csv"))
