library(Seurat)
ROOT <- "/rds/general/user/ts3225/ephemeral/R/QC4"
primary_rds <- file.path(ROOT, "PDAC_primary_filtered.rds")
met_rds <- file.path(ROOT, "PDAC_met_filtered.rds")

PDAC_primary <- readRDS(primary_rds)
PDAC_met <- readRDS(met_rds)

required_meta <- c(
  "dataset_GSE",
  "tumor",
  "patient",
  "subject_ID",
  "celltype_annotation",
  "DEcelltype_annotation",
  "n_cells_global",
  "n_cells_celltype"
)

PDAC_combined <- merge(
PDAC_primary,
y = PDAC_met,
add.cell.ids = c("primary", "met")
)
rm(PDAC_primary, PDAC_met)
DefaultAssay(PDAC_combined) <- "decontXcounts"


PDAC_combined <- JoinLayers(
PDAC_combined,
assay = "decontXsymbols"
)

global_info <- unique(
  PDAC_combined@meta.data[
    ,
    c(
      "dataset_GSE",
      "tumor",
      "patient",
      "subject_ID",
      "n_cells_global"
    )
  ]
)

DEcelltype_info <- unique(
  PDAC_combined@meta.data[
    !is.na(PDAC_combined$DEcelltype_annotation),
    c(
      "dataset_GSE",
      "tumor",
      "patient",
      "subject_ID",
      "DEcelltype_annotation",
      "n_cells_celltype"
    )
  ]
)

pseudo_PDAC_all <- AggregateExpression(
PDAC_combined,
assays = "decontXcounts",
return.seurat = TRUE,
group.by = c("dataset_GSE","tumor", "patient")
)

pseudo_PDAC_celltype <- AggregateExpression(
PDAC_combined,
assays = "decontXcounts",
return.seurat = TRUE,
group.by = c("dataset_GSE", "tumor", "patient", "celltype_annotation")
)

pseudo_PDAC_DEcelltype <- AggregateExpression(
PDAC_combined,
assays = "decontXcounts",
return.seurat = TRUE,
group.by = c("dataset_GSE", "tumor", "patient", "DEcelltype_annotation")
)

global_original_key <- paste(
  global_info$dataset_GSE,
  global_info$tumor,
  global_info$patient,
  sep = "__"
)

global_pb_key <- paste(
  pseudo_PDAC_all$dataset_GSE,
  pseudo_PDAC_all$tumor,
  pseudo_PDAC_all$patient,
  sep = "__"
)

global_match <- match(
  global_pb_key,
  global_original_key
)

stopifnot(!anyNA(global_match))

pseudo_PDAC_all$subject_ID <-
  global_info$subject_ID[global_match]

pseudo_PDAC_all$n_cells_global <-
  global_info$n_cells_global[global_match]

head(
  pseudo_PDAC_all@meta.data[
    ,
    c(
      "dataset_GSE",
      "tumor",
      "patient",
      "subject_ID",
      "n_cells_global"
    )
  ]
)

DE_original_key <- paste(
  DEcelltype_info$dataset_GSE,
  DEcelltype_info$tumor,
  DEcelltype_info$patient,
  DEcelltype_info$DEcelltype_annotation,
  sep = "__"
)

DE_pb_key <- paste(
  pseudo_PDAC_DEcelltype$dataset_GSE,
  pseudo_PDAC_DEcelltype$tumor,
  pseudo_PDAC_DEcelltype$patient,
  pseudo_PDAC_DEcelltype$DEcelltype_annotation,
  sep = "__"
)

DE_match <- match(
  DE_pb_key,
  DE_original_key
)

stopifnot(!anyNA(DE_match))


pseudo_PDAC_DEcelltype$subject_ID <-
  DEcelltype_info$subject_ID[DE_match]

pseudo_PDAC_DEcelltype$n_cells_celltype <-
  DEcelltype_info$n_cells_celltype[DE_match]


head(
  pseudo_PDAC_DEcelltype@meta.data[
    ,
    c(
      "dataset_GSE",
      "tumor",
      "patient",
      "subject_ID",
      "DEcelltype_annotation",
      "n_cells_celltype"
    )
  ]
)

detailed_info <- unique(
  PDAC_combined@meta.data[
    ,
    c(
      "dataset_GSE",
      "tumor",
      "patient",
      "subject_ID",
      "celltype_annotation"
    )
  ]
)

detailed_original_key <- paste(
  detailed_info$dataset_GSE,
  detailed_info$tumor,
  detailed_info$patient,
  detailed_info$celltype_annotation,
  sep = "__"
)

detailed_pb_key <- paste(
  pseudo_PDAC_celltype$dataset_GSE,
  pseudo_PDAC_celltype$tumor,
  pseudo_PDAC_celltype$patient,
  pseudo_PDAC_celltype$celltype_annotation,
  sep = "__"
)

detailed_match <- match(
  detailed_pb_key,
  detailed_original_key
)

stopifnot(!anyNA(detailed_match))

pseudo_PDAC_celltype$subject_ID <-
  detailed_info$subject_ID[detailed_match]




rm(PDAC_combined)

saveRDS(
  pseudo_PDAC_all,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/pseudo_PDAC_all.rds",
  compress = FALSE)

saveRDS(
  pseudo_PDAC_celltype,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/pseudo_PDAC_celltype.rds",
  compress = FALSE)

saveRDS(
  pseudo_PDAC_DEcelltype,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/pseudo_PDAC_DEcelltype.rds",
  compress = FALSE)


