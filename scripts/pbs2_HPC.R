library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(SingleCellExperiment)
library(celda)
library(decontX)
library(scDblFinder)
library(scater)
library(scran)
library(patchwork)
library(dplyr)
library(harmony)
library(SingleR)
library(celldex)
library(scRNAseq)
library(scater)
library(dplyr)
library(broom)
library(DESeq2)
library(ggrepel)
library(data.table)
library(Matrix)
library(BiocParallel)
out_dir <- "/rds/general/user/ts3225/ephemeral/R/QC4"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

LIN_DIR <- "/rds/general/user/ts3225/ephemeral/R/Lin"

P01 <- Read10X(file.path(LIN_DIR, "P1"))
P02 <- Read10X(file.path(LIN_DIR, "P2"))
P03 <- Read10X(file.path(LIN_DIR, "P3"))
P04 <- Read10X(file.path(LIN_DIR, "P4"))
P05 <- Read10X(file.path(LIN_DIR, "P5"))
P06 <- Read10X(file.path(LIN_DIR, "P6"))
P07 <- Read10X(file.path(LIN_DIR, "P7"))
P08 <- Read10X(file.path(LIN_DIR, "P8"))
P09 <- Read10X(file.path(LIN_DIR, "P9"))
P10 <- Read10X(file.path(LIN_DIR, "P10"))

M01 <- Read10X(file.path(LIN_DIR, "M1"))
M02 <- Read10X(file.path(LIN_DIR, "M2"))
M03 <- Read10X(file.path(LIN_DIR, "M4"))
M04 <- Read10X(file.path(LIN_DIR, "M5"))
M05 <- Read10X(file.path(LIN_DIR, "M6"))

P01 <- CreateSeuratObject(counts = P01, project = "P01", min.cells = 3, min.features = 200)
P02 <- CreateSeuratObject(counts = P02, project = "P02", min.cells = 3, min.features = 200)
P03 <- CreateSeuratObject(counts = P03, project = "P03", min.cells = 3, min.features = 200)
P04 <- CreateSeuratObject(counts = P04, project = "P04", min.cells = 3, min.features = 200)
P05 <- CreateSeuratObject(counts = P05, project = "P05", min.cells = 3, min.features = 200)
P06 <- CreateSeuratObject(counts = P06, project = "P06", min.cells = 3, min.features = 200)
P07 <- CreateSeuratObject(counts = P07, project = "P07", min.cells = 3, min.features = 200)
P08 <- CreateSeuratObject(counts = P08, project = "P08", min.cells = 3, min.features = 200)
P09 <- CreateSeuratObject(counts = P09, project = "P09", min.cells = 3, min.features = 200)
P10 <- CreateSeuratObject(counts = P10, project = "P10", min.cells = 3, min.features = 200)
M01 <- CreateSeuratObject(counts = M01, project = "M01", min.cells = 3, min.features = 200)
M02 <- CreateSeuratObject(counts = M02, project = "M02", min.cells = 3, min.features = 200)
M03 <- CreateSeuratObject(counts = M03, project = "M03", min.cells = 3, min.features = 200)
M04 <- CreateSeuratObject(counts = M04, project = "M04", min.cells = 3, min.features = 200)
M05 <- CreateSeuratObject(counts = M05, project = "M05", min.cells = 3, min.features = 200)


P01$tumor <- "primary"
P02$tumor <- "primary"
P03$tumor <- "primary"
P04$tumor <- "primary"
P05$tumor <- "primary"
P06$tumor <- "primary"
P07$tumor <- "primary"
P08$tumor <- "primary"
P09$tumor <- "primary"
P10$tumor <- "primary"


P01$stage <- "localised"
P02$stage <- "localised"
P03$stage <- "localised"
P04$stage <- "localised"
P05$stage <- "localised"
P06$stage <- "localised"
P07$stage <- "localised"
P08$stage <- "localised"
P09$stage <- "localised"
P10$stage <- "localised"




M01$tumor <- "metastasis"
M02$tumor <- "metastasis"
M03$tumor <- "metastasis"
M04$tumor <- "metastasis"
M05$tumor <- "metastasis"


M01$stage <- "metastatic"
M02$stage <- "metastatic"
M03$stage <- "metastatic"
M04$stage <- "metastatic"
M05$stage <- "metastatic"


P01$dataset <- "one"
P02$dataset <- "one"
P03$dataset <- "one"
P04$dataset <- "one"
P05$dataset <- "one"
P06$dataset <- "one"
P07$dataset <- "one"
P08$dataset <- "one"
P09$dataset <- "one"
P10$dataset <- "one"

M01$dataset <- "one"
M02$dataset <- "one"
M03$dataset <- "one"
M04$dataset <- "one"
M05$dataset <- "one"




P01$location <- "pancreas"
P02$location <- "pancreas"
P03$location <- "pancreas"
P04$location <- "pancreas"
P05$location <- "pancreas"
P06$location <- "pancreas"
P07$location <- "pancreas"
P08$location <- "pancreas"
P09$location <- "pancreas"
P10$location <- "pancreas"

M01$location <- "liver"
M02$location <-  "liver"
M03$location <-  "liver"
M04$location <-  "liver"
M05$location <-  "liver"





WERBA_DIR <- "/rds/general/user/ts3225/ephemeral/R/Werba"

P11 <- Read10X(file.path(WERBA_DIR, "P4"))
P12 <- Read10X(file.path(WERBA_DIR, "P5"))
P13 <- Read10X(file.path(WERBA_DIR, "P7"))
P14 <- Read10X(file.path(WERBA_DIR, "P9"))
P15 <- Read10X(file.path(WERBA_DIR, "P13"))
P16 <- Read10X(file.path(WERBA_DIR, "P15"))
P17 <- Read10X(file.path(WERBA_DIR, "P19"))
P18 <- Read10X(file.path(WERBA_DIR, "P20"))
P19 <- Read10X(file.path(WERBA_DIR, "P22"))
P20 <- Read10X(file.path(WERBA_DIR, "P23"))
P21 <- Read10X(file.path(WERBA_DIR, "P26"))

M06 <- Read10X(file.path(WERBA_DIR, "P1"))
M07 <- Read10X(file.path(WERBA_DIR, "P2"))
M08 <- Read10X(file.path(WERBA_DIR, "P11"))
M09 <- Read10X(file.path(WERBA_DIR, "P16"))
M10 <- Read10X(file.path(WERBA_DIR, "P18"))
M11 <- Read10X(file.path(WERBA_DIR, "P21"))
M12 <- Read10X(file.path(WERBA_DIR, "P24"))
M13 <- Read10X(file.path(WERBA_DIR, "P25"))
M14 <- Read10X(file.path(WERBA_DIR, "P27"))


P11 <- CreateSeuratObject(counts = P11, project = "P11", min.cells = 3, min.features = 200)
P12 <- CreateSeuratObject(counts = P12, project = "P12", min.cells = 3, min.features = 200)
P13 <- CreateSeuratObject(counts = P13, project = "P13", min.cells = 3, min.features = 200)
P14 <- CreateSeuratObject(counts = P14, project = "P14", min.cells = 3, min.features = 200)
P15 <- CreateSeuratObject(counts = P15, project = "P15", min.cells = 3, min.features = 200)
P16 <- CreateSeuratObject(counts = P16, project = "P16", min.cells = 3, min.features = 200)
P17 <- CreateSeuratObject(counts = P17, project = "P17", min.cells = 3, min.features = 200)
P18 <- CreateSeuratObject(counts = P18, project = "P18", min.cells = 3, min.features = 200)
P19 <- CreateSeuratObject(counts = P19, project = "P19", min.cells = 3, min.features = 200)
P20 <- CreateSeuratObject(counts = P20, project = "P20", min.cells = 3, min.features = 200)
P21 <- CreateSeuratObject(counts = P21, project = "P21", min.cells = 3, min.features = 200)

M06 <- CreateSeuratObject(counts = M06, project = "M06", min.cells = 3, min.features = 200)
M07 <- CreateSeuratObject(counts = M07, project = "M07", min.cells = 3, min.features = 200)
M08 <- CreateSeuratObject(counts = M08, project = "M08", min.cells = 3, min.features = 200)
M09 <- CreateSeuratObject(counts = M09, project = "M09", min.cells = 3, min.features = 200)
M10 <- CreateSeuratObject(counts = M10, project = "M10", min.cells = 3, min.features = 200)
M11 <- CreateSeuratObject(counts = M11, project = "M11", min.cells = 3, min.features = 200)
M12 <- CreateSeuratObject(counts = M12, project = "M12", min.cells = 3, min.features = 200)
M13 <- CreateSeuratObject(counts = M13, project = "M13", min.cells = 3, min.features = 200)
M14 <- CreateSeuratObject(counts = M14, project = "M14", min.cells = 3, min.features = 200)


P11$tumor <- "primary"
P12$tumor <- "primary"
P13$tumor <- "primary"
P14$tumor <- "primary"
P15$tumor <- "primary"
P16$tumor <- "primary"
P17$tumor <- "primary"
P18$tumor <- "primary"
P19$tumor <- "primary"
P20$tumor <- "primary"
P21$tumor <- "primary"

M06$tumor <- "metastasis"
M07$tumor <- "metastasis"
M08$tumor <- "metastasis"
M09$tumor <- "metastasis"
M10$tumor <- "metastasis"
M11$tumor <- "metastasis"
M12$tumor <- "metastasis"
M13$tumor <- "metastasis"
M14$tumor <- "metastasis"



P11$stage <- "localised"
P12$stage <- "localised"
P13$stage <- "localised"
P14$stage <- "metastatic"
P15$stage <- "locally advanced"
P16$stage <- "localised"
P17$stage <- "localised"
P18$stage <- "metastatic"
P19$stage <- "locally advanced"
P20$stage <- "localised"
P21$stage <- "metastatic"

M06$stage <- "metastatic"
M07$stage <- "metastatic"
M08$stage <- "metastatic"
M09$stage <- "metastatic"
M10$stage <- "metastatic"
M11$stage <- "metastatic"
M12$stage <- "metastatic"
M13$stage <- "metastatic"
M14$stage <- "metastatic"



P11$dataset <- "two"
P12$dataset <- "two"
P13$dataset <- "two"
P14$dataset <- "two"
P15$dataset <- "two"
P16$dataset <- "two"
P17$dataset <- "two"
P18$dataset <- "two"
P19$dataset <- "two"
P20$dataset <- "two"
P21$dataset <- "two"

M06$dataset <- "two"
M07$dataset <- "two"
M08$dataset <- "two"
M09$dataset <- "two"
M10$dataset <- "two"
M11$dataset <- "two"
M12$dataset <- "two"
M13$dataset <- "two"
M14$dataset <- "two"



P11$location <- "pancreas"
P12$location  <- "pancreas"
P13$location  <- "pancreas"
P14$location  <- "pancreas"
P15$location  <- "pancreas"
P16$location  <- "pancreas"
P17$location  <- "pancreas"
P18$location  <- "pancreas"
P19$location  <- "pancreas"
P20$location  <- "pancreas"
P21$location  <- "pancreas"

M06$location <- "liver"
M07$location  <- "liver"
M08$location  <- "liver"
M09$location  <- "liver"
M10$location  <- "liver"
M11$location  <- "liver"
M12$location  <- "liver"
M13$location  <- "liver"
M14$location  <- "liver"

infile <- "/rds/general/user/ts3225/ephemeral/R/GSE263733_Raw_counts.txt.gz"


con <- gzfile(infile, "rt")
cell_names <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
close(con)

sample_id <- sub("_.*$", "", cell_names)

sample_counts <- sort(table(sample_id), decreasing = TRUE)
sample_counts



samples_to_use <- c( "PB2032-PAAD-T", "PB2032-PAAD-T2", "PB2151-PAAD-T", 
                     "PB2151-PAAD-T1", "PB2155-PAAD-T", "PB2191-PAAD-T", 
                     "PB2203-PAAD-T", "PB2218-PAAD-T", "PB2219-PAAD-T", 
                     "PB2256-PAAD-T", "PB2264-PAAD-T", "PB2265-PAAD-T",  
                     "PB2266-PAAD-T", "PB2268-PAAD-T", "PB2281-PAAD-T", 
                     "PB2286-PAAD-T", "PB2287-PAAD-T", "PB2311-PAAD-T",
                     "PB2341-PAAD-T", "PB2349-PAAD-T", "PB2366-PAAD-T", 
                     "PB2409-PAAD-T", "PB2410-PAAD-T",
                     "PB2155-Livermeta-T", "PB2191-Livermeta-T", "PB2264-Livermeta-T", 
                     "PB2311-Livermeta-T", "PB2349-Livermeta-T", "PB2409-Livermeta-T",
                     "PB2410-Livermeta-T")


setdiff(samples_to_use, names(sample_counts))
cells_by_sample <- split(seq_along(cell_names), sample_id)

make_seurat_for_sample <- function(sample_name) {
  if (!sample_name %in% names(cells_by_sample)) {
    stop("Sample not found: ", sample_name)
  }
  
  cols <- cells_by_sample[[sample_name]]
  
  if (length(cols) == 0) {
    stop("No cells found for sample: ", sample_name)
  }
  
  dt <- fread(
    cmd = paste("gzip -cd", shQuote(infile)),
    skip = 1,
    header = FALSE,
    select = c(1, cols + 1),
    data.table = FALSE
  )
  
  genes <- make.unique(as.character(dt[[1]]))
  counts <- dt[, -1, drop = FALSE]
  
  count_matrix <- as.matrix(counts)
  storage.mode(count_matrix) <- "numeric"
  
  mat <- Matrix(count_matrix, sparse = TRUE)
  rownames(mat) <- genes
  colnames(mat) <- cell_names[cols]
  
  mat <- mat[rownames(mat) != "" & !is.na(rownames(mat)), , drop = FALSE]
  
  CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
}


seurat_list <- lapply(samples_to_use, make_seurat_for_sample)
names(seurat_list) <- samples_to_use


P22 <- seurat_list[["PB2032-PAAD-T"]]
P23 <- seurat_list[["PB2032-PAAD-T2"]]
P24 <- seurat_list[["PB2151-PAAD-T"]]
P25 <- seurat_list[["PB2151-PAAD-T1"]]
P26 <- seurat_list[["PB2155-PAAD-T"]]
P27 <- seurat_list[["PB2191-PAAD-T"]]
P28 <- seurat_list[["PB2203-PAAD-T"]]
P29 <- seurat_list[["PB2218-PAAD-T"]]
P30 <- seurat_list[["PB2219-PAAD-T"]]
P31 <- seurat_list[["PB2256-PAAD-T"]]
P32 <- seurat_list[["PB2264-PAAD-T"]]
P33 <- seurat_list[["PB2265-PAAD-T"]]
P34 <- seurat_list[["PB2266-PAAD-T"]]
P35 <- seurat_list[["PB2268-PAAD-T"]]
P36 <- seurat_list[["PB2281-PAAD-T"]]
P37 <- seurat_list[["PB2286-PAAD-T"]]
P38 <- seurat_list[["PB2287-PAAD-T"]]
P39 <- seurat_list[["PB2311-PAAD-T"]]
P40 <- seurat_list[["PB2341-PAAD-T"]]
P41 <- seurat_list[["PB2349-PAAD-T"]]
P42 <- seurat_list[["PB2366-PAAD-T"]]
P43 <- seurat_list[["PB2409-PAAD-T"]]
P44 <- seurat_list[["PB2410-PAAD-T"]]


M15 <- seurat_list[["PB2155-Livermeta-T"]]
M16 <- seurat_list[["PB2191-Livermeta-T"]]
M17 <- seurat_list[["PB2264-Livermeta-T"]]
M18 <- seurat_list[["PB2311-Livermeta-T"]]
M19 <- seurat_list[["PB2349-Livermeta-T"]]
M20 <- seurat_list[["PB2409-Livermeta-T"]]
M21 <- seurat_list[["PB2410-Livermeta-T"]]



P22$tumor <- "primary"
P23$tumor <- "primary"
P24$tumor <- "primary"
P25$tumor <- "primary"
P26$tumor <- "primary"
P27$tumor <- "primary"
P28$tumor <- "primary"
P29$tumor <- "primary"
P30$tumor <- "primary"
P31$tumor <- "primary"
P32$tumor <- "primary"
P33$tumor <- "primary"
P34$tumor <- "primary"
P35$tumor <- "primary"
P36$tumor <- "primary"
P37$tumor <- "primary"
P38$tumor <- "primary"
P39$tumor <- "primary"
P40$tumor <- "primary"
P41$tumor <- "primary"
P42$tumor <- "primary"
P43$tumor <- "primary"
P44$tumor <- "primary"

M15$tumor <- "metastasis"
M16$tumor <- "metastasis"
M17$tumor <- "metastasis"
M18$tumor <- "metastasis"
M19$tumor <- "metastasis"
M20$tumor <- "metastasis"
M21$tumor <- "metastasis"



P22$stage <- "locally advanced"
P23$stage <- "locally advanced"
P24$stage <- "metastatic"
P25$stage <- "metastatic"
P26$stage <- "metastatic"
P27$stage <- "metastatic"
P28$stage <- "metastatic"
P29$stage <- "metastatic"
P30$stage <- "metastatic"
P31$stage <- "locally advanced"
P32$stage <- "metastatic"
P33$stage <- "metastatic"
P34$stage <- "locally advanced"
P35$stage <- "locally advanced"
P36$stage <- "metastatic"
P37$stage <- "metastatic"
P38$stage <- "locally advanced"
P39$stage <- "metastatic"
P40$stage <- "metastatic"
P41$stage <- "metastatic"
P42$stage <- "locally advanced"
P43$stage <- "metastatic"
P44$stage <- "metastatic"

M15$stage <- "metastatic"
M16$stage <- "metastatic"
M17$stage <- "metastatic"
M18$stage <- "metastatic"
M19$stage <- "metastatic"
M20$stage <- "metastatic"
M21$stage <- "metastatic"


P22$dataset <- "three"
P23$dataset <- "three"
P24$dataset <- "three"
P25$dataset <- "three"
P26$dataset <- "three"
P27$dataset <- "three"
P28$dataset <- "three"
P29$dataset <- "three"
P30$dataset <- "three"
P31$dataset <- "three"
P32$dataset <- "three"
P33$dataset <- "three"
P34$dataset <- "three"
P35$dataset <- "three"
P36$dataset <- "three"
P37$dataset <- "three"
P38$dataset <- "three"
P39$dataset <- "three"
P40$dataset <- "three"
P41$dataset <- "three"
P42$dataset <- "three"
P43$dataset <- "three"
P44$dataset <- "three"

M15$dataset <- "three"
M16$dataset <- "three"
M17$dataset <- "three"
M18$dataset <- "three"
M19$dataset <- "three"
M20$dataset <- "three"
M21$dataset <- "three"




P22$location <- "pancreas"
P23$location <- "pancreas"
P24$location <- "pancreas"
P25$location <- "pancreas"
P26$location <- "pancreas"
P27$location <- "pancreas"
P28$location <- "pancreas"
P29$location <- "pancreas"
P30$location <- "pancreas"
P31$location <- "pancreas"
P32$location <- "pancreas"
P33$location <- "pancreas"
P34$location <- "pancreas"
P35$location <- "pancreas"
P36$location <- "pancreas"
P37$location <- "pancreas"
P38$location <- "pancreas"
P39$location <- "pancreas"
P40$location <- "pancreas"
P41$location <- "pancreas"
P42$location <- "pancreas"
P43$location <- "pancreas"
P44$location <- "pancreas"

M15$location <- "liver"
M16$location <- "liver"
M17$location <- "liver"
M18$location <- "liver"
M19$location <- "liver"
M20$location <- "liver"
M21$location <- "liver"



base_dir <- "/rds/general/user/ts3225/ephemeral/R/CG"

P45 <- Read10X(file.path(base_dir, "42_V3", "filtered_feature_bc_matrix"))
P46 <- Read10X(file.path(base_dir, "45", "filtered_feature_bc_matrix"))
P47 <- Read10X(file.path(base_dir, "47", "filtered_feature_bc_matrix"))


M22 <- Read10X(file.path(base_dir, "43", "filtered_feature_bc_matrix"))
M23 <- Read10X(file.path(base_dir, "46", "filtered_feature_bc_matrix"))
M24 <- Read10X(file.path(base_dir, "48", "filtered_feature_bc_matrix"))
M25 <- Read10X(file.path(base_dir, "49", "filtered_feature_bc_matrix"))


P45 <- CreateSeuratObject(counts = P45, project = "P45", min.cells = 3, min.features = 200)
P46 <- CreateSeuratObject(counts = P46, project = "P46", min.cells = 3, min.features = 200)
P47 <- CreateSeuratObject(counts = P47, project = "P47", min.cells = 3, min.features = 200)

M22 <- CreateSeuratObject(counts = M22, project = "M22", min.cells = 3, min.features = 200)
M23 <- CreateSeuratObject(counts = M23, project = "M23", min.cells = 3, min.features = 200)
M24 <- CreateSeuratObject(counts = M24, project = "M24", min.cells = 3, min.features = 200)
M25 <- CreateSeuratObject(counts = M25, project = "M25", min.cells = 3, min.features = 200)


P45$tumor <- "primary"
P46$tumor <- "primary"
P47$tumor <- "primary"

M22$tumor <- "metastasis"
M23$tumor <- "metastasis"
M24$tumor <- "metastasis"
M25$tumor <- "metastasis"





P45$stage <- "metastatic"
P46$stage <- "metastatic"
P47$stage <- "metastatic"


M22$stage <- "metastatic"
M23$stage <- "metastatic"
M24$stage <- "metastatic"
M25$stage <- "metastatic"

P45$dataset <- "four"
P46$dataset <- "four"
P47$dataset <- "four"

M22$dataset <- "four"
M23$dataset <- "four"
M24$dataset <- "four"
M25$dataset <- "four"



P45$location <- "pancreas"
P46$location <- "pancreas"
P47$location <- "pancreas"

M22$location <- "liver"
M23$location <- "liver"
M24$location <- "liver"
M25$location <- "liver"




P01$orig.ident <- "P01"
P02$orig.ident <- "P02"
P03$orig.ident <- "P03"
P04$orig.ident <- "P04"
P05$orig.ident <- "P05"
P06$orig.ident <- "P06"
P07$orig.ident <- "P07"
P08$orig.ident <- "P08"
P09$orig.ident <- "P09"
P10$orig.ident <- "P10"
P11$orig.ident <- "P11"
P12$orig.ident <- "P12"
P13$orig.ident <- "P13"
P14$orig.ident <- "P14"
P15$orig.ident <- "P15"
P16$orig.ident <- "P16"
P17$orig.ident <- "P17"
P18$orig.ident <- "P18"
P19$orig.ident <- "P19"
P20$orig.ident <- "P20"
P21$orig.ident <- "P21"
P22$orig.ident <- "P22"
P23$orig.ident <- "P23"
P24$orig.ident <- "P24"
P25$orig.ident <- "P25"
P26$orig.ident <- "P26"
P27$orig.ident <- "P27"
P28$orig.ident <- "P28"
P29$orig.ident <- "P29"
P30$orig.ident <- "P30"
P31$orig.ident <- "P31"
P32$orig.ident <- "P32"
P33$orig.ident <- "P33"
P34$orig.ident <- "P34"
P35$orig.ident <- "P35"
P36$orig.ident <- "P36"
P37$orig.ident <- "P37"
P38$orig.ident <- "P38"
P39$orig.ident <- "P39"
P40$orig.ident <- "P40"
P41$orig.ident <- "P41"
P42$orig.ident <- "P42"
P43$orig.ident <- "P43"
P44$orig.ident <- "P44"
P45$orig.ident <- "P45"
P46$orig.ident <- "P46"
P47$orig.ident <- "P47"

M01$orig.ident <- "M01"
M02$orig.ident <- "M02"
M03$orig.ident <- "M03"
M04$orig.ident <- "M04"
M05$orig.ident <- "M05"
M06$orig.ident <- "M06"
M07$orig.ident <- "M07"
M08$orig.ident <- "M08"
M09$orig.ident <- "M09"
M10$orig.ident <- "M10"
M11$orig.ident <- "M11"
M12$orig.ident <- "M12"
M13$orig.ident <- "M13"
M14$orig.ident <- "M14"
M15$orig.ident <- "M15"
M16$orig.ident <- "M16"
M17$orig.ident <- "M17"
M18$orig.ident <- "M18"
M19$orig.ident <- "M19"
M20$orig.ident <- "M20"
M21$orig.ident <- "M21"
M22$orig.ident <- "M22"
M23$orig.ident <- "M23"
M24$orig.ident <- "M24"
M25$orig.ident <- "M25"



P01$patient <- "P01"
P02$patient <- "P02"
P03$patient <- "P03"
P04$patient <- "P04"
P05$patient <- "P05"
P06$patient <- "P06"
P07$patient <- "P07"
P08$patient <- "P08"
P09$patient <- "P09"
P10$patient <- "P10"
P11$patient <- "P11"
P12$patient <- "P12"
P13$patient <- "P13"
P14$patient <- "P14"
P15$patient <- "P15"
P16$patient <- "P16"
P17$patient <- "P17"
P18$patient <- "P18"
P19$patient <- "P19"
P20$patient <- "P20"
P21$patient <- "P21"


P22$patient <- "P22"
P23$patient <- "P22"
P24$patient <- "P23"
P25$patient <- "P23"
P26$patient <- "P24"
P27$patient <- "P25"
P28$patient <- "P26"
P29$patient <- "P27"
P30$patient <- "P28"
P31$patient <- "P29"
P32$patient <- "P30"
P33$patient <- "P31"
P34$patient <- "P32"
P35$patient <- "P33"
P36$patient <- "P34"
P37$patient <- "P35"
P38$patient <- "P36"
P39$patient <- "P37"
P40$patient <- "P38"

P41$patient <- "P39"
P42$patient <- "P40"
P43$patient <- "P41"
P44$patient <- "P42"
P45$patient <- "P43"
P46$patient <- "P44"
P47$patient <- "P45"

M01$patient <- "M01"
M02$patient <- "M02"
M03$patient <- "M03"
M04$patient <- "M04"
M05$patient <- "M05"
M06$patient <- "M06"
M07$patient <- "M07"
M08$patient <- "M08"
M09$patient <- "M09"
M10$patient <- "M10"
M11$patient <- "M11"
M12$patient <- "M12"
M13$patient <- "M13"
M14$patient <- "M14"
M15$patient <- "M15"
M16$patient <- "M16"
M17$patient <- "M17"
M18$patient <- "M18"
M19$patient <- "M19"
M20$patient <- "M20"
M21$patient <- "M21"
M22$patient <- "M22"
M23$patient <- "M23"
M24$patient <- "M24"
M25$patient <- "M25"




PDAC <- merge( P01,
               y = c(P02, P03, P04, P05, P06, P07, P08, P09, P10, P11,
                     P12, P13, P14, P15, P16, P17, P18, P19, P20, P21,
                     P22, P23, P24, P25, P26, P27, P28, P29, P30, P31,
                     P32, P33, P34, P35, P36, P37, P38, P39, P40, P41,
                     P42, P43, P44, P45, P46, P47, M01, M02, M03, M04, M05,
                     M06, M07, M08, M09, M10, M11, M12, M13, M14, M15,
                     M16, M17, M18, M19, M20, M21, M22, M23, M24, M25
               ),
               add.cell.ids = c ("P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10",
                                 "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20",
                                 "P21", "P22", "P23", "P24", "P25", "P26", "P27", "P28", "P29", "P30",
                                 "P31", "P32", "P33", "P34", "P35", "P36", "P37", "P38", "P39", "P40",
                                 "P41", "P42", "P43", "P44", "P45", "P46", "P47", "M01", "M02", "M03", "M04",
                                 "M05", "M06", "M07", "M08", "M09", "M10", "M11", "M12", "M13", "M14",
                                 "M15", "M16", "M17", "M18", "M19", "M20", "M21", "M22", "M23", "M24",
                                 "M25"
               ),
               project = "PDAC")

PDAC_joined <- JoinLayers(PDAC, assay = "RNA")


saveRDS(PDAC_joined, file = "/rds/general/user/ts3225/ephemeral/R/QC4/PDAC_joined.rds")

rm(
  list = c(
    sprintf("P%02d", 1:47),
    sprintf("M%02d", 1:25)))

sample_metadata <- unique(PDAC_joined@meta.data[, c("orig.ident", "tumor", "location", "stage", "dataset")])
rownames(sample_metadata) <- sample_metadata$orig.ident

PDAC_bulk <- AggregateExpression(
  object = PDAC_joined,
  assays = "RNA",
  group.by = "orig.ident",
  return.seurat = TRUE)


PDAC_bulk$orig.ident <- colnames(PDAC_bulk)
PDAC_bulk$tumor <- sample_metadata[PDAC_bulk$orig.ident, "tumor"]
PDAC_bulk$location <- sample_metadata[PDAC_bulk$orig.ident, "location"]
PDAC_bulk$stage <- sample_metadata[PDAC_bulk$orig.ident, "stage"]
PDAC_bulk$dataset <- sample_metadata[PDAC_bulk$orig.ident, "dataset"]


PDAC_bulk <- NormalizeData(PDAC_bulk)
PDAC_bulk <- FindVariableFeatures(PDAC_bulk, selection.method = "vst", nfeatures = 3000)
PDAC_bulk <- ScaleData(PDAC_bulk, features = VariableFeatures(PDAC_bulk))

set.seed(1234)
PDAC_bulk <- RunPCA(PDAC_bulk, features = VariableFeatures(PDAC_bulk), npcs = 30, seed.use = 1234)
pca_var <- PDAC_bulk[["pca"]]@stdev^2
pca_var <- pca_var / sum(pca_var) * 100

pca_df <- as.data.frame(Embeddings(PDAC_bulk, "pca")[, 1:2])
pca_df$orig.ident <- rownames(pca_df)
pca_df$tumor <- PDAC_bulk$tumor
pca_df$location <- PDAC_bulk$location
pca_df$stage <- PDAC_bulk$stage
pca_df$dataset <- PDAC_bulk$dataset

PCA1 <- ggplot(pca_df, aes(PC_1, PC_2, color = stage, shape = location, label = orig.ident)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5) +
  theme_classic() +
  labs(
    title = "PCA before QC",
    x = paste0("PC1: ", round(pca_var[1], 1), "% variance"),
    y = paste0("PC2: ", round(pca_var[2], 1), "% variance"),
    colour = "Disease Stage",
    shape = "Tumor Location"
  )


saveRDS(
  PCA1,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/PCA1.rds"
)
ggsave(file.path(out_dir, "PCA1.pdf"), plot = PCA1, width = 10, height = 10)

set.seed(1234)
sce <- as.SingleCellExperiment(PDAC_joined, assay = "RNA")
sce$orig.ident <- PDAC_joined$orig.ident

bp <- SerialParam(RNGseed = 1234)
sce <- scDblFinder(sce, samples = "orig.ident", BPPARAM = bp)
PDAC_joined$scDblFinder.score <- sce$scDblFinder.score
PDAC_joined$scDblFinder.class <- sce$scDblFinder.class

table(PDAC_joined$scDblFinder.class)

before <- table(PDAC_joined$orig.ident)

after <- table(PDAC_joined$orig.ident[PDAC_joined$scDblFinder.class == "singlet"])

after <- after[names(before)]
after[is.na(after)] <- 0


scDblFinder_numbers <- data.frame(
  sample = names(before),
  before = as.numeric(before),
  after = as.numeric(after[names(before)])
)


scDblFinder_numbers <- scDblFinder_numbers %>%
  mutate(
    removed = before - after,
    percent_removed = round(100 * removed / before, 2)
  )

write.csv(scDblFinder_numbers, file = "/rds/general/user/ts3225/ephemeral/R/QC4/scDblFinder_numbers.csv")

scDblFinder_long <- scDblFinder_numbers %>%
  select(sample, before, after) %>%
  pivot_longer(
    cols = c(before, after),
    names_to = "scdbl_step",
    values_to = "cells"
  )

write.csv(scDblFinder_long, file = "/rds/general/user/ts3225/ephemeral/R/QC4/scDblFinder_long.csv")

scDblFinder_long$scdbl_step <- factor(
  scDblFinder_long$scdbl_step,
  levels = c("before", "after"),
  labels = c("Before", "After")
)

scdbl1 <- ggplot(scDblFinder_long, aes(x = sample, y = cells, fill = scdbl_step)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(
    aes(label = cells),
    position = position_dodge(width = 0.9),
    vjust = -0.1,
    size = 2
  ) +   
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Cell counts before and after scDblFinder",
    x = "Samples",
    y = "Cell counts",
    fill = ""
  )


saveRDS(
  scdbl1,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/scdbl1.rds"
)
ggsave(file.path(out_dir, "scdbl1.pdf"), plot = scdbl1, width = 20, height = 10)

scdbl2 <- ggplot(scDblFinder_numbers, aes(x = sample, y = removed)) +
  geom_col(fill = "firebrick") +
  geom_text(
    aes(label = removed),
    vjust = -0.3,
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Cells removed by scDblFinder",
    x = "Samples",
    y = "Number of predicted doublets"
  )

saveRDS(
  scdbl2,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/scdbl2.rds"
)
ggsave(file.path(out_dir, "scdbl2.pdf"),  plot = scdbl2, width = 20, height = 10)



PDAC_no_doublets <- subset(PDAC_joined,
                           subset = scDblFinder.class == "singlet")

sample_metadata_scdbl <- unique(PDAC_no_doublets@meta.data[, c("orig.ident", "tumor", "location", "stage", "dataset")])
rownames(sample_metadata_scdbl) <- sample_metadata_scdbl$orig.ident

PDAC_scdbl_bulk <- AggregateExpression(
  object = PDAC_no_doublets,
  assays = "RNA",
  group.by = "orig.ident",
  return.seurat = TRUE
)

PDAC_scdbl_bulk$orig.ident <- colnames(PDAC_scdbl_bulk)
PDAC_scdbl_bulk$tumor <- sample_metadata_scdbl[PDAC_scdbl_bulk$orig.ident, "tumor"]
PDAC_scdbl_bulk$location <- sample_metadata_scdbl[PDAC_scdbl_bulk$orig.ident, "location"]
PDAC_scdbl_bulk$stage <- sample_metadata_scdbl[PDAC_scdbl_bulk$orig.ident, "stage"]
PDAC_scdbl_bulk$dataset <- sample_metadata_scdbl[PDAC_scdbl_bulk$orig.ident, "dataset"]

PDAC_scdbl_bulk <- NormalizeData(PDAC_scdbl_bulk)
PDAC_scdbl_bulk <- FindVariableFeatures(PDAC_scdbl_bulk, selection.method = "vst", nfeatures = 3000)
PDAC_scdbl_bulk <- ScaleData(PDAC_scdbl_bulk, features = VariableFeatures(PDAC_scdbl_bulk))

set.seed(1234)
PDAC_scdbl_bulk <- RunPCA(PDAC_scdbl_bulk, features = VariableFeatures(PDAC_scdbl_bulk), npcs = 30, seed.use = 1234)

pca_var_scdbl <- PDAC_scdbl_bulk[["pca"]]@stdev^2
pca_var_scdbl <- pca_var_scdbl / sum(pca_var_scdbl) * 100

pca_df_scdbl <- as.data.frame(Embeddings(PDAC_scdbl_bulk, "pca")[, 1:2])
pca_df_scdbl$orig.ident <- rownames(pca_df_scdbl)
pca_df_scdbl$tumor <- PDAC_scdbl_bulk$tumor
pca_df_scdbl$location <- PDAC_scdbl_bulk$location
pca_df_scdbl$stage <- PDAC_scdbl_bulk$stage
pca_df_scdbl$dataset <- PDAC_scdbl_bulk$dataset

PCA2 <- ggplot(pca_df_scdbl, aes(PC_1, PC_2, color = stage, shape = location, label = orig.ident)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5) +
  theme_classic() +
  labs(
    title = "PCA after scDblFinder",
    x = paste0("PC1: ", round(pca_var_scdbl[1], 1), "% variance"),
    y = paste0("PC2: ", round(pca_var_scdbl[2], 1), "% variance"),
    colour = "Disease Stage",
    shape = "Tumor Location"
  )

saveRDS(
  PCA2,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/PCA2.rds"
)
ggsave(file.path(out_dir, "PCA2.pdf"), plot = PCA2, width = 10, height = 10)


sce_singlets <- as.SingleCellExperiment(PDAC_no_doublets, assay = "RNA")
sce_singlets$orig.ident <- PDAC_no_doublets$orig.ident
set.seed(1234)
sce_singlets <- decontX(
  sce_singlets,
  assayName = "counts",
  batch = sce_singlets$orig.ident,
  seed = 1234
)

PDAC_no_doublets$decontX_contamination <- sce_singlets$decontX_contamination

PDAC_no_doublets[["decontXcounts"]] <- CreateAssay5Object(
  counts = round(decontXcounts(sce_singlets))
)

PDAC_decontX <- PDAC_no_doublets
DefaultAssay(PDAC_decontX) <- "decontXcounts"
raw <- GetAssayData(PDAC_decontX, assay = "RNA", layer = "counts")
dx  <- GetAssayData(PDAC_decontX, assay = "decontXcounts", layer = "counts")

common_cells <- intersect(colnames(raw), colnames(dx))
common_genes <- intersect(rownames(raw), rownames(dx))

raw <- raw[common_genes, common_cells, drop = FALSE]
dx  <- dx[common_genes, common_cells, drop = FALSE]

decontX_qc <- data.frame(
  RNA_counts = Matrix::colSums(raw),
  decontX_counts = Matrix::colSums(dx),
  RNA_features = Matrix::colSums(raw > 0),
  decontX_features = Matrix::colSums(dx > 0),
  row.names = common_cells
)

decontX_qc$counts_subtracted_by_decontX <- 
  decontX_qc$RNA_counts - decontX_qc$decontX_counts

decontX_qc$percentage_of_counts_subtracted_by_decontX <- 
  100 * decontX_qc$counts_subtracted_by_decontX / decontX_qc$RNA_counts

PDAC_decontX <- AddMetaData(PDAC_decontX, decontX_qc)


decontX_summary_table <- PDAC_decontX@meta.data %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    median_contamination = median(decontX_contamination, na.rm = TRUE),
    mean_contamination = mean(decontX_contamination, na.rm = TRUE),
    median_percent_of_counts_subtracted = median(percentage_of_counts_subtracted_by_decontX, na.rm = TRUE),
    mean_percent_of_counts_subtracted = mean(percentage_of_counts_subtracted_by_decontX, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

write.csv(decontX_summary_table, file = "/rds/general/user/ts3225/ephemeral/R/QC4/decontX_summary_table.csv")

decontX_summary_table_clean <- decontX_summary_table %>%
  dplyr::mutate(
    median_contamination = round(median_contamination, 3),
    mean_contamination = round(mean_contamination, 3),
    median_percent_of_counts_subtracted = round(median_percent_of_counts_subtracted, 2),
    mean_percent_of_counts_subtracted = round(mean_percent_of_counts_subtracted, 2)
  )


write.csv(decontX_summary_table_clean, file = "/rds/general/user/ts3225/ephemeral/R/QC4/PdecontX_summary_table_clean.csv")

decontX_summary_plot <- decontX_summary_table_clean %>%
  dplyr::select(
    orig.ident,
    median_contamination,
    median_percent_of_counts_subtracted ,
    n_cells
  ) %>%
  tidyr::pivot_longer(
    cols = c(median_contamination, median_percent_of_counts_subtracted),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      median_contamination = "Median estimated contamination",
      median_percent_of_counts_subtracted = "Median % of UMIs subtracted"
    ),
    plot_value = dplyr::if_else(
      metric == "Median estimated contamination",
      value * 100,
      value
    )
  )
write.csv(decontX_summary_plot, file = "/rds/general/user/ts3225/ephemeral/R/QC4/decontX_summary_plot.csv")

decontX_summary_plot$metric <- factor(
  decontX_summary_plot$metric,
  levels = c(
    "Median estimated contamination",
    "Median % of UMIs subtracted"
  )
)

decontX1 <- ggplot(decontX_summary_plot, aes(x = orig.ident, y = plot_value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(
    aes(label = paste0(round(plot_value, 1), "%")),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "DecontX correction summary per sample",
    x = "Samples",
    y = "Percentage of correction",
    fill = ""
  )

saveRDS(
  decontX1,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/decontX1.rds"
)
ggsave(file.path(out_dir, "decontX1.pdf"), plot = decontX1, width = 30, height = 10)


markers <- list(
  "T cells" = c("CD3D", "CD3E"),
  "NK cells" = c("GNLY", "NKG7"),
  "B cells" = c("CD79A", "CD79B", "MS4A1"),
  "Epithelial cells" = c("EPCAM", "KRT18", "KRT8"),
  "Fibroblasts" = c("COL1A1", "DCN"),
  "Endothelial cells" = c("PECAM1", "VWF"),
  "Macrophages" = c("CD68", "CD163", "LYZ"))


PDAC_decontX <- NormalizeData(PDAC_decontX, assay = "RNA")
PDAC_decontX <- NormalizeData(PDAC_decontX, assay = "decontXcounts")

p_before <- DotPlot(
  PDAC_decontX,
  features = markers,
  assay = "RNA",
  group.by = "orig.ident"
) + RotatedAxis() + ggtitle("Before DecontX") +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 14)
  )

p_after <- DotPlot(
  PDAC_decontX,
  features = markers,
  assay = "decontXcounts",
  group.by = "orig.ident"
) + RotatedAxis() + ggtitle("After DecontX: corrected counts")+
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 14)
  )


decontX3 <- p_before / p_after
saveRDS(
  decontX3,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/decontX3.rds"
)

ggsave(file.path(out_dir, "decontX3.pdf"), plot = decontX3, width = 30, height = 20)


sample_metadata_decontX <- unique(PDAC_decontX@meta.data[, c("orig.ident", "tumor", "location", "stage", "dataset")])
rownames(sample_metadata_decontX) <- sample_metadata_decontX$orig.ident

PDAC_decontX_bulk <- AggregateExpression(
  object = PDAC_decontX,
  assays = "decontXcounts",
  group.by = "orig.ident",
  return.seurat = TRUE
)

PDAC_decontX_bulk$orig.ident <- colnames(PDAC_decontX_bulk)
PDAC_decontX_bulk$tumor <- sample_metadata_decontX[PDAC_decontX_bulk$orig.ident, "tumor"]
PDAC_decontX_bulk$location <- sample_metadata_decontX[PDAC_decontX_bulk$orig.ident, "location"]
PDAC_decontX_bulk$stage <- sample_metadata_decontX[PDAC_decontX_bulk$orig.ident, "stage"]
PDAC_decontX_bulk$dataset <- sample_metadata_decontX[PDAC_decontX_bulk$orig.ident, "dataset"]

PDAC_decontX_bulk <- NormalizeData(PDAC_decontX_bulk)
PDAC_decontX_bulk <- FindVariableFeatures(PDAC_decontX_bulk, selection.method = "vst", nfeatures = 3000)
PDAC_decontX_bulk <- ScaleData(PDAC_decontX_bulk, features = VariableFeatures(PDAC_decontX_bulk))

set.seed(1234)
PDAC_decontX_bulk <- RunPCA(PDAC_decontX_bulk, features = VariableFeatures(PDAC_decontX_bulk), npcs = 30, seed.use = 1234)

pca_var_decontX <- PDAC_decontX_bulk[["pca"]]@stdev^2
pca_var_decontX <- pca_var_decontX / sum(pca_var_decontX) * 100

pca_df_decontX <- as.data.frame(Embeddings(PDAC_decontX_bulk, "pca")[, 1:2])
pca_df_decontX$orig.ident <- rownames(pca_df_decontX)
pca_df_decontX$tumor <- PDAC_decontX_bulk$tumor
pca_df_decontX$location <- PDAC_decontX_bulk$location
pca_df_decontX$stage <- PDAC_decontX_bulk$stage
pca_df_decontX$dataset <- PDAC_decontX_bulk$dataset

PCA3 <- ggplot(pca_df_decontX, aes(PC_1, PC_2, color = stage, shape = location, label = orig.ident)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5) +
  theme_classic() +
  labs(
    title = "PCA after scDblFinder + DecontX",
    x = paste0("PC1: ", round(pca_var_decontX[1], 1), "% variance"),
    y = paste0("PC2: ", round(pca_var_decontX[2], 1), "% variance"),
    colour = "Disease Stage",
    shape = "Tumor Location"
  )

saveRDS(
  PCA3,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/PCA3.rds"
)
ggsave(file.path(out_dir, "PCA3.pdf"), plot = PCA3, width = 10, height = 10)


DefaultAssay(PDAC_decontX) <- "RNA"

PDAC_decontX[["percent.mt"]] <- PercentageFeatureSet(PDAC_decontX, pattern = "^MT")
qc_violin <- VlnPlot(PDAC_decontX, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),  group.by = "orig.ident", ncol = 3, pt.size = 0, raster = FALSE) &
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10))

saveRDS(
  qc_violin,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/qc_violin.rds"
)
ggsave(file.path(out_dir, "qc_violin.pdf"), plot = qc_violin, width = 30, height = 20)


plot1 <- FeatureScatter(PDAC_decontX, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot2 <- FeatureScatter(PDAC_decontX, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
Combined_plot <- plot1 + plot2
saveRDS(Combined_plot, file.path(out_dir, "Combined_plot.rds"))
ggsave(file.path(out_dir, "Combined_plot.pdf"), plot = Combined_plot, width = 20, height = 20)

QC1 <- ggplot(PDAC_decontX@meta.data, aes(color=orig.ident, x= nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()
saveRDS(QC1 , file.path(out_dir, "QC1.rds"))
ggsave(file.path(out_dir, "QC1.pdf"), plot = QC1, width = 15, height = 10)


QC2 <- ggplot(PDAC_decontX@meta.data, aes(x=percent.mt,fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  scale_x_log10()+
  theme_classic() +
  geom_vline(xintercept = 25,color="red",linetype="dotted")
saveRDS(QC2 , file.path(out_dir, "QC2.rds"))
ggsave(file.path(out_dir, "QC2.pdf"), plot = QC2, width = 15, height = 10)


PDAC_clean <- subset(PDAC_decontX, subset = percent.mt < 25)
DefaultAssay(PDAC_clean) <- "decontXcounts"


sample_metadata_QC <- unique(PDAC_clean@meta.data[, c("orig.ident", "tumor", "location", "stage", "dataset")])
rownames(sample_metadata_QC) <- sample_metadata_QC$orig.ident

PDAC_QC <- AggregateExpression(
  object = PDAC_clean,
  assays = "decontXcounts",
  group.by = "orig.ident",
  return.seurat = TRUE
)

PDAC_QC$orig.ident <- colnames(PDAC_QC)
PDAC_QC$tumor <- sample_metadata_QC[PDAC_QC$orig.ident, "tumor"]
PDAC_QC$location <- sample_metadata_QC[PDAC_QC$orig.ident, "location"]
PDAC_QC$stage <- sample_metadata_QC[PDAC_QC$orig.ident, "stage"]
PDAC_QC$dataset <- sample_metadata_QC[PDAC_QC$orig.ident, "dataset"]

PDAC_QC <- NormalizeData(PDAC_QC)
PDAC_QC <- FindVariableFeatures(PDAC_QC, selection.method = "vst", nfeatures = 3000)
PDAC_QC <- ScaleData(PDAC_QC, features = VariableFeatures(PDAC_QC))

set.seed(1234)
PDAC_QC <- RunPCA(
  PDAC_QC,
  features = VariableFeatures(PDAC_QC),
  npcs = 30, seed.use = 1234
)

pca_var_QC <- PDAC_QC[["pca"]]@stdev^2
pca_var_QC <- pca_var_QC / sum(pca_var_QC) * 100

pca_df_QC <- as.data.frame(Embeddings(PDAC_QC, "pca")[, 1:2])
pca_df_QC$orig.ident <- rownames(pca_df_QC)
pca_df_QC$tumor <- PDAC_QC$tumor
pca_df_QC$location <- PDAC_QC$location
pca_df_QC$stage <- PDAC_QC$stage
pca_df_QC$dataset <- PDAC_QC$dataset

PCA4 <- ggplot(pca_df_QC, aes(PC_1, PC_2, color = stage, shape = location, label = orig.ident)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3.5) +
  theme_classic() +
  labs(
    title = "PCA after QC filtering",
    x = paste0("PC1: ", round(pca_var_QC[1], 1), "% variance"),
    y = paste0("PC2: ", round(pca_var_QC[2], 1), "% variance"),
    colour = "Disease Stage",
    shape = "Tumor Location"
  )

saveRDS(
  PCA4,
  file = "/rds/general/user/ts3225/ephemeral/R/QC4/PCA4.rds")

ggsave(file.path(out_dir, "PCA4.pdf"), plot = PCA4, width = 10, height = 10)

saveRDS(PDAC_clean, file = "/rds/general/user/ts3225/ephemeral/R/QC4/PDAC_clean.rds")

PDAC_primary <- subset(PDAC_clean, subset = tumor == "primary")
PDAC_met <- subset(PDAC_clean, subset = tumor == "metastasis")

DefaultAssay(PDAC_primary) <- "decontXcounts"
PDAC_primary <- JoinLayers(PDAC_primary, assay = "decontXcounts")


counts_mat <- LayerData(
  PDAC_primary,
  assay = "decontXcounts",
  layer = "counts"
)

PDAC_primary$nCount_decontX <- Matrix::colSums(counts_mat)
PDAC_primary$nFeature_decontX <- Matrix::colSums(counts_mat > 0)

PDAC_primary <- subset(
  PDAC_primary,
  subset = nCount_decontX > 0 & nFeature_decontX > 200
)

DefaultAssay(PDAC_primary) <- "decontXcounts"
PDAC_primary[["decontXcounts"]] <- split(
  PDAC_primary[["decontXcounts"]],
  f = PDAC_primary$orig.ident
)

PDAC_primary <- NormalizeData(PDAC_primary)
PDAC_primary <- FindVariableFeatures(PDAC_primary, selection.method = "vst", nfeatures = 3000)
top20primary <- head(VariableFeatures(PDAC_primary), 20)


primaryplot1 <- VariableFeaturePlot(PDAC_primary)
primaryplot2 <- LabelPoints(plot = primaryplot1, points = top20primary, repel = TRUE)
featuresprimary <- primaryplot1 + primaryplot2
saveRDS(featuresprimary, file.path(out_dir, "featuresprimary.rds"))
ggsave(file.path(out_dir, "featuresprimary.pdf"), plot = featuresprimary, width = 20, height = 20)


PDAC_primary <- ScaleData(PDAC_primary, features = VariableFeatures(PDAC_primary))

set.seed(1234)
PDAC_primary <- RunPCA(PDAC_primary, features = VariableFeatures(PDAC_primary), npcs = 50, seed.use = 1234)
pca.emb <- Embeddings(PDAC_primary, "pca")[, 1:30]

sum(!is.finite(as.matrix(pca.emb)))
sum(duplicated(as.data.frame(pca.emb)))

ElbowPlotprimary <- ElbowPlot(PDAC_primary, ndims = 50)
ggsave(file.path(out_dir, "ElbowPlotprimary.pdf"), plot = ElbowPlotprimary, width = 10, height = 10)



PDAC_primary <- FindNeighbors(PDAC_primary, reduction = "pca", dims = 1:30)
set.seed(1234)
PDAC_primary <- FindClusters(PDAC_primary, resolution = 0.5, random.seed = 1234)

set.seed(1234)
PDAC_primary <- RunUMAP(
  PDAC_primary,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap.unintegrated",
  seed.use = 1234
)

DimPlot1 <- DimPlot(PDAC_primary, reduction = "umap.unintegrated", group.by = "dataset", raster = FALSE)
saveRDS(DimPlot1, file.path(out_dir, "DimPlot1.rds"))
ggsave(file.path(out_dir, "DimPlot1.pdf"), plot = DimPlot1, width = 20, height = 10)
DimPlot2 <- DimPlot(PDAC_primary, reduction = "umap.unintegrated", group.by = "orig.ident", raster = FALSE)
saveRDS(DimPlot2, file.path(out_dir, "DimPlot2.rds"))
ggsave(file.path(out_dir, "DimPlot2.pdf"), plot = DimPlot2, width = 20, height = 10)



set.seed(1234)
PDAC_primary <- IntegrateLayers(
  object = PDAC_primary,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

PDAC_primary <- FindNeighbors(PDAC_primary, reduction = "harmony", dims = 1:30)
set.seed(1234)
PDAC_primary <- FindClusters(PDAC_primary, resolution = 0.5, random.seed = 1234)

set.seed(1234)
PDAC_primary <- RunUMAP(
  PDAC_primary,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap_harmony",
  seed.use = 1234
)

DimPlot3 <- DimPlot(PDAC_primary, reduction = "umap_harmony", group.by = "dataset", raster = FALSE)
DimPlot4 <- DimPlot(PDAC_primary, reduction = "umap_harmony", group.by = "orig.ident", raster = FALSE)
saveRDS(DimPlot3, file.path(out_dir, "DimPlot3.rds"))
saveRDS(DimPlot4, file.path(out_dir, "DimPlot4.rds"))
ggsave(file.path(out_dir, "DimPlot3.pdf"), plot = DimPlot3, width = 20, height = 10)
ggsave(file.path(out_dir, "DimPlot4.pdf"), plot = DimPlot4, width = 20, height = 10)


PDAC_primary <- JoinLayers(PDAC_primary, assay = "decontXcounts")
DefaultAssay(PDAC_primary) <- "decontXcounts"
primaryseuratmarkers0.5 <- FindAllMarkers(PDAC_primary,
                                          only.pos = TRUE,
                                          min.pct = 0.25,
                                          logfc.threshold = 0.25)

write.csv(primaryseuratmarkers0.5, file = "/rds/general/user/ts3225/ephemeral/R/QC4/primaryseuratmarkers0.5.csv")



primarytop50_markers <- primaryseuratmarkers0.5 %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50)

write.csv(primarytop50_markers, file = "/rds/general/user/ts3225/ephemeral/R/QC4/primarytop50_markers.csv")

saveRDS(PDAC_primary, file = "/rds/general/user/ts3225/ephemeral/R/QC4/PDAC_primary.rds")


DefaultAssay(PDAC_met) <- "decontXcounts"
PDAC_met <- JoinLayers(PDAC_met, assay = "decontXcounts")


counts_mat2 <- LayerData(
  PDAC_met,
  assay = "decontXcounts",
  layer = "counts"
)

PDAC_met$nCount_decontX <- Matrix::colSums(counts_mat2)
PDAC_met$nFeature_decontX <- Matrix::colSums(counts_mat2 > 0)

PDAC_met <- subset(
  PDAC_met,
  subset = nCount_decontX > 0 & nFeature_decontX > 200
)

DefaultAssay(PDAC_met) <- "decontXcounts"
PDAC_met[["decontXcounts"]] <- split(
  PDAC_met[["decontXcounts"]],
  f = PDAC_met$orig.ident
)
PDAC_met <- NormalizeData(PDAC_met)
PDAC_met <- FindVariableFeatures(PDAC_met, selection.method = "vst", nfeatures = 3000)
top20met <- head(VariableFeatures(PDAC_met), 20)


metplot1 <- VariableFeaturePlot(PDAC_met)
metplot2 <- LabelPoints(plot = metplot1, points = top20met, repel = TRUE)
featuresmet <- metplot1 + metplot2
saveRDS(featuresmet, file.path(out_dir, "featuresmet.rds"))
ggsave(file.path(out_dir, "featuresmet.pdf"), plot = featuresmet, width = 20, height = 20)

metallgenes <- rownames(PDAC_met)
PDAC_met <- ScaleData(PDAC_met, features = VariableFeatures(PDAC_met))

set.seed(1234)
PDAC_met <- RunPCA(PDAC_met, features = VariableFeatures(PDAC_met), npcs = 50, seed.use = 1234)
pca.emb2 <- Embeddings(PDAC_met, "pca")[, 1:30]

sum(!is.finite(as.matrix(pca.emb2)))
sum(duplicated(as.data.frame(pca.emb2)))

ElbowPlotmet <- ElbowPlot(PDAC_met, ndims = 50)
ggsave(file.path(out_dir, "ElbowPlotmet.pdf"), plot = ElbowPlotmet, width = 10, height = 10)



PDAC_met <- FindNeighbors(PDAC_met, reduction = "pca", dims = 1:30)
set.seed(1234)
PDAC_met <- FindClusters(PDAC_met, resolution = 0.5, random.seed = 1234)

set.seed(1234)
PDAC_met <- RunUMAP(
  PDAC_met,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap.unintegrated",
  seed.use = 1234
)

DimPlot5 <- DimPlot(PDAC_met, reduction = "umap.unintegrated", group.by = "dataset", raster = FALSE)
saveRDS(DimPlot5, file.path(out_dir, "DimPlot5.rds"))
ggsave(file.path(out_dir, "DimPlot5.pdf"), plot = DimPlot5, width = 20, height = 10)
DimPlot6 <- DimPlot(PDAC_met, reduction = "umap.unintegrated", group.by = "orig.ident", raster = FALSE)
saveRDS(DimPlot6, file.path(out_dir, "DimPlot6.rds"))
ggsave(file.path(out_dir, "DimPlot6.pdf"), plot = DimPlot6, width = 20, height = 10)




set.seed(1234)
PDAC_met <- IntegrateLayers(
  object = PDAC_met,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

PDAC_met <- FindNeighbors(PDAC_met, reduction = "harmony", dims = 1:30)
set.seed(1234)
PDAC_met <- FindClusters(PDAC_met, resolution = 0.5, random.seed = 1234)

set.seed(1234)
PDAC_met <- RunUMAP(
  PDAC_met,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap_harmony",
  seed.use = 1234
)

DimPlot7 <- DimPlot(PDAC_met, reduction = "umap_harmony", group.by = "dataset", raster = FALSE)
DimPlot8 <- DimPlot(PDAC_met, reduction = "umap_harmony", group.by = "orig.ident", raster = FALSE)
saveRDS(DimPlot7, file.path(out_dir, "DimPlot7.rds"))
saveRDS(DimPlot8, file.path(out_dir, "DimPlot8.rds"))
ggsave(file.path(out_dir, "DimPlot7.pdf"), plot = DimPlot7, width = 20, height = 10)
ggsave(file.path(out_dir, "DimPlot8.pdf"), plot = DimPlot8, width = 20, height = 10)


PDAC_met <- JoinLayers(PDAC_met, assay = "decontXcounts")
DefaultAssay(PDAC_met) <- "decontXcounts"
metseuratmarkers0.5 <- FindAllMarkers(PDAC_met,
                                      only.pos = TRUE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25)

write.csv(metseuratmarkers0.5, file = "/rds/general/user/ts3225/ephemeral/R/QC4/metseuratmarkers0.5.csv")



mettop50_markers <- metseuratmarkers0.5 %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50)

write.csv(mettop50_markers, file = "/rds/general/user/ts3225/ephemeral/R/QC4/mettop50_markers.csv")

saveRDS(PDAC_met, file = "/rds/general/user/ts3225/ephemeral/R/QC4/PDAC_met.rds")

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

packages_used <- c(
  "Seurat", "SeuratObject", "tidyverse", "ggplot2",
  "SingleCellExperiment", "celda", "decontX", "scDblFinder",
  "scater", "scran", "patchwork", "dplyr", "harmony",
  "SingleR", "celldex", "scRNAseq", "broom", "DESeq2",
  "ggrepel", "data.table", "Matrix", "BiocParallel"
)

package_versions <- data.frame(
  package = packages_used,
  version = sapply(packages_used, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      as.character(packageVersion(pkg))
    } else {
      NA_character_
    }
  }),
  row.names = NULL
)

write.csv(
  package_versions,
  file.path(out_dir, "package_versions.csv"),
  row.names = FALSE
)



