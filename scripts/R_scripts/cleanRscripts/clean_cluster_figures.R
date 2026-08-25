library(Seurat)
library(tidyverse)
library(ggplot2)
library(grid)

#figure for composition#####making figures
PDAC_primary_filtered <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/PDAC_primary_filtered.rds")
PDAC_met_filtered <- readRDS(file = "/Users/thirisantracy/Desktop/thesis/PDAC_met_filtered.rds")

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


PDAC_primary_filtered$celltype_annotation <- factor(
  PDAC_primary_filtered$celltype_annotation,
  levels = all_celltypes
)

PDAC_met_filtered$celltype_annotation <- factor(
  PDAC_met_filtered$celltype_annotation,
  levels = all_celltypes
)

table(PDAC_primary_filtered$celltype_annotation)
table(PDAC_met_filtered$celltype_annotation)

primary_composition <- PDAC_primary_filtered@meta.data %>%
  mutate(
    patient = as.character(patient),
    celltype_annotation = as.character(celltype_annotation)
  ) %>%
  dplyr::count(
    dataset_GSE,
    patient,
    subject_ID,
    celltype_annotation,
    name = "n_cells"
  ) %>%
  group_by(dataset_GSE, patient, subject_ID) %>%
  mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells
  ) %>%
  ungroup()

met_composition <- PDAC_met_filtered@meta.data %>%
  mutate(
    patient = as.character(patient),
    celltype_annotation = as.character(celltype_annotation)
  ) %>%
  dplyr::count(
    dataset_GSE,
    patient,
    subject_ID,
    celltype_annotation,
    name = "n_cells"
  ) %>%
  group_by(dataset_GSE, patient, subject_ID) %>%
  mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells
  ) %>%
  ungroup()


primary_composition$celltype_annotation <- factor(
  primary_composition$celltype_annotation,
  levels = all_celltypes
)

met_composition$celltype_annotation <- factor(
  met_composition$celltype_annotation,
  levels = all_celltypes
)

########################################################
# ORDER PRIMARY SAMPLES NUMERICALLY
########################################################

primary_composition <- primary_composition %>%
  dplyr::mutate(
    patient = as.character(patient),
    sample_number = readr::parse_number(patient)
  )

primary_dataset_order <- primary_composition %>%
  dplyr::group_by(dataset_GSE) %>%
  dplyr::summarise(
    first_sample = min(sample_number),
    .groups = "drop"
  ) %>%
  dplyr::arrange(first_sample) %>%
  dplyr::pull(dataset_GSE)

primary_composition <- primary_composition %>%
  dplyr::arrange(sample_number) %>%
  dplyr::mutate(
    patient = factor(
      patient,
      levels = unique(patient)
    ),
    dataset_GSE = factor(
      dataset_GSE,
      levels = primary_dataset_order
    )
  )


########################################################
# ORDER METASTATIC SAMPLES NUMERICALLY
########################################################

met_composition <- met_composition %>%
  dplyr::mutate(
    patient = as.character(patient),
    sample_number = readr::parse_number(patient)
  )

met_dataset_order <- met_composition %>%
  dplyr::group_by(dataset_GSE) %>%
  dplyr::summarise(
    first_sample = min(sample_number),
    .groups = "drop"
  ) %>%
  dplyr::arrange(first_sample) %>%
  dplyr::pull(dataset_GSE)

met_composition <- met_composition %>%
  dplyr::arrange(sample_number) %>%
  dplyr::mutate(
    patient = factor(
      patient,
      levels = unique(patient)
    ),
    dataset_GSE = factor(
      dataset_GSE,
      levels = met_dataset_order
    )
  )


levels(primary_composition$patient)
levels(primary_composition$dataset_GSE)

levels(met_composition$patient)
levels(met_composition$dataset_GSE)

primary_composition_plot <- ggplot(
  primary_composition,
  aes(
    x = patient,
    y = proportion * 100,
    fill = celltype_annotation
  )
) +
  
  geom_col(
    width = 0.85
  ) +
  
  facet_grid(
    ~ dataset_GSE,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_fill_manual(
    values = celltype_colours,
    drop = TRUE
  ) +
  
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = "Sample",
    y = "Composition of the cell types (%)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 11,
      face = "plain",
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      size = 11,
      face = "plain",
      margin = margin(r = 10)
    ),
    
    strip.background = element_blank(),
    
    strip.text.x = element_text(
      size = 8,
      face = "bold",
      margin = margin(
        b = 7,
        l = 4,
        r = 4
      )
    ),
    
    panel.spacing.x = unit(
      1.2,
      "lines"
    ),
    
    legend.position = "right",
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.key.height = unit(
      0.38,
      "cm"
    ),
    
    legend.key.width = unit(
      0.38,
      "cm"
    ),
    
    plot.margin = margin(
      8, 8, 8, 8
    )
  )

primary_composition_plot

ggsave(
  "/Users/thirisantracy/Desktop/thesis/images/primary_composition_final.pdf",
  primary_composition_plot,
  width = 16,
  height = 6,
  units = "in",
  bg = "white"
)

met_composition_plot <- ggplot(
  met_composition,
  aes(
    x = patient,
    y = proportion * 100,
    fill = celltype_annotation
  )
) +
  
  geom_col(
    width = 0.85
  ) +
  
  facet_grid(
    ~ dataset_GSE,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_fill_manual(
    values = celltype_colours,
    drop = TRUE
  ) +
  
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = "Sample",
    y = "Composition of the cell types (%)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 11,
      face = "plain",
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      size = 11,
      face = "plain",
      margin = margin(r = 10)
    ),
    
    strip.background = element_blank(),
    
    strip.text.x = element_text(
      size = 8,
      face = "bold",
      margin = margin(
        b = 7,
        l = 4,
        r = 4
      )
    ),
    
    panel.spacing.x = unit(
      1.2,
      "lines"
    ),
    
    legend.position = "right",
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.key.height = unit(
      0.38,
      "cm"
    ),
    
    legend.key.width = unit(
      0.38,
      "cm"
    ),
    
    plot.margin = margin(
      8, 8, 8, 8
    )
  )

ggsave(
  "/Users/thirisantracy/Desktop/thesis/images/met_composition_final.pdf",
  met_composition_plot,
  width = 11,
  height = 6,
  units = "in",
  bg = "white"
)
##########composition for supplementary figure

primary_composition <- PDAC_primary_filtered@meta.data %>%
  dplyr::mutate(
    dataset_GSE = as.character(dataset_GSE),
    patient = as.character(patient),
    subject_ID = as.character(subject_ID),
    celltype_annotation = as.character(celltype_annotation)
  ) %>%
  dplyr::filter(
    !is.na(celltype_annotation)
  ) %>%
  dplyr::count(
    dataset_GSE,
    patient,
    subject_ID,
    celltype_annotation,
    name = "n_cells"
  ) %>%
  dplyr::group_by(
    dataset_GSE,
    patient,
    subject_ID
  ) %>%
  dplyr::mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells
  ) %>%
  dplyr::ungroup()

met_composition <- PDAC_met_filtered@meta.data %>%
  dplyr::mutate(
    dataset_GSE = as.character(dataset_GSE),
    patient = as.character(patient),
    subject_ID = as.character(subject_ID),
    celltype_annotation = as.character(celltype_annotation)
  ) %>%
  dplyr::filter(
    !is.na(celltype_annotation)
  ) %>%
  dplyr::count(
    dataset_GSE,
    patient,
    subject_ID,
    celltype_annotation,
    name = "n_cells"
  ) %>%
  dplyr::group_by(
    dataset_GSE,
    patient,
    subject_ID
  ) %>%
  dplyr::mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells
  ) %>%
  dplyr::ungroup()


primary_samples <- primary_composition %>%
  dplyr::distinct(
    dataset_GSE,
    patient,
    subject_ID
  )

met_samples <- met_composition %>%
  dplyr::distinct(
    dataset_GSE,
    patient,
    subject_ID
  )

nrow(primary_samples)
nrow(met_samples)

primary_sample_info <- primary_composition %>%
  dplyr::distinct(
    dataset_GSE,
    patient,
    subject_ID,
    total_cells
  )

primary_composition_complete <- tidyr::crossing(
  primary_sample_info,
  celltype_annotation = primary_types
) %>%
  dplyr::left_join(
    primary_composition %>%
      dplyr::select(
        dataset_GSE,
        patient,
        subject_ID,
        celltype_annotation,
        n_cells,
        proportion
      ),
    by = c(
      "dataset_GSE",
      "patient",
      "subject_ID",
      "celltype_annotation"
    )
  ) %>%
  dplyr::mutate(
    n_cells = tidyr::replace_na(
      n_cells,
      0L
    ),
    proportion = tidyr::replace_na(
      proportion,
      0
    ),
    tumor = "Primary"
  )

met_sample_info <- met_composition %>%
  dplyr::distinct(
    dataset_GSE,
    patient,
    subject_ID,
    total_cells
  )

met_composition_complete <- tidyr::crossing(
  met_sample_info,
  celltype_annotation = met_types
) %>%
  dplyr::left_join(
    met_composition %>%
      dplyr::select(
        dataset_GSE,
        patient,
        subject_ID,
        celltype_annotation,
        n_cells,
        proportion
      ),
    by = c(
      "dataset_GSE",
      "patient",
      "subject_ID",
      "celltype_annotation"
    )
  ) %>%
  dplyr::mutate(
    n_cells = tidyr::replace_na(
      n_cells,
      0L
    ),
    proportion = tidyr::replace_na(
      proportion,
      0
    ),
    tumor = "Metastasis"
  )


primary_composition_complete %>%
  dplyr::count(
    celltype_annotation
  )

met_composition_complete %>%
  dplyr::count(
    celltype_annotation)

primary_composition_summary <- primary_composition_complete %>%
  dplyr::mutate(
    percentage = proportion * 100
  ) %>%
  dplyr::group_by(
    celltype_annotation
  ) %>%
  dplyr::summarise(
    
    samples_detected = sum(
      n_cells > 0
    ),
    
    total_samples = dplyr::n(),
    
    detection_percent =
      100 * samples_detected / total_samples,
    
    median_percent = median(
      percentage,
      na.rm = TRUE
    ),
    
    Q1_percent = quantile(
      percentage,
      0.25,
      na.rm = TRUE
    ),
    
    Q3_percent = quantile(
      percentage,
      0.75,
      na.rm = TRUE
    ),
    
    minimum_percent = min(
      percentage,
      na.rm = TRUE
    ),
    
    maximum_percent = max(
      percentage,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

primary_composition_summary_formatted <-
  primary_composition_summary %>%
  dplyr::mutate(
    
    `Samples detected, n/N (%)` = paste0(
      samples_detected,
      "/",
      total_samples,
      " (",
      round(detection_percent, 1),
      "%)"
    ),
    
    `Median (%)` = round(
      median_percent,
      2
    ),
    
    `IQR (%)` = paste0(
      round(Q1_percent, 2),
      "–",
      round(Q3_percent, 2)
    ),
    
    `Min–max (%)` = paste0(
      round(minimum_percent, 2),
      "–",
      round(maximum_percent, 2)
    )
    
  ) %>%
  dplyr::select(
    `Cell type` = celltype_annotation,
    `Samples detected, n/N (%)`,
    `Median (%)`,
    `IQR (%)`,
    `Min–max (%)`
  )

View(primary_composition_summary_formatted)

met_composition_summary <- met_composition_complete %>%
  dplyr::mutate(
    percentage = proportion * 100
  ) %>%
  dplyr::group_by(
    celltype_annotation
  ) %>%
  dplyr::summarise(
    
    samples_detected = sum(
      n_cells > 0
    ),
    
    total_samples = dplyr::n(),
    
    detection_percent =
      100 * samples_detected / total_samples,
    
    median_percent = median(
      percentage,
      na.rm = TRUE
    ),
    
    Q1_percent = quantile(
      percentage,
      0.25,
      na.rm = TRUE
    ),
    
    Q3_percent = quantile(
      percentage,
      0.75,
      na.rm = TRUE
    ),
    
    minimum_percent = min(
      percentage,
      na.rm = TRUE
    ),
    
    maximum_percent = max(
      percentage,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

met_composition_summary_formatted <-
  met_composition_summary %>%
  dplyr::mutate(
    
    `Samples detected, n/N (%)` = paste0(
      samples_detected,
      "/",
      total_samples,
      " (",
      round(detection_percent, 1),
      "%)"
    ),
    
    `Median (%)` = round(
      median_percent,
      2
    ),
    
    `IQR (%)` = paste0(
      round(Q1_percent, 2),
      "–",
      round(Q3_percent, 2)
    ),
    
    `Min–max (%)` = paste0(
      round(minimum_percent, 2),
      "–",
      round(maximum_percent, 2)
    )
    
  ) %>%
  dplyr::select(
    `Cell type` = celltype_annotation,
    `Samples detected, n/N (%)`,
    `Median (%)`,
    `IQR (%)`,
    `Min–max (%)`
  )

View(met_composition_summary_formatted)


selected_celltypes <- c(
  "CAFs",
  "Tumour epithelial cells",
  "TAMs"
)

composition_selected <- dplyr::bind_rows(
  
  primary_composition_complete %>%
    dplyr::filter(
      celltype_annotation %in%
        selected_celltypes
    ),
  
  met_composition_complete %>%
    dplyr::filter(
      celltype_annotation %in%
        selected_celltypes
    )
)

composition_selected %>%
  dplyr::count(
    celltype_annotation,
    tumor
  )

selected_colours <- c(
  "CAFs" = "#667BC6",
  "Tumour epithelial cells" = "#4C78A8",
  "TAMs" = "#D95F4C"
)

selected_composition_plot <- ggplot(
  composition_selected,
  aes(
    x = tumor,
    y = proportion * 100,
    fill = celltype_annotation
  )
) +
  
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.75
  ) +
  
  geom_jitter(
    width = 0.12,
    size = 1.5,
    alpha = 0.65,
    colour = "grey45"
  ) +
  
  facet_wrap(
    ~ celltype_annotation,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = selected_colours,
    breaks = selected_celltypes
  ) +
  
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    x = "Tumour",
    y = "Composition of cell types (%)",
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    
    axis.text.x = element_text(
      size = 9
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    axis.title.x = element_text(
      size = 11,
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      size = 11,
      margin = margin(r = 8)
    ),
    
    panel.spacing.x = unit(
      0.6,
      "lines"
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
    )
  )

selected_composition_plot



ggsave(
  "/Users/thirisantracy/Desktop/thesis/images/selected_celltype_composition_final.pdf",
  selected_composition_plot,
  width = 10,
  height = 5.5,
  units = "in",
  bg = "white"
)
