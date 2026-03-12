#!/usr/bin/env Rscript
# Make reference coverage grid: databases (rows) x cancer types (columns)
# Run from package root: Rscript inst/figures/make_reference_coverage.R

devtools::load_all(quiet = TRUE)

precog <- PhenoMapR:::get_data("precog")
tcga <- PhenoMapR:::get_data("tcga")
pediatric <- PhenoMapR:::get_data("pediatric")
ici <- PhenoMapR:::get_data("ici")

# PRECOG full name -> canonical code (TCGA or PRECOG-only)
precog_to_tcga <- c(
  Adrenocortical = "ACC", Appendix = "Appendix", Bladder = "BLCA",
  Brain_Astrocytoma = "LGG", Brain_Glioblastoma = "GBM", Brain_Glioma = "LGG",
  Brain_Medulloblastoma = "MB", Brain_Meningioma = "Meningioma",
  Brain_Neuroblastoma = "NB", Breast = "BRCA", Breast_Metastasis = "BRCA_Met",
  Cervical = "CESC", Colon = "COAD", Colon_Metastasis = "COAD_Met",
  Colon_Rectal = "READ", Desmoid = "Desmoid", Gastric = "STAD",
  Germ_cell_tumors = "TGCT", Head_and_neck = "HNSC",
  Head_and_neck_Larynx = "HNSC", Head_and_neck_Oral_SCC = "HNSC",
  Hematopoietic_AML = "LAML", Hematopoietic_B.ALL = "ALL",
  Hematopoietic_Burkitt_lymphoma = "Burkitt", Hematopoietic_CLL = "CLL",
  Hematopoietic_DLBCL = "DLBC", Hematopoietic_DLBCL_CHOP.treated = "DLBC",
  Hematopoietic_DLBCL_RCHOP.treated = "DLBC", Hematopoietic_FL = "FL",
  Hematopoietic_Mantle_cell_lymphoma = "MCL",
  Hematopoietic_Multiple_myeloma = "MM",
  Hematopoietic_Peripheral_T.cell_Lymphoma = "PTCL",
  Kidney = "KIRC", Liver = "LIHC", Liver_Primary = "LIHC",
  Lung_ADENO = "LUAD", Lung_LCC = "Lung_LCC", Lung_SCC = "LUSC",
  Lung_SCLC = "SCLC", Melanoma = "SKCM", Melanoma_Metastasis = "SKCM_Met",
  Mesothelioma = "MESO", Ovarian = "OV", Ovarian_Epithelial = "OV",
  Pancreatic = "PAAD", Pancreatic_Metastasis = "PAAD_Met", Prostate = "PRAD",
  Sarcoma_Ewing_sarcoma = "ESFT", Sarcoma_Liposarcoma = "SARC",
  Sarcoma_Osteosarcoma = "OSTEO", Sarcoma_Uterine = "UCS"
)

# ICI prefix -> TCGA code
ici_to_tcga <- c(PDAC = "PAAD", CRC = "COAD", ESCC = "ESCA", OCCC = "OV")

# Pediatric: keep native codes
tcga_types <- setdiff(colnames(tcga), c("PRECOG_metaZ", "TCGA_metaZ"))
pediatric_types <- setdiff(colnames(pediatric), "meta Z (all tumors)")
ici_prefixes <- unique(sub("_.*", "", colnames(ici)))
ici_prefixes <- setdiff(ici_prefixes, "Mixed")

# Canonical columns: TCGA + PRECOG-only + Pediatric + ICI
all_canon <- unique(c(
  tcga_types,
  precog_to_tcga,
  pediatric_types,
  ifelse(ici_prefixes %in% names(ici_to_tcga), ici_to_tcga[ici_prefixes], ici_prefixes)
))
canonical <- sort(all_canon)

# Build presence matrix
db_names <- c("Adult PRECOG", "Pediatric PRECOG", "TCGA", "ICI PRECOG")
mat <- matrix(0, nrow = 4, ncol = length(canonical))
rownames(mat) <- db_names
colnames(mat) <- canonical

# Adult PRECOG
for (pt in colnames(precog)) {
  canon <- if (pt %in% names(precog_to_tcga)) precog_to_tcga[pt] else pt
  if (canon %in% canonical) mat["Adult PRECOG", canon] <- 1
}

# Pediatric PRECOG
for (pt in pediatric_types) {
  if (pt %in% canonical) mat["Pediatric PRECOG", pt] <- 1
}

# TCGA
for (t in tcga_types) {
  if (t %in% canonical) mat["TCGA", t] <- 1
}

# ICI PRECOG
for (ip in ici_prefixes) {
  canon <- if (ip %in% names(ici_to_tcga)) ici_to_tcga[ip] else ip
  if (canon %in% canonical) mat["ICI PRECOG", canon] <- 1
}

# Drop empty columns
mat <- mat[, colSums(mat) > 0, drop = FALSE]

# Order columns: TCGA, then PRECOG metastasis, Pediatric, ICI
tcga_cols <- intersect(tcga_types, colnames(mat))
met_cols <- c("BRCA_Met", "PAAD_Met", "COAD_Met", "SKCM_Met")
met_cols <- intersect(met_cols, colnames(mat))
ped_cols <- intersect(pediatric_types, colnames(mat))
ici_cols <- setdiff(colnames(mat), c(tcga_cols, met_cols, ped_cols))
col_order <- c(sort(tcga_cols), sort(met_cols), sort(ped_cols), sort(ici_cols))
mat <- mat[, col_order, drop = FALSE]

# Color scheme per database
color_scheme <- c(
  "PRECOG" = "#44536a",
  "TCGA" = "#9f9f9f",
  "ICI" = "#95c277",
  "Pediatric" = "#d58457"
)
# Map row names to scheme keys
db_to_key <- c(
  "Adult PRECOG" = "PRECOG",
  "Pediatric PRECOG" = "Pediatric",
  "TCGA" = "TCGA",
  "ICI PRECOG" = "ICI"
)

# Long format for ggplot2: database, cancer_type, available (0/1)
library(ggplot2)
df <- expand.grid(
  database = db_names,
  cancer_type = colnames(mat),
  stringsAsFactors = FALSE
)
df$available <- as.vector(mat[cbind(
  match(df$database, db_names),
  match(df$cancer_type, colnames(mat))
)])

# Order y-axis (database) by number of cancers available (most first)
count_per_db <- setNames(rowSums(mat), db_names)
db_order <- names(sort(count_per_db, decreasing = TRUE))
df$database <- factor(df$database, levels = db_order)

# Fill: factor with explicit levels so colors map correctly
fill_levels <- c("Not available", "PRECOG", "TCGA", "ICI", "Pediatric")
db_key_vec <- setNames(
  c("PRECOG", "Pediatric", "TCGA", "ICI"),
  c("Adult PRECOG", "Pediatric PRECOG", "TCGA", "ICI PRECOG")
)
df$fill <- ifelse(df$available == 1L, db_key_vec[as.character(df$database)], "Not available")
df$fill <- factor(df$fill, levels = fill_levels)

fill_colors <- c(
  "Not available" = "#F7F7F7",
  "PRECOG" = unname(color_scheme["PRECOG"]),
  "TCGA" = unname(color_scheme["TCGA"]),
  "ICI" = unname(color_scheme["ICI"]),
  "Pediatric" = unname(color_scheme["Pediatric"])
)

# ggplot2: minimal theme, square tiles, large font, minimal margins
out <- file.path("inst", "figures", "reference_coverage.png")
n_cancer <- ncol(mat)
n_db <- 4L
p <- ggplot(df, aes(x = .data$cancer_type, y = .data$database, fill = .data$fill)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(
    values = fill_colors,
    name = NULL,
    breaks = fill_levels,
    drop = FALSE
  ) +
  scale_x_discrete(expand = c(0, 0), guide = guide_axis(angle = 90)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(
    x = NULL,
    y = NULL,
    title = "Reference coverage: cancer types available for scoring"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "plain", size = 16),
    legend.text = element_text(size = 12),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    plot.margin = margin(0, 0, 0, 0),
    axis.ticks.length = unit(0, "pt"),
    axis.ticks = element_blank()
  )

# Save with a tightly sized PNG device to avoid extra whitespace from aspect constraints.
# We size the device to the tile grid plus room for axes/title/legend.
tile_px <- 22L
panel_w <- n_cancer * tile_px
panel_h <- n_db * tile_px

pad_left <- 220L   # y-axis labels
pad_right <- 30L
pad_top <- 120L    # legend + title
pad_bottom <- 260L # rotated x labels

png(
  filename = out,
  width = pad_left + panel_w + pad_right,
  height = pad_top + panel_h + pad_bottom,
  res = 150,
  bg = "white"
)
print(p)
dev.off()
message("Saved: ", out)
