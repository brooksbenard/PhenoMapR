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

# Plot
out <- file.path("inst", "figures", "reference_coverage.png")
png(out, width = 1400, height = 340, res = 120, bg = "white")
par(mar = c(10, 6, 3, 2))
image(1:ncol(mat), 1:4, t(mat[nrow(mat):1, , drop = FALSE]),
  col = c("#F7F7F7", "#2166AC"), axes = FALSE, xlab = "", ylab = "")
axis(1, at = 1:ncol(mat), labels = colnames(mat), las = 2, cex.axis = 0.75)
axis(2, at = 1:4, labels = rev(db_names), las = 1, cex.axis = 1)
title(main = "Reference coverage: cancer types available for scoring")
legend("topright", legend = c("Available", "Not available"),
  fill = c("#2166AC", "#F7F7F7"), bty = "n", cex = 0.9)
dev.off()

message("Saved: ", out)
