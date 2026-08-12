#!/usr/bin/env Rscript
# Build a compact, bundled copy of the cell-type-resolved prognosis-related
# gene (crPRG) lists from Tang et al., Cancer Cell 2026 (ctPANDA), used by the
# methodology ctPANDA orthogonal validation
# (manuscript/benchmarks/methodology/test_ctpanda_enrichment.R).
#
# Source (open access supplemental tables, Cell Press / Elsevier CDN):
#   Tang R, Jin L, Chen C, et al. Mapping cell type-resolved transcriptomic
#   profiles to patient survival in pancreatic cancer. Cancer Cell (2026).
#   https://doi.org/10.1016/j.ccell.2026.05.012
#     Table S2A  -> mmc3.xlsx  (cell type-resolved PRGs; Worse=SS/Better=LS/NA)
#     Table S4B  -> mmc5.xlsx  (conserved LTS/STS DEGs with HR and log2FC)
#
# We redistribute only the small derived gene lists (gene symbol + cell type +
# prognostic direction, plus HR/log2FC for the 830 conserved features) so the
# methodology harness runs deterministically without network access.
#
# Usage:
#   Rscript tools/make_ctpanda_prgs.R            # downloads mmc files to /tmp
#   Rscript tools/make_ctpanda_prgs.R mmc3.xlsx mmc5.xlsx   # use local files

suppressPackageStartupMessages({
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Install 'readxl' to build the ctPANDA data.", call. = FALSE)
  }
})

args <- commandArgs(trailingOnly = TRUE)
pii <- "S1535610826002576"
base_url <- sprintf("https://ars.els-cdn.com/content/image/1-s2.0-%s-", pii)

fetch <- function(mmc, local) {
  if (!is.null(local) && file.exists(local)) return(local)
  dest <- file.path(tempdir(), mmc)
  if (!file.exists(dest)) {
    url <- paste0(base_url, mmc)
    message("Downloading ", url)
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

s2_file <- fetch("mmc3.xlsx", if (length(args) >= 1) args[[1]] else NULL)
s4_file <- fetch("mmc5.xlsx", if (length(args) >= 2) args[[2]] else NULL)

## ---- Table S2A: gene x cell-type matrix of Worse/Better/NA ----------------
s2a <- suppressMessages(readxl::read_excel(s2_file, sheet = "Table S2A", skip = 2))
s2a <- as.data.frame(s2a, stringsAsFactors = FALSE)
gene_col <- names(s2a)[1]
ct_cols <- setdiff(names(s2a), gene_col)
s2a_long <- do.call(rbind, lapply(ct_cols, function(ct) {
  v <- s2a[[ct]]
  keep <- !is.na(v) & v %in% c("Worse", "Better")
  if (!any(keep)) return(NULL)
  data.frame(
    gene = as.character(s2a[[gene_col]][keep]),
    cell_type = ct,
    direction = ifelse(v[keep] == "Worse", "SS", "LS"),
    stringsAsFactors = FALSE
  )
}))
rownames(s2a_long) <- NULL

## ---- Table S4B: conserved LTS/STS DEGs with HR / log2FC -------------------
s4b <- suppressMessages(readxl::read_excel(s4_file, sheet = "Table S4B", skip = 2))
s4b <- as.data.frame(s4b, stringsAsFactors = FALSE)
names(s4b) <- c("feature", "cell_type", "gene", "prognosis", "HR", "adj_p", "log2FC")[seq_along(names(s4b))]
s4b$direction <- ifelse(s4b$prognosis == "Worse", "SS", "LS")
s4b$HR <- suppressWarnings(as.numeric(s4b$HR))
s4b$adj_p <- suppressWarnings(as.numeric(s4b$adj_p))
s4b$log2FC <- suppressWarnings(as.numeric(s4b$log2FC))

obj <- list(
  s2a = s2a_long,
  s4b = s4b,
  source = list(
    citation = paste0(
      "Tang R, Jin L, Chen C, et al. Mapping cell type-resolved transcriptomic ",
      "profiles to patient survival in pancreatic cancer. Cancer Cell (2026)."
    ),
    doi = "10.1016/j.ccell.2026.05.012",
    tables = c("Table S2A (mmc3.xlsx)", "Table S4B (mmc5.xlsx)"),
    built = as.character(Sys.Date())
  )
)

outdir <- "inst/extdata/ctpanda"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- file.path(outdir, "ctpanda_prgs.rds")
saveRDS(obj, out, compress = "xz")

message("Wrote ", out, " (", round(file.info(out)$size / 1024, 1), " KB)")
message("S2A cell-type-gene calls: ", nrow(s2a_long),
        " | Ductal_PDAC SS=", sum(s2a_long$cell_type == "Ductal_PDAC" & s2a_long$direction == "SS"),
        " LS=", sum(s2a_long$cell_type == "Ductal_PDAC" & s2a_long$direction == "LS"))
message("S4B conserved features: ", nrow(s4b))
