#!/usr/bin/env Rscript
#
# Incremental TCGA PRECOG bulk survival stratification:
# - download one TCGA TPM gz at a time
# - score via PhenoMapR (built-in PRECOG meta-z signatures)
# - stratify by median PhenoMapR score (High vs Low) *within each sample stratum*
# - compute Cox PH HR + CI and log-rank p-value per TCGA type
# - append sample-level results to a single consolidated table
# - delete the local TPM file after each cancer type is processed
#
# Expected local inputs:
# - Liu et al. outcomes (Cell 2018): vignettes/Liu_Cell_2018_TCGA_Outcomes.xlsx
#   Sheet "TCGA-CDR": OS/OS.time for most cancers; PFI/PFI.time for DLBC, TGCT, READ.
# - Optional: data/tcga/clinical_PANCAN_patient_with_followup.tsv (limits which TCGA types run)
#
# Outputs:
# - results/tcga_precog_forest_results.tsv — stratum-level Cox/log-rank (merged with forest-plot fields)
# - results/tcga_precog_forest_plot_data.tsv — exact rows/columns used for the forest figure (HR, CI, cor, etc.)
# - results/tcga_precog_forest_panel.tsv — sample-level scores + survival (minimal, for tertile script)
# - results/tcga_precog_forest_samples_annotated.tsv — per-sample PhenoMapR scores + Liu TCGA-CDR + PANCAN clinical
# - results/tcga_outcomes_methods_note.txt — why Liu vs PANCAN OS can differ
# - results/tcga_precog_forest.png — forest + patient n + PRECOG gene n
# Optional flags: --precog_full, --tcga_full (unshrunk meta-z for PRECOG–TCGA r in the plot)

suppressPackageStartupMessages({
  library(data.table)
  library(PhenoMapR)
  library(survival)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) == 0) return(default)
  if (isTRUE(default == TRUE)) return(TRUE)
  if (i == length(args)) return(TRUE)
  args[[i + 1]]
}

z_score_cutoff <- as.numeric(get_flag("--z_score_cutoff", "2"))
reference <- get_flag("--reference", "precog")

drive_folder_id <- get_flag("--drive_folder_id", "1zGT4k1y1KA9R78p6OMp9nmpQiDKupzJm")
drive_html_url <- paste0(
  "https://drive.google.com/drive/folders/",
  drive_folder_id,
  "?usp=drive_link"
)

tpm_dir <- get_flag("--tpm_dir", "data/tcga")
clinical_file <- get_flag("--clinical_file", file.path(tpm_dir, "clinical_PANCAN_patient_with_followup.tsv"))
liu_outcomes_xlsx <- get_flag("--liu_outcomes", "vignettes/Liu_Cell_2018_TCGA_Outcomes.xlsx")
out_dir <- get_flag("--out_dir", "results")
# Unshrunk meta-z tables for PRECOG–TCGA Pearson r (forest point sizes); optional, repo root by default
precog_full_path <- get_flag("--precog_full", "PRECOG_V2_Cancer_Type_Meta_Zscores_Final.rds")
tcga_full_path <- get_flag("--tcga_full", "TCGA_metaz.csv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

reference_mat <- PhenoMapR:::get_data(reference)

if (!file.exists(liu_outcomes_xlsx)) {
  stop("Missing Liu outcomes file: ", liu_outcomes_xlsx)
}

# DLBC, TGCT, READ: PFI time/event; all other TCGA types: OS time/event (Liu Cell 2018 TCGA-CDR).
USE_PFI_CANCERS <- c("DLBC", "TGCT", "READ")

load_liu_cdr <- function(path) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Install readxl to read Liu outcomes: install.packages(\"readxl\")")
  }
  df <- readxl::read_excel(path, sheet = "TCGA-CDR")
  data.table::as.data.table(df)
}

cat("Loading Liu outcomes:", liu_outcomes_xlsx, "\n")
liu_dt <- load_liu_cdr(liu_outcomes_xlsx)

clin <- NULL
if (file.exists(clinical_file)) {
  cat("Loading clinical (optional filter):", clinical_file, "\n")
  clin <- fread(clinical_file, showProgress = FALSE)
} else {
  cat("Clinical file not found (optional); using all Liu cancer types that map to PRECOG.\n")
}

# Adult PRECOG full-name -> TCGA code mapping (from inst/figures/make_reference_coverage.R)
precog_to_tcga <- c(
  Adrenocortical = "ACC", Appendix = "Appendix", Bladder = "BLCA", Brain_Astrocytoma = "LGG",
  Brain_Glioblastoma = "GBM", Brain_Glioma = "LGG", Brain_Medulloblastoma = "MB",
  Brain_Meningioma = "Meningioma", Brain_Neuroblastoma = "NB", Breast = "BRCA",
  Breast_Metastasis = "BRCA_Met", Cervical = "CESC", Colon = "COAD", Colon_Metastasis = "COAD_Met",
  Colon_Rectal = "READ", Desmoid = "Desmoid", Gastric = "STAD", Germ_cell_tumors = "TGCT",
  Head_and_neck = "HNSC", Head_and_neck_Larynx = "HNSC", Head_and_neck_Oral_SCC = "HNSC",
  Hematopoietic_AML = "LAML", Hematopoietic_B.ALL = "ALL",
  Hematopoietic_Burkitt_lymphoma = "Burkitt", Hematopoietic_CLL = "CLL",
  Hematopoietic_DLBCL = "DLBC", Hematopoietic_DLBCL_CHOP.treated = "DLBC",
  Hematopoietic_DLBCL_RCHOP.treated = "DLBC", Hematopoietic_FL = "FL",
  Hematopoietic_Mantle_cell_lymphoma = "MCL", Hematopoietic_Multiple_myeloma = "MM",
  Hematopoietic_Peripheral_T.cell_Lymphoma = "PTCL", Kidney = "KIRC", Liver = "LIHC",
  Liver_Primary = "LIHC", Lung_ADENO = "LUAD", Lung_LCC = "Lung_LCC", Lung_SCC = "LUSC",
  Lung_SCLC = "SCLC", Melanoma = "SKCM", Melanoma_Metastasis = "SKCM_Met", Mesothelioma = "MESO",
  Ovarian = "OV", Ovarian_Epithelial = "OV", Pancreatic = "PAAD",
  Pancreatic_Metastasis = "PAAD_Met", Prostate = "PRAD", Sarcoma_Ewing_sarcoma = "ESFT",
  Sarcoma_Liposarcoma = "SARC", Sarcoma_Osteosarcoma = "OSTEO", Sarcoma_Uterine = "UCS"
)

concordant_precog_label <- function(tcga_code) {
  labels <- names(precog_to_tcga)[precog_to_tcga == tcga_code]
  if (length(labels) == 0) return(NA_character_)

  # Prefer non-metastatic, then non-_Primary variants.
  non_met <- labels[!grepl("Metastasis", labels, fixed = TRUE)]
  labels <- if (length(non_met) > 0) non_met else labels

  non_primary <- labels[!grepl("_Primary", labels, fixed = TRUE)]
  if (length(non_primary) > 0) labels <- non_primary

  labels[[1]]
}

metastatic_precog_label <- function(tcga_code) {
  # Return the Adult PRECOG metastatic label (e.g. Breast_Metastasis) for a TCGA
  # code when available; else NA.
  canon_met <- paste0(tcga_code, "_Met")
  labels <- names(precog_to_tcga)[precog_to_tcga == canon_met]
  if (length(labels) == 0) return(NA_character_)
  # Prefer *_Metastasis labels if multiple
  met <- labels[grepl("Metastasis", labels, fixed = TRUE)]
  if (length(met) > 0) return(met[[1]])
  labels[[1]]
}

extract_tcga_code_from_filename <- function(path) {
  # e.g. TCGA_ACC_tpm.fullIDs.remapped.tsv.gz -> ACC
  sub("^TCGA_([^_]+)_tpm.*$", "\\1", basename(path))
}

sample_to_patient <- function(sample_ids) {
  parts <- strsplit(sample_ids, ".", fixed = TRUE)
  vapply(parts, function(x) {
    if (length(x) >= 3) paste(x[1], x[2], x[3], sep = "-") else NA_character_
  }, character(1))
}

sample_to_type_code <- function(sample_ids) {
  # TCGA sample type is embedded in the 4th dot-separated field, e.g.
  # "TCGA.OR.A5JG.01A.11R.A29S.07" -> "01" (Primary Tumor)
  parts <- strsplit(sample_ids, ".", fixed = TRUE)
  vapply(parts, function(x) {
    if (length(x) >= 4) substr(x[4], 1, 2) else NA_character_
  }, character(1))
}

sample_stratum <- function(sample_type) {
  # Map TCGA sample_type codes to analysis strata.
  # We keep normals separate from tumor samples on the forest plot.
  # 01/02/03/... are treated as "Tumor" except 06 metastatic and 11 normal.
  ifelse(sample_type == "11", "Normal",
    ifelse(sample_type == "06", "Metastatic", "Tumor")
  )
}

build_survival_liu <- function(liu_sub) {
  if (nrow(liu_sub) == 0L) {
    return(data.table(
      patient = character(),
      os_time = numeric(),
      os_event = integer()
    ))
  }
  tcga_code <- unique(as.character(liu_sub$type))
  if (length(tcga_code) != 1L) stop("Expected exactly one cancer type in liu_sub.")
  use_pfi <- tcga_code %in% USE_PFI_CANCERS
  if (use_pfi) {
    tt <- suppressWarnings(as.numeric(liu_sub$PFI.time))
    ev <- suppressWarnings(as.integer(round(liu_sub$PFI)))
  } else {
    tt <- suppressWarnings(as.numeric(liu_sub$OS.time))
    ev <- suppressWarnings(as.integer(round(liu_sub$OS)))
  }
  dt <- data.table(
    patient = liu_sub$bcr_patient_barcode,
    os_time = tt,
    os_event = ev
  )
  dt <- dt[is.finite(os_time) & !is.na(os_event)]
  dt <- dt[order(os_time), .SD[1], by = patient]
  dt
}

#' Per-sample panel + Liu TCGA-CDR columns (liu_*) + PANCAN clinical (pancan_*) + derived PANCAN OS for comparison.
enrich_panel_samples <- function(panel, liu_dt, clin) {
  p <- data.table::copy(panel)
  p[, survival_endpoint := ifelse(tcga_code %in% USE_PFI_CANCERS, "PFI", "OS")]
  p[, survival_time_used_days := os_time]
  p[, survival_event_used := os_event]

  l <- data.table::copy(liu_dt)
  id_l <- c("bcr_patient_barcode", "type")
  if (!all(id_l %in% names(l))) {
    stop("Liu TCGA-CDR sheet must include bcr_patient_barcode and type.")
  }
  rest_l <- setdiff(names(l), id_l)
  if (length(rest_l)) {
    new_l <- paste0("liu_", rest_l)
    if (any(duplicated(new_l))) new_l <- make.unique(new_l, sep = "_")
    data.table::setnames(l, rest_l, new_l)
  }
  data.table::setnames(l, "bcr_patient_barcode", "patient")
  data.table::setnames(l, "type", "tcga_code")
  p <- merge(p, l, by = c("patient", "tcga_code"), all.x = TRUE)

  if (!is.null(clin) && nrow(clin) > 0L) {
    c <- data.table::copy(clin)
    id_c <- c("bcr_patient_barcode", "acronym")
    if (!all(id_c %in% names(c))) {
      stop("PANCAN clinical must include bcr_patient_barcode and acronym.")
    }
    rest_c <- setdiff(names(c), id_c)
    if (length(rest_c)) {
      new_c <- paste0("pancan_", rest_c)
      if (any(duplicated(new_c))) new_c <- make.unique(new_c, sep = "_")
      data.table::setnames(c, rest_c, new_c)
    }
    data.table::setnames(c, "bcr_patient_barcode", "patient")
    data.table::setnames(c, "acronym", "tcga_code")
    p <- merge(p, c, by = c("patient", "tcga_code"), all.x = TRUE)
  }

  if ("pancan_vital_status" %in% names(p)) {
    p[, pancan_os_time := ifelse(
      pancan_vital_status == "Dead",
      suppressWarnings(as.numeric(pancan_days_to_death)),
      ifelse(
        pancan_vital_status == "Alive",
        suppressWarnings(as.numeric(pancan_days_to_last_followup)),
        NA_real_
      )
    )]
    p[, pancan_os_event := NA_integer_]
    p[pancan_vital_status == "Dead", pancan_os_event := 1L]
    p[pancan_vital_status == "Alive", pancan_os_event := 0L]
  } else {
    p[, pancan_os_time := NA_real_]
    p[, pancan_os_event := NA_integer_]
  }

  p[, liu_minus_pancan_os_time_days := ifelse(
    survival_endpoint == "OS" & is.finite(survival_time_used_days) & is.finite(pancan_os_time),
    survival_time_used_days - pancan_os_time,
    NA_real_
  )]
  p[, outcomes_table := "Liu_Cell_2018_TCGA_CDR"]
  p
}

write_outcomes_note <- function(path) {
  txt <- c(
    "TCGA PhenoMapR forest analysis — outcomes sources and hazard ratios",
    "=================================================================",
    "",
    "Primary survival endpoint for Cox / log-rank:",
    "  Liu et al. (Cell 2018) supplementary table: Liu_Cell_2018_TCGA_Outcomes.xlsx, sheet \"TCGA-CDR\".",
    "",
    "Endpoint rules (same as script):",
    "  - Most cancers: OS and OS.time (days).",
    "  - DLBC, TGCT, READ: PFI and PFI.time.",
    "",
    "Why hazard ratios can differ from an analysis using only clinical_PANCAN_patient_with_followup.tsv:",
    "  - Liu TCGA-CDR is a curated, pan-cancer harmonized endpoint set with publication-specific censoring,",
    "    redaction (see column Redaction), and inclusion rules.",
    "  - PANCAN uses GDC fields vital_status, days_to_death, and days_to_last_followup from the clinical",
    "    snapshot; follow-up dates and death/censor definitions need not match Liu row-for-row.",
    "  - PhenoMapR scores are merged to samples; survival uses patient-level Liu times after merge with",
    "    expression, so the analyzed risk set follows Liu, not raw PANCAN OS alone.",
    "",
    "The file tcga_precog_forest_samples_annotated.tsv includes:",
    "  - PhenoMapR score columns (score, precog_label_used, stratum, ...)",
    "  - Liu columns prefixed with liu_ (including liu_OS, liu_OS.time, liu_PFI, liu_PFI.time, ...)",
    "  - PANCAN columns prefixed with pancan_ plus pancan_os_time / pancan_os_event derived like the legacy",
    "    build_os() helper (Dead -> days_to_death, Alive -> days_to_last_followup).",
    "  - liu_minus_pancan_os_time_days compares Liu vs PANCAN-derived OS time when survival_endpoint == OS.",
    ""
  )
  writeLines(txt, path)
}

probe_first_bytes_gzip <- function(file_id) {
  # Check gzip magic bytes (0x1f 0x8b) from Google Drive "uc?export=download"
  # Read only first 2 bytes for speed.
  u <- paste0("https://drive.google.com/uc?export=download&id=", file_id)
  con <- try(url(u, open = "rb"), silent = TRUE)
  if (inherits(con, "try-error")) return(FALSE)
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  b <- try(readBin(con, what = "raw", n = 2), silent = TRUE)
  if (inherits(b, "try-error") || length(b) < 2) return(FALSE)
  identical(b[1], as.raw(0x1f)) && identical(b[2], as.raw(0x8b))
}

gzip_original_name <- function(file_id) {
  # Parse gzip original filename (FNAME) from the gzip header when present.
  # Properly handles optional fields per RFC 1952: FEXTRA/FNAME/FCOMMENT/FHCRC.
  u <- paste0("https://drive.google.com/uc?export=download&id=", file_id)
  con <- try(url(u, open = "rb"), silent = TRUE)
  if (inherits(con, "try-error")) return(NA_character_)
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  b <- try(readBin(con, what = "raw", n = 1024), silent = TRUE)
  if (inherits(b, "try-error") || length(b) < 10) return(NA_character_)
  if (!(identical(b[1], as.raw(0x1f)) && identical(b[2], as.raw(0x8b)))) return(NA_character_)

  flg <- as.integer(b[4])
  pos <- 11L  # 1-indexed position after the fixed 10-byte header

  # FEXTRA: 2-byte little-endian length, then that many bytes
  if (bitwAnd(flg, 0x04) != 0) {
    if (length(b) < pos + 1L) return(NA_character_)
    xlen <- as.integer(b[pos]) + 256L * as.integer(b[pos + 1L])
    pos <- pos + 2L + xlen
    if (pos > length(b)) return(NA_character_)
  }

  # FNAME: NUL-terminated string
  if (bitwAnd(flg, 0x08) == 0) return(NA_character_)
  if (pos > length(b)) return(NA_character_)
  nul_rel <- which(b[pos:length(b)] == as.raw(0x00))
  if (length(nul_rel) == 0) return(NA_character_)
  end <- pos + nul_rel[1] - 2L
  if (end < pos) return(NA_character_)
  name_raw <- b[pos:end]
  paste(rawToChar(name_raw, multiple = TRUE), collapse = "")
}

gzip_original_name_from_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 12) return(NA_character_)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  b <- readBin(con, what = "raw", n = 1024)
  if (length(b) < 10) return(NA_character_)
  if (!(identical(b[1], as.raw(0x1f)) && identical(b[2], as.raw(0x8b)))) return(NA_character_)
  flg <- as.integer(b[4])
  pos <- 11L
  if (bitwAnd(flg, 0x04) != 0) {
    if (length(b) < pos + 1L) return(NA_character_)
    xlen <- as.integer(b[pos]) + 256L * as.integer(b[pos + 1L])
    pos <- pos + 2L + xlen
    if (pos > length(b)) return(NA_character_)
  }
  if (bitwAnd(flg, 0x08) == 0) return(NA_character_)
  if (pos > length(b)) return(NA_character_)
  nul_rel <- which(b[pos:length(b)] == as.raw(0x00))
  if (length(nul_rel) == 0) return(NA_character_)
  end <- pos + nul_rel[1] - 2L
  if (end < pos) return(NA_character_)
  name_raw <- b[pos:end]
  paste(rawToChar(name_raw, multiple = TRUE), collapse = "")
}

is_probably_html_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(TRUE)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  b <- readBin(con, what = "raw", n = 64)
  if (length(b) == 0) return(TRUE)
  txt <- tolower(paste(rawToChar(b, multiple = TRUE), collapse = ""))
  grepl("<!doctype", txt, fixed = TRUE) || grepl("<html", txt, fixed = TRUE)
}

file_has_gzip_magic <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 2) return(FALSE)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  b <- readBin(con, what = "raw", n = 2)
  length(b) == 2 && identical(b[1], as.raw(0x1f)) && identical(b[2], as.raw(0x8b))
}

# Download a Google Drive file ID, handling the "virus scan warning" HTML flow.
# Uses the R 'curl' package (avoids calling system curl).
drive_download <- function(file_id, destfile) {
  dir.create(dirname(destfile), showWarnings = FALSE, recursive = TRUE)
  tmpdir <- file.path(dirname(destfile), "_tmp_drive")
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
  tmp <- tempfile(fileext = ".bin", tmpdir = tmpdir)
  on.exit({
    unlink(tmp)
  }, add = TRUE)

  if (!requireNamespace("curl", quietly = TRUE)) {
    return(FALSE)
  }

  url1 <- paste0("https://drive.google.com/uc?export=download&id=", file_id)

  # Attempt direct download first (works for smaller files).
  ok1 <- tryCatch({
    curl::curl_download(url1, destfile = tmp, quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok1) return(FALSE)
  if (file_has_gzip_magic(tmp)) {
    file.copy(tmp, destfile, overwrite = TRUE)
    return(TRUE)
  }
  if (!is_probably_html_file(tmp)) return(FALSE)

  html <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  # Virus-scan warning page: extract form action + hidden inputs
  action <- NA_character_
  m_act <- regexpr("(?s)<form[^>]+action=\\\"([^\\\"]+)\\\"", html, perl = TRUE)
  if (m_act[1] > 0) {
    action <- regmatches(html, m_act)
    action <- sub("(?s)<form[^>]+action=\\\"([^\\\"]+)\\\".*", "\\1", action, perl = TRUE)
  }
  if (is.na(action) || !nzchar(action)) action <- "https://drive.usercontent.google.com/download"
  # extract hidden inputs: name/value pairs
  inp_re <- "<input[^>]+type=\\\"hidden\\\"[^>]+name=\\\"([^\\\"]+)\\\"[^>]+value=\\\"([^\\\"]*)\\\""
  m_inp <- gregexpr(inp_re, html, perl = TRUE)
  hits <- regmatches(html, m_inp)[[1]]
  if (length(hits) == 0) return(FALSE)
  names_v <- sub(inp_re, "\\1", hits, perl = TRUE)
  vals_v <- sub(inp_re, "\\2", hits, perl = TRUE)
  params <- stats::setNames(vals_v, names_v)
  # Ensure required params exist
  if (is.null(params[["id"]])) params[["id"]] <- file_id
  if (is.null(params[["export"]])) params[["export"]] <- "download"
  # Keep only the fields that Drive's warning form expects (avoid 400s).
  keep_keys <- intersect(names(params), c("id", "export", "confirm", "uuid"))
  params <- params[keep_keys]
  if (is.null(params[["id"]])) params[["id"]] <- file_id
  if (is.null(params[["export"]])) params[["export"]] <- "download"
  # Build query string
  # Build query string. Values here are simple (ids, uuid, 't'); avoid over-encoding.
  q <- paste0(
    names(params),
    "=",
    unname(params),
    collapse = "&"
  )
  url2 <- paste0(action, "?", q)

  # Download from the computed usercontent URL.
  ok2 <- tryCatch({
    curl::curl_download(url2, destfile = tmp, quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Drive confirm download failed: ", conditionMessage(e))
    FALSE
  })
  if (!ok2) return(FALSE)
  if (!file_has_gzip_magic(tmp)) return(FALSE)

  file.copy(tmp, destfile, overwrite = TRUE)
  TRUE
}

extract_drive_tpm_id_map <- function(html) {
  # Build a mapping from TCGA_*.tsv.gz filename -> Drive file_id from the folder HTML.
  # We look for a nearby "application/x-gzip" block and then a TCGA filename.
  q <- "(?:&quot;|\\\")"
  re <- paste0(
    "\\[\\[null,", q, "([A-Za-z0-9_-]{20,})", q,
    "\\],null,null,null,", q, "application/x-gzip", q,
    "[\\s\\S]{0,260}?", q, "(TCGA_[^&\\\"<>]+?_tpm\\.fullIDs\\.remapped\\.tsv\\.gz)", q
  )
  m <- gregexpr(re, html, perl = TRUE)
  hits <- regmatches(html, m)[[1]]
  if (length(hits) == 0) return(setNames(character(0), character(0)))
  ids <- sub(re, "\\1", hits, perl = TRUE)
  fns <- sub(re, "\\2", hits, perl = TRUE)
  # Prefer first occurrence per filename
  out <- ids[!duplicated(fns)]
  names(out) <- fns[!duplicated(fns)]
  out
}

find_drive_file_id <- function(filename, tpm_id_map) {
  if (length(tpm_id_map) == 0) return(NA_character_)
  if (!filename %in% names(tpm_id_map)) return(NA_character_)
  fid <- unname(tpm_id_map[filename])
  if (is.na(fid) || !nzchar(fid)) NA_character_ else fid
}

download_one_tpm <- function(tcga_code, precog_label, tpm_id_map) {
  file_name <- paste0("TCGA_", tcga_code, "_tpm.fullIDs.remapped.tsv.gz")
  local_path <- file.path(tpm_dir, file_name)
  if (file.exists(local_path) && file.info(local_path)$size > 1024 * 1024) {
    cat("Using existing TPM:", local_path, "\n")
    return(local_path)
  }

  cat("Downloading:", file_name, "...\n")
  fid <- find_drive_file_id(file_name, tpm_id_map)
  if (is.na(fid)) {
    cat("  Skipping; could not find downloadable Drive file id for:", file_name, "\n", sep = "")
    return(NA_character_)
  }

  ok <- tryCatch({
    drive_download(fid, local_path)
  }, error = function(e) FALSE)

  if (!ok || !file.exists(local_path) || file.info(local_path)$size < 1024 * 1024) {
    cat("  Download failed or produced too-small file for:", file_name, "\n", sep = "")
    return(NA_character_)
  }

  # Validate we downloaded the correct TCGA file (not a different gzip).
  want <- sub("\\.gz$", "", file_name, ignore.case = TRUE)
  got <- gzip_original_name_from_file(local_path)
  if (is.na(got) || !identical(got, want)) {
    cat("  Downloaded gzip does not match expected file. Expected=", want, " got=", got, "\n", sep = "")
    try(unlink(local_path), silent = TRUE)
    return(NA_character_)
  }

  local_path
}

# Genes used in PhenoMap scoring: |PRECOG z| > cutoff and gene in TPM (same rule as weighted_sum_scoring)
count_precog_genes_scoring <- function(refmat, tpm_genes, precog_label, z_score_cutoff) {
  if (!precog_label %in% colnames(refmat)) return(NA_integer_)
  z <- suppressWarnings(as.numeric(as.vector(refmat[, precog_label, drop = TRUE])))
  names(z) <- rownames(refmat)
  ok <- !is.na(z) & abs(z) > z_score_cutoff
  length(intersect(names(z)[ok], tpm_genes))
}

score_tpm_to_sample_scores <- function(tpm_file,
                                       tcga_code,
                                       precog_label_primary,
                                       precog_label_metastatic = NA_character_,
                                       z_score_cutoff,
                                       reference_mat) {
  # fread can read gz directly.
  dt <- fread(tpm_file, showProgress = FALSE, fill = TRUE)
  if (ncol(dt) < 2) {
    stop("Expression file parsed with <2 columns (Gene + samples).")
  }
  # Ensure unique column names (some TCGA files can contain duplicate sample IDs)
  cn <- names(dt)
  cn_u <- make.unique(cn, sep = "_dup")
  if (!identical(cn, cn_u)) {
    setnames(dt, cn_u)
  }
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- dt[[1]]
  tpm_genes <- rownames(mat)
  n_genes_primary <- count_precog_genes_scoring(reference_mat, tpm_genes, precog_label_primary, z_score_cutoff)
  n_genes_met <- if (!is.na(precog_label_metastatic)) {
    count_precog_genes_scoring(reference_mat, tpm_genes, precog_label_metastatic, z_score_cutoff)
  } else {
    NA_integer_
  }

  # Score all samples with primary reference
  scores_primary <- PhenoMap(
    expression = mat,
    reference = reference,
    cancer_type = precog_label_primary,
    z_score_cutoff = z_score_cutoff,
    verbose = FALSE
  )
  score_col_primary <- colnames(scores_primary)[1]

  # Optional: score all samples with metastatic reference, if available
  scores_met <- NULL
  score_col_met <- NA_character_
  if (!is.na(precog_label_metastatic)) {
    scores_met <- PhenoMap(
      expression = mat,
      reference = reference,
      cancer_type = precog_label_metastatic,
      z_score_cutoff = z_score_cutoff,
      verbose = FALSE
    )
    score_col_met <- colnames(scores_met)[1]
  }

  sample_ids <- rownames(scores_primary)
  sample_type <- sample_to_type_code(sample_ids)
  stratum <- sample_stratum(sample_type)

  score_primary <- as.numeric(scores_primary[[score_col_primary]])
  score_used <- score_primary
  precog_used <- rep(precog_label_primary, length(sample_ids))

  if (!is.null(scores_met)) {
    is_met_sample <- stratum == "Metastatic"
    score_met <- as.numeric(scores_met[[score_col_met]])
    score_used[is_met_sample] <- score_met[is_met_sample]
    precog_used[is_met_sample] <- precog_label_metastatic
  }

  list(
    sample_dt = data.table(
      tcga_code = tcga_code,
      stratum = stratum,
      patient = sample_to_patient(sample_ids),
      sample_id = sample_ids,
      sample_type = sample_type,
      precog_label_used = precog_used,
      score = score_used
    )[!is.na(patient)],
    n_genes_primary = n_genes_primary,
    n_genes_met = n_genes_met
  )
}

tcga_types <- sort(intersect(unique(as.character(liu_dt$type)), unique(unname(precog_to_tcga))))
if (!is.null(clin)) {
  tcga_types <- sort(intersect(tcga_types, unique(clin$acronym)))
}
tcga_only <- get_flag("--tcga_only", "")
if (is.character(tcga_only) && nzchar(tcga_only)) {
  tcga_types <- intersect(tcga_types, strsplit(tcga_only, ",", fixed = TRUE)[[1]])
  cat("Limiting to tcga_only:", paste(tcga_types, collapse = ", "), "\n")
}

cat("Fetching Drive folder HTML once for id lookup...\n")
dir.create(tpm_dir, showWarnings = FALSE, recursive = TRUE)
html_tmp <- tempfile(fileext = ".html")
system2("curl", c("-Ls", drive_html_url, "-o", html_tmp), stdout = NULL, stderr = NULL)
html <- paste(readLines(html_tmp, warn = FALSE), collapse = "\n")
unlink(html_tmp)
tpm_id_map <- extract_drive_tpm_id_map(html)
cat("Found TCGA TPM entries in folder HTML:", length(tpm_id_map), "\n")

panel_all <- list()
type_results <- list()

for (tcga_code in tcga_types) {
  precog_label_primary <- concordant_precog_label(tcga_code)
  if (is.na(precog_label_primary)) {
    cat("Skipping", tcga_code, ": no concordant PRECOG label.\n")
    next
  }
  precog_label_met <- metastatic_precog_label(tcga_code)
  cat("\n=== Processing TCGA:", tcga_code,
      "PRECOG(primary):", precog_label_primary,
      if (!is.na(precog_label_met)) paste0(" PRECOG(met): ", precog_label_met) else "",
      "===\n", sep = " ")
  liu_sub <- liu_dt[type == tcga_code]
  os_dt <- build_survival_liu(liu_sub)
  endpoint_lab <- if (tcga_code %in% USE_PFI_CANCERS) "PFI" else "OS"
  if (nrow(os_dt) < 20) {
    cat("  Skipping; too few ", endpoint_lab, " patients: ", nrow(os_dt), "\n", sep = "")
    next
  }

  tpm_local <- download_one_tpm(tcga_code, precog_label_primary, tpm_id_map)
  if (is.na(tpm_local)) next

  panel_dt <- tryCatch({
    scr <- score_tpm_to_sample_scores(
      tpm_file = tpm_local,
      tcga_code = tcga_code,
      precog_label_primary = precog_label_primary,
      precog_label_metastatic = precog_label_met,
      z_score_cutoff = z_score_cutoff,
      reference_mat = reference_mat
    )
    score_dt <- scr$sample_dt

    dat <- merge(os_dt, score_dt, by = "patient", all = FALSE)
    if (nrow(dat) < 20) NULL else {

    # Per-stratum median split + survival tests; keep strata separate (e.g. ACC_Normal)
    out_panels <- list()
    for (st in sort(unique(dat$stratum))) {
      sub <- dat[stratum == st]
      if (nrow(sub) < 20) next

      med <- median(sub$score, na.rm = TRUE)
      sub[, score_grp := factor(ifelse(score >= med, "High", "Low"), levels = c("Low", "High"))]
      sub[, median_score := med]
      sub[, forest_label := ifelse(st == "Tumor", tcga_code, paste0(tcga_code, "_", st))]

      fit <- coxph(Surv(os_time, os_event) ~ score_grp, data = sub)
      coef_hi <- coef(fit)[["score_grpHigh"]]
      hr <- exp(coef_hi)
      ci <- exp(confint(fit)[1, ])

      lr <- survdiff(Surv(os_time, os_event) ~ score_grp, data = sub)
      p_lr <- 1 - pchisq(lr$chisq, 1)

      n_precog_genes <- if (st == "Metastatic" && !is.na(precog_label_met) && is.finite(scr$n_genes_met)) {
        as.integer(scr$n_genes_met)
      } else {
        as.integer(scr$n_genes_primary)
      }

      type_results[[paste(tcga_code, st, sep = "__")]] <- data.table(
        forest_label = unique(sub$forest_label),
        tcga_code = tcga_code,
        stratum = st,
        survival_endpoint = endpoint_lab,
        precog_label_used = if (st == "Metastatic" && !is.na(precog_label_met)) precog_label_met else precog_label_primary,
        n_samples = nrow(sub),
        n_patients = uniqueN(sub$patient),
        n_precog_genes = n_precog_genes,
        median_score = med,
        hr = hr,
        ci_low = ci[1],
        ci_high = ci[2],
        p_logrank = p_lr
      )

      out_panels[[st]] <- sub
    }

      if (length(out_panels) == 0) NULL else {
        rbindlist(out_panels, use.names = TRUE, fill = TRUE)
      }
    }
  }, error = function(e) {
    cat("  ERROR processing ", tcga_code, ": ", conditionMessage(e), "\n", sep = "")
    NULL
  })

  # Always delete local tpm after processing this tcga_code
  if (file.exists(tpm_local)) {
    cat("  Deleting local TPM:", tpm_local, "\n")
    try(unlink(tpm_local), silent = TRUE)
  }

  if (!is.null(panel_dt)) panel_all[[tcga_code]] <- panel_dt
}

if (length(panel_all) == 0 || length(type_results) == 0) {
  stop("No successful TCGA types processed. Check downloads/concordance mapping.")
}

panel_dt_all <- rbindlist(panel_all, use.names = TRUE, fill = TRUE)
type_dt <- rbindlist(type_results, use.names = TRUE, fill = TRUE)
type_dt <- type_dt[order(hr, decreasing = TRUE)]

# Forest plot
plot_dt <- type_dt[
  stratum != "Normal" &
    is.finite(hr) & !is.na(hr) &
    is.finite(ci_low) & !is.na(ci_low) &
    is.finite(ci_high) & !is.na(ci_high)
]
plot_dt[, sig := ifelse(is.finite(p_logrank) & !is.na(p_logrank) & p_logrank < 0.05, "Significant", "Not significant")]
plot_dt$sig <- factor(plot_dt$sig, levels = c("Significant", "Not significant"))

# Point size: Pearson r between PRECOG and TCGA meta-z (same cancer mapping), using
# unshrunk files when --precog_full / --tcga_full exist (else bundled data).
load_full_meta_z_for_cor <- function(precog_path, tcga_path) {
  pf <- NULL
  tf <- NULL
  if (file.exists(precog_path)) {
    pf <- readRDS(precog_path)
    cat("Forest plot: using full PRECOG meta-z for correlations:", precog_path, "\n")
  } else {
    cat("Forest plot: full PRECOG file not found:", precog_path, "(using bundled precog for r)\n")
  }
  if (file.exists(tcga_path)) {
    df <- utils::read.csv(tcga_path, check.names = FALSE, stringsAsFactors = FALSE)
    gn <- df[[1L]]
    rownames(df) <- gn
    df[[1L]] <- NULL
    for (nm in c("PRECOG_metaZ", "TCGA_metaZ")) {
      if (nm %in% names(df)) df[[nm]] <- NULL
    }
    tf <- as.matrix(df)
    storage.mode(tf) <- "double"
    cat("Forest plot: using full TCGA meta-z CSV for correlations:", tcga_path, "\n")
  } else {
    cat("Forest plot: full TCGA CSV not found:", tcga_path, "(using bundled tcga for r)\n")
  }
  list(precog = pf, tcga = tf)
}

full_cor <- load_full_meta_z_for_cor(precog_full_path, tcga_full_path)
precog_cor <- full_cor$precog
tcga_cor <- full_cor$tcga
precog_ref_b <- PhenoMapR:::get_data("precog")
tcga_ref_b <- PhenoMapR:::get_data("tcga")

calc_precog_tcga_cor_one <- function(pf, tf, precog_label, tcga_code, min_genes = 50L) {
  if (is.null(pf) || is.null(tf)) return(NA_real_)
  if (is.na(precog_label) || is.na(tcga_code)) return(NA_real_)
  if (!precog_label %in% colnames(pf)) return(NA_real_)
  if (!tcga_code %in% colnames(tf)) return(NA_real_)
  g <- intersect(rownames(pf), rownames(tf))
  if (length(g) < min_genes) return(NA_real_)
  v1 <- suppressWarnings(as.numeric(as.vector(pf[g, precog_label, drop = TRUE])))
  v2 <- suppressWarnings(as.numeric(as.vector(tf[g, tcga_code, drop = TRUE])))
  if (length(v1) != length(g) || length(v2) != length(g)) return(NA_real_)
  ok <- is.finite(v1) & is.finite(v2)
  if (sum(ok) < min_genes) return(NA_real_)
  suppressWarnings(stats::cor(v1[ok], v2[ok], method = "pearson"))
}

calc_precog_tcga_cor <- function(precog_label, tcga_code) {
  r <- calc_precog_tcga_cor_one(precog_cor, tcga_cor, precog_label, tcga_code)
  if (is.finite(r)) return(r)
  # Fallback: shrunk package objects
  calc_precog_tcga_cor_one(precog_ref_b, tcga_ref_b, precog_label, tcga_code)
}

plot_dt[, precog_tcga_cor := mapply(calc_precog_tcga_cor, precog_label_used, tcga_code)]
plot_dt[, precog_tcga_cor_abs := abs(precog_tcga_cor)]
# Impute missing |cor| for sizing only (keeps all strata on plot; rare when gene overlap is low)
med_cor_abs <- stats::median(plot_dt$precog_tcga_cor_abs, na.rm = TRUE)
if (!is.finite(med_cor_abs)) med_cor_abs <- 0.5
plot_dt[, precog_tcga_cor_abs_plot := precog_tcga_cor_abs]
plot_dt[!is.finite(precog_tcga_cor_abs_plot), precog_tcga_cor_abs_plot := med_cor_abs]

# Point outline: hollow for negative or missing PRECOG–TCGA r; filled for r >= 0
plot_dt[, cor_shape := fifelse(
  is.finite(precog_tcga_cor) & precog_tcga_cor >= 0,
  "r_pos",
  "r_neg"
)]
plot_dt$cor_shape <- factor(plot_dt$cor_shape, levels = c("r_pos", "r_neg"))

# Shared y order for forest + patient-count bar (factor so both panels align)
plot_dt[, y_ord := stats::reorder(forest_label, hr)]

plot_dt[, z_score_cutoff := z_score_cutoff]
plot_dt[, reference := reference]

type_dt[, included_in_forest_plot := forest_label %in% plot_dt$forest_label]
plot_extras <- c(
  "sig", "precog_tcga_cor", "precog_tcga_cor_abs", "precog_tcga_cor_abs_plot",
  "cor_shape", "y_ord"
)
type_dt_export <- merge(
  type_dt,
  plot_dt[, c("forest_label", plot_extras), with = FALSE],
  by = "forest_label",
  all.x = TRUE
)
type_dt_export[, z_score_cutoff := z_score_cutoff]
type_dt_export[, reference := reference]

panel_out <- file.path(out_dir, "tcga_precog_forest_panel.tsv")
type_out <- file.path(out_dir, "tcga_precog_forest_results.tsv")
plot_data_out <- file.path(out_dir, "tcga_precog_forest_plot_data.tsv")
samples_ann_out <- file.path(out_dir, "tcga_precog_forest_samples_annotated.tsv")
note_out <- file.path(out_dir, "tcga_outcomes_methods_note.txt")

fwrite(panel_dt_all, panel_out, sep = "\t")
fwrite(type_dt_export, type_out, sep = "\t")
fwrite(plot_dt, plot_data_out, sep = "\t")

panel_enriched <- enrich_panel_samples(panel_dt_all, liu_dt, clin)
panel_enriched[, z_score_cutoff := z_score_cutoff]
panel_enriched[, reference := reference]
fwrite(panel_enriched, samples_ann_out, sep = "\t")

write_outcomes_note(note_out)

cat("\nWrote sample panel (minimal, tertile-compatible):", panel_out, "\n")
cat("Wrote stratum-level results + forest-plot fields:", type_out, "\n")
cat("Wrote exact forest figure data:", plot_data_out, "\n")
cat("Wrote annotated samples (Liu + PANCAN + PhenoMapR):", samples_ann_out, "\n")
cat("Wrote outcomes note (Liu vs PANCAN):", note_out, "\n")

# Point size: |cor| mapped to [0, 1] on the scale; visual point size range still scales with n rows
cor_lim <- c(0, 1)
n_forest <- nrow(plot_dt)
size_rng <- c(
  max(1.1, 1.35 + 0.02 * n_forest),
  min(7, 2.05 + 0.14 * n_forest)
)

# Larger text throughout the figure (axis, legend, y category labels)
plot_base_size <- 16
plot_title_size <- 18
legend_pt <- 4.5

p_forest <- ggplot(plot_dt, aes(y = y_ord, x = hr, color = sig, size = precog_tcga_cor_abs_plot, shape = cor_shape)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.6) +
  geom_point(data = ~ .x[.x$cor_shape == "r_pos", ], stroke = 0) +
  geom_point(data = ~ .x[.x$cor_shape == "r_neg", ], stroke = 0.85) +
  scale_x_log10() +
  scale_shape_manual(
    values = c(r_pos = 16, r_neg = 1),
    labels = c(r_pos = "r \u2265 0", r_neg = "r < 0"),
    name = "PRECOG vs. TCGA\nMeta-z correlation"
  ) +
  scale_size_continuous(
    range = size_rng,
    limits = cor_lim,
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    name = "|cor(PRECOG, TCGA)|"
  ) +
  scale_color_manual(
    values = c("Significant" = "#B2182B", "Not significant" = "#4D4D4D"),
    name = "Log-rank Survival"
  ) +
  labs(
    x = "Hazard Ratio\n(High vs. Low PhenoMapR score)",
    y = "TCGA cancer"
  ) +
  theme_minimal(base_size = plot_base_size) +
  guides(
    color = ggplot2::guide_legend(override.aes = list(size = legend_pt, stroke = 0)),
    shape = ggplot2::guide_legend(override.aes = list(size = legend_pt)),
    # Size legend: show |cor| 0–1 as a gradient of point sizes (do not fix size in override.aes)
    size = ggplot2::guide_legend(override.aes = list(shape = 16))
  )

p_counts <- ggplot(plot_dt, aes(y = y_ord, x = n_patients)) +
  geom_col(width = 0.72, fill = "#5C5C5C", alpha = 0.88) +
  labs(x = "# of TCGA\nPatients", y = NULL) +
  theme_minimal(base_size = plot_base_size) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

if (!"n_precog_genes" %in% names(plot_dt)) {
  plot_dt[, n_precog_genes := NA_integer_]
}
plot_dt[, n_precog_genes_plot := ifelse(is.finite(n_precog_genes), as.numeric(n_precog_genes), NA_real_)]

p_genes <- ggplot(plot_dt, aes(y = y_ord, x = n_precog_genes_plot)) +
  geom_col(width = 0.72, fill = "#4A6FA5", alpha = 0.88) +
  scale_x_continuous(
    labels = function(x) tolower(scales::label_number(scale_cut = scales::cut_short_scale())(x))
  ) +
  labs(x = "# of genes in\nPRECOG signature", y = NULL) +
  theme_minimal(base_size = plot_base_size) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

p <- p_forest + p_counts + p_genes +
  plot_layout(widths = c(3.4, 1, 1), guides = "collect") +
  plot_annotation(
    title = "TCGA survival by median PhenoMapR PRECOG score",
    subtitle = "Outcomes: Liu Cell 2018 TCGA-CDR (OS; PFI for DLBC, TGCT, READ)",
    theme = theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = plot_title_size),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = plot_base_size - 2)
    )
  )

forest_png <- file.path(out_dir, "tcga_precog_forest.png")
ggsave(forest_png, p, width = 12.5, height = max(4.8, 0.24 * n_forest + 2.2), dpi = 200)
cat("Wrote forest plot:", forest_png, "\n")

cat("\nTop results (sorted by HR):\n")
print(type_dt[, .(
  forest_label, tcga_code, stratum, survival_endpoint, precog_label_used,
  n_samples, n_patients, n_precog_genes, hr, ci_low, ci_high, p_logrank
)])

