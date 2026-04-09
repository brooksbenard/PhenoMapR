#!/usr/bin/env Rscript
#
# Incremental TCGA PRECOG bulk survival stratification:
# - download one TCGA TPM gz at a time
# - score via PhenoMapR (built-in PRECOG meta-z signatures)
# - stratify by median, mean, and max-selected log-rank (optimal cutpoint) PhenoMapR score *within each stratum*
# - compute Cox PH HR + CI; median/mean use log-rank p at split; max-stat uses HL global p + HR at optimal split
# - append sample-level results to a single consolidated table
# - delete the local TPM file after each cancer type is processed
#
# Expected local inputs:
# - PANCAN clinical follow-up: data/tcga/clinical_PANCAN_patient_with_followup.tsv (or --clinical_file)
#   Outcomes: overall survival (OS) only — vital_status (Dead/Alive), days_to_death,
#   days_to_last_followup (GDC-style fields). Progression fields in this PANCAN snapshot are
#   empty for DLBC, TGCT, and READ, so only OS is available for those types.
#
# Outputs:
# - results/tcga_precog_forest_results.tsv — median split: stratum-level Cox/log-rank + forest-plot fields
# - results/tcga_precog_forest_plot_data.tsv — median split: exact forest figure rows
# - results/tcga_precog_forest_mean_results.tsv — mean split (same structure; threshold column mean_score)
# - results/tcga_precog_forest_mean_plot_data.tsv — mean split forest figure rows
# - results/tcga_precog_forest_maxstat_results.tsv — max-stat optimal split (p_maxstat = HL-adjusted global test)
# - results/tcga_precog_forest_maxstat_plot_data.tsv — max-stat forest figure rows
# - results/tcga_precog_forest_panel.tsv — patient-level scores + survival + median/mean/maxstat group labels
# - results/tcga_precog_forest_samples_annotated.tsv — per-sample PhenoMapR scores + PANCAN clinical (pancan_*)
# - results/tcga_outcomes_methods_note.txt — outcomes source (PANCAN OS)
# - results/tcga_precog_forest.png — median split forest + patient n + PRECOG gene n
# - results/tcga_precog_forest_mean.png — mean split forest (same layout)
# - results/tcga_precog_forest_maxstat.png — max-stat optimal split forest (same layout)
# Optional flags: --precog_full, --tcga_full (unshrunk meta-z for PRECOG–TCGA r in the plot)

suppressPackageStartupMessages({
  library(data.table)
  library(PhenoMapR)
  library(survival)
  library(maxstat)
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
out_dir <- get_flag("--out_dir", "results")
# Unshrunk meta-z tables for PRECOG–TCGA Pearson r (forest point sizes); optional, repo root by default
precog_full_path <- get_flag("--precog_full", "PRECOG_V2_Cancer_Type_Meta_Zscores_Final.rds")
tcga_full_path <- get_flag("--tcga_full", "TCGA_metaz.csv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

reference_mat <- PhenoMapR:::get_data(reference)

if (!file.exists(clinical_file)) {
  stop(
    "Missing PANCAN clinical file: ", clinical_file,
    "\nPlace clinical_PANCAN_patient_with_followup.tsv under data/tcga/ or pass --clinical_file <path>."
  )
}
cat("Loading PANCAN clinical:", clinical_file, "\n")
clin <- fread(clinical_file, showProgress = FALSE)
if (!all(c("bcr_patient_barcode", "acronym", "vital_status", "days_to_death", "days_to_last_followup") %in% names(clin))) {
  stop("PANCAN clinical TSV must include bcr_patient_barcode, acronym, vital_status, days_to_death, days_to_last_followup.")
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

na_pancan_num <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "[Not Available]", "[Not Applicable]", "NA", "Unknown")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

#' One TCGA acronym from PANCAN clinical: OS time (days) and event (1 = dead, 0 = alive).
build_survival_pancan <- function(clin_sub) {
  if (nrow(clin_sub) == 0L) {
    return(data.table(
      patient = character(),
      os_time = numeric(),
      os_event = integer()
    ))
  }
  vs <- tolower(trimws(as.character(clin_sub$vital_status)))
  dtd <- na_pancan_num(clin_sub$days_to_death)
  dtlfu <- na_pancan_num(clin_sub$days_to_last_followup)

  dt <- data.table(
    patient = as.character(clin_sub$bcr_patient_barcode),
    vs = vs,
    dtd = dtd,
    dtlfu = dtlfu
  )

  dt <- dt[vs %in% c("dead", "alive")]
  if (nrow(dt) == 0L) {
    return(data.table(patient = character(), os_time = numeric(), os_event = integer()))
  }

  # Patient-level aggregation:
  # - event: 1 if any record says dead, else 0
  # - time: if dead -> max(days_to_death); if alive -> max(days_to_last_followup)
  dt_dead <- dt[vs == "dead" & is.finite(dtd), .(os_event = 1L, os_time = max(dtd, na.rm = TRUE)), by = patient]
  dt_alive <- dt[vs == "alive" & is.finite(dtlfu), .(os_event = 0L, os_time = max(dtlfu, na.rm = TRUE)), by = patient]
  out <- rbindlist(list(dt_dead, dt_alive), use.names = TRUE, fill = TRUE)
  out <- out[is.finite(os_time) & !is.na(os_event)]
  if (nrow(out) == 0L) return(out[, .(patient, os_time, os_event)])
  # If both dead and alive rows exist for a patient, keep dead.
  out <- out[order(-os_event, -os_time), .SD[1], by = patient]
  out[, .(patient, os_time, os_event)]
}

#' Per-sample panel + PANCAN clinical (pancan_*) + derived pancan_os_* (same OS definition as survival fit).
enrich_panel_samples <- function(panel, clin, outcomes_label) {
  p <- data.table::copy(panel)
  p[, survival_endpoint := "OS"]
  p[, survival_time_used_days := os_time]
  p[, survival_event_used := os_event]

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
    vs <- tolower(trimws(as.character(p$pancan_vital_status)))
    p[, pancan_os_time := ifelse(
      vs == "dead",
      na_pancan_num(pancan_days_to_death),
      ifelse(vs == "alive", na_pancan_num(pancan_days_to_last_followup), NA_real_)
    )]
    p[, pancan_os_event := NA_integer_]
    p[vs == "dead", pancan_os_event := 1L]
    p[vs == "alive", pancan_os_event := 0L]
  } else {
    p[, pancan_os_time := NA_real_]
    p[, pancan_os_event := NA_integer_]
  }

  p[, outcomes_table := outcomes_label]
  p
}

write_outcomes_note <- function(path, outcomes_label) {
  txt <- c(
    "TCGA PhenoMapR forest analysis — outcomes sources and hazard ratios",
    "=================================================================",
    "",
    "Primary survival endpoint for Cox / log-rank: overall survival (OS) from PANCAN clinical follow-up.",
    paste0("  File (basename): ", outcomes_label),
    "",
    "Definition (same as script):",
    "  - vital_status Dead: time = days_to_death, event = 1.",
    "  - vital_status Alive: time = days_to_last_followup, event = 0 (censored).",
    "  - Other vital_status values are excluded. OS is used for all TCGA types.",
    "",
    "DLBC, TGCT, READ:",
    "  This PANCAN snapshot has days_to_patient_progression_free / patient_progression_status",
    "  effectively empty for these cohorts, so progression-free endpoints cannot be derived from",
    "  this file. The forest plot uses OS only for every cancer type.",
    "",
    "The file tcga_precog_forest_samples_annotated.tsv includes:",
    "  - PhenoMapR score columns (score, precog_label_used, stratum, ...)",
    "  - PANCAN columns prefixed with pancan_ plus pancan_os_time / pancan_os_event (aligned with the fit).",
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

tcga_types <- sort(intersect(unique(as.character(clin$acronym)), unique(unname(precog_to_tcga))))
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
type_results_median <- list()
type_results_mean <- list()
type_results_maxstat <- list()

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
  clin_sub <- clin[acronym == tcga_code]
  if (tcga_code %in% c("DLBC", "TGCT", "READ")) {
    n_pfi_time <- sum(is.finite(na_pancan_num(clin_sub$days_to_patient_progression_free)))
    cat(
      "  Outcomes note: ", tcga_code, " — PANCAN progression fields non-missing: ",
      n_pfi_time, " / ", nrow(clin_sub), " (days_to_patient_progression_free); using OS only.\n",
      sep = ""
    )
  }
  os_dt <- build_survival_pancan(clin_sub)
  endpoint_lab <- "OS"
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

    # IMPORTANT: survival is at patient-level. TCGA has multiple samples per patient, so we
    # aggregate scores to one row per patient *within each stratum* before fitting Cox/log-rank.
    score_dt_patient <- score_dt[
      is.finite(score) & !is.na(score),
      .(
        score = mean(score, na.rm = TRUE),
        precog_label_used = precog_label_used[1],
        tcga_code = tcga_code[1]
      ),
      by = .(patient, stratum)
    ]

    dat <- merge(os_dt, score_dt_patient, by = "patient", all = FALSE)
    if (nrow(dat) < 20) NULL else {

    # Per-stratum median and mean splits + survival tests; keep strata separate (e.g. ACC_Normal)
    out_panels <- list()
    for (st in sort(unique(dat$stratum))) {
      sub <- dat[stratum == st]
      if (nrow(sub) < 20) next

      med <- median(sub$score, na.rm = TRUE)
      mn <- mean(sub$score, na.rm = TRUE)
      forest_label_st <- ifelse(st == "Tumor", tcga_code, paste0(tcga_code, "_", st))
      precog_used_st <- if (st == "Metastatic" && !is.na(precog_label_met)) precog_label_met else precog_label_primary

      n_precog_genes <- if (st == "Metastatic" && !is.na(precog_label_met) && is.finite(scr$n_genes_met)) {
        as.integer(scr$n_genes_met)
      } else {
        as.integer(scr$n_genes_primary)
      }

      # Median split (score >= median -> High)
      sub_med <- data.table::copy(sub)
      sub_med[, score_grp := factor(ifelse(score >= med, "High", "Low"), levels = c("Low", "High"))]
      fit <- coxph(Surv(os_time, os_event) ~ score_grp, data = sub_med)
      coef_hi <- coef(fit)[["score_grpHigh"]]
      hr <- exp(coef_hi)
      ci <- exp(confint(fit)[1, ])
      lr <- survdiff(Surv(os_time, os_event) ~ score_grp, data = sub_med)
      p_lr <- 1 - pchisq(lr$chisq, 1)

      type_results_median[[paste(tcga_code, st, sep = "__")]] <- data.table(
        forest_label = forest_label_st,
        tcga_code = tcga_code,
        stratum = st,
        survival_endpoint = endpoint_lab,
        precog_label_used = precog_used_st,
        n_samples = nrow(sub_med),
        n_patients = uniqueN(sub_med$patient),
        n_precog_genes = n_precog_genes,
        median_score = med,
        hr = hr,
        ci_low = ci[1],
        ci_high = ci[2],
        p_logrank = p_lr
      )

      # Mean split (score >= mean -> High)
      sub_mean <- data.table::copy(sub)
      sub_mean[, score_grp := factor(ifelse(score >= mn, "High", "Low"), levels = c("Low", "High"))]
      fit_m <- coxph(Surv(os_time, os_event) ~ score_grp, data = sub_mean)
      coef_hi_m <- coef(fit_m)[["score_grpHigh"]]
      hr_m <- exp(coef_hi_m)
      ci_m <- exp(confint(fit_m)[1, ])
      lr_m <- survdiff(Surv(os_time, os_event) ~ score_grp, data = sub_mean)
      p_lr_m <- 1 - pchisq(lr_m$chisq, 1)

      type_results_mean[[paste(tcga_code, st, sep = "__")]] <- data.table(
        forest_label = forest_label_st,
        tcga_code = tcga_code,
        stratum = st,
        survival_endpoint = endpoint_lab,
        precog_label_used = precog_used_st,
        n_samples = nrow(sub_mean),
        n_patients = uniqueN(sub_mean$patient),
        n_precog_genes = n_precog_genes,
        mean_score = mn,
        hr = hr_m,
        ci_low = ci_m[1],
        ci_high = ci_m[2],
        p_logrank = p_lr_m
      )

      # Max-selected log-rank split (optimal cutpoint; global p from Hothorn–Lausen HL approximation)
      sub_ms_ok <- FALSE
      sub_ms <- NULL
      cp_ms <- NA_real_
      ms <- tryCatch(
        maxstat::maxstat.test(
          survival::Surv(os_time, os_event) ~ score,
          data = sub,
          smethod = "LogRank",
          pmethod = "HL",
          minprop = 0.1,
          maxprop = 0.9
        ),
        error = function(e) NULL
      )
      if (!is.null(ms)) {
        cp_candidate <- unname(ms$estimate["estimated cutpoint"])
        if (length(cp_candidate) == 1L && is.finite(cp_candidate)) {
          sub_ms <- data.table::copy(sub)
          sub_ms[, score_grp := factor(fifelse(score > cp_candidate, "High", "Low"), levels = c("Low", "High"))]
          if (uniqueN(sub_ms$score_grp) >= 2L) {
            fit_ms <- tryCatch(
              coxph(Surv(os_time, os_event) ~ score_grp, data = sub_ms),
              error = function(e) NULL
            )
            if (!is.null(fit_ms) && "score_grpHigh" %in% names(coef(fit_ms))) {
              coef_hi_ms <- coef(fit_ms)[["score_grpHigh"]]
              if (is.finite(coef_hi_ms)) {
                hr_ms <- exp(coef_hi_ms)
                ci_ms_x <- exp(confint(fit_ms)[1, ])
                lr_ms <- survdiff(Surv(os_time, os_event) ~ score_grp, data = sub_ms)
                p_lr_split <- 1 - pchisq(lr_ms$chisq, 1)
                cp_ms <- as.numeric(cp_candidate)
                p_ms <- as.numeric(ms$p.value)
                type_results_maxstat[[paste(tcga_code, st, sep = "__")]] <- data.table(
                  forest_label = forest_label_st,
                  tcga_code = tcga_code,
                  stratum = st,
                  survival_endpoint = endpoint_lab,
                  precog_label_used = precog_used_st,
                  n_samples = nrow(sub_ms),
                  n_patients = uniqueN(sub_ms$patient),
                  n_precog_genes = n_precog_genes,
                  maxstat_cutpoint = cp_ms,
                  p_maxstat = p_ms,
                  p_logrank_at_split = p_lr_split,
                  hr = hr_ms,
                  ci_low = ci_ms_x[1],
                  ci_high = ci_ms_x[2],
                  p_logrank = p_ms
                )
                sub_ms_ok <- TRUE
              }
            }
          }
        }
      }

      # Panel: patient-level row; score_grp / median_score = median; score_grp_mean / mean_score = mean; max-stat if fit OK
      sub_panel <- data.table::copy(sub)
      sub_panel[, forest_label := forest_label_st]
      sub_panel[, score_grp := factor(ifelse(score >= med, "High", "Low"), levels = c("Low", "High"))]
      sub_panel[, median_score := med]
      sub_panel[, score_grp_mean := factor(ifelse(score >= mn, "High", "Low"), levels = c("Low", "High"))]
      sub_panel[, mean_score := mn]
      if (isTRUE(sub_ms_ok)) {
        sub_panel[, score_grp_maxstat := as.character(sub_ms$score_grp)]
        sub_panel[, maxstat_cutpoint := cp_ms]
      } else {
        sub_panel[, score_grp_maxstat := NA_character_]
        sub_panel[, maxstat_cutpoint := NA_real_]
      }
      out_panels[[st]] <- sub_panel
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

if (length(panel_all) == 0 || length(type_results_median) == 0 || length(type_results_mean) == 0) {
  stop("No successful TCGA types processed. Check downloads/concordance mapping.")
}

panel_dt_all <- rbindlist(panel_all, use.names = TRUE, fill = TRUE)
type_dt_median <- rbindlist(type_results_median, use.names = TRUE, fill = TRUE)
type_dt_median <- type_dt_median[order(hr, decreasing = TRUE)]
type_dt_mean <- rbindlist(type_results_mean, use.names = TRUE, fill = TRUE)
type_dt_mean <- type_dt_mean[order(hr, decreasing = TRUE)]

has_maxstat <- length(type_results_maxstat) > 0L
if (!has_maxstat) {
  warning("No max-stat strata produced results; skipping max-stat forest outputs.")
}
type_dt_maxstat <- if (has_maxstat) {
  x <- rbindlist(type_results_maxstat, use.names = TRUE, fill = TRUE)
  x[order(hr, decreasing = TRUE)]
} else {
  NULL
}

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

finalize_forest_plot_dt <- function(type_dt) {
  plot_dt <- type_dt[
    stratum != "Normal" &
      is.finite(hr) & !is.na(hr) &
      is.finite(ci_low) & !is.na(ci_low) &
      is.finite(ci_high) & !is.na(ci_high)
  ]
  # Median/mean: p_logrank = log-rank at fixed split. Max-stat: p_logrank duplicates p_maxstat (HL global p).
  plot_dt[, sig := ifelse(is.finite(p_logrank) & !is.na(p_logrank) & p_logrank < 0.05, "Significant", "Not significant")]
  plot_dt$sig <- factor(plot_dt$sig, levels = c("Significant", "Not significant"))
  plot_dt[, precog_tcga_cor := mapply(calc_precog_tcga_cor, precog_label_used, tcga_code)]
  plot_dt[, precog_tcga_cor_abs := abs(precog_tcga_cor)]
  med_cor_abs <- stats::median(plot_dt$precog_tcga_cor_abs, na.rm = TRUE)
  if (!is.finite(med_cor_abs)) med_cor_abs <- 0.5
  plot_dt[, precog_tcga_cor_abs_plot := precog_tcga_cor_abs]
  plot_dt[!is.finite(precog_tcga_cor_abs_plot), precog_tcga_cor_abs_plot := med_cor_abs]
  plot_dt[, cor_shape := fifelse(
    is.finite(precog_tcga_cor) & precog_tcga_cor >= 0,
    "r_pos",
    "r_neg"
  )]
  plot_dt$cor_shape <- factor(plot_dt$cor_shape, levels = c("r_pos", "r_neg"))
  plot_dt[, y_ord := stats::reorder(forest_label, hr)]
  plot_dt[, z_score_cutoff := z_score_cutoff]
  plot_dt[, reference := reference]
  plot_dt
}

plot_dt_median <- finalize_forest_plot_dt(type_dt_median)
plot_dt_mean <- finalize_forest_plot_dt(type_dt_mean)
plot_dt_maxstat <- if (has_maxstat) finalize_forest_plot_dt(type_dt_maxstat) else NULL

merge_type_with_plot_extras <- function(type_dt, plot_dt) {
  x <- data.table::copy(type_dt)
  x[, included_in_forest_plot := forest_label %in% plot_dt$forest_label]
  plot_extras <- c(
    "sig", "precog_tcga_cor", "precog_tcga_cor_abs", "precog_tcga_cor_abs_plot",
    "cor_shape", "y_ord"
  )
  out <- merge(
    x,
    plot_dt[, c("forest_label", plot_extras), with = FALSE],
    by = "forest_label",
    all.x = TRUE
  )
  out[, z_score_cutoff := z_score_cutoff]
  out[, reference := reference]
  out
}

type_dt_export_median <- merge_type_with_plot_extras(type_dt_median, plot_dt_median)
type_dt_export_mean <- merge_type_with_plot_extras(type_dt_mean, plot_dt_mean)
type_dt_export_maxstat <- if (has_maxstat) {
  merge_type_with_plot_extras(type_dt_maxstat, plot_dt_maxstat)
} else {
  NULL
}

panel_out <- file.path(out_dir, "tcga_precog_forest_panel.tsv")
type_out_median <- file.path(out_dir, "tcga_precog_forest_results.tsv")
plot_data_out_median <- file.path(out_dir, "tcga_precog_forest_plot_data.tsv")
type_out_mean <- file.path(out_dir, "tcga_precog_forest_mean_results.tsv")
plot_data_out_mean <- file.path(out_dir, "tcga_precog_forest_mean_plot_data.tsv")
type_out_maxstat <- file.path(out_dir, "tcga_precog_forest_maxstat_results.tsv")
plot_data_out_maxstat <- file.path(out_dir, "tcga_precog_forest_maxstat_plot_data.tsv")
samples_ann_out <- file.path(out_dir, "tcga_precog_forest_samples_annotated.tsv")
note_out <- file.path(out_dir, "tcga_outcomes_methods_note.txt")

fwrite(panel_dt_all, panel_out, sep = "\t")
fwrite(type_dt_export_median, type_out_median, sep = "\t")
fwrite(plot_dt_median, plot_data_out_median, sep = "\t")
fwrite(type_dt_export_mean, type_out_mean, sep = "\t")
fwrite(plot_dt_mean, plot_data_out_mean, sep = "\t")
if (has_maxstat) {
  fwrite(type_dt_export_maxstat, type_out_maxstat, sep = "\t")
  fwrite(plot_dt_maxstat, plot_data_out_maxstat, sep = "\t")
}

outcomes_label <- basename(clinical_file)
panel_enriched <- enrich_panel_samples(panel_dt_all, clin, outcomes_label)
panel_enriched[, z_score_cutoff := z_score_cutoff]
panel_enriched[, reference := reference]
fwrite(panel_enriched, samples_ann_out, sep = "\t")

write_outcomes_note(note_out, outcomes_label)

cat("\nWrote patient-level panel (median, mean, max-stat group columns):", panel_out, "\n")
cat("Wrote median-split results + forest fields:", type_out_median, "\n")
cat("Wrote median-split forest figure data:", plot_data_out_median, "\n")
cat("Wrote mean-split results + forest fields:", type_out_mean, "\n")
cat("Wrote mean-split forest figure data:", plot_data_out_mean, "\n")
if (has_maxstat) {
  cat("Wrote max-stat results + forest fields:", type_out_maxstat, "\n")
  cat("Wrote max-stat forest figure data:", plot_data_out_maxstat, "\n")
}
cat("Wrote annotated samples (PANCAN + PhenoMapR):", samples_ann_out, "\n")
cat("Wrote outcomes note (PANCAN OS):", note_out, "\n")

build_tcga_forest_multipanel <- function(plot_dt, title, forest_png, color_legend_name = NULL) {
  legend_color_name <- if (is.null(color_legend_name)) "Log-rank Survival" else color_legend_name
  pd <- data.table::copy(plot_dt)
  cor_lim <- c(0, 1)
  n_forest <- nrow(pd)
  size_rng <- c(
    max(1.1, 1.35 + 0.02 * n_forest),
    min(7, 2.05 + 0.14 * n_forest)
  )
  plot_base_size <- 16
  plot_title_size <- 18
  legend_pt <- 4.5

  if (!"n_precog_genes" %in% names(pd)) {
    pd[, n_precog_genes := NA_integer_]
  }
  pd[, n_precog_genes_plot := ifelse(is.finite(n_precog_genes), as.numeric(n_precog_genes), NA_real_)]

  p_forest <- ggplot(pd, aes(y = y_ord, x = hr, color = sig, size = precog_tcga_cor_abs_plot, shape = cor_shape)) +
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
      name = legend_color_name
    ) +
    labs(
      x = "Hazard Ratio\n(High vs. Low PhenoMapR score)",
      y = "TCGA cancer"
    ) +
    theme_minimal(base_size = plot_base_size) +
    guides(
      color = ggplot2::guide_legend(override.aes = list(size = legend_pt, stroke = 0)),
      shape = ggplot2::guide_legend(override.aes = list(size = legend_pt)),
      size = ggplot2::guide_legend(override.aes = list(shape = 16))
    )

  p_counts <- ggplot(pd, aes(y = y_ord, x = n_patients)) +
    geom_col(width = 0.72, fill = "#5C5C5C", alpha = 0.88) +
    labs(x = "# of TCGA\nPatients", y = NULL) +
    theme_minimal(base_size = plot_base_size) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()
    )

  p_genes <- ggplot(pd, aes(y = y_ord, x = n_precog_genes_plot)) +
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
      title = title,
      theme = theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = plot_title_size)
      )
    )

  ggsave(forest_png, p, width = 12.5, height = max(4.8, 0.24 * n_forest + 2.2), dpi = 200)
  cat("Wrote forest plot:", forest_png, "\n")
  invisible(NULL)
}

build_tcga_forest_multipanel(
  plot_dt_median,
  "TCGA survival by median PhenoMapR PRECOG score",
  file.path(out_dir, "tcga_precog_forest.png")
)
build_tcga_forest_multipanel(
  plot_dt_mean,
  "TCGA survival by mean-split PhenoMapR PRECOG score",
  file.path(out_dir, "tcga_precog_forest_mean.png")
)
if (has_maxstat) {
  build_tcga_forest_multipanel(
    plot_dt_maxstat,
    "TCGA survival by max-stat optimal PhenoMapR PRECOG split",
    file.path(out_dir, "tcga_precog_forest_maxstat.png"),
    color_legend_name = "Max-selected log-rank\n(HL global p)"
  )
}

cat("\nMedian split — top results (sorted by HR):\n")
print(type_dt_median[, .(
  forest_label, tcga_code, stratum, survival_endpoint, precog_label_used,
  n_samples, n_patients, n_precog_genes, median_score, hr, ci_low, ci_high, p_logrank
)])
cat("\nMean split — top results (sorted by HR):\n")
print(type_dt_mean[, .(
  forest_label, tcga_code, stratum, survival_endpoint, precog_label_used,
  n_samples, n_patients, n_precog_genes, mean_score, hr, ci_low, ci_high, p_logrank
)])
if (has_maxstat) {
  cat("\nMax-stat split — top results (sorted by HR):\n")
  print(type_dt_maxstat[, .(
    forest_label, tcga_code, stratum, survival_endpoint, precog_label_used,
    n_samples, n_patients, n_precog_genes, maxstat_cutpoint, p_maxstat, p_logrank_at_split,
    hr, ci_low, ci_high
  )])
}
