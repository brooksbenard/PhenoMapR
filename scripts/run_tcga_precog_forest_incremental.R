#!/usr/bin/env Rscript
#
# Incremental TCGA PRECOG bulk survival stratification:
# - download one TCGA TPM gz at a time
# - score via PhenoMapR (built-in PRECOG meta-z signatures)
# - stratify by median and mean PhenoMapR score *within each stratum*
# - compute Cox PH HR + CI; median/mean use log-rank p at split
# - append sample-level results to a single consolidated table
# - delete the local TPM file after each cancer type is processed
#
# Expected local inputs:
# - TCGA Pan-Cancer Clinical Data Resource, Supplemental Table S1 sheet 1 (Excel):
#   default path data/tcga/TCGA-CDR-SupplementalTableS1.xlsx (or --cdr_file).
#   Uses columns type, bcr_patient_barcode, OS/OS.time, DSS/DSS.time, DFI/DFI.time, PFI/PFI.time.
#
# Outputs (one set per requested endpoint, suffix _{OS|DSS|DFI|PFI}):
# - results/tcga_precog_forest_{ep}_results.tsv — median split
# - results/tcga_precog_forest_{ep}_plot_data.tsv
# - results/tcga_precog_forest_{ep}_mean_results.tsv — mean split
# - results/tcga_precog_forest_{ep}_mean_plot_data.tsv
# - results/tcga_precog_forest_{ep}_panel.tsv — patient-level scores + survival + median/mean labels
# - results/tcga_precog_forest_{ep}_samples_annotated.tsv — scores + CDR covariates (cdr_*)
# - results/tcga_outcomes_methods_note.txt — CDR endpoint definitions
# - results/tcga_precog_forest_{ep}.png / tcga_precog_forest_{ep}_mean.png — forests
# - results/tcga_precog_forest_combined_dual.png — OS (PFI for PRAD, DLBC, LGG, READ, TGCT); PRECOG + TCGA HRs (dodged)
# - results/tcga_precog_forest_combined_dual_{precog_only,tcga_only}.png — same cohorts as combined dual; OS-style 3-panel forests (one signature each)
# - results/tcga_precog_forest_combined_dual_{results,plot_data}.tsv
# - results/tcga_metaz_signature_forest_{ep}* — same endpoints/strata as PRECOG forests; PhenoMapR reference = built-in TCGA meta-z (cancer_type = TCGA code)
# - results/tcga_metaz_signature_forest_combined.png (+ _results.tsv, _plot_data.tsv) — Tumor-only pan-cancer; PFI for BRCA, DLBC, KICH, LGG, PRAD, READ, TGCT, THCA, THYM else OS; CDR\u2229Drive TPM cohorts; min 10 patients. _results.tsv = raw Cox; _plot_data.tsv caps non-finite or HR/CI > 100 for drawing (plot_hr_ci_capped).
# Optional flags: --precog_full, --tcga_full; --endpoints OS,DSS,DFI,PFI (subset)

suppressPackageStartupMessages({
  library(data.table)
  library(PhenoMapR)
  library(survival)
  library(readxl)
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
cdr_file <- get_flag("--cdr_file", file.path(tpm_dir, "TCGA-CDR-SupplementalTableS1.xlsx"))
out_dir <- get_flag("--out_dir", "results")
# Unshrunk meta-z tables for PRECOG–TCGA Pearson r (forest point sizes); optional, repo root by default
precog_full_path <- get_flag("--precog_full", "PRECOG_V2_Cancer_Type_Meta_Zscores_Final.rds")
tcga_full_path <- get_flag("--tcga_full", "TCGA_metaz.csv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

reference_mat <- PhenoMapR:::get_data(reference)

if (!file.exists(cdr_file)) {
  stop(
    "Missing TCGA-CDR workbook: ", cdr_file,
    "\nPlace TCGA-CDR-SupplementalTableS1.xlsx under data/tcga/ or pass --cdr_file <path>."
  )
}
cat("Loading TCGA-CDR sheet 1:", cdr_file, "\n")
clin <- suppressWarnings(as.data.table(readxl::read_excel(cdr_file, sheet = 1L)))
if (ncol(clin) > 0L && grepl("^\\.\\.\\.", names(clin)[1L])) {
  clin <- clin[, -1L, with = FALSE]
}
need_cdr <- c(
  "bcr_patient_barcode", "type", "OS", "OS.time", "DSS", "DSS.time",
  "DFI", "DFI.time", "PFI", "PFI.time"
)
if (!all(need_cdr %in% names(clin))) {
  stop("TCGA-CDR sheet 1 must include: ", paste(need_cdr, collapse = ", "))
}
data.table::setnames(clin, "type", "acronym", skip_absent = TRUE)

cdr_endpoints_all <- list(
  OS = list(time = "OS.time", event = "OS", label = "OS"),
  DSS = list(time = "DSS.time", event = "DSS", label = "DSS"),
  DFI = list(time = "DFI.time", event = "DFI", label = "DFI"),
  PFI = list(time = "PFI.time", event = "PFI", label = "PFI")
)
endpoints_arg <- get_flag("--endpoints", "OS,DSS,DFI,PFI")
endpoint_names <- trimws(strsplit(endpoints_arg, ",", fixed = TRUE)[[1L]])
endpoint_names <- endpoint_names[nzchar(endpoint_names)]
unknown_ep <- setdiff(endpoint_names, names(cdr_endpoints_all))
if (length(unknown_ep)) {
  stop("Unknown --endpoints labels: ", paste(unknown_ep, collapse = ", "))
}
if (length(endpoint_names) == 0L) stop("No endpoints after parsing --endpoints.")
cdr_endpoints <- cdr_endpoints_all[endpoint_names]
cat("Endpoints:", paste(names(cdr_endpoints), collapse = ", "), "\n")

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

#' Patient-level survival from TCGA-CDR columns (time in days, event 1/0).
build_survival_cdr <- function(cdr_sub, time_col, event_col) {
  if (nrow(cdr_sub) == 0L) {
    return(data.table(
      patient = character(),
      surv_time = numeric(),
      surv_event = integer()
    ))
  }
  tm <- suppressWarnings(as.numeric(cdr_sub[[time_col]]))
  ev <- suppressWarnings(as.integer(cdr_sub[[event_col]]))
  out <- data.table(
    patient = as.character(cdr_sub$bcr_patient_barcode),
    surv_time = tm,
    surv_event = ev
  )
  out[is.finite(surv_time) & surv_time > 0 & !is.na(surv_event)]
}

#' Per-sample panel + full CDR row (cdr_* prefix) for annotation export.
enrich_panel_samples_cdr <- function(panel, cdr_clin, outcomes_basename, endpoint_lab) {
  p <- data.table::copy(panel)
  p[, survival_endpoint := endpoint_lab]
  p[, survival_time_used_days := surv_time]
  p[, survival_event_used := surv_event]

  cm <- data.table::copy(cdr_clin)
  if (!all(c("bcr_patient_barcode", "acronym") %in% names(cm))) {
    stop("CDR table must include bcr_patient_barcode and acronym (type).")
  }
  rest <- setdiff(names(cm), c("bcr_patient_barcode", "acronym"))
  if (length(rest)) {
    nn <- paste0("cdr_", rest)
    if (any(duplicated(nn))) nn <- make.unique(nn, sep = "_")
    data.table::setnames(cm, rest, nn)
  }
  data.table::setnames(cm, "bcr_patient_barcode", "patient")
  data.table::setnames(cm, "acronym", "tcga_code")
  p <- merge(p, cm, by = c("patient", "tcga_code"), all.x = TRUE)
  p[, outcomes_table := outcomes_basename]
  p
}

write_outcomes_note <- function(path, cdr_basename) {
  txt <- c(
    "TCGA PhenoMapR forest analysis — outcomes (TCGA-CDR)",
    "=====================================================",
    "",
    "Survival endpoints come from the TCGA Pan-Cancer Clinical Data Resource (CDR),",
    paste0("Supplemental Table S1, first sheet, file: ", cdr_basename, "."),
    "",
    "Endpoints (as defined in the supplement):",
    "  OS  — OS.time (days), OS (event indicator).",
    "  DSS — DSS.time, DSS.",
    "  DFI — DFI.time, DFI.",
    "  PFI — PFI.time, PFI.",
    "",
    "The script includes patients with finite strictly positive time and non-missing event (0 = censored, 1 = event).",
    "Separate forest tables and figures are written per endpoint (filename suffix _OS, _DSS, _DFI, _PFI).",
    "",
    "Annotated sample tables add CDR covariates with prefix cdr_.",
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
drive_download <- function(file_id, destfile, n_tries = 6L, sleep_sec = 5) {
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

  # Use a more tolerant handle: large files can intermittently reset/timeout.
  handle <- curl::new_handle(
    connecttimeout = 60,
    timeout = 3600,
    low_speed_time = 600,
    low_speed_limit = 1
  )

  # Attempt direct download first; if we get HTML warning, follow confirm flow.
  attempt_once <- function() {
    ok1 <- tryCatch({
      curl::curl_download(url1, destfile = tmp, quiet = TRUE, handle = handle)
      TRUE
    }, error = function(e) {
      message("Drive direct download failed: ", conditionMessage(e))
      FALSE
    })
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
    # Build query string. Values here are simple (ids, uuid, 't'); avoid over-encoding.
    q <- paste0(names(params), "=", unname(params), collapse = "&")
    url2 <- paste0(action, "?", q)

    ok2 <- tryCatch({
      curl::curl_download(url2, destfile = tmp, quiet = TRUE, handle = handle)
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

  for (i in seq_len(max(1L, as.integer(n_tries)))) {
    ok <- attempt_once()
    if (isTRUE(ok)) return(TRUE)
    if (i < n_tries) {
      Sys.sleep(sleep_sec * i)
    }
  }
  FALSE
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
    drive_download(fid, local_path, n_tries = 6L, sleep_sec = 5)
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

count_tcga_genes_scoring <- function(refmat, tpm_genes, tcga_code, z_score_cutoff) {
  if (!tcga_code %in% colnames(refmat)) return(NA_integer_)
  z <- suppressWarnings(as.numeric(as.vector(refmat[, tcga_code, drop = TRUE])))
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
  tcga_ref_mat <- PhenoMapR:::get_data("tcga")
  n_genes_tcga <- count_tcga_genes_scoring(tcga_ref_mat, tpm_genes, tcga_code, z_score_cutoff)

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

  scores_tcga <- PhenoMap(
    expression = mat,
    reference = "tcga",
    cancer_type = tcga_code,
    z_score_cutoff = z_score_cutoff,
    verbose = FALSE
  )
  score_col_tcga <- colnames(scores_tcga)[1L]
  score_tcga_vec <- as.numeric(scores_tcga[[score_col_tcga]])

  list(
    sample_dt = data.table(
      tcga_code = tcga_code,
      stratum = stratum,
      patient = sample_to_patient(sample_ids),
      sample_id = sample_ids,
      sample_type = sample_type,
      precog_label_used = precog_used,
      score = score_used,
      score_tcga = score_tcga_vec
    )[!is.na(patient)],
    n_genes_primary = n_genes_primary,
    n_genes_met = n_genes_met,
    n_genes_tcga = n_genes_tcga
  )
}

#' PhenoMapR TCGA meta-z scores only (no PRECOG reference). For cohorts without a PRECOG label map.
score_tpm_tcga_metaz_only <- function(tpm_file, tcga_code, z_score_cutoff) {
  dt <- fread(tpm_file, showProgress = FALSE, fill = TRUE)
  if (ncol(dt) < 2) {
    stop("Expression file parsed with <2 columns (Gene + samples).")
  }
  cn <- names(dt)
  cn_u <- make.unique(cn, sep = "_dup")
  if (!identical(cn, cn_u)) {
    setnames(dt, cn_u)
  }
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- dt[[1]]
  tpm_genes <- rownames(mat)
  tcga_ref_mat <- PhenoMapR:::get_data("tcga")
  n_genes_tcga <- count_tcga_genes_scoring(tcga_ref_mat, tpm_genes, tcga_code, z_score_cutoff)

  scores_tcga <- PhenoMap(
    expression = mat,
    reference = "tcga",
    cancer_type = tcga_code,
    z_score_cutoff = z_score_cutoff,
    verbose = FALSE
  )
  score_col_tcga <- colnames(scores_tcga)[1L]
  sample_ids <- rownames(scores_tcga)
  sample_type <- sample_to_type_code(sample_ids)
  stratum <- sample_stratum(sample_type)
  score_tcga_vec <- as.numeric(scores_tcga[[score_col_tcga]])

  list(
    sample_dt = data.table(
      tcga_code = tcga_code,
      stratum = stratum,
      patient = sample_to_patient(sample_ids),
      sample_id = sample_ids,
      sample_type = sample_type,
      precog_label_used = NA_character_,
      score = NA_real_,
      score_tcga = score_tcga_vec
    )[!is.na(patient)],
    n_genes_primary = NA_integer_,
    n_genes_met = NA_integer_,
    n_genes_tcga = n_genes_tcga
  )
}

tcga_only <- get_flag("--tcga_only", "")

cat("Fetching Drive folder HTML once for id lookup...\n")
dir.create(tpm_dir, showWarnings = FALSE, recursive = TRUE)
html_tmp <- tempfile(fileext = ".html")
system2("curl", c("-Ls", drive_html_url, "-o", html_tmp), stdout = NULL, stderr = NULL)
html <- paste(readLines(html_tmp, warn = FALSE), collapse = "\n")
unlink(html_tmp)
tpm_id_map <- extract_drive_tpm_id_map(html)
cat("Found TCGA TPM entries in folder HTML:", length(tpm_id_map), "\n")

tcga_codes_with_tpm <- sort(unique(vapply(names(tpm_id_map), extract_tcga_code_from_filename, character(1))))
tcga_types_all <- sort(intersect(unique(as.character(clin$acronym)), tcga_codes_with_tpm))
if (is.character(tcga_only) && nzchar(tcga_only)) {
  only <- strsplit(tcga_only, ",", fixed = TRUE)[[1]]
  tcga_types_all <- intersect(tcga_types_all, only)
  cat("Limiting to tcga_only:", paste(tcga_types_all, collapse = ", "), "\n")
}
tcga_types_precog <- tcga_types_all[vapply(tcga_types_all, function(x) {
  !is.na(concordant_precog_label(x))
}, logical(1))]
cat(
  "TCGA cohorts (CDR \u2229 TPM in Drive):", length(tcga_types_all),
  "| with PRECOG label map:", length(tcga_types_precog), "\n"
)

type_results_median <- stats::setNames(vector("list", length(cdr_endpoints)), names(cdr_endpoints))
type_results_mean <- stats::setNames(vector("list", length(cdr_endpoints)), names(cdr_endpoints))
panel_by_ep <- stats::setNames(vector("list", length(cdr_endpoints)), names(cdr_endpoints))
type_results_median_tcga_meta <- stats::setNames(vector("list", length(cdr_endpoints)), names(cdr_endpoints))
type_results_mean_tcga_meta <- stats::setNames(vector("list", length(cdr_endpoints)), names(cdr_endpoints))
panel_by_ep_tcga_meta <- stats::setNames(vector("list", length(cdr_endpoints)), names(cdr_endpoints))
for (ep in names(cdr_endpoints)) {
  type_results_median[[ep]] <- list()
  type_results_mean[[ep]] <- list()
  panel_by_ep[[ep]] <- list()
  type_results_median_tcga_meta[[ep]] <- list()
  type_results_mean_tcga_meta[[ep]] <- list()
  panel_by_ep_tcga_meta[[ep]] <- list()
}

# Summary forest: OS everywhere except PFI for these TCGA acronyms (CDR endpoints).
pfi_outcome_tcga <- c("PRAD", "DLBC", "LGG", "READ", "TGCT")
# Combined TCGA meta-z forest only: PFI for these types; OS for all others (independent of dual plot rule).
pfi_outcome_tcga_metaz_combined <- c(
  "BRCA", "DLBC", "KICH", "LGG", "PRAD", "READ", "TGCT", "THCA", "THYM"
)
type_results_combined_dual <- list()
type_results_combined_tcga_metaz <- list()

#' Median split on a numeric score vector; Cox HR + log-rank (requires surv_time, surv_event in tb).
cox_median_hr_one <- function(tb, score_vec) {
  if (nrow(tb) < 10L) return(NULL)
  med <- stats::median(score_vec, na.rm = TRUE)
  dt <- data.table::copy(tb)
  dt[, score_split := score_vec]
  dt[, score_grp := factor(fifelse(score_split >= med, "High", "Low"), levels = c("Low", "High"))]
  if (uniqueN(dt$score_grp) < 2L) return(NULL)
  fit <- tryCatch(
    coxph(Surv(surv_time, surv_event) ~ score_grp, data = dt),
    error = function(e) NULL
  )
  if (is.null(fit) || !("score_grpHigh" %in% names(coef(fit)))) return(NULL)
  ch <- coef(fit)[["score_grpHigh"]]
  if (!is.finite(ch)) return(NULL)
  hr <- exp(ch)
  ci <- tryCatch(exp(confint(fit)[1, ]), error = function(e) NULL)
  if (is.null(ci) || length(ci) < 2L) return(NULL)
  lr <- tryCatch(survdiff(Surv(surv_time, surv_event) ~ score_grp, data = dt), error = function(e) NULL)
  if (is.null(lr)) return(NULL)
  p_lr <- suppressWarnings(as.numeric(1 - pchisq(lr$chisq, 1)))
  list(
    hr = hr,
    ci_low = ci[1],
    ci_high = ci[2],
    p_logrank = p_lr,
    median_threshold = med,
    n_samples = nrow(dt),
    n_patients = uniqueN(dt$patient)
  )
}

for (tcga_code in tcga_types_precog) {
  precog_label_primary <- concordant_precog_label(tcga_code)
  precog_label_met <- metastatic_precog_label(tcga_code)
  cat("\n=== Processing TCGA:", tcga_code,
      "PRECOG(primary):", precog_label_primary,
      if (!is.na(precog_label_met)) paste0(" PRECOG(met): ", precog_label_met) else "",
      "===\n", sep = " ")
  tpm_local <- download_one_tpm(tcga_code, precog_label_primary, tpm_id_map)
  if (is.na(tpm_local)) next

  tryCatch(
    {
      scr <- score_tpm_to_sample_scores(
        tpm_file = tpm_local,
        tcga_code = tcga_code,
        precog_label_primary = precog_label_primary,
        precog_label_metastatic = precog_label_met,
        z_score_cutoff = z_score_cutoff,
        reference_mat = reference_mat
      )
      score_dt <- scr$sample_dt

      score_dt_patient <- score_dt[
        is.finite(score) & is.finite(score_tcga) & !is.na(score) & !is.na(score_tcga),
        .(
          score = mean(score, na.rm = TRUE),
          score_tcga = mean(score_tcga, na.rm = TRUE),
          precog_label_used = precog_label_used[1],
          tcga_code = tcga_code[1]
        ),
        by = .(patient, stratum)
      ]

      score_dt_patient_tcga_meta <- score_dt[
        is.finite(score_tcga) & !is.na(score_tcga),
        .(
          score_tcga = mean(score_tcga, na.rm = TRUE),
          precog_label_used = precog_label_used[1],
          tcga_code = tcga_code[1]
        ),
        by = .(patient, stratum)
      ]

      ep_comb <- if (tcga_code %in% pfi_outcome_tcga) {
        list(time = "PFI.time", event = "PFI", label = "PFI")
      } else {
        list(time = "OS.time", event = "OS", label = "OS")
      }
      surv_dt_comb <- build_survival_cdr(clin[acronym == tcga_code], ep_comb$time, ep_comb$event)
      dat_comb <- merge(surv_dt_comb, score_dt_patient, by = "patient", all = FALSE)
      if (nrow(dat_comb) >= 20) {
        for (st in sort(unique(dat_comb$stratum))) {
          sub_c <- dat_comb[stratum == st]
          if (nrow(sub_c) < 20) next
          forest_label_st <- ifelse(st == "Tumor", tcga_code, paste0(tcga_code, "_", st))
          precog_used_st <- if (st == "Metastatic" && !is.na(precog_label_met)) precog_label_met else precog_label_primary
          n_precog_genes_st <- if (st == "Metastatic" && !is.na(precog_label_met) && is.finite(scr$n_genes_met)) {
            as.integer(scr$n_genes_met)
          } else {
            as.integer(scr$n_genes_primary)
          }
          n_tcga_genes_st <- as.integer(scr$n_genes_tcga)

          rp <- cox_median_hr_one(sub_c, sub_c$score)
          rt <- cox_median_hr_one(sub_c, sub_c$score_tcga)
          rows_dual <- list()
          if (!is.null(rp)) {
            rows_dual[[length(rows_dual) + 1L]] <- data.table(
              forest_label = forest_label_st,
              tcga_code = tcga_code,
              stratum = st,
              survival_endpoint = ep_comb$label,
              outcome_rule = ifelse(tcga_code %in% pfi_outcome_tcga, "PFI_cohort", "OS_cohort"),
              precog_label_used = precog_used_st,
              signature_source = "PRECOG",
              n_samples = rp$n_samples,
              n_patients = rp$n_patients,
              n_precog_genes = n_precog_genes_st,
              n_tcga_genes = n_tcga_genes_st,
              median_threshold = rp$median_threshold,
              hr = rp$hr,
              ci_low = rp$ci_low,
              ci_high = rp$ci_high,
              p_logrank = rp$p_logrank
            )
          }
          if (!is.null(rt)) {
            rows_dual[[length(rows_dual) + 1L]] <- data.table(
              forest_label = forest_label_st,
              tcga_code = tcga_code,
              stratum = st,
              survival_endpoint = ep_comb$label,
              outcome_rule = ifelse(tcga_code %in% pfi_outcome_tcga, "PFI_cohort", "OS_cohort"),
              precog_label_used = precog_used_st,
              signature_source = "TCGA",
              n_samples = rt$n_samples,
              n_patients = rt$n_patients,
              n_precog_genes = n_precog_genes_st,
              n_tcga_genes = n_tcga_genes_st,
              median_threshold = rt$median_threshold,
              hr = rt$hr,
              ci_low = rt$ci_low,
              ci_high = rt$ci_high,
              p_logrank = rt$p_logrank
            )
          }
          if (length(rows_dual) > 0L) {
            type_results_combined_dual[[paste(tcga_code, st, sep = "__")]] <- rbindlist(rows_dual, use.names = TRUE, fill = TRUE)
          }
        }
      }

      # Pan-cancer combined TCGA meta-z: one row per TCGA code (primary Tumor only); OS vs PFI from
      # `pfi_outcome_tcga_metaz_combined` (not the dual-plot rule). Merge on TCGA score only.
      ep_tm <- if (tcga_code %in% pfi_outcome_tcga_metaz_combined) {
        list(time = "PFI.time", event = "PFI", label = "PFI")
      } else {
        list(time = "OS.time", event = "OS", label = "OS")
      }
      surv_dt_tm <- build_survival_cdr(clin[acronym == tcga_code], ep_tm$time, ep_tm$event)
      dat_comb_tm <- merge(surv_dt_tm, score_dt_patient_tcga_meta, by = "patient", all = FALSE)
      sub_tumor_tm <- dat_comb_tm[stratum == "Tumor"]
      if (nrow(sub_tumor_tm) >= 10L) {
        rt <- cox_median_hr_one(sub_tumor_tm, sub_tumor_tm$score_tcga)
        if (!is.null(rt)) {
          type_results_combined_tcga_metaz[[paste(tcga_code, "Tumor", sep = "__")]] <- data.table(
            forest_label = tcga_code,
            tcga_code = tcga_code,
            stratum = "Tumor",
            survival_endpoint = ep_tm$label,
            outcome_rule = ifelse(tcga_code %in% pfi_outcome_tcga_metaz_combined, "PFI_cohort", "OS_cohort"),
            precog_label_used = precog_label_primary,
            n_samples = rt$n_samples,
            n_patients = rt$n_patients,
            n_precog_genes = as.integer(scr$n_genes_tcga),
            median_score = rt$median_threshold,
            hr = rt$hr,
            ci_low = rt$ci_low,
            ci_high = rt$ci_high,
            p_logrank = rt$p_logrank
          )
        }
      }

      for (ep in names(cdr_endpoints)) {
        ep_cfg <- cdr_endpoints[[ep]]
        ep_lab <- ep_cfg$label
        surv_dt <- build_survival_cdr(clin[acronym == tcga_code], ep_cfg$time, ep_cfg$event)
        if (nrow(surv_dt) < 20) {
          cat(
            "  Skipping ", ep, " (", ep_lab, "): too few CDR patients for ", tcga_code,
            " (n=", nrow(surv_dt), ")\n",
            sep = ""
          )
          next
        }
        dat <- merge(surv_dt, score_dt_patient, by = "patient", all = FALSE)
        if (nrow(dat) < 20) {
          cat(
            "  Skipping ", ep, " (", ep_lab, "): PRECOG+TCGA score overlap < 20 for ", tcga_code,
            " (n=", nrow(dat), ")\n",
            sep = ""
          )
        }

        if (nrow(dat) >= 20) {
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

          sub_med <- data.table::copy(sub)
          sub_med[, score_grp := factor(ifelse(score >= med, "High", "Low"), levels = c("Low", "High"))]
          fit <- coxph(Surv(surv_time, surv_event) ~ score_grp, data = sub_med)
          coef_hi <- coef(fit)[["score_grpHigh"]]
          hr <- exp(coef_hi)
          ci <- exp(confint(fit)[1, ])
          lr <- survdiff(Surv(surv_time, surv_event) ~ score_grp, data = sub_med)
          p_lr <- 1 - pchisq(lr$chisq, 1)

          type_results_median[[ep]][[paste(tcga_code, st, sep = "__")]] <- data.table(
            forest_label = forest_label_st,
            tcga_code = tcga_code,
            stratum = st,
            survival_endpoint = ep_lab,
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

          sub_mean <- data.table::copy(sub)
          sub_mean[, score_grp := factor(ifelse(score >= mn, "High", "Low"), levels = c("Low", "High"))]
          fit_m <- coxph(Surv(surv_time, surv_event) ~ score_grp, data = sub_mean)
          coef_hi_m <- coef(fit_m)[["score_grpHigh"]]
          hr_m <- exp(coef_hi_m)
          ci_m <- exp(confint(fit_m)[1, ])
          lr_m <- survdiff(Surv(surv_time, surv_event) ~ score_grp, data = sub_mean)
          p_lr_m <- 1 - pchisq(lr_m$chisq, 1)

          type_results_mean[[ep]][[paste(tcga_code, st, sep = "__")]] <- data.table(
            forest_label = forest_label_st,
            tcga_code = tcga_code,
            stratum = st,
            survival_endpoint = ep_lab,
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

          sub_panel <- data.table::copy(sub)
          sub_panel[, forest_label := forest_label_st]
          sub_panel[, score_grp := factor(ifelse(score >= med, "High", "Low"), levels = c("Low", "High"))]
          sub_panel[, median_score := med]
          sub_panel[, score_grp_mean := factor(ifelse(score >= mn, "High", "Low"), levels = c("Low", "High"))]
          sub_panel[, mean_score := mn]
          out_panels[[st]] <- sub_panel
        }

        if (length(out_panels) > 0L) {
          panel_by_ep[[ep]][[tcga_code]] <- rbindlist(out_panels, use.names = TRUE, fill = TRUE)
        }
        }

        dat_tm <- merge(surv_dt, score_dt_patient_tcga_meta, by = "patient", all = FALSE)
        if (nrow(dat_tm) < 20) {
          cat(
            "  Skipping ", ep, " (", ep_lab, "): TCGA meta-z score overlap < 20 for ", tcga_code,
            " (n=", nrow(dat_tm), ")\n",
            sep = ""
          )
        } else {
          out_panels_tm <- list()
          n_tcga_sig_genes <- as.integer(scr$n_genes_tcga)
          for (st in sort(unique(dat_tm$stratum))) {
            sub_tm <- dat_tm[stratum == st]
            if (nrow(sub_tm) < 20) next

            med_t <- median(sub_tm$score_tcga, na.rm = TRUE)
            mn_t <- mean(sub_tm$score_tcga, na.rm = TRUE)
            forest_label_st <- ifelse(st == "Tumor", tcga_code, paste0(tcga_code, "_", st))
            precog_used_st <- if (st == "Metastatic" && !is.na(precog_label_met)) {
              precog_label_met
            } else {
              precog_label_primary
            }

            sub_med_t <- data.table::copy(sub_tm)
            sub_med_t[, score_grp := factor(ifelse(score_tcga >= med_t, "High", "Low"), levels = c("Low", "High"))]
            fit_t <- coxph(Surv(surv_time, surv_event) ~ score_grp, data = sub_med_t)
            coef_hi_t <- coef(fit_t)[["score_grpHigh"]]
            hr_t <- exp(coef_hi_t)
            ci_t <- exp(confint(fit_t)[1, ])
            lr_t <- survdiff(Surv(surv_time, surv_event) ~ score_grp, data = sub_med_t)
            p_lr_t <- 1 - pchisq(lr_t$chisq, 1)

            type_results_median_tcga_meta[[ep]][[paste(tcga_code, st, sep = "__")]] <- data.table(
              forest_label = forest_label_st,
              tcga_code = tcga_code,
              stratum = st,
              survival_endpoint = ep_lab,
              precog_label_used = precog_used_st,
              n_samples = nrow(sub_med_t),
              n_patients = uniqueN(sub_med_t$patient),
              n_precog_genes = n_tcga_sig_genes,
              median_score = med_t,
              hr = hr_t,
              ci_low = ci_t[1],
              ci_high = ci_t[2],
              p_logrank = p_lr_t
            )

            sub_mean_t <- data.table::copy(sub_tm)
            sub_mean_t[, score_grp := factor(ifelse(score_tcga >= mn_t, "High", "Low"), levels = c("Low", "High"))]
            fit_mt <- coxph(Surv(surv_time, surv_event) ~ score_grp, data = sub_mean_t)
            coef_hi_mt <- coef(fit_mt)[["score_grpHigh"]]
            hr_mt <- exp(coef_hi_mt)
            ci_mt <- exp(confint(fit_mt)[1, ])
            lr_mt <- survdiff(Surv(surv_time, surv_event) ~ score_grp, data = sub_mean_t)
            p_lr_mt <- 1 - pchisq(lr_mt$chisq, 1)

            type_results_mean_tcga_meta[[ep]][[paste(tcga_code, st, sep = "__")]] <- data.table(
              forest_label = forest_label_st,
              tcga_code = tcga_code,
              stratum = st,
              survival_endpoint = ep_lab,
              precog_label_used = precog_used_st,
              n_samples = nrow(sub_mean_t),
              n_patients = uniqueN(sub_mean_t$patient),
              n_precog_genes = n_tcga_sig_genes,
              mean_score = mn_t,
              hr = hr_mt,
              ci_low = ci_mt[1],
              ci_high = ci_mt[2],
              p_logrank = p_lr_mt
            )

            sub_panel_t <- data.table::copy(sub_tm)
            sub_panel_t[, forest_label := forest_label_st]
            sub_panel_t[, score := score_tcga]
            sub_panel_t[, score_grp := factor(ifelse(score_tcga >= med_t, "High", "Low"), levels = c("Low", "High"))]
            sub_panel_t[, median_score := med_t]
            sub_panel_t[, score_grp_mean := factor(ifelse(score_tcga >= mn_t, "High", "Low"), levels = c("Low", "High"))]
            sub_panel_t[, mean_score := mn_t]
            out_panels_tm[[st]] <- sub_panel_t
          }

          if (length(out_panels_tm) > 0L) {
            panel_by_ep_tcga_meta[[ep]][[tcga_code]] <- rbindlist(out_panels_tm, use.names = TRUE, fill = TRUE)
          }
        }
      }
    },
    error = function(e) {
      cat("  ERROR processing ", tcga_code, ": ", conditionMessage(e), "\n", sep = "")
    }
  )

  if (file.exists(tpm_local)) {
    cat("  Deleting local TPM:", tpm_local, "\n")
    try(unlink(tpm_local), silent = TRUE)
  }
}

tcga_metaz_only_codes <- setdiff(tcga_types_all, tcga_types_precog)
if (length(tcga_metaz_only_codes) > 0L) {
  cat(
    "\nAdditional cohorts (TCGA meta-z only, no PRECOG map): ",
    paste(tcga_metaz_only_codes, collapse = ", "), "\n",
    sep = ""
  )
}
for (tcga_code in tcga_metaz_only_codes) {
  cat("\n=== Processing TCGA:", tcga_code, " (TCGA meta-z only; no PRECOG label) ===\n", sep = "")
  tpm_local <- download_one_tpm(tcga_code, tcga_code, tpm_id_map)
  if (is.na(tpm_local)) next

  tryCatch(
    {
      scr <- score_tpm_tcga_metaz_only(tpm_local, tcga_code, z_score_cutoff = z_score_cutoff)
      score_dt <- scr$sample_dt
      score_dt_patient_tcga_meta <- score_dt[
        is.finite(score_tcga) & !is.na(score_tcga),
        .(
          score_tcga = mean(score_tcga, na.rm = TRUE),
          precog_label_used = precog_label_used[1],
          tcga_code = tcga_code[1]
        ),
        by = .(patient, stratum)
      ]

      ep_tm <- if (tcga_code %in% pfi_outcome_tcga_metaz_combined) {
        list(time = "PFI.time", event = "PFI", label = "PFI")
      } else {
        list(time = "OS.time", event = "OS", label = "OS")
      }
      surv_dt_tm <- build_survival_cdr(clin[acronym == tcga_code], ep_tm$time, ep_tm$event)
      dat_comb_tm <- merge(surv_dt_tm, score_dt_patient_tcga_meta, by = "patient", all = FALSE)
      sub_tumor_tm <- dat_comb_tm[stratum == "Tumor"]
      if (nrow(sub_tumor_tm) >= 10L) {
        rt <- cox_median_hr_one(sub_tumor_tm, sub_tumor_tm$score_tcga)
        if (!is.null(rt)) {
          type_results_combined_tcga_metaz[[paste(tcga_code, "Tumor", sep = "__")]] <- data.table(
            forest_label = tcga_code,
            tcga_code = tcga_code,
            stratum = "Tumor",
            survival_endpoint = ep_tm$label,
            outcome_rule = ifelse(tcga_code %in% pfi_outcome_tcga_metaz_combined, "PFI_cohort", "OS_cohort"),
            precog_label_used = NA_character_,
            n_samples = rt$n_samples,
            n_patients = rt$n_patients,
            n_precog_genes = as.integer(scr$n_genes_tcga),
            median_score = rt$median_threshold,
            hr = rt$hr,
            ci_low = rt$ci_low,
            ci_high = rt$ci_high,
            p_logrank = rt$p_logrank
          )
        }
      }
    },
    error = function(e) {
      cat("  ERROR processing ", tcga_code, ": ", conditionMessage(e), "\n", sep = "")
    }
  )

  if (file.exists(tpm_local)) {
    cat("  Deleting local TPM:", tpm_local, "\n")
    try(unlink(tpm_local), silent = TRUE)
  }
}

any_results <- FALSE
for (ep in names(cdr_endpoints)) {
  if (length(type_results_median[[ep]]) > 0L || length(type_results_median_tcga_meta[[ep]]) > 0L) {
    any_results <- TRUE
  }
}
if (length(type_results_combined_dual) > 0L || length(type_results_combined_tcga_metaz) > 0L) {
  any_results <- TRUE
}
if (!any_results) {
  stop("No successful TCGA strata for any endpoint. Check downloads, CDR overlap, and --endpoints.")
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

#' Like `finalize_forest_plot_dt`, but keeps rows with non-finite Cox HR/CI (e.g. separation).
#' Caps HR/CI for log-scale forest: (1) non-finite path; (2) any HR or CI > 100 scaled into (0, 100].
#' Sets `plot_hr_ci_capped` TRUE when display values differ from raw Cox (raw stays in _results.tsv).
finalize_forest_plot_dt_tcga_metaz_combined <- function(type_dt) {
  plot_dt <- data.table::copy(type_dt[stratum != "Normal"])
  if (nrow(plot_dt) == 0L) {
    return(plot_dt)
  }

  plot_dt[, plot_hr_ci_capped := FALSE]
  plot_dt[, bad := !is.finite(hr) | !is.finite(ci_low) | !is.finite(ci_high) |
    ci_low <= 0 | ci_high <= ci_low | hr <= 0]

  cap_hr <- 100
  cap_ci_hi <- 100
  floor_ci_lo <- 0.05

  if (plot_dt[, any(bad)]) {
    plot_dt[bad == TRUE, plot_hr_ci_capped := TRUE]
    plot_dt[bad == TRUE, hr_eff := fifelse(is.finite(hr) & hr > 0, pmin(hr, cap_hr), cap_hr * 0.65)]
    plot_dt[bad == TRUE, ci_low_eff := fifelse(
      is.finite(ci_low) & ci_low > 0,
      pmax(ci_low, floor_ci_lo),
      pmax(floor_ci_lo, hr_eff / 200)
    )]
    plot_dt[bad == TRUE, ci_high_eff := fifelse(
      is.finite(ci_high) & ci_high < Inf,
      pmin(ci_high, cap_ci_hi),
      pmin(cap_ci_hi, pmax(hr_eff * 80, ci_low_eff * 25))
    )]
    plot_dt[bad == TRUE & ci_high_eff <= ci_low_eff, ci_high_eff := pmin(cap_ci_hi, ci_low_eff * 10)]
    plot_dt[bad == TRUE & hr_eff <= ci_low_eff, hr_eff := sqrt(ci_low_eff * ci_high_eff)]
    plot_dt[bad == TRUE & hr_eff >= ci_high_eff, hr_eff := sqrt(ci_low_eff * ci_high_eff)]
    plot_dt[bad == TRUE, `:=`(hr = hr_eff, ci_low = ci_low_eff, ci_high = ci_high_eff)]
    plot_dt[, c("hr_eff", "ci_low_eff", "ci_high_eff") := NULL]
  }

  plot_dt[, bad := NULL]
  plot_dt <- plot_dt[
    is.finite(hr) & is.finite(ci_low) & is.finite(ci_high) &
      ci_low > 0 & ci_high > ci_low & ci_high > hr
  ]

  # Finite Cox but HR or CI above 100: scale proportionally so max is 100 (forest x-axis cap).
  plot_dt[, big100 := is.finite(hr) & is.finite(ci_low) & is.finite(ci_high) & pmax(hr, ci_low, ci_high) > 100]
  if (plot_dt[, any(big100)]) {
    plot_dt[big100 == TRUE, plot_hr_ci_capped := TRUE]
    plot_dt[big100 == TRUE, m100 := pmax(hr, ci_high, ci_low)]
    plot_dt[big100 == TRUE, f100 := 100 / m100]
    plot_dt[big100 == TRUE, `:=`(hr = hr * f100, ci_low = ci_low * f100, ci_high = ci_high * f100)]
    plot_dt[big100 == TRUE, ci_high := pmin(ci_high, 100)]
    plot_dt[big100 == TRUE, ci_low := pmax(ci_low, floor_ci_lo)]
    plot_dt[big100 == TRUE, hr := pmin(pmax(hr, ci_low * 1.02), ci_high / 1.02)]
    plot_dt[, c("m100", "f100", "big100") := NULL]
  }

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
  plot_dt[]
}

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

build_tcga_forest_multipanel <- function(
  plot_dt,
  title,
  forest_png,
  color_legend_name = NULL,
  hr_x_label = NULL,
  genes_x_label = NULL,
  simple_forest_points = FALSE
) {
  legend_color_name <- if (is.null(color_legend_name)) "Log-rank Outcome" else color_legend_name
  hr_x_label <- if (is.null(hr_x_label)) {
    "Hazard Ratio\n(High vs. Low PhenoMapR score)"
  } else {
    hr_x_label
  }
  genes_x_label <- if (is.null(genes_x_label)) {
    "# of genes in\nPRECOG signature"
  } else {
    genes_x_label
  }
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

  if (isTRUE(simple_forest_points)) {
    pt_sz <- mean(size_rng)
    p_forest <- ggplot(pd, aes(y = y_ord, x = hr, color = sig)) +
      geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
      geom_errorbar(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.6) +
      geom_point(shape = 16L, size = pt_sz, stroke = 0) +
      scale_x_log10() +
      scale_color_manual(
        values = c("Significant" = "#B2182B", "Not significant" = "#4D4D4D"),
        name = legend_color_name
      ) +
      labs(
        x = hr_x_label,
        y = "TCGA cancer"
      ) +
      theme_minimal(base_size = plot_base_size) +
      guides(
        color = ggplot2::guide_legend(override.aes = list(size = legend_pt, stroke = 0))
      )
  } else {
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
        x = hr_x_label,
        y = "TCGA cancer"
      ) +
      theme_minimal(base_size = plot_base_size) +
      guides(
        color = ggplot2::guide_legend(override.aes = list(size = legend_pt, stroke = 0)),
        shape = ggplot2::guide_legend(override.aes = list(size = legend_pt)),
        size = ggplot2::guide_legend(override.aes = list(shape = 16))
      )
  }

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
    labs(x = genes_x_label, y = NULL) +
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

finalize_dual_forest_plot_dt <- function(dt_long) {
  plot_dt <- dt_long[
    stratum != "Normal" &
      is.finite(hr) & !is.na(hr) &
      is.finite(ci_low) & !is.na(ci_low) &
      is.finite(ci_high) & !is.na(ci_high)
  ]
  plot_dt[, sig := ifelse(is.finite(p_logrank) & !is.na(p_logrank) & p_logrank < 0.05, "Significant", "Not significant")]
  plot_dt$sig <- factor(plot_dt$sig, levels = c("Significant", "Not significant"))
  plot_dt[, precog_tcga_cor := mapply(calc_precog_tcga_cor, precog_label_used, tcga_code)]
  plot_dt[, precog_tcga_cor_abs := abs(precog_tcga_cor)]
  med_cor_abs <- stats::median(plot_dt$precog_tcga_cor_abs, na.rm = TRUE)
  if (!is.finite(med_cor_abs)) med_cor_abs <- 0.5
  plot_dt[, precog_tcga_cor_abs_plot := precog_tcga_cor_abs]
  plot_dt[!is.finite(precog_tcga_cor_abs_plot), precog_tcga_cor_abs_plot := med_cor_abs]

  plot_dt[, cor_shape := fifelse(
    signature_source == "TCGA",
    NA_character_,
    fifelse(is.finite(precog_tcga_cor) & precog_tcga_cor >= 0, "r_pos", "r_neg")
  )]
  plot_dt[, shape_plot := fifelse(signature_source == "TCGA", "TCGA", cor_shape)]
  plot_dt[, shape_plot := factor(shape_plot, levels = c("r_pos", "r_neg", "TCGA"))]

  ord <- plot_dt[signature_source == "PRECOG", .(forest_label, hro = hr)]
  if (nrow(ord) < uniqueN(plot_dt$forest_label)) {
    ord_t <- plot_dt[signature_source == "TCGA", .(forest_label, hro = hr)]
    miss <- setdiff(unique(plot_dt$forest_label), ord$forest_label)
    if (length(miss) > 0L) {
      ord <- rbind(ord, ord_t[forest_label %in% miss], use.names = TRUE, fill = TRUE)
      ord <- unique(ord, by = "forest_label")
    }
  }
  ord <- unique(ord, by = "forest_label")
  plot_dt <- merge(plot_dt, ord, by = "forest_label", all.x = TRUE)
  plot_dt[, y_ord := stats::reorder(forest_label, hro)]
  plot_dt[, z_score_cutoff := z_score_cutoff]
  plot_dt[, reference := reference]
  plot_dt[]
}

#' Rows from `finalize_dual_forest_plot_dt` for one signature; reorder y by that signature's HR.
#' Gene-count panel uses `n_precog_genes` or `n_tcga_genes` via `n_precog_genes` for `build_tcga_forest_multipanel`.
prep_combined_dual_slice_os_style <- function(plot_dual, signature_source_keep) {
  x <- data.table::copy(plot_dual[signature_source == signature_source_keep])
  if (nrow(x) == 0L) return(x)
  x[, y_ord := stats::reorder(forest_label, hr)]
  x[, cor_shape := fifelse(
    is.finite(precog_tcga_cor) & precog_tcga_cor >= 0,
    "r_pos",
    "r_neg"
  )]
  x[, cor_shape := factor(cor_shape, levels = c("r_pos", "r_neg"))]
  if (signature_source_keep == "TCGA") {
    x[, n_precog_genes := n_tcga_genes]
  }
  x[]
}

build_tcga_forest_dual_multipanel <- function(plot_dt, title, forest_png) {
  n_y <- data.table::uniqueN(plot_dt$forest_label)
  size_rng <- c(
    max(1.1, 1.35 + 0.02 * n_y),
    min(7, 2.05 + 0.14 * n_y)
  )
  plot_base_size <- 16
  plot_title_size <- 18
  legend_pt <- 4.5

  # Numeric y + small offset (PRECOG vs TCGA) — works on all ggplot2 versions (no position_dodge height=).
  uy <- levels(plot_dt$y_ord)
  plot_dt[, y_base := match(as.character(forest_label), uy)]
  plot_dt[, y_pos := y_base + fifelse(signature_source == "PRECOG", -0.22, 0.22)]

  p_forest <- ggplot(plot_dt, aes(
    y = y_pos,
    x = hr,
    colour = sig,
    size = precog_tcga_cor_abs_plot,
    group = signature_source
  )) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      linewidth = 0.55,
      width = 0.06
    ) +
    geom_point(
      data = plot_dt[shape_plot == "r_pos"],
      aes(shape = shape_plot),
      stroke = 0
    ) +
    geom_point(
      data = plot_dt[shape_plot == "r_neg"],
      aes(shape = shape_plot),
      stroke = 0.85
    ) +
    geom_point(
      data = plot_dt[shape_plot == "TCGA"],
      aes(shape = shape_plot),
      stroke = 0
    ) +
    scale_y_continuous(breaks = seq_along(uy), labels = uy, expand = expansion(add = 0.35)) +
    scale_x_log10() +
    scale_shape_manual(
      values = c(r_pos = 16L, r_neg = 1L, TCGA = 15L),
      breaks = c("r_pos", "r_neg", "TCGA"),
      labels = c("PRECOG (r \u2265 0)", "PRECOG (r < 0)", "TCGA meta-z"),
      name = "Score"
    ) +
    scale_size_continuous(
      range = size_rng,
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "|cor(PRECOG, TCGA)|"
    ) +
    scale_color_manual(
      values = c("Significant" = "#B2182B", "Not significant" = "#4D4D4D"),
      name = "Log-rank Outcome"
    ) +
    labs(
      x = "Hazard ratio (high vs low, median split)",
      y = NULL
    ) +
    theme_minimal(base_size = plot_base_size) +
    guides(
      color = ggplot2::guide_legend(order = 1L, override.aes = list(size = legend_pt, stroke = 0)),
      shape = ggplot2::guide_legend(order = 2L, override.aes = list(size = legend_pt)),
      size = ggplot2::guide_legend(order = 3L, override.aes = list(shape = 16L))
    )

  pu <- unique(plot_dt[, .(forest_label, y_base, n_patients)])
  # orientation = "y": horizontal bars (x = bar length). Numeric y alone would draw vertical cols at x=n_patients.
  p_counts <- ggplot(pu, aes(x = n_patients, y = y_base)) +
    geom_col(orientation = "y", width = 0.72, fill = "#5C5C5C", alpha = 0.88) +
    scale_y_continuous(breaks = seq_along(uy), labels = uy, expand = expansion(add = 0.35)) +
    labs(x = "# of TCGA\nPatients", y = NULL) +
    theme_minimal(base_size = plot_base_size) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()
    )

  pu_g <- unique(plot_dt[, .(forest_label, y_base, n_precog_genes, n_tcga_genes)])
  pg_long <- rbind(
    pu_g[, .(forest_label, y_base, genes = n_precog_genes, sig_lab = "PRECOG")],
    pu_g[, .(forest_label, y_base, genes = n_tcga_genes, sig_lab = "TCGA")]
  )
  pg_long[, y_pos := y_base + fifelse(sig_lab == "PRECOG", -0.21, 0.21)]
  p_genes <- ggplot(pg_long, aes(x = genes, y = y_pos, fill = sig_lab)) +
    geom_col(orientation = "y", width = 0.32, alpha = 0.9) +
    scale_y_continuous(breaks = seq_along(uy), labels = uy, expand = expansion(add = 0.35)) +
    scale_fill_manual(values = c(PRECOG = "#4A6FA5", TCGA = "#8B6914"), name = "Genes in\nsignature") +
    scale_x_continuous(
      labels = function(x) tolower(scales::label_number(scale_cut = scales::cut_short_scale())(x))
    ) +
    labs(x = "# of genes (|z| cutoff)", y = NULL) +
    theme_minimal(base_size = plot_base_size) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()
    )

  p <- p_forest + p_counts + p_genes +
    patchwork::plot_layout(widths = c(3.4, 1, 1.15), guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = ggplot2::element_text(hjust = 0.5, size = plot_title_size))
    )

  ggplot2::ggsave(forest_png, p, width = 13.2, height = max(4.8, 0.24 * n_y + 2.4), dpi = 200)
  cat("Wrote combined dual-signature forest plot:", forest_png, "\n")
  invisible(NULL)
}

outcomes_basename <- basename(cdr_file)
note_out <- file.path(out_dir, "tcga_outcomes_methods_note.txt")
write_outcomes_note(note_out, outcomes_basename)
cat("Wrote outcomes note (TCGA-CDR):", note_out, "\n")

hr_x_tcga_metaz <- "Hazard Ratio\n(High vs. Low TCGA meta-z score)"
genes_x_tcga_metaz <- "# of genes in\nTCGA meta-z signature"

for (ep in names(cdr_endpoints)) {
  ep_cfg <- cdr_endpoints[[ep]]
  ep_lab <- ep_cfg$label

  if (length(type_results_median[[ep]]) > 0L) {
    type_dt_median <- rbindlist(type_results_median[[ep]], use.names = TRUE, fill = TRUE)
    type_dt_median <- type_dt_median[order(hr, decreasing = TRUE)]
    type_dt_mean <- rbindlist(type_results_mean[[ep]], use.names = TRUE, fill = TRUE)
    type_dt_mean <- type_dt_mean[order(hr, decreasing = TRUE)]

    plot_dt_median <- finalize_forest_plot_dt(type_dt_median)
    plot_dt_mean <- finalize_forest_plot_dt(type_dt_mean)

    type_dt_export_median <- merge_type_with_plot_extras(type_dt_median, plot_dt_median)
    type_dt_export_mean <- merge_type_with_plot_extras(type_dt_mean, plot_dt_mean)

    base <- file.path(out_dir, paste0("tcga_precog_forest_", ep))
    panel_out <- paste0(base, "_panel.tsv")
    type_out_median <- paste0(base, "_results.tsv")
    plot_data_out_median <- paste0(base, "_plot_data.tsv")
    type_out_mean <- paste0(base, "_mean_results.tsv")
    plot_data_out_mean <- paste0(base, "_mean_plot_data.tsv")
    samples_ann_out <- paste0(base, "_samples_annotated.tsv")

    panel_dt_ep <- rbindlist(panel_by_ep[[ep]], use.names = TRUE, fill = TRUE)
    fwrite(panel_dt_ep, panel_out, sep = "\t")
    fwrite(type_dt_export_median, type_out_median, sep = "\t")
    fwrite(plot_dt_median, plot_data_out_median, sep = "\t")
    fwrite(type_dt_export_mean, type_out_mean, sep = "\t")
    fwrite(plot_dt_mean, plot_data_out_mean, sep = "\t")

    panel_enriched <- enrich_panel_samples_cdr(panel_dt_ep, clin, outcomes_basename, ep_lab)
    panel_enriched[, z_score_cutoff := z_score_cutoff]
    panel_enriched[, reference := reference]
    fwrite(panel_enriched, samples_ann_out, sep = "\t")

    cat("\n--- Endpoint ", ep_lab, " (PRECOG reference) ---\n", sep = "")
    cat("Wrote panel:", panel_out, "\n")
    cat("Wrote median-split results:", type_out_median, "\n")
    cat("Wrote median-split plot data:", plot_data_out_median, "\n")
    cat("Wrote mean-split results:", type_out_mean, "\n")
    cat("Wrote mean-split plot data:", plot_data_out_mean, "\n")
    cat("Wrote annotated samples:", samples_ann_out, "\n")

    build_tcga_forest_multipanel(
      plot_dt_median,
      paste0("TCGA ", ep_lab, " — median PhenoMapR PRECOG score"),
      paste0(base, ".png")
    )
    build_tcga_forest_multipanel(
      plot_dt_mean,
      paste0("TCGA ", ep_lab, " — mean-split PhenoMapR PRECOG score"),
      paste0(base, "_mean.png")
    )

    cat("\nMedian split — top results (", ep_lab, ", PRECOG, sorted by HR):\n", sep = "")
    print(type_dt_median[, .(
      forest_label, tcga_code, stratum, survival_endpoint, precog_label_used,
      n_samples, n_patients, n_precog_genes, median_score, hr, ci_low, ci_high, p_logrank
    )])
    cat("\nMean split — top results (", ep_lab, ", PRECOG, sorted by HR):\n", sep = "")
    print(type_dt_mean[, .(
      forest_label, tcga_code, stratum, survival_endpoint, precog_label_used,
      n_samples, n_patients, n_precog_genes, mean_score, hr, ci_low, ci_high, p_logrank
    )])
  } else {
    warning("No PRECOG forest strata for endpoint ", ep, " — skipping PRECOG outputs for this endpoint.")
  }

  if (length(type_results_median_tcga_meta[[ep]]) > 0L) {
    type_tm_med <- rbindlist(type_results_median_tcga_meta[[ep]], use.names = TRUE, fill = TRUE)
    type_tm_med <- type_tm_med[order(hr, decreasing = TRUE)]
    type_tm_mean <- rbindlist(type_results_mean_tcga_meta[[ep]], use.names = TRUE, fill = TRUE)
    type_tm_mean <- type_tm_mean[order(hr, decreasing = TRUE)]

    plot_tm_med <- finalize_forest_plot_dt(type_tm_med)
    plot_tm_mean <- finalize_forest_plot_dt(type_tm_mean)

    type_tm_export_med <- merge_type_with_plot_extras(type_tm_med, plot_tm_med)
    type_tm_export_mean <- merge_type_with_plot_extras(type_tm_mean, plot_tm_mean)
    type_tm_export_med[, reference := "tcga"]
    type_tm_export_mean[, reference := "tcga"]
    plot_tm_med[, reference := "tcga"]
    plot_tm_mean[, reference := "tcga"]

    base_tm <- file.path(out_dir, paste0("tcga_metaz_signature_forest_", ep))
    panel_tm_out <- paste0(base_tm, "_panel.tsv")
    type_tm_out_med <- paste0(base_tm, "_results.tsv")
    plot_tm_out_med <- paste0(base_tm, "_plot_data.tsv")
    type_tm_out_mean <- paste0(base_tm, "_mean_results.tsv")
    plot_tm_out_mean <- paste0(base_tm, "_mean_plot_data.tsv")
    samples_tm_out <- paste0(base_tm, "_samples_annotated.tsv")

    panel_tm_ep <- rbindlist(panel_by_ep_tcga_meta[[ep]], use.names = TRUE, fill = TRUE)
    fwrite(panel_tm_ep, panel_tm_out, sep = "\t")
    fwrite(type_tm_export_med, type_tm_out_med, sep = "\t")
    fwrite(plot_tm_med, plot_tm_out_med, sep = "\t")
    fwrite(type_tm_export_mean, type_tm_out_mean, sep = "\t")
    fwrite(plot_tm_mean, plot_tm_out_mean, sep = "\t")

    panel_tm_enr <- enrich_panel_samples_cdr(panel_tm_ep, clin, outcomes_basename, ep_lab)
    panel_tm_enr[, z_score_cutoff := z_score_cutoff]
    panel_tm_enr[, reference := "tcga"]
    fwrite(panel_tm_enr, samples_tm_out, sep = "\t")

    cat("\n--- Endpoint ", ep_lab, " (TCGA meta-z reference) ---\n", sep = "")
    cat("Wrote panel:", panel_tm_out, "\n")
    cat("Wrote median-split results:", type_tm_out_med, "\n")
    cat("Wrote median-split plot data:", plot_tm_out_med, "\n")
    cat("Wrote mean-split results:", type_tm_out_mean, "\n")
    cat("Wrote mean-split plot data:", plot_tm_out_mean, "\n")
    cat("Wrote annotated samples:", samples_tm_out, "\n")

    build_tcga_forest_multipanel(
      plot_tm_med,
      paste0("TCGA ", ep_lab, " — median PhenoMapR TCGA meta-z signature"),
      paste0(base_tm, ".png"),
      hr_x_label = hr_x_tcga_metaz,
      genes_x_label = genes_x_tcga_metaz
    )
    build_tcga_forest_multipanel(
      plot_tm_mean,
      paste0("TCGA ", ep_lab, " — mean-split PhenoMapR TCGA meta-z signature"),
      paste0(base_tm, "_mean.png"),
      hr_x_label = hr_x_tcga_metaz,
      genes_x_label = genes_x_tcga_metaz
    )

    cat("\nMedian split — top results (", ep_lab, ", TCGA meta-z, sorted by HR):\n", sep = "")
    print(type_tm_med[, .(
      forest_label, tcga_code, stratum, survival_endpoint, precog_label_used,
      n_samples, n_patients, n_precog_genes, median_score, hr, ci_low, ci_high, p_logrank
    )])
    cat("\nMean split — top results (", ep_lab, ", TCGA meta-z, sorted by HR):\n", sep = "")
    print(type_tm_mean[, .(
      forest_label, tcga_code, stratum, survival_endpoint, precog_label_used,
      n_samples, n_patients, n_precog_genes, mean_score, hr, ci_low, ci_high, p_logrank
    )])
  } else {
    warning("No TCGA meta-z forest strata for endpoint ", ep, " — skipping tcga_metaz outputs for this endpoint.")
  }
}

# Combined summary: OS for most cancers; PFI for PRAD, DLBC, LGG, READ, TGCT.
# Two HRs per stratum (PRECOG vs TCGA meta-z median splits).
type_dt_dual <- rbindlist(type_results_combined_dual, use.names = TRUE, fill = TRUE)
if (nrow(type_dt_dual) > 0L) {
  plot_dual <- finalize_dual_forest_plot_dt(type_dt_dual)
  fwrite(
    type_dt_dual,
    file.path(out_dir, "tcga_precog_forest_combined_dual_results.tsv"),
    sep = "\t"
  )
  fwrite(
    plot_dual,
    file.path(out_dir, "tcga_precog_forest_combined_dual_plot_data.tsv"),
    sep = "\t"
  )
  cat("\nWrote combined dual-signature results:", file.path(out_dir, "tcga_precog_forest_combined_dual_results.tsv"), "\n")
  cat("Wrote combined dual-signature plot data:", file.path(out_dir, "tcga_precog_forest_combined_dual_plot_data.tsv"), "\n")
  build_tcga_forest_dual_multipanel(
    plot_dual,
    paste0(
      "TCGA summary — OS (PFI for PRAD, DLBC, LGG, READ, TGCT); ",
      "median-split HR: PRECOG vs TCGA meta-z signature"
    ),
    file.path(out_dir, "tcga_precog_forest_combined_dual.png")
  )

  plot_precog_only <- prep_combined_dual_slice_os_style(plot_dual, "PRECOG")
  if (nrow(plot_precog_only) > 0L) {
    build_tcga_forest_multipanel(
      plot_precog_only,
      paste0(
        "TCGA summary — OS (PFI for PRAD, DLBC, LGG, READ, TGCT); ",
        "median-split HR (PRECOG signature)"
      ),
      file.path(out_dir, "tcga_precog_forest_combined_dual_precog_only.png")
    )
  }
  plot_tcga_only <- prep_combined_dual_slice_os_style(plot_dual, "TCGA")
  if (nrow(plot_tcga_only) > 0L) {
    build_tcga_forest_multipanel(
      plot_tcga_only,
      paste0(
        "TCGA summary — OS (PFI for PRAD, DLBC, LGG, READ, TGCT); ",
        "median-split HR (TCGA meta-z signature)"
      ),
      file.path(out_dir, "tcga_precog_forest_combined_dual_tcga_only.png"),
      hr_x_label = "Hazard Ratio\n(High vs. Low TCGA meta-z score)",
      genes_x_label = "# of genes in\nTCGA signature"
    )
  }
} else {
  warning("No rows for combined OS/PFI dual-signature forest (check CDR overlap and scoring).")
}

type_dt_comb_tm <- rbindlist(type_results_combined_tcga_metaz, use.names = TRUE, fill = TRUE)
if (nrow(type_dt_comb_tm) > 0L) {
  plot_comb_tm <- finalize_forest_plot_dt_tcga_metaz_combined(type_dt_comb_tm)
  plot_comb_tm[, reference := "tcga"]
  comb_tm_res <- file.path(out_dir, "tcga_metaz_signature_forest_combined_results.tsv")
  comb_tm_pd <- file.path(out_dir, "tcga_metaz_signature_forest_combined_plot_data.tsv")
  comb_tm_png <- file.path(out_dir, "tcga_metaz_signature_forest_combined.png")
  fwrite(type_dt_comb_tm, comb_tm_res, sep = "\t")
  fwrite(plot_comb_tm, comb_tm_pd, sep = "\t")
  cat(
    "\nWrote combined TCGA meta-z signature results (",
    uniqueN(type_dt_comb_tm$tcga_code), " TCGA codes, ",
    nrow(type_dt_comb_tm), " Tumor rows; CDR\u2229Drive cohorts available: ",
    length(tcga_types_all), "): ", comb_tm_res, "\n",
    sep = ""
  )
  cat("Wrote combined TCGA meta-z signature plot data:", comb_tm_pd, "\n")
  build_tcga_forest_multipanel(
    plot_comb_tm,
    paste0(
      "TCGA summary — PFI for BRCA, DLBC, KICH, LGG, PRAD, READ, TGCT, THCA, THYM; OS otherwise; ",
      "median-split HR (PhenoMapR TCGA meta-z, Tumor; plot HR/CI capped at 100 when needed)"
    ),
    comb_tm_png,
    hr_x_label = hr_x_tcga_metaz,
    genes_x_label = genes_x_tcga_metaz,
    simple_forest_points = TRUE
  )
} else {
  warning("No rows for combined TCGA meta-z signature forest (check CDR overlap and TCGA scoring).")
}
