#!/usr/bin/env Rscript
## Patch anatomical + GSE242912 validation caches without full GEO rebuild.
repo_root <- Sys.getenv("PHENOMAPR_REPO", unset = getwd())
out_dir <- file.path(repo_root, "vignettes")
bulk_val_cohorts <- c("GSE115650", "GSE56787", "GSE26852")

bulk_path <- file.path(out_dir, "FSHD_validation_bulk_plot_df.rds")
if (file.exists(bulk_path)) {
  bulk <- readRDS(bulk_path)
  bulk <- bulk[bulk$cohort %in% bulk_val_cohorts, , drop = FALSE]
  bulk$score_scaled <- ave(
    bulk$score,
    bulk$cohort,
    FUN = function(x) as.numeric(scale(x))
  )
  saveRDS(bulk, bulk_path)
  message("Patched bulk validation cache: ", nrow(bulk), " rows (3 cohorts)")
}

summary_path <- file.path(out_dir, "FSHD_validation_summary_table.rds")
if (file.exists(summary_path)) {
  sm <- readRDS(summary_path)
  bulk_val_roles <- c(
    "Bulk validation", "Bulk validation (mechanistic)", "Bulk validation (time course)"
  )
  sm <- sm[
    !(sm$Role %in% bulk_val_roles & !sm$GEO %in% bulk_val_cohorts),
    ,
    drop = FALSE
  ]
  saveRDS(sm, summary_path)
  message("Patched validation summary table: ", nrow(sm), " rows")
}

load_wellstone_comprehensive_df <- function(cache_dir) {
  local_path <- file.path(cache_dir, "FSHD_wellstone_comprehensive_df.rda")
  url <- paste0(
    "https://github.com/fredhutch/Wellstone_BiLateral_Biopsy/raw/master/data/",
    "comprehensive_df.rda"
  )
  if (!file.exists(local_path)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    utils::download.file(url, local_path, mode = "wb", quiet = TRUE)
  }
  env <- new.env()
  load(local_path, envir = env)
  env$comprehensive_df
}

anat_path <- file.path(out_dir, "FSHD_validation_anatomical_df.rds")
if (file.exists(anat_path)) {
  anatomical_df <- readRDS(anat_path)
  anatomical_df$score_scaled <- as.numeric(scale(anatomical_df$score))
  saveRDS(anatomical_df, anat_path)
  message("Patched anatomical cache: ", nrow(anatomical_df), " rows")
}

gse_path <- file.path(out_dir, "FSHD_validation_gse242912_df.rds")
lr_path <- file.path(out_dir, "FSHD_validation_gse242912_lr_df.rds")
if (file.exists(gse_path)) {
  gse242912_df <- readRDS(gse_path)
  comprehensive_df <- load_wellstone_comprehensive_df(out_dir)
  clin_meta <- comprehensive_df[, c(
    "sample_id", "STIR_RATING", "STIR_status", "FAT_FRACTION", "Cumulative Score"
  )]
  names(clin_meta) <- c(
    "sample_id", "STIR_RATING", "STIR_status", "fat_fraction", "path_score"
  )
  gse242912_df <- merge(
    gse242912_df[, setdiff(
      names(gse242912_df),
      c("STIR_RATING", "STIR_status", "fat_fraction", "path_score")
    ), drop = FALSE],
    clin_meta,
    by = "sample_id",
    all.x = TRUE
  )
  saveRDS(gse242912_df, gse_path)

  suffix_lr_cols <- function(df, sfx) {
    nm <- names(df)
    nm[nm != "subject"] <- paste0(nm[nm != "subject"], "_", sfx)
    nm <- sub("^STIR_RATING_", "stir_rating_", nm)
    nm <- sub("^STIR_status_", "stir_status_", nm)
    names(df) <- nm
    df
  }
  lr_cols <- c(
    "subject", "score", "score_scaled",
    "STIR_RATING", "STIR_status", "fat_fraction", "path_score"
  )
  lr_cols <- intersect(lr_cols, names(gse242912_df))
  left <- suffix_lr_cols(gse242912_df[gse242912_df$location == "L", lr_cols, drop = FALSE], "L")
  right <- suffix_lr_cols(gse242912_df[gse242912_df$location == "R", lr_cols, drop = FALSE], "R")
  gse242912_lr_df <- merge(left, right, by = "subject")
  gse242912_lr_df$score_diff <- abs(gse242912_lr_df$score_L - gse242912_lr_df$score_R)
  gse242912_lr_df$score_scaled_diff <- abs(
    gse242912_lr_df$score_scaled_L - gse242912_lr_df$score_scaled_R
  )
  med_diff <- stats::median(gse242912_lr_df$score_diff, na.rm = TRUE)
  gse242912_lr_df$score_discordant <- gse242912_lr_df$score_diff > med_diff
  gse242912_lr_df$stir_discordant <- with(
    gse242912_lr_df,
    !is.na(stir_status_L) & !is.na(stir_status_R) & stir_status_L != stir_status_R
  )
  gse242912_lr_df$score_higher_on_stir_pos <- with(
    gse242912_lr_df,
    ifelse(
      stir_status_L == "STIR+" & stir_status_R == "STIR-",
      score_scaled_L > score_scaled_R,
      ifelse(
        stir_status_R == "STIR+" & stir_status_L == "STIR-",
        score_scaled_R > score_scaled_L,
        NA
      )
    )
  )
  gse242912_lr_df$signed_score_stir_pos_minus_neg <- with(
    gse242912_lr_df,
    ifelse(
      stir_status_L == "STIR+" & stir_status_R == "STIR-",
      score_scaled_L - score_scaled_R,
      ifelse(
        stir_status_R == "STIR+" & stir_status_L == "STIR-",
        score_scaled_R - score_scaled_L,
        NA_real_
      )
    )
  )
  gse242912_lr_df$score_stir_alignment <- with(
    gse242912_lr_df,
    ifelse(
      score_scaled_L > score_scaled_R &
        stir_status_L == "STIR+" & stir_status_R == "STIR-",
      "Higher score on STIR+ side",
      ifelse(
        score_scaled_R > score_scaled_L &
          stir_status_R == "STIR+" & stir_status_L == "STIR-",
        "Higher score on STIR+ side",
        ifelse(
          stir_status_L == stir_status_R,
          paste0("Both ", stir_status_L),
          "STIR discordant; higher score on STIR- side"
        )
      )
    )
  )
  if (all(c("fat_fraction_L", "fat_fraction_R") %in% names(gse242912_lr_df))) {
    gse242912_lr_df$signed_fat_high_minus_low <- with(
      gse242912_lr_df,
      ifelse(
        score_scaled_L > score_scaled_R,
        fat_fraction_L - fat_fraction_R,
        fat_fraction_R - fat_fraction_L
      )
    )
  }
  if (all(c("path_score_L", "path_score_R") %in% names(gse242912_lr_df))) {
    gse242912_lr_df$signed_path_high_minus_low <- with(
      gse242912_lr_df,
      ifelse(
        score_scaled_L > score_scaled_R,
        path_score_L - path_score_R,
        path_score_R - path_score_L
      )
    )
  }
  if (file.exists(lr_path)) {
    old_lr <- readRDS(lr_path)
    old_lr <- old_lr[, !duplicated(names(old_lr)), drop = FALSE]
    clin_cols <- c("dux4_rlogsum_v2", "fat_fraction_v2", "stir_rating_v2", "pathology_v2")
    clin_cols <- intersect(clin_cols, names(old_lr))
    if (length(clin_cols) > 0L) {
      clin_df <- old_lr[, c("subject", clin_cols), drop = FALSE]
      clin_df <- clin_df[, !duplicated(names(clin_df)), drop = FALSE]
      gse242912_lr_df <- merge(gse242912_lr_df, clin_df, by = "subject", all.x = TRUE)
    }
  }
  gse242912_lr_df <- gse242912_lr_df[, !duplicated(names(gse242912_lr_df)), drop = FALSE]
  saveRDS(gse242912_lr_df, lr_path)
  message("Patched GSE242912 cache: ", nrow(gse242912_df), " biopsies; ",
          sum(gse242912_lr_df$stir_discordant, na.rm = TRUE), " STIR-discordant pairs")
}

long_path <- file.path(out_dir, "FSHD_validation_longitudinal_df.rds")
if (file.exists(long_path) && requireNamespace("GEOquery", quietly = TRUE)) {
  extract_char_field <- function(pd, field) {
    ch_cols <- grep("^characteristics", names(pd), value = TRUE)
    if (length(ch_cols) == 0) return(rep(NA_character_, nrow(pd)))
    ch <- apply(pd[, ch_cols, drop = FALSE], 1, function(x) paste(x, collapse = " | "))
    pat <- paste0("(?i)", field, ":\\s*([^|]+)")
    out <- sub(paste0(".*", pat, ".*"), "\\1", ch, perl = TRUE)
    out[!grepl(pat, ch, perl = TRUE)] <- NA_character_
    trimws(out)
  }
  long <- readRDS(long_path)
  common <- long$patient_id
  g115 <- suppressMessages(GEOquery::getGEO("GSE115650", GSEMatrix = TRUE)[[1]])
  pd115 <- Biobase::pData(g115)
  meta115 <- tolower(apply(pd115[, grep("^characteristics", names(pd115), value = TRUE), drop = FALSE], 1, paste, collapse = " | "))
  pd115 <- pd115[grepl("fshd|facioscapulohumeral", meta115), , drop = FALSE]
  pd115$patient_id <- extract_char_field(pd115, "sample_id")
  g140 <- suppressMessages(GEOquery::getGEO("GSE140261", GSEMatrix = TRUE)[[1]])
  pd140 <- Biobase::pData(g140)
  pd140$title <- as.character(pd140$title)
  match_v2_num <- function(field) {
    suppressWarnings(as.numeric(extract_char_field(pd140, field)))[
      match(paste0(common, "b"), pd140$title)
    ]
  }
  mean_v1_num <- function(field) {
    vapply(common, function(pid) {
      idx <- which(pd115$patient_id == pid)
      if (length(idx) == 0L) return(NA_real_)
      mean(
        suppressWarnings(as.numeric(extract_char_field(pd115[idx, , drop = FALSE], field))),
        na.rm = TRUE
      )
    }, numeric(1))
  }
  long$path_score_v1 <- mean_v1_num("path.score")
  long$inflam_v1 <- mean_v1_num("inflam")
  long$active_v1 <- mean_v1_num("active")
  long$stir_v1 <- mean_v1_num("stir")
  long$t1_v1 <- mean_v1_num("t1")
  long$fat_fraction_v1 <- mean_v1_num("fat.fraction")
  long$rnaseq_score_v1 <- mean_v1_num("rnaseq.score")
  long$path_score_v2 <- match_v2_num("pathology.score")
  long$inflam_v2 <- match_v2_num("inflammation")
  long$t1_rating_v2 <- match_v2_num("t1_rating")
  long$fat_fraction_v2 <- match_v2_num("fat_fraction")
  long$stir_rating_v2 <- match_v2_num("stir_rating")
  long$dux4_rlogsum_v2 <- match_v2_num("dux4.rlogsum")
  if (!"score_visit1_scaled" %in% names(long)) {
    long$score_visit1_scaled <- as.numeric(scale(long$score_visit1))
    long$score_visit2_scaled <- as.numeric(scale(long$score_visit2))
  }
  long$score_delta <- long$score_visit2 - long$score_visit1
  long$score_delta_scaled <- long$score_visit2_scaled - long$score_visit1_scaled
  long$delta_fat_fraction <- long$fat_fraction_v2 - long$fat_fraction_v1
  long$delta_stir <- long$stir_rating_v2 - long$stir_v1
  long$delta_path_score <- long$path_score_v2 - long$path_score_v1
  long$delta_inflam <- long$inflam_v2 - long$inflam_v1
  long$delta_t1 <- long$t1_rating_v2 - long$t1_v1
  long$score_increased <- long$score_delta > 0
  saveRDS(long, long_path)
  message("Patched longitudinal cache: ", nrow(long), " paired subjects with clinical metadata")
}
