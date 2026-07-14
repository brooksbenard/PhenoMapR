# Shared CytoSPACE spatial map helpers (row/col Visium grid from CytoSPACE metadata).
# Sourced by render_spatial_pair_spatial_maps.R and render_spatial_cytospace_location_plots.R.

.spatial_point_size_mode <- function(mode = NULL) {
  mode <- tolower(mode %||% Sys.getenv("SPATIAL_POINT_SIZE", unset = "scaled"))
  if (!mode %in% c("scaled", "uniform")) {
    stop("point_size_mode must be 'scaled' or 'uniform', not: ", mode)
  }
  mode
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.spatial_load_cytospace_seurat <- function(rds_path, gd_id = NULL, env_var = "PHENOMAPR_SPATIAL_RDS_URL") {
  if (!file.exists(rds_path)) {
    url <- Sys.getenv(env_var, unset = "")
    if (nzchar(url)) {
      message("Downloading CytoSPACE RDS from ", env_var, " ...")
      utils::download.file(url, rds_path, mode = "wb", quiet = TRUE)
    } else if (!is.null(gd_id) && requireNamespace("googledrive", quietly = TRUE)) {
      options(googledrive_quiet = TRUE)
      googledrive::drive_deauth()
      message("Downloading CytoSPACE RDS from Google Drive ...")
      googledrive::drive_download(googledrive::as_id(gd_id), rds_path, overwrite = TRUE)
    }
  }
  if (!file.exists(rds_path)) {
    stop("CytoSPACE RDS not found: ", rds_path)
  }
  readRDS(rds_path)
}

.spatial_pick_score_col <- function(meta_names) {
  hit <- grep("weighted_sum_score", meta_names, value = TRUE)
  if (length(hit)) return(hit[1])
  stop(
    "No weighted_sum_score column in CytoSPACE metadata. ",
    "Expected a column matching weighted_sum_score*; found: ",
    paste(head(meta_names, 12), collapse = ", ")
  )
}

.spatial_cytospace_cell_df <- function(seurat, score_col = NULL) {
  if (is.null(score_col)) {
    score_col <- .spatial_pick_score_col(names(seurat@meta.data))
  }
  if (is.na(score_col) || !nzchar(score_col)) {
    stop("No PhenoMapR score column in CytoSPACE metadata")
  }
  if (!all(c("row", "col", "CellType") %in% names(seurat@meta.data))) {
    stop("CytoSPACE metadata must include row, col, and CellType")
  }

  df <- as.data.frame(seurat@meta.data, stringsAsFactors = FALSE)
  df$Cell_id <- rownames(seurat@meta.data)
  df$Cell <- if ("Cell" %in% names(df)) as.character(df$Cell) else df$Cell_id
  df <- df[, c("Cell", "row", "col", "CellType", score_col), drop = FALSE]
  names(df)[5] <- score_col

  tab_loc <- table(paste(df$row, df$col, sep = "_"))
  df$points_per_location <- as.integer(tab_loc[paste(df$row, df$col, sep = "_")])
  df$percentile <- dplyr::percent_rank(df[[score_col]])
  df$prognostic_group <- dplyr::case_when(
    df$percentile < 0.05 ~ "Favorable",
    df$percentile >= 0.95 ~ "Adverse",
    TRUE ~ "Other"
  )
  df$pheno_z <- as.numeric(scale(df[[score_col]]))

  df <- df[
    stats::complete.cases(df[, c("Cell", "row", "col", "CellType", "prognostic_group")]),
    ,
    drop = FALSE
  ]
  df$Cell <- as.character(df$Cell)
  df$row <- as.numeric(df$row)
  df$col <- as.numeric(df$col)
  df$CellType <- as.character(df$CellType)
  df$prognostic_group <- as.character(df$prognostic_group)
  df <- df[!is.na(df$CellType) & nzchar(trimws(df$CellType)), , drop = FALSE]
  df$ct_pg <- paste0(df$CellType, "_", df$prognostic_group)

  list(cell_df = df, score_col = score_col)
}

.spatial_jitter_params <- function(cell_df, point_range = c(0.35, 1.4), uniform_size = 0.55) {
  rng_row <- diff(range(cell_df$row, na.rm = TRUE))
  rng_col <- diff(range(-cell_df$col, na.rm = TRUE))
  list(
    w = max(0.25, if (rng_row > 0) rng_row * 0.025 else 0.15),
    h = max(0.25, if (rng_col > 0) rng_col * 0.025 else 0.15),
    point_range = point_range,
    uniform_size = uniform_size
  )
}

.spatial_order_by_abs <- function(data, metric_col, extra_cols = character(0)) {
  if (!nrow(data)) return(data)
  cols <- unique(c(metric_col, extra_cols))
  cols <- intersect(cols, names(data))
  if (!length(cols)) return(data)
  ord <- do.call(
    order,
    c(
      lapply(cols, function(col) abs(as.numeric(data[[col]]))),
      list(na.last = TRUE)
    )
  )
  data[ord, , drop = FALSE]
}

.spatial_geom_jitter <- function(
    data,
    mapping = NULL,
    order_col = NULL,
    jitter_params,
    point_size_mode = "scaled",
    fixed_size = NULL,
    alpha = 0.85,
    ...) {
  if (!is.null(order_col) && order_col %in% names(data)) {
    data <- .spatial_order_by_abs(data, order_col)
  }
  if (is.null(mapping)) {
    mapping <- ggplot2::aes(x = .data$row, y = -.data$col)
  }

  mode <- .spatial_point_size_mode(point_size_mode)
  dots <- list(...)
  explicit_size <- !is.null(dots$size)

  if (mode == "uniform" && !explicit_size) {
    dots$size <- jitter_params$uniform_size
    mapping <- mapping # size comes from fixed arg, not aes
  }

  do.call(
    ggplot2::geom_jitter,
    c(
      list(
        data = data,
        mapping = mapping,
        alpha = alpha,
        width = jitter_params$w,
        height = jitter_params$h,
        shape = 16
      ),
      dots
    )
  )
}

.spatial_size_scale <- function(
    point_size_mode = "scaled",
    jitter_params,
    guide = "none",
    name = "Cells per spot") {
  if (.spatial_point_size_mode(point_size_mode) == "uniform") {
    return(list())
  }
  list(
    ggplot2::scale_size_continuous(
      range = jitter_params$point_range,
      trans = "reverse",
      guide = guide,
      name = name
    )
  )
}

.spatial_add_size_scale <- function(plot, point_size_mode, jitter_params, guide = "none", name = "Cells per spot") {
  scales <- .spatial_size_scale(point_size_mode, jitter_params, guide = guide, name = name)
  if (length(scales)) plot + scales[[1]] else plot
}

.spatial_size_mapping <- function(point_size_mode, jitter_params) {
  if (.spatial_point_size_mode(point_size_mode) == "uniform") {
    return(list())
  }
  list(size = ggplot2::aes(size = .data$points_per_location))
}

.spatial_map_theme <- function(base_size = 11) {
  ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#444444", size = base_size - 1),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 2),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

.spatial_output_suffix <- function(point_size_mode) {
  if (.spatial_point_size_mode(point_size_mode) == "uniform") "_uniform" else ""
}
