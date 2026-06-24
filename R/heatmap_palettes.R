#' List built-in marker heatmap color palette options
#'
#' Returns RColorBrewer ([ColorBrewer2](https://colorbrewer2.org/)) palette
#' names grouped by type and viridis palette names from the
#' \href{https://github.com/sjmgarnier/viridis}{\pkg{viridis}} family.
#'
#' @return A list with elements \code{brewer} (list of \code{sequential},
#'   \code{diverging}, \code{qualitative} character vectors) and
#'   \code{viridis} (character vector of palette names).
#' @export
list_marker_heatmap_color_palettes <- function() {
  brewer <- list(
    sequential = character(0),
    diverging = character(0),
    qualitative = character(0)
  )
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    info <- RColorBrewer::brewer.pal.info
    brewer$sequential <- rownames(info[info$category == "seq", , drop = FALSE])
    brewer$diverging <- rownames(info[info$category == "div", , drop = FALSE])
    brewer$qualitative <- rownames(info[info$category == "qual", , drop = FALSE])
  }
  list(
    brewer = brewer,
    viridis = .viridis_palette_names()
  )
}


#' Default phenotype-group colors for marker heatmaps
#'
#' @return Named character vector (\code{Most Favorable}, \code{Other},
#'   \code{Most Adverse}).
#' @keywords internal
.default_phenotype_palette <- function() {
  c(`Most Adverse` = "#B2182B", Other = "#F7F7F7", `Most Favorable` = "#2166AC")
}


#' Default PhenoMapR score bar colors (low, mid, high)
#' @keywords internal
.default_score_colors <- function() {
  c("#2166AC", "#FFFFFF", "#B2182B")
}


#' ColorBrewer diverging \code{RdGy} palette, 11 classes (standard order)
#' @keywords internal
.rdgy11_brewer <- c(
  "#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF",
  "#E0E0E0", "#BABABA", "#878787", "#4D4D4D", "#1A1A1A"
)


#' Default scaled-expression heatmap colors (high = red, low = black)
#' @keywords internal
.default_expression_colors <- function() {
  rev(.rdgy11_brewer)
}


#' @keywords internal
.viridis_palette_names <- function() {
  c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")
}


#' @keywords internal
.normalize_palette_spec <- function(spec) {
  if (is.null(spec) || identical(spec, "default")) {
    return(list(source = "default"))
  }
  if (is.character(spec) && length(spec) == 1L && !is.na(spec)) {
    if (identical(spec, "default")) {
      return(list(source = "default"))
    }
    if (grepl("^brewer:", spec, ignore.case = TRUE)) {
      return(list(source = "brewer", name = sub("^brewer:", "", spec, ignore.case = TRUE)))
    }
    if (grepl("^viridis:", spec, ignore.case = TRUE)) {
      return(list(source = "viridis", name = sub("^viridis:", "", spec, ignore.case = TRUE)))
    }
    stop("Unrecognized palette spec character: ", spec)
  }
  if (is.list(spec)) {
    src <- spec$source %||% "default"
    if (!is.character(src) || length(src) != 1L) {
      stop("'source' in palette spec must be a single character string.")
    }
    out <- list(source = src)
    if (!is.null(spec$name)) out$name <- as.character(spec$name)[1L]
    if (!is.null(spec$colors)) out$colors <- as.character(spec$colors)
    return(out)
  }
  stop("Palette spec must be NULL, \"default\", a string like \"brewer:RdBu\", or a list.")
}


#' @keywords internal
.brewer_palette_colors <- function(name, n) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Install 'RColorBrewer' to use ColorBrewer palettes.")
  }
  n <- as.integer(n)[1L]
  if (is.na(n) || n < 1L) stop("'n' must be >= 1 for ColorBrewer palettes.")
  info <- RColorBrewer::brewer.pal.info
  if (!name %in% rownames(info)) {
    stop("Unknown ColorBrewer palette: ", name)
  }
  max_n <- info[name, "maxcolors"]
  if (n == 1L) {
    return(RColorBrewer::brewer.pal(max_n, name)[1L])
  }
  if (n <= max_n) {
    return(RColorBrewer::brewer.pal(n, name))
  }
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(max_n, name))(n)
}


#' @keywords internal
.viridis_palette_colors <- function(name, n) {
  n <- as.integer(n)[1L]
  if (is.na(n) || n < 1L) stop("'n' must be >= 1 for viridis palettes.")
  name <- as.character(name)[1L]
  if (!name %in% .viridis_palette_names()) {
    stop("Unknown viridis palette: ", name)
  }
  pkg <- if (requireNamespace("viridisLite", quietly = TRUE)) {
    "viridisLite"
  } else if (requireNamespace("viridis", quietly = TRUE)) {
    "viridis"
  } else {
    stop("Install 'viridisLite' or 'viridis' to use viridis palettes.")
  }
  fn <- get(name, envir = asNamespace(pkg), mode = "function")
  fn(n)
}


#' @keywords internal
.interpolate_diverging_colors <- function(colors, n) {
  colors <- as.character(colors)
  if (length(colors) < 2L) {
    stop("Need at least two colors to interpolate a diverging palette.")
  }
  if (length(colors) == n) return(colors)
  grDevices::colorRampPalette(colors)(n)
}


#' Resolve phenotype-group colors for marker heatmaps
#' @keywords internal
.resolve_phenotype_palette <- function(spec) {
  spec <- .normalize_palette_spec(spec)
  if (identical(spec$source, "default")) {
    return(.default_phenotype_palette())
  }
  if (identical(spec$source, "manual")) {
    cols <- spec$colors
    if (is.null(cols) || length(cols) < 3L) {
      stop("Manual phenotype palette requires three colors (favorable, other, adverse).")
    }
    names(cols) <- NULL
    return(c(
      `Most Favorable` = cols[1L],
      Other = cols[2L],
      `Most Adverse` = cols[3L]
    ))
  }
  if (identical(spec$source, "brewer")) {
    name <- spec$name %||% "RdBu"
    cols <- .brewer_palette_colors(name, 3L)
    return(c(
      `Most Favorable` = cols[1L],
      Other = cols[2L],
      `Most Adverse` = cols[3L]
    ))
  }
  if (identical(spec$source, "viridis")) {
    cols <- .viridis_palette_colors(spec$name %||% "viridis", 3L)
    return(c(
      `Most Favorable` = cols[1L],
      Other = cols[2L],
      `Most Adverse` = cols[3L]
    ))
  }
  stop("Unsupported phenotype palette source: ", spec$source)
}


#' Resolve PhenoMapR score-bar colors (low, mid, high)
#' @keywords internal
.resolve_score_colors <- function(spec) {
  spec <- .normalize_palette_spec(spec)
  if (identical(spec$source, "default")) {
    return(.default_score_colors())
  }
  if (identical(spec$source, "manual")) {
    cols <- spec$colors
    if (is.null(cols) || length(cols) < 2L) {
      stop("Manual score palette requires at least two colors (low and high).")
    }
    if (length(cols) == 2L) {
      return(c(cols[1L], "#FFFFFF", cols[2L]))
    }
    return(cols[seq_len(min(3L, length(cols)))])
  }
  if (identical(spec$source, "brewer")) {
    return(.brewer_palette_colors(spec$name %||% "RdBu", 3L))
  }
  if (identical(spec$source, "viridis")) {
    cols <- .viridis_palette_colors(spec$name %||% "viridis", 3L)
    return(c(cols[1L], "#FFFFFF", cols[3L]))
  }
  stop("Unsupported score palette source: ", spec$source)
}


#' Resolve scaled-expression heatmap colors
#' @keywords internal
.resolve_expression_colors <- function(spec, n = 11L) {
  spec <- .normalize_palette_spec(spec)
  n <- as.integer(n)[1L]
  if (is.na(n) || n < 2L) n <- 11L
  if (identical(spec$source, "default")) {
    cols <- .default_expression_colors()
    if (length(cols) == n) return(cols)
    mid_i <- (length(cols) + 1L) %/% 2L
    return(.interpolate_diverging_colors(cols[c(1L, mid_i, length(cols))], n))
  }
  if (identical(spec$source, "manual")) {
    cols <- spec$colors
    if (is.null(cols) || length(cols) < 2L) {
      stop("Manual expression palette requires at least two colors.")
    }
    if (length(cols) == 2L) {
      cols <- c(cols[1L], grDevices::colorRampPalette(c(cols[1L], "#FFFFFF", cols[2L]))(3L)[2L], cols[2L])
    }
    return(.interpolate_diverging_colors(cols, n))
  }
  if (identical(spec$source, "brewer")) {
    return(.brewer_palette_colors(spec$name %||% "RdGy", n))
  }
  if (identical(spec$source, "viridis")) {
    return(.viridis_palette_colors(spec$name %||% "viridis", n))
  }
  stop("Unsupported expression palette source: ", spec$source)
}


#' Resolve cell-type colors for marker heatmaps
#' @keywords internal
.resolve_celltype_palette_from_spec <- function(spec, levels, override = NULL) {
  levels <- unique(as.character(levels))
  levels <- levels[nzchar(levels)]
  if (length(levels) == 0L) return(character(0))
  if (!is.null(override)) {
    pal <- override[levels]
    pal[is.na(pal)] <- "#BBBBBB"
    return(pal)
  }
  spec <- .normalize_palette_spec(spec)
  n <- length(levels)
  if (identical(spec$source, "default")) {
    return(get_celltype_palette(levels))
  }
  if (identical(spec$source, "manual")) {
    cols <- spec$colors
    if (is.null(cols) || !length(cols)) {
      stop("Manual cell-type palette requires at least one color.")
    }
    cols <- rep(cols, length.out = n)
    return(stats::setNames(cols, levels))
  }
  if (identical(spec$source, "brewer")) {
    cols <- .brewer_palette_colors(spec$name %||% "Set2", max(3L, n))
    return(stats::setNames(cols[seq_len(n)], levels))
  }
  if (identical(spec$source, "viridis")) {
    cols <- .viridis_palette_colors(spec$name %||% "viridis", n)
    return(stats::setNames(cols, levels))
  }
  stop("Unsupported cell-type palette source: ", spec$source)
}


#' Build marker heatmap color scheme list with defaults filled in
#' @keywords internal
.normalize_marker_heatmap_color_schemes <- function(color_schemes) {
  defaults <- list(
    phenotype = "default",
    score = "default",
    expression = "default",
    celltype = "default"
  )
  if (is.null(color_schemes)) return(defaults)
  if (!is.list(color_schemes)) {
    stop("'color_schemes' must be NULL or a named list.")
  }
  for (nm in names(color_schemes)) {
    if (nm %in% names(defaults)) defaults[[nm]] <- color_schemes[[nm]]
  }
  defaults
}


#' PhenoMapR score color function from resolved low/mid/high colors
#' @keywords internal
.phenomap_score_col_fun_from_colors <- function(smin, smax, score_colors) {
  score_colors <- as.character(score_colors)
  if (length(score_colors) == 2L) {
    score_colors <- c(score_colors[1L], "#FFFFFF", score_colors[2L])
  }
  low <- score_colors[1L]
  mid <- if (length(score_colors) >= 3L) score_colors[2L] else "#FFFFFF"
  high <- score_colors[length(score_colors)]
  if (!is.finite(smin) || !is.finite(smax) || smin == smax) {
    smin <- -1
    smax <- 1
  }
  if (smin < 0 && smax > 0) {
    circlize::colorRamp2(c(smin, 0, smax), c(low, mid, high))
  } else if (smax <= 0) {
    circlize::colorRamp2(c(smin, 0), c(low, mid))
  } else if (smin >= 0) {
    circlize::colorRamp2(c(0, smax), c(mid, high))
  } else {
    circlize::colorRamp2(c(smin, smax), c(low, high))
  }
}


#' Scaled-expression color function from resolved palette colors
#' @keywords internal
.scaled_expr_col_fun_from_colors <- function(scale_clip, expr_colors) {
  n <- length(expr_colors)
  if (n < 2L) {
    expr_colors <- .default_expression_colors()
    n <- length(expr_colors)
  }
  breaks <- seq(scale_clip[1], scale_clip[2], length.out = n)
  circlize::colorRamp2(breaks, expr_colors)
}
