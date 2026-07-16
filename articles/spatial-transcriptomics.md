# Scoring spatial transcriptomics data with PhenoMapR

### Overview

This vignette demonstrates applying
[PhenoMapR](https://brooksbenard.github.io/PhenoMapR) to spatial
transcriptomics data in two parts. **Part 1** uses **spot-level** 10X
Visium data from an HTAN PDAC sample to show broad tissue regions of
phenotypic interest. **Part 2** uses the same sample after we mapped
paired single-cell data to the spots using
[**CytoSPACE**](https://www.nature.com/articles/s41587-023-01697-9); we
then mirror the single-cell workflow (score by cell type, identify
marker genes, and plot heatmaps) at **cell** resolution and assess
**co-localization** of prognostic cell types.

The sample (`HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test`) is from a pancreatic
cancer dataset available from [HTAN](https://humantumoratlas.org/).

### Download data and vizualize H&E

``` r

suppressPackageStartupMessages({
  library(PhenoMapR)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggchicklet2)
})

knitr::opts_chunk$set(fig.width = 12, out.width = "100%", warning = FALSE)
theme_set(theme_minimal(base_size = 14))

options(googledrive_quiet = TRUE)
googledrive::drive_deauth()

.is_drive_quota_json <- function(path) {
  info <- file.info(path)
  if (is.null(path) || !isTRUE(file.exists(path)) ||
      is.na(info$size) || info$size < 10 || info$size > 5000) {
    return(FALSE)
  }
  txt <- tryCatch(
    paste(readLines(path, n = 30L, warn = FALSE), collapse = "\n"),
    error = function(e) ""
  )
  grepl("downloadQuotaExceeded|The download quota for this file", txt)
}
.is_valid_rds_blob <- function(path) {
  if (is.null(path) || !isTRUE(file.exists(path))) return(FALSE)
  if (.is_drive_quota_json(path)) {
    unlink(path)
    return(FALSE)
  }
  info <- file.info(path)
  if (is.na(info$size) || info$size < 1e5) return(FALSE)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 6L)
  is_gzip <- length(magic) >= 2L && identical(magic[1:2], as.raw(c(0x1f, 0x8b)))
  is_xz <- length(magic) >= 6L && identical(magic[1:6], as.raw(c(0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00)))
  is_rds <- length(magic) >= 2L && identical(rawToChar(magic[1:2]), "RD")
  is_gzip || is_xz || is_rds
}
.resolve_vignette_rds <- function(filename, drive_id, env_var = NULL) {
  candidates <- unique(c(
    filename,
    file.path("vignettes", filename),
    file.path("..", filename),
    file.path("..", "vignettes", filename)
  ))
  hit <- candidates[file.exists(candidates) &
                      vapply(candidates, .is_valid_rds_blob, logical(1))][1]
  if (!is.na(hit)) return(hit)
  dest <- filename
  if (!is.null(env_var) && nzchar(env_var)) {
    url <- Sys.getenv(env_var, unset = "")
    if (nzchar(url)) {
      message("Downloading ", filename, " from ", env_var)
      ok <- tryCatch({
        utils::download.file(url, dest, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e) FALSE)
      if (isTRUE(ok) && isTRUE(.is_valid_rds_blob(dest))) return(dest)
      if (file.exists(dest)) unlink(dest)
    }
  }
  message("Downloading ", filename, " from Google Drive")
  ok <- tryCatch({
    googledrive::drive_download(
      googledrive::as_id(drive_id), dest, overwrite = TRUE
    )
    TRUE
  }, error = function(e) FALSE)
  if (isTRUE(ok) && isTRUE(.is_valid_rds_blob(dest))) return(dest)
  if (file.exists(dest)) unlink(dest)
  stop(
    "Could not obtain ", filename,
    if (!is.null(env_var) && nzchar(env_var)) {
      paste0(". Place it under vignettes/ or set ", env_var, ".")
    } else {
      ". Place it under vignettes/."
    },
    call. = FALSE
  )
}

gd_id_spot <- "1OkIr7ksAWxKVjtdlGqYHMidvHZZsySEE"
rds_spot <- .resolve_vignette_rds(
  "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test.hgnc.rds",
  gd_id_spot,
  "PHENOMAPR_SPATIAL_RDS_URL"
)

seurat_spot <- readRDS(rds_spot)

h_and_e <- SpatialPlot(object = seurat_spot,
  features = NULL,
  image.alpha = 1,
  pt.size.factor = 0) +
  guides(fill = "none") 

print(h_and_e)
```

![](spatial-transcriptomics_files/figure-html/load-setup-1.png)

### Score spots with PhenoMapR

``` r

# Score spots
scores_spot <- PhenoMapR::PhenoMap(
    expression = seurat_spot,
    reference = "precog",
    cancer_type = "Pancreatic",
    assay = "Spatial",
    layer = "data",
    verbose = FALSE
  )

# Add scores to Seurat object
seurat_spot <- PhenoMapR::add_scores_to_seurat(seurat_spot, scores_spot)

score_col_spot <- grep("^PhenoMapR_", names(scores_spot), value = TRUE)[1]
## Z-scaled PhenoMapR score on spots (for hex / metadata plots)
sc_sp <- suppressWarnings(as.numeric(seurat_spot@meta.data[[score_col_spot]]))
seurat_spot$phenomapr_scaled_spot <- as.numeric(scale(sc_sp))

plot_score_distribution(
  seurat_spot$phenomapr_scaled_spot,
  main = "PRECOG Pancreatic score distribution (Visium spots)",
  base_size = 14
)
```

![](spatial-transcriptomics_files/figure-html/score-spots-part1-1.png)

### PhenoMapR score distribution across spots

``` r

p <- SpatialPlot(
  object = seurat_spot,
  features = "PhenoMapR_Pancreatic", image.alpha = 0
)
cols_pheno <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

df <- p$data  # extract data from SpatialPlot

spot_hex_bins <- max(10, min(50L, as.integer(round(sqrt(nrow(df))))))

hex_phenomapr <- ggplot(df, aes(
  x = x,
  y = -y,
  z = scale(PhenoMapR_Pancreatic)
)) +
  stat_summary_hex(fun = mean, bins = spot_hex_bins) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    trans = scales::pseudo_log_trans(sigma = 0.2),
    # limits = c(-4, 4),
    name = "PhenoMapR\nz-score"
  )+
  coord_fixed() +
  theme_void()


(h_and_e | hex_phenomapr) + 
  plot_layout(guides = "collect") & 
  theme(aspect.ratio = 1)
```

![](spatial-transcriptomics_files/figure-html/plot-spot-scores-1.png) We
can see that there are clear morphological associations with where
adverse and favorable signals are being mapped to.

Next, let’s define the most adverse and favorable groups by selecting
the top and bottom 5th percentile of spots based on PhenoMapR scores.

``` r

df <- SpatialPlot(
  object = seurat_spot,
  features = "PhenoMapR_Pancreatic", image.alpha = 0
)$data %>%
  as.data.frame()

groups_spot <- PhenoMapR::define_phenotype_groups(df, score_columns = "PhenoMapR_Pancreatic", percentile = 0.05)

df <- left_join(df, groups_spot, by = c("cell" = "cell_id"))

# Helper function for mode
mode_val <- function(z) {
  z <- z[is.finite(z)]
  if (!length(z)) return(NA_real_)
  as.numeric(names(sort(table(z), decreasing = TRUE)[1L]))
}

hex_phenomapr_groups <- ggplot(df, aes(
  x = x,
  y = -y,
  z = as.numeric(factor(
    phenotype_group_PhenoMapR_Pancreatic,
    levels = c("Most Favorable", "Other", "Most Adverse")
  ))
)) +
  stat_summary_hex(
    aes(
      fill = after_stat(factor(
        value,
        levels = 1:3,
        labels = c("Most Favorable", "Other", "Most Adverse")
      ))
    ),
    bins = 50,
    fun = mode_val,
    colour = NA
  ) +
  scale_fill_manual(
    values = c(
      "Most Favorable" = "#2166AC",
      "Other"            = "#F7F7F7",
      "Most Adverse"     = "#B2182B"
    ),
    drop = FALSE,
    name = "PhenoMapR\nGroup"
  ) +
  coord_fixed() +
  theme_void()


(hex_phenomapr | hex_phenomapr_groups) + 
  plot_layout(guides = "collect") & 
  theme(aspect.ratio = 1)
```

![](spatial-transcriptomics_files/figure-html/unnamed-chunk-1-1.png)
Here, we can still see that the most favorable spots seem to cluster
together, while a subset of the most adverse also tend to co-localize.

## Part 2: CytoSPACE-mapped cells

This sample from HTAN contains a paired single-cell sample from the same
tissue block. In order to leverage the increased resolution of the
single-cell data in a spatial context, we use CytoSPACE to map the
single-cells to their spots.

### Load and score sample with PhenoMapR (CytoSPACE)

Here, we load the pre-processed Seurat object with the cells already
maped.

``` r

# Part 1 spot object is no longer needed; free memory before loading CytoSPACE cells.
for (nm in c("seurat_spot", "scores_spot", "h_and_e", "hex_phenomapr",
             "hex_phenomapr_groups", "df", "groups_spot", "p")) {
  if (exists(nm, inherits = FALSE)) rm(list = nm, inherits = FALSE)
}
gc(verbose = FALSE)
```

    ##            used  (Mb) gc trigger  (Mb)  max used  (Mb)
    ## Ncells  4401432 235.1    7907672 422.4   7907672 422.4
    ## Vcells 11321600  86.4  126673958 966.5 126319018 963.8

``` r

# Load the **CytoSPACE** object (single cells placed on Visium coordinates).
rds_cyto <- .resolve_vignette_rds(
  "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds",
  "1gcOyLriW9bIFNbDuQN6Vi1UsrMGKDxll",
  "PHENOMAPR_SPATIAL_CYTO_RDS_URL"
)

seurat <- readRDS(rds_cyto)

# Score CytoSPACE-mapped cells. Prefer counts: many CytoSPACE exports keep an
# empty Assay5 `data` layer (PhenoMapR falls back to counts automatically).
scores_spatial <- PhenoMapR::PhenoMap(
    expression = seurat,
    reference = "precog",
    cancer_type = "Pancreatic",
    assay = "Spatial",
    layer = "counts",
    verbose = FALSE
  )
  seurat <- PhenoMapR::add_scores_to_seurat(seurat, scores_spatial)
score_col <- grep("^PhenoMapR_", names(scores_spatial), value = TRUE)[1]

plot_score_distribution(
  seurat@meta.data[[score_col]],
  main = "PRECOG Pancreatic score distribution (CytoSPACE cells)",
  base_size = 14
)
```

![](spatial-transcriptomics_files/figure-html/score-spots-1.png)

``` r

# Always build the cell_locations table used by location maps and
# co-localization. Pre-rendered location PNGs must not skip this:
# downstream chunks (scores / tails / coloc) depend on it.
cell_locations <- seurat@meta.data %>%
  as.data.frame() %>%
  dplyr::select("Cell", "row", "col", "CellType", dplyr::all_of(score_col)) %>%
  dplyr::group_by(.data$row, .data$col) %>%
  dplyr::mutate(points_per_location = dplyr::n()) %>%
  dplyr::ungroup()
names(cell_locations)[names(cell_locations) == score_col] <- "PhenoMapR_Pancreatic"
cell_locations$CellType <- as.factor(cell_locations$CellType)
rng_row <- diff(range(cell_locations$row, na.rm = TRUE))
rng_col <- diff(range(-cell_locations$col, na.rm = TRUE))
spatial_jitter_w <- max(0.25, if (rng_row > 0) rng_row * 0.025 else 0.15)
spatial_jitter_h <- max(0.25, if (rng_col > 0) rng_col * 0.025 else 0.15)
spatial_point_range <- c(0.5, 1.6)
```

### Score by cell type

Before looking at any spatial information, we plot the **PhenoMapR**
score by cell type to see cell types most enriched in the adverse and
favorable prognostic groups.

``` r

spatial_celltype_pal <- NULL
spatial_celltype_col <- NULL
meta_names <- names(seurat@meta.data)
celltype_col <- "CellType"

df <- seurat@meta.data
n_meta <- nrow(df)
ann_vec <- NULL
if (!is.null(celltype_col)) {
  raw <- df[[celltype_col]]
  if (is.list(raw)) {
    ann_vec <- vapply(raw, function(x) if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)[1], character(1))
  } else {
    ann_vec <- as.vector(raw)
  }
  if (length(ann_vec) != n_meta) ann_vec <- NULL
}

  df$annotation <- factor(ann_vec, exclude = NULL)

  pal <- PhenoMapR::get_celltype_palette(levels(df$annotation))
print(ggplot(df, aes(
    x = reorder(.data$annotation, .data[[score_col]], FUN = median),
    y = .data[[score_col]],
    fill = .data$annotation
  )) +
        # geom_boxplot(outlier.alpha = 0.5, median.linewidth = 0.5,
        #          outlier.fill = NULL,
        #          outlier.color = NULL, outlier.shape = 21) +
      ggchicklet2::geom_chicklet_boxplot(radius = grid::unit(1, "pt"),
                                     outlier.alpha = 0.5, median.linewidth = 0.5,
                 outlier.fill = NULL,
                 outlier.color = NULL, outlier.shape = 21) +
    scale_fill_manual(values = pal, name = celltype_col) +
    coord_cartesian(ylim = c(-10000, 15000)) +
    guides(fill = guide_legend(ncol = 2)) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(y = "PRECOG Pancreatic score", x = celltype_col, title = "PhenoMapR Score by Cell Type"))
```

![](spatial-transcriptomics_files/figure-html/score-by-annotation-1.png)
As seen in the other pancreatic cancer vignettes, the ductal cell type
is the most associated with the adverse prognostic signal, while the
beta and delta cells are most associated with favorable signal.
Interestingly, plasma cells show quite a wide range of score and contain
some of the most favorably prognostic cells across the sample.

### Where are the different cell types?

CytoSPACE maps single cells onto the Visium tissue grid; every analysis
below uses **`row` / `col`** from the CytoSPACE metadata (and
**`imagerow` / `imagecol`** tissue coordinates for spatial CellChat
distance weighting on the same mapped cells). Point size is **scaled
inversely by cells per spot** when multiple cells share a Visium
location (default), or **uniform** for easier visual comparison —
regenerate both with
`Rscript scripts/render_spatial_cytospace_location_plots.R`.

![](../reference/figures/spatial_cytospace_locations_scaled.png)

![](../reference/figures/spatial_cytospace_locations_uniform.png)

The code below mirrors the location render script when pre-rendered
figures are absent (scaled points only):

``` r

spatial_celltype_pal <- PhenoMapR::get_celltype_palette(levels(cell_locations$CellType))
ct_freq <- sort(table(cell_locations$CellType, useNA = "no"), decreasing = TRUE)
ct_order <- names(ct_freq)
cell_locations$celltype_zorder <- as.numeric(
  factor(as.character(cell_locations$CellType), levels = ct_order)
)
ct_pal <- if (!is.null(spatial_celltype_pal)) {
  spatial_celltype_pal
} else {
  PhenoMapR::get_celltype_palette(levels(cell_locations$CellType))
}

cytospace_loc <- ggplot(
  cell_locations,
  aes(
    x = .data$row, y = -.data$col, color = .data$CellType,
    size = points_per_location, zorder = .data$celltype_zorder
  )
) +
  geom_jitter(alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16) +
  scale_color_manual(values = ct_pal, name = "Cell Type", na.value = "grey90") +
  scale_size_continuous(range = spatial_point_range, trans = "reverse", name = "Cells per spot") +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  coord_fixed(ratio = 0.6) +
  theme_void()

print(cytospace_loc)
```

Interestingly, it appears that several cell types tend to co-localize
with themselves, some of which seem to be localized to the areas where
we saw the greatest adverse and favorable scores mapping to at the spot
level.

#### Where raw PhenoMapR scores are

Spatial map of PhenoMapR score (z-scaled for color gradient). Blue =
more favorable, red = more adverse. When pre-rendered figures are
available (pkgdown), the combined location panel above already includes
this view; the stand-alone score map is shown next.

![](../reference/figures/spatial_cytospace_phenomapr_scaled.png)

``` r

cell_locations <- cell_locations[order(abs(cell_locations$PhenoMapR_Pancreatic)), ]

sc_phenomapr <- ggplot(cell_locations, aes(
  x = .data$row,
  y = -.data$col,
  color = scale(PhenoMapR_Pancreatic),
  size = .data$points_per_location
)) +
  geom_jitter(alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16) +
  scale_size_continuous(range = spatial_point_range, trans = "reverse", name = "Cells per spot") +
  scale_color_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0,
    trans = scales::pseudo_log_trans(sigma = 0.2),
    name = "PhenoMapR\nz-score"
  ) +
  coord_fixed(ratio = 0.6) +
  theme_void()

(cytospace_loc | sc_phenomapr) +
  theme(aspect.ratio = 1)
```

#### Where 5th percentile cells are

In order to make the visualization more apparent, we restrict the
spatial map of prognostic groups to: top 5% (Most Adverse), bottom 5%
(Most Favorable), and the rest (Other).

![](../reference/figures/spatial_cytospace_prognostic_tails_scaled.png)

``` r

cell_locations <- cell_locations %>%
  mutate(percentile = percent_rank(PhenoMapR_Pancreatic)) %>%
  mutate(prognostic_group = case_when(
    percentile < 0.05 ~ "Favorable",
    percentile >= 0.95 ~ "Adverse",
    TRUE ~ "Other"
  ))

df_other <- cell_locations %>%
  dplyr::filter(prognostic_group == "Other")

df_extreme <- cell_locations %>%
  dplyr::filter(prognostic_group != "Other")

sc_phenomapr_5 <- ggplot() +
  geom_jitter(
    data = df_other,
    aes(x = .data$row, y = -.data$col, color = .data$prognostic_group,
        size = points_per_location),
    alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16
  ) +
  geom_jitter(
    data = df_extreme,
    aes(x = .data$row, y = -.data$col, color = .data$prognostic_group,
        size = points_per_location),
    alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16
  ) +
  scale_color_manual(
    values = c(`Adverse` = "#B2182B", Other = "#f7f7f7", `Favorable` = "#2166AC"),
    name = "Prognostic group",
    na.value = "grey90",
    drop = FALSE
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_size_continuous(range = spatial_point_range, trans = "reverse", name = "Cells per spot") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  coord_fixed(ratio = 0.6) +
  theme_void()

(sc_phenomapr | sc_phenomapr_5) +
  theme(aspect.ratio = 1)
```

``` r

# Ensure prognostic_group exists for co-localization even when live
# location / tails plotting was skipped (prerendered figures present).
if (!"prognostic_group" %in% names(cell_locations)) {
  cell_locations <- cell_locations %>%
    mutate(percentile = dplyr::percent_rank(PhenoMapR_Pancreatic)) %>%
    mutate(prognostic_group = dplyr::case_when(
      percentile < 0.05 ~ "Favorable",
      percentile >= 0.95 ~ "Adverse",
      TRUE ~ "Other"
    ))
}
```

Visually, adverse and favorable cells can appear to segregate on the
tissue. The next section tests three **CytoSPACE-grounded** spatial
hypotheses on prognostic **`CellType_PrognosticGroup`** labels:

1.  **Uniquely spatial structure between prognostic cells** — do adverse
    vs favorable groups of the same or different cell types co-occur in
    neighborhoods more than expected by chance? (`nhood_enrichment`
    z-scores; local `cooccur_local` scores)
2.  **PhenoMapR tracks co-localization context** — among cells in a
    reference prognostic group, does PhenoMapR correlate with local
    co-localization to a target group, or with the PhenoMapR scores of
    co-localized target neighbors? (Spearman ρ heatmaps)
3.  **Co-localization enriches for spatial CellChat L-R pairs** — among
    significantly co-localized sender → receiver pairs, is there
    enriched ligand–receptor communication in **spatial CellChat** mode
    using CytoSPACE-mapped Visium coordinates? (integrated
    `dual_spatial_lr` pairs)

We quantify (1)–(2) with
**[`spatialCooccur`](https://github.com/juninamo/spatialCooccur)**
(Inamo *et al.*, [medRxiv
2025](https://doi.org/10.1101/2025.08.05.25332835)) and (3) with
**[CellChat v2](https://github.com/jinworks/CellChat)** spatial
inference on the same CytoSPACE object.

### Co-localization of prognostic cells (CytoSPACE)

We use **CytoSPACE-mapped Visium `row` / `col`** coordinates (same as
the location maps above), assembled as a **`data.frame`** with **`x`**,
**`y`**, and **`CellType_PrognosticGroup`** labels
(e.g. `Ductal_Adverse`). This lets us ask whether, for example,
**adverse ductal cells** co-localize differently than **favorable ductal
cells**, and whether those spatial patterns track PhenoMapR scores and
spatial CellChat L-R support.

Install **`spatialCooccur`** from GitHub if needed:
`remotes::install_github("juninamo/spatialCooccur")`.

**Question 1 — spatial determinants:**
**`spatialCooccur::nhood_enrichment()`** builds a **kNN graph** on
CytoSPACE coordinates, aggregates **co-occurrence** between prognostic
label pairs, and compares observed counts to **permutation** nulls
(z-scores; BH FDR on two-sided *p*). **`cooccur_local()`** scores each
cell by whether its **radius** neighborhood contains both a reference
and target prognostic label; pairwise heatmaps summarize mean non-zero
scores per reference group.

**Question 2 — PhenoMapR vs co-localization:** for each reference group
**G** and target **T**, we compute Spearman ρ between (a) each cell’s
PhenoMapR score and its local **`cooccur_local(G, T)`** score, and (b)
each cell’s PhenoMapR score and the **mean PhenoMapR score of T
neighbors** in the same radius graph.

When the pre-rendered PNGs in `inst/figures/` are present (as in the
pkgdown build), the vignette **shows but does not execute** the analysis
code below and **embeds those static figures** instead. Delete a PNG to
re-run that analysis live on the full CytoSPACE cell set (no
subsampling). Regenerate all figures without knitting:
`Rscript scripts/render_spatial_colocalization_heatmap.R`.

#### Neighborhood enrichment (`nhood_enrichment`) — Question 1

``` r

if (!requireNamespace("spatialCooccur", quietly = TRUE)) {
  stop("Install spatialCooccur: remotes::install_github(\"juninamo/spatialCooccur\")")
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
    !requireNamespace("circlize", quietly = TRUE)) {
  stop("Install ComplexHeatmap and circlize for the co-localization heatmap")
}
if (!exists("cell_locations") ||
    !all(c("Cell", "row", "col", "CellType", "prognostic_group") %in% names(cell_locations))) {
  stop("Expected `cell_locations` with Cell, row, col, CellType, prognostic_group")
}

dfb <- cell_locations[
  stats::complete.cases(cell_locations[, c("Cell", "row", "col", "CellType", "prognostic_group")]),
  ,
  drop = FALSE
]
dfb$Cell <- as.character(dfb$Cell)
dfb$row <- as.numeric(dfb$row)
dfb$col <- as.numeric(dfb$col)
dfb$CellType <- as.character(dfb$CellType)
dfb$prognostic_group <- as.character(dfb$prognostic_group)
dfb$prognostic_group <- dplyr::recode(
  dfb$prognostic_group,
  `Most Adverse` = "Adverse",
  `Most Favorable` = "Favorable"
)
dfb <- dfb[dfb$prognostic_group %in% c("Adverse", "Favorable", "Other"), , drop = FALSE]
dfb <- dfb[!is.na(dfb$CellType) & nzchar(trimws(dfb$CellType)), , drop = FALSE]
dfb$ct_pg <- paste0(dfb$CellType, "_", dfb$prognostic_group)

scoc_df <- data.frame(
  x = dfb$row,
  y = dfb$col,
  ct_pg = dfb$ct_pg,
  stringsAsFactors = FALSE
)
rownames(scoc_df) <- dfb$Cell
scoc_df$cell_type <- scoc_df$ct_pg

xy_mat <- as.matrix(scoc_df[, c("x", "y"), drop = FALSE])
if (any(duplicated(xy_mat) | duplicated(xy_mat, fromLast = TRUE))) {
  rng <- max(
    diff(range(xy_mat[, 1], na.rm = TRUE)),
    diff(range(xy_mat[, 2], na.rm = TRUE)),
    na.rm = TRUE
  )
  eps <- if (is.finite(rng) && rng > 0) rng * 1e-5 else 1e-6
  set.seed(4L)
  scoc_df$x <- scoc_df$x + stats::runif(nrow(scoc_df), 0, eps)
  scoc_df$y <- scoc_df$y + stats::runif(nrow(scoc_df), 0, eps)
}

k_sc <- min(20L, max(5L, nrow(scoc_df) - 1L))
n_perm_sc <- 100L
n_grp <- length(unique(scoc_df$ct_pg))
if (nrow(scoc_df) < (k_sc + 3L) || n_grp < 2L) {
  stop("Not enough cells or label groups for nhood_enrichment")
}

scoc_nhood <- spatialCooccur::nhood_enrichment(
  scoc_df,
  cluster_key = "ct_pg",
  neighbors.k = k_sc,
  connectivity_key = "nn",
  transformation = TRUE,
  n_perms = n_perm_sc,
  seed = 42L,
  n_jobs = 1L
)
if (is.null(scoc_nhood$zscore)) {
  stop("spatialCooccur::nhood_enrichment returned no zscore matrix")
}

zm <- as.matrix(scoc_nhood$zscore)
lab_clean <- function(nm) gsub("^Cluster", "", nm)
rn <- lab_clean(rownames(zm))
cn <- lab_clean(colnames(zm))

suffix_pg <- c("Adverse", "Favorable", "Other")
.parse_ct_pg <- function(lab) {
  for (s in suffix_pg) {
    suff <- paste0("_", s)
    if (nzchar(lab) && endsWith(lab, suff)) {
      ct <- substr(lab, 1L, nchar(lab) - nchar(suff))
      return(c(ct = ct, pg = s))
    }
  }
  c(ct = lab, pg = NA_character_)
}

all_labs <- sort(unique(c(rn, cn)))
meta_labs <- as.data.frame(
  do.call(rbind, lapply(all_labs, function(l) {
    p <- .parse_ct_pg(l)
    data.frame(label = l, CellType = p[["ct"]], Prognostic = p[["pg"]], stringsAsFactors = FALSE)
  })),
  stringsAsFactors = FALSE
)
meta_labs$Prognostic <- factor(meta_labs$Prognostic, levels = suffix_pg)
meta_labs <- meta_labs[!is.na(meta_labs$Prognostic), , drop = FALSE]
ord <- order(meta_labs$CellType, meta_labs$Prognostic)
labs_ord <- meta_labs$label[ord]
labs_ord <- labs_ord[labs_ord %in% rn & labs_ord %in% cn]
if (length(labs_ord) < 2L) {
  stop("Not enough overlapping labels in zscore matrix for heatmap")
}

mat <- zm[match(labs_ord, rn), match(labs_ord, cn), drop = FALSE]
dimnames(mat) <- list(labs_ord, labs_ord)

row_ct <- as.character(meta_labs$CellType[match(rownames(mat), meta_labs$label)])
row_pg <- as.character(meta_labs$Prognostic[match(rownames(mat), meta_labs$label)])
col_ct <- as.character(meta_labs$CellType[match(colnames(mat), meta_labs$label)])
col_pg <- as.character(meta_labs$Prognostic[match(colnames(mat), meta_labs$label)])

pal_pg <- c(Adverse = "#B2182B", Favorable = "#2166AC", Other = "#f7f7f7")
uct <- sort(unique(c(row_ct, col_ct)))
pal_ct <- PhenoMapR::get_celltype_palette(uct)

col_ncells <- rep(0L, ncol(mat))
names(col_ncells) <- colnames(mat)
tab_ct <- table(scoc_df$ct_pg)
hit <- names(col_ncells) %in% names(tab_ct)
col_ncells[hit] <- as.integer(tab_ct[names(col_ncells)[hit]])

.spatial_anno_mm <- 2.5
.spatial_bar_mm <- 10
.spatial_hm_cell_mm <- 2.8

.spatial_hm_wh <- function(mat) {
  list(
    width = grid::unit(ncol(mat) * .spatial_hm_cell_mm, "mm"),
    height = grid::unit(nrow(mat) * .spatial_hm_cell_mm, "mm")
  )
}

.spatial_png_inches <- function(mat, pad_w = 3.5, pad_h = 3) {
  wh <- .spatial_hm_wh(mat)
  top_anno_mm <- .spatial_bar_mm + 2 * .spatial_anno_mm
  left_anno_mm <- 2 * .spatial_anno_mm
  c(
    width = as.numeric(grid::convertWidth(wh$width, "in", valueOnly = TRUE)) +
      as.numeric(grid::convertWidth(grid::unit(left_anno_mm, "mm"), "in", valueOnly = TRUE)) + pad_w,
    height = as.numeric(grid::convertHeight(wh$height, "in", valueOnly = TRUE)) +
      as.numeric(grid::convertHeight(grid::unit(top_anno_mm, "mm"), "in", valueOnly = TRUE)) + pad_h
  )
}

.spatial_bar_axis_at <- function(x) {
  mx <- max(x, 1L, na.rm = TRUE)
  if (mx <= 1000) {
    stats::pretty(c(0, mx), n = 3)
  } else {
    c(0, round(mx / 2), mx)
  }
}

.spatial_make_ha <- function(mat) {
  row_ct_ha <- as.character(meta_labs$CellType[match(rownames(mat), meta_labs$label)])
  row_pg_ha <- as.character(meta_labs$Prognostic[match(rownames(mat), meta_labs$label)])
  col_ct_ha <- as.character(meta_labs$CellType[match(colnames(mat), meta_labs$label)])
  col_pg_ha <- as.character(meta_labs$Prognostic[match(colnames(mat), meta_labs$label)])
  col_nc <- col_ncells[colnames(mat)]
  col_nc[is.na(col_nc)] <- 0L
  anno_u <- grid::unit(.spatial_anno_mm, "mm")
  bar_u <- grid::unit(.spatial_bar_mm, "mm")
  bar_axis_at <- .spatial_bar_axis_at(col_nc)
  list(
    row = ComplexHeatmap::rowAnnotation(
      CellType = row_ct_ha,
      `Prognostic Group` = row_pg_ha,
      col = list(CellType = pal_ct, `Prognostic Group` = pal_pg),
      show_annotation_name = FALSE,
      simple_anno_size = anno_u,
      show_legend = c(FALSE)
    ),
    col = ComplexHeatmap::HeatmapAnnotation(
      `# cells` = ComplexHeatmap::anno_barplot(
        col_nc,
        gp = grid::gpar(fill = "#666666"),
        border = FALSE,
        height = bar_u,
        ylim = c(0, max(col_nc, 1L, na.rm = TRUE)),
        axis_param = list(
          at = bar_axis_at,
          gp = grid::gpar(fontsize = 6),
          side = "left"
        )
      ),
      CellType = col_ct_ha,
      `Prognostic Group` = col_pg_ha,
      col = list(CellType = pal_ct, `Prognostic Group` = pal_pg),
      annotation_name_side = "right",
      annotation_name_gp = grid::gpar(fontsize = 7),
      annotation_height = grid::unit(c(.spatial_bar_mm, .spatial_anno_mm, .spatial_anno_mm), "mm")
    )
  )
}

mat_p <- 2 * stats::pnorm(-abs(mat))
mat_p[!is.finite(mat)] <- NA_real_
mat_p_adj <- mat_p
pv_flat <- as.vector(mat_p)
ok_flat <- is.finite(pv_flat)
mat_p_adj[] <- NA_real_
mat_p_adj[ok_flat] <- stats::p.adjust(pv_flat[ok_flat], method = "BH")

mabs <- suppressWarnings(max(abs(as.numeric(mat)), na.rm = TRUE))
if (!is.finite(mabs) || mabs <= 0) mabs <- 1
col_fun <- circlize::colorRamp2(c(-mabs, 0, mabs), c("#7F312FFF", "#f7f7f7", "#005C55FF"))

ha_nhood <- .spatial_make_ha(mat)
row_ha <- ha_nhood$row
col_ha <- ha_nhood$col
hm_wh_nhood <- .spatial_hm_wh(mat)

ht <- ComplexHeatmap::Heatmap(
  mat,
  name = "Z",
  col = col_fun,
  width = hm_wh_nhood$width,
  height = hm_wh_nhood$height,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  left_annotation = row_ha,
  top_annotation = col_ha,
  row_title = "Reference Cell Type",
  column_title = "Neighborhood Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  heatmap_legend_param = list(title = "Z-score"),
  border = TRUE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    qv <- mat_p_adj[i, j]
    if (!is.finite(qv)) return(invisible(NULL))
    sym <- if (qv < 0.001) "***" else if (qv < 0.01) "**" else if (qv < 0.05) "*" else return(invisible(NULL))
    grid::grid.text(sym, x, y, gp = grid::gpar(fontsize = 9))
  }
)

ComplexHeatmap::draw(
  ht,
  column_title = "Neighborhood Enrichment of Prognostic Cell Types",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

ht_nhood_cluster <- ComplexHeatmap::Heatmap(
  mat,
  name = "Z",
  col = col_fun,
  width = hm_wh_nhood$width,
  height = hm_wh_nhood$height,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  left_annotation = row_ha,
  top_annotation = col_ha,
  row_title = "Reference Cell Type",
  column_title = "Neighborhood Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  heatmap_legend_param = list(title = "Z-score"),
  border = TRUE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    qv <- mat_p_adj[i, j]
    if (!is.finite(qv)) return(invisible(NULL))
    sym <- if (qv < 0.001) "***" else if (qv < 0.01) "**" else if (qv < 0.05) "*" else return(invisible(NULL))
    grid::grid.text(sym, x, y, gp = grid::gpar(fontsize = 9))
  }
)

ComplexHeatmap::draw(
  ht_nhood_cluster,
  column_title = "Neighborhood Enrichment of Prognostic Cell Types (clustered)",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# --- Local co-localization (spatialCooccur::cooccur_local) ---
rad_sc <- max(
  3,
  0.04 * max(
    diff(range(scoc_df$x, na.rm = TRUE)),
    diff(range(scoc_df$y, na.rm = TRUE)),
    na.rm = TRUE
  )
)
local_maxnsteps <- 3L
local_min_cells <- 3L

# cooccur_local() rebuilds neighborhoods on every call; precompute once and reuse
# the same radius + diffusion logic for all reference/target pairs.
coords_local <- as.matrix(scoc_df[, c("x", "y"), drop = FALSE])
neighbors_local <- Seurat::FindNeighbors(coords_local, k.param = k_sc, verbose = FALSE)
adj_local <- neighbors_local$nn
rownames(adj_local) <- colnames(adj_local) <- rownames(scoc_df)
degrees_local <- Matrix::colSums(adj_local) + 1
res_local <- RANN::nn2(
  data = coords_local,
  query = coords_local,
  searchtype = "radius",
  radius = rad_sc,
  k = k_sc
)
label_vec <- as.character(scoc_df$cell_type)
neighbor_labels <- vector("list", nrow(scoc_df))
for (i in seq_len(nrow(scoc_df))) {
  ni <- res_local$nn.idx[i, ]
  ni <- ni[ni > 0 & ni != i]
  neighbor_labels[[i]] <- if (length(ni)) label_vec[ni] else character(0)
}

.spatial_cooccur_local_scores <- function(cluster_x, cluster_y) {
  binary <- vapply(
    neighbor_labels,
    function(nl) {
      if (!length(nl)) return(0)
      as.numeric(any(nl == cluster_x) & any(nl == cluster_y))
    },
    numeric(1)
  )
  if (local_maxnsteps <= 0L) {
    return(binary)
  }
  s_norm <- matrix(binary) / degrees_local
  as.numeric((adj_local %*% s_norm) + s_norm)
}

all_ref_labels <- labs_ord
mat_local <- matrix(
  NA_real_,
  nrow = length(all_ref_labels),
  ncol = length(all_ref_labels),
  dimnames = list(all_ref_labels, all_ref_labels)
)
ref_idx_by_label <- split(seq_len(nrow(scoc_df)), scoc_df$ct_pg)

for (ref_label in all_ref_labels) {
  ref_idx <- ref_idx_by_label[[ref_label]]
  if (is.null(ref_idx) || length(ref_idx) < local_min_cells) next

  for (target_label in all_ref_labels) {
    # Equivalent to:
    # loc <- spatialCooccur::cooccur_local(
    #   scoc_df,
    #   cluster_x = ref_label,
    #   cluster_y = target_label,
    #   connectivity_key = "nn",
    #   neighbors.k = k_sc,
    #   radius = rad_sc,
    #   maxnsteps = local_maxnsteps
    # )
    # scores <- loc[[1]]
    scores <- .spatial_cooccur_local_scores(ref_label, target_label)
    s_ref <- scores[ref_idx]
    s_nz <- s_ref[is.finite(s_ref) & s_ref > 0]
    mat_local[ref_label, target_label] <- if (length(s_nz)) mean(s_nz) else 0
  }
}

if (all(!is.finite(mat_local))) {
  stop("spatialCooccur::cooccur_local returned no scores for any reference group")
}

mabs_local <- suppressWarnings(max(as.numeric(mat_local), na.rm = TRUE))
if (!is.finite(mabs_local) || mabs_local <= 0) mabs_local <- 1
col_fun_local <- circlize::colorRamp2(c(0, mabs_local), c("#f7f7f7", "#005C55FF"))

ha_local <- .spatial_make_ha(mat_local)
row_ha_local <- ha_local$row
col_ha_local <- ha_local$col
hm_wh_local <- .spatial_hm_wh(mat_local)

ht_local <- ComplexHeatmap::Heatmap(
  mat_local,
  name = "Mean score",
  col = col_fun_local,
  width = hm_wh_local$width,
  height = hm_wh_local$height,
  na_col = "#eeeeee",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  left_annotation = row_ha_local,
  top_annotation = col_ha_local,
  row_title = "Reference Cell Type",
  column_title = "Target Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  border = TRUE,
  heatmap_legend_param = list(title = "Mean non-zero\nlocal score")
)

ComplexHeatmap::draw(
  ht_local,
  column_title = "Local Co-localization of Prognostic Cell Types",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

mat_local_cluster <- mat_local
mat_local_cluster[!is.finite(mat_local_cluster)] <- 0

ht_local_cluster <- ComplexHeatmap::Heatmap(
  mat_local_cluster,
  name = "Mean score",
  col = col_fun_local,
  width = hm_wh_local$width,
  height = hm_wh_local$height,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  left_annotation = row_ha_local,
  top_annotation = col_ha_local,
  row_title = "Reference Cell Type",
  column_title = "Target Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  border = TRUE,
  heatmap_legend_param = list(title = "Mean non-zero\nlocal score")
)
grDevices::pdf(NULL)
ht_local_cluster_drawn <- ComplexHeatmap::draw(ht_local_cluster)
local_row_ord <- ComplexHeatmap::row_order(ht_local_cluster_drawn)
local_col_ord <- ComplexHeatmap::column_order(ht_local_cluster_drawn)
grDevices::dev.off()

local_outline_threshold <- 0.75
ht_local_cluster_outline <- ComplexHeatmap::Heatmap(
  mat_local_cluster,
  name = "Mean score",
  col = col_fun_local,
  width = hm_wh_local$width,
  height = hm_wh_local$height,
  row_order = local_row_ord,
  column_order = local_col_ord,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  left_annotation = row_ha_local,
  top_annotation = col_ha_local,
  row_title = "Reference Cell Type",
  column_title = "Target Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  border = TRUE,
  heatmap_legend_param = list(title = "Mean non-zero\nlocal score"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    val <- mat_local_cluster[i, j]
    if (is.finite(val) && val >= local_outline_threshold) {
      grid::grid.rect(
        x, y, w, h,
        gp = grid::gpar(col = "black", lwd = 1.5, fill = NA)
      )
    }
  }
)

ComplexHeatmap::draw(
  ht_local_cluster_drawn,
  column_title = "Local Co-localization of Prognostic Cell Types (clustered)",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

ComplexHeatmap::draw(
  ht_local_cluster_outline,
  column_title = paste0(
    "Local Co-localization of Prognostic Cell Types (clustered; score \u2265 ",
    local_outline_threshold, " outlined)"
  ),
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# --- PhenoMapR score vs local co-localization ---
cor_min_cells <- 10L
pheno_vec <- stats::setNames(as.numeric(dfb[[score_col]]), dfb$Cell)
pheno_vec <- pheno_vec[rownames(scoc_df)]

neighbor_idx <- vector("list", nrow(scoc_df))
for (i in seq_len(nrow(scoc_df))) {
  ni <- res_local$nn.idx[i, ]
  neighbor_idx[[i]] <- ni[ni > 0 & ni != i]
}

mat_nbr_pheno <- matrix(
  NA_real_,
  nrow = nrow(scoc_df),
  ncol = length(all_ref_labels),
  dimnames = list(rownames(scoc_df), all_ref_labels)
)
for (i in seq_len(nrow(scoc_df))) {
  ni <- neighbor_idx[[i]]
  if (!length(ni)) next
  nl <- label_vec[ni]
  pv <- pheno_vec[ni]
  for (target_label in unique(nl)) {
    hit <- nl == target_label
    if (any(hit)) {
      mat_nbr_pheno[i, target_label] <- mean(pv[hit], na.rm = TRUE)
    }
  }
}

mat_cor_pheno_cooccur <- matrix(
  NA_real_,
  nrow = length(all_ref_labels),
  ncol = length(all_ref_labels),
  dimnames = list(all_ref_labels, all_ref_labels)
)
mat_cor_pheno_neighbor <- matrix(
  NA_real_,
  nrow = length(all_ref_labels),
  ncol = length(all_ref_labels),
  dimnames = list(all_ref_labels, all_ref_labels)
)

.spatial_spearman_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < cor_min_cells) return(NA_real_)
  x_ok <- x[ok]
  y_ok <- y[ok]
  if (stats::sd(x_ok) == 0 || stats::sd(y_ok) == 0) return(NA_real_)
  stats::cor(x_ok, y_ok, method = "spearman")
}

for (ref_label in all_ref_labels) {
  ref_idx <- ref_idx_by_label[[ref_label]]
  if (is.null(ref_idx) || length(ref_idx) < cor_min_cells) next
  pheno_ref <- pheno_vec[ref_idx]
  for (target_label in all_ref_labels) {
    cooc_ref <- .spatial_cooccur_local_scores(ref_label, target_label)[ref_idx]
    mat_cor_pheno_cooccur[ref_label, target_label] <- .spatial_spearman_cor(pheno_ref, cooc_ref)
    nbr_ref <- mat_nbr_pheno[ref_idx, target_label]
    mat_cor_pheno_neighbor[ref_label, target_label] <- .spatial_spearman_cor(pheno_ref, nbr_ref)
  }
}

.spatial_spearman_palette <- function() {
  c(
    "#283C93", "#4960C1", "#6987D6", "#91ABE4", "#C5D5F1", "#FDF9FA",
    "#E9A5B8", "#D5708F", "#BC5270", "#A83E5D", "#8D2842"
  )
}

.spatial_spearman_col_fun <- function() {
  cols <- .spatial_spearman_palette()
  breaks <- seq(-1, 1, length.out = length(cols))
  circlize::colorRamp2(breaks, cols)
}

col_fun_cor <- .spatial_spearman_col_fun()

.spatial_cor_dist <- function(x) {
  x[!is.finite(x)] <- 0
  stats::dist(x)
}

.draw_cor_heatmap <- function(mat_cor, title, clustered = FALSE, outline_abs_threshold = NULL) {
  ha_cor <- .spatial_make_ha(mat_cor)
  hm_wh <- .spatial_hm_wh(mat_cor)
  row_order <- NULL
  col_order <- NULL
  if (isTRUE(clustered) && !is.null(outline_abs_threshold)) {
    ht_cluster <- ComplexHeatmap::Heatmap(
      mat_cor,
      name = "Spearman",
      col = col_fun_cor,
      width = hm_wh$width,
      height = hm_wh$height,
      na_col = "#eeeeee",
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      clustering_distance_rows = .spatial_cor_dist,
      clustering_distance_columns = .spatial_cor_dist,
      show_row_names = FALSE,
      show_column_names = FALSE,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      left_annotation = ha_cor$row,
      top_annotation = ha_cor$col,
      border = TRUE
    )
    grDevices::pdf(NULL)
    ht_cluster_drawn <- ComplexHeatmap::draw(ht_cluster)
    row_order <- ComplexHeatmap::row_order(ht_cluster_drawn)
    col_order <- ComplexHeatmap::column_order(ht_cluster_drawn)
    grDevices::dev.off()
    clustered <- FALSE
  }
  cell_fun_outline <- NULL
  if (!is.null(outline_abs_threshold)) {
    thr <- outline_abs_threshold
    cell_fun_outline <- function(j, i, x, y, w, h, fill) {
      val <- mat_cor[i, j]
      if (is.finite(val) && abs(val) > thr) {
        grid::grid.rect(
          x, y, w, h,
          gp = grid::gpar(col = "black", lwd = 1.5, fill = NA)
        )
      }
    }
  }
  ComplexHeatmap::Heatmap(
    mat_cor,
    name = "Spearman",
    col = col_fun_cor,
    width = hm_wh$width,
    height = hm_wh$height,
    na_col = "#eeeeee",
    cluster_rows = clustered,
    cluster_columns = clustered,
    row_order = row_order,
    column_order = col_order,
    clustering_distance_rows = if (clustered) .spatial_cor_dist else "euclidean",
    clustering_distance_columns = if (clustered) .spatial_cor_dist else "euclidean",
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 7),
    left_annotation = ha_cor$row,
    top_annotation = ha_cor$col,
    row_title = "Reference Cell Type",
    column_title = "Target Cell Type",
    row_title_gp = grid::gpar(fontsize = 12),
    column_title_gp = grid::gpar(fontsize = 12),
    row_title_side = "left",
    column_title_side = "bottom",
    border = TRUE,
    heatmap_legend_param = list(title = "Spearman rho"),
    cell_fun = cell_fun_outline
  ) -> ht_cor
  ComplexHeatmap::draw(
    ht_cor,
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
}

cor_outline_threshold <- 0.5

.draw_cor_heatmap(
  mat_cor_pheno_cooccur,
  "PhenoMapR Score vs Local Co-localization (Spearman rho)",
  clustered = FALSE
)
.draw_cor_heatmap(
  mat_cor_pheno_cooccur,
  "PhenoMapR Score vs Local Co-localization (Spearman rho, clustered)",
  clustered = TRUE
)
.draw_cor_heatmap(
  mat_cor_pheno_cooccur,
  paste0(
    "PhenoMapR Score vs Local Co-localization (Spearman rho, clustered; |rho| > ",
    cor_outline_threshold, " outlined)"
  ),
  clustered = TRUE,
  outline_abs_threshold = cor_outline_threshold
)
.draw_cor_heatmap(
  mat_cor_pheno_neighbor,
  "PhenoMapR Score vs Co-occurring Neighbor PhenoMapR (Spearman rho)",
  clustered = FALSE
)
.draw_cor_heatmap(
  mat_cor_pheno_neighbor,
  "PhenoMapR Score vs Co-occurring Neighbor PhenoMapR (Spearman rho, clustered)",
  clustered = TRUE
)
.draw_cor_heatmap(
  mat_cor_pheno_neighbor,
  paste0(
    "PhenoMapR Score vs Co-occurring Neighbor PhenoMapR (Spearman rho, clustered; |rho| > ",
    cor_outline_threshold, " outlined)"
  ),
  clustered = TRUE,
  outline_abs_threshold = cor_outline_threshold
)
```

![](../reference/figures/spatial_colocalization_nhood_enrichment.png)

Clustering rows and columns (without dendrograms) reorders the
neighborhood enrichment matrix to highlight groups with similar
co-occurrence profiles.

![](../reference/figures/spatial_colocalization_nhood_enrichment_clustered.png)

#### Local co-localization scores (`cooccur_local`) — Question 1

The heatmap below shows **pairwise local co-localization scores** across
all **`CellType_PrognosticGroup`** labels (not neighborhood enrichment
z-scores). Each cell is the **mean non-zero `cooccur_local()` score**
for reference group (row) versus target group (column).

![](../reference/figures/spatial_colocalization_colocalization_scores.png)

![](../reference/figures/spatial_colocalization_colocalization_scores_clustered.png)

Cells with mean non-zero local scores **≥ 0.75** are outlined in black
on the clustered heatmap below.

![](../reference/figures/spatial_colocalization_colocalization_scores_clustered_outlined.png)

#### PhenoMapR score vs local co-localization — Question 2

We next ask whether a cell’s **PhenoMapR score** is associated with its
**local co-localization** context. For each reference group **G** (rows)
and target group **T** (columns):

1.  **PhenoMapR vs co-localization score:** among cells in **G**,
    Spearman correlation between each cell’s PhenoMapR score and its
    **`cooccur_local(G, T)`** score (same radius and diffusion as
    above).
2.  **PhenoMapR vs neighbor PhenoMapR:** among cells in **G** with at
    least one **T** neighbor in the radius, Spearman correlation between
    the cell’s PhenoMapR score and the **mean PhenoMapR score of T
    neighbors** (e.g., for adverse ductal cells, is co-localization with
    favorable plasma cells associated with both the ductal cell’s score
    and the plasma neighbors’ scores?).

![](../reference/figures/spatial_colocalization_pheno_vs_cooccur_spearman.png)

![](../reference/figures/spatial_colocalization_pheno_vs_cooccur_spearman_clustered.png)

Cells with **\|Spearman rho\| \> 0.5** are outlined in **black** on the
clustered heatmaps below (from
`scripts/render_spatial_colocalization_heatmap.R`). After integrating
**spatial CellChat**, dual-positive sender → receiver pairs receive an
additional **purple** outline on separate figures
(`*_clustered_outlined_integrated.png`; see CellChat section below).

![](../reference/figures/spatial_colocalization_pheno_vs_cooccur_spearman_clustered_outlined.png)

![](../reference/figures/spatial_colocalization_pheno_vs_neighbor_pheno_spearman.png)

Clustered versions reorder row and column annotations with the heatmap
(no dendrograms). Spearman **rho** uses a custom blue–white–pink
diverging palette (ρ = −1 navy blue, ρ = 0 off-white, ρ = +1 burgundy
pink).

![](../reference/figures/spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered.png)

![](../reference/figures/spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered_outlined.png)

### Co-localization meets cell-cell signaling (CellChat) — Question 3

The co-localization and PhenoMapR correlation analyses above identify
**which prognostic cell groups neighbor one another** (Question 1) and
whether **PhenoMapR scores track local mixing** (Question 2). Question 3
asks whether those **significantly co-localized** sender → receiver
pairs are enriched for **spatial CellChat** ligand–receptor (L-R)
interactions on the same CytoSPACE-mapped cells.

We integrate the spatial metrics with **[CellChat
v2](https://github.com/jinworks/CellChat)** ([Jin *et al.*, *Nat.
Protoc.* 2024](https://doi.org/10.1038/s41596-024-01045-4)) in **spatial
mode**, using CytoSPACE-mapped **Visium tissue coordinates**
(`GetTissueCoordinates`) and distance constraints from the CellChat
spatial FAQ (65 µm spot size, `contact.range = 100` µm for Visium,
`interaction.range = 250` µm for secreted signaling):

1.  **Pairwise spatial summary** — for every `CellType_PrognosticGroup`
    sender → receiver pair, combine **neighborhood enrichment** z-score,
    **local co-localization** score (`cooccur_local`), and Spearman
    **ρ** (PhenoMapR vs co-localization / neighbor PhenoMapR).
2.  **Spatial CellChat inference** — run on cells in the **Adverse** and
    **Favorable** tails (≥10 cells per group), using normalized
    **Spatial** assay expression, CytoSPACE **imagerow/imagecol**
    coordinates, and the same combined labels (`Ductal_Adverse`,
    `Plasma_Favorable`, …). Communication probabilities use
    `datatype = "spatial"` with `distance.use = TRUE` so L-R pairs are
    weighted by physical cell–cell distance on tissue.
3.  **Interaction mode** — split CellChat probabilities into
    **physical** (Cell-Cell Contact + ECM-Receptor) vs **secreted**
    (Secreted Signaling) using `CellChatDB` annotations.
4.  **Integrated filters** — flag pairs with neighborhood enrichment
    (*q* \< 0.05), local co-localization (score ≥ 0.5), and/or spatial
    CellChat L-R support; test whether physical vs secreted modes
    associate differently with spatial proximity and PhenoMapR
    correlations.

Install CellChat if needed:
`remotes::install_github("jinworks/CellChat")` (requires
**BiocNeighbors** from Bioconductor). Pre-rendered figures and the pair
summary table are in `inst/figures/` and `inst/data/`; regenerate with:

`Rscript scripts/render_spatial_colocalization_cellchat.R`

The code chunk below mirrors the render script; it is not executed
during knitting (use the pre-rendered outputs or run the script
locally).

``` r

if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("Install CellChat: remotes::install_github(\"jinworks/CellChat\")")
}
if (!exists("scoc_df") || !exists("mat_z") || !exists("mat_local")) {
  stop("Run the co-localization chunk first (or use pre-rendered figures)")
}

# Spatial CellChat on CytoSPACE prognostic tail cells (Visium coordinates + distance constraints)
meta_cc <- seurat@meta.data
meta_cc$percentile <- dplyr::percent_rank(meta_cc[[score_col]])
meta_cc$prognostic_group <- dplyr::case_when(
  meta_cc$percentile < 0.05 ~ "Favorable",
  meta_cc$percentile >= 0.95 ~ "Adverse",
  TRUE ~ "Other"
)
meta_cc$ct_pg <- paste0(meta_cc$CellType, "_", meta_cc$prognostic_group)
cells_cc <- rownames(meta_cc)[meta_cc$prognostic_group %in% c("Adverse", "Favorable")]
expr_counts <- Seurat::GetAssayData(seurat, assay = "Spatial", layer = "counts")[, cells_cc, drop = FALSE]
cs <- Matrix::colSums(expr_counts)
expr_norm <- log1p(Matrix::t(Matrix::t(expr_counts) / pmax(cs, 1) * 10000))
meta_cc <- meta_cc[cells_cc, , drop = FALSE]
meta_cc$labels <- meta_cc$ct_pg
meta_cc$samples <- factor(meta_cc$orig.ident)
keep_labs <- names(table(meta_cc$labels))[table(meta_cc$labels) >= 10L]
cells_cc <- rownames(meta_cc)[meta_cc$labels %in% keep_labs]

spatial.locs <- Seurat::GetTissueCoordinates(
  seurat, scale = NULL, cols = c("imagerow", "imagecol")
)[cells_cc, , drop = FALSE]
spot.size <- 65
spatial.factors <- data.frame(
  ratio = spot.size / seurat@images[[Seurat::Images(seurat)[1]]]@scale.factors$spot,
  tol = spot.size / 2
)
d_um <- CellChat::computeCellDistance(
  coordinates = spatial.locs,
  ratio = spatial.factors$ratio,
  tol = spatial.factors$tol
)
min_d <- min(d_um[d_um > 0], na.rm = TRUE)
scale.distance <- min(1, 1.1 / min_d)

cellchat <- CellChat::createCellChat(
  object = expr_norm[, cells_cc, drop = FALSE],
  meta = meta_cc[cells_cc, , drop = FALSE],
  group.by = "labels",
  datatype = "spatial",
  coordinates = spatial.locs,
  spatial.factors = spatial.factors
)
cellchat@DB <- CellChat::CellChatDB.human
cellchat <- CellChat::subsetData(cellchat)
cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat::computeCommunProb(
  cellchat,
  type = "truncatedMean",
  trim = 0.1,
  distance.use = TRUE,
  interaction.range = 250,
  scale.distance = scale.distance,
  contact.dependent = TRUE,
  contact.range = 100
)
cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
cellchat <- CellChat::computeCommunProbPathway(cellchat)
cellchat <- CellChat::aggregateNet(cellchat)

cc_df <- CellChat::subsetCommunication(cellchat)
cc_agg <- cc_df %>%
  dplyr::group_by(source, target) %>%
  dplyr::summarise(
    cellchat_prob = sum(.data$prob, na.rm = TRUE),
    cellchat_n_lr = dplyr::n(),
    .groups = "drop"
  )

# Merge spatial + CellChat pair metrics (full pipeline in render script)
pairs_path <- .prerender_data_path(.spatial_cellchat_pairs_rds, must_exist = FALSE)
if (!is.na(pairs_path)) {
  pair_prognostic <- readRDS(pairs_path)
}
```

The bubble plot below compares **neighborhood enrichment** (x) with
**spatial CellChat** aggregated L-R probability (y), using CytoSPACE
cell locations and Visium distance constraints. Point size reflects
**local co-localization**; red points are pairs with both spatial
co-localization and predicted signaling.

![](../reference/figures/spatial_coloc_cellchat_bubble.png)

The integrated heatmap clusters prognostic sender → receiver pairs by a
combined score (enrichment z-score, local co-localization, and spatial
CellChat probability).

![](../reference/figures/spatial_coloc_cellchat_dual_heatmap.png)

![](../reference/figures/spatial_coloc_cellchat_chord.png)

Top pathways for the strongest spatial + signaling pairs
(e.g. **Fibroblast_Adverse → Ductal_Adverse**):

![](../reference/figures/spatial_coloc_cellchat_lr_bubble.png)

#### Physical vs secreted CellChat interactions

CellChat annotates each ligand–receptor pair as **Cell-Cell Contact**,
**ECM-Receptor** (grouped here as **physical**), or **Secreted
Signaling** (**secreted**). We ask whether physical interactions track
**neighborhood enrichment** and **local co-localization** more strongly
than secreted interactions, and whether **physical interaction
fraction** correlates with PhenoMapR score associations.

![](../reference/figures/spatial_coloc_cellchat_type_spatial.png)

![](../reference/figures/spatial_coloc_cellchat_type_pheno.png)

![](../reference/figures/spatial_coloc_cellchat_type_assoc.png)

#### Integrated spatial + PhenoMapR + CellChat heatmaps

The four-panel heatmap below uses the same row/column order as the
PhenoMapR neighbor Spearman matrix: **(1)** neighbor PhenoMapR ρ with
dual-evidence outlines, **(2)** log1p spatial CellChat probability,
**(3)** integrated score (z-scaled spatial, PhenoMapR ρ, and CellChat),
and **(4)** dual-positive flag.

![](../reference/figures/spatial_coloc_integrated_four_panel.png)

![](../reference/figures/spatial_coloc_cellchat_prob_heatmap_clustered.png)

![](../reference/figures/spatial_coloc_integration_evidence.png)

Spearman heatmaps with **black** outlines for **\|ρ\| \> 0.5** and
**purple** outlines for pairs with both spatial co-localization and
spatial CellChat L-R support:

![](../reference/figures/spatial_colocalization_pheno_vs_cooccur_spearman_clustered_outlined_integrated.png)

![](../reference/figures/spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered_outlined_integrated.png)

#### Top co-localized expression axes

The plots below summarize the strongest **spatially co-localized L-R
expression axes** among dual-positive cross-type pairs
(CellChat-inspired topic composition and spatial axis maps).

![](../reference/figures/spatial_coloc_lr_topic_composition.png)

![](../reference/figures/spatial_coloc_pathway_topic_composition.png)

![](../reference/figures/spatial_coloc_expression_axis_heatmap.png)

![](../reference/figures/spatial_coloc_top_lr_spatial_axes.png)

#### Spatial pair evidence maps (CytoSPACE locations)

For the top **dual spatial + CellChat** cross-type pairs, seven-panel
maps show prognostic tail cells, sender/receiver roles, PhenoMapR score,
local co-localization, CellChat sender L-R potential, integrated
hotspots, and target-neighbor fraction — all on **CytoSPACE
`row`/`col`** coordinates. Each pair is saved in **scaled** (point size
∝ 1/cells per spot) and **uniform** (fixed point size) versions under
`inst/figures/spatial_pair_maps/`; regenerate with
`Rscript scripts/render_spatial_pair_spatial_maps.R` (after the CellChat
render script).

**Interpretation.** Many high-scoring pairs are **within-group**
(e.g. adverse ductal cells enriched near other adverse ductal cells),
consistent with spatial clustering of extreme phenotypes. Cross-type
pairs such as **Fibroblast_Adverse → Ductal_Adverse** combine moderate
neighborhood enrichment, local mixing, and spatial CellChat L-R
probability—candidate **stromal–epithelial** interactions among adverse
prognostic cells that are both co-localized and within signaling range
on tissue.

**Physical vs secreted.** The association table above summarizes how
**physical** (contact/ECM) vs **secreted** CellChat probabilities relate
to neighborhood enrichment, local co-localization, and PhenoMapR spatial
correlations among prognostic pairs with signaling support. In general,
**physical** interaction probability tracks immediate co-localization
more strongly than **secreted** signaling, which can remain nonzero over
broader neighborhoods where neighbor PhenoMapR scores covary.

### Prognostic markers

[`find_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/find_phenotype_markers.md)
supports the same **marker scopes** on spatial data (cells mapped to
spots) as in the single-cell vignette:

- **Cell type agnostic** — `marker_scope = "phenotype_groups"`: **Most
  Adverse** vs all other labeled cells and **Most Favorable** vs all
  other labeled cells, **ignoring** cell type.
- **Cell type × phenotype vs all opposite-tail cells** —
  `marker_scope = "cell_type_specific"`,
  `celltype_contrast = "vs_opposite_tail"`: for each cell type, that
  type in one tail vs **every cell in the opposite tail** (works when a
  type exists in only one tail).
- **Cell type specific (within cell type)** —
  `marker_scope = "cell_type_specific"`, default
  `celltype_contrast = "within_cell_type"`: that type in one tail vs
  **other phenotype bins within the same type**.

Below we run all three contrasts and draw
**[`plot_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/plot_phenotype_markers.md)**
heatmaps with **five genes per block** and **five `anno_mark` labels per
block** (`top_n_markers = 5`, `n_mark_labels = 5`). Gene labels are
colored by cell type (`color_mark_labels_by_celltype = TRUE`) on the
cell-type-specific heatmaps. On CI/pkgdown builds, only the **global**
contrast is rendered to stay within memory limits; run locally for all
three heatmaps.

#### Step 1: Phenotype groups and marker discovery

``` r

markers_global <- NULL
markers_ct_within <- NULL
markers_ct_opp <- NULL

cells <- colnames(seurat)
scores_df_markers <- seurat@meta.data[cells, score_col, drop = FALSE]
rownames(scores_df_markers) <- cells
groups_markers <- PhenoMapR::define_phenotype_groups(
  scores_df_markers,
  percentile = 0.05,
  score_columns = score_col
)
group_col <- grep("phenotype_group", names(groups_markers), value = TRUE)[1]
seurat@meta.data[[group_col]] <-
  groups_markers[rownames(seurat@meta.data), group_col, drop = TRUE]

ct_col_markers <- NULL
if (exists("spatial_df_celltype_col") && !is.null(spatial_df_celltype_col) &&
    spatial_df_celltype_col %in% names(seurat@meta.data)) {
  ct_col_markers <- spatial_df_celltype_col
} else if (exists("celltype_col") && !is.null(celltype_col) &&
           celltype_col %in% names(seurat@meta.data)) {
  ct_col_markers <- celltype_col
}

group_vec <- seurat@meta.data[cells, group_col]
group_df <- data.frame(
  cell_id = cells,
  phenotype_group = as.character(group_vec),
  stringsAsFactors = FALSE
)
if (!is.null(ct_col_markers)) {
  group_df$cell_type <- as.character(seurat@meta.data[cells, ct_col_markers])
}

max_cells_markers <- if (.is_vignette_ci()) 2500L else 5000L

markers_global <- PhenoMapR::find_phenotype_markers(
  seurat,
  group_labels = group_df,
  group_column = "phenotype_group",
  cell_id_column = "cell_id",
  marker_scope = "phenotype_groups",
  assay = "Spatial",
  layer = "counts",
  max_cells_per_ident = max_cells_markers,
  verbose = FALSE
)

if (.spatial_run_full_markers && !is.null(ct_col_markers)) {
  markers_ct_opp <- PhenoMapR::find_phenotype_markers(
    seurat,
    group_labels = group_df,
    group_column = "phenotype_group",
    cell_id_column = "cell_id",
    cell_type_column = "cell_type",
    marker_scope = "cell_type_specific",
    celltype_contrast = "vs_opposite_tail",
    assay = "Spatial",
    layer = "counts",
    max_cells_per_ident = max_cells_markers,
    verbose = FALSE
  )
  markers_ct_within <- PhenoMapR::find_phenotype_markers(
    seurat,
    group_labels = group_df,
    group_column = "phenotype_group",
    cell_id_column = "cell_id",
    cell_type_column = "cell_type",
    marker_scope = "cell_type_specific",
    celltype_contrast = "within_cell_type",
    assay = "Spatial",
    layer = "counts",
    max_cells_per_ident = max_cells_markers,
    verbose = FALSE
  )
}
```

#### Step 2: Heatmap prep (expression matrix + metadata)

``` r

expr_pm <- NULL
meta_pm <- NULL
group_col_hm <- NULL
ct_col_hm <- NULL
pal_ct <- NULL

.get_seurat_layer <- function(obj, assay, layer) {
  mat <- tryCatch(
    Seurat::GetAssayData(obj, layer = layer, assay = assay),
    error = function(e) tryCatch(Seurat::GetAssayData(obj, slot = layer, assay = assay), error = function(e2) NULL)
  )
  if (is.null(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
    mat <- tryCatch(
      SeuratObject::LayerData(obj, layer = layer, assay = assay),
      error = function(e) NULL
    )
  }
  if (is.null(mat) || nrow(mat) == 0L || ncol(mat) == 0L) NULL else mat
}

if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)) {
  genes_hm <- character(0)
  for (m in list(markers_global, markers_ct_opp, markers_ct_within)) {
    if (is.null(m)) next
    for (df in list(m$adverse_markers, m$favorable_markers)) {
      if (!is.null(df) && nrow(df) > 0L && "gene" %in% names(df)) {
        genes_hm <- union(genes_hm, df$gene)
      }
    }
  }

  expr <- NULL
  layer_used <- "data"
  assay_use_hm <- if (exists("assay_use") && !is.null(assay_use)) as.character(assay_use)[1] else "Spatial"
  assay_order <- unique(c(assay_use_hm, "Spatial", "RNA", "SCT"))
  for (a in assay_order) {
    if (!a %in% names(seurat@assays)) next
    expr <- .get_seurat_layer(seurat, a, "data")
    layer_used <- "data"
    if (is.null(expr)) {
      expr <- .get_seurat_layer(seurat, a, "counts")
      layer_used <- "counts"
    }
    if (!is.null(expr)) break
  }

  if (!is.null(expr) && ncol(expr) > 0) {
    cells_expr <- colnames(expr)
    if (is.null(cells_expr)) cells_expr <- character(0)
    cells_obj <- colnames(seurat)
    cells_meta <- rownames(seurat@meta.data)
    cells_use <- intersect(cells_expr, cells_obj)
    if (length(cells_use) == 0 && length(cells_meta) > 0) {
      cells_use <- intersect(cells_expr, cells_meta)
    }
    if (length(cells_use) == 0 && length(cells_expr) > 0 && any(grepl("_", cells_expr, fixed = TRUE))) {
      stripped <- sub("^[^_]+_", "", cells_expr)
      cells_use <- intersect(stripped, cells_obj)
      if (length(cells_use) == 0) cells_use <- intersect(stripped, cells_meta)
      if (length(cells_use) > 0) {
        expr_cols <- vapply(cells_use, function(c) cells_expr[which(stripped == c)[1]], character(1))
        expr <- expr[, expr_cols, drop = FALSE]
        colnames(expr) <- cells_use
      }
    }
    if (length(cells_use) == 0) {
      message(
        "No overlapping cells between expression and metadata; skipping marker heatmaps. ",
        "expr ncol=", length(cells_expr), ", obj ncol=", length(cells_obj)
      )
    } else {
      max_cells_hm <- if (.is_vignette_ci()) 2000L else 4000L
      if (length(cells_use) > max_cells_hm) {
        set.seed(5L)
        cells_use <- sample(cells_use, max_cells_hm)
      }
      expr_sub <- expr[, cells_use, drop = FALSE]
      if (length(genes_hm) > 0) {
        genes_use <- intersect(genes_hm, rownames(expr_sub))
        if (length(genes_use) > 0) {
          expr_sub <- expr_sub[genes_use, , drop = FALSE]
        }
      }
      if (identical(layer_used, "counts")) {
        if (requireNamespace("Matrix", quietly = TRUE) &&
            (inherits(expr_sub, "Matrix") || inherits(expr_sub, "sparseMatrix"))) {
          cs <- Matrix::colSums(expr_sub)
          expr_sub <- Matrix::t(Matrix::t(expr_sub) / pmax(cs, 1) * 10000)
          expr_sub <- log1p(expr_sub)
        } else {
          cs <- colSums(expr_sub)
          expr_sub <- t(t(expr_sub) / pmax(cs, 1) * 10000)
          expr_sub <- log1p(expr_sub)
        }
      }
      expr_pm <- as.matrix(expr_sub)
      if (exists("expr", inherits = FALSE)) rm(expr)
      gc(verbose = FALSE)
      meta_pm <- seurat@meta.data[cells_use, , drop = FALSE]
      meta_pm$cell_id_plot <- cells_use
      group_col_hm <- paste0("phenotype_group_", score_col)
      if (!group_col_hm %in% names(meta_pm)) group_col_hm <- group_col
      ct_col_hm <- ct_col_markers
      if (is.null(ct_col_hm) || !ct_col_hm %in% names(meta_pm)) {
        meta_pm$celltype_plot <- factor(rep("Cell", nrow(meta_pm)))
        ct_col_hm <- "celltype_plot"
      } else {
        meta_pm[[ct_col_hm]] <- factor(as.character(meta_pm[[ct_col_hm]]))
      }
      pal_ct <- PhenoMapR::get_celltype_palette(levels(meta_pm[[ct_col_hm]]))
    }
  } else {
    message("Could not extract a non-empty expression matrix for marker heatmaps.")
  }
}
```

#### Step 3: Cell-type-agnostic marker heatmap

``` r

if (!is.null(expr_pm) && !is.null(markers_global)) {
  adv_ok <- !is.null(markers_global$adverse_markers) && nrow(markers_global$adverse_markers) > 0L
  fav_ok <- !is.null(markers_global$favorable_markers) && nrow(markers_global$favorable_markers) > 0L
  if (adv_ok && fav_ok) {
    PhenoMapR::plot_phenotype_markers(
      markers = markers_global,
      expr_mat = expr_pm,
      meta = meta_pm,
      cell_id_col = "cell_id_plot",
      group_col = group_col_hm,
      score_col = score_col,
      celltype_col = ct_col_hm,
      celltype_palette = pal_ct,
      heatmap_type = "global",
      top_n_markers = 5L,
      n_mark_labels = 5L,
      p_adj_threshold = 0.05,
      column_title = "Global phenotype marker genes (CytoSPACE cells)"
    )
  }
}
```

#### Step 4: Cell-type × phenotype vs all opposite-tail cells

``` r

if (.spatial_run_full_markers && !is.null(expr_pm) && !is.null(markers_ct_opp)) {
  PhenoMapR::plot_phenotype_markers(
    markers = markers_ct_opp,
    expr_mat = expr_pm,
    meta = meta_pm,
    cell_id_col = "cell_id_plot",
    group_col = group_col_hm,
    score_col = score_col,
    celltype_col = ct_col_hm,
    celltype_palette = pal_ct,
    heatmap_type = "cell_type_specific",
    top_n_markers = 5L,
    n_mark_labels = 5L,
    color_mark_labels_by_celltype = TRUE,
    outline_marker_blocks = TRUE,
    block_outline_color = "black",
    p_adj_threshold = 0.05,
    column_title = "Cell-type markers vs all opposite-tail cells (CytoSPACE cells)"
  )
}
```

#### Step 5: Cell-type-specific markers (within cell type)

``` r

if (.spatial_run_full_markers && !is.null(expr_pm) && !is.null(markers_ct_within)) {
  PhenoMapR::plot_phenotype_markers(
    markers = markers_ct_within,
    expr_mat = expr_pm,
    meta = meta_pm,
    cell_id_col = "cell_id_plot",
    group_col = group_col_hm,
    score_col = score_col,
    celltype_col = ct_col_hm,
    celltype_palette = pal_ct,
    heatmap_type = "cell_type_specific",
    top_n_markers = 5L,
    n_mark_labels = 5L,
    color_mark_labels_by_celltype = TRUE,
    outline_marker_blocks = TRUE,
    block_outline_color = "black",
    p_adj_threshold = 0.05,
    column_title = "Cell-type-specific phenotype marker genes (CytoSPACE cells)"
  )
}
```

#### Stacked barplot: adverse vs. favorable by cell type

Each column is a cell type; the height is filled by the number of
adverse (5th percentile) or favorable (5th percentile) cells.

``` r

if (!is.null(group_col)) {
  ct_col_bar <- NULL
  if (exists("spatial_df_celltype_col") && !is.null(spatial_df_celltype_col) && spatial_df_celltype_col %in% names(seurat@meta.data)) {
    ct_col_bar <- spatial_df_celltype_col
  } else if (exists("celltype_col") && !is.null(celltype_col) && celltype_col %in% names(seurat@meta.data)) {
    ct_col_bar <- celltype_col
  } else {
    for (c in c("CellType", "Celltype..major.lineage.", "cell_type", "celltype", "annotation", "seurat_clusters")) {
      if (c %in% names(seurat@meta.data) && length(unique(na.omit(seurat@meta.data[[c]]))) >= 2) {
        ct_col_bar <- c
        break
      }
    }
  }
  if (!is.null(ct_col_bar)) {
    meta <- seurat@meta.data
    meta$pg <- meta[[group_col]]
    meta$ct <- as.character(meta[[ct_col_bar]])
    meta$ct_ok <- !is.na(meta$ct) & nzchar(trimws(meta$ct))
    idx_extreme <- meta$pg %in% c("Most Adverse", "Most Favorable") & meta$ct_ok
    df_bar <- as.data.frame(table(
      CellType = meta$ct[idx_extreme],
      Prognostic_group = meta$pg[idx_extreme],
      useNA = "no"
    ))
    if (nrow(df_bar) > 0) {
      df_bar$Prognostic_group <- factor(df_bar$Prognostic_group, levels = c("Most Favorable", "Most Adverse"))
      pal_bar <- c(`Most Adverse` = "#B2182B", `Most Favorable` = "#2166AC")
      adverse <- df_bar[df_bar$Prognostic_group == "Most Adverse", c("CellType", "Freq")]
      fav <- df_bar[df_bar$Prognostic_group == "Most Favorable", c("CellType", "Freq")]
      names(adverse)[2] <- "n_adverse"
      names(fav)[2] <- "n_fav"
      df_labels <- merge(adverse, fav, by = "CellType", all = TRUE)
      df_labels$n_adverse[is.na(df_labels$n_adverse)] <- 0
      df_labels$n_fav[is.na(df_labels$n_fav)] <- 0
      df_labels$label <- paste0(df_labels$n_adverse, "/", df_labels$n_fav)
      df_labels$total <- df_labels$n_adverse + df_labels$n_fav
      ct_ord <- levels(reorder(df_bar$CellType, df_bar$Freq, function(x) -sum(x)))
      df_labels$CellType <- factor(df_labels$CellType, levels = ct_ord)
      p_bar <- ggplot(df_bar, aes(x = reorder(.data$CellType, .data$Freq, function(x) -sum(x)), y = .data$Freq, fill = .data$Prognostic_group)) +
        geom_col(position = "stack") +
        # ggchicklet2::geom_chicklet_bar(radius = grid::unit(1.5, "pt"), position = "stack") +
        geom_text(data = df_labels, aes(x = .data$CellType, y = .data$total, label = .data$n_adverse),
                  inherit.aes = FALSE, vjust = -0.3, size = 3.5, color = "#B2182B", hjust = 1.05) +
        geom_text(data = df_labels, aes(x = .data$CellType, y = .data$total, label = paste0("/", .data$n_fav)),
                  inherit.aes = FALSE, vjust = -0.3, size = 3.5, color = "#2166AC", hjust = -0.05) +
        scale_fill_manual(values = pal_bar) +
        labs(x = "Cell type", y = "Number of cells", fill = "Prognostic group",
             title = "Adverse vs. favorable cells (5th percentile) by cell type") +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), legend.position = "top")
      print(p_bar)
    }
  }
}
```

![](spatial-transcriptomics_files/figure-html/stacked-bar-adverse-favorable-by-celltype-1.png)
In line with the other Pancreatic cancer vignettes, PhenoMapR identifies
ductal cells (and some fibroblasts) as the most associated cell type
with an adverse PhenoMapR prognostic score. Alpha, plasma, and beta
cells are the most associated with the more favorable prognostic signal.
Interestingly, fibroblasts seems to comprise both adverse and favorable
phenotypes.

### Summary

- **Part 1 (spots)**: Spot-level object
  `HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test.hgnc.rds`; **H&E** via
  `SpatialPlot`; PhenoMapR scoring; **hex** summary of z-scaled scores
  plus **point** map of 5th-percentile groups (**same `coord_fixed` /
  limits as H&E**; no jitter).
- **Part 2 (CytoSPACE)**:
  `HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds`; same PRECOG
  **Pancreatic** reference at **cell** resolution; score distribution,
  score by cell type, prognostic groups, and spatial maps (jittered
  cells per spot).
- **Co-localization**: pre-rendered
  **`spatialCooccur::nhood_enrichment()`** and
  **`spatialCooccur::cooccur_local()`** heatmaps (ordered and
  clustered), plus **Spearman correlations** of PhenoMapR scores with
  local co-localization scores and with co-occurring neighbors’
  PhenoMapR scores; regenerate with
  `scripts/render_spatial_colocalization_heatmap.R`.
- **Spatial CellChat integration**: co-localization metrics merged with
  distance-constrained **[CellChat
  v2](https://github.com/jinworks/CellChat)** L-R inference on
  CytoSPACE-mapped cells (`datatype = "spatial"`); integrated four-panel
  heatmaps, dual-evidence Spearman outlines (black = \|ρ\| \> 0.5,
  purple = spatial + CellChat), and evidence-tier summary; regenerate
  with `scripts/render_spatial_colocalization_cellchat.R` (or
  `scripts/render_spatial_integrated_heatmaps.R` from cached
  `inst/data/spatial_coloc_cellchat_pairs.rds`).
- **Markers**: three
  **[`plot_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/plot_phenotype_markers.md)**
  heatmaps on CytoSPACE cells — (1) **cell-type agnostic**
  (`heatmap_type = "global"`), (2) **cell type × phenotype vs all
  opposite-tail cells** (`celltype_contrast = "vs_opposite_tail"`), (3)
  **within cell type** (`celltype_contrast = "within_cell_type"`); each
  uses **five genes and five `anno_mark` labels per block**.

### References

**\[1\]** Benard, B. A. et al. PRECOG update: an augmented resource of
clinical outcome associations with gene expression for adult, pediatric,
and immunotherapy cohorts. Nucleic Acids Res. 54, D1579–D1589 (2026).

### Session Info

``` r

sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggchicklet2_0.7.0  patchwork_1.3.2    dplyr_1.2.1        ggplot2_4.0.3     
    ## [5] Seurat_5.5.1       SeuratObject_5.4.0 sp_2.2-1           PhenoMapR_0.1.0   
    ## [9] testthat_3.3.2    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3     shape_1.4.6.1          jsonlite_2.0.0        
    ##   [4] magrittr_2.0.5         spatstat.utils_3.2-3   farver_2.1.2          
    ##   [7] rmarkdown_2.31         GlobalOptions_0.1.4    fs_2.1.0              
    ##  [10] ragg_1.5.2             vctrs_0.7.3            ROCR_1.0-12           
    ##  [13] spatstat.explore_3.8-1 htmltools_0.5.9        progress_1.2.3        
    ##  [16] curl_7.1.0             sass_0.4.10            sctransform_0.4.3     
    ##  [19] parallelly_1.48.0      KernSmooth_2.23-26     bslib_0.11.0          
    ##  [22] htmlwidgets_1.6.4      desc_1.4.3             ica_1.0-3             
    ##  [25] plyr_1.8.9             plotly_4.12.0          zoo_1.8-15            
    ##  [28] cachem_1.1.0           igraph_2.3.3           iterators_1.0.14      
    ##  [31] mime_0.13              lifecycle_1.0.5        pkgconfig_2.0.3       
    ##  [34] Matrix_1.7-5           R6_2.6.1               fastmap_1.2.0         
    ##  [37] clue_0.3-68            fitdistrplus_1.2-6     future_1.70.0         
    ##  [40] shiny_1.14.0           digest_0.6.39          colorspace_2.1-3      
    ##  [43] S4Vectors_0.50.1       rprojroot_2.1.1        tensor_1.5.1          
    ##  [46] RSpectra_0.16-2        irlba_2.3.7            pkgload_1.5.3         
    ##  [49] textshaping_1.0.5      labeling_0.4.3         progressr_1.0.0       
    ##  [52] spatstat.sparse_3.2-0  httr_1.4.8             polyclip_1.10-7       
    ##  [55] abind_1.4-8            compiler_4.6.1         gargle_1.6.1          
    ##  [58] doParallel_1.0.17      withr_3.0.3            S7_0.2.2              
    ##  [61] fastDummies_1.7.6      hexbin_1.28.5          pkgbuild_1.4.8        
    ##  [64] MASS_7.3-65            rjson_0.2.23           tools_4.6.1           
    ##  [67] lmtest_0.9-40          otel_0.2.0             googledrive_2.1.2     
    ##  [70] httpuv_1.6.17          future.apply_1.20.2    goftest_1.2-3         
    ##  [73] glue_1.8.1             nlme_3.1-169           promises_1.5.0        
    ##  [76] grid_4.6.1             Rtsne_0.17             cluster_2.1.8.2       
    ##  [79] reshape2_1.4.5         generics_0.1.4         gtable_0.3.6          
    ##  [82] spatstat.data_3.1-9    tidyr_1.3.2            data.table_1.18.4     
    ##  [85] hms_1.1.4              BiocGenerics_0.58.1    spatstat.geom_3.8-1   
    ##  [88] RcppAnnoy_0.0.23       foreach_1.5.2          ggrepel_0.9.8         
    ##  [91] RANN_2.6.2             pillar_1.11.1          stringr_1.6.0         
    ##  [94] spam_2.11-4            RcppHNSW_0.7.0         later_1.4.8           
    ##  [97] circlize_0.4.18        splines_4.6.1          lattice_0.22-9        
    ## [100] survival_3.8-6         deldir_2.0-4           tidyselect_1.2.1      
    ## [103] ComplexHeatmap_2.28.0  miniUI_0.1.2           pbapply_1.7-4         
    ## [106] knitr_1.51             gridExtra_2.3.1        IRanges_2.46.0        
    ## [109] scattermore_1.2        stats4_4.6.1           xfun_0.60             
    ## [112] brio_1.1.5             matrixStats_1.5.0      stringi_1.8.7         
    ## [115] lazyeval_0.2.3         yaml_2.3.12            evaluate_1.0.5        
    ## [118] codetools_0.2-20       tibble_3.3.1           cli_3.6.6             
    ## [121] uwot_0.2.4             xtable_1.8-8           reticulate_1.46.0     
    ## [124] systemfonts_1.3.2      jquerylib_0.1.4        Rcpp_1.1.2            
    ## [127] globals_0.19.1         spatstat.random_3.5-0  png_0.1-9             
    ## [130] spatstat.univar_3.2-0  parallel_4.6.1         pkgdown_2.2.1         
    ## [133] prettyunits_1.2.0      dotCall64_1.2          listenv_1.0.0         
    ## [136] viridisLite_0.4.3      scales_1.4.0           ggridges_0.5.7        
    ## [139] purrr_1.2.2            crayon_1.5.3           GetoptLong_1.1.1      
    ## [142] rlang_1.3.0            cowplot_1.2.0
