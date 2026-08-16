#!/usr/bin/env Rscript
# Compare DepMap 2D PAAD GSEA NES when phenotype tails are 5% vs 10%.
# Malignant-ductal (Ductal type 2) adverse/favorable markers; same FC/p filters.
#
# Usage (repo root):
#   Rscript manuscript/benchmarks/methodology/compare_depmap_tail_5_vs_10.R

suppressPackageStartupMessages({
  args_cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_cmd, value = TRUE)
  script_path <- if (length(file_arg)) {
    normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
  } else {
    normalizePath(
      file.path(
        "manuscript", "benchmarks", "methodology",
        "compare_depmap_tail_5_vs_10.R"
      ),
      winslash = "/",
      mustWork = FALSE
    )
  }
  meth_dir <- dirname(script_path)
  int_dir <- normalizePath(file.path(meth_dir, "..", "integration"), winslash = "/")
  source(file.path(int_dir, "config.R"), local = TRUE)
  source(file.path(meth_dir, "helpers_methodology.R"), local = TRUE)
  source(file.path(meth_dir, "depmap_helpers.R"), local = TRUE)
  root <- integration_repo_root()
  methodology_init(options = list(tier = "full", download = FALSE), root = root)
  pkgload::load_all(root, quiet = TRUE, helpers = FALSE, export_all = FALSE)
})

if (!requireNamespace("fgsea", quietly = TRUE)) stop("fgsea required")

rep_dir <- file.path(meth_dir, "reports")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

# Match manuscript DepMap GSEA panel (FC sweep: fav n.s. at log2FC >= 0.75)
marker_p_adj <- 0.05
marker_min_log2fc <- 0.75
nperm <- 10000L

message("Loading CRA001160 (full cohort, including non-tumor — manuscript DepMap setting)...")
dat <- methodology_load_cra001160(tier = "full", tumor_only = FALSE)
malignant_ductal <- methodology_malignant_ductal_celltypes()

message("PhenoMap scoring...")
scores <- PhenoMap(
  expression = dat$expression,
  reference = methodology_get_option("reference"),
  cancer_type = methodology_get_option("cancer_type"),
  verbose = FALSE
)

message("Loading DepMap PAAD gene scores...")
gene_scores <- methodology_depmap_paad_gene_scores()
universe <- rownames(dat$expression)

.gsea_row <- function(gsea_res, pathway) {
  tab <- gsea_res$table
  hit <- tab[as.character(tab$pathway) == pathway, , drop = FALSE]
  if (!nrow(hit)) {
    return(data.frame(
      NES = NA_real_, pval = NA_real_, padj = NA_real_,
      ES = NA_real_, size = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    NES = as.numeric(hit$NES[1]),
    pval = as.numeric(hit$pval[1]),
    padj = as.numeric(hit$padj[1]),
    ES = as.numeric(hit$ES[1]),
    size = as.integer(hit$size[1]),
    stringsAsFactors = FALSE
  )
}

.run_one_percentile <- function(percentile) {
  message(sprintf("\n=== Phenotype tails: %.0f%% ===", 100 * percentile))
  groups <- define_phenotype_groups(scores, percentile = percentile)
  gcol <- grep("^phenotype_group_", colnames(groups), value = TRUE)[1]
  id_col <- if (".cell_id" %in% names(dat$metadata)) ".cell_id" else "cell_id"
  groups$cell_type <- dat$metadata$celltype_original[
    match(groups$cell_id, dat$metadata[[id_col]])
  ]

  message("find_phenotype_markers (within cell type)...")
  mk <- find_phenotype_markers(
    expression = dat$expression,
    group_labels = groups,
    group_column = gcol,
    cell_id_column = "cell_id",
    cell_type_column = "cell_type",
    marker_scope = "cell_type_specific",
    celltype_contrast = "within_cell_type",
    celltype_unique_genes = FALSE,
    verbose = FALSE,
    max_cells_per_ident = 5000L
  )

  adv_raw <- methodology_filter_malignant_ductal_markers(mk$adverse_markers)
  fav_raw <- methodology_filter_malignant_ductal_markers(mk$favorable_markers)
  adv <- methodology_significant_positive_markers(
    adv_raw, p_adj_threshold = marker_p_adj, min_log2fc = marker_min_log2fc
  )
  fav <- methodology_significant_positive_markers(
    fav_raw, p_adj_threshold = marker_p_adj, min_log2fc = marker_min_log2fc
  )
  adv_genes <- unique(as.character(adv$gene))
  fav_genes <- unique(as.character(fav$gene))
  message(sprintf(
    "Malignant-ductal markers (p_adj<=%.2f, log2FC>=%.2f): adverse=%d genes, favorable=%d genes",
    marker_p_adj, marker_min_log2fc, length(adv_genes), length(fav_genes)
  ))
  message(sprintf(
    "  pre-filter raw: adverse=%d, favorable=%d",
    nrow(adv_raw), nrow(fav_raw)
  ))

  empty_gsea <- function() {
    list(
      table = data.frame(
        pathway = character(), NES = numeric(), pval = numeric(),
        padj = numeric(), ES = numeric(), size = integer(),
        stringsAsFactors = FALSE
      )
    )
  }

  if (length(adv_genes) < 10L && length(fav_genes) < 10L) {
    message("Skipping GSEA: both gene sets below minSize=10 after filters")
    gsea_eff <- empty_gsea()
    gsea_dep <- empty_gsea()
    dep_val <- list(rank_enrichment = NULL)
  } else {
    # Gene-effect essentiality GSEA (2D Chronos)
    gsea_eff <- tryCatch(
      methodology_depmap_gsea(
        adverse_genes = adv_genes,
        favorable_genes = fav_genes,
        gene_scores = gene_scores,
        universe = universe,
        nperm = nperm
      ),
      error = function(e) {
        message("gene_effect GSEA failed: ", conditionMessage(e))
        empty_gsea()
      }
    )

    # Dependency-probability GSEA via metric remapping
    dep_val <- tryCatch(
      methodology_depmap_metric_validation(
        adverse_genes = adv_genes,
        favorable_genes = fav_genes,
        gene_scores = gene_scores,
        universe = universe,
        metric = "dependency",
        dep_cutoff = 0.5
      ),
      error = function(e) {
        message("dependency validation failed: ", conditionMessage(e))
        list(gsea = empty_gsea(), rank_enrichment = NULL)
      }
    )
    gsea_dep <- if (!is.null(dep_val$gsea)) dep_val$gsea else empty_gsea()
  }

  overlap_adv <- length(intersect(adv_genes, gene_scores$gene))
  overlap_fav <- length(intersect(fav_genes, gene_scores$gene))

  ae <- .gsea_row(gsea_eff, "Adverse")
  fe <- .gsea_row(gsea_eff, "Favorable")
  ad <- .gsea_row(gsea_dep, "Adverse")
  fd <- .gsea_row(gsea_dep, "Favorable")

  rbind(
    data.frame(
      percentile = percentile,
      metric = "gene_effect",
      min_log2fc = marker_min_log2fc,
      n_adverse_genes = length(adv_genes),
      n_favorable_genes = length(fav_genes),
      n_adverse_raw = nrow(adv_raw),
      n_favorable_raw = nrow(fav_raw),
      n_adverse_in_depmap = overlap_adv,
      n_favorable_in_depmap = overlap_fav,
      adverse_NES = ae$NES,
      adverse_pval = ae$pval,
      adverse_padj = ae$padj,
      adverse_size = ae$size,
      favorable_NES = fe$NES,
      favorable_pval = fe$pval,
      favorable_padj = fe$padj,
      favorable_size = fe$size,
      stringsAsFactors = FALSE
    ),
    data.frame(
      percentile = percentile,
      metric = "dependency",
      min_log2fc = marker_min_log2fc,
      n_adverse_genes = length(adv_genes),
      n_favorable_genes = length(fav_genes),
      n_adverse_raw = nrow(adv_raw),
      n_favorable_raw = nrow(fav_raw),
      n_adverse_in_depmap = if (!is.null(dep_val$rank_enrichment)) {
        dep_val$rank_enrichment$n_adverse_in_depmap
      } else {
        overlap_adv
      },
      n_favorable_in_depmap = if (!is.null(dep_val$rank_enrichment)) {
        dep_val$rank_enrichment$n_favorable_in_depmap
      } else {
        overlap_fav
      },
      adverse_NES = ad$NES,
      adverse_pval = ad$pval,
      adverse_padj = ad$padj,
      adverse_size = ad$size,
      favorable_NES = fd$NES,
      favorable_pval = fd$pval,
      favorable_padj = fd$padj,
      favorable_size = fd$size,
      stringsAsFactors = FALSE
    )
  )
}

res <- rbind(.run_one_percentile(0.10), .run_one_percentile(0.05))
out <- file.path(rep_dir, "depmap_2d_gsea_nes_tail_5_vs_10_fc075.csv")
utils::write.csv(res, out, row.names = FALSE)

# Wide delta table
wide <- do.call(rbind, lapply(unique(res$metric), function(m) {
  d10 <- res[res$metric == m & res$percentile == 0.10, , drop = FALSE]
  d05 <- res[res$metric == m & res$percentile == 0.05, , drop = FALSE]
  data.frame(
    metric = m,
    min_log2fc = marker_min_log2fc,
    n_adv_10 = d10$n_adverse_genes,
    n_adv_5 = d05$n_adverse_genes,
    n_fav_10 = d10$n_favorable_genes,
    n_fav_5 = d05$n_favorable_genes,
    NES_adv_10 = d10$adverse_NES,
    NES_adv_5 = d05$adverse_NES,
    delta_NES_adv = d05$adverse_NES - d10$adverse_NES,
    pval_adv_10 = d10$adverse_pval,
    pval_adv_5 = d05$adverse_pval,
    padj_adv_10 = d10$adverse_padj,
    padj_adv_5 = d05$adverse_padj,
    NES_fav_10 = d10$favorable_NES,
    NES_fav_5 = d05$favorable_NES,
    delta_NES_fav = d05$favorable_NES - d10$favorable_NES,
    pval_fav_10 = d10$favorable_pval,
    pval_fav_5 = d05$favorable_pval,
    padj_fav_10 = d10$favorable_padj,
    padj_fav_5 = d05$favorable_padj,
    stringsAsFactors = FALSE
  )
}))
out_w <- file.path(rep_dir, "depmap_2d_gsea_nes_tail_5_vs_10_fc075_delta.csv")
utils::write.csv(wide, out_w, row.names = FALSE)

message("\n=== 5% vs 10% DepMap 2D GSEA NES (min_log2fc=", marker_min_log2fc, ") ===")
print(wide, row.names = FALSE)
message("Wrote ", out)
message("Wrote ", out_w)
invisible(res)
