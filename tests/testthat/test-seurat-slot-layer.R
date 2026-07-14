# Seurat v4 slot / v5 layer compatibility for PhenoMapR public APIs

test_that(".resolve_seurat_layer_or_slot prefers layer and detects conflicts", {
  expect_equal(.resolve_seurat_layer_or_slot(layer = "data", slot = NULL), "data")
  expect_equal(.resolve_seurat_layer_or_slot(layer = "counts", slot = NULL), "counts")
  # Legacy slot= alone still works when layer is default
  expect_equal(.resolve_seurat_layer_or_slot(layer = "data", slot = "counts"), "counts")
  expect_equal(.resolve_seurat_layer_or_slot(layer = "counts", slot = "counts"), "counts")
  expect_error(
    .resolve_seurat_layer_or_slot(layer = "counts", slot = "data"),
    "Conflicting Seurat matrix selection"
  )
})

test_that(".get_assay_data_compat reads counts/data on installed Seurat", {
  skip_if_not_installed("Seurat")
  set.seed(1)
  counts <- matrix(
    rpois(40, 2),
    nrow = 10,
    ncol = 4,
    dimnames = list(paste0("g", 1:10), paste0("c", 1:4))
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  m_counts <- .get_assay_data_compat(obj, assay = "RNA", slot = "counts")
  expect_equal(dim(m_counts), dim(counts))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  m_data <- .get_assay_data_compat(obj, assay = "RNA", slot = "data")
  expect_equal(dim(m_data), dim(counts))
})

test_that(".get_assay_data_compat falls back from empty data to counts", {
  skip_if_not_installed("Seurat")
  set.seed(2)
  counts <- matrix(
    rpois(40, 2),
    nrow = 10,
    ncol = 4,
    dimnames = list(paste0("g", 1:10), paste0("c", 1:4))
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  # Fresh Assays often have empty 'data' until NormalizeData; compat should
  # still return a usable counts matrix when 'data' is requested.
  expect_warning(
    m <- .get_assay_data_compat(obj, assay = "RNA", slot = "data"),
    "using 'counts' instead"
  )
  expect_equal(dim(m), dim(counts))
  expect_true(!is.null(rownames(m)))
})

test_that("PhenoMap accepts layer= (preferred) and slot= alias on Seurat objects", {
  skip_if_not_installed("Seurat")
  data(precog, package = "PhenoMapR", envir = environment())
  genes <- rownames(precog)[seq_len(min(25, nrow(precog)))]
  n_cells <- 6
  counts <- matrix(
    pmax(0, rnorm(length(genes) * n_cells)),
    nrow = length(genes),
    ncol = n_cells,
    dimnames = list(genes, paste0("C", seq_len(n_cells)))
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)

  scores_layer <- PhenoMap(
    expression = obj,
    reference = "precog",
    cancer_type = "Breast",
    assay = "RNA",
    layer = "data",
    verbose = FALSE
  )
  scores_slot <- PhenoMap(
    expression = obj,
    reference = "precog",
    cancer_type = "Breast",
    assay = "RNA",
    slot = "data",
    verbose = FALSE
  )
  expect_equal(scores_layer, scores_slot)
})

test_that("PhenoMap spatial Seurat respects layer/slot (not forced to counts)", {
  skip_if_not_installed("Seurat")
  data(precog, package = "PhenoMapR", envir = environment())
  genes <- rownames(precog)[seq_len(min(25, nrow(precog)))]
  n_cells <- 5
  counts <- matrix(
    pmax(0, rnorm(length(genes) * n_cells)),
    nrow = length(genes),
    ncol = n_cells,
    dimnames = list(genes, paste0("C", seq_len(n_cells)))
  )
  # Mimic a spatial object by having an assay named Spatial
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "Spatial"))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  expect_identical(detect_input_type(obj), "seurat_spatial")

  scores <- PhenoMap(
    expression = obj,
    reference = "precog",
    cancer_type = "Breast",
    assay = "Spatial",
    layer = "data",
    verbose = FALSE
  )
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), n_cells)
})
