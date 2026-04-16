#' Replace Visium slide images using 10X `spatial` outputs
#'
#' Rebuilds each image slot with [Seurat::Read10X_Image()], which reads
#' `tissue_positions_list.csv` or `tissue_positions.csv`, `scalefactors_json.json`,
#' and the tissue PNG. This matches **spot centroids to the H&E** as Space Ranger
#' defines them—useful when a saved `.rds` has drifted or was built with
#' inconsistent coordinate logic.
#'
#' Set `spatial_dir` to the **`spatial` folder** inside the Visium run (the directory
#' that contains the CSV, JSON, and `tissue_lowres_image.png`).
#'
#' **Barcodes:** only spots/cells whose names appear in both `object` and the 10X
#' positions file are kept. If your object uses renamed barcodes (e.g. some
#' deconvolution outputs), intersection may be empty and the image slot is left
#' unchanged (a warning is issued).
#'
#' @param object A [Seurat::Seurat-class] object with at least one Visium image.
#' @param spatial_dir Path to the 10X **`spatial`** directory.
#' @param filter.matrix Passed to [Seurat::Read10X_Image()]. Use `FALSE` (default)
#'   to load all tissue positions, then keep only [colnames()] of `object`.
#' @param image.name PNG file name in `spatial_dir` (default lowres, as in Seurat).
#'
#' @return `object` with image slot(s) replaced for cells that overlap the new image.
#' @export
#'
#' @examples
#' \dontrun{
#' ## After Load10X_Spatial / download, if spots look misaligned:
#' obj <- sync_visium_spatial_image_from_10x_dir(obj, spatial_dir = "/path/to/spatial")
#' }
sync_visium_spatial_image_from_10x_dir <- function(
    object,
    spatial_dir,
    filter.matrix = FALSE,
    image.name = "tissue_lowres_image.png") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required.", call. = FALSE)
  }
  if (length(object@images) == 0L) {
    return(object)
  }
  if (!dir.exists(spatial_dir)) {
    stop("spatial_dir does not exist: ", spatial_dir, call. = FALSE)
  }
  spatial_dir <- normalizePath(spatial_dir, winslash = "/", mustWork = TRUE)
  pos_glob <- Sys.glob(file.path(spatial_dir, "*tissue_positions*"))
  if (length(pos_glob) == 0L) {
    warning(
      "No tissue_positions* file in ",
      spatial_dir,
      "; object unchanged.",
      call. = FALSE
    )
    return(object)
  }
  assay <- Seurat::DefaultAssay(object)
  inames <- names(object@images)
  for (iname in inames) {
    old_img <- object@images[[iname]]
    image_type <- if (inherits(old_img, "VisiumV1")) {
      "VisiumV1"
    } else {
      "VisiumV2"
    }
    img_new <- tryCatch(
      Seurat::Read10X_Image(
        image.dir = spatial_dir,
        image.name = image.name,
        assay = assay,
        slice = iname,
        filter.matrix = filter.matrix,
        image.type = image_type
      ),
      error = function(e) {
        warning(
          "Read10X_Image failed (",
          iname,
          "): ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
    if (is.null(img_new)) {
      next
    }
    keep <- intersect(SeuratObject::Cells(img_new), colnames(object))
    if (length(keep) == 0L) {
      warning(
        "No overlapping barcodes between object and Read10X_Image for ",
        iname,
        call. = FALSE
      )
      next
    }
    object[[iname]] <- img_new[keep]
  }
  object
}
