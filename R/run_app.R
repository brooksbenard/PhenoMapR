#' Launch the PhenoMapR Shiny app
#'
#' Starts a local Shiny session that walks users through the full PhenoMapR
#' workflow: upload data, choose a reference, score, tag phenotype groups,
#' and find marker genes. The app is a thin UI over the same exported
#' functions you would call from a script (\code{PhenoMap},
#' \code{define_phenotype_groups}, \code{find_phenotype_markers},
#' \code{derive_reference_from_bulk}, etc.).
#'
#' On a remote workstation or HPC node, run the app with
#' \code{run_app(host = "0.0.0.0", port = 3838, launch.browser = FALSE)} and
#' point a browser at \code{http://<server>:<port>}. No data leaves the
#' machine the app runs on.
#'
#' @param host Character. Network interface to bind to (default \code{"127.0.0.1"};
#'   use \code{"0.0.0.0"} to expose on all interfaces for remote use).
#' @param port Integer or \code{NULL} (random free port).
#' @param launch.browser Logical. Open the app in the default browser
#'   (\code{TRUE}); on headless servers set \code{FALSE}.
#' @param max_upload_size_mb Numeric. Per-file upload limit in megabytes.
#'   Defaults to \code{Inf} (no limit), so users on a workstation or HPC node
#'   with plenty of memory can upload very large \code{.rds} / \code{.h5ad}
#'   objects without hitting Shiny's default 5 MB cap. Set to a finite number
#'   only if you want to cap uploads (e.g. on a shared web server).
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @return Invisibly, the return value of \code{shiny::runApp()}.
#'
#' @examples
#' \dontrun{
#' # Local laptop use:
#' PhenoMapR::run_app()
#'
#' # Remote workstation, listen on all interfaces:
#' PhenoMapR::run_app(host = "0.0.0.0", port = 3838, launch.browser = FALSE)
#' }
#'
#' @export
run_app <- function(host = "127.0.0.1",
                    port = NULL,
                    launch.browser = interactive(),
                    max_upload_size_mb = Inf,
                    ...) {

  required <- c("shiny", "bslib", "DT", "ggplot2", "dplyr")
  missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "PhenoMapR::run_app() needs these packages installed: ",
      paste(missing, collapse = ", "),
      ". Install with: install.packages(c(\"",
      paste(missing, collapse = "\", \""), "\"))."
    )
  }

  app_dir <- system.file("shiny", package = "PhenoMapR")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    # devtools::load_all() / source clone fallback: look next to this file.
    here <- tryCatch({
      this_ofile <- sys.frame(1)$ofile
      if (is.null(this_ofile)) "" else dirname(this_ofile)
    }, error = function(e) "")
    candidates <- c(
      file.path(here, "..", "inst", "shiny"),
      file.path(getwd(), "inst", "shiny")
    )
    for (cand in candidates) {
      cand <- normalizePath(cand, mustWork = FALSE)
      if (dir.exists(cand) && file.exists(file.path(cand, "app.R"))) {
        app_dir <- cand
        break
      }
    }
  }
  if (!nzchar(app_dir) || !dir.exists(app_dir) ||
      !file.exists(file.path(app_dir, "app.R"))) {
    stop(
      "Shiny app directory not found. If you are running from a clone, ",
      "set the working directory to the repo root and call ",
      "shiny::runApp(\"inst/shiny\"). If you installed PhenoMapR from ",
      "GitHub, make sure the install includes inst/shiny (re-install with ",
      "remotes::install_github(\"brooksbenard/PhenoMapR\"))."
    )
  }

  # Make sure the app's local files (helpers.R, www/) can be resolved when
  # Shiny's working directory differs from the package install location.
  old_env <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = NA_character_)
  Sys.setenv(PHENOMAPR_SHINY_DIR = app_dir)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("PHENOMAPR_SHINY_DIR")
    else Sys.setenv(PHENOMAPR_SHINY_DIR = old_env)
  }, add = TRUE)

  # Default behaviour is "no limit": Shiny treats the option as bytes and
  # compares actual upload size against it with `<`, so Inf effectively
  # disables the cap. Pass a finite max_upload_size_mb to cap explicitly.
  max_size_bytes <- if (is.finite(max_upload_size_mb)) {
    max_upload_size_mb * 1024 * 1024
  } else {
    Inf
  }
  old_opt <- getOption("shiny.maxRequestSize")
  options(shiny.maxRequestSize = max_size_bytes)
  on.exit(options(shiny.maxRequestSize = old_opt), add = TRUE)

  shiny::runApp(
    appDir = app_dir,
    host = host,
    port = port,
    launch.browser = launch.browser,
    ...
  )
}
