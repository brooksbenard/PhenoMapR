# =============================================================================
# PhenoMapR — Shiny app
#
# Launch from the package: PhenoMapR::run_app()
# Launch from a clone:      shiny::runApp("inst/shiny")
# Launch on a server:       PhenoMapR::run_app(host = "0.0.0.0", port = 3838)
#
# This single-file app guides users through the full PhenoMapR workflow:
#   1. Upload expression data (matrix / Seurat / SCE / AnnData)
#   2. Pick a phenotype signature (built-in PRECOG/TCGA/Pediatric/ICI or custom)
#      — also includes reference diagnostics: gene coverage, top prognostic
#        genes, reference signature heatmap.
#   3. Compute weighted-sum scores with PhenoMap()
#   4. Tag adverse / favorable phenotype groups (merged into the Score tab)
#   5. Find phenotype marker genes (global or cell-type-specific)
# =============================================================================

suppressPackageStartupMessages({
  required <- c("shiny", "bslib", "DT", "ggplot2", "dplyr", "colourpicker")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "The PhenoMapR Shiny app needs these R packages installed: ",
      paste(missing, collapse = ", "),
      ". Install them with: install.packages(c(\"",
      paste(missing, collapse = "\", \""), "\"))."
    )
  }
  # Attach the shiny / bslib / DT / ggplot2 / dplyr packages, but
  # SKIP any that are already on the search path. This sidesteps a
  # nasty environment-specific reload trap reported in the wild:
  #
  #   "Loading required package: shiny
  #    Failed with error: 'Package 'shiny' version 1.8.1.1 cannot
  #    be unloaded: namespace 'shiny' is imported by 'miniUI',
  #    'Seurat', 'AUCell' so cannot be unloaded'"
  #
  # When the user has already attached shiny indirectly (via Seurat,
  # AUCell, miniUI...) and an on-disk version mismatch causes
  # library(shiny) to attempt a reload, that reload chain fails
  # because shiny is locked by its importers. Simply not re-attaching
  # what is already attached avoids the reload entirely. The
  # packages we need are already loaded into the search path by the
  # earlier requireNamespace() calls above and remain fully usable.
  for (pkg in required) {
    if (!paste0("package:", pkg) %in% search()) {
      do.call(library, list(pkg, character.only = TRUE))
    }
  }
})

if (!requireNamespace("PhenoMapR", quietly = TRUE)) {
  stop(
    "PhenoMapR is not installed in this R library. ",
    "Install from GitHub:\n  remotes::install_github(\"brooksbenard/PhenoMapR\")\n",
    "or from a local clone:\n  remotes::install_local(\"/path/to/PhenoMapR\")."
  )
}

# Source local helpers (parsers, validators, %||%)
local_dir <- if (nzchar(Sys.getenv("PHENOMAPR_SHINY_DIR"))) {
  Sys.getenv("PHENOMAPR_SHINY_DIR")
} else if (file.exists("helpers.R")) {
  getwd()
} else {
  system.file("shiny", package = "PhenoMapR")
}
source(file.path(local_dir, "helpers.R"), local = TRUE)
if (!nzchar(Sys.getenv("PHENOMAPR_SHINY_DEMO_RDS", unset = ""))) {
  demo_rds <- normalizePath(
    file.path(local_dir, "..", "extdata", "shiny", "PAAD_CRA001160_demo_5000.rds"),
    winslash = "/",
    mustWork = FALSE
  )
  if (file.exists(demo_rds)) {
    Sys.setenv(PHENOMAPR_SHINY_DEMO_RDS = demo_rds)
  }
}

# Allow uploads of any size by default. Shiny's default 5 MB cap is too small
# for sc / spatial RDS / h5ad files, and PhenoMapR users typically run this
# app on a workstation or HPC node where the OS / RAM is the only real limit.
# Override with options(shiny.maxRequestSize = <bytes>) before launching if
# you want to cap uploads (e.g. on a shared web server).
if (is.null(getOption("shiny.maxRequestSize"))) {
  options(shiny.maxRequestSize = Inf)
}

# ============================================================================
# Reference catalogue + cancer-type pre-load (cheap)
# ============================================================================
REFERENCE_CHOICES <- c(
  "PRECOG"           = "precog",
  "TCGA"             = "tcga",
  "Pediatric PRECOG" = "pediatric_precog",
  "ICI-PRECOG"       = "ici_precog",
  "Custom"           = "_custom"
)
.cancer_types_cache <- new.env(parent = emptyenv())
get_cancer_types <- function(ref) {
  if (!exists(ref, envir = .cancer_types_cache, inherits = FALSE)) {
    assign(ref, safe_list_cancer_types(ref), envir = .cancer_types_cache)
  }
  get(ref, envir = .cancer_types_cache, inherits = FALSE)
}

# ============================================================================
# UI
# ============================================================================

# Reusable help icons for parameters that need tooltips.
help_icon <- function(text) {
  tags$span(
    class = "help-tip",
    title = text,
    tags$i(class = "bi bi-info-circle", style = "color:#56688a;")
  )
}

ui <- page_navbar(
  title = tagList(
    # Brand: large "PhenoMapR" wordmark followed by the hex logo. The brand
    # is wrapped in a clickable <a> that bumps `input$nav_to_welcome` via
    # Shiny.setInputValue; the server-side observer calls
    # `bslib::nav_select()` to switch the navbar to the Welcome panel.
    # `return false` prevents the default `href="#"` from scrolling the
    # page. The Welcome tab itself is hidden from the navbar tab strip
    # (see styles.css `.nav-welcome-hidden`), so the brand is the only
    # way to get back to it.
    tags$a(
      class = "navbar-brand-link",
      href = "#",
      title = "Back to overview",
      onclick = paste(
        "Shiny.setInputValue('nav_to_welcome', Math.random(),",
        "{priority: 'event'}); return false;"
      ),
      tags$span("PhenoMapR", class = "brand-title"),
      tags$img(src = "PhenoMapR_logo.png",
               class = "brand-logo",
               alt = "PhenoMapR hex logo")
    )
  ),
  id = "main_nav",
  fillable = FALSE,
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#264653",
    secondary = "#2A9D8F",
    base_font = font_google("Source Sans Pro"),
    heading_font = font_google("Source Sans Pro")
  ),
  header = tags$head(
    # Cache-bust styles.css on every app boot so users don't have to
    # hard-refresh after we ship CSS-only tweaks (button positioning,
    # colors, etc). The query string is the file's mtime in ms --
    # changes whenever styles.css does, so the browser fetches a
    # fresh copy after each edit but caches normally between sessions
    # when the file is unchanged.
    tags$link(rel = "stylesheet",
              href = local({
                css_path <- file.path(
                  if (nzchar(Sys.getenv("PHENOMAPR_SHINY_DIR"))) {
                    Sys.getenv("PHENOMAPR_SHINY_DIR")
                  } else if (file.exists("www/styles.css")) {
                    getwd()
                  } else {
                    system.file("shiny", package = "PhenoMapR")
                  },
                  "www", "styles.css")
                v <- tryCatch(
                  as.integer(file.info(css_path)$mtime),
                  error = function(e) as.integer(Sys.time())
                )
                if (is.na(v)) v <- as.integer(Sys.time())
                sprintf("styles.css?v=%d", v)
              })),
    tags$link(rel = "stylesheet",
              href = "https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.css"),
    # Inject the centered busy-overlay markup + custom-message handlers
    # used by phenomapr_busy_show()/phenomapr_busy_hide().
    phenomapr_busy_assets()
  ),

  # -------------------------------------------------------------------------
  # Welcome
  # -------------------------------------------------------------------------
  nav_panel(
    # The .nav-welcome-hidden marker class lets styles.css hide this nav
    # link from the navbar tab strip without disabling the panel itself
    # (the brand-link above is the only way back to this panel now).
    title = tagList(tags$span(class = "nav-welcome-hidden",
                              icon("house"), " Welcome")),
    value = "welcome",
    card(
      card_header(tags$strong("Overview of PhenoMapR Shiny App")),
      card_body(
        layout_columns(
          col_widths = c(7, 5),
          tagList(
            tags$div(
              class = "welcome-intro",
              tags$p(
                tags$strong("PhenoMapR"),
                " is a scaleable approach to map sample phenotypes associated ",
                "with bulk expression onto single-cell, spatial, bulk, and ",
                "pseudobulk transcriptomics data. PhenoMapR comes pre-loaded ",
                "with the largest database of pan-cancer prognostic meta ",
                "z-scores (",
                tags$strong(class = "ref-name ref-name-precog", "PRECOG"),
                ", ",
                tags$strong(class = "ref-name ref-name-pediatric",
                            "Pediatric PRECOG"),
                ", ",
                tags$strong(class = "ref-name ref-name-ici", "ICI PRECOG"),
                ", & ",
                tags$strong(class = "ref-name ref-name-tcga", "TCGA"),
                "). ",
                tags$a(
                  class = "welcome-intro-ref-link",
                  href = "https://precog.stanford.edu/",
                  target = "_blank",
                  rel = "noopener noreferrer",
                  icon("link"),
                  " PRECOG database (precog.stanford.edu)"
                ),
                " ",
                tags$a(
                  class = "welcome-intro-ref-link",
                  href = paste0(
                    "https://academic.oup.com/nar/article/54/D1/D1579/",
                    "8324954"
                  ),
                  target = "_blank",
                  rel = "noopener noreferrer",
                  icon("book-open"),
                  " PRECOG paper (NAR 2026)"
                )
              ),
              tags$p(
                "If you are interested in mapping non-cancer phenotypes ",
                "(e.g. Alzheimer's risk), PhenoMapR will take a custom bulk ",
                "expression reference you supply and generate a phenotype ",
                "signature. This signature can then be used to score bulk, ",
                "single-cell, or spatial transcriptomics data."
              )
            ),

            # "How to use this app" — card-style panel with a numbered
            # stepper. Each step's headline is large + brand-primary so it
            # reads as the section's "main part," while the body text stays
            # in a softer color for visual hierarchy.
            tags$div(
              class = "welcome-section welcome-howto",
              tags$div(
                class = "welcome-section-header",
                icon("compass"),
                tags$span("How to use this app")
              ),
              tags$div(
                class = "welcome-section-body",
                tags$ol(
                  class = "welcome-steps",
                  tags$li(
                    class = "welcome-step",
                    tags$div(class = "welcome-step-num", "1"),
                    tags$div(
                      class = "welcome-step-text",
                      tags$div(class = "welcome-step-title", "Upload data"),
                      tags$div(
                        class = "welcome-step-desc",
                        "PhenoMapR can process most of the common gene ",
                        "expression file formats (e.g. RDS / Seurat object, ",
                        "SingleCellExperiment, AnnData / h5ad, 10x HDF5 ",
                        "(.h5), data.frame, matrix etc.). Additional tabular ",
                        "cell metadata can be loaded as well (helpful when an ",
                        "expression format does not contain, for example, ",
                        "cell type annotations)."
                      )
                    )
                  ),
                  tags$li(
                    class = "welcome-step",
                    tags$div(class = "welcome-step-num", "2"),
                    tags$div(
                      class = "welcome-step-text",
                      tags$div(class = "welcome-step-title",
                               "Choose a phenotype signature"),
                      tags$div(
                        class = "welcome-step-desc",
                        "Pick a built-in cohort meta-z (",
                        tags$strong(class = "ref-name ref-name-precog", "PRECOG"),
                        " / ",
                        tags$strong(class = "ref-name ref-name-tcga", "TCGA"),
                        " / ",
                        tags$strong(class = "ref-name ref-name-pediatric",
                                    "Pediatric PRECOG"),
                        " / ",
                        tags$strong(class = "ref-name ref-name-ici", "ICI PRECOG"),
                        "), upload a precomputed signature file, or ",
                        tags$em("derive"),
                        " a custom signature from your own bulk expression ",
                        "+ phenotype on the ",
                        tags$strong("Phenotype"),
                        " tab (binary, continuous, or survival outcomes)."
                      )
                    )
                  ),
                  tags$li(
                    class = "welcome-step",
                    tags$div(class = "welcome-step-num", "3"),
                    tags$div(
                      class = "welcome-step-text",
                      tags$div(class = "welcome-step-title",
                               "Score samples, cells, or locations"),
                      tags$div(
                        class = "welcome-step-desc",
                        "Compute PhenoMapR scores for each cell / sample / ",
                        "location and view high-level score distribution ",
                        "statistics."
                      )
                    )
                  ),
                  tags$li(
                    class = "welcome-step",
                    tags$div(class = "welcome-step-num", "4"),
                    tags$div(
                      class = "welcome-step-text",
                      tags$div(class = "welcome-step-title",
                               "Visualization"),
                      tags$div(
                        class = "welcome-step-desc",
                        "Overlay PhenoMapR score, cell type, source, or ",
                        "phenotype-group label onto any 2D embedding stored ",
                        "on your object (UMAP, tSNE, PCA, or tissue / spot ",
                        "coordinates for spatial inputs)."
                      )
                    )
                  ),
                  tags$li(
                    class = "welcome-step",
                    tags$div(class = "welcome-step-num", "5"),
                    tags$div(
                      class = "welcome-step-text",
                      tags$div(class = "welcome-step-title", "Marker genes"),
                      tags$div(
                        class = "welcome-step-desc",
                        "Run differential expression for the ",
                        tags$span(class = "phenomapr-phenotype-minus",
                                  HTML("&ldquo;most phenotype &minus;&rdquo;")),
                        " vs. ",
                        tags$span(class = "phenomapr-phenotype-plus",
                                  HTML("&ldquo;most phenotype +&rdquo;")),
                        " tails (globally or per cell type) and draw a ",
                        "ComplexHeatmap of the top markers per group."
                      )
                    )
                  )
                )
              )
            )
            # NOTE: the former "Server / remote use" collapsible card
            # was removed at the user\'s request -- the
            # PhenoMapR::run_app(host = "0.0.0.0", port = 3838) recipe
            # remains documented in the package README and the
            # ?run_app help page, which are the appropriate homes for
            # deployment instructions. Keeping the welcome page focused
            # on the in-app workflow (How it works) makes the
            # landing screen tidier for the typical local-laptop user.
          ),
          tags$div(
            class = "welcome-image",
            tags$img(
              src = "PhenoMapR_stacked_modalities.png",
              alt = "PhenoMapR scoring across bulk, single-cell, and spatial modalities",
              style = "max-width: 100%; height: auto; border-radius: 6px; box-shadow: 0 2px 12px rgba(0,0,0,0.08);"
            ),
            tags$div(
              class = "welcome-image-caption",
              "PhenoMapR maps bulk-derived prognostic signatures onto matched single-cell and spatial data."
            )
          )
        )
      )
    ),
    div(class = "version-note",
        sprintf("PhenoMapR %s · Shiny %s",
                as.character(packageVersion("PhenoMapR")),
                as.character(packageVersion("shiny"))))
  ),

  # -------------------------------------------------------------------------
  # 1. Upload data
  # -------------------------------------------------------------------------
  nav_panel(
    title = tagList(
      tags$span(class = "nav-step", icon("upload"), " 1. Data"),
      tags$span(class = "nav-step-arrow", "\u2192")
    ),
    value = "upload",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        # ---- Expression input sub-section ----
        # Wrapped in `.phenomapr-compact-stack expression-compact-stack`
        # so the file picker, "Use a tiny demo matrix instead"
        # actionLink, and the matrix-diagnostics renderUI sit close
        # together (the renderUI is often empty for object inputs and
        # would otherwise leave a wide blank gap under the
        # actionLink before the <hr/> divider).
        tags$div(
          class = "phenomapr-compact-stack expression-compact-stack",
          h4("Expression input"),
          helpText(
            "Upload your bulk, single-cell, or spatial dataset that you want scored."
          ),
          phenomapr_file_input(
            "expr_file",
            label = NULL,
            accept = c(".rds", ".h5", ".h5ad", ".tsv", ".csv", ".txt"),
            width = "100%"
          ),
          actionLink("use_demo", "Use demo dataset instead",
                     icon = icon("flask")),
          uiOutput("demo_source_panel"),
          # Matrix-only diagnostics panel: gene-ID style (HUGO / ENSG /
          # mixed), expression format (raw counts / CPM/TPM / log-norm /
          # z-scaled), and a single-cell vs bulk guess via sparsity. The
          # `uiOutput` is rendered only when state$expr_summary$kind is
          # "matrix" (see the server-side renderUI further below). When
          # cleanup is recommended the same panel also exposes
          # checkboxes + a "Clean & normalize" button that runs
          # clean_matrix_input() in place on state$expression.
          uiOutput("expr_matrix_diagnostics")
        ),

        hr(),
        # ---- Cell metadata sub-section ----
        # Wrapped in `.phenomapr-compact-stack` so the three column
        # dropdowns (cell-ID / cell-type / source) and the surrounding
        # helper UI (status, upload <details>, file picker) all share
        # the tightened sidebar rhythm. Without this wrapper the
        # default Shiny .form-group spacing leaves a large vertical
        # gap between dropdowns that pushes the "About input data"
        # card below the fold on shorter viewports.
        tags$div(
          class = "phenomapr-compact-stack metadata-compact-stack",
          h4("Cell metadata"),
          uiOutput("metadata_status"),
          # Metadata upload section -- single static <details> wrapper
          # housing the dynamic <summary> label, a short helpText, and
          # the file picker. Stable DOM identity is preserved by
          # NEVER re-rendering the wrapper or the file picker via
          # renderUI: only the <summary> *content* is server-rendered
          # through `metadata_upload_summary`, which is replaced
          # surgically without disturbing the file picker below it.
          # The `open` attribute (auto-open when no metadata is
          # present, auto-close once metadata arrives) is toggled via
          # a custom message handler ("phenomapr-set-details-open")
          # registered in helpers.R.
          tags$details(
            id    = "metadata_upload_details",
            class = "metadata-upload-details",
            # Custom container = tags$summary so the disclosure
            # trigger is a real HTML5 <summary> as a direct child of
            # the <details>. The renderUI body supplies just the
            # summary's text/icon children; the <summary> itself (id
            # + class) is stable across re-renders, which means the
            # collapsed/open state controlled by the wrapper's `open`
            # attribute is never disturbed.
            uiOutput(
              "metadata_upload_summary",
              container = function(...)
                tags$summary(class = "metadata-upload-details-summary", ...)
            ),
            tags$div(
              class = "metadata-upload-help",
              helpText(
                "Optional \u2014 only needed for matrix uploads, or to ",
                "attach extra annotations on top of an object's built-in ",
                "metadata. A cell-ID column must match the expression ",
                "matrix columns."
              )
            ),
            tags$div(
              class = "metadata-upload-file-host",
              phenomapr_file_input(
                "meta_file",
                label  = NULL,
                accept = c(".rds", ".tsv", ".csv", ".txt"),
                width  = "100%"
              )
            )
          ),
          selectInput("meta_cell_id_col", "Cell ID column", choices = NULL),
          selectInput("meta_cell_type_col", "Cell type column (optional)",
                      choices = NULL),
          selectInput("meta_source_col", "Source / group column (optional)",
                      choices = NULL)
        )
      ),
      card(
        fill = FALSE,
        card_body(
          padding = 0,
          tags$details(
            class = "phenomapr-about-markers-details",
            tags$summary(
              icon("circle-info"),
              tags$span(
                class = "phenomapr-about-markers-summary-label",
                " About input data"
              ),
              tags$span(
                class = "phenomapr-about-markers-summary-hint",
                " \u2014 click to expand"
              )
            ),
            layout_columns(
              col_widths = c(4, 4, 4),
              gap = "1rem",
              tags$div(
                class = "data-input-note",
                tags$div(class = "data-input-note-title",
                         icon("file-import"), " PhenoMapR accepts"),
                markdown(
                  "* **Matrix / data.frame** \u2014 rows = genes, columns = cells/samples.
                * **Seurat** (v4 / v5; assays `RNA`, `Spatial`).
                * **SingleCellExperiment** / **SpatialExperiment**.
                * **AnnData** (via `reticulate`)."
                )
              ),
              tags$div(
                class = "data-input-note",
                tags$div(class = "data-input-note-title",
                         icon("dna"), " Gene IDs"),
                markdown(
                  "For matrix uploads, gene IDs must be HUGO symbols
                (e.g. `TP53`). ENSG IDs are detected and flagged;
                convert with `biomaRt` / `AnnotationDbi` before
                uploading."
                )
              ),
              tags$div(
                class = "data-input-note",
                tags$div(class = "data-input-note-title",
                         icon("server"), " Upload size"),
                markdown(
                  "By default this app accepts uploads of any size \u2014
                analyses should work as long as your machine has
                enough RAM to process them. Cap the limit with
                `options(shiny.maxRequestSize = <bytes>)` before
                launching if needed."
                )
              )
            )
          )
        )
      ),
      # ---- Dataset overview --------------------------------------------------
      card(
        card_header(icon("table-list"), " Dataset overview"),
        card_body(
          uiOutput("dataset_overview_summary"),
          layout_columns(
            col_widths = c(6, 6),
            card(
              phenomapr_card_header_modal_dl(
                tags$strong("Cells per cell type"),
                panel_id = "celltype_count_plot"
              ),
              card_body(plotOutput("celltype_count_plot", height = "260px"))
            ),
            card(
              phenomapr_card_header_modal_dl(
                tags$strong("Cells per source / group"),
                panel_id = "source_count_plot"
              ),
              card_body(plotOutput("source_count_plot", height = "260px"))
            )
          ),
          layout_columns(
            col_widths = c(6, 6),
            card(
              phenomapr_card_header_modal_dl(
                tags$strong("Cell type × source composition"),
                panel_id = "celltype_source_plot"
              ),
              card_body(plotOutput("celltype_source_plot", height = "260px"))
            ),
            card(
              phenomapr_card_header_dl(
                tags$strong("Metadata columns"),
                download_id = "metadata_columns_tbl_download"
              ),
              card_body(DTOutput("metadata_columns_tbl"))
            )
          )
        )
      )
    )
  ),

  # -------------------------------------------------------------------------
  # 2. Phenotype
  # -------------------------------------------------------------------------
  nav_panel(
    title = tagList(
      tags$span(class = "nav-step", icon("pen-nib"), " 2. Phenotype"),
      tags$span(class = "nav-step-arrow", "\u2192")
    ),
    value = "reference",  # internal value preserved for stability
    layout_sidebar(
      # fillable = FALSE prevents bslib from stretching cards to fill the
      # main column, so the small "Phenotype signature status" / "Gene
      # coverage" cards hug their content and the signature plot sits right
      # below them instead of being pushed down by trailing whitespace.
      fillable = FALSE,
      sidebar = sidebar(
        width = 360,
        # Top-level signature source. The radio group is wrapped in a
        # `.reference-pills` div so styles.css can re-skin each option as a
        # colored click-target widget (the underlying input is still a
        # shiny radio so server-side `input$reference_choice` semantics are
        # unchanged). The label uses a plain `h4` to match the sizing of
        # the other sidebar section headers ("Expression input", "Cell
        # metadata", "PhenoMap() parameters", …).
        h4("Phenotype signature source"),
        tags$div(
          class = "reference-pills",
          radioButtons("reference_choice", label = NULL,
                       choices = REFERENCE_CHOICES, selected = "precog")
        ),
        conditionalPanel(
          "input.reference_choice != '_custom'",
          selectInput("cancer_type", "Select Cancer / tissue type",
                      choices = NULL)
        ),
        # Signature |z| cutoff slider. Lives OUTSIDE the
        # built-in-only conditionalPanel so it remains visible (and
        # therefore wired up) for custom phenotype signatures as
        # well; the Phenotype-signature histogram below uses
        # input$z_score_cutoff to draw live vertical lines on either
        # side of zero, mirroring the cutoff overlay the score
        # distribution histogram already does for percentile
        # thresholds. PhenoMap() honours this cutoff for every
        # reference (not just built-ins), so showing the slider for
        # custom signatures matches what the scoring step does
        # downstream.
        sliderInput("z_score_cutoff",
                    "Set threshold for Signature |z| cutoff",
                    min = 0, max = 5, value = 2, step = 0.1),
        conditionalPanel(
          "input.reference_choice != '_custom'",
          # Pre-filter note. Built-in signatures ship inside the
          # package already filtered to keep the install lightweight,
          # so the slider can only tighten the cutoff further (values
          # below the pre-filter are no-ops on the data that ships
          # with PhenoMapR). The pre-filter is per-reference:
          #   PRECOG / TCGA / Pediatric PRECOG --> |z| >= 2
          #   ICI PRECOG                       --> |z| >= 1
          # ICI is relaxed to 1 because several ICI columns drop to
          # zero finite genes once trimmed at 2 (the package would
          # otherwise be unable to score those signatures at all);
          # the looser cutoff keeps every ICI column usable while
          # only adding ~1 MB to the install. Two nested
          # conditionalPanels keep the message reference-aware
          # without burying the class string ("phenomapr-prefilter-
          # note") inside a renderUI -- the smoke test still finds
          # it inline between the slider and the custom-source
          # conditionalPanel.
          tags$div(
            class = "phenomapr-prefilter-note",
            icon("info-circle"),
            conditionalPanel(
              "input.reference_choice == 'ici_precog'",
              HTML(paste0(
                " ICI PRECOG ships pre-filtered to ",
                "<strong>|z| &ge; 1</strong> (every other built-in ",
                "ships at <strong>|z| &ge; 2</strong>). Selecting ICI ",
                "auto-sets the slider to 1 because several cohorts have ",
                "no genes above 2. The slider can only tighten this ",
                "cutoff further; values below 1 are no-ops on the ",
                "shipped ICI matrix."
              ))
            ),
            conditionalPanel(
              "input.reference_choice != 'ici_precog'",
              HTML(paste0(
                " This built-in phenotype signature ships ",
                "pre-filtered to <strong>|z| &ge; 2</strong> to keep ",
                "the package lightweight. The slider can only ",
                "tighten this cutoff further; values below 2 are ",
                "no-ops. (ICI PRECOG is the exception and ships at ",
                "<strong>|z| &ge; 1</strong>)"
              ))
            )
          )
        ),
        conditionalPanel(
          "input.reference_choice == '_custom'",
          # Surface the two custom paths as a clearly labelled radio. The
          # "Derive…" branch is what generates a custom phenotype signature
          # right inside the app from a user's own bulk + phenotype data.
          radioButtons(
            "custom_source", "Custom signature method",
            choices = c(
              "Derive a phenotype signature from bulk + phenotype" = "derive",
              "Upload a precomputed signature file"                = "upload"
            ),
            selected = "derive"
          ),
          conditionalPanel(
            "input.custom_source == 'upload'",
            phenomapr_file_input(
              "custom_ref_file",
              "Signature file (.rds / .tsv / .csv)",
              accept = c(".rds", ".tsv", ".csv", ".txt"), width = "100%"
            ),
            helpText(
              "Gene symbols as row names, a single numeric column of z-scores."
            )
          ),
          conditionalPanel(
            "input.custom_source == 'derive'",
            # The inputs in this block are stacked tightly via the
            # `phenomapr-compact-stack` class -- ordinary Shiny `.form-group`
            # bottom-margins push the file pickers, ID / outcome dropdowns,
            # Phenotype-type and Binary-positive radios way apart, which
            # buried the actual "Derive" button below the fold on smaller
            # sidebars. See `phenomapr-compact-stack` in styles.css.
            helpText(
              "Upload bulk expression (genes × samples) and a sample-level ",
              "phenotype table. PhenoMapR will compute per-gene z-scores ",
              "(Cox for survival, logistic for binary, correlation for ",
              "continuous) that you can immediately use to score your ",
              "single-cell / spatial data."
            ),
            tags$div(
              class = "phenomapr-compact-stack",
              phenomapr_file_input("derive_bulk_file", "Bulk expression (genes × samples)",
                                   accept = c(".rds", ".tsv", ".csv", ".txt")),
              uiOutput("derive_bulk_summary"),
              radioButtons(
                "derive_platform", "Bulk expression platform",
                choices = c("RNA-seq" = "rnaseq", "Microarray" = "microarray"),
                selected = "rnaseq",
                inline = TRUE
              ),
              conditionalPanel(
                "input.derive_platform == 'microarray'",
                phenomapr_file_input(
                  "derive_probe_file",
                  "GPL probe annotation (probes → genes)",
                  accept = c(".tsv", ".csv", ".txt", ".rds")
                ),
                helpText(
                  "Required when row names are platform probe IDs. The table ",
                  "should include a probe ID column and a gene symbol column ",
                  "(GEO GPL format). Expression is quantile-normalized and ",
                  "gene-scaled following the original PRECOG microarray pipeline."
                )
              ),
              checkboxInput(
                "derive_meta_z_mode",
                "Combine multiple bulk cohorts into one meta-z signature",
                value = FALSE
              ),
              conditionalPanel(
                "input.derive_meta_z_mode",
                helpText(
                  "Upload bulk + phenotype for each cohort, then click ",
                  HTML('"Add cohort to meta-z list"'),
                  " before uploading the next. When all cohorts are added, ",
                  "click Derive phenotype signature to run a Stouffer ",
                  "meta-analysis across studies (like built-in PRECOG meta-z ",
                  "signatures)."
                ),
                actionButton(
                  "derive_add_study", "Add cohort to meta-z list",
                  icon = icon("plus"), class = "btn-outline-secondary"
                ),
                actionButton(
                  "derive_clear_studies", "Clear meta-z cohort list",
                  class = "btn-outline-secondary"
                ),
                verbatimTextOutput("derive_meta_z_study_status")
              ),
              phenomapr_file_input("derive_phen_file", "Phenotype table (rows = samples)",
                                   accept = c(".rds", ".tsv", ".csv", ".txt")),
              selectInput("derive_id_col", "Sample ID column", choices = NULL),
              selectInput("derive_pheno_col", "Outcome column", choices = NULL),
              selectInput(
                "derive_type", "Phenotype type",
                choices = c("auto", "binary", "continuous", "survival"),
                selected = "auto"
              ),
              uiOutput("derive_pheno_preview"),
              conditionalPanel(
                "input.derive_type == 'survival'",
                selectInput("derive_time_col", "Time column", choices = NULL),
                selectInput("derive_event_col", "Event column", choices = NULL)
              ),
              conditionalPanel(
                "input.derive_type == 'binary' || input.derive_type == 'auto'",
                radioButtons(
                  "derive_binary_positive", "Binary positive level",
                  choices = c("Second factor level (default)" = "second",
                              "First factor level"            = "first"),
                  selected = "second"
                ),
                helpText(
                  "Controls which level of the binary outcome is coded as ",
                  "y = 1 in the logistic regression — i.e. which class ",
                  "positive z-scores point to. Switch to *First* when the ",
                  "first level of your factor (e.g. \"Responder\", ",
                  "\"Mutated\") is the phenotype you care about."
                )
              )
            ),
            selectInput(
              "derive_hugo_species", "Gene symbol species",
              choices = c("human", "mouse"), selected = "human"
            ),
            checkboxInput(
              "derive_normalize",
              "Normalize counts to log2(CPM + 1) if input looks raw",
              value = TRUE
            ),
            actionButton("derive_run", "Derive phenotype signature",
                         icon = icon("flask"), class = "btn-primary mt-2"),
            br(), br(),
            downloadButton("download_derived_reference",
                           "Download phenotype signature (RDS)",
                           class = "btn-outline-primary")
          )
        )
        # NOTE: the former "Reference diagnostics" sidebar block (top-N
        # + direction inputs) has been moved into the "Top prognostic
        # genes" card body in the main area, so all reference-related
        # controls live with the table they affect instead of being
        # split between the sidebar and the body.
      ),
      card(
        # The .signature-card-body class lets us position the
        # "Further reading" box (.phenotype-signature-refs) absolutely in
        # the upper-right corner of the card body without it stealing
        # space from the bullet list below. Per-database name labels are
        # color-coded to match the sidebar's reference-source palette so
        # the bullet list visually maps onto the picker pills above.
        card_header(icon("circle-info"), " Choosing a phenotype signature"),
        card_body(
          class = "signature-card-body",
          # Compact "Further reading" callout pinned to the upper-right of
          # the card. Provides the canonical PRECOG paper + website links
          # without breaking the flow of the explainer text.
          tags$div(
            class = "phenotype-signature-refs",
            tags$div(class = "phenotype-signature-refs-title",
                     "Further reading"),
            tags$a(
              href = "https://academic.oup.com/nar/article/54/D1/D1579/8324954",
              target = "_blank", rel = "noopener noreferrer",
              icon("book-open"), " PRECOG paper (NAR 2026)"
            ),
            tags$a(
              href = "https://precog.stanford.edu/",
              target = "_blank", rel = "noopener noreferrer",
              icon("link"), " precog.stanford.edu"
            )
          ),
          tags$p(tags$strong("Built-in signatures"),
                 " are meta-z-scores across many cohorts:"),
          tags$ul(
            class = "phenotype-signature-list",
            tags$li(
              tags$strong(class = "ref-name ref-name-precog", "PRECOG"),
              " - ", tags$strong("166 studies"),
              ", ~", tags$strong("18,000 adult patients"),
              " across ", tags$strong("39 cancer histologies"), "."
            ),
            tags$li(
              tags$strong(class = "ref-name ref-name-tcga", "TCGA"),
              " - single-cohort prognostic z, ~",
              tags$strong("11,000 patients"),
              " across ", tags$strong("33 adult cancer types"), "."
            ),
            tags$li(
              tags$strong(class = "ref-name ref-name-pediatric",
                          "Pediatric PRECOG"),
              " - ", tags$strong("32 studies"),
              ", ~", tags$strong("4,000 pediatric patients"),
              " across ", tags$strong("12 cancers"), "."
            ),
            tags$li(
              tags$strong(class = "ref-name ref-name-ici", "ICI PRECOG"),
              " - ", tags$strong("51 studies"),
              ", ~", tags$strong("4,000 immunotherapy-treated patients"),
              " across ", tags$strong("20 cancer subtypes"), "."
            )
          ),
          tags$p(
            tags$strong("Custom phenotype signature"),
            " is a ", tags$code("data.frame"),
            " with gene rownames and a single z-score column. Generate one ",
            tags$em("inside this app"),
            " from your own bulk expression + phenotype data via ",
            tags$em("Derive a phenotype signature"),
            ", or upload a precomputed file you've built elsewhere with ",
            tags$code("derive_reference_from_bulk()"),
            " (binary, continuous, or survival outcomes)."
          )
        )
      ),
      # Combined: phenotype signature status (left) + gene coverage
      # (right) inside a single card. The two halves used to live as
      # separate cards stacked side-by-side; merging them under one
      # card_header trims the redundant chrome (two card titles +
      # two card bodies + their borders) and lets the user read
      # "what reference is loaded, how restrictive is the cutoff, and
      # how much of my expression matrix overlaps it" as a single
      # status block. The download icon in the upper-right of the
      # header is wired to gene_coverage_tbl_download (same TSV
      # export as before) -- the signature-status side has no
      # download because it is computed live from the inputs.
      card(
        fill = FALSE,
        class = "compact-coverage-card phenomapr-signature-coverage-card",
        phenomapr_card_header_dl(
          tags$strong("Phenotype signature & gene coverage"),
          download_id = "gene_coverage_tbl_download",
          tooltip     = "Download gene coverage (TSV)"
        ),
        card_body(
          class = "compact-coverage-body",
          layout_columns(
            col_widths = c(6, 6),
            tags$div(
              class = "phenomapr-signature-coverage-half",
              tags$h6(
                "Phenotype signature status",
                class = "phenomapr-signature-coverage-subhead"
              ),
              uiOutput("reference_status")
            ),
            tags$div(
              class = "phenomapr-signature-coverage-half",
              tags$h6(
                "Gene coverage",
                class = "phenomapr-signature-coverage-subhead"
              ),
              helpText(
                "Fraction of your expression genes that overlap the chosen signature."
              ),
              DTOutput("gene_coverage_tbl")
            )
          )
        )
      ),
      # ---- Phenotype signature (2/3) + Top prognostic genes (1/3) ----------
      # The signature plot sits on the left at 2/3 width; the top
      # prognostic genes table is the right-hand companion at 1/3
      # width with its own top-N + direction controls embedded in the
      # card body (no more sidebar round-trip).
      layout_columns(
        col_widths = c(8, 4),
        # Both cards share an explicit height so the row reads as a
        # matched pair regardless of how DT lays out its pagination
        # chrome at runtime. The plot is fixed at 200 px and the DT
        # paginates to ~6 rows, both of which fit comfortably in the
        # 320 px envelope. The plot card gets fill = TRUE so its body
        # stretches to fill the card.
        card(
          height = "320px",
          phenomapr_card_header_modal_dl(
            tags$strong("Phenotype signature"),
            panel_id = "reference_signature_plot"
          ),
          card_body(
            class = "phenotype-signature-body",
            plotOutput("reference_signature_plot", height = "200px")
          )
        ),
        card(
          height = "320px",
          phenomapr_card_header_dl(
            tags$strong("Top prognostic genes"),
            download_id = "top_genes_tbl_download"
          ),
          card_body(
            class = "top-prognostic-genes-body",
            # Direction-only filtering was removed at the user's
            # request -- the DT shows the entire signature ordered by
            # |z| (direction = "both") and the user filters via the
            # built-in search box.
            DTOutput("top_genes_tbl")
          )
        )
      ),
      # ---- Available cancer types (full-width row beneath) -----------------
      # Lists the cancer / tissue types that ship with the currently-
      # selected built-in signature; populated server-side from
      # output$cancer_types_list. Full-width so long lists (PRECOG has
      # 39 cancers) wrap naturally instead of forcing horizontal scroll.
      card(
        card_header(
          icon("clipboard-list"),
          tags$strong(" Available cancer types"),
          tags$small(
            class = "text-muted",
            style = "margin-left: 0.5rem; font-weight: 400;",
            "Cancer / tissue types shipped with the selected built-in signature."
          )
        ),
        card_body(verbatimTextOutput("cancer_types_list"))
      ),
      # ---- Derived signature detail panel ------------------------------------
      conditionalPanel(
        "input.reference_choice == '_custom' && input.custom_source == 'derive'",
        card(
          card_header(icon("flask"), tags$strong(" Derived phenotype signature details")),
          card_body(
            helpText(
              "Appears after deriving a custom phenotype signature from your ",
              "bulk + phenotype upload. Lists the top-ranked positive and ",
              "negative z-score genes and lets you save the signature for ",
              "re-use."
            ),
            uiOutput("derived_reference_summary"),
            layout_columns(
              col_widths = c(6, 6),
              card(
                phenomapr_card_header_dl(
                  tagList(
                    tags$strong("Top positive ("),
                    phenomapr_phenotype_plus("phenotype +"),
                    tags$strong(") genes")
                  ),
                  download_id = "derived_top_pos_tbl_download"
                ),
                card_body(DTOutput("derived_top_pos_tbl"))
              ),
              card(
                phenomapr_card_header_dl(
                  tagList(
                    tags$strong("Top negative ("),
                    phenomapr_phenotype_minus("phenotype \u2212"),
                    tags$strong(") genes")
                  ),
                  download_id = "derived_top_neg_tbl_download"
                ),
                card_body(DTOutput("derived_top_neg_tbl"))
              )
            )
          )
        )
      )
    )
  ),

  # -------------------------------------------------------------------------
  # 3. Score
  # -------------------------------------------------------------------------
  nav_panel(
    title = tagList(
      tags$span(class = "nav-step", icon("chart-line"), " 3. Score"),
      tags$span(class = "nav-step-arrow", "\u2192")
    ),
    value = "score",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        # ---- Scoring summary -------------------------------------------------
        # Top-of-sidebar three-line context block: what kind of data the
        # user uploaded in step 1, which phenotype signature was picked in
        # step 2, and which cancer / tissue type that signature is being
        # interpreted against. Lets users confirm they are about to score
        # the right (data, signature, cancer-type) triple without having
        # to bounce back to earlier tabs.
        h4(tagList(icon("compass"), " Scoring summary")),
        uiOutput("score_data_status"),
        hr(),
        # ---- PhenoMap() parameters block --------------------------------
        # Wrapped in a `phenomapr-compact-stack score-params-compact-stack`
        # div so all children (heading, slot/assay controls, pseudobulk
        # toggle, action / download buttons) share the same tight
        # vertical rhythm and do not accrue Bootstrap's default
        # form-group margins around every input.
        tags$div(
          class = "phenomapr-compact-stack score-params-compact-stack",
          h4("Scoring parameters"),
          # Slot / assay controls. Only rendered for object-typed inputs
          # (Seurat / SCE / SpatialExperiment / AnnData), since plain
          # matrices and data.frames have no concept of an assay or a
          # named layer -- their values are scored directly. The
          # `score_show_slot_block` flag is driven by state$expr_summary
          # below (TRUE iff kind != "matrix"); using a conditionalPanel
          # rather than a uiOutput keeps the inputs in the DOM with
          # stable IDs so the observer further down can still call
          # updateRadioButtons / updateTextInput on them when an object
          # is loaded next.
          conditionalPanel(
            condition = "output.score_show_slot_block",
            # Hide the slot/assay rationale behind a collapsible
            # <details> chip labeled "Details" so the sidebar stays
            # uncluttered for users who already know what to score
            # against. The chevron (>) state is driven by HTML5; no
            # JavaScript needed.
            tags$details(
              class = "score-slot-details",
              tags$summary("Details"),
              helpText(
                icon("info-circle"),
                " PhenoMapR computes a weighted sum of expression x reference ",
                "z-score across signature genes. It expects ",
                tags$strong("log-normalized expression"),
                " (the 'data' / log-counts layer); choose 'counts (raw)' only ",
                "if your object lacks a normalized layer."
              )
            ),
            radioButtons(
              "score_slot", "Layer / slot to score against",
              choices = c("data (log-normalized)" = "data",
                          "counts (raw)"          = "counts"),
              selected = "data"
            ),
            textInput("score_assay",
                      "Assay (Seurat / SCE; blank = auto)",
                      value = "")
          ),
          conditionalPanel(
            condition = "output.score_allow_pseudobulk",
            helpText(
              icon("info-circle"),
              " For large single-cell or spatial datasets, enable ",
              tags$strong("Aggregate to pseudobulk"),
              " to sum expression within each patient, sample, cluster, or ",
              "tissue section before scoring (one score per group instead of ",
              "per cell/spot)."
            ),
            checkboxInput("pseudobulk", "Aggregate to pseudobulk", value = FALSE),
            conditionalPanel(
              "input.pseudobulk == true",
              selectInput("pseudobulk_group_by",
                          "Group cells by (metadata column)",
                          choices = NULL),
              helpText(
                "Pick the column whose values define each pseudobulk sample ",
                "(e.g. patient / donor / sample / cluster / tissue core / ",
                "spatial slide). Cells sharing that value are summed into one ",
                "pseudobulk profile before scoring."
              )
            )
          ),
          conditionalPanel(
            condition = "output.score_is_bulk_matrix",
            helpText(
              icon("info-circle"),
              " Bulk expression is already one profile per sample; ",
              "pseudobulking is not applicable."
            )
          ),
          actionButton("run_score", "Compute PhenoMapR scores",
                       icon = icon("play"), class = "btn-primary"),
          downloadButton("download_scores", "Download score table (TSV)",
                         class = "btn-outline-primary")
        ),
        hr(),
        h4("Phenotype groups"),
        helpText(
          HTML(paste0(
            "After scoring, the top / bottom percentile of cells are ",
            "automatically tagged as ",
            as.character(phenomapr_phenotype_plus()),
            " and ",
            as.character(phenomapr_phenotype_minus()),
            " as you drag the slider. These labels feed the marker-finding ",
            "step and downstream visualizations."
          ))
        ),
        sliderInput("percentile", "Tail percentile (top and bottom)",
                    min = 0.01, max = 0.40, value = 0.05, step = 0.01),
        downloadButton("download_groups", "Download labels (TSV)",
                       class = "btn-outline-primary")
      ),
      card(
        # `fill = FALSE` so the card sizes to its prose body instead of
        # stretching to fill the column height when its taller siblings
        # below (the Score-distribution + Cells-ordered cards) push the
        # layout. Without it, bslib's default flex behavior leaves a
        # large empty band beneath the markdown blurb -- the same fix
        # we already use on the Markers tab's About-markers card.
        fill = FALSE,
        card_header(icon("circle-info"), " About scoring"),
        card_body(
          markdown(
            "Each cell / sample gets a weighted sum of expression × reference
            z-score across genes that pass `|z| > z_score_cutoff`. Higher
            score = stronger phenotype-of-interest signal (the *direction*
            depends on the reference you chose). Use **Aggregate to
            pseudobulk** for very large single-cell datasets to score per
            sample / cluster / patient instead of per cell."
          )
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          phenomapr_card_header_modal_dl(
            tags$strong("Score distribution"),
            panel_id = "score_dist_plot"
          ),
          card_body(
            # Score scale toggle drives both the histogram data and its
            # title. "PhenoMapR score" uses the raw weighted-sum score
            # straight out of PhenoMap(); "Z-score" uses the same value
            # standardized (mean 0, sd 1) across the cells / samples so
            # users can see the distribution in a unit-free scale.
            radioButtons(
              "dist_score_scale", NULL,
              choices = c("PhenoMapR score" = "raw",
                          "Z-score"         = "scaled"),
              selected = "raw",
              inline = TRUE
            ),
            plotOutput("score_dist_plot", height = "300px")
          )
        ),
        card(
          phenomapr_card_header_modal_dl(
            tags$strong("Cells ordered by PhenoMapR score"),
            panel_id = "score_rank_plot"
          ),
          card_body(
            radioButtons(
              "rank_color_by", NULL,
              choices = c("Score"     = "score",
                          "Cell type" = "cell_type",
                          "Source"    = "source"),
              inline = TRUE, selected = "score"
            ),
            plotOutput("score_rank_plot", height = "320px")
          )
        )
      ),
      # ---- Phenotype groups (merged from the old standalone tab) ----------
      # NOTE: the two PLOTS that used to live here ("Score by cell type
      # and source" and "Per-cell-type group enrichment") have been
      # moved to the 4. Visualization tab so the Score tab can stay
      # focused on the score table + group-size summary, and the
      # visualization tab gathers all of the plot panels in one place.
      # The underlying server-side renderPlot / register_plot_download
      # wiring is unchanged -- only the UI placement moved.
      tags$div(
        class = "score-section-divider",
        tags$h3(tagList(icon("layer-group"), " Phenotype groups"))
      ),
      # About-tails + Group-sizes share a single row at 50/50 width so
      # the conceptual explanation and the live group counts sit next
      # to each other rather than stacking.
      layout_columns(
        col_widths = c(6, 6),
        card(
          # `fill = FALSE` -- prose body is short, so size to content
          # instead of stretching to whatever the taller sibling card
          # would be. See same comment on the "About scoring" card above.
          fill = FALSE,
          card_header(icon("circle-info"), " About phenotype tails"),
          card_body(
            tags$p(
              "PhenoMapR partitions cells into ",
              phenomapr_phenotype_plus(),
              " (top ", tags$em("N"), " %), ",
              phenomapr_phenotype_minus(),
              " (bottom ", tags$em("N"), " %), and ",
              tags$strong("Other"),
              " based on the chosen score column. These labels feed the ",
              "marker-finding step and any downstream visualizations."
            ),
            tags$p(
              "Lower the percentile to make the tails tighter (more ",
              "discriminative but fewer cells); raise it to include more cells ",
              "(statistically more robust but biologically fuzzier)."
            )
          )
        ),
        card(
          # Same rationale -- the Group-sizes table is exactly N+1 rows
          # (one per phenotype group), so we don't want it stretched
          # vertically. Pinning fill = FALSE here also matches the
          # About-tails card width-mate so both peers align cleanly.
          fill = FALSE,
          card_header(tags$strong("Group sizes")),
          card_body(uiOutput("group_summary"))
        )
      )
    )
  ),

  # -------------------------------------------------------------------------
  # 4. Visualization
  #
  # Generalized from "UMAP / embedding" to also handle spatial-transcriptomics
  # inputs. Spatial coordinates (Seurat `@images`, SpatialExperiment
  # `spatialCoords()`, AnnData `obsm["spatial"]`) surface as a synthetic
  # "spatial" reduction in the same dropdown, and the plot below applies
  # `coord_fixed()` + a reversed Y so spots are drawn in their tissue frame.
  # -------------------------------------------------------------------------
  nav_panel(
    title = tagList(
      tags$span(class = "nav-step", icon("braille"), " 4. Visualization"),
      tags$span(class = "nav-step-arrow", "\u2192")
    ),
    value = "umap",  # value kept for backward-compat with bookmarks/tests
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        # ---- Embedding selector + help sub-section ----
        # Wrapped in `.phenomapr-compact-stack embedding-compact-stack`
        # so the section header, Reduction dropdown, reduction-status
        # chip, spatial-sample selector, and the two collapsible
        # <details> blocks sit tight together above the <hr/> that
        # separates this group from the "Color cells by" controls
        # below. Without this wrapper the default Shiny .form-group
        # spacing + the <details> default margins leave a large
        # vertical gap pushing the Color/Point-size controls below
        # the fold on shorter viewports.
        tags$div(
          class = "phenomapr-compact-stack embedding-compact-stack",
          h4("Embedding"),
          # Collapse the long "how detection works" copy into a <details>
          # block so the sidebar isn't dominated by a wall of text by
          # default. The Reduction dropdown is the primary control and now
          # sits right under the section header.
          selectInput("umap_reduction", "Reduction", choices = NULL),
          uiOutput("umap_reduction_status"),
          # Tissue-section / core / library switcher. Surfaces ONLY when
          # the currently-selected embedding is spatial AND the loaded
          # object stamps >1 unique section label on its spots (Seurat
          # @images named entries, SpatialExperiment colData$sample_id,
          # AnnData obs[library_id]/etc). Server-side renderUI keeps the
          # control hidden otherwise so the sidebar isn't cluttered for
          # single-section datasets.
          uiOutput("spatial_sample_selector"),
          tags$details(
            class = "embedding-help-details",
            tags$summary("About reductions & auto-detection"),
            helpText(
              "Reductions stored on a Seurat / SingleCellExperiment / ",
              "AnnData object are picked up automatically. Spatial inputs ",
              "(Seurat with @images, SpatialExperiment, AnnData with ",
              "obsm['spatial']) also surface tissue coordinates as a ",
              "\"spatial\" reduction in the dropdown above. Imaging-based ",
              "spatial datasets that carry a cell-segmentation mask ",
              "(Seurat 5 FOV objects from Xenium / CosMx / MERSCOPE / ",
              "Visium HD; AnnData with obsm['segmentation']) additionally ",
              "expose a \"segmentation\" reduction \u2014 the per-cell ",
              "segmentation centroids drawn in their tissue frame, which ",
              "lines up 1:1 with the segmented cells in the underlying ",
              "image (whereas \"spatial\" plots the upstream pipeline's ",
              "canonical xy, often a spot grid). CosMx exports that only ",
              "store FOV pixel centroids in obs (e.g. x_FOV_px / y_FOV_px) ",
              "also surface \"segmentation\" automatically. When polygon ",
              "boundaries are present (Seurat 5 FOV @boundaries, or a ",
              "phenomapr_polygons table on a matrix), choose the ",
              "\"segmentation\" reduction and switch \"Segmentation render ",
              "style\" to Cell masks (polygons). For matrix uploads, the ",
              "app additionally auto-detects 2D coordinate column pairs ",
              "(e.g. UMAP_1/UMAP_2, tSNE_1/tSNE_2, X_umap_0/X_umap_1) on ",
              "the metadata table from the Data tab \u2014 no second ",
              "upload needed. Use the file picker below to supply a ",
              "separate embedding file if you don't have one in your ",
              "metadata."
            )
          ),
          tags$details(
            class = "embedding-help-details embedding-upload-details",
            tags$summary("Upload a separate embedding file"),
            helpText(
              "Only needed when the metadata table on the Data tab does not ",
              "already contain coordinate columns. Expects either an RDS / ",
              "TSV / CSV with a cell-ID column and two numeric coordinate ",
              "columns."
            ),
            phenomapr_file_input(
              "umap_upload",
              label = NULL,
              accept = c(".tsv", ".csv", ".txt", ".rds"), width = "100%"
            )
          )
        ),
        hr(),
        radioButtons(
          "umap_color_by", "Color cells by",
          choices = c("PhenoMapR score" = "score",
                      "Cell type"       = "cell_type",
                      "Source / group"  = "source",
                      "Phenotype group" = "group"),
          selected = "score"
        ),
        # Score-scale toggle: only relevant when the embedding is colored
        # by score. "PhenoMapR score" plots the raw weighted-sum score
        # straight out of PhenoMap(); "Z-score" plots the same value
        # standardized (mean 0, sd 1) across all cells / samples so a
        # 0-centered diverging palette matches the unit-free magnitudes.
        # Mirrors the same control on the Score tab so users can compare
        # the same statistic across tabs without re-deriving it.
        conditionalPanel(
          "input.umap_color_by == 'score'",
          radioButtons(
            "umap_score_scale", "Score scale",
            choices = c("PhenoMapR score" = "raw",
                        "Z-score"         = "scaled"),
            selected = "raw",
            inline = TRUE
          )
        ),
        sliderInput("umap_point_size", "Point size", min = 0.2, max = 3,
                    value = 0.8, step = 0.1),
        sliderInput("umap_point_alpha", "Point opacity", min = 0.1, max = 1,
                    value = 0.75, step = 0.05),
        conditionalPanel(
          "output.spatial_polygon_masks_available",
          radioButtons(
            "spatial_render_style", "Segmentation render style",
            choices = c("Points (centroids)" = "points",
                        "Cell masks (polygons)" = "masks"),
            selected = "points",
            inline = FALSE
          ),
          sliderInput("spatial_mask_alpha", "Mask fill opacity", min = 0.1, max = 1,
                      value = 0.85, step = 0.05)
        ),
        checkboxInput("umap_facet_source", "Facet by source", value = FALSE),
        downloadButton("download_umap_table", "Download embedding (TSV)",
                       class = "btn-outline-primary"),
        hr(),
        # ---- "Score by cell type and source" plot controls ----------------
        # These two controls live in the sidebar (rather than inline in
        # the boxplot card_body) so they share the same visual rhythm
        # as every other plot knob on this tab and don't eat vertical
        # space inside the card. They feed the second card on this tab
        # (`output$score_box_source_plot`):
        #
        #   * Source palette         -- Color (brand) vs Greyscale (B/W ramp)
        #   * Cell-type legend ncol  -- 1..8, spreads the legend across
        #     N columns so it stays on-canvas even with many cell types
        #     (e.g. 30+ in a PDAC atlas).
        tags$div(
          class = "phenomapr-score-box-controls",
          tags$h6("Score by cell type and source"),
          radioButtons(
            "score_box_source_palette",
            "Source palette",
            choices = c("Color" = "color", "Greyscale" = "greyscale"),
            selected = "color",
            inline = TRUE
          ),
          numericInput(
            "score_box_celltype_legend_ncol",
            "Cell-type legend columns",
            value = 1L, min = 1L, max = 8L, step = 1L
          )
        )
      ),
      card(
        phenomapr_card_header_modal_dl(
          tags$strong("Embedding"),
          panel_id = "umap_plot"
        ),
        card_body(plotOutput("umap_plot", height = "560px"))
      ),
      # ---- Score by cell type and source (moved from 3. Score tab) ----
      #
      # The Score-tab boxplot was relocated here so all plot panels live
      # in the Visualization tab. The renderPlot / phenomapr_register_
      # plot_download wiring keyed on "score_box_source_plot" is
      # unchanged; only this UI card moved.
      #
      # Gets BOTH a "Download plot" (modal) button AND a "Download plot
      # data" button. The data button exposes the same per-cell score
      # table that the standalone "Score table" panel used to surface
      # -- we removed that panel earlier (it duplicated the sidebar
      # "Download score table (TSV)" button and pushed the more useful
      # boxplot down below the fold), but kept the underlying data
      # download accessible from the boxplot header so users can grab
      # the rows the plot was built from without bouncing back to the
      # sidebar.
      card(
        phenomapr_card_header_modal_dl(
          tags$strong("Score by cell type and source"),
          panel_id = "score_box_source_plot",
          plot_label = "Download plot",
          tooltip = "Download plot (with options)",
          data_download_id = "score_table_download",
          data_label = "Download plot data",
          data_tooltip = "Download plot data (TSV, includes phenotype group labels)"
        ),
        card_body(
          helpText(
            "Boxplot of (scaled) PhenoMapR scores per cell type, ordered ",
            "from lowest to highest median score. With no source / group ",
            "column mapped, the plot annotates a single global ANOVA ",
            "F-test across cell types. With exactly 2 sources mapped, ",
            "Wilcoxon brackets compare sources within each cell type. ",
            "With 3+ sources mapped, a one-way ANOVA is run within each ",
            "cell type. Only significant (p < 0.05) brackets are drawn. ",
            "Source palette and cell-type legend columns are configured ",
            "in the sidebar."
          ),
          plotOutput("score_box_source_plot", height = "440px")
        )
      ),
      # ---- Per-cell-type group enrichment (moved from 3. Score tab) ----
      # Same migration rationale as the boxplot above. The downstream
      # download / render path is unchanged.
      card(
        phenomapr_card_header_modal_dl(
          tags$strong("Per-cell-type group enrichment"),
          panel_id = "group_by_celltype_plot"
        ),
        card_body(
          plotOutput("group_by_celltype_plot", height = "320px")
        )
      )
    )
  ),

  # -------------------------------------------------------------------------
  # 5. Markers
  # -------------------------------------------------------------------------
  nav_panel(
    title = tagList(
      tags$span(class = "nav-step nav-step-last", icon("dna"), " 5. Markers")
    ),
    value = "markers",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        # Sidebar heading -- mirrors the h4 banners we use in the
        # other tabs ("Cell metadata", "Phenotype groups", etc.) so
        # the sidebar reads as a labeled control surface rather than
        # a wall of widgets.
        h4("Marker parameters"),
        radioButtons(
          "marker_scope", "Scope",
          # Labels track the customer-facing phenotype rebrand:
          # "adverse / favorable" -> "phenotype + / phenotype \u2212"
          # mirrors the legend remap in score_rank_plot /
          # group_by_celltype_plot and the home-page hero blurb.
          #
          # Three modes, listed in order of *increasing stringency* so
          # the user reads them top-to-bottom from "broadest" to
          # "narrowest":
          #   1. "phenotype_groups": cohort-wide adverse-vs-favorable
          #      (cell type agnostic).
          #   2. "cell_type_vs_opposite": for each cell type in a tail,
          #      compare against ALL cells in the opposite tail
          #      (regardless of cell type). Mid-stringency -- still
          #      surfaces phenotype-driven markers for cell types that
          #      exist in only one tail (e.g. only adverse ductal, no
          #      favorable ductal), where the within-cell-type
          #      contrast below would return empty.
          #   3. "cell_type_specific": for each cell type, adverse vs
          #      favorable cells of the SAME cell type only. Most
          #      stringent. Returns nothing for a cell type when the
          #      same cell type does not exist in the opposite tail.
          choiceNames = list(
            tagList(
              "Cohort-wide (",
              phenomapr_phenotype_plus("phenotype +"),
              " vs ",
              phenomapr_phenotype_minus("phenotype \u2212"),
              ")"
            ),
            "Cell-type \u00D7 phenotype vs all opposite-tail cells",
            "Cell-type specific (within phenotype groups)"
          ),
          choiceValues = c(
            "phenotype_groups",
            "cell_type_vs_opposite",
            "cell_type_specific"
          ),
          selected = "phenotype_groups"
        ),
        sliderInput("marker_min_pct", "Min. fraction expressing (min.pct)",
                    min = 0, max = 0.5, value = 0.1, step = 0.05),
        sliderInput("marker_logfc", "Min. |log2FC|",
                    min = 0, max = 2, value = 0.25, step = 0.05),
        sliderInput("marker_pval", "p-value threshold",
                    min = 0.001, max = 0.1, value = 0.05, step = 0.001),
        numericInput("marker_maxcells", "Max cells per group",
                     value = 5000L, min = 100L, step = 500L),
        actionButton("run_markers", "Find markers",
                     icon = icon("microscope"), class = "btn-primary"),
        downloadButton("download_markers", "Download marker tables (RDS)",
                       class = "btn-outline-primary"),
        # ---- Heatmap drawing controls (moved into the sidebar) ----------
        # Wrapped in `.phenomapr-compact-stack marker-heatmap-compact-stack`
        # so the heatmap-specific controls (helpText, two numericInputs,
        # the Draw button) read as a tight grouping right under the
        # divider, instead of inheriting the default Shiny .form-group
        # spacing that would push the Draw button below the fold.
        hr(),
        tags$div(
          class = "phenomapr-compact-stack marker-heatmap-compact-stack",
          h4("Marker-gene heatmap"),
          helpText(
            "After running marker discovery, configure and draw a ",
            "ComplexHeatmap of the top markers per block."
          ),
          numericInput("hm_top_n", "Top genes per block", value = 20,
                       min = 5, max = 200, step = 5),
          numericInput("hm_n_labels", "Gene labels per block", value = 5,
                       min = 0, max = 50, step = 1),
          checkboxInput("hm_color_mark_labels",
                        "Color gene labels by cell type",
                        value = TRUE),
          # Block-outline controls. Only meaningful for
          # heatmap_type = "cell_type_specific" (the per-cell-type
          # heatmap is already split into one row slice per
          # cell-type * phenotype-bin block, with column extents
          # computed from the meta_df group/celltype columns -- see
          # marker_heatmap_args -> plot_phenotype_markers()). For the
          # cohort-wide heatmap this checkbox is a no-op.
          # Default OFF so the heatmap reads cleanly without
          # additional rules; flipping ON immediately re-renders the
          # heatmap with thick black outlines around each block,
          # which is what users actually want when they need to
          # explain "this block here".
          checkboxInput("hm_block_borders",
                        "Black borders around cell-type \u00D7 phenotype blocks",
                        value = FALSE),
          tags$details(
            class = "phenomapr-hm-colors-details",
            tags$summary(
              tags$strong("Color schemes"),
              tags$span(class = "phenomapr-hm-colors-summary-hint",
                        " \u2014 phenotype, score, cell types, heatmap")
            ),
            marker_heatmap_color_controls(
              "phenotype",
              "Phenotype groups",
              brewer_kind = "both",
              default_brewer = "RdBu",
              manual_ui = tagList(
                colourInput("hm_colors_phenotype_fav",
                            label = phenomapr_phenotype_minus("Most Phenotype -"),
                            "#2166AC"),
                colourInput("hm_colors_phenotype_other", "Other", "#F7F7F7"),
                colourInput("hm_colors_phenotype_adv",
                            label = phenomapr_phenotype_plus("Most Phenotype +"),
                            "#B2182B")
              )
            ),
            marker_heatmap_color_controls(
              "score",
              "PhenoMapR score bar",
              brewer_kind = "diverging",
              default_brewer = "RdBu",
              manual_ui = tagList(
                colourInput("hm_colors_score_low", "Low score", "#2166AC"),
                colourInput("hm_colors_score_mid", "Zero / midpoint", "#FFFFFF"),
                colourInput("hm_colors_score_high", "High score", "#B2182B")
              )
            ),
            marker_heatmap_color_controls(
              "celltype",
              "Cell types",
              brewer_kind = "qualitative",
              default_brewer = "Set2",
              manual_ui = textInput(
                "hm_colors_celltype_manual",
                "Cell type colors (comma-separated hex)",
                value = "",
                placeholder = "#2A9D8F, #E76F51, #264653, ..."
              )
            ),
            marker_heatmap_color_controls(
              "expr",
              "Scaled expression heatmap",
              brewer_kind = "diverging",
              default_brewer = "RdGy",
              manual_ui = tagList(
                colourInput("hm_colors_expr_low", "Low scaled expression",
                            "#1A1A1A"),
                colourInput("hm_colors_expr_mid", "Midpoint", "#FFFFFF"),
                colourInput("hm_colors_expr_high", "High scaled expression",
                            "#67001F")
              )
            )
          ),
          actionButton("draw_marker_heatmap", "Draw heatmap",
                       icon = icon("th"), class = "btn-primary")
        )
      ),
      card(
        # `fill = FALSE` keeps this short explanatory card from being
        # stretched to fill the column height when its taller siblings
        # below (the marker tables and the heatmap card) push the
        # layout. Without it, bslib's default flex behavior leaves a
        # large empty band beneath the markdown blurb.
        fill = FALSE,
        # Wrap the entire card body in <details>/<summary> so the user
        # gets the chip-style "About markers" header to click into when
        # they want a refresher, but the marker tables / heatmap below
        # take precedence in the default vertical real estate. Class
        # mirrors `.embedding-help-details` (sidebar's spiritual cousin)
        # so the chevron / spacing already feel native.
        card_body(
          padding = 0,
          tags$details(
            class = "phenomapr-about-markers-details",
            tags$summary(
              icon("circle-info"),
              tags$span(class = "phenomapr-about-markers-summary-label",
                        " About markers"),
              tags$span(
                class = "phenomapr-about-markers-summary-hint",
                " \u2014 click to expand"
              )
            ),
            layout_columns(
              col_widths = c(7, 5),
              gap = "1.25rem",
              # ---- Left column: simplified copy mapped to the schematic
              tags$div(
                class = "phenomapr-about-markers-text",
                tags$p(
                  "Wilcoxon DE (via ",
                  tags$code("presto"),
                  " when installed, otherwise base R; ",
                  tags$code("Seurat::FindMarkers"),
                  " for Seurat inputs in cohort-wide mode) between the ",
                  phenomapr_phenotype_plus(),
                  " and ",
                  phenomapr_phenotype_minus(),
                  " tails, listed here from broadest to most ",
                  "stringent."
                ),
                tags$p(
                  tags$strong("(1) Cohort-wide."),
                  " Whole ",
                  phenomapr_phenotype_plus(),
                  " tail vs whole ",
                  phenomapr_phenotype_minus(),
                  " tail \u2014 no cell-type partitioning."
                ),
                tags$p(
                  tags$strong("(2) Cell-type \u00D7 phenotype vs all opposite-tail cells."),
                  " One cell type in one tail vs ",
                  tags$em("every"),
                  " cell in the opposite tail. Still works when a cell type ",
                  "exists in only one tail (where mode 3 returns empty)."
                ),
                tags$p(
                  tags$strong("(3) Cell-type specific (within phenotype groups)."),
                  " One cell type in the ",
                  phenomapr_phenotype_plus(),
                  " tail vs the ",
                  tags$em("same"),
                  " cell type in the ",
                  phenomapr_phenotype_minus(),
                  " tail. Requires at least ",
                  tags$strong("5 cells"),
                  " of that type in ",
                  tags$strong("both"),
                  " phenotype tails."
                )
              ),
              # ---- Right column: matched schematic figure
              # Sourced from inst/figures/PhenoMapR_marker_schematic.png
              # and copied into inst/shiny/www/ on package build so Shiny
              # serves it with no addResourcePath() plumbing.
              tags$div(
                class = "phenomapr-about-markers-figure",
                tags$img(
                  src = "PhenoMapR_marker_schematic.png",
                  alt = paste(
                    "Schematic of the three marker-gene contrast scopes:",
                    "(1) cohort-wide Most Phenotype + vs Most Phenotype -,",
                    "(2) cell-type x phenotype vs all opposite-tail cells,",
                    "(3) cell-type specific within phenotype groups."
                  ),
                  loading = "lazy",
                  class = "phenomapr-about-markers-img"
                )
              )
            )
          )
        )
      ),
      navset_tab(
        nav_panel(
          title = tags$span(
            class = "phenomapr-phenotype-plus",
            "Phenotype + markers"
          ),
          card_body(
            phenomapr_panel_banner_dl("adverse_markers_tbl_download"),
            DTOutput("adverse_markers_tbl")
          )
        ),
        nav_panel(
          title = tags$span(
            class = "phenomapr-phenotype-minus",
            "Phenotype \u2212 markers"
          ),
          card_body(
            phenomapr_panel_banner_dl("favorable_markers_tbl_download"),
            DTOutput("favorable_markers_tbl")
          )
        )
      ),
      card(
        phenomapr_card_header_modal_dl(
          tags$strong("Marker-gene heatmap"),
          panel_id = "marker_heatmap"
        ),
        card_body(
          helpText(
            "Heatmap of the top markers from the ",
            phenomapr_phenotype_plus("phenotype +"),
            " and ",
            phenomapr_phenotype_minus("phenotype \u2212"),
            " tails. When markers were computed with cell-type-specific ",
            "scope, the heatmap shows one column slice per ",
            "phenotype-group \u00D7 cell-type block; otherwise it shows a ",
            "cohort-wide view. Use the controls in the sidebar to set the ",
            "number of genes per block, the number of labels, and to redraw."
          ),
          # The imageOutput is only mounted *after* the user clicks "Draw
          # heatmap" in the sidebar. Mounting it lazily prevents Shiny /
          # browser placeholders from rendering inside the empty
          # `.shiny-image-output` container before any PNG exists, which
          # was surfacing as a stray small teal "1" badge in some
          # browsers when the card sat empty. It also keeps the card
          # visually compact until the user actually requests a heatmap.
          #
          # ComplexHeatmap + anno_mark: rendered to a PNG file via
          # renderImage() rather than renderPlot() because Shiny replays
          # plot expressions on each device-size change, and anno_mark
          # labels do not survive those replays (lines redraw, but the
          # text glyphs come out blank). renderImage() serves the PNG
          # we draw once at a fixed DPI, which keeps gene labels intact.
          conditionalPanel(
            "input.draw_marker_heatmap > 0",
            imageOutput("marker_heatmap", height = "640px")
          ),
          conditionalPanel(
            "!input.draw_marker_heatmap",
            tags$div(
              class = "marker-heatmap-empty",
              tags$em(
                "Click ",
                tags$strong("Draw heatmap"),
                " in the sidebar to render the marker-gene heatmap once",
                " marker discovery has completed."
              )
            )
          )
        )
      )
    )
  ),

  nav_spacer(),
  # NOTE: the former "Plot appearance" navbar popover (corner radius
  # slider) has been folded into the plot-download modal alongside
  # the other export-time aesthetic controls (base font size, width,
  # height, DPI, file format). On-screen plots use the default
  # radius; per-export tweaks happen in the modal so users can iterate
  # against the live preview without juggling two UI surfaces.
  nav_item(
    tags$a(
      href = "https://github.com/brooksbenard/PhenoMapR",
      target = "_blank",
      tagList(icon("github"), " GitHub")
    )
  ),
  # End-session control: opens a confirmation modal that calls stopApp() so
  # the R process can release memory and return control to the console / shut
  # down the deployed app cleanly. See the matching observers below.
  nav_item(
    tags$button(
      id = "end_session",
      type = "button",
      class = "btn btn-sm btn-outline-danger end-session-btn action-button",
      title = "End the PhenoMapR session and stop the app",
      tagList(icon("power-off"), " End session")
    )
  )
)


# ============================================================================
# Server
# ============================================================================
server <- function(input, output, session) {

  # ------------------------------------------------------------------------
  # Brand click -> jump to Welcome panel
  # ------------------------------------------------------------------------
  # The PhenoMapR wordmark + hex logo in the navbar is a clickable link
  # that bumps `input$nav_to_welcome` (see the `title = tagList(...)` in
  # `page_navbar` above). The Welcome nav-link itself is hidden from the
  # tab strip, so this observer is the only path back to the overview.
  observeEvent(input$nav_to_welcome, {
    bslib::nav_select(id = "main_nav", selected = "welcome",
                      session = session)
  })

  # ------------------------------------------------------------------------
  # End-session control (navbar "End session" button)
  # ------------------------------------------------------------------------
  # Two-step shutdown: the navbar button opens a confirmation modal, and
  # only after the user clicks the modal's "End session" button do we call
  # shiny::stopApp(). This frees the R process so memory used by the loaded
  # Seurat / SCE / AnnData object is reclaimed, and on local installs it
  # returns control to the console.
  observeEvent(input$end_session, {
    showModal(modalDialog(
      title = tagList(icon("power-off"), " End the PhenoMapR session?"),
      tags$p(
        "This will stop the Shiny app and release any memory used by the ",
        "currently loaded expression object, scoring results, and marker ",
        "outputs."
      ),
      tags$p(
        tags$strong("Anything not yet exported will be lost."),
        " You can relaunch the app at any time with ",
        tags$code("PhenoMapR::run_app()"), "."
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("end_session_confirm", "End session",
                     icon = icon("power-off"),
                     class = "btn-danger")
      ),
      easyClose = TRUE,
      size = "s"
    ))
  })

  observeEvent(input$end_session_confirm, {
    removeModal()
    showNotification(
      "Ending PhenoMapR session\u2026",
      type = "message", duration = 2, id = "end_session_msg"
    )
    # Close the user's WebSocket connection cleanly first so the browser
    # gets a "Disconnected from server" overlay instead of a hard error,
    # then stop the R-side app event loop. stopApp() returns invisible()
    # to the caller of runApp() (typically the console).
    session$onSessionEnded(function() {
      shiny::stopApp(returnValue = invisible(NULL))
    })
    session$close()
  })

  # ------------------------------------------------------------------------
  # Reactive state container
  # ------------------------------------------------------------------------
  state <- reactiveValues(
    expression = NULL,         # Seurat / SCE / matrix / AnnData
    expr_summary = NULL,       # list(kind, n_genes, n_samples, ...)
    metadata = NULL,           # data.frame with cell IDs
    meta_columns = character(0),
    metadata_source = "(none)", # "object" | "upload" | "demo" | "(none)"
    demo_source = NULL,        # list describing CRA001160 subsample when demo loaded
    reference = NULL,          # "precog" / "tcga" / ... OR data.frame
    reference_label = "(none)",
    derive_bulk = NULL,
    derive_phen = NULL,
    derive_probe_annot = NULL,
    derive_meta_z_studies = list(),
    scores = NULL,
    groups = NULL,
    markers = NULL
  )

  # Captured plot / table objects for the panel download buttons.
  # Each render function below populates its slot at draw time
  # (panel_objects$<output_id> <- p or df); the matching
  # downloadHandler reads back from this store on click. See
  # phenomapr_register_plot_download() / *_table_download() in helpers.R.
  panel_objects <- reactiveValues()

  # On-screen plots use the default ggchicklet2 radius (3 pt). The
  # per-export radius lives in the plot-download modal and is applied
  # at preview / save time via .apply_chicklet_radius(), so no global
  # option observe() is needed any more.
  options(phenomapr.plot_radius_pt = 3)

  # Base font size for ggplot theme_minimal() across every panel.
  # Bumped from the previous mix of 12 / 13 to a single 14 for
  # easier readability at default plot sizes; downstream renderers
  # read this through `.theme_base_size()`.
  .theme_base_size <- function() 14L

  # -----------------------------------------------------------------------
  # Plot-download modal: per-panel observers + shared downloadHandler.
  #
  # Every plot panel's card header carries an actionButton whose
  # inputId is paste0(<panel_id>, "_modal_btn"). When clicked, we
  # stash the panel id in `active_plot_dl()` and pop up the
  # `phenomapr_plot_download_modal()` dialog. The dialog contains a
  # real downloadButton wired to `output$plot_dl_action`, which
  # reads the stashed id and the modal inputs (width / height / DPI /
  # format / base font size) to produce the file.
  #
  # Tables keep the direct downloadButton flow registered above via
  # `phenomapr_register_table_download()`.
  # -----------------------------------------------------------------------

  # Display labels + sensible default export dimensions per plot
  # panel. Defaults are seeded into the modal so users get a
  # reasonable starting point per panel.
  plot_panels <- list(
    celltype_count_plot      = list(label = "Cells per cell type",
                                     defaults = list(width = 8,  height = 5)),
    source_count_plot        = list(label = "Cells per source / group",
                                     defaults = list(width = 8,  height = 5)),
    celltype_source_plot     = list(label = "Cell type x source composition",
                                     defaults = list(width = 9,  height = 5)),
    reference_signature_plot = list(label = "Phenotype signature",
                                     defaults = list(width = 9,  height = 5)),
    score_dist_plot          = list(label = "Score distribution",
                                     defaults = list(width = 8,  height = 5)),
    score_rank_plot          = list(label = "Cells ordered by score",
                                     defaults = list(width = 9,  height = 5)),
    score_box_source_plot    = list(label = "Score by cell type and source",
                                     defaults = list(width = 10, height = 6)),
    umap_plot                = list(label = "Embedding",
                                     defaults = list(width = 10, height = 8)),
    group_by_celltype_plot   = list(label = "Per-cell-type group enrichment",
                                     defaults = list(width = 9,  height = 5)),
    marker_heatmap           = list(label = "Marker-gene heatmap",
                                     defaults = list(width = 12, height = 7))
  )

  active_plot_dl <- reactiveVal(NULL)

  # Persistent per-session "appearance" preferences (base font size +
  # corner radius). showModal() rebuilds the modal UI every time, so
  # without these we would snap back to the factory defaults on every
  # open. Width / height / DPI / format are NOT persisted because they
  # are panel-specific (e.g. UMAP defaults to 10x8, count plots default
  # to 8x5) -- each panel's own defaults win on first open.
  user_appearance_prefs <- reactiveValues(
    base_size = NULL,
    radius_pt = NULL
  )

  # Spin up one observer per plot panel. Each one watches its own
  # *_modal_btn input and opens the shared modal with that panel's
  # label + default dimensions. `local()` captures the loop variable
  # so each observer closes over the right `pid`.
  for (pid in names(plot_panels)) {
    local({
      panel_id <- pid
      meta     <- plot_panels[[panel_id]]
      observeEvent(input[[paste0(panel_id, "_modal_btn")]], {
        active_plot_dl(panel_id)
        # Seed the modal with this panel's per-panel defaults, then
        # overlay the user's persisted appearance tweaks (if any) so
        # they carry from one modal open to the next.
        defs <- meta$defaults
        prev_base   <- isolate(user_appearance_prefs$base_size)
        prev_radius <- isolate(user_appearance_prefs$radius_pt)
        if (!is.null(prev_base))   defs$base_size  <- prev_base
        if (!is.null(prev_radius)) defs$radius_pt  <- prev_radius
        showModal(phenomapr_plot_download_modal(
          panel_label = meta$label,
          defaults    = defs
        ))
      }, ignoreInit = TRUE)
    })
  }

  # Capture the user's appearance tweaks as they change so they
  # persist across modal opens (and across panels). We only track the
  # appearance knobs here -- width / height / dpi / format are
  # intentionally not persisted (see comment above).
  observe({
    v <- input$plot_dl_base_size
    if (!is.null(v) && is.finite(v)) user_appearance_prefs$base_size <- v
  })
  observe({
    v <- input$plot_dl_radius_pt
    if (!is.null(v) && is.finite(v)) user_appearance_prefs$radius_pt <- v
  })

  # ----------------------------------------------------------------------
  # Live preview inside the download modal. Mirrors what the saved
  # file will look like at the user's chosen width / height / base
  # font size. The renderPlot is dynamically sized via its `width` /
  # `height` arguments (functions of `input$plot_dl_width` /
  # `input$plot_dl_height`) so the on-screen aspect ratio matches
  # what the file will produce. The base font size is applied as a
  # theme override (same code path used by the export dispatcher) so
  # the preview reflects the actual font sizing of the exported file.
  #
  # DPI and file format intentionally do NOT affect the preview --
  # the screen is fixed-DPI and PDF/SVG vs PNG/JPEG/TIFF only changes
  # the on-disk encoding, not the rendering. That caveat is called
  # out in the modal's helpText so users do not get confused.
  #
  # marker_heatmap is special-cased: it is not a ggplot, the export
  # path re-runs plot_phenotype_markers() into a grDevices device,
  # and trying to render it twice (once for the main panel via
  # imageOutput, once for the preview via renderPlot) is both slow
  # and brittle. Instead the preview shows a friendly placeholder.

  # ---- Unit conversion helpers ----
  # The modal exposes `plot_dl_units = "in" | "cm"`. width / height
  # inputs carry the value in whatever unit the user picked; every
  # downstream consumer (renderPlot preview sizing, ggsave / device
  # calls in the downloadHandler) needs the value in INCHES because
  # that is what grDevices and ggsave expect.
  .CM_PER_INCH <- 2.54
  .to_inches <- function(v, units) {
    if (identical(units, "cm")) v / .CM_PER_INCH else v
  }
  # State: remember the previous unit so we know which direction to
  # rescale when the user toggles the radio. Default "in".
  prev_dl_units <- reactiveVal("in")

  # When the user toggles the units radio, rescale width / height to
  # KEEP THE PHYSICAL SIZE CONSTANT (8 in -> 20.32 cm rather than
  # 8 cm), update the labels to carry the new unit, and bump max +
  # step accordingly so the input keeps reasonable bounds.
  observeEvent(input$plot_dl_units, {
    cur <- input$plot_dl_units
    prv <- isolate(prev_dl_units())
    if (identical(cur, prv)) return()

    factor <- if (identical(prv, "in") && identical(cur, "cm")) {
      .CM_PER_INCH
    } else if (identical(prv, "cm") && identical(cur, "in")) {
      1 / .CM_PER_INCH
    } else {
      1
    }
    w <- suppressWarnings(as.numeric(isolate(input$plot_dl_width)))
    h <- suppressWarnings(as.numeric(isolate(input$plot_dl_height)))
    if (!is.finite(w) || w <= 0) w <- 8 * (if (identical(cur, "cm")) .CM_PER_INCH else 1)
    if (!is.finite(h) || h <= 0) h <- 6 * (if (identical(cur, "cm")) .CM_PER_INCH else 1)
    if (factor != 1) {
      w <- round(w * factor, 2)
      h <- round(h * factor, 2)
    }
    is_cm <- identical(cur, "cm")
    new_max  <- if (is_cm) 80 else 30
    new_step <- if (is_cm) 1  else 0.5
    new_w_label <- sprintf("Width (%s)",  if (is_cm) "cm" else "inches")
    new_h_label <- sprintf("Height (%s)", if (is_cm) "cm" else "inches")
    updateNumericInput(session, "plot_dl_width",
                       label = new_w_label,
                       value = w, max = new_max, step = new_step)
    updateNumericInput(session, "plot_dl_height",
                       label = new_h_label,
                       value = h, max = new_max, step = new_step)
    prev_dl_units(cur)
  }, ignoreInit = TRUE)

  # Reset the "previous units" memory each time a new modal opens so
  # the first rescale event always references the modal's seeded unit
  # rather than whatever the last modal left behind.
  observeEvent(active_plot_dl(), {
    prev_dl_units(input$plot_dl_units %||% "in")
  }, ignoreInit = FALSE)

  .plot_dl_preview_dims <- function(input,
                                    max_w = 640, max_h = 320,
                                    fallback = list(w = 8, h = 6)) {
    w <- suppressWarnings(as.numeric(input$plot_dl_width))
    h <- suppressWarnings(as.numeric(input$plot_dl_height))
    units <- input$plot_dl_units %||% "in"
    if (!is.finite(w) || w <= 0) w <- fallback$w * (if (identical(units, "cm")) .CM_PER_INCH else 1)
    if (!is.finite(h) || h <= 0) h <- fallback$h * (if (identical(units, "cm")) .CM_PER_INCH else 1)
    # Convert to inches before sizing the preview pixel box so the
    # screen aspect ratio matches what ggsave will produce on disk.
    w_in <- .to_inches(w, units)
    h_in <- .to_inches(h, units)
    # Fit (w_in x h_in) inches into a max_w x max_h pixel box
    # preserving aspect ratio. Whichever dimension hits the cap first
    # decides the scale factor.
    sw <- max_w / w_in
    sh <- max_h / h_in
    s  <- min(sw, sh)
    list(width  = max(120L, round(w_in * s)),
         height = max(120L, round(h_in * s)))
  }

  # Friendly empty / error state for the preview pane. Returns a
  # ggplot so renderPlot's auto-print works.
  .plot_dl_placeholder <- function(message) {
    ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0, y = 0,
        label = paste(strwrap(message, width = 42), collapse = "\n"),
        size = 4.2, colour = "#6c757d"
      ) +
      ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
      ggplot2::theme_void()
  }

  output$plot_dl_preview <- renderPlot({
    pid <- active_plot_dl()
    if (is.null(pid)) {
      return(.plot_dl_placeholder("Open a plot's download icon to see a preview."))
    }
    bsize  <- suppressWarnings(as.numeric(input$plot_dl_base_size))
    radius <- suppressWarnings(as.numeric(input$plot_dl_radius_pt))
    if (!is.finite(bsize)  || bsize  <= 0) bsize  <- 14
    if (!is.finite(radius) || radius <  0) radius <- 3

    if (identical(pid, "marker_heatmap")) {
      return(.plot_dl_placeholder(
        paste("Live preview is disabled for the marker heatmap because",
              "it requires a full ComplexHeatmap re-render. The saved",
              "file is generated by re-running plot_phenotype_markers()",
              "at the chosen size, DPI and format.")
      ))
    }

    p <- panel_objects[[pid]]
    if (is.null(p)) {
      return(.plot_dl_placeholder(
        paste0("Plot '", pid, "' has not been rendered yet. ",
               "Close this dialog, scroll to the panel so it renders, ",
               "and re-open the download.")
      ))
    }

    if (inherits(p, c("ggplot", "patchwork"))) {
      p <- p + ggplot2::theme(text = ggplot2::element_text(size = bsize))
      p <- .apply_chicklet_radius(p, radius)
      return(p)
    }
    # Non-ggplot fallback (e.g. ComplexHeatmap, recordedplot). We have
    # to print() these manually because renderPlot's auto-print is
    # tuned for ggplot / lattice.
    print(p)
    invisible(NULL)
  },
  width  = function() .plot_dl_preview_dims(input)$width,
  height = function() .plot_dl_preview_dims(input)$height,
  res    = 96
  )

  # Shared downloadHandler: produces the file for whichever plot is
  # currently active. `marker_heatmap` is special-cased because its
  # underlying object is not a ggplot -- it has to be re-rendered
  # from the cached args via plot_phenotype_markers(). All other
  # panels stash a ggplot/patchwork object in `panel_objects[[pid]]`
  # at render time, so they share the .phenomapr_save_plot() path.
  output$plot_dl_action <- downloadHandler(
    filename = function() {
      pid <- isolate(active_plot_dl()) %||% "plot"
      fmt <- isolate(input$plot_dl_format) %||% "png"
      phenomapr_dl_filename(pid, fmt)
    },
    content = function(file) {
      pid <- isolate(active_plot_dl())
      if (is.null(pid)) {
        .phenomapr_write_placeholder_png(
          file, "No active plot for download.")
        removeModal()
        return(invisible(NULL))
      }
      fmt    <- isolate(input$plot_dl_format)   %||% "png"
      units  <- isolate(input$plot_dl_units)    %||% "in"
      w      <- as.numeric(isolate(input$plot_dl_width)  %||% 8)
      h      <- as.numeric(isolate(input$plot_dl_height) %||% 6)
      dpi    <- as.numeric(isolate(input$plot_dl_dpi)    %||% 300)
      bsize  <- as.numeric(isolate(input$plot_dl_base_size) %||% 14)
      radius <- as.numeric(isolate(input$plot_dl_radius_pt) %||% 3)

      # Validate numeric inputs; fall back to defaults if NA / <= 0.
      if (!is.finite(w)      || w      <= 0) w      <- 8
      if (!is.finite(h)      || h      <= 0) h      <- 6
      if (!is.finite(dpi)    || dpi    <= 0) dpi    <- 300
      if (!is.finite(radius) || radius < 0)  radius <- 3

      # ggsave + grDevices::png/jpeg/tiff/pdf/svg all expect width and
      # height in inches. Convert here so the rest of the export path
      # is unit-agnostic. .to_inches is a no-op for "in".
      w <- .to_inches(w, units)
      h <- .to_inches(h, units)

      if (identical(pid, "marker_heatmap")) {
        # Heatmap path: re-render via plot_phenotype_markers() into a
        # device matching the user's format selection. Width/height
        # are interpreted as inches; for raster formats we multiply by
        # dpi to get pixel dimensions.
        args <- isolate(marker_heatmap_args())
        if (is.null(args)) {
          .phenomapr_write_placeholder_png(
            file, "Draw the heatmap first, then re-try the download."
          )
          removeModal()
          return(invisible(NULL))
        }
        # Mirror the live block-border toggle so the export matches
        # what the user is currently looking at.
        borders_on <- isTRUE(isolate(input$hm_block_borders))
        args$outline_marker_blocks <- borders_on
        args$block_outline_color   <- if (borders_on) "black" else "white"
        # Use a slightly slimmer line (~1.5pt) so the borders read as
        # clean dividers rather than thick chrome around each block.
        args$block_outline_lwd     <- if (borders_on) 1.5 else 1
        args$color_mark_labels_by_celltype <- isTRUE(isolate(input$hm_color_mark_labels))
        args$color_schemes <- marker_heatmap_color_schemes_from_input(isolate(input))
        fmt_lc <- tolower(fmt)
        ok <- tryCatch({
          switch(fmt_lc,
            png  = grDevices::png(file, width = w * dpi, height = h * dpi, res = dpi),
            jpeg = grDevices::jpeg(file, width = w * dpi, height = h * dpi, res = dpi, quality = 95),
            tiff = grDevices::tiff(file, width = w * dpi, height = h * dpi, res = dpi),
            pdf  = grDevices::pdf(file, width = w, height = h),
            svg  = grDevices::svg(file, width = w, height = h),
            stop("Unsupported format: ", fmt_lc)
          )
          on.exit(if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(),
                                                     silent = TRUE),
                  add = TRUE)
          do.call(PhenoMapR::plot_phenotype_markers,
                  c(args, list(draw = TRUE)))
          TRUE
        }, error = function(e) {
          if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
          FALSE
        })
        if (!isTRUE(ok)) {
          .phenomapr_write_placeholder_png(
            file, "plot_phenotype_markers() failed for this download."
          )
        }
        removeModal()
        return(invisible(NULL))
      }

      plot_obj <- isolate(panel_objects[[pid]])
      if (is.null(plot_obj)) {
        .phenomapr_write_placeholder_png(
          file, paste0("Plot '", pid, "' has not been rendered yet.")
        )
        removeModal()
        return(invisible(NULL))
      }

      ok <- tryCatch({
        .phenomapr_save_plot(
          file      = file,
          plot_obj  = plot_obj,
          format    = fmt,
          width     = w,
          height    = h,
          dpi       = dpi,
          base_size = bsize,
          radius_pt = radius
        )
        TRUE
      }, error = function(e) {
        message("plot_dl_action failed for ", pid, ": ",
                conditionMessage(e))
        FALSE
      })
      if (!isTRUE(ok)) {
        .phenomapr_write_placeholder_png(
          file, paste0("Could not export plot '", pid, "'.")
        )
      }
      removeModal()
      invisible(NULL)
    }
  )

  # ------------------------------------------------------------------------
  # Hybrid file pickers (browser upload + server filesystem browse).
  # Each `*_pick` reactive returns NULL or list(datapath, name, source).
  # See phenomapr_file_input() / phenomapr_file_pick() in helpers.R.
  # ------------------------------------------------------------------------
  shiny_file_roots <- phenomapr_app_server_roots()
  expr_file_pick       <- phenomapr_file_pick("expr_file",       input, output, session, shiny_file_roots,
                                              accept = c(".rds", ".h5", ".h5ad", ".tsv", ".csv", ".txt"))
  meta_file_pick       <- phenomapr_file_pick("meta_file",       input, output, session, shiny_file_roots,
                                              accept = c(".rds", ".tsv", ".csv", ".txt"))
  custom_ref_file_pick <- phenomapr_file_pick("custom_ref_file", input, output, session, shiny_file_roots,
                                              accept = c(".rds", ".tsv", ".csv", ".txt"))
  derive_bulk_file_pick <- phenomapr_file_pick("derive_bulk_file", input, output, session, shiny_file_roots,
                                               accept = c(".rds", ".tsv", ".csv", ".txt"))
  derive_phen_file_pick <- phenomapr_file_pick("derive_phen_file", input, output, session, shiny_file_roots,
                                               accept = c(".rds", ".tsv", ".csv", ".txt"))
  derive_probe_file_pick <- phenomapr_file_pick("derive_probe_file", input, output, session, shiny_file_roots,
                                                accept = c(".rds", ".tsv", ".csv", ".txt"))
  umap_upload_pick     <- phenomapr_file_pick("umap_upload",     input, output, session, shiny_file_roots,
                                              accept = c(".tsv", ".csv", ".txt", ".rds"))

  # ------------------------------------------------------------------------
  # 1. Expression upload
  # ------------------------------------------------------------------------
  observeEvent(expr_file_pick(), {
    pick <- expr_file_pick()
    # NULL transition: the user clicked the "remove" affordance on the
    # file picker. Treat this as a full data reset so the rest of the
    # app reflects the cleared state. The metadata reset is gated on
    # the source: only data that came WITH the expression object (or
    # was synthesized from it) is dropped -- a separately uploaded
    # metadata file (`source = "upload"`) is preserved so it does not
    # silently disappear with the expression.
    if (is.null(pick)) {
      state$expression <- NULL
      state$expr_summary <- NULL
      state$demo_source <- NULL
      if (!identical(state$metadata_source %||% "", "upload")) {
        state$metadata <- NULL
        state$meta_columns <- character(0)
        state$metadata_source <- "(none)"
      }
      return()
    }
    phenomapr_busy_show("Reading expression file...", pick$name)
    # CRITICAL: defer the heavy synchronous file read into a
    # later::later() callback. Without this, the observer body would
    # block libuv\'s I/O loop for the duration of the read, so the
    # phenomapr-busy-show custom message we just queued never gets
    # transmitted to the browser until AFTER the read finishes -- by
    # which point we\'ve already queued the busy-hide right behind it
    # and the user sees no popup at all. Yielding via later() lets
    # this observer\'s flushReact complete first (which actually
    # delivers the show message), and the heavy read then runs in
    # a fresh tick where it\'s free to block libuv without losing
    # the show signal. See the busy-overlay file-comment in
    # helpers.R for the full rationale.
    pick_local <- pick
    sess <- session
    later::later(function() {
      res <- tryCatch(
        parse_expression_upload(pick_local$datapath, pick_local$name),
        error = function(e) {
          showNotification(paste0("Upload failed: ", conditionMessage(e)),
                           type = "error", duration = 8, session = sess)
          NULL
        }
      )
      # Hide the popup AS SOON AS the heavy file read is done, BEFORE
      # we touch state$. State assignments invalidate downstream
      # reactives; deferring them with another later() ensures the
      # busy-hide custom message lands in its own (tiny) flush ahead
      # of the heavy renderUI/renderPlot cascade triggered by the
      # state mutations.
      phenomapr_busy_hide(session = sess)
      if (is.null(res)) return()

      res_obj <- res
      later::later(function() {
        tryCatch({
          state$expression <- res_obj$object
          state$expr_summary <- res_obj
          state$demo_source <- NULL
          md <- extract_object_metadata(res_obj$object)
          state$metadata <- md
          state$meta_columns <- if (!is.null(md)) colnames(md) else character(0)
          state$metadata_source <- if (!is.null(md)) "object" else "(none)"
          shiny::showNotification(
            paste(res_obj$notes, collapse = " "),
            type = "message", duration = 5,
            session = sess
          )

        # When loading an AnnData (.h5ad), surface a clear message if the
        # in-object metadata could not be auto-extracted. Some round-trip
        # paths (e.g., anndataR -> Python anndata) can produce obs columns
        # with pandas extension dtypes that the standard
        # reticulate::py_to_r() pipeline can't decode in one pass. Our
        # extractor tries three layered fallbacks (decategorize, columnwise);
        # if all of them came back empty we tell the user explicitly so they
        # can either fix the source object or upload metadata separately.
        if (identical(res_obj$kind %||% "", "anndata")) {
          diag_msgs <- if (is.null(md)) {
            tryCatch(
              phenomapr_anndata_obs_df(res_obj$object)$warnings,
              error = function(e) character(0)
            )
          } else {
            attr(md, "anndata_obs_warnings", exact = TRUE) %||% character(0)
          }
          if (is.null(md)) {
            diag <- if (length(diag_msgs)) {
              paste(utils::head(diag_msgs, 3L), collapse = " | ")
            } else {
              "AnnData.obs returned an empty data.frame."
            }
            shiny::showNotification(
              paste0(
                "AnnData metadata could not be auto-detected. ",
                "You can still score the cells, but cell-type-aware steps ",
                "and groupings need a metadata table. Upload your cell ",
                "metadata in the next section, or re-export the .h5ad ",
                "with plain object/string columns. Diagnostic: ", diag
              ),
              type = "warning", duration = 12, id = "anndata_meta_warn",
              session = sess
            )
          } else if (length(diag_msgs)) {
            shiny::showNotification(
              sprintf(
                "AnnData metadata loaded with %d cell%s x %d column%s, but %d non-fatal warning%s during extraction (e.g. %s).",
                nrow(md), ifelse(nrow(md) == 1L, "", "s"),
                ncol(md) - 1L, ifelse(ncol(md) - 1L == 1L, "", "s"),
                length(diag_msgs), ifelse(length(diag_msgs) == 1L, "", "s"),
                substr(diag_msgs[[1L]], 1L, 140L)
              ),
              type = "message", duration = 10, id = "anndata_meta_partial",
              session = sess
            )
          } else if (!is.null(md)) {
            shiny::showNotification(
              sprintf("Detected %d cell%s x %d metadata column%s from AnnData.obs.",
                      nrow(md), ifelse(nrow(md) == 1L, "", "s"),
                      ncol(md) - 1L, ifelse(ncol(md) - 1L == 1L, "", "s")),
              type = "message", duration = 5, id = "anndata_meta_ok",
              session = sess
            )
          }
        }
      }, error = function(e) {
        shiny::showNotification(
          paste0("Internal error while applying upload: ", conditionMessage(e)),
          type = "error", duration = 10, session = sess
        )
      })
      }, delay = 0)  # close inner later::later() (state propagation)
    }, delay = 0)    # close outer later::later() (heavy file read)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$use_demo, {
    demo <- make_shiny_demo_dataset()
    m <- demo$expression
    state$expression <- m
    state$expr_summary <- summarize_expression_object(m)
    state$expr_summary$notes <- sprintf(
      "Demo subsample (%d genes \u00d7 %d cells) from %s with UMAP coordinates.",
      nrow(m), ncol(m), demo$source_info$accession %||% "CRA001160"
    )
    state$expr_summary$is_demo <- TRUE
    state$metadata <- demo$metadata
    state$meta_columns <- colnames(state$metadata)
    state$metadata_source <- "demo"
    state$demo_source <- demo$source_info
    updateSliderInput(session, "percentile", value = 0.10)
    updateRadioButtons(session, "reference_choice", selected = "precog")
    demo_cts <- get_cancer_types("precog")
    demo_panc <- if ("Pancreatic" %in% demo_cts) "Pancreatic" else demo_cts[1L]
    updateSelectInput(session, "cancer_type", choices = demo_cts, selected = demo_panc)
    n_fmt <- format(demo$source_info$n_cells_sampled %||% ncol(m), big.mark = ",")
    showNotification(
      sprintf(
        "Demo loaded: %s pre-selected cells from %s (PRECOG Pancreatic selected).",
        n_fmt, demo$source_info$accession %||% "CRA001160"
      ),
      type = "message", duration = 6
    )
  })

  output$demo_source_panel <- renderUI({
    if (!identical(state$metadata_source %||% "", "demo")) return(NULL)
    info <- state$demo_source
    if (is.null(info)) return(NULL)

    fmt_int <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
    tags$details(
      class = "phenomapr-about-markers-details phenomapr-demo-source-details",
      tags$summary(
        icon("database"),
        tags$span(
          class = "phenomapr-about-markers-summary-label",
          sprintf(" Demo source: %s", info$accession %||% "CRA001160")
        ),
        tags$span(
          class = "phenomapr-about-markers-summary-hint",
          " \u2014 click to expand"
        )
      ),
      tags$div(
        class = "phenomapr-demo-source-body",
        tags$p(
          tags$strong("Dataset: "),
          sprintf("%s \u2014 %s", info$accession, info$title)
        ),
        tags$p(
          tags$strong("Cancer type: "),
          info$cancer_type
        ),
        tags$p(
          tags$strong("Cohort: "),
          sprintf(
            "%s cells from %s patients (%s tumors, %s normal pancreases)",
            fmt_int(info$n_cells_total),
            fmt_int(info$n_patients),
            fmt_int(info$n_tumors),
            fmt_int(info$n_controls)
          )
        ),
        tags$p(
          tags$strong("Publication: "),
          tags$a(
            href = info$paper_url,
            target = "_blank",
            rel = "noopener noreferrer",
            info$paper_label
          )
        ),
        tags$p(
          tags$strong("Processed data: "),
          tags$a(
            href = info$data_source_url,
            target = "_blank",
            rel = "noopener noreferrer",
            info$data_source
          )
        ),
        tags$p(
          tags$strong("This load: "),
          sprintf(
            "Pre-selected subsample of %s cells \u00d7 %s genes (fixed demo subset). Default signature: %s.",
            fmt_int(info$n_cells_sampled),
            fmt_int(info$n_genes_sampled),
            info$reference_signature
          )
        )
      )
    )
  })

  # Optional metadata upload
  observeEvent(meta_file_pick(), {
    pick <- meta_file_pick()
    # Cleared by the user via the remove affordance. Drop any
    # uploaded metadata, falling back to whatever the expression
    # object carried (if anything). We deliberately do NOT touch
    # state$metadata when the source was "object" -- removing the
    # metadata upload should not also strip metadata that came with
    # the loaded Seurat/SCE/AnnData.
    if (is.null(pick)) {
      if (identical(state$metadata_source %||% "", "upload")) {
        state$metadata <- NULL
        state$meta_columns <- character(0)
        state$metadata_source <- "(none)"
        # If the expression object still has in-object metadata, restore it.
        if (!is.null(state$expression)) {
          md_obj <- tryCatch(extract_object_metadata(state$expression),
                             error = function(e) NULL)
          if (!is.null(md_obj)) {
            state$metadata <- md_obj
            state$meta_columns <- colnames(md_obj)
            state$metadata_source <- "object"
          }
        }
      }
      return()
    }
    phenomapr_busy_show("Loading cell metadata...", pick$name)
    # See the expression-upload observer above for the deferred-work
    # rationale: heavy reads MUST run inside later::later() so the
    # busy_show custom message can be transmitted to the browser
    # before libuv gets blocked by the read.
    pick_local <- pick
    sess <- session
    later::later(function() {
      md <- tryCatch(
        parse_metadata_upload(pick_local$datapath, pick_local$name),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 8,
                           session = sess)
          NULL
        }
      )
      phenomapr_busy_hide(session = sess)
      if (is.null(md)) return()

      md_local <- md
      later::later(function() {
        tryCatch({
          state$metadata <- md_local
          state$meta_columns <- colnames(md_local)
          state$metadata_source <- "upload"
          shiny::showNotification("Metadata loaded.",
                                  type = "message", duration = 4,
                                  session = sess)
        }, error = function(e) {
          shiny::showNotification(
            paste0("Internal error while applying metadata: ",
                   conditionMessage(e)),
            type = "error", duration = 10, session = sess
          )
        })
      }, delay = 0)
    }, delay = 0)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  output$meta_columns_available <- reactive({ length(state$meta_columns) > 0 })
  outputOptions(output, "meta_columns_available", suspendWhenHidden = FALSE)

  # Metadata-upload <summary> body + open-state toggle.
  #
  # The <details> WRAPPER is rendered statically in the sidebar UI
  # together with the file picker, so the picker's DOM identity stays
  # stable across state$metadata transitions (see the sidebar UI
  # comment for the rationale). Here we only:
  #
  #   1. server-render the SUMMARY content, which morphs between a
  #      "No metadata detected -- upload a tabular metadata file"
  #      prompt (when expression is loaded but no metadata was found)
  #      and a more neutral "Override / supplement with a tabular
  #      metadata file" label (when metadata is already in hand).
  #
  #   2. toggle the `open` attribute on the static <details> via a
  #      custom JS message ("phenomapr-set-details-open"), so the
  #      panel auto-opens when no metadata is present and auto-closes
  #      once metadata arrives. The user can always click the chevron
  #      to override.
  # Children of the static <summary>. Only these inner nodes change
  # across state$metadata transitions; the <summary> wrapper itself
  # (created by the uiOutput container in the sidebar UI) keeps a
  # stable DOM identity, which is what guarantees the file picker
  # below it never gets torn down.
  output$metadata_upload_summary <- renderUI({
    has_expr  <- !is.null(state$expression)
    no_meta   <- is.null(state$metadata) || !length(state$meta_columns)
    prompt    <- isTRUE(has_expr && no_meta)

    if (prompt) {
      tagList(
        tags$span(class = "mu-badge", icon("upload")),
        tags$strong(
          "No metadata detected \u2014 upload a tabular metadata file"
        )
      )
    } else {
      "Override / supplement with a tabular metadata file"
    }
  })

  # Auto-open / auto-close the static <details> AND tag it with the
  # `mu-prompt` class so CSS can apply the salmon-tinted prompt look
  # to the summary. We push the desired state via a custom message;
  # the JS handler in helpers.R sets or removes both the `open`
  # attribute and the class.
  observe({
    has_expr <- !is.null(state$expression)
    no_meta  <- is.null(state$metadata) || !length(state$meta_columns)
    prompt   <- isTRUE(has_expr && no_meta)
    session$sendCustomMessage(
      "phenomapr-set-details-state",
      list(id     = "metadata_upload_details",
           open   = prompt,
           prompt = prompt)
    )
  })

  # Status banner above the metadata dropdowns: explains where the metadata
  # currently in use came from (object / upload / demo). When no metadata
  # has been loaded yet we render NOTHING here -- the collapsible
  # "No metadata detected -- upload a tabular metadata file" header on the
  # details/summary block immediately below already conveys the empty
  # state, so the previous "No metadata loaded yet ..." em-paragraph was
  # redundant.
  output$metadata_status <- renderUI({
    if (is.null(state$metadata) || !length(state$meta_columns)) {
      return(NULL)
    }
    src <- state$metadata_source %||% "(unknown)"
    kind <- state$expr_summary$kind %||% "loaded"
    label <- switch(
      src,
      object = sprintf("Auto-detected from %s object", kind),
      upload = "Loaded from uploaded file",
      demo   = "Demo dataset metadata",
      sprintf("Loaded (%s)", src)
    )
    tags$div(
      class = "metadata-status metadata-status-ok",
      tags$div(class = "ms-title",
               tags$span(class = "ms-badge", icon("check-circle")),
               tags$strong(label)),
      tags$div(class = "ms-detail",
               sprintf("%s cells × %d columns",
                       .fmt_int(nrow(state$metadata)),
                       ncol(state$metadata)))
    )
  })

  # Pick a sensible default column from `cols` using a list of regex patterns
  # (case-insensitive). The first column matching the first pattern wins;
  # if none match, the next pattern is tried; falls back to "(none)".
  .pick_default <- function(cols, patterns, exclude = NULL) {
    available <- setdiff(cols, c("(none)", exclude))
    for (p in patterns) {
      hits <- grep(p, available, ignore.case = TRUE, value = TRUE, perl = TRUE)
      if (length(hits)) return(hits[1L])
    }
    "(none)"
  }

  observe({
    cols <- c("(none)", state$meta_columns)

    # Cell ID: prefer the `.cell_id` column we always inject in
    # extract_object_metadata(); fall back to common barcode names.
    cell_default <- if (".cell_id" %in% cols) ".cell_id" else {
      .pick_default(cols, c("^cell_?id$", "^cell$", "^barcode$", "^CellID$"))
    }
    if (cell_default == "(none)" && length(state$meta_columns)) {
      cell_default <- state$meta_columns[1L]
    }

    # Cell type: strong matches before generic "type"; exclude `orig.ident`
    # and `feature_type` which incidentally contain the substring.
    ct_default <- .pick_default(
      cols,
      patterns = c(
        "^cell_type$",
        "cell_type_minor|minor[._\\s-]?lineage",
        "^(cell[._\\s]?type|celltype)$",
        "annotated_celltype|predicted\\.celltype|leiden_celltype|seurat_clusters_annot",
        "^Annotation$|^annotation$",
        "^type$"
      ),
      exclude = c(
        "orig.ident", "feature_type",
        "cell_type_major", "cell_type_original"
      )
    )

    # Source / group: cohort / sample / condition / donor / patient / tissue.
    src_default <- .pick_default(
      cols,
      patterns = c(
        "^source$|^Source$|^sample_source$",
        "^group$|^Group$|^cohort$|^Cohort$",
        "^condition$|^Condition$|^treatment$",
        "^donor$|^patient$|^subject$",
        "^tissue$|^region$|^organ$",
        "^sample$|^Sample$|^orig\\.ident$"
      )
    )

    updateSelectInput(session, "meta_cell_id_col",   choices = cols, selected = cell_default)
    updateSelectInput(session, "meta_cell_type_col", choices = cols, selected = ct_default)
    updateSelectInput(session, "meta_source_col",    choices = cols, selected = src_default)

    # Pseudobulk grouping: prefer patient / donor / sample / orig.ident /
    # cluster / region / slide / core columns over the first-available
    # metadata column. Falls back to the first column otherwise so the
    # input still has a meaningful selection.
    pb_default <- .pick_default(
      cols,
      patterns = c(
        "^(patient|donor|subject)(_id)?$",
        "^(sample|sample_id|library|orig\\.ident)$",
        "^(seurat_clusters|cluster|leiden|louvain)(_id)?$",
        "^(slide|core|spot|fov|section|capture_area)(_id)?$",
        "^(tissue|region|condition|treatment)$"
      )
    )
    if (pb_default == "(none)" && length(state$meta_columns)) {
      pb_default <- state$meta_columns[1L]
    }
    updateSelectInput(session, "pseudobulk_group_by",
                      choices = setdiff(cols, "(none)"),
                      selected = pb_default)
  })

  # ------------------------------------------------------------------------
  # 3. Score sidebar — sync with what was detected in 1. Data
  # ------------------------------------------------------------------------
  # `score_data_status` is the top-of-sidebar "Scoring summary" block.
  # It surfaces the three pieces of context users need to confirm before
  # clicking "Compute PhenoMapR scores":
  #
  #   1. Input data type detected in 1. Data (Seurat / SCE / AnnData /
  #      Spatial / Matrix), plus the assay name and available layers.
  #   2. Phenotype signature source chosen in 2. Phenotype (PRECOG /
  #      TCGA / Pediatric PRECOG / ICI PRECOG / custom).
  #   3. Cancer / tissue type the signature is being interpreted
  #      against (only applies to built-in signatures).
  #
  # Each row uses a small teal icon + label/value pair so the user can
  # scan the (data, signature, cancer-type) triple at a glance.
  output$score_data_status <- renderUI({
    s <- state$expr_summary

    # ---- Row 1: input data type ------------------------------------------
    if (is.null(s)) {
      data_row <- tags$div(
        class = "sds-row sds-row-empty",
        tags$div(class = "sds-icon", icon("circle-info")),
        tags$div(class = "sds-content",
                 tags$div(class = "sds-label", "Input data"),
                 tags$em("Load an expression dataset in 1. Data first."))
      )
    } else {
      kind_label <- switch(
        s$kind %||% "loaded",
        seurat  = "Seurat object",
        sce     = "SingleCellExperiment",
        spatial = "Spatial dataset",
        matrix  = "Expression matrix",
        anndata = "AnnData",
        sprintf("%s object", s$kind %||% "loaded")
      )
      bits <- character(0)
      if (!is.na(s$default_assay %||% NA_character_) &&
          nzchar(s$default_assay)) {
        bits <- c(bits, sprintf("assay: %s", s$default_assay))
      }
      if (length(s$layers_avail %||% character(0))) {
        bits <- c(bits, sprintf("layers: %s",
                                paste(s$layers_avail, collapse = ", ")))
      }
      # For plain matrix / data.frame inputs there is no "assay" or
      # "layer" concept to report, so we instead surface the heuristic
      # single-cell vs bulk classification from
      # detect_expression_format() (cached in expr_diagnostics()). This
      # gives users a one-glance confirmation that the file they
      # uploaded is being interpreted the way they expect before
      # clicking "Compute PhenoMapR scores".
      if (identical(s$kind %||% "", "matrix")) {
        diag <- tryCatch(expr_diagnostics(), error = function(e) NULL)
        if (!is.null(diag) && length(diag$sc_or_bulk %||% "")) {
          sb_label <- switch(
            diag$sc_or_bulk,
            single_cell = "predicted single-cell",
            bulk        = "predicted bulk",
            unclear     = "single-cell vs bulk unclear",
            diag$sc_or_bulk
          )
          if (nzchar(sb_label)) bits <- c(bits, sb_label)
        }
      }
      data_row <- tags$div(
        class = "sds-row sds-row-ok",
        tags$div(class = "sds-icon", icon("check-circle")),
        tags$div(class = "sds-content",
                 tags$div(class = "sds-label", "Input data"),
                 tags$div(class = "sds-value", kind_label),
                 if (length(bits))
                   tags$div(class = "sds-detail",
                            paste(bits, collapse = "  \u00B7  ")))
      )
    }

    # ---- Row 2: phenotype signature source --------------------------------
    rc <- input$reference_choice %||% "(none)"
    phen_label <- switch(
      rc,
      "precog"           = "PRECOG (built-in)",
      "tcga"             = "TCGA (built-in)",
      "pediatric_precog" = "Pediatric PRECOG (built-in)",
      "pediatric"        = "Pediatric PRECOG (built-in)",
      "ici_precog"       = "ICI PRECOG (built-in)",
      "ici"              = "ICI PRECOG (built-in)",
      "_custom"          = {
        cs <- input$custom_source %||% "upload"
        if (identical(cs, "derive"))
          "Custom signature (derived from your bulk + phenotype)"
        else
          "Custom signature (uploaded)"
      },
      sprintf("%s", rc)
    )
    phen_icon <- if (identical(rc, "_custom")) "flask" else "dna"
    phen_row <- tags$div(
      class = "sds-row sds-row-ok",
      tags$div(class = "sds-icon", icon(phen_icon)),
      tags$div(class = "sds-content",
               tags$div(class = "sds-label", "Phenotype source"),
               tags$div(class = "sds-value", phen_label))
    )

    # ---- Row 3: cancer / tissue type --------------------------------------
    if (identical(rc, "_custom")) {
      ct_row <- tags$div(
        class = "sds-row sds-row-na",
        tags$div(class = "sds-icon", icon("minus")),
        tags$div(class = "sds-content",
                 tags$div(class = "sds-label", "Cancer / tissue type"),
                 tags$em("Not applicable for custom signatures"))
      )
    } else {
      ct <- input$cancer_type
      if (is.null(ct) || !nzchar(ct)) {
        ct_row <- tags$div(
          class = "sds-row sds-row-empty",
          tags$div(class = "sds-icon", icon("circle-question")),
          tags$div(class = "sds-content",
                   tags$div(class = "sds-label", "Cancer / tissue type"),
                   tags$em("Select one in 2. Phenotype."))
        )
      } else {
        ct_row <- tags$div(
          class = "sds-row sds-row-ok",
          tags$div(class = "sds-icon", icon("ribbon")),
          tags$div(class = "sds-content",
                   tags$div(class = "sds-label", "Cancer / tissue type"),
                   tags$div(class = "sds-value", ct))
        )
      }
    }

    tags$div(class = "score-data-status", data_row, phen_row, ct_row)
  })

  # JS-readable flags that gate the slot/assay UI in the 3. Score
  # sidebar. Computed here (not inline in renderUI) so they update
  # the moment state$expr_summary changes and stay live even while
  # the user is on a different tab.
  #
  #   * score_show_slot_block = TRUE iff the loaded expression input
  #     has a meaningful concept of an assay / layer (Seurat / SCE /
  #     SpatialExperiment / AnnData). Plain matrices and data.frames
  #     return FALSE -- the slot radio + assay textInput are then
  #     hidden because they have no effect on PhenoMap()'s matrix
  #     path. NULL state$expr_summary (no upload yet) also returns
  #     FALSE so the empty-state row in score_data_status is the
  #     ONLY message users see before they load anything.
  #
  #   * score_have_matrix is the dual flag, TRUE iff something IS
  #     loaded AND it's matrix-class. Used by the small "matrices
  #     scored directly" note that sits in place of the hidden slot
  #     block, so the message only appears for matrix uploads -- not
  #     in the pre-upload empty state.
  output$score_show_slot_block <- reactive({
    s <- state$expr_summary
    !is.null(s) && !identical(s$kind %||% "", "matrix")
  })
  outputOptions(output, "score_show_slot_block", suspendWhenHidden = FALSE)

  output$score_have_matrix <- reactive({
    s <- state$expr_summary
    !is.null(s) && identical(s$kind %||% "", "matrix")
  })
  outputOptions(output, "score_have_matrix", suspendWhenHidden = FALSE)

  # Pseudobulk is only meaningful for single-cell / spatial inputs (including
  # matrix uploads classified as single-cell-like). Bulk matrices are already
  # sample-level profiles.
  score_allow_pseudobulk <- reactive({
    s <- state$expr_summary
    if (is.null(s)) return(FALSE)
    kind <- s$kind %||% ""
    if (kind %in% c("seurat", "sce", "spatial", "anndata")) return(TRUE)
    if (isTRUE(s$is_demo)) return(TRUE)
    if (!identical(kind, "matrix")) return(FALSE)
    d <- expr_diagnostics()
    if (is.null(d)) return(FALSE)
    identical(d$sc_or_bulk, "single_cell")
  })

  output$score_allow_pseudobulk <- score_allow_pseudobulk
  outputOptions(output, "score_allow_pseudobulk", suspendWhenHidden = FALSE)

  output$score_is_bulk_matrix <- reactive({
    s <- state$expr_summary
    if (is.null(s) || !identical(s$kind %||% "", "matrix")) return(FALSE)
    if (isTRUE(s$is_demo)) return(FALSE)
    d <- expr_diagnostics()
    !is.null(d) && identical(d$sc_or_bulk, "bulk")
  })
  outputOptions(output, "score_is_bulk_matrix", suspendWhenHidden = FALSE)

  observeEvent(state$expr_summary, {
    if (!isTRUE(score_allow_pseudobulk())) {
      updateCheckboxInput(session, "pseudobulk", value = FALSE)
    }
  }, ignoreNULL = FALSE)

  # Whenever a fresh object lands in state$expr_summary, refresh:
  #   * score_slot choices  -> only layers we actually have (drop scale.data
  #     entirely; if the object lacks a normalized "data" layer, default to
  #     "counts" so PhenoMap() doesn't error out).
  #   * score_assay value   -> the detected default assay name (Seurat/SCE).
  # For plain matrices we leave the defaults alone (data + blank assay) —
  # PhenoMap() handles matrices without consulting either field.
  observeEvent(state$expr_summary, {
    s <- state$expr_summary
    if (is.null(s)) return()

    nice <- c("data" = "data (log-normalized)",
              "counts" = "counts (raw)")

    if (identical(s$kind, "seurat") || identical(s$kind, "spatial")) {
      avail <- intersect(c("data", "counts"), s$layers_avail %||% character(0))
      if (!length(avail)) {
        # Couldn't introspect layers (very old Seurat object etc.) -- expose
        # both options and let the user choose. Default to "data".
        avail <- c("data", "counts")
      }
      choices <- setNames(avail, nice[avail])
      sel <- if ("data" %in% avail) "data" else "counts"
      updateRadioButtons(session, "score_slot",
                         choices = choices, selected = sel)
    } else if (identical(s$kind, "sce")) {
      # SCE: layers ARE the assay matrices. Prefer logcounts -> "data" map.
      la <- s$layers_avail %||% character(0)
      has_lognorm <- any(grepl("(?i)^(logcounts|normcounts|data|lognorm)$", la))
      avail <- c(if (has_lognorm) "data", if (any(la == "counts")) "counts")
      if (!length(avail)) avail <- c("data", "counts")
      choices <- setNames(avail, nice[avail])
      sel <- if ("data" %in% avail) "data" else "counts"
      updateRadioButtons(session, "score_slot",
                         choices = choices, selected = sel)
    } else {
      # AnnData / plain matrix — slot radio is informational only for the
      # matrix case; show both choices so the AnnData path can pick.
      avail <- c("data", "counts")
      choices <- setNames(avail, nice[avail])
      updateRadioButtons(session, "score_slot",
                         choices = choices, selected = "data")
    }

    # Carry over the detected assay name into the text input so users don't
    # have to retype it (still editable; blank still means "auto").
    da <- s$default_assay %||% NA_character_
    if (!is.na(da) && nzchar(da)) {
      updateTextInput(session, "score_assay", value = da)
    }
  }, ignoreInit = TRUE)

  # Per-cell unified table (score + group + cell_type + source) shared by
  # multiple downstream tabs (UMAP, score-by-source, marker heatmap, …).
  cell_table <- reactive({
    if (is.null(state$scores)) return(NULL)
    build_cell_table(
      scores = state$scores,
      groups = state$groups,
      metadata = state$metadata,
      cell_id_col = input$meta_cell_id_col,
      cell_type_col = input$meta_cell_type_col,
      source_col = input$meta_source_col,
      score_column = active_score_column(state$scores)
    )
  })

  # ------------------------------------------------------------------------
  # Dataset overview (Data tab additions)
  # ------------------------------------------------------------------------
  output$dataset_overview_summary <- renderUI({
    s <- state$expr_summary
    if (is.null(s)) return(NULL)
    md <- state$metadata
    md_cols <- colnames(md %||% data.frame())

    ct_col <- if (input$meta_cell_type_col %in% md_cols) {
      input$meta_cell_type_col
    } else NA_character_
    src_col <- if (input$meta_source_col %in% md_cols) {
      input$meta_source_col
    } else NA_character_
    n_ct <- if (!is.na(ct_col)) length(unique(na.omit(md[[ct_col]]))) else NA_integer_
    n_src <- if (!is.na(src_col)) length(unique(na.omit(md[[src_col]]))) else NA_integer_

    # Auto-detect patient / sample columns from the metadata so the user
    # doesn't have to wire up another dropdown for these counts.
    patient_col <- detect_patient_column(md)
    sample_col  <- detect_sample_column(md)
    n_patients  <- count_distinct_meta(md, patient_col)
    n_samples   <- count_distinct_meta(md, sample_col)

    emb_avail <- list_available_embeddings(state$expression)
    if (!length(emb_avail) && !is.null(state$metadata)) {
      meta_emb_names <- names(detect_metadata_embeddings(state$metadata))
      if (length(meta_emb_names)) {
        emb_avail <- paste0(meta_emb_names, " (metadata)")
      }
    }

    # Show which column was used for patient / sample tile so the user can
    # sanity-check the auto-detection.
    patient_caption <- if (is.na(n_patients)) NULL
                       else tags$div(class = "stat-caption",
                                     sprintf("from %s", patient_col))
    sample_caption  <- if (is.na(n_samples))  NULL
                       else tags$div(class = "stat-caption",
                                     sprintf("from %s", sample_col))

    # Stat tiles flow in a single CSS-grid row so all tiles share the full
    # width of the overview card equally. Tiles with potentially longer
    # text content (Input kind, Available embeddings) get `stat-box-wide`,
    # which makes them span 2 grid cells so their text has room without
    # forcing the count tiles to be wider than they need to be.
    tags$div(
      class = "stat-row",
      # "Input: …" (e.g. Seurat / SingleCellExperiment / matrix / AnnData)
      # was formerly shown in the "Loaded data" panel — surfaced here so
      # the Dataset overview is a single, complete summary of the input.
      div(class = "stat-box stat-box-wide",
          tags$div(class = "stat-label", "Input"),
          tags$div(class = "stat-value-sm stat-value-input",
                   s$kind %||% "—")),
      div(class = "stat-box",
          tags$div(class = "stat-label", "Cells"),
          tags$div(class = "stat-value", .fmt_int(s$n_samples %||% 0))),
      div(class = "stat-box",
          tags$div(class = "stat-label", "Genes"),
          tags$div(class = "stat-value", .fmt_int(s$n_genes %||% 0))),
      div(class = "stat-box",
          tags$div(class = "stat-label", "Patients"),
          tags$div(class = "stat-value",
                   if (is.na(n_patients)) "—" else .fmt_int(n_patients)),
          patient_caption),
      div(class = "stat-box",
          tags$div(class = "stat-label", "Samples"),
          tags$div(class = "stat-value",
                   if (is.na(n_samples)) "—" else .fmt_int(n_samples)),
          sample_caption),
      div(class = "stat-box",
          tags$div(class = "stat-label", "Cell types"),
          tags$div(class = "stat-value",
                   if (is.na(n_ct)) "—" else .fmt_int(n_ct))),
      div(class = "stat-box",
          tags$div(class = "stat-label", "Sources / groups"),
          tags$div(class = "stat-value",
                   if (is.na(n_src)) "—" else .fmt_int(n_src))),
      div(class = "stat-box stat-box-wide",
          tags$div(class = "stat-label", "Available embeddings"),
          tags$div(class = "stat-value-sm",
                   if (length(emb_avail)) paste(emb_avail, collapse = ", ")
                   else "(none detected on this object)"))
    )
  })

  output$celltype_count_plot <- renderPlot({
    md <- state$metadata
    ct_col <- input$meta_cell_type_col
    req(!is.null(md), nzchar(ct_col), ct_col != "(none)", ct_col %in% colnames(md))
    df <- as.data.frame(table(md[[ct_col]], useNA = "no"))
    colnames(df) <- c("cell_type", "n")
    df <- df[order(-df$n), ]
    df$cell_type <- factor(df$cell_type, levels = df$cell_type)
    pal <- tryCatch(PhenoMapR::get_celltype_palette(as.character(df$cell_type)),
                    error = function(e) NULL)
    p <- ggplot(df, aes(x = cell_type, y = n, fill = cell_type)) +
      .geom_rounded_col() +
      labs(x = NULL, y = "Cells", fill = "Cell type") +
      theme_minimal(base_size = .theme_base_size()) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1),
            legend.position = "none")
    if (!is.null(pal)) p <- p + scale_fill_manual(values = pal)
    panel_objects$celltype_count_plot <- p
    p
  })
  phenomapr_register_plot_download(output, "celltype_count_plot",
    function() isolate(panel_objects$celltype_count_plot),
    width = 8, height = 5)

  output$source_count_plot <- renderPlot({
    md <- state$metadata
    src_col <- input$meta_source_col
    req(!is.null(md), nzchar(src_col), src_col != "(none)", src_col %in% colnames(md))
    # as.character() first: integer FOV/core IDs would otherwise stay numeric
    # and trip discrete_scale ("Continuous value supplied to a discrete scale").
    df <- as.data.frame(table(as.character(md[[src_col]]), useNA = "no"),
                        stringsAsFactors = FALSE)
    colnames(df) <- c("source", "n")
    df <- df[order(-df$n), ]
    df$source <- factor(df$source, levels = df$source)
    p <- ggplot(df, aes(x = source, y = n, fill = source)) +
      .geom_rounded_col() +
      scale_fill_phenomapr_d() +
      labs(x = NULL, y = "Cells", fill = "Source") +
      theme_minimal(base_size = .theme_base_size()) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1),
            legend.position = "none")
    panel_objects$source_count_plot <- p
    p
  })
  phenomapr_register_plot_download(output, "source_count_plot",
    function() isolate(panel_objects$source_count_plot),
    width = 8, height = 5)

  output$celltype_source_plot <- renderPlot({
    md <- state$metadata
    ct_col <- input$meta_cell_type_col
    src_col <- input$meta_source_col
    req(!is.null(md),
        nzchar(ct_col), ct_col != "(none)", ct_col %in% colnames(md),
        nzchar(src_col), src_col != "(none)", src_col %in% colnames(md))
    df <- data.frame(
      cell_type = as.character(md[[ct_col]]),
      source    = as.character(md[[src_col]]),
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$cell_type) & !is.na(df$source), ]
    if (!nrow(df)) {
      panel_objects$celltype_source_plot <- NULL
      return(NULL)
    }
    df_count <- df %>% dplyr::count(cell_type, source)
    p <- ggplot(df_count, aes(x = cell_type, y = n, fill = source)) +
      .geom_rounded_stack() +
      scale_fill_phenomapr_d() +
      labs(x = NULL, y = "Cells", fill = "Source") +
      theme_minimal(base_size = .theme_base_size()) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
    panel_objects$celltype_source_plot <- p
    p
  })
  phenomapr_register_plot_download(output, "celltype_source_plot",
    function() isolate(panel_objects$celltype_source_plot),
    width = 9, height = 5)

  output$metadata_columns_tbl <- renderDT({
    md <- state$metadata
    req(!is.null(md))
    df <- data.frame(
      column = colnames(md),
      type   = vapply(md, function(x) class(x)[1L], character(1L)),
      unique = vapply(md, function(x) length(unique(x)), integer(1L)),
      n_na   = vapply(md, function(x) sum(is.na(x)), integer(1L)),
      stringsAsFactors = FALSE
    )
    panel_objects$metadata_columns_tbl <- df
    datatable(df, rownames = FALSE,
              options = list(pageLength = 8, dom = "tip"))
  })
  phenomapr_register_table_download(output, "metadata_columns_tbl",
    function() isolate(panel_objects$metadata_columns_tbl))

  # ----- Matrix diagnostics + cleanup ----------------------------------------
  # When the uploaded expression is a plain matrix / data.frame (i.e.
  # kind == "matrix"), run detect_expression_format() to surface gene-ID
  # style (HUGO / ENSG / mixed), expression scale (raw counts / CPM/TPM
  # / log-normalized / z-scaled), and a single-cell vs bulk guess driven
  # by sparsity. The result is cached in expr_diagnostics() so the
  # downstream UI and the "Clean & normalize" observer share one
  # detection pass.
  expr_diagnostics <- reactive({
    s <- state$expr_summary
    if (is.null(s) || !identical(s$kind, "matrix")) return(NULL)
    expr <- state$expression
    if (is.null(expr)) return(NULL)
    # Suppress the "NAs introduced by coercion" warning that some
    # S4 .Internal( ) paths emit when vapply()-coercing 1x1 sparse
    # returns to numeric(1). The detector itself does the right
    # thing now (samples @x directly for sparseMatrix inputs); the
    # withCallingHandlers wrapper is belt-and-braces so a stray
    # warning never bubbles into the Shiny log.
    out <- tryCatch(
      withCallingHandlers(
        detect_expression_format(expr, verbose = FALSE),
        warning = function(w) {
          if (grepl("NAs introduced by coercion",
                    conditionMessage(w), fixed = TRUE)) {
            invokeRestart("muffleWarning")
          }
        }
      ),
      error = function(e) {
        showNotification(
          paste0("Could not auto-detect expression format: ",
                 conditionMessage(e)),
          type = "warning", duration = 8
        )
        NULL
      }
    )
    out
  })

  output$expr_matrix_diagnostics <- renderUI({
    d <- expr_diagnostics()
    if (is.null(d)) return(NULL)

    fmt_class <- switch(
      d$format,
      raw_counts = "diag-warn",
      cpm_or_tpm = "diag-warn",
      z_scaled   = "diag-warn",
      log_normalized = "diag-ok",
      "diag-info"
    )
    id_class <- switch(
      d$gene_id_kind,
      hugo = "diag-ok",
      ensembl = "diag-warn",
      mixed = "diag-warn",
      "diag-info"
    )
    id_label <- switch(
      d$gene_id_kind,
      hugo = sprintf("HUGO symbols (%d/%d)",
                     d$n_hugo_like, d$n_genes),
      ensembl = sprintf("Ensembl IDs (%d/%d ENSG-prefixed)",
                        d$n_ensembl, d$n_genes),
      mixed = sprintf("Mixed (%d ENSG + %d HUGO-like)",
                      d$n_ensembl, d$n_hugo_like),
      "Unknown gene-ID style"
    )
    sc_class <- switch(
      d$sc_or_bulk,
      single_cell = "diag-info",
      bulk        = "diag-info",
      "diag-info"
    )

    diag_row <- function(class_, label, value, detail = NULL) {
      tags$div(
        class = paste("diag-row", class_),
        tags$div(class = "diag-label", label),
        tags$div(class = "diag-value", value),
        if (!is.null(detail))
          tags$div(class = "diag-detail", detail)
      )
    }

    stats_line <- if (all(is.finite(d$stats[c("min", "max", "mean")]))) {
      sprintf(
        "min %.2f - max %.2f - mean %.2f - %.0f%% integer - %.0f%% zeros",
        d$stats[["min"]], d$stats[["max"]], d$stats[["mean"]],
        d$stats[["frac_integer"]] * 100, d$stats[["frac_zero"]] * 100
      )
    } else NULL

    cleanup_block <- NULL
    needs_cleanup <- d$gene_id_kind %in% c("ensembl", "mixed") ||
      d$n_dup > 0L ||
      d$format %in% c("raw_counts", "cpm_or_tpm")
    if (needs_cleanup) {
      default_mode <- if (identical(d$sc_or_bulk, "single_cell"))
        "single_cell" else "bulk"
      cleanup_block <- tags$div(
        class = "diag-cleanup",
        tags$div(class = "diag-cleanup-title",
                 icon("wand-magic-sparkles"),
                 tags$strong(" Cleanup options")),
        checkboxInput(
          "diag_do_hugo",
          "Clean gene IDs to approved HUGO symbols (HGNChelper)",
          value = d$gene_id_kind %in% c("ensembl", "mixed") ||
            d$n_dup > 0L
        ),
        checkboxInput(
          "diag_do_collapse",
          "Collapse duplicate gene rows by mean",
          value = d$n_dup > 0L
        ),
        checkboxInput(
          "diag_do_normalize",
          "Log-normalize expression",
          value = d$format %in% c("raw_counts", "cpm_or_tpm")
        ),
        radioButtons(
          "diag_norm_mode",
          "Normalization mode",
          choices = c("Auto-detect (recommended)" = "auto",
                      "Single cell (Seurat log1p of library-size-scaled counts)" = "single_cell",
                      "Bulk (log2(CPM + 1))" = "bulk"),
          selected = if (identical(d$sc_or_bulk, "unclear"))
            "auto" else default_mode,
          inline = FALSE
        ),
        actionButton(
          "diag_run_clean",
          tagList(icon("broom"), " Clean & normalize"),
          class = "btn-primary btn-sm",
          width = "100%"
        )
      )
    } else {
      cleanup_block <- tags$div(
        class = "diag-cleanup diag-cleanup-clean",
        icon("check-circle"),
        tags$span(" Matrix looks ready -- no cleanup recommended.")
      )
    }

    tags$div(
      class = "expr-matrix-diag",
      tags$div(class = "diag-title",
               icon("microscope"),
               tags$strong(" Matrix diagnostics")),
      diag_row(id_class, "Gene IDs",
               id_label,
               detail = if (d$n_dup > 0L)
                 sprintf("%d duplicate gene ID(s)%s",
                         d$n_dup,
                         if (length(d$dup_examples))
                           paste0(" (e.g. ",
                                  paste(d$dup_examples, collapse = ", "), ")")
                         else "")),
      diag_row(fmt_class, "Expression format",
               sprintf("%s (%s confidence)",
                       d$format_label, d$format_confidence),
               detail = stats_line),
      diag_row(sc_class, "Data type",
               d$sc_or_bulk_label),
      cleanup_block,
      if (length(d$recommendations))
        tags$ul(class = "diag-recs",
                lapply(d$recommendations, tags$li))
    )
  })

  # When the user clicks "Clean & normalize", run clean_matrix_input() on
  # state$expression with the chosen options and update state in place.
  # Recomputes state$expr_summary so the rest of the UI (Score-tab
  # status, embedding, gene-coverage) sees the cleaned object.
  observeEvent(input$diag_run_clean, {
    req(state$expression)
    if (!identical(state$expr_summary$kind %||% "", "matrix")) {
      showNotification("Cleanup only applies to matrix / data.frame uploads.",
                       type = "warning", duration = 5)
      return()
    }
    phenomapr_busy_show("Cleaning matrix...",
                        "Running HUGO cleanup and / or normalization")
    # See the expression-upload observer for the deferred-work
    # rationale: heavy synchronous work MUST run inside later::later()
    # so the busy_show custom message reaches the browser BEFORE
    # libuv is blocked by the work.
    expr_snapshot <- state$expression
    do_hugo          <- isTRUE(input$diag_do_hugo)
    do_collapse_dups <- isTRUE(input$diag_do_collapse)
    do_normalize     <- isTRUE(input$diag_do_normalize)
    norm_mode        <- input$diag_norm_mode %||% "auto"
    sess <- session
    later::later(function() {
      res <- tryCatch(
        clean_matrix_input(
          expr_snapshot,
          do_hugo          = do_hugo,
          do_collapse_dups = do_collapse_dups,
          do_normalize     = do_normalize,
          mode             = norm_mode,
          hugo_species     = "human",
          verbose          = FALSE
        ),
        error = function(e) {
          showNotification(paste0("Cleanup failed: ",
                                  conditionMessage(e)),
                           type = "error", duration = 8, session = sess)
          NULL
        }
      )
      phenomapr_busy_hide(session = sess)
      if (is.null(res)) return()
      later::later(function() {
        tryCatch({
          state$expression <- res$matrix
          state$expr_summary <- summarize_expression_object(res$matrix)
          state$expr_summary$notes <- paste(
            "Matrix cleaned via clean_matrix_input():",
            paste(res$steps, collapse = "; ")
          )
          shiny::showNotification(
            paste0(
              "Cleanup complete -- ",
              if (length(res$steps)) paste(res$steps, collapse = "; ")
              else "no changes were necessary."
            ),
            type = "message", duration = 8, session = sess
          )
        }, error = function(e) {
          shiny::showNotification(
            paste0("Internal error after cleanup: ", conditionMessage(e)),
            type = "error", duration = 10, session = sess
          )
        })
      }, delay = 0)
    }, delay = 0)
  })

  # ------------------------------------------------------------------------
  # 2. Reference
  # ------------------------------------------------------------------------
  # Populate cancer_type dropdown when reference changes.
  # ICI PRECOG ships pre-filtered to |z| >= 1 and several cohorts have
  # zero genes above 2; auto-set the Signature |z| slider to 1 so
  # "Compute PhenoMapR scores" does not return an empty score table.
  observeEvent(input$reference_choice, {
    if (input$reference_choice == "_custom") return()
    cts <- get_cancer_types(input$reference_choice)
    updateSelectInput(session, "cancer_type", choices = cts, selected = cts[1L])
    if (identical(input$reference_choice, "ici_precog")) {
      updateSliderInput(session, "z_score_cutoff", value = 1)
    } else {
      # Restore the package default when leaving ICI, but only if the
      # user is still sitting on the ICI-recommended value of 1.
      cut_now <- suppressWarnings(as.numeric(input$z_score_cutoff))
      if (is.finite(cut_now) && isTRUE(all.equal(cut_now, 1))) {
        updateSliderInput(session, "z_score_cutoff", value = 2)
      }
    }
  })

  # Custom: upload file
  observeEvent(custom_ref_file_pick(), {
    pick <- custom_ref_file_pick()
    if (is.null(pick)) {
      # Only clear custom-signature state; leave PRECOG/TCGA selections alone.
      if (identical(state$reference_source %||% "", "custom") ||
          grepl("^Custom \\(", state$reference_label %||% "")) {
        state$reference <- NULL
        state$reference_label <- ""
      }
      return()
    }
    phenomapr_busy_show("Loading custom signature...", pick$name)
    # Defer the heavy synchronous read so busy_show reaches the
    # browser before libuv is blocked.
    pick_local <- pick
    sess <- session
    later::later(function() {
      ref <- tryCatch(
        parse_reference_upload(pick_local$datapath, pick_local$name),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 8,
                           session = sess)
          NULL
        }
      )
      phenomapr_busy_hide(session = sess)
      if (is.null(ref)) return()
      later::later(function() {
        state$reference <- ref
        state$reference_label <- paste0("Custom (", pick_local$name, ")")
        showNotification(
          sprintf("Custom reference loaded (%d genes).", nrow(ref)),
          type = "message", duration = 4, session = sess
        )
      }, delay = 0)
    }, delay = 0)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  # Custom: derive from bulk + phenotype
  observeEvent(derive_bulk_file_pick(), {
    pick <- derive_bulk_file_pick()
    if (is.null(pick)) {
      state$derive_bulk <- NULL
      return()
    }
    phenomapr_busy_show("Loading bulk expression...", pick$name)
    pick_local <- pick
    sess <- session
    later::later(function() {
      res <- tryCatch(
        parse_expression_upload(pick_local$datapath, pick_local$name),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 8,
                           session = sess); NULL
        }
      )
      phenomapr_busy_hide(session = sess)
      if (is.null(res)) return()
      later::later(function() {
        state$derive_bulk <- res$object
        showNotification(
          sprintf("Bulk expression: %s genes \u00d7 %s samples.",
                  .fmt_int(res$n_genes %||% 0L), .fmt_int(res$n_samples %||% 0L)),
          type = "message", duration = 5, session = sess
        )
      }, delay = 0)
    }, delay = 0)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  output$derive_bulk_summary <- renderUI({
    obj <- state$derive_bulk
    if (is.null(obj)) return(NULL)
    summ <- summarize_expression_object(obj)
    tags$div(
      class = "help-block",
      tags$strong("Bulk loaded: "),
      sprintf("%s genes × %s samples (%s).",
              .fmt_int(summ$n_genes %||% 0L),
              .fmt_int(summ$n_samples %||% 0L),
              summ$kind)
    )
  })

  output$derive_pheno_preview <- renderUI({
    df <- state$derive_phen
    col <- input$derive_pheno_col
    if (is.null(df) || is.null(col) || !nzchar(col) || !(col %in% colnames(df))) return(NULL)
    vals <- df[[col]]
    tbl <- if (is.numeric(vals)) {
      qs <- stats::quantile(vals, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
      tags$div(
        tags$p(tags$strong("Outcome column: "), col, " (numeric)"),
        tags$p(tags$strong("min · Q1 · median · Q3 · max: "),
               paste(.fmt_num(qs), collapse = " · "))
      )
    } else {
      tbl <- as.data.frame(table(as.character(vals), useNA = "ifany"))
      colnames(tbl) <- c("level", "n")
      tbl$n <- .fmt_int(tbl$n)
      tags$div(
        tags$p(tags$strong("Outcome column: "), col, " (categorical)"),
        tag_table(tbl)
      )
    }
    tags$div(class = "help-block", tbl)
  })

  observeEvent(derive_phen_file_pick(), {
    pick <- derive_phen_file_pick()
    if (is.null(pick)) {
      state$derive_phen <- NULL
      updateSelectInput(session, "derive_id_col", choices = character(0))
      updateSelectInput(session, "derive_pheno_col", choices = character(0))
      updateSelectInput(session, "derive_time_col", choices = character(0))
      updateSelectInput(session, "derive_event_col", choices = character(0))
      return()
    }
    phenomapr_busy_show("Loading phenotype table...", pick$name)
    pick_local <- pick
    sess <- session
    later::later(function() {
      df <- tryCatch(
        parse_metadata_upload(pick_local$datapath, pick_local$name),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 8,
                           session = sess); NULL
        }
      )
      phenomapr_busy_hide(session = sess)
      if (is.null(df)) return()
      later::later(function() {
        state$derive_phen <- df
        cols <- colnames(df)
        updateSelectInput(sess, "derive_id_col", choices = cols, selected = cols[1L])
        updateSelectInput(sess, "derive_pheno_col", choices = cols,
                          selected = if (length(cols) > 1L) cols[2L] else cols[1L])
        updateSelectInput(sess, "derive_time_col", choices = cols)
        updateSelectInput(sess, "derive_event_col", choices = cols)
      }, delay = 0)
    }, delay = 0)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(derive_probe_file_pick(), {
    pick <- derive_probe_file_pick()
    if (is.null(pick)) {
      state$derive_probe_annot <- NULL
      return()
    }
    phenomapr_busy_show("Loading GPL probe annotation...", pick$name)
    pick_local <- pick
    sess <- session
    later::later(function() {
      df <- tryCatch(
        parse_metadata_upload(pick_local$datapath, pick_local$name),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 8,
                           session = sess); NULL
        }
      )
      phenomapr_busy_hide(session = sess)
      if (is.null(df)) return()
      later::later(function() {
        state$derive_probe_annot <- df
        showNotification(
          sprintf("Probe annotation loaded (%d rows).", nrow(df)),
          type = "message", duration = 4, session = sess
        )
      }, delay = 0)
    }, delay = 0)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  output$derive_meta_z_study_status <- renderText({
    n <- length(state$derive_meta_z_studies %||% list())
    if (n == 0L) {
      "Meta-z cohort list: (empty — add at least one bulk + phenotype cohort)"
    } else {
      sprintf("Meta-z cohort list: %d cohort(s) queued.", n)
    }
  })

  observeEvent(input$derive_add_study, {
    req(state$derive_bulk, state$derive_phen)
    if (input$derive_platform == "microarray" && is.null(state$derive_probe_annot)) {
      showNotification(
        "Microarray cohorts need a GPL probe annotation before adding to the meta-z list.",
        type = "warning", duration = 8
      )
      return()
    }
    n <- length(state$derive_meta_z_studies %||% list()) + 1L
    label <- sprintf("cohort_%02d", n)
    state$derive_meta_z_studies[[label]] <- list(
      bulk_expression = state$derive_bulk,
      phenotype = state$derive_phen,
      sample_id_column = input$derive_id_col,
      phenotype_column = if (input$derive_type == "survival") NULL else input$derive_pheno_col,
      phenotype_type = input$derive_type,
      survival_time = if (input$derive_type == "survival") input$derive_time_col else NULL,
      survival_event = if (input$derive_type == "survival") input$derive_event_col else NULL,
      normalize = isTRUE(input$derive_normalize),
      platform = input$derive_platform %||% "rnaseq",
      probe_annotation = state$derive_probe_annot
    )
    showNotification(
      sprintf("Added %s to meta-z list (%d cohort(s) total).", label, n),
      type = "message", duration = 4
    )
  })

  observeEvent(input$derive_clear_studies, {
    state$derive_meta_z_studies <- list()
    showNotification("Meta-z cohort list cleared.", type = "message", duration = 3)
  })

  observeEvent(input$derive_run, {
    meta_z_mode <- isTRUE(input$derive_meta_z_mode)
    if (meta_z_mode) {
      req(length(state$derive_meta_z_studies %||% list()) >= 1L)
    } else {
      req(state$derive_bulk, state$derive_phen)
    }
    phenomapr_busy_show(
      if (meta_z_mode) "Deriving meta-z signature..." else "Deriving phenotype signature...",
      sprintf("Bulk + phenotype | %s outcome", input$derive_type %||% "binary")
    )
    bin_pos <- if (input$derive_binary_positive %in% c("first", "second")) {
      input$derive_binary_positive
    } else "second"
    bulk_snapshot  <- state$derive_bulk
    phen_snapshot  <- state$derive_phen
    probe_snapshot <- state$derive_probe_annot
    studies_snapshot <- state$derive_meta_z_studies
    derive_args <- list(
      sample_id_column = input$derive_id_col,
      phenotype_column = if (input$derive_type == "survival") NULL else input$derive_pheno_col,
      phenotype_type   = input$derive_type,
      survival_time    = if (input$derive_type == "survival") input$derive_time_col else NULL,
      survival_event   = if (input$derive_type == "survival") input$derive_event_col else NULL,
      normalize        = isTRUE(input$derive_normalize),
      platform         = input$derive_platform %||% "rnaseq",
      hugo_species     = input$derive_hugo_species %||% "human",
      binary_positive_reference = bin_pos,
      meta_z_mode      = meta_z_mode
    )
    sess <- session
    later::later(function() {
      ref <- tryCatch(
        if (isTRUE(derive_args$meta_z_mode)) {
          studies <- studies_snapshot
          if (!length(studies) && !is.null(bulk_snapshot) && !is.null(phen_snapshot)) {
            studies <- list(cohort_01 = list(
              bulk_expression = bulk_snapshot,
              phenotype = phen_snapshot,
              sample_id_column = derive_args$sample_id_column,
              phenotype_column = derive_args$phenotype_column,
              phenotype_type = derive_args$phenotype_type,
              survival_time = derive_args$survival_time,
              survival_event = derive_args$survival_event,
              normalize = derive_args$normalize,
              platform = derive_args$platform,
              probe_annotation = probe_snapshot
            ))
          }
          PhenoMapR::derive_meta_z_from_bulk_studies(
            studies = studies,
            meta_z_label = "custom_meta_z",
            hugo_species = derive_args$hugo_species,
            binary_positive_reference = derive_args$binary_positive_reference,
            verbose = TRUE
          )
        } else {
          if (derive_args$platform == "microarray" && is.null(probe_snapshot) &&
              !is.null(bulk_snapshot)) {
            rn <- rownames(bulk_snapshot)
            if (length(rn) && mean(grepl("^(AFFX-|ILMN_)", rn, ignore.case = TRUE)) >= 0.5) {
              stop("Microarray probe-level bulk data requires a GPL probe annotation file.")
            }
          }
          PhenoMapR::derive_reference_from_bulk(
            bulk_expression  = bulk_snapshot,
            phenotype        = phen_snapshot,
            sample_id_column = derive_args$sample_id_column,
            phenotype_column = derive_args$phenotype_column,
            phenotype_type   = derive_args$phenotype_type,
            survival_time    = derive_args$survival_time,
            survival_event   = derive_args$survival_event,
            normalize        = derive_args$normalize,
            platform         = derive_args$platform,
            probe_annotation = probe_snapshot,
            hugo_species     = derive_args$hugo_species,
            binary_positive_reference = derive_args$binary_positive_reference,
            verbose          = TRUE
          )
        },
        error = function(e) {
          showNotification(paste0("Reference derivation failed: ",
                                  conditionMessage(e)),
                           type = "error", duration = 10, session = sess); NULL
        }
      )
      phenomapr_busy_hide(session = sess)
      if (is.null(ref)) {
        return()
      }
      derive_type_local <- derive_args$phenotype_type
      later::later(function() {
        state$reference <- ref
        state$reference_label <- if (isTRUE(derive_args$meta_z_mode)) {
          sprintf(
            "Custom meta-z (%d cohorts)",
            attr(ref, "n_studies") %||% length(studies_snapshot)
          )
        } else {
          sprintf(
            "Custom (derived from bulk + phenotype, type = %s, platform = %s)",
            attr(ref, "phenotype_type") %||% derive_type_local,
            attr(ref, "platform") %||% derive_args$platform
          )
        }
        showNotification(
          sprintf("Reference derived (%s genes).", .fmt_int(nrow(ref))),
          type = "message", duration = 5, session = sess
        )
      }, delay = 0)
    }, delay = 0)
  })

  # ---- Derived reference preview ------------------------------------------
  output$derived_reference_summary <- renderUI({
    ref <- state$reference
    if (is.null(ref) || !is.data.frame(ref)) {
      return(tags$p(tags$em(
        "Upload bulk expression + phenotype, then click \"Derive reference\" to populate this section."
      )))
    }
    z <- as.numeric(ref[[1L]])
    n_pos <- sum(z > 2, na.rm = TRUE)
    n_neg <- sum(z < -2, na.rm = TRUE)
    rng <- range(z, na.rm = TRUE)
    tagList(
      tags$p(tags$strong("Reference label: "), state$reference_label %||% "Custom"),
      tags$p(tags$strong("Genes (rows): "), .fmt_int(nrow(ref))),
      tags$p(tags$strong("z-score column: "), colnames(ref)[1L]),
      tags$p(tags$strong("|z| > 2 genes: "),
             sprintf("%s positive, %s negative", .fmt_int(n_pos), .fmt_int(n_neg))),
      tags$p(tags$strong("z range: "),
             sprintf("[%s, %s]", .fmt_num(rng[1L]), .fmt_num(rng[2L])))
    )
  })

  output$derived_top_pos_tbl <- renderDT({
    ref <- state$reference
    req(is.data.frame(ref))
    df <- data.frame(
      gene = rownames(ref),
      z_score = as.numeric(ref[[1L]]),
      stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$z_score) & df$z_score > 0, , drop = FALSE]
    df <- df[order(-df$z_score), , drop = FALSE]
    df <- head(df, 50L)
    df$z_score <- round(df$z_score, 3)
    panel_objects$derived_top_pos_tbl <- df
    datatable(df, rownames = FALSE,
              options = list(pageLength = 10, dom = "tip", scrollX = TRUE))
  })
  phenomapr_register_table_download(output, "derived_top_pos_tbl",
    function() isolate(panel_objects$derived_top_pos_tbl))

  output$derived_top_neg_tbl <- renderDT({
    ref <- state$reference
    req(is.data.frame(ref))
    df <- data.frame(
      gene = rownames(ref),
      z_score = as.numeric(ref[[1L]]),
      stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$z_score) & df$z_score < 0, , drop = FALSE]
    df <- df[order(df$z_score), , drop = FALSE]
    df <- head(df, 50L)
    df$z_score <- round(df$z_score, 3)
    panel_objects$derived_top_neg_tbl <- df
    datatable(df, rownames = FALSE,
              options = list(pageLength = 10, dom = "tip", scrollX = TRUE))
  })
  phenomapr_register_table_download(output, "derived_top_neg_tbl",
    function() isolate(panel_objects$derived_top_neg_tbl))

  output$download_derived_reference <- downloadHandler(
    filename = function() {
      sprintf("phenomapr_derived_reference_%s.rds",
              format(Sys.time(), "%Y%m%d_%H%M%S"))
    },
    content = function(file) {
      req(state$reference)
      ref <- state$reference
      if (!is.data.frame(ref)) {
        stop("No derived reference data.frame available. Derive one first.")
      }
      saveRDS(ref, file)
    }
  )

  # Sync built-in reference selection (no upload required)
  observe({
    rc <- input$reference_choice
    if (is.null(rc) || !nzchar(rc) || rc == "_custom") return()
    state$reference <- rc
    state$reference_label <- sprintf("%s · %s", rc, input$cancer_type %||% "?")
  })

  output$reference_status <- renderUI({
    ref <- state$reference
    if (is.null(ref)) return(tags$p(tags$em("No reference selected.")))
    # Resolve the live |z| threshold once so both the "|z| cutoff"
    # row and the "Significant genes" row stay in sync as the user
    # scrubs the slider. NA / non-finite values fall back to "(n/a)".
    cut_raw <- suppressWarnings(as.numeric(input$z_score_cutoff))
    cut_val <- if (length(cut_raw) == 1L && is.finite(cut_raw)) {
      cut_raw
    } else {
      NA_real_
    }
    cut_lbl <- if (is.finite(cut_val)) format(cut_val) else "(n/a)"
    # Count signature genes that pass |z| >= cutoff. We surface this
    # alongside the cutoff itself so the user can immediately see how
    # restrictive the current threshold is for the chosen reference;
    # bumping the slider one notch and watching this number react is
    # often the fastest way to settle on a sensible value.
    n_sig <- NA_integer_
    if (is.character(ref)) {
      ct <- input$cancer_type %||% NA_character_
      valid_cts <- tryCatch(get_cancer_types(ref), error = function(e) character())
      if (is.finite(cut_val) && !is.na(ct) && nzchar(ct) && ct %in% valid_cts) {
        # `n = Inf` returns the entire signature; we then count
        # entries with |z_score| >= cutoff. This is the same gene
        # set that `phenomap()` will pre-filter on for scoring, so
        # the number is exactly "how many genes are about to drive
        # the score" -- not a paraphrase.
        n_sig <- tryCatch({
          tg <- PhenoMapR::get_top_prognostic_genes(
            reference   = ref,
            cancer_type = ct,
            n           = Inf,
            direction   = "both"
          )
          z_vec <- suppressWarnings(as.numeric(tg$z_score))
          sum(is.finite(z_vec) & abs(z_vec) >= cut_val)
        }, error = function(e) NA_integer_)
      }
      tagList(
        tags$p(tags$strong("Built-in: "), ref),
        tags$p(tags$strong("Cancer type: "), input$cancer_type %||% "(none)"),
        tags$p(tags$strong("|z| cutoff: "), cut_lbl),
        tags$p(
          tags$strong("Significant genes: "),
          if (is.na(n_sig)) tags$em("(set cancer type & cutoff)")
          else .fmt_int(n_sig)
        )
      )
    } else {
      # Custom signature: count rows in the user's reference whose
      # absolute z-score is at or above the cutoff.
      n_sig <- if (is.finite(cut_val)) {
        z_col <- which(vapply(as.data.frame(ref), is.numeric, logical(1)))[1L]
        if (is.na(z_col)) NA_integer_ else {
          z_vec <- suppressWarnings(as.numeric(as.data.frame(ref)[[z_col]]))
          sum(is.finite(z_vec) & abs(z_vec) >= cut_val)
        }
      } else {
        NA_integer_
      }
      tagList(
        tags$p(tags$strong("Source: "), state$reference_label),
        tags$p(tags$strong("Genes: "), .fmt_int(nrow(ref))),
        tags$p(tags$strong("z column: "), colnames(ref)[1L]),
        tags$p(tags$strong("|z| cutoff: "), cut_lbl),
        tags$p(
          tags$strong("Significant genes: "),
          if (is.na(n_sig)) tags$em("(set a numeric z column & cutoff)")
          else .fmt_int(n_sig)
        )
      )
    }
  })

  output$gene_coverage_tbl <- renderDT({
    req(state$expr_summary, state$reference)
    genes <- state$expr_summary$gene_ids
    if (is.character(state$reference)) {
      cov <- PhenoMapR::get_gene_coverage(genes, reference = state$reference)
    } else {
      cov <- data.frame(
        reference = "custom",
        total_genes = length(genes),
        covered_genes = length(intersect(genes, rownames(state$reference))),
        coverage_pct = round(100 * length(intersect(genes, rownames(state$reference))) /
                               max(1L, length(genes)), 2)
      )
    }
    panel_objects$gene_coverage_tbl <- cov
    datatable(cov, rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })
  phenomapr_register_table_download(output, "gene_coverage_tbl",
    function() isolate(panel_objects$gene_coverage_tbl))

  # Phenotype-signature histogram.
  #
  # In both built-in (PRECOG / TCGA / Pediatric PRECOG / ICI PRECOG)
  # and custom (data.frame / matrix returned by
  # derive_reference_from_bulk()) modes, the panel now shows the
  # SAME thing: a histogram of every gene\'s reference z-score in the
  # currently selected signature. Two material differences:
  #
  #   * Built-in references ship pre-filtered to keep the install
  #     lightweight (PRECOG / TCGA / Pediatric PRECOG -> |z| >= 2;
  #     ICI PRECOG -> |z| >= 1), so the histogram has nothing to
  #     show in the |z| < cutoff interior. We facet by sign and let
  #     each tail use its own x scale, which produces a "broken
  #     axis" layout that omits the empty middle entirely (the user
  #     asked us to "simplify" the built-in histogram by not showing
  #     the empty interior -- this is the cleanest way to do that
  #     without hand-rolling ggh4x dependencies).
  #
  #   * Custom signatures often span the full distribution
  #     including |z| < 2 (e.g. genes with weak survival
  #     associations), so we render them as a single, continuous
  #     histogram with no x-axis trimming.
  #
  # In both cases, the input$z_score_cutoff slider in the sidebar
  # drives a pair of dashed vertical lines at +cut and -cut
  # superimposed on the histogram, mirroring the percentile-cutoff
  # overlay the Score-distribution panel uses. Updating the slider
  # reactively redraws the lines without recomputing the histogram
  # data, so the user can scrub through cutoff values and see the
  # corresponding gene partition immediately.
  output$reference_signature_plot <- renderPlot({
    req(state$reference)
    cut_raw <- suppressWarnings(as.numeric(input$z_score_cutoff))
    cut_val <- if (length(cut_raw) == 1L && is.finite(cut_raw)) {
      cut_raw
    } else {
      NA_real_
    }
    cut_label <- if (is.finite(cut_val)) {
      sprintf("|z| = %g", cut_val)
    } else {
      NA_character_
    }
    base_size <- .theme_base_size()

    if (is.data.frame(state$reference) || is.matrix(state$reference)) {
      # Custom signature: full continuous histogram of every gene\'s
      # z-score. The first numeric column is the z-score column by
      # convention (derive_reference_from_bulk() / the upload path
      # both emit a single-column data.frame).
      ref_df <- as.data.frame(state$reference)
      z_col <- which(vapply(ref_df, is.numeric, logical(1)))[1L]
      if (is.na(z_col)) {
        panel_objects$reference_signature_plot <- NULL
        return(NULL)
      }
      z_vec <- as.numeric(ref_df[[z_col]])
      z_vec <- z_vec[is.finite(z_vec)]
      if (!length(z_vec)) {
        panel_objects$reference_signature_plot <- NULL
        return(NULL)
      }
      df <- data.frame(z = z_vec)
      ttl <- sprintf("Custom phenotype signature (%s gene%s)",
                     .fmt_int(length(z_vec)),
                     if (length(z_vec) == 1L) "" else "s")
      p <- ggplot(df, aes(x = z)) +
        .geom_rounded_histogram(bins = 50, fill = "#2A9D8F", color = "white") +
        labs(x = "Reference z-score", y = "Count", title = ttl) +
        theme_minimal(base_size = base_size)
      if (is.finite(cut_val) && cut_val > 0) {
        p <- p +
          ggplot2::geom_vline(
            xintercept = c(-cut_val, cut_val),
            linetype = "dashed", colour = "#264653", linewidth = 0.7
          ) +
          ggplot2::annotate(
            "text", x = -cut_val, y = Inf, label = cut_label,
            hjust = 1.05, vjust = 1.3, size = 3.2, colour = "#264653"
          ) +
          ggplot2::annotate(
            "text", x = cut_val, y = Inf, label = cut_label,
            hjust = -0.05, vjust = 1.3, size = 3.2, colour = "#264653"
          )
      }
      panel_objects$reference_signature_plot <- p
      return(p)
    }

    # Built-in: gate on `cancer_type` actually being a valid cancer
    # type for the currently selected reference. When the user
    # switches reference (e.g. PRECOG -> TCGA) Shiny fires
    # `input$reference_choice` immediately but `input$cancer_type`
    # only updates asynchronously via `updateSelectInput()`. Without
    # this guard, the renderer transiently sees a stale
    # (reference, cancer_type) pair like ("tcga", "Adrenocortical")
    # and `get_top_prognostic_genes()` throws.
    valid_cts <- get_cancer_types(state$reference)
    req(input$cancer_type, input$cancer_type %in% valid_cts)
    tg <- PhenoMapR::get_top_prognostic_genes(
      reference = state$reference, cancer_type = input$cancer_type,
      # Inf -> head(scores, Inf) returns the full sorted table.
      n = Inf, direction = "both"
    )
    z_vec <- as.numeric(tg$z_score)
    z_vec <- z_vec[is.finite(z_vec)]
    n_genes <- length(z_vec)
    if (!n_genes) {
      panel_objects$reference_signature_plot <- NULL
      return(NULL)
    }
    df <- data.frame(z = z_vec)
    df$side <- factor(
      ifelse(df$z < 0, "Negative tail (z < 0)", "Positive tail (z > 0)"),
      levels = c("Negative tail (z < 0)", "Positive tail (z > 0)")
    )
    ttl <- sprintf("Phenotype signature distribution (%s \u00B7 %s, %s genes)",
                   state$reference, input$cancer_type %||% "",
                   .fmt_int(n_genes))
    # facet_wrap(scales = "free_x") gives each tail its own
    # x-scale, producing the "broken axis" layout the user asked for
    # (built-in references ship pre-filtered: |z| >= 2 for PRECOG /
    # TCGA / Pediatric PRECOG, |z| >= 1 for ICI PRECOG, so the
    # interior is empty in either case).
    p <- ggplot(df, aes(x = z)) +
      .geom_rounded_histogram(bins = 30, fill = "#2A9D8F", color = "white") +
      ggplot2::facet_wrap(~ side, scales = "free_x", nrow = 1) +
      labs(x = "Reference z-score", y = "Count", title = ttl) +
      theme_minimal(base_size = base_size) +
      theme(strip.background = element_rect(fill = "#F4F1DE",
                                             colour = NA),
            strip.text = element_text(face = "bold"),
            panel.spacing.x = grid::unit(1.2, "lines"))
    # Live |z| cutoff vlines: -cut on the negative facet, +cut on
    # the positive facet. We pass a per-facet data.frame so each
    # facet only draws its own line (and the label sits inside the
    # right facet, not duplicated across both).
    if (is.finite(cut_val) && cut_val > 0) {
      vline_df <- data.frame(
        side = factor(
          c("Negative tail (z < 0)", "Positive tail (z > 0)"),
          levels = levels(df$side)
        ),
        x = c(-cut_val, cut_val),
        hjust = c(1.05, -0.05)
      )
      p <- p +
        ggplot2::geom_vline(
          data = vline_df,
          ggplot2::aes(xintercept = x),
          linetype = "dashed", colour = "#264653", linewidth = 0.7,
          inherit.aes = FALSE
        ) +
        ggplot2::geom_text(
          data = vline_df,
          ggplot2::aes(x = x, y = Inf, label = cut_label, hjust = hjust),
          vjust = 1.3, size = 3.2, colour = "#264653",
          inherit.aes = FALSE
        )
    }
    panel_objects$reference_signature_plot <- p
    p
  })
  phenomapr_register_plot_download(output, "reference_signature_plot",
    function() isolate(panel_objects$reference_signature_plot),
    width = 9, height = 5)

  # ------------------------------------------------------------------------
  # 3. Score
  # ------------------------------------------------------------------------
  observeEvent(input$run_score, {
    req(state$expression, state$reference)
    ref_label <- if (is.character(state$reference)) {
      sprintf("%s | %s", toupper(state$reference),
              input$cancer_type %||% "")
    } else "custom signature"
    phenomapr_busy_show("Computing PhenoMap scores...", ref_label)
    # CRITICAL: defer the heavy PhenoMap() computation into a
    # later::later() callback so the busy_show custom message can
    # actually reach the browser before libuv is blocked. See the
    # busy-overlay file-comment in helpers.R for the full
    # rationale.
    #
    # Direction radio was removed from the sidebar; PhenoMap() now
    # always uses the built-in signature's native direction
    # (reference_sign = 1L: higher score = worse prognosis), which
    # is the convention used by PRECOG / TCGA / Pediatric PRECOG /
    # ICI PRECOG signatures.
    expr_snapshot <- state$expression
    if (isTRUE(input$pseudobulk) && isTRUE(score_allow_pseudobulk()) &&
        identical(state$expr_summary$kind %||% "", "matrix") &&
        !is.null(state$metadata)) {
      cid <- input$meta_cell_id_col %||% ".cell_id"
      if (identical(cid, "(none)")) cid <- ".cell_id"
      expr_snapshot <- attach_matrix_coldata(
        expr_snapshot, state$metadata, cell_id_col = cid
      )
    }
    ref_snapshot  <- state$reference
    score_args <- list(
      cancer_type    = input$cancer_type,
      z_score_cutoff = input$z_score_cutoff,
      pseudobulk     = isTRUE(input$pseudobulk) && isTRUE(score_allow_pseudobulk()),
      group_by       = input$pseudobulk_group_by,
      assay          = input$score_assay,
      slot           = input$score_slot
    )
    sess <- session
    later::later(function() {
      scores <- tryCatch(
        run_phenomap_with_progress(
          expression = expr_snapshot,
          reference = ref_snapshot,
          cancer_type = score_args$cancer_type,
          z_score_cutoff = score_args$z_score_cutoff,
          pseudobulk = score_args$pseudobulk,
          group_by = score_args$group_by,
          assay = score_args$assay,
          slot = score_args$slot,
          reference_sign = 1L
        ),
        error = function(e) {
          showNotification(paste0("PhenoMap failed: ", conditionMessage(e)),
                           type = "error", duration = 12, session = sess); NULL
        }
      )
      # Guard against legacy empty score frames (character-only) so the
      # phenotype-group observer cannot surface the cryptic
      # "No numeric score columns found in 'scores'" message.
      if (!is.null(scores)) {
        n_num <- sum(vapply(scores, is.numeric, logical(1)))
        if (n_num < 1L) {
          showNotification(
            paste0(
              "PhenoMap returned no numeric score columns. ",
              "For ICI PRECOG, set Signature |z| cutoff to 1 ",
              "(many ICI cohorts have no genes with |z| > 2)."
            ),
            type = "error", duration = 12, session = sess
          )
          scores <- NULL
        }
      }
      # Hide the popup BEFORE assigning to state$scores so the hide
      # custom message lands in its own (tiny) flush ahead of the
      # ~half-dozen renderPlot/renderUI outputs that re-render off
      # state$scores (histogram, rank plot, box-by-source, group
      # enrichment, score table, status, etc.).
      phenomapr_busy_hide(session = sess)
      if (is.null(scores)) return()
      later::later(function() {
        tryCatch({
          state$scores <- scores
          state$groups <- NULL  # invalidate downstream
          state$markers <- NULL
          shiny::showNotification(
            sprintf("Scored %s (%.1fs).",
                    .fmt_n_units(nrow(scores), "row"),
                    attr(scores, "elapsed_s") %||% NA_real_),
            type = "message", duration = 5,
            session = sess
          )
        }, error = function(e) {
          shiny::showNotification(
            paste0("Internal error while applying scores: ",
                   conditionMessage(e)),
            type = "error", duration = 10, session = sess
          )
        })
      }, delay = 0)
    }, delay = 0)
  })

  # Score distribution histogram. The "Score scale" radio in the card
  # header lets the user flip between the raw weighted-sum PhenoMapR score
  # and its standardized (mean 0, sd 1) counterpart. Both modes hand the
  # appropriate vector to `plot_score_distribution()` along with a matching
  # plot title so the y-axis context is always obvious.
  output$score_dist_plot <- renderPlot({
    req(state$scores)
    s <- state$scores
    cn <- colnames(s)[1L]
    scale_choice <- input$dist_score_scale %||% "raw"
    if (identical(scale_choice, "scaled")) {
      v <- as.numeric(scale(s[[cn]]))
      df <- data.frame(score = v)
      ttl <- sprintf("Z-score distribution (%s)", cn)
    } else {
      v <- as.numeric(s[[cn]])
      df <- data.frame(score = v)
      ttl <- sprintf("PhenoMapR score distribution (%s)", cn)
    }
    p <- PhenoMapR::plot_score_distribution(df, score_column = "score",
                                            main = ttl,
                                            base_size = .theme_base_size())

    # Overlay the current Phenotype-groups tail thresholds as dashed
    # vertical lines on the histogram so users can see, in real time,
    # exactly where the slider is splitting the distribution. The
    # cutoffs are recomputed from the same vector that drives the
    # histogram (`v`) so they line up with the bars regardless of the
    # PhenoMapR / Z-score scale toggle above.
    pct <- input$percentile %||% 0.05
    pct <- max(min(pct, 0.49), 0.001)
    vf  <- v[is.finite(v)]
    if (length(vf) >= 2L) {
      q_lo <- stats::quantile(vf, pct,       na.rm = TRUE, names = FALSE)
      q_hi <- stats::quantile(vf, 1 - pct,   na.rm = TRUE, names = FALSE)
      tail_label <- sprintf("Tail = %.0f%%", pct * 100)
      p <- p +
        ggplot2::geom_vline(
          xintercept = c(q_lo, q_hi),
          linetype   = "dashed",
          colour     = "#264653",
          linewidth  = 0.7
        ) +
        ggplot2::annotate(
          "text", x = q_lo, y = Inf, label = tail_label,
          hjust = 1.05, vjust = 1.3, size = 3.2, colour = "#264653"
        ) +
        ggplot2::annotate(
          "text", x = q_hi, y = Inf, label = tail_label,
          hjust = -0.05, vjust = 1.3, size = 3.2, colour = "#264653"
        )
    }
    panel_objects$score_dist_plot <- p
    p
  })
  phenomapr_register_plot_download(output, "score_dist_plot",
    function() isolate(panel_objects$score_dist_plot),
    width = 8, height = 5)

  # Score table data. Used to be a DT::datatable panel that lived in
  # its own card right below the "Score by cell type and source"
  # boxplot, but the visual panel was redundant with the sidebar's
  # "Download score table (TSV)" affordance and just pushed the more
  # useful Phenotype-groups section below the fold. We removed the
  # panel and reduced this to a plain reactive: the boxplot card
  # header now hosts a "Download plot data (TSV)" button whose
  # downloadHandler reads this reactive directly, so users still get
  # the per-cell score table without needing it rendered as a table.
  #
  # Schema:
  #   cell_id                                    -- cell identifier
  #   <score_col>                                -- raw PhenoMapR score
  #   scaled_<score_col>                         -- z-scaled score
  #   phenotype_group_<score_col>                -- "Most Adverse" /
  #     "Most Favorable" / "Other" (only present when scores have
  #     been thresholded into groups, i.e. state$groups is set).
  #
  # The phenotype-group merge means the "Download plot data" TSV is
  # self-contained: the user can recolour / refilter the boxplot
  # offline using the canonical group labels without also having to
  # download the sidebar "Download labels (TSV)" file separately.
  # We keep the CANONICAL strings ("Most Adverse" / "Most Favorable")
  # in the export rather than the displayed "Most Phenotype +/-"
  # labels so the file is interoperable with PhenoMapR R functions
  # (find_phenotype_markers, plot_phenotype_markers, the
  # prognostic_analysis tests + vignettes) that pattern-match those
  # literals.
  score_table_data <- reactive({
    req(state$scores)
    s <- state$scores
    s$cell_id <- rownames(s)
    score_cols <- setdiff(colnames(s), "cell_id")
    numeric_cols <- score_cols[vapply(s[score_cols], is.numeric, logical(1))]
    for (col in numeric_cols) {
      v <- as.numeric(scale(s[[col]]))
      s[[paste0("scaled_", col)]] <- v
    }

    # Optionally splice phenotype_group_<score_col> columns from
    # state$groups. state$groups carries cell IDs in BOTH `cell_id`
    # column AND rownames; we join on rownames so the alignment
    # survives even if some scoring runs strip the cell_id column
    # downstream.
    grp_cols <- character(0)
    if (!is.null(state$groups)) {
      grp_cols <- grep("^phenotype_group_", colnames(state$groups),
                       value = TRUE)
      if (length(grp_cols)) {
        join_idx <- match(rownames(s), rownames(state$groups))
        for (gc in grp_cols) {
          s[[gc]] <- state$groups[[gc]][join_idx]
        }
      }
    }

    ordered_cols <- c(
      "cell_id",
      unlist(lapply(score_cols, function(col) {
        if (col %in% numeric_cols) c(col, paste0("scaled_", col)) else col
      }), use.names = FALSE),
      grp_cols
    )
    s[, ordered_cols, drop = FALSE]
  })
  phenomapr_register_table_download(output, "score_table",
                                    function() score_table_data())

  output$download_scores <- downloadHandler(
    filename = function() {
      sprintf("phenomapr_scores_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S"))
    },
    content = function(file) {
      req(state$scores)
      df <- state$scores; df$cell_id <- rownames(df)
      df <- df[, c("cell_id", setdiff(colnames(df), "cell_id"))]
      utils::write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  # ------------------------------------------------------------------------
  # Score rank plot
  # ------------------------------------------------------------------------
  output$score_rank_plot <- renderPlot({
    df <- cell_table()
    req(df, "score" %in% colnames(df))
    d <- df[is.finite(df$score), ]
    if (!nrow(d)) {
      panel_objects$score_rank_plot <- NULL
      return(NULL)
    }
    d <- d[order(d$score), ]
    d$rank <- seq_len(nrow(d))
    color_by <- input$rank_color_by %||% "score"
    # `legend_pt_override` enlarges the small per-cell dots (0.7pt) in the
    # legend to a comfortable readable size (4pt) without bloating the
    # plotted points themselves. Cell type and Source legends both get
    # this treatment so users can identify colors at a glance.
    legend_pt_override <- ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 4, alpha = 1)
      )
    )
    # Coerce categorical aesthetics before ggplot() captures `d`. Mutating
    # after ggplot(d) does not update plot data, so integer source IDs
    # would still hit scale_color_phenomapr_d as continuous values.
    if (color_by == "cell_type" && "cell_type" %in% colnames(d)) {
      d$cell_type <- as.character(d$cell_type)
    }
    if (color_by == "source" && "source" %in% colnames(d)) {
      d$source <- factor(as.character(d$source),
                         levels = unique(as.character(d$source)))
    }
    base <- ggplot(d, aes(x = rank, y = score)) +
      labs(x = "Rank by PhenoMapR score", y = "PhenoMapR score") +
      theme_minimal(base_size = .theme_base_size())
    if (color_by == "cell_type" && "cell_type" %in% colnames(d)) {
      pal <- tryCatch(
        PhenoMapR::get_celltype_palette(as.character(unique(d$cell_type))),
        error = function(e) NULL
      )
      p <- base +
        geom_point(aes(color = cell_type), size = 0.7, alpha = 0.85) +
        labs(color = "Cell type") +
        legend_pt_override
      if (!is.null(pal)) p <- p + scale_color_manual(values = pal)
      panel_objects$score_rank_plot <- p
      return(p)
    }
    if (color_by == "source" && "source" %in% colnames(d)) {
      # Reuse the PhenoMapR brand palette so the "Source" coloring in this
      # plot matches the "Cells per source / group" and "Cell type ×
      # source composition" plots on the Data tab. Sort levels by their
      # appearance order in the input metadata so the color mapping is
      # stable across refreshes.
      p <- base +
        geom_point(aes(color = source), size = 0.7, alpha = 0.85) +
        scale_color_phenomapr_d() +
        labs(color = "Source") +
        legend_pt_override
      panel_objects$score_rank_plot <- p
      return(p)
    }
    p <- base +
      geom_point(aes(color = score), size = 0.7, alpha = 0.85) +
      scale_color_gradient2(
        low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
        midpoint = 0,
        name = "Score"
      )
    panel_objects$score_rank_plot <- p
    p
  })
  phenomapr_register_plot_download(output, "score_rank_plot",
    function() isolate(panel_objects$score_rank_plot),
    width = 9, height = 5)

  # ------------------------------------------------------------------------
  # Per-cell-type × source boxplot with Wilcoxon brackets
  # ------------------------------------------------------------------------
  output$score_box_source_plot <- renderPlot({
    df <- cell_table()
    req(df,
        "score" %in% colnames(df),
        "cell_type" %in% colnames(df))
    has_source <- "source" %in% colnames(df) &&
      length(unique(stats::na.omit(df$source))) >= 2L

    keep <- is.finite(df$score) & !is.na(df$cell_type)
    if (has_source) keep <- keep & !is.na(df$source)
    d <- df[keep, , drop = FALSE]
    if (!nrow(d)) {
      panel_objects$score_box_source_plot <- NULL
      return(NULL)
    }

    dl <- data.frame(
      Score     = d$score,
      Cell_type = as.character(d$cell_type),
      stringsAsFactors = FALSE
    )
    if (has_source) {
      dl$Source <- as.character(d$source)
    }

    # Always order cell types from lowest to highest median PhenoMapR score
    # (computed on the unscaled score; identical to the centred / scaled
    # ranking, since `scale()` is a monotone transform across all cells).
    med <- stats::aggregate(Score ~ Cell_type, data = dl,
                            FUN = function(x) stats::median(x, na.rm = TRUE))
    med <- med[order(med$Score), , drop = FALSE]
    cell_levels <- as.character(med$Cell_type)

    # Choose the appropriate significance test for the bracket annotation:
    #   - No source mapped (only ordered cell types): single global
    #     one-way ANOVA across cell types (F-test). Returns a single
    #     bracket spanning the full x-axis carrying the ANOVA p-value
    #     -- "is there ANY difference between cell types?".
    #   - Source mapped with exactly 2 levels: per-cell-type Wilcoxon
    #     contrast between the two sources (the classic 2-sample test).
    #   - Source mapped with 3+ levels: per-cell-type one-way ANOVA
    #     across sources (Score ~ Source) within each cell type, since
    #     pairwise Wilcoxon is no longer the right comparison for >2
    #     groups.
    # All three helpers honour the median-based `cell_levels` ordering
    # so the bracket xmin/xmax positions line up with the plot axis.
    n_src_levels <- if (has_source) {
      length(unique(stats::na.omit(dl$Source)))
    } else 0L
    if (!has_source) {
      pval_df <- celltype_anova_pvalue(
        dl, "Score", "Cell_type", cell_levels = cell_levels
      )
    } else if (n_src_levels == 2L) {
      pval_df <- celltype_source_pvalues(
        dl, "Score", "Cell_type", "Source",
        cell_levels = cell_levels, significant_only = TRUE
      )
    } else {
      pval_df <- celltype_source_anova_pvalues(
        dl, "Score", "Cell_type", "Source",
        cell_levels = cell_levels, significant_only = TRUE
      )
    }

    dl$Cell_type <- factor(dl$Cell_type, levels = cell_levels)

    pal_cells <- tryCatch(
      PhenoMapR::get_celltype_palette(cell_levels),
      error = function(e) NULL
    )

    # Cell-type legend columns -- defensively coerce to a positive
    # integer so a NULL / blank / non-finite numericInput doesn't break
    # the plot. Bounded to [1, 8] to match the input's max.
    legend_ncol_raw <- suppressWarnings(
      as.integer(input$score_box_celltype_legend_ncol)
    )
    legend_ncol <- if (length(legend_ncol_raw) == 1L &&
                      is.finite(legend_ncol_raw) &&
                      legend_ncol_raw >= 1L) {
      min(8L, legend_ncol_raw)
    } else {
      1L
    }

    if (has_source) {
      # Source levels are colored with the PhenoMapR brand palette by
      # default so this plot lines up 1:1 with the other source-keyed
      # plots ("Cells per source / group", "Cell type x source
      # composition", score rank plot, UMAP). The "Greyscale" radio
      # toggle swaps in a black-to-light-grey ramp instead -- useful
      # for publications that print in black/white, or when the user
      # wants all the visual emphasis on the cell-type fill colour
      # and only a secondary monochrome cue on Source. The level
      # order is frozen via factor() so palette -> level mapping is
      # stable in either branch.
      src_levels <- sort(unique(dl$Source))
      dl$Source <- factor(dl$Source, levels = src_levels)
      use_greyscale <- identical(input$score_box_source_palette, "greyscale")
      src_pal <- if (use_greyscale) {
        # gray.colors() with a slight off-black start (0.15) and an
        # off-white end (0.75) keeps adjacent shades visually
        # distinguishable on box outlines without flat-out white,
        # which would disappear into the panel background.
        setNames(grDevices::gray.colors(
          length(src_levels), start = 0.15, end = 0.75
        ), src_levels)
      } else {
        setNames(pm_brand_palette(length(src_levels)), src_levels)
      }
      p <- ggplot(dl, aes(x = Cell_type, y = scale(Score),
                          fill = Cell_type, color = Source)) +
        .geom_rounded_boxplot(
          outlier.alpha    = 0.5,
          median.linewidth = 0.5,
          outlier.fill     = NULL,
          outlier.color    = NULL,
          outlier.shape    = 21
        ) +
        scale_color_manual(values = src_pal, name = "Source") +
        guides(
          fill  = guide_legend(order = 1, ncol = legend_ncol, byrow = TRUE),
          color = guide_legend(order = 2)
        )
      title <- "PhenoMapR Score distribution by Source and Cell Type"
    } else {
      p <- ggplot(dl, aes(x = Cell_type, y = scale(Score), fill = Cell_type)) +
        .geom_rounded_boxplot(
          outlier.alpha    = 0.5,
          median.linewidth = 0.5,
          outlier.shape    = 21
        ) +
        guides(fill = guide_legend(ncol = legend_ncol, byrow = TRUE))
      title <- "PhenoMapR Score distribution by Cell Type"
    }

    # Subtitle. When per-cell-type bracket annotations are drawn (the
    # Wilcoxon-2-source path or the per-CT ANOVA path), prepend a
    # significance-key line so users know what * / ** / *** mean
    # without scanning a separate legend. The bare-ANOVA case
    # (annotation render_kind, no stars on the plot) doesn't get the
    # key because the label there literally reads "ANOVA (p = X.XXX)".
    render_kind  <- if (!is.null(pval_df))
      (attr(pval_df, "render_kind") %||% "bracket") else "bracket"
    has_brackets <- !is.null(pval_df) && nrow(pval_df) > 0L &&
      !identical(render_kind, "annotation")
    plot_subtitle <- if (has_brackets) {
      "Significance: *** p < 0.001, ** p < 0.01, * p < 0.05"
    } else NULL

    p <- p +
      labs(x = NULL, y = "PhenoMapR score (scaled)",
           title = title, subtitle = plot_subtitle) +
      theme_minimal(base_size = .theme_base_size()) +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        plot.title      = element_text(hjust = 0.5),
        plot.subtitle   = element_text(
          hjust = 0.5, size = .theme_base_size() * 0.85,
          colour = "#3a3a3a", margin = ggplot2::margin(t = 1, b = 6)
        )
      )
    if (!is.null(pal_cells)) {
      p <- p + scale_fill_manual(values = pal_cells, name = "Cell type")
    } else {
      p <- p + labs(fill = "Cell type")
    }

    if (!is.null(pval_df) && nrow(pval_df) > 0) {
      # Two annotation modes:
      #   - "annotation": single centred text label, no bracket or
      #     connector line (currently the global-ANOVA case from
      #     celltype_anova_pvalue()).
      #   - default      : per-cell-type significance brackets drawn
      #     manually via geom_segment + geom_text. We do NOT use
      #     ggsignif::geom_signif() here because in some
      #     ggsignif/ggplot2 version combos the manual-mode renderer
      #     drops short labels ("*", "**") while keeping long ones
      #     ("***"), producing the empty-bracket look the user
      #     reported. Drawing brackets ourselves with explicit
      #     geom_segment + geom_text guarantees every bracket carries
      #     its label regardless of length.
      if (identical(render_kind, "annotation")) {
        # Draw the ANOVA p-value as a single centred text annotation
        # above the boxes. No bracket, no significance stars -- the
        # plot stays uncluttered and the bare "ANOVA (p = ...)"
        # label tells the story.
        x_mid <- pval_df$x_mid[1L] %||%
          ((pval_df$xmin[1L] + pval_df$xmax[1L]) / 2)
        p <- p +
          ggplot2::annotate(
            "text",
            x = x_mid,
            y = pval_df$y_pos[1L],
            label = pval_df$label[1L],
            size = 4.2,
            fontface = "plain",
            colour = "black"
          )
      } else {
        # Manual bracket drawing. Each row of pval_df contributes
        # three primitives:
        #   1. A horizontal segment from xmin to xmax at y_pos
        #      (the bracket spine).
        #   2. Two short vertical segments at xmin / xmax dropping
        #      down a small amount (the bracket "tips").
        #   3. A geom_text label centred above the spine, raised by
        #      a constant offset that scales with the y range so the
        #      label never overlaps the spine.
        # The spine y is taken straight from helpers.R's pval_df
        # (which uses `y_top * 1.05`); for the label we add a small
        # constant in scaled-score units. This keeps every label
        # visible and at the same height regardless of star length.
        y_top <- attr(pval_df, "y_top") %||% max(pval_df$y_pos, na.rm = TRUE)
        if (!is.finite(y_top)) y_top <- 0
        # Tip length in y units. Using a small fraction of y_top
        # keeps the bracket tips proportionate to the data range.
        tip_h <- max(0.08, abs(y_top) * 0.05)
        # Label offset above the spine. A bit larger than tip_h so
        # the text never collides with the bracket spine.
        lab_offset <- max(0.18, abs(y_top) * 0.10)

        seg_df <- data.frame(
          x    = pval_df$xmin,
          xend = pval_df$xmax,
          y    = pval_df$y_pos,
          yend = pval_df$y_pos,
          stringsAsFactors = FALSE
        )
        tick_df <- data.frame(
          x    = c(pval_df$xmin, pval_df$xmax),
          xend = c(pval_df$xmin, pval_df$xmax),
          y    = c(pval_df$y_pos, pval_df$y_pos),
          yend = c(pval_df$y_pos - tip_h, pval_df$y_pos - tip_h),
          stringsAsFactors = FALSE
        )
        text_df <- data.frame(
          x     = (pval_df$xmin + pval_df$xmax) / 2,
          y     = pval_df$y_pos + lab_offset,
          label = as.character(pval_df$label),
          stringsAsFactors = FALSE
        )

        p <- p +
          geom_segment(
            data = seg_df,
            aes(x = x, xend = xend, y = y, yend = yend),
            inherit.aes = FALSE, color = "black", linewidth = 0.5
          ) +
          geom_segment(
            data = tick_df,
            aes(x = x, xend = xend, y = y, yend = yend),
            inherit.aes = FALSE, color = "black", linewidth = 0.5
          ) +
          geom_text(
            data = text_df,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 5, fontface = "bold", colour = "black",
            vjust = 0
          )

        # Expand the y axis upward so the labels (especially the
        # bold "*" / "**" / "***" text at y_pos + lab_offset) are
        # never clipped by the panel boundary. Without this, ggplot's
        # auto-scaling can crop the top row of stars on smaller
        # devices, which is exactly the symptom the user reported.
        y_max <- max(text_df$y, na.rm = TRUE) + lab_offset * 0.6
        p <- p + ggplot2::expand_limits(y = y_max)
      }
    }

    panel_objects$score_box_source_plot <- p
    p
  })
  phenomapr_register_plot_download(output, "score_box_source_plot",
    function() isolate(panel_objects$score_box_source_plot),
    width = 11, height = 6)

  # ------------------------------------------------------------------------
  # UMAP / embedding tab
  # ------------------------------------------------------------------------
  # Custom embedding upload (overrides object-resident reductions).
  uploaded_embedding <- reactiveVal(NULL)

  observeEvent(umap_upload_pick(), {
    pick <- umap_upload_pick()
    if (is.null(pick)) {
      uploaded_embedding(NULL)
      return()
    }
    ext <- tolower(tools::file_ext(pick$name))
    obj <- tryCatch({
      if (ext == "rds") {
        readRDS(pick$datapath)
      } else if (ext == "tsv" || ext == "txt") {
        utils::read.delim(pick$datapath, sep = "\t",
                          stringsAsFactors = FALSE, check.names = FALSE)
      } else if (ext == "csv") {
        utils::read.csv(pick$datapath,
                        stringsAsFactors = FALSE, check.names = FALSE)
      } else {
        stop("Unsupported embedding file type: ", ext)
      }
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 8); NULL
    })
    if (is.null(obj)) return()
    df <- if (is.matrix(obj) || is.data.frame(obj)) as.data.frame(obj) else NULL
    if (is.null(df) || ncol(df) < 2L) {
      showNotification("Embedding file must have at least 2 columns.",
                       type = "error", duration = 6); return()
    }
    # Pick the cell-id column heuristically.
    id_col <- intersect(colnames(df), c("cell_id", "Cell", "barcode", "CellID"))[1L]
    if (is.na(id_col)) {
      id_col <- colnames(df)[1L]
    }
    # Numeric coord columns
    num_cols <- vapply(df, is.numeric, logical(1L))
    coord_cols <- colnames(df)[num_cols]
    if (length(coord_cols) < 2L) {
      showNotification("Embedding file needs two numeric coordinate columns.",
                       type = "error", duration = 6); return()
    }
    emb_df <- data.frame(
      cell_id = as.character(df[[id_col]]),
      dim1 = as.numeric(df[[coord_cols[1L]]]),
      dim2 = as.numeric(df[[coord_cols[2L]]]),
      dim1_name = coord_cols[1L],
      dim2_name = coord_cols[2L],
      stringsAsFactors = FALSE
    )
    uploaded_embedding(emb_df)
    showNotification(
      sprintf("Custom embedding loaded (%s cells).", .fmt_int(nrow(emb_df))),
      type = "message", duration = 4
    )
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  # Populate the reduction selector with:
  #  1. embeddings stored on the loaded object (Seurat / SCE / AnnData),
  #  2. embeddings auto-detected as coordinate column pairs on the metadata
  #     table when the expression input is a plain matrix / data.frame
  #     (so users don't have to upload the same TSV twice), and
  #  3. an explicitly uploaded embedding file if one was supplied.
  detected_meta_embeddings <- reactive({
    obj_avail <- list_available_embeddings(state$expression)
    # Only mine the metadata for coordinate pairs when the loaded object
    # has none of its own — otherwise the on-object reductions take priority.
    if (length(obj_avail) > 0L) return(list())
    detect_metadata_embeddings(state$metadata)
  })

  observe({
    obj_avail  <- list_available_embeddings(state$expression)
    meta_embs  <- detected_meta_embeddings()
    has_upload <- !is.null(uploaded_embedding())

    choices <- character(0)
    if (has_upload) {
      choices <- c(choices, setNames("(uploaded)", "(uploaded file)"))
    }
    for (a in obj_avail) {
      choices <- c(choices, setNames(a, a))
    }
    for (s in names(meta_embs)) {
      choices <- c(
        choices,
        setNames(paste0("__meta__:", s), paste0(s, "  (from metadata)"))
      )
    }
    if (!length(choices)) {
      updateSelectInput(session, "umap_reduction",
                        choices = c("(none available)" = "(none available)"),
                        selected = "(none available)")
    } else {
      updateSelectInput(session, "umap_reduction",
                        choices = choices, selected = unname(choices)[1L])
    }
  })

  # "Tissue section" / "Core" / "Library" switcher for multi-sample
  # spatial objects. Renders ONLY when (a) current embedding is spatial
  # and (b) the embedding df spans >1 section.
  #
  # Two distinct UX modes share the same selectInput value:
  #   * "__all__" -> every section in one panel. Slide-mm / global coords
  #     overlay naturally; per-FOV pixel coords switch to the object's
  #     global spatial frame (see apply_global_spatial_coords_for_combined()).
  #   * any specific section name -> zoom to one tissue. The adjacent
  #     adjacent Prev / Next buttons walk through the section list one
  #     step at a time so users don't have to scroll a long dropdown.
  #     This is the "click through" mode -- intentionally separate from
  #     the facet mode so users get the panoramic view AND the focused
  #     view without one mode hiding the other.
  output$spatial_sample_selector <- renderUI({
    emb <- current_embedding()
    if (is.null(emb)) return(NULL)
    if (!isTRUE(any(emb$is_spatial))) return(NULL)
    sections <- unique(as.character(emb$sample))
    sections <- sections[!is.na(sections) & nzchar(sections)]
    if (length(sections) < 2L) return(NULL)
    # Stable, alphabetical ordering so the dropdown doesn't reshuffle as
    # data lands. "All sections (facet view)" always stays at the top.
    sections <- sort(sections)
    choices <- c(
      "All sections (combined)" = "__all__",
      setNames(sections, sections)
    )
    tagList(
      selectInput(
        "spatial_sample",
        label = tagList(icon("layer-group"), " Tissue section / core"),
        choices = choices,
        selected = "__all__"
      ),
      # Prev / Next stepper -- meaningful only when the user has chosen
      # an individual section. Hidden under "All sections (facet view)"
      # via conditionalPanel so the arrows aren't misleading in that
      # mode (there's nothing to step to when every section is shown).
      conditionalPanel(
        condition = "input.spatial_sample != '__all__'",
        tags$div(
          class = "spatial-sample-stepper",
          actionButton(
            "spatial_sample_prev",
            label = NULL,
            icon = icon("chevron-left"),
            title = "Previous section",
            class = "btn btn-sm btn-outline-secondary spatial-sample-step-btn"
          ),
          tags$span(
            class = "spatial-sample-counter",
            uiOutput("spatial_sample_counter", inline = TRUE)
          ),
          actionButton(
            "spatial_sample_next",
            label = NULL,
            icon = icon("chevron-right"),
            title = "Next section",
            class = "btn btn-sm btn-outline-secondary spatial-sample-step-btn"
          )
        )
      )
    )
  })

  # "42 / 60" counter shown between the Prev / Next arrows so users
  # know how far through the section list they are. Returns NULL (which
  # renders nothing) when we're in facet view or on a single-section
  # dataset.
  output$spatial_sample_counter <- renderUI({
    emb <- current_embedding()
    if (is.null(emb) || !isTRUE(any(emb$is_spatial))) return(NULL)
    sections <- sort(unique(as.character(emb$sample)))
    sections <- sections[!is.na(sections) & nzchar(sections)]
    if (length(sections) < 2L) return(NULL)
    cur <- input$spatial_sample %||% sections[1L]
    if (identical(cur, "__all__")) return(NULL)
    idx <- match(cur, sections)
    if (is.na(idx)) return(NULL)
    sprintf("%d / %d", idx, length(sections))
  })

  # Helper: re-derive the ordered section list from the current
  # embedding. Called by both Prev and Next observers so we don't
  # duplicate the filter / sort logic. Returns character(0) when the
  # current embedding has fewer than two sections.
  .current_spatial_sections <- function() {
    emb <- current_embedding()
    if (is.null(emb) || !isTRUE(any(emb$is_spatial))) return(character(0))
    sections <- sort(unique(as.character(emb$sample)))
    sections <- sections[!is.na(sections) & nzchar(sections)]
    if (length(sections) < 2L) return(character(0))
    sections
  }

  # Step backward through the section list. From "All sections" we
  # wrap to the last section -- so a user who lands in facet view and
  # immediately clicks Prev gets shown the final tissue (the visual
  # opposite of clicking Next, which would show the first).
  observeEvent(input$spatial_sample_prev, {
    sections <- .current_spatial_sections()
    if (!length(sections)) return()
    cur <- input$spatial_sample %||% sections[1L]
    new_val <- if (identical(cur, "__all__")) {
      sections[length(sections)]
    } else {
      idx <- match(cur, sections)
      new_idx <- if (is.na(idx) || idx <= 1L) length(sections) else idx - 1L
      sections[new_idx]
    }
    updateSelectInput(session, "spatial_sample", selected = new_val)
  })

  # Step forward through the section list. Wraps from the last section
  # back to the first. From "All sections", clicking Next exits facet
  # view by selecting the first individual section.
  observeEvent(input$spatial_sample_next, {
    sections <- .current_spatial_sections()
    if (!length(sections)) return()
    cur <- input$spatial_sample %||% sections[1L]
    new_val <- if (identical(cur, "__all__")) {
      sections[1L]
    } else {
      idx <- match(cur, sections)
      new_idx <- if (is.na(idx) || idx >= length(sections)) 1L else idx + 1L
      sections[new_idx]
    }
    updateSelectInput(session, "spatial_sample", selected = new_val)
  })

  output$umap_reduction_status <- renderUI({
    obj_avail <- list_available_embeddings(state$expression)
    meta_embs <- detected_meta_embeddings()
    if (length(obj_avail) > 0L) {
      tags$div(
        class = "umap-status umap-status-ok",
        tags$strong(sprintf("%d reduction%s on the loaded object.",
                            length(obj_avail),
                            if (length(obj_avail) == 1L) "" else "s"))
      )
    } else if (length(meta_embs) > 0L) {
      tags$div(
        class = "umap-status umap-status-ok",
        tags$strong(sprintf("Auto-detected %d embedding%s in the metadata table:",
                            length(meta_embs),
                            if (length(meta_embs) == 1L) "" else "s")),
        tags$br(),
        tags$small(paste(
          vapply(names(meta_embs), function(s) {
            sprintf("%s = %s / %s", s,
                    meta_embs[[s]]$dim1, meta_embs[[s]]$dim2)
          }, character(1L)),
          collapse = "; "
        ))
      )
    } else {
      tags$div(
        class = "umap-status umap-status-empty",
        tags$em("No embedding found on the loaded object or metadata table."),
        tags$br(),
        tags$small("Upload a separate embedding file below, or add ",
                   tags$code("UMAP_1 / UMAP_2"),
                   " (or similar) coordinate columns to your metadata.")
      )
    }
  })

  current_embedding <- reactive({
    sel <- input$umap_reduction
    if (is.null(sel) || !nzchar(sel) || sel == "(none available)") return(NULL)
    if (sel == "(uploaded)") return(uploaded_embedding())
    if (startsWith(sel, "__meta__:")) {
      nm <- sub("^__meta__:", "", sel)
      meta_embs <- detected_meta_embeddings()
      if (!nm %in% names(meta_embs)) return(NULL)
      p <- meta_embs[[nm]]
      return(extract_embedding_from_metadata(
        state$metadata,
        dim1_col = p$dim1,
        dim2_col = p$dim2,
        cell_id_col = input$meta_cell_id_col
      ))
    }
    extract_embedding(state$expression, sel)
  })

  spatial_polygon_masks_available <- reactive({
    obj <- state$expression
    if (is.null(obj) || !spatial_polygons_available(obj)) return(FALSE)
    sel <- input$umap_reduction %||% ""
    if (!identical(sel, "segmentation")) return(FALSE)
    emb <- current_embedding()
    !is.null(emb) && isTRUE(any(emb$is_spatial))
  })
  output$spatial_polygon_masks_available <- spatial_polygon_masks_available
  outputOptions(output, "spatial_polygon_masks_available", suspendWhenHidden = FALSE)

  output$umap_plot <- renderPlot({
    emb <- current_embedding()
    req(emb)
    color_by <- input$umap_color_by %||% "score"
    ct <- cell_table()
    df <- emb
    df <- dplyr::left_join(df, ct %||% data.frame(cell_id = character(0)), by = "cell_id")
    score_col <- active_score_column(state$scores)
    if (!is.null(state$scores) && nzchar(score_col) && score_col %in% colnames(state$scores)) {
      df$score <- as.numeric(state$scores[[score_col]][match(df$cell_id, rownames(state$scores))])
      grp_col <- resolve_phenotype_group_column(state$groups, score_col)
      if (!is.na(grp_col) && !is.null(state$groups)) {
        df$group <- as.character(state$groups[[grp_col]][match(df$cell_id, state$groups$cell_id)])
      }
    }
    pt_size <- input$umap_point_size %||% 0.8
    pt_alpha <- input$umap_point_alpha %||% 0.75
    mask_alpha <- input$spatial_mask_alpha %||% 0.85
    # Detect spatial-frame embeddings (set by `extract_embedding()` when
    # the user picks the synthetic "spatial" reduction). Spatial plots
    # need (a) equal aspect so tissue isn't squashed and (b) a reversed
    # y-axis since image-space coordinates have origin at top-left.
    is_spatial <- isTRUE(any(df$is_spatial))
    spatial_pick <- input$spatial_sample %||% "__all__"
    all_sections <- is_spatial && identical(spatial_pick, "__all__")

    # Multi-section spatial objects: when the user picks a single
    # section in the sidebar, restrict the plot to that section so
    # different tissues don't get drawn on top of each other in the
    # same coord_fixed() frame. "__all__" (or no selector at all)
    # keeps every spot -- ggplot's coord_fixed() puts all sections on
    # the same scale, which is what users want when they're scanning.
    if (is_spatial && "sample" %in% colnames(df)) {
      pick <- input$spatial_sample %||% "__all__"
      if (nzchar(pick) && !identical(pick, "__all__")) {
        df <- df[!is.na(df$sample) & df$sample == pick, , drop = FALSE]
        req(nrow(df) > 0L)
      }
    }

    use_masks <- is_spatial && isTRUE(spatial_polygon_masks_available()) &&
      identical(input$spatial_render_style %||% "points", "masks")
    poly_df <- NULL
    if (use_masks) {
      section_arg <- NULL
      if (is_spatial && "sample" %in% colnames(df) &&
          !identical(spatial_pick, "__all__")) {
        section_arg <- spatial_pick
      }
      poly_df <- extract_spatial_polygons(state$expression, section = section_arg)
      req(!is.null(poly_df) && nrow(poly_df) > 0L)
      poly_df <- poly_df[poly_df$cell_id %in% df$cell_id, , drop = FALSE]
      req(nrow(poly_df) > 0L)
    }

    if (all_sections && .spatial_coords_are_fov_local(df)) {
      df <- apply_global_spatial_coords_for_combined(state$expression, df)
    }

    # Discrete aesthetics must be character/factor BEFORE ggplot() captures
    # `df`. Mutating after ggplot(df) leaves integer FOV/core IDs continuous
    # in the plot data and breaks scale_*_phenomapr_d().
    if ("cell_type" %in% colnames(df)) {
      df$cell_type <- as.character(df$cell_type)
    }
    if ("source" %in% colnames(df)) {
      df$source <- factor(as.character(df$source),
                          levels = unique(as.character(df$source)))
    }
    if ("group" %in% colnames(df)) {
      df$group <- as.character(df$group)
    }

    base <- if (use_masks) {
      ggplot(poly_df, aes(x = x, y = y)) +
        labs(x = unique(df$dim1_name)[1L] %||% "x",
             y = unique(df$dim2_name)[1L] %||% "y") +
        theme_minimal(base_size = .theme_base_size()) +
        theme(panel.grid.minor = element_blank())
    } else {
      ggplot(df, aes(x = dim1, y = dim2)) +
        labs(x = unique(df$dim1_name)[1L] %||% "dim1",
             y = unique(df$dim2_name)[1L] %||% "dim2") +
        theme_minimal(base_size = .theme_base_size()) +
        theme(panel.grid.minor = element_blank())
    }
    if (is_spatial) {
      base <- base +
        coord_fixed() +
        scale_y_reverse() +
        labs(x = NULL, y = NULL) +
        theme(
          axis.text  = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()
        )
    }

    if (color_by == "score" && "score" %in% colnames(df) && any(is.finite(df$score))) {
      # Honor the sidebar "Score scale" toggle: when "Z-score" is
      # selected we standardize the raw PhenoMapR score (mean 0, sd 1)
      # before plotting so cells are colored by their unit-free z-value
      # instead of the raw weighted-sum magnitude. This mirrors the
      # behavior of the Score-tab distribution plot toggle.
      score_scale <- input$umap_score_scale %||% "raw"
      if (identical(score_scale, "scaled")) {
        finite_idx <- is.finite(df$score)
        z <- rep(NA_real_, nrow(df))
        if (any(finite_idx)) z[finite_idx] <- as.numeric(scale(df$score[finite_idx]))
        df$score_to_plot <- z
        legend_name <- "Z-score"
      } else {
        df$score_to_plot <- df$score
        legend_name <- "PhenoMapR\nscore"
      }
      lim <- max(abs(df$score_to_plot), na.rm = TRUE)
      if (!is.finite(lim) || lim == 0) lim <- 1
      if (use_masks) {
        poly_col <- join_spatial_polygon_colors(poly_df, df, "score_to_plot")
        req(!is.null(poly_col) && nrow(poly_col) > 0L)
        p <- base +
          geom_polygon(aes(group = cell_id, fill = score_to_plot),
                       data = poly_col, color = NA, alpha = mask_alpha) +
          scale_fill_gradient2(
            low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
            midpoint = 0, limits = c(-lim, lim),
            oob = scales::squish, name = legend_name
          )
      } else {
        p <- base +
          geom_point(data = df, aes(color = score_to_plot),
                     size = pt_size, alpha = pt_alpha) +
          scale_color_gradient2(
            low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
            midpoint = 0, limits = c(-lim, lim),
            oob = scales::squish, name = legend_name
          )
      }
    } else if (color_by == "cell_type" && "cell_type" %in% colnames(df)) {
      pal <- tryCatch(
        PhenoMapR::get_celltype_palette(as.character(unique(df$cell_type))),
        error = function(e) NULL
      )
      if (use_masks) {
        poly_col <- join_spatial_polygon_colors(poly_df, df, "cell_type")
        req(!is.null(poly_col) && nrow(poly_col) > 0L)
        p <- base +
          geom_polygon(aes(group = cell_id, fill = cell_type),
                       data = poly_col, color = NA, alpha = mask_alpha) +
          labs(fill = "Cell type")
        if (!is.null(pal)) p <- p + scale_fill_manual(values = pal)
      } else {
        p <- base +
          geom_point(aes(color = cell_type), size = pt_size, alpha = pt_alpha) +
          labs(color = "Cell type") +
          ggplot2::guides(color = ggplot2::guide_legend(
            override.aes = list(size = 4, alpha = 1)
          ))
        if (!is.null(pal)) p <- p + scale_color_manual(values = pal)
      }
    } else if (color_by == "source" && "source" %in% colnames(df)) {
      # Apply the PhenoMapR brand palette so the Source coloring here
      # matches the Source coloring on the Data tab and Score tab plots.
      # (source levels already factorized above, before ggplot() capture.)
      if (use_masks) {
        poly_col <- join_spatial_polygon_colors(poly_df, df, "source")
        req(!is.null(poly_col) && nrow(poly_col) > 0L)
        p <- base +
          geom_polygon(aes(group = cell_id, fill = source),
                       data = poly_col, color = NA, alpha = mask_alpha) +
          scale_fill_phenomapr_d() +
          labs(fill = "Source")
      } else {
        p <- base +
          geom_point(aes(color = source), size = pt_size, alpha = pt_alpha) +
          scale_color_phenomapr_d() +
          labs(color = "Source") +
          ggplot2::guides(color = ggplot2::guide_legend(
            override.aes = list(size = 4, alpha = 1)
          ))
      }
    } else if (color_by == "group" && "group" %in% colnames(df)) {
      # Display-only legend remap: the underlying group factor levels
      # stay "Most Adverse" / "Most Favorable" / "Other" (those strings
      # are the canonical output of define_phenotype_groups() and are
      # consumed verbatim by find_phenotype_markers() and friends), but
      # the LEGEND shows the user-facing "Most Phenotype +/-" labels.
      group_vals <- c(
        "Most Adverse"   = "#B2182B",
        "Most Favorable" = "#2166AC",
        "Other"          = "#BBBBBB"
      )
      group_labs <- c(
        "Most Adverse"   = "Most Phenotype +",
        "Most Favorable" = "Most Phenotype \u2212",
        "Other"          = "Other"
      )
      if (use_masks) {
        poly_col <- join_spatial_polygon_colors(poly_df, df, "group")
        req(!is.null(poly_col) && nrow(poly_col) > 0L)
        p <- base +
          geom_polygon(aes(group = cell_id, fill = group),
                       data = poly_col, color = NA, alpha = mask_alpha) +
          scale_fill_manual(
            values = group_vals, labels = group_labs,
            na.value = "#E0E0E0", name = "Group"
          )
      } else {
        p <- base +
          geom_point(aes(color = group), size = pt_size, alpha = pt_alpha) +
          scale_color_manual(
            values = group_vals, labels = group_labs,
            na.value = "#E0E0E0", name = "Group"
          ) +
          ggplot2::guides(color = ggplot2::guide_legend(
            override.aes = list(size = 4, alpha = 1)
          ))
      }
    } else {
      if (use_masks) {
        p <- base +
          geom_polygon(aes(group = cell_id), fill = "#555555",
                       color = NA, alpha = mask_alpha)
      } else {
        p <- base + geom_point(size = pt_size, alpha = pt_alpha, color = "#555555")
      }
    }

    if (isTRUE(input$umap_facet_source) && "source" %in% colnames(df)) {
      p <- p + facet_wrap(~ source)
    }
    panel_objects$umap_plot <- p
    p
  })
  phenomapr_register_plot_download(output, "umap_plot",
    function() isolate(panel_objects$umap_plot),
    width = 10, height = 8)

  output$download_umap_table <- downloadHandler(
    filename = function() {
      sprintf("phenomapr_embedding_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S"))
    },
    content = function(file) {
      emb <- current_embedding(); req(emb)
      df <- dplyr::left_join(emb, cell_table() %||% data.frame(cell_id = character(0)),
                             by = "cell_id")
      utils::write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  # ------------------------------------------------------------------------
  # Phenotype groups (merged into the Score tab UI)
  #
  # Tagging is now fully reactive: any time the user drags the percentile
  # slider, switches the score column, or re-computes scores, groups are
  # automatically re-derived. No explicit "Tag phenotype groups" click is
  # required. We also invalidate `state$markers` whenever groups change so
  # the Markers tab doesn't surface stale results.
  # ------------------------------------------------------------------------
  observe({
    scores <- state$scores
    if (is.null(scores)) {
      state$groups  <- NULL
      state$markers <- NULL
      return()
    }
    req(input$percentile)
    groups <- tryCatch(
      PhenoMapR::define_phenotype_groups(
        scores = scores,
        percentile = input$percentile
      ),
      error = function(e) {
        showNotification(paste0("define_phenotype_groups failed: ",
                                conditionMessage(e)),
                         type = "error", duration = 8); NULL
      }
    )
    if (is.null(groups)) return()

    # Merge with metadata when we have it, so downstream marker analysis can
    # use cell type / source labels. Treat NULL / "" / "(none)" all as
    # "not selected", and fall back to `.cell_id` when no column is mapped.
    id_col <- input$meta_cell_id_col
    if (is.null(id_col) || !nzchar(id_col) || id_col == "(none)") {
      id_col <- if (!is.null(state$metadata) && ".cell_id" %in% colnames(state$metadata)) {
        ".cell_id"
      } else if (!is.null(state$metadata) && ncol(state$metadata) >= 1L) {
        colnames(state$metadata)[1L]
      } else NA_character_
    }

    if (!is.null(state$metadata) && !is.na(id_col) &&
        id_col %in% colnames(state$metadata)) {
      md <- state$metadata
      mapped <- c(input$meta_cell_type_col, input$meta_source_col)
      mapped <- mapped[!vapply(mapped, function(x) is.null(x) || !nzchar(x) || x == "(none)",
                               logical(1))]
      keep_cols <- intersect(mapped, colnames(md))
      md$.cell_id <- as.character(md[[id_col]])
      if (length(keep_cols)) {
        md <- md[, c(".cell_id", keep_cols), drop = FALSE]
        groups$.cell_id <- groups$cell_id
        groups <- dplyr::left_join(groups, md, by = ".cell_id")
        groups$.cell_id <- NULL
      }
    }

    state$groups  <- groups
    state$markers <- NULL
  })

  output$group_summary <- renderUI({
    if (is.null(state$groups)) return(tags$p(tags$em("Compute scores to see group counts here.")))
    g <- state$groups
    score_col <- active_score_column(state$scores)
    grp_col <- resolve_phenotype_group_column(g, score_col)
    if (is.na(grp_col)) return(tags$p(tags$em("No phenotype group column found.")))
    tbl <- as.data.frame(table(g[[grp_col]], useNA = "ifany"),
                         stringsAsFactors = FALSE)
    colnames(tbl) <- c("Phenotype group", "Cells")
    tbl$Cells <- .fmt_int(tbl$Cells)
    # Display remap: the underlying R data layer ALWAYS uses the
    # canonical strings emitted by define_phenotype_groups()
    # ("Most Adverse" / "Most Favorable" / "Other") -- we must NOT
    # rewrite state$groups itself or we break downstream marker
    # finders, vignettes, and tests that rely on those literals.
    # Here we substitute the user-facing labels at render time only,
    # using HTML so we can bold the +/- markers without bringing in
    # a renderer (tag_table emits a plain <table>).
    display_map <- c(
      "Most Adverse"   = as.character(phenomapr_phenotype_plus()),
      "Most Favorable" = as.character(phenomapr_phenotype_minus()),
      "Other"          = "Other"
    )
    raw_vals <- as.character(tbl[["Phenotype group"]])
    mapped <- display_map[raw_vals]
    mapped[is.na(mapped)] <- raw_vals[is.na(mapped)]
    # tag_table() expects a data.frame with character cells; lapply
    # HTML(...) below tells Shiny to render the cell contents
    # verbatim (so <strong> and &minus; survive).
    tbl[["Phenotype group"]] <- lapply(mapped, HTML)
    tagList(
      tags$p(tags$strong("Group column: "), grp_col),
      tag_table(tbl)
    )
  })

  output$group_by_celltype_plot <- renderPlot({
    req(state$groups)
    g <- state$groups
    score_col <- active_score_column(state$scores)
    grp_col <- resolve_phenotype_group_column(g, score_col)
    if (is.na(grp_col)) {
      panel_objects$group_by_celltype_plot <- NULL
      return(NULL)
    }
    ct_col <- input$meta_cell_type_col
    if (is.null(ct_col) || ct_col == "(none)" || !ct_col %in% colnames(g)) {
      panel_objects$group_by_celltype_plot <- NULL
      return(NULL)
    }
    df <- g[, c(grp_col, ct_col)]
    colnames(df) <- c("group", "cell_type")
    df <- df[!is.na(df$group) & !is.na(df$cell_type), ]
    df_count <- df %>%
      dplyr::count(cell_type, group) %>%
      dplyr::group_by(cell_type) %>%
      dplyr::mutate(frac = n / sum(n)) %>%
      dplyr::ungroup()

    # Order cell types by ascending median PhenoMapR score so this
    # enrichment plot matches the x-axis order of the
    # `Score by cell type and source` boxplot directly above it. We
    # pull scores from `cell_table()` (the same reactive that drives
    # the boxplot); when scores aren't available yet we fall back to
    # alphabetical ordering on the cell types observed here. The
    # fallback also covers the rare case where the cell-type column
    # was set BEFORE scores landed.
    ct_levels <- {
      cell_tbl <- tryCatch(cell_table(), error = function(e) NULL)
      have_scores <- !is.null(cell_tbl) &&
        all(c("score", "cell_type") %in% colnames(cell_tbl)) &&
        any(is.finite(cell_tbl$score))
      if (have_scores) {
        med <- stats::aggregate(
          score ~ cell_type,
          data  = data.frame(
            score     = as.numeric(cell_tbl$score),
            cell_type = as.character(cell_tbl$cell_type),
            stringsAsFactors = FALSE
          ),
          FUN   = function(x) stats::median(x, na.rm = TRUE)
        )
        med <- med[order(med$score), , drop = FALSE]
        present <- as.character(unique(df_count$cell_type))
        ordered <- as.character(med$cell_type[
          as.character(med$cell_type) %in% present
        ])
        # Tail any cell types that show up here but not in the score
        # table (e.g. NA/0-score rows filtered out of cell_table)
        # alphabetically so they still appear instead of disappearing.
        tail <- sort(setdiff(present, ordered))
        c(ordered, tail)
      } else {
        sort(as.character(unique(df_count$cell_type)))
      }
    }
    df_count$cell_type <- factor(as.character(df_count$cell_type),
                                 levels = ct_levels)

    p <- ggplot(df_count, aes(x = cell_type, y = frac, fill = group)) +
      .geom_rounded_stack(color = NA) +
      scale_y_continuous(labels = scales::percent_format()) +
      # See score_rank_plot for rationale: data values stay as the
      # canonical "Most Adverse" / "Most Favorable" strings; only the
      # legend labels surface the user-facing "Most Phenotype +/-"
      # variants.
      scale_fill_manual(
        values = c(
          "Most Adverse"   = "#B2182B",
          "Most Favorable" = "#2166AC",
          "Other"          = "#BBBBBB"
        ),
        labels = c(
          "Most Adverse"   = "Most Phenotype +",
          "Most Favorable" = "Most Phenotype \u2212",
          "Other"          = "Other"
        )
      ) +
      labs(x = "Cell type", y = "Fraction of cells", fill = "Group") +
      theme_minimal(base_size = .theme_base_size()) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    panel_objects$group_by_celltype_plot <- p
    p
  })
  phenomapr_register_plot_download(output, "group_by_celltype_plot",
    function() isolate(panel_objects$group_by_celltype_plot),
    width = 9, height = 5)

  output$download_groups <- downloadHandler(
    filename = function() {
      sprintf("phenomapr_groups_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S"))
    },
    content = function(file) {
      req(state$groups)
      utils::write.table(state$groups, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  # ------------------------------------------------------------------------
  # 5. Markers
  # ------------------------------------------------------------------------
  observeEvent(input$run_markers, {
    req(state$expression, state$groups)
    phenomapr_busy_show(
      "Finding marker genes...",
      sprintf("Scope: %s", input$marker_scope %||% "across-all-cells")
    )

    score_col <- active_score_column(state$scores)
    grp_col <- resolve_phenotype_group_column(state$groups, score_col)
    if (is.na(grp_col)) {
      phenomapr_busy_hide()
      showNotification("Could not find a phenotype_group_* column for the selected score.",
                       type = "error", duration = 6); return()
    }
    ct_col_use <- input$meta_cell_type_col
    if (is.null(ct_col_use) || ct_col_use == "(none)") {
      ct_col_use <- NULL
    } else if (!ct_col_use %in% colnames(state$groups) &&
               !is.null(state$metadata) && ct_col_use %in% colnames(state$metadata)) {
      # groups$cell_id merge may use a differently named column; keep the
      # metadata column name for find_phenotype_markers() regardless.
    } else if (!ct_col_use %in% colnames(state$groups)) {
      ct_col_use <- NULL
    }

    # CRITICAL: defer the heavy find_phenotype_markers() call into
    # later::later() so the busy_show custom message can reach the
    # browser before libuv is blocked. See helpers.R\'s busy-overlay
    # file-comment for the full rationale.
    expr_snapshot   <- state$expression
    groups_snapshot <- state$groups

    # Decode the single sidebar radio into the (marker_scope,
    # celltype_contrast) tuple expected by find_phenotype_markers().
    # Three UI choices:
    #   "phenotype_groups"        -> cohort-wide
    #   "cell_type_specific"      -> per cell type, within-tail contrast
    #   "cell_type_vs_opposite"   -> per cell type, vs ALL cells in the
    #                                opposite phenotype tail
    scope_ui <- input$marker_scope %||% "phenotype_groups"
    if (identical(scope_ui, "cell_type_vs_opposite")) {
      marker_scope_use      <- "cell_type_specific"
      celltype_contrast_use <- "vs_opposite_tail"
    } else if (identical(scope_ui, "cell_type_specific")) {
      marker_scope_use      <- "cell_type_specific"
      celltype_contrast_use <- "within_cell_type"
    } else {
      marker_scope_use      <- "phenotype_groups"
      celltype_contrast_use <- "within_cell_type"
    }

    marker_args <- list(
      group_column        = grp_col,
      marker_scope        = marker_scope_use,
      celltype_contrast   = celltype_contrast_use,
      cell_type_column    = ct_col_use,
      min.pct             = input$marker_min_pct,
      logfc.threshold     = input$marker_logfc,
      pval_threshold      = input$marker_pval,
      max_cells_per_ident = input$marker_maxcells
    )
    sess <- session
    later::later(function() {
      markers <- tryCatch(
        PhenoMapR::find_phenotype_markers(
          expression = expr_snapshot,
          group_labels = groups_snapshot,
          group_column = marker_args$group_column,
          cell_id_column = "cell_id",
          marker_scope = marker_args$marker_scope,
          celltype_contrast = marker_args$celltype_contrast,
          cell_type_column = marker_args$cell_type_column,
          min.pct = marker_args$min.pct,
          logfc.threshold = marker_args$logfc.threshold,
          pval_threshold = marker_args$pval_threshold,
          max_cells_per_ident = marker_args$max_cells_per_ident,
          verbose = TRUE
        ),
        error = function(e) {
          showNotification(paste0("find_phenotype_markers failed: ",
                                  conditionMessage(e)),
                           type = "error", duration = 10, session = sess); NULL
        }
      )
      # Hide BEFORE assigning state$markers so the hide custom
      # message lands in a flush ahead of the heavy marker-table
      # renders (renderDT * 2, plus the marker heatmap UI).
      phenomapr_busy_hide(session = sess)
      if (is.null(markers)) return()
      later::later(function() {
        tryCatch({
          state$markers <- markers
          shiny::showNotification(
            tagList(
              "Found ",
              .fmt_int(nrow(markers$adverse_markers)),
              " ",
              phenomapr_phenotype_plus("phenotype +"),
              " and ",
              .fmt_int(nrow(markers$favorable_markers)),
              " ",
              phenomapr_phenotype_minus("phenotype \u2212"),
              " markers."
            ),
            type = "message", duration = 5, session = sess
          )
        }, error = function(e) {
          shiny::showNotification(
            paste0("Internal error while applying markers: ",
                   conditionMessage(e)),
            type = "error", duration = 10, session = sess
          )
        })
      }, delay = 0)
    }, delay = 0)
  })

  # The adverse / favorable marker tables can run into the hundreds of
  # rows for cell-type-specific or large-cohort panels. Without a
  # vertical scroll cap, the DT widget renders at full content height
  # and pushes through (visually overlapping) the Marker-gene heatmap
  # card directly below it. Bound each table to ~480px of vertical
  # space and let DataTables' native scroller handle scrolling within
  # that window -- pagination + filter / search controls remain in the
  # toolbar above so users can still page through long result sets.
  .markers_dt_options <- list(
    pageLength = 15,
    scrollX = TRUE,
    scrollY = "480px",
    scrollCollapse = TRUE
  )

  output$adverse_markers_tbl <- renderDT({
    req(state$markers)
    df <- state$markers$adverse_markers
    panel_objects$adverse_markers_tbl <- df
    if (!nrow(df)) return(datatable(df, options = list(dom = "t")))
    datatable(df, rownames = FALSE, options = .markers_dt_options)
  })
  phenomapr_register_table_download(output, "adverse_markers_tbl",
    function() isolate(panel_objects$adverse_markers_tbl))
  output$favorable_markers_tbl <- renderDT({
    req(state$markers)
    df <- state$markers$favorable_markers
    panel_objects$favorable_markers_tbl <- df
    if (!nrow(df)) return(datatable(df, options = list(dom = "t")))
    datatable(df, rownames = FALSE, options = .markers_dt_options)
  })
  phenomapr_register_table_download(output, "favorable_markers_tbl",
    function() isolate(panel_objects$favorable_markers_tbl))

  output$download_markers <- downloadHandler(
    filename = function() sprintf("phenomapr_markers_%s.rds", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      req(state$markers); saveRDS(state$markers, file)
    }
  )

  # ------------------------------------------------------------------------
  # Marker-gene heatmap (uses plot_phenotype_markers + ComplexHeatmap)
  # ------------------------------------------------------------------------
  # We cache the *inputs* to plot_phenotype_markers() — not the returned
  # Heatmap object — so the actual draw happens inside the function with its
  # full set of merged legends and right-side padding. Calling
  # `ComplexHeatmap::draw()` directly on the returned object skipped the
  # custom annotation legends (HeatmapAnnotation uses show_legend = FALSE)
  # and the wide padding needed for left-side anno_mark labels, which was
  # clipping both the gene labels and the colour-bar legends in the app.
  marker_heatmap_args <- reactiveVal(NULL)

  observeEvent(input$draw_marker_heatmap, {
    req(state$markers, state$expression, state$scores)
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
        !requireNamespace("circlize", quietly = TRUE)) {
      showNotification(
        "Install ComplexHeatmap + circlize for marker heatmaps.",
        type = "warning", duration = 8
      ); return()
    }
    phenomapr_busy_show(
      "Building marker heatmap...",
      "Assembling expression matrix and per-cell metadata"
    )
    # Defer the heavy assembly work (extract_expression_matrix +
    # joins) into later::later() so the busy_show custom message
    # can reach the browser before libuv gets blocked. See
    # helpers.R\'s busy-overlay file-comment for the rationale.
    # The popup will stay up through the renderImage draw too,
    # because that happens inside this same logical "busy"
    # window from the user\'s perspective.
    expr_obj_snapshot   <- state$expression
    markers_snapshot    <- state$markers
    scores_snapshot     <- state$scores
    groups_snapshot     <- state$groups
    metadata_snapshot   <- state$metadata
    score_assay_input   <- input$score_assay %||% ""
    score_slot_input    <- input$score_slot   %||% "data"
    marker_scope_input  <- input$marker_scope
    meta_ct_col_input   <- input$meta_cell_type_col
    meta_id_col_input   <- input$meta_cell_id_col
    hm_top_n_input      <- input$hm_top_n     %||% 20L
    hm_n_labels_input   <- input$hm_n_labels  %||% 5L
    hm_block_borders_in <- isTRUE(input$hm_block_borders)
    hm_color_mark_labels_in <- isTRUE(input$hm_color_mark_labels)
    hm_color_schemes_in <- marker_heatmap_color_schemes_from_input(input)
    score_col_active <- active_score_column(scores_snapshot)
    marker_pval_input   <- input$marker_pval %||% 0.05
    sess <- session
    later::later(function() {
      # For AnnData inputs, the marker heatmap only needs the marker genes,
      # we pass that list straight to extract_expression_matrix() so that
      # subsetting happens on the Python side and we don't materialise the
      # full expression matrix in R.
      marker_genes <- unique(c(
        if (!is.null(markers_snapshot$adverse_markers))   markers_snapshot$adverse_markers$gene,
        if (!is.null(markers_snapshot$favorable_markers)) markers_snapshot$favorable_markers$gene
      ))
      marker_genes <- marker_genes[nzchar(marker_genes)]

      expr_mat <- tryCatch(
        extract_expression_matrix(
          expr_obj_snapshot,
          assay = if (nzchar(score_assay_input)) score_assay_input else NULL,
          slot = score_slot_input,
          gene_subset = if (length(marker_genes)) marker_genes else NULL
        ),
        error = function(e) {
          showNotification(paste0("Could not extract expression matrix: ",
                                  conditionMessage(e)),
                           type = "error", duration = 8, session = sess); NULL
        }
      )
      if (is.null(expr_mat)) {
        phenomapr_busy_hide(session = sess)
        return()
      }

      grp_col <- resolve_phenotype_group_column(groups_snapshot, score_col_active)
      if (is.na(grp_col)) {
        phenomapr_busy_hide(session = sess)
        showNotification("No phenotype_group_* column for the selected score \u2014 re-tag groups first.",
                         type = "error", duration = 6, session = sess); return()
      }
      score_name <- score_col_active
      if (!nzchar(score_name) || !(score_name %in% colnames(scores_snapshot))) {
        score_name <- colnames(scores_snapshot)[1L]
      }
      scores_df <- data.frame(
        cell_id = rownames(scores_snapshot),
        .score = as.numeric(scores_snapshot[[score_name]]),
        stringsAsFactors = FALSE
      )
      colnames(scores_df)[2L] <- score_name
      meta_df <- data.frame(
        cell_id = as.character(groups_snapshot$cell_id),
        stringsAsFactors = FALSE
      )
      meta_df$.group <- groups_snapshot[[grp_col]]
      colnames(meta_df)[2L] <- grp_col
      meta_df <- dplyr::left_join(meta_df, scores_df, by = "cell_id")

      ct_col_use <- meta_ct_col_input
      has_ct <- !is.null(ct_col_use) && nzchar(ct_col_use) && ct_col_use != "(none)" &&
        !is.null(metadata_snapshot) && ct_col_use %in% colnames(metadata_snapshot)
      if (has_ct) {
        md <- metadata_snapshot
        id_col <- meta_id_col_input
        if (!nzchar(id_col) || !(id_col %in% colnames(md))) id_col <- ".cell_id"
        if (!(id_col %in% colnames(md))) id_col <- colnames(md)[1L]
        md$.cell_id <- as.character(md[[id_col]])
        ct_df <- data.frame(
          cell_id = md$.cell_id,
          .ct = as.character(md[[ct_col_use]]),
          stringsAsFactors = FALSE
        )
        colnames(ct_df)[2L] <- ct_col_use
        meta_df <- dplyr::left_join(meta_df, ct_df, by = "cell_id")
      }

      # Both per-cell-type contrasts ("cell_type_specific" /
      # "cell_type_vs_opposite") render as a cell-type-specific
      # heatmap; the underlying marker tables already carry a
      # cell_type column for either contrast.
      is_ct <- (marker_scope_input %in%
                  c("cell_type_specific", "cell_type_vs_opposite")) &&
               has_ct
      heatmap_type <- if (is_ct) "cell_type_specific" else "global"
      if (identical(marker_scope_input, "cell_type_vs_opposite")) {
        celltype_contrast_hm <- "vs_opposite_tail"
      } else if (identical(marker_scope_input, "cell_type_specific")) {
        celltype_contrast_hm <- "within_cell_type"
      } else {
        celltype_contrast_hm <- "vs_opposite_tail"
      }

      phenomapr_busy_hide(session = sess)
      later::later(function() {
        marker_heatmap_args(list(
          markers = markers_snapshot,
          expr_mat = expr_mat,
          meta = meta_df,
          cell_id_col = "cell_id",
          group_col = grp_col,
          score_col = score_name,
          celltype_col = if (has_ct) ct_col_use else NULL,
          heatmap_type = heatmap_type,
          celltype_contrast = celltype_contrast_hm,
          top_n_markers = hm_top_n_input,
          n_mark_labels = hm_n_labels_input,
          p_adj_threshold = marker_pval_input,
          color_mark_labels_by_celltype = hm_color_mark_labels_in,
          # Block-outline parameters captured at draw time. The
          # renderImage handler also consults the live
          # `input$hm_block_borders` value, so toggling the
          # checkbox after the heatmap exists triggers a
          # cheap re-render without rebuilding the cached
          # expression matrix or metadata.
          outline_marker_blocks = hm_block_borders_in,
          block_outline_color = if (hm_block_borders_in) "black" else "white",
          # 1.5pt is slightly slimmer than the previous 2pt default
          # so the borders read as clean dividers instead of heavy
          # chrome around each block.
          block_outline_lwd = if (hm_block_borders_in) 1.5 else 1,
          color_schemes = hm_color_schemes_in,
          use_raster = FALSE
        ))
        showNotification("Heatmap ready -- drawing...", type = "message",
                         duration = 3, session = sess)
      }, delay = 0)
    }, delay = 0)
  })

  output$marker_heatmap <- renderImage(
    {
      args <- marker_heatmap_args(); req(args)

      # Live block-border toggle: the checkbox feeds renderImage
      # directly so flipping it re-renders the cached heatmap with
      # the new outline color without rebuilding the (expensive)
      # expression matrix or metadata join. The cached `args` was
      # populated at "Draw heatmap" time -- we just override the
      # outline knobs here from the current input value.
      borders_on <- isTRUE(input$hm_block_borders)
      args$outline_marker_blocks <- borders_on
      args$block_outline_color   <- if (borders_on) "black" else "white"
      # Slightly slimmer line (1.5pt) than the original 2pt default.
      args$block_outline_lwd     <- if (borders_on) 1.5 else 1
      args$color_mark_labels_by_celltype <- isTRUE(input$hm_color_mark_labels)
      args$color_schemes <- marker_heatmap_color_schemes_from_input(input)

      # Render to a real PNG at fixed pixel size + DPI. ComplexHeatmap +
      # anno_mark are notoriously fragile under Shiny's renderPlot()
      # replay (the labels render the first time but get blanked out
      # on resize / re-render). renderImage() bypasses replays entirely
      # and the explicit res = 110 keeps font metrics matched between
      # the layout pass and the text grobs that draw the labels.
      width  <- session$clientData$output_marker_heatmap_width  %||% 1400
      height <- session$clientData$output_marker_heatmap_height %||% 640
      width  <- max(800,  as.integer(width))
      height <- max(480,  as.integer(height))

      tmp <- tempfile(fileext = ".png")
      ok <- tryCatch({
        grDevices::png(tmp, width = width, height = height, res = 110)
        on.exit(if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(),
                                                   silent = TRUE),
                add = TRUE)
        do.call(
          PhenoMapR::plot_phenotype_markers,
          c(args, list(draw = TRUE))
        )
        TRUE
      }, error = function(e) {
        showNotification(paste0("plot_phenotype_markers failed: ",
                                conditionMessage(e)),
                         type = "error", duration = 10)
        FALSE
      })
      if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
      if (!isTRUE(ok)) {
        return(list(src = "", contentType = "image/png",
                    width = width, height = height,
                    alt = "Marker heatmap unavailable."))
      }
      panel_objects$marker_heatmap <- list(width = width, height = height)
      list(
        src         = tmp,
        contentType = "image/png",
        width       = width,
        height      = height,
        alt         = "Marker gene heatmap"
      )
    },
    deleteFile = TRUE
  )

  # Heatmap download: re-render to a PNG at the current panel size from
  # the cached `marker_heatmap_args`. Going through plot_phenotype_markers()
  # again (rather than reusing the live renderImage tempfile) is the
  # only reliable path -- the imageOutput is mounted with
  # deleteFile = TRUE, so the temp PNG is removed as soon as the
  # browser fetches it.
  output$marker_heatmap_download <- downloadHandler(
    filename = function() phenomapr_dl_filename("marker_heatmap", "png"),
    content = function(file) {
      args <- isolate(marker_heatmap_args())
      dims <- isolate(panel_objects$marker_heatmap)
      width  <- (dims$width  %||% 1400)
      height <- (dims$height %||% 640)
      if (is.null(args)) {
        .phenomapr_write_placeholder_png(
          file, "Draw the heatmap first, then re-try the download."
        )
        return(invisible(NULL))
      }
      # Mirror the live block-border toggle so the downloaded PNG
      # matches what the user currently sees on screen.
      borders_on <- isTRUE(isolate(input$hm_block_borders))
      args$outline_marker_blocks <- borders_on
      args$block_outline_color   <- if (borders_on) "black" else "white"
      # Slightly slimmer line (1.5pt) than the original 2pt default.
      args$block_outline_lwd     <- if (borders_on) 1.5 else 1
      args$color_mark_labels_by_celltype <- isTRUE(isolate(input$hm_color_mark_labels))
      args$color_schemes <- marker_heatmap_color_schemes_from_input(isolate(input))
      ok <- tryCatch({
        grDevices::png(file, width = width, height = height, res = 110)
        on.exit(if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(),
                                                   silent = TRUE),
                add = TRUE)
        do.call(PhenoMapR::plot_phenotype_markers,
                c(args, list(draw = TRUE)))
        TRUE
      }, error = function(e) {
        if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
        FALSE
      })
      if (!isTRUE(ok)) {
        .phenomapr_write_placeholder_png(
          file, "plot_phenotype_markers() failed for this download."
        )
      }
      invisible(NULL)
    }
  )

  # ------------------------------------------------------------------------
  # Reference diagnostics (rendered at the bottom of the Phenotype tab —
  # was a standalone "6. Diagnostics" nav panel before). The reference and
  # cancer type follow the scoring-signature picker at the top of the
  # sidebar (`reference_choice` + `cancer_type`); only the diagnostic
  # knobs (top-N / direction) are local. This way the diagnostics auto-
  # propagate whenever the user changes the phenotype signature — no
  # explicit "Show top prognostic genes" button is needed.
  # ------------------------------------------------------------------------
  output$top_genes_tbl <- renderDT({
    rc <- input$reference_choice
    # Custom signatures don't ship with cancer-type-resolved gene tables.
    # Show a friendly placeholder instead of failing.
    if (is.null(rc) || identical(rc, "_custom")) {
      placeholder <- data.frame(
        note = paste0(
          "Reference diagnostics are only available for the built-in ",
          "signatures (PRECOG, TCGA, Pediatric PRECOG, ICI PRECOG). ",
          "Switch the 'Phenotype signature source' above to view top ",
          "prognostic genes."
        )
      )
      panel_objects$top_genes_tbl <- placeholder
      return(datatable(
        placeholder,
        rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE,
                       searching = FALSE)
      ))
    }
    # Guard against the same async-update race the signature plot hits:
    # `input$reference_choice` flips synchronously when the user clicks a
    # pill, but `input$cancer_type` only catches up after the round-trip
    # update of its select choices. Skip rendering until the two are in
    # sync, so we never call `get_top_prognostic_genes(reference="tcga",
    # cancer_type="Adrenocortical")` and similar.
    valid_cts <- get_cancer_types(rc)
    req(input$cancer_type, input$cancer_type %in% valid_cts)
    tg <- tryCatch(
      PhenoMapR::get_top_prognostic_genes(
        reference   = rc,
        cancer_type = input$cancer_type,
        # `n = Inf` returns the entire signature -- the per-table search
        # box (enabled below) is the new way users find specific genes,
        # so we no longer cap with a sidebar top-N input.
        n           = Inf,
        # Direction radio was removed at the user's request; show all
        # genes (both positive and negative z-scores) ordered by |z|.
        # Users filter by sign in DT via the z_score column header
        # sort + search.
        direction   = "both"
      ),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8)
        NULL
      }
    )
    req(!is.null(tg))
    panel_objects$top_genes_tbl <- tg
    datatable(
      tg,
      rownames = FALSE,
      # No pagination -- show every gene in the signature as one
      # scrollable list contained inside the card body. dom = "ft"
      # drops the length menu ("Show N entries"), the info row
      # ("Showing N of N entries"), and the pagination buttons; only
      # the search box at the top and the table itself are rendered.
      # The fixed scrollY keeps the scroll affordance inside the
      # card so the row layout stays exactly 320 px tall and the
      # rest of the page does not move when the user scrolls the
      # gene list.
      options = list(
        paging          = FALSE,
        dom             = "ft",
        scrollY         = "210px",
        scrollCollapse  = TRUE,
        scrollX         = TRUE,
        searching       = TRUE,
        searchHighlight = TRUE
      )
    )
  })
  phenomapr_register_table_download(output, "top_genes_tbl",
    function() isolate(panel_objects$top_genes_tbl))

  output$cancer_types_list <- renderPrint({
    rc <- input$reference_choice
    if (is.null(rc) || identical(rc, "_custom")) {
      cat("Custom phenotype signature selected.\n",
          "Built-in cancer-type lists are not applicable for custom",
          "signatures.\n", sep = "")
      return(invisible(NULL))
    }
    cts <- get_cancer_types(rc)
    cat(sprintf("Reference: %s\n  %d cancer types available\n\n",
                rc, length(cts)))
    print(cts)
  })

  # ------------------------------------------------------------------------
  # Session cleanup
  # ------------------------------------------------------------------------
  session$onSessionEnded(function() {
    gc(verbose = FALSE)
  })
}

shinyApp(ui = ui, server = server)
