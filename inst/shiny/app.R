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
  required <- c("shiny", "bslib", "DT", "ggplot2", "dplyr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "The PhenoMapR Shiny app needs these R packages installed: ",
      paste(missing, collapse = ", "),
      ". Install them with: install.packages(c(\"",
      paste(missing, collapse = "\", \""), "\"))."
    )
  }
  library(shiny)
  library(bslib)
  library(DT)
  library(ggplot2)
  library(dplyr)
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
    tags$link(rel = "stylesheet", href = "styles.css"),
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
      card_header("Overview of PhenoMapR Shiny App"),
      card_body(
        layout_columns(
          col_widths = c(7, 5),
          tagList(
            markdown(
              "**PhenoMapR** maps phenotypic signal from bulk expression onto
              single-cell, spatial, and pseudobulk transcriptomics data using
              weighted-sum scoring against curated prognostic z-scores
              (PRECOG, TCGA, Pediatric PRECOG, ICI PRECOG) or a custom
              reference you supply."
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
                        HTML("RDS (Seurat / SCE / matrix), TSV/CSV ",
                             "(genes &times; samples), or 10X HDF5 / ",
                             "AnnData <code>.h5ad</code>. Add optional ",
                             "cell metadata for cell-type-specific ",
                             "analyses.")
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
                        HTML("Pick a built-in cohort meta-z ",
                             "(PRECOG / TCGA / Pediatric / ICI), upload ",
                             "a precomputed signature file, or ",
                             "<em>derive</em> a custom signature from your ",
                             "own bulk expression + phenotype on the ",
                             "<strong>Phenotype</strong> tab (binary, ",
                             "continuous, or survival outcomes).")
                      )
                    )
                  ),
                  tags$li(
                    class = "welcome-step",
                    tags$div(class = "welcome-step-num", "3"),
                    tags$div(
                      class = "welcome-step-text",
                      tags$div(class = "welcome-step-title", "Score"),
                      tags$div(
                        class = "welcome-step-desc",
                        HTML("Compute weighted-sum PhenoMapR scores for each ",
                             "cell / sample / spot, view histograms, ",
                             "score-rank plots, and per-cell-type &times; ",
                             "source Wilcoxon comparisons, then immediately ",
                             "tag the top / bottom <em>N</em>&nbsp;% of cells ",
                             "as <strong>Most Adverse</strong> / ",
                             "<strong>Most Favorable</strong> phenotype ",
                             "groups.")
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
                        HTML("Overlay PhenoMapR score, cell type, source, ",
                             "or phenotype-group label onto any 2D ",
                             "embedding stored on your object (UMAP, tSNE, ",
                             "PCA, or tissue / spot coordinates for ",
                             "spatial inputs) or an uploaded coordinate ",
                             "file.")
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
                        HTML("Run differential expression for the adverse ",
                             "vs. favorable tails (globally or per cell ",
                             "type) and draw a ComplexHeatmap of the top ",
                             "markers per block.")
                      )
                    )
                  )
                )
              )
            ),

            # "Server / remote use" — companion card matching the stepper
            # chrome, with the code block visually called out so users can
            # spot the remote-launch command at a glance.
            tags$div(
              class = "welcome-section welcome-server",
              tags$div(
                class = "welcome-section-header",
                icon("server"),
                tags$span("Server / remote use")
              ),
              tags$div(
                class = "welcome-section-body",
                tags$p(
                  HTML("This app is just a thin UI over the regular ",
                       "PhenoMapR API. Run it on a remote workstation with:")
                ),
                # Emit <pre><code>…</code></pre> as a single HTML string so
                # htmltools' default pretty-printing doesn't insert newlines
                # between the <pre> and <code> tags — those newlines are
                # preserved by `white-space: pre` on .welcome-code and would
                # render as a stray blank line below the command.
                HTML(paste0(
                  '<pre class="welcome-code"><code>',
                  'PhenoMapR::run_app(host = "0.0.0.0", port = 3838)',
                  '</code></pre>'
                )),
                tags$p(
                  HTML("and open ",
                       "<code>http://&lt;server&gt;:3838</code> from your ",
                       "browser. <em>No data leaves the machine the app ",
                       "runs on.</em>")
                )
              )
            )
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
        h4("Expression input"),
        phenomapr_file_input(
          "expr_file",
          label = NULL,
          accept = c(".rds", ".h5", ".h5ad", ".tsv", ".csv", ".txt"),
          width = "100%"
        ),
        actionLink("use_demo", "Use a tiny demo matrix instead",
                   icon = icon("flask")),

        hr(),
        h4("Cell metadata"),
        helpText(
          "Metadata is read automatically from Seurat / SingleCellExperiment / ",
          "AnnData objects. PhenoMapR guesses sensible defaults for the cell-ID, ",
          "cell-type, and source columns — adjust them if needed."
        ),
        uiOutput("metadata_status"),
        # The metadata-upload panel auto-expands when an expression file has
        # been loaded but no metadata could be extracted from it, so the user
        # is prompted to attach a tabular metadata file before picking the
        # cell-ID / cell-type / source columns below.
        uiOutput("metadata_upload_panel"),
        selectInput("meta_cell_id_col", "Cell ID column", choices = NULL),
        selectInput("meta_cell_type_col", "Cell type column (optional)",
                    choices = NULL),
        selectInput("meta_source_col", "Source / group column (optional)",
                    choices = NULL)
      ),
      card(
        card_header(icon("circle-info"), " About data input"),
        card_body(
          markdown(
            "PhenoMapR accepts:

            * **Matrix / data.frame** — rows = genes, columns = cells/samples.
            * **Seurat** (v4 / v5; assays `RNA`, `Spatial`).
            * **SingleCellExperiment** / **SpatialExperiment**.
            * **AnnData** (via `reticulate`).

            For matrix uploads, gene IDs must be HUGO symbols (e.g. `TP53`).
            ENSG IDs are detected and flagged; convert with `biomaRt` /
            `AnnotationDbi` before uploading.

            By default this app accepts uploads of any size — Shiny's 5 MB
            cap is removed so very large `.rds` / `.h5ad` objects work as
            long as your machine has enough RAM to hold them. Cap the limit
            with `options(shiny.maxRequestSize = <bytes>)` before launching
            if needed."
          )
        )
      ),
      uiOutput("expr_axes_warning"),

      # ---- Dataset overview --------------------------------------------------
      card(
        card_header(icon("table-list"), " Dataset overview"),
        card_body(
          uiOutput("dataset_overview_summary"),
          layout_columns(
            col_widths = c(6, 6),
            card(
              card_header(tags$strong("Cells per cell type")),
              card_body(plotOutput("celltype_count_plot", height = "260px"))
            ),
            card(
              card_header(tags$strong("Cells per source / group")),
              card_body(plotOutput("source_count_plot", height = "260px"))
            )
          ),
          layout_columns(
            col_widths = c(6, 6),
            card(
              card_header(tags$strong("Cell type × source composition")),
              card_body(plotOutput("celltype_source_plot", height = "260px"))
            ),
            card(
              card_header(tags$strong("Metadata columns")),
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
          selectInput("cancer_type", "Cancer / tissue type", choices = NULL),
          sliderInput("z_score_cutoff", "Signature |z| cutoff",
                      min = 0, max = 5, value = 2, step = 0.1),
          radioButtons("reference_sign", "Direction",
                       choices = c("Higher score = worse prognosis" = 1L,
                                   "Higher score = better prognosis" = -1L),
                       selected = 1L)
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
            helpText(
              "Upload bulk expression (genes × samples) and a sample-level ",
              "phenotype table. PhenoMapR will compute per-gene z-scores ",
              "(Cox for survival, logistic for binary, correlation for ",
              "continuous) that you can immediately use to score your ",
              "single-cell / spatial data."
            ),
            phenomapr_file_input("derive_bulk_file", "Bulk expression (genes × samples)",
                                 accept = c(".rds", ".tsv", ".csv", ".txt")),
            uiOutput("derive_bulk_summary"),
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
        ),
        # ---- Reference diagnostics inputs (moved from the old "6. Diagnostics"
        # tab so all reference-signature info lives in one place). The
        # reference + cancer type are NOT re-asked here — the diagnostics
        # below automatically follow the "Phenotype signature source" /
        # "Cancer / tissue type" selection at the top of this sidebar. Only
        # the diagnostic-specific knobs (top-N and direction) live here.
        hr(),
        h4(tagList(icon("stethoscope"), " Reference diagnostics")),
        helpText(
          "Top prognostic genes for the phenotype signature selected ",
          "above. Adjust how many genes to list and the direction filter."
        ),
        numericInput("top_genes_n", "Top N genes",
                     value = 50, min = 5, max = 1000),
        radioButtons(
          "top_genes_dir", "Direction",
          choices = c("both" = "both", "positive" = "positive",
                      "negative" = "negative"),
          selected = "both",
          inline = TRUE
        )
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
              tags$span(class = "ref-name ref-name-precog", "PRECOG"),
              " - ", tags$strong("166 studies"),
              ", ~", tags$strong("18,000 adult patients"),
              " across ", tags$strong("39 cancer histologies"), "."
            ),
            tags$li(
              tags$span(class = "ref-name ref-name-tcga", "TCGA"),
              " - single-cohort prognostic z, ~",
              tags$strong("11,000 patients"),
              " across ", tags$strong("33 adult cancer types"), "."
            ),
            tags$li(
              tags$span(class = "ref-name ref-name-pediatric",
                        "Pediatric PRECOG"),
              " - ", tags$strong("32 studies"),
              ", ~", tags$strong("4,000 pediatric patients"),
              " across ", tags$strong("12 cancers"), "."
            ),
            tags$li(
              tags$span(class = "ref-name ref-name-ici", "ICI PRECOG"),
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
      # Side-by-side: signature status + gene coverage. Both are small
      # status-style cards, so showing them on one row saves vertical
      # space and the signature plot below sits closer to the top.
      layout_columns(
        col_widths = c(6, 6),
        card(
          fill = FALSE,
          card_header(tags$strong("Phenotype signature status")),
          card_body(uiOutput("reference_status"))
        ),
        card(
          fill = FALSE,
          class = "compact-coverage-card",
          card_header(tags$strong("Gene coverage")),
          card_body(
            class = "compact-coverage-body",
            helpText(
              "Fraction of your expression genes that overlap the chosen signature."
            ),
            DTOutput("gene_coverage_tbl")
          )
        )
      ),
      card(
        card_header(tags$strong("Phenotype signature")),
        card_body(
          helpText("Top / bottom z-score genes for the selected signature column."),
          plotOutput("reference_signature_plot", height = "260px")
        )
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
                card_header(tags$strong("Top positive (adverse-direction) genes")),
                card_body(DTOutput("derived_top_pos_tbl"))
              ),
              card(
                card_header(tags$strong("Top negative (favorable-direction) genes")),
                card_body(DTOutput("derived_top_neg_tbl"))
              )
            )
          )
        )
      ),
      # ---- Reference diagnostics (folded in from the old "6. Diagnostics") -
      # Lets users inspect any built-in reference's top prognostic genes /
      # cancer-type list right next to the signature they're using to score,
      # without bouncing to a separate tab.
      tags$div(
        class = "phenotype-section-divider",
        tags$h3(tagList(icon("stethoscope"), " Reference diagnostics"))
      ),
      card(
        card_header(icon("circle-info"), " About reference diagnostics"),
        card_body(
          markdown(
            "These tables follow the **Phenotype signature source** and
            **Cancer / tissue type** you picked at the top of the sidebar
            — they update automatically whenever you change the signature.
            Use the *Top N genes* and *Direction* knobs in the sidebar to
            tune the listing."
          )
        )
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(
          card_header(tags$strong("Top prognostic genes")),
          card_body(
            helpText(
              "Auto-updates from the phenotype signature selected above. ",
              "Adjust top-N and direction in the sidebar."
            ),
            DTOutput("top_genes_tbl")
          )
        ),
        card(
          card_header(tags$strong("Available cancer types")),
          card_body(verbatimTextOutput("cancer_types_list"))
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
        h4("PhenoMap() parameters"),
        radioButtons(
          "score_slot", "Seurat / SCE layer",
          choices = c("data (log-normalized)" = "data",
                      "counts (raw)" = "counts",
                      "scale.data" = "scale.data"),
          selected = "data"
        ),
        textInput("score_assay", "Assay (Seurat / SCE; blank = auto)", value = ""),
        checkboxInput("pseudobulk", "Aggregate to pseudobulk", value = FALSE),
        conditionalPanel(
          "input.pseudobulk == true",
          selectInput("pseudobulk_group_by",
                      "Group cells by (column in metadata)",
                      choices = NULL)
        ),
        actionButton("run_score", "Compute PhenoMapR scores",
                     icon = icon("play"), class = "btn-primary"),
        downloadButton("download_scores", "Download score table (TSV)",
                       class = "btn-outline-primary"),
        hr(),
        h4("Phenotype groups"),
        helpText(
          "After scoring, the top / bottom percentile of cells are ",
          "automatically tagged as Most Phenotype + and Most Phenotype - ",
          "as you drag the slider. These labels feed the marker-finding ",
          "step and downstream visualizations."
        ),
        sliderInput("percentile", "Tail percentile (top and bottom)",
                    min = 0.01, max = 0.40, value = 0.05, step = 0.01),
        selectInput("groups_score_col",
                    "Score column to tag (if multiple)",
                    choices = NULL),
        downloadButton("download_groups", "Download labels (TSV)",
                       class = "btn-outline-primary")
      ),
      card(
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
          card_header(tags$strong("Score distribution")),
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
          card_header(tags$strong("Cells ordered by PhenoMapR score")),
          card_body(
            helpText(
              "Each cell is shown as a point at its rank along the score axis. ",
              "Use the colour control to highlight cell types or source / group ",
              "annotations (if cell metadata has been provided)."
            ),
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
      card(
        card_header(tags$strong("Score by cell type and source")),
        card_body(
          helpText(
            "Boxplot of (scaled) PhenoMapR scores per cell type, ordered ",
            "from lowest to highest median score. When a source / group ",
            "column is also mapped on the Data tab, boxes are split by ",
            "source within each cell type and Wilcoxon brackets compare ",
            "sources within each cell type. When only a cell-type column ",
            "is mapped, the plot shows one box per cell type and annotates ",
            "pairwise Wilcoxon tests between cell types. Only significant ",
            "(p < 0.05) brackets are drawn — brackets stacked by p-value."
          ),
          plotOutput("score_box_source_plot", height = "440px")
        )
      ),
      card(
        card_header(tags$strong("Score table")),
        card_body(DTOutput("score_table"))
      ),

      # ---- Phenotype groups (merged from the old standalone tab) ----------
      tags$div(
        class = "score-section-divider",
        tags$h3(tagList(icon("layer-group"), " Phenotype groups"))
      ),
      card(
        card_header(icon("circle-info"), " About phenotype tails"),
        card_body(
          markdown(
            "PhenoMapR partitions cells into **Most Adverse** (top *N* %),
            **Most Favorable** (bottom *N* %), and **Other** based on the
            chosen score column. These labels feed the marker-finding step
            and any downstream visualizations.

            Lower the percentile to make the tails tighter (more
            discriminative but fewer cells); raise it to include more cells
            (statistically more robust but biologically fuzzier)."
          )
        )
      ),
      card(
        card_header(tags$strong("Group sizes")),
        card_body(uiOutput("group_summary"))
      ),
      card(
        card_header(tags$strong("Per-cell-type group enrichment")),
        card_body(
          helpText("Only shown when a cell-type column has been selected above."),
          plotOutput("group_by_celltype_plot", height = "320px")
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
        h4("Embedding"),
        # Collapse the long "how detection works" copy into a <details>
        # block so the sidebar isn't dominated by a wall of text by
        # default. The Reduction dropdown is the primary control and now
        # sits right under the section header.
        selectInput("umap_reduction", "Reduction", choices = NULL),
        uiOutput("umap_reduction_status"),
        tags$details(
          class = "embedding-help-details",
          tags$summary("About reductions & auto-detection"),
          helpText(
            "Reductions stored on a Seurat / SingleCellExperiment / ",
            "AnnData object are picked up automatically. Spatial inputs ",
            "(Seurat with @images, SpatialExperiment, AnnData with ",
            "obsm['spatial']) also surface tissue coordinates as a ",
            "\"spatial\" reduction in the dropdown above. For matrix ",
            "uploads, the app additionally auto-detects 2D coordinate ",
            "column pairs (e.g. UMAP_1/UMAP_2, tSNE_1/tSNE_2, ",
            "X_umap_0/X_umap_1) on the metadata table from the Data tab ",
            "\u2014 no second upload needed. Use the file picker below ",
            "to supply a separate embedding file if you don't have one ",
            "in your metadata."
          )
        ),
        tags$details(
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
        checkboxInput("umap_facet_source", "Facet by source", value = FALSE),
        downloadButton("download_umap_table", "Download embedding (TSV)",
                       class = "btn-outline-primary")
      ),
      card(
        card_header(tags$strong("Embedding")),
        card_body(plotOutput("umap_plot", height = "560px"))
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
        radioButtons(
          "marker_scope", "Scope",
          choices = c("Cohort-wide (adverse vs favorable)" = "phenotype_groups",
                      "Cell-type specific"                 = "cell_type_specific"),
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
        hr(),
        h4("Marker-gene heatmap"),
        helpText(
          "After running marker discovery, configure and draw a ",
          "ComplexHeatmap of the top markers per block."
        ),
        numericInput("hm_top_n", "Top genes per block", value = 20,
                     min = 5, max = 200, step = 5),
        numericInput("hm_n_labels", "Gene labels per block", value = 5,
                     min = 0, max = 50, step = 1),
        actionButton("draw_marker_heatmap", "Draw heatmap",
                     icon = icon("th"), class = "btn-primary")
      ),
      card(
        # `fill = FALSE` keeps this short explanatory card from being
        # stretched to fill the column height when its taller siblings
        # below (the marker tables and the heatmap card) push the
        # layout. Without it, bslib's default flex behavior leaves a
        # large empty band beneath the markdown blurb.
        fill = FALSE,
        card_header(icon("circle-info"), " About markers"),
        card_body(
          markdown(
            "Wilcoxon-based differential expression between the adverse /
            favorable tails and the rest of the data (via `presto` when
            installed, base R otherwise; Seurat's `FindMarkers` is used for
            Seurat inputs in cohort-wide mode).

            **Cell-type-specific** contrasts each phenotype tail × cell type
            against other phenotype groups in the *same* cell type by
            default."
          )
        )
      ),
      navset_tab(
        nav_panel(
          "Adverse markers",
          card_body(DTOutput("adverse_markers_tbl"))
        ),
        nav_panel(
          "Favorable markers",
          card_body(DTOutput("favorable_markers_tbl"))
        )
      ),
      card(
        card_header(tags$strong("Marker-gene heatmap")),
        card_body(
          helpText(
            "Heatmap of the top markers from the adverse and favorable tails. ",
            "When markers were computed with cell-type-specific scope, the ",
            "heatmap shows one column slice per phenotype-group × cell-type ",
            "block; otherwise it shows a cohort-wide view. Use the controls ",
            "in the sidebar to set the number of genes per block, the number ",
            "of labels, and to redraw."
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
    reference = NULL,          # "precog" / "tcga" / ... OR data.frame
    reference_label = "(none)",
    derive_bulk = NULL,
    derive_phen = NULL,
    scores = NULL,
    groups = NULL,
    markers = NULL
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
  umap_upload_pick     <- phenomapr_file_pick("umap_upload",     input, output, session, shiny_file_roots,
                                              accept = c(".tsv", ".csv", ".txt", ".rds"))

  # ------------------------------------------------------------------------
  # 1. Expression upload
  # ------------------------------------------------------------------------
  observeEvent(expr_file_pick(), {
    pick <- req(expr_file_pick())
    phenomapr_busy_show("Reading expression file...", pick$name)
    on.exit(phenomapr_busy_hide(), add = TRUE)
    res <- tryCatch(
      parse_expression_upload(pick$datapath, pick$name),
      error = function(e) {
        phenomapr_busy_hide()
        showNotification(paste0("Upload failed: ", conditionMessage(e)),
                         type = "error", duration = 8)
        NULL
      }
    )
    if (is.null(res)) return()
    state$expression <- res$object
    state$expr_summary <- res
    state$metadata <- extract_object_metadata(res$object)
    state$meta_columns <- if (!is.null(state$metadata)) colnames(state$metadata) else character(0)
    state$metadata_source <- if (!is.null(state$metadata)) "object" else "(none)"
    showNotification(paste(res$notes, collapse = " "), type = "message", duration = 5)
  })

  observeEvent(input$use_demo, {
    set.seed(7)
    n_genes <- 200L; n_cells <- 60L
    m <- matrix(rpois(n_genes * n_cells, lambda = 3), nrow = n_genes, ncol = n_cells)
    rownames(m) <- c("TP53","MYC","EGFR","BRCA1","CDKN2A","KRAS","PIK3CA","PTEN","MKI67","CD68",
                     paste0("GENE", seq_len(n_genes - 10L)))
    colnames(m) <- paste0("Cell_", seq_len(n_cells))
    state$expression <- m
    state$expr_summary <- summarize_expression_object(m)
    state$expr_summary$notes <- "Demo matrix (200 genes × 60 cells) loaded."
    state$metadata <- data.frame(
      .cell_id = colnames(m),
      cell_type = sample(c("Acinar","Ductal","Macrophage","CD8T","Fibroblast"), n_cells, TRUE),
      Source = sample(c("Tumor","Normal"), n_cells, TRUE),
      stringsAsFactors = FALSE
    )
    state$meta_columns <- colnames(state$metadata)
    state$metadata_source <- "demo"
    showNotification("Demo matrix loaded.", type = "message", duration = 4)
  })

  # Optional metadata upload
  observeEvent(meta_file_pick(), {
    pick <- req(meta_file_pick())
    phenomapr_busy_show("Loading cell metadata...", pick$name)
    on.exit(phenomapr_busy_hide(), add = TRUE)
    md <- tryCatch(
      parse_metadata_upload(pick$datapath, pick$name),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8)
        NULL
      }
    )
    if (is.null(md)) return()
    state$metadata <- md
    state$meta_columns <- colnames(md)
    state$metadata_source <- "upload"
    showNotification("Metadata loaded.", type = "message", duration = 4)
  })

  output$meta_columns_available <- reactive({ length(state$meta_columns) > 0 })
  outputOptions(output, "meta_columns_available", suspendWhenHidden = FALSE)

  # Metadata-upload panel (a <details> block above the cell-ID / cell-type /
  # source dropdowns). When an expression file has been loaded but no
  # metadata was detected, the panel is rendered with `open` so the user is
  # prompted to attach a tabular metadata file. Once metadata is available
  # — from the object, a demo, or an uploaded file — the panel collapses to
  # keep the sidebar tidy. The user can still toggle it manually.
  output$metadata_upload_panel <- renderUI({
    has_expr <- !is.null(state$expression)
    no_meta  <- is.null(state$metadata) || !length(state$meta_columns)
    auto_open <- has_expr && no_meta

    summary_label <- if (auto_open) {
      "No metadata detected — upload a tabular metadata file"
    } else {
      "Override / supplement with a tabular metadata file"
    }

    body <- tagList(
      helpText(
        "Optional — only needed for matrix uploads, or to attach extra ",
        "annotations on top of an object's built-in metadata. A cell-ID ",
        "column must match the expression matrix columns."
      ),
      phenomapr_file_input(
        "meta_file",
        label = NULL,
        accept = c(".rds", ".tsv", ".csv", ".txt"),
        width = "100%"
      )
    )

    if (auto_open) {
      tags$details(
        class = "metadata-upload-details metadata-upload-details-prompt",
        open = "open",
        tags$summary(
          tags$span(class = "mu-badge", icon("upload")),
          tags$strong(summary_label)
        ),
        body
      )
    } else {
      tags$details(
        class = "metadata-upload-details",
        tags$summary(summary_label),
        body
      )
    }
  })

  # Status banner above the metadata dropdowns: explains where the metadata
  # currently in use came from (object / upload / demo / nothing yet).
  output$metadata_status <- renderUI({
    if (is.null(state$metadata) || !length(state$meta_columns)) {
      return(tags$div(
        class = "metadata-status metadata-status-empty",
        tags$em("No metadata loaded yet — load a Seurat / SCE / AnnData object, ",
                "or upload a tabular metadata file below.")
      ))
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
        "^(cell[._\\s]?type|celltype|cell_type_[a-z]+)$",
        "annotated_celltype|predicted\\.celltype|leiden_celltype|seurat_clusters_annot",
        "^Annotation$|^annotation$",
        "^type$"
      ),
      exclude = c("orig.ident", "feature_type")
    )

    # Source / group: cohort / sample / condition / donor / patient / tissue.
    src_default <- .pick_default(
      cols,
      patterns = c(
        "^source$|^Source$|^sample_source$",
        "^condition$|^Condition$|^treatment$",
        "^donor$|^patient$|^subject$",
        "^tissue$|^region$|^organ$",
        "^sample$|^Sample$|^orig\\.ident$"
      )
    )

    updateSelectInput(session, "meta_cell_id_col",   choices = cols, selected = cell_default)
    updateSelectInput(session, "meta_cell_type_col", choices = cols, selected = ct_default)
    updateSelectInput(session, "meta_source_col",    choices = cols, selected = src_default)
    updateSelectInput(session, "pseudobulk_group_by",
                      choices = setdiff(cols, "(none)"),
                      selected = state$meta_columns[1L])
  })

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
      source_col = input$meta_source_col
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
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1),
            legend.position = "none")
    if (!is.null(pal)) p <- p + scale_fill_manual(values = pal)
    p
  })

  output$source_count_plot <- renderPlot({
    md <- state$metadata
    src_col <- input$meta_source_col
    req(!is.null(md), nzchar(src_col), src_col != "(none)", src_col %in% colnames(md))
    df <- as.data.frame(table(md[[src_col]], useNA = "no"))
    colnames(df) <- c("source", "n")
    df <- df[order(-df$n), ]
    df$source <- factor(df$source, levels = df$source)
    ggplot(df, aes(x = source, y = n, fill = source)) +
      .geom_rounded_col() +
      scale_fill_phenomapr_d() +
      labs(x = NULL, y = "Cells", fill = "Source") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1),
            legend.position = "none")
  })

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
    if (!nrow(df)) return(NULL)
    df_count <- df %>% dplyr::count(cell_type, source)
    ggplot(df_count, aes(x = cell_type, y = n, fill = source)) +
      .geom_rounded_stack() +
      scale_fill_phenomapr_d() +
      labs(x = NULL, y = "Cells", fill = "Source") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
  })

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
    datatable(df, rownames = FALSE,
              options = list(pageLength = 8, dom = "tip"))
  })

  output$expr_axes_warning <- renderUI({
    s <- state$expr_summary
    if (is.null(s) || s$kind != "matrix") return(NULL)
    msgs <- character(0)
    if (!is.null(s$gene_ids) && length(s$gene_ids)) {
      if (mean(grepl("^ENSG\\d", s$gene_ids)) > 0.5) {
        msgs <- c(msgs, paste(
          "Gene names look like ENSG IDs; convert to HUGO symbols before scoring."
        ))
      }
    }
    if (!is.na(s$n_samples) && !is.na(s$n_genes) && s$n_samples > 10 * s$n_genes) {
      msgs <- c(msgs, paste(
        "Matrix has many more columns than rows — the loader will transpose to genes × samples."
      ))
    }
    if (!length(msgs)) return(NULL)
    div(class = "alert alert-warning",
        tags$strong("Heads up:"), tags$ul(lapply(msgs, tags$li)))
  })

  # ------------------------------------------------------------------------
  # 2. Reference
  # ------------------------------------------------------------------------
  # Populate cancer_type dropdown when reference changes.
  observeEvent(input$reference_choice, {
    if (input$reference_choice == "_custom") return()
    cts <- get_cancer_types(input$reference_choice)
    updateSelectInput(session, "cancer_type", choices = cts, selected = cts[1L])
  })

  # Custom: upload file
  observeEvent(custom_ref_file_pick(), {
    pick <- req(custom_ref_file_pick())
    phenomapr_busy_show("Loading custom signature...", pick$name)
    on.exit(phenomapr_busy_hide(), add = TRUE)
    ref <- tryCatch(
      parse_reference_upload(pick$datapath, pick$name),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8)
        NULL
      }
    )
    if (is.null(ref)) return()
    state$reference <- ref
    state$reference_label <- paste0("Custom (", pick$name, ")")
    showNotification(
      sprintf("Custom reference loaded (%d genes).", nrow(ref)),
      type = "message", duration = 4
    )
  })

  # Custom: derive from bulk + phenotype
  observeEvent(derive_bulk_file_pick(), {
    pick <- req(derive_bulk_file_pick())
    phenomapr_busy_show("Loading bulk expression...", pick$name)
    on.exit(phenomapr_busy_hide(), add = TRUE)
    res <- tryCatch(
      parse_expression_upload(pick$datapath, pick$name),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8); NULL
      }
    )
    if (!is.null(res)) {
      state$derive_bulk <- res$object
      showNotification(
        sprintf("Bulk expression: %s genes × %s samples.",
                .fmt_int(res$n_genes %||% 0L), .fmt_int(res$n_samples %||% 0L)),
        type = "message", duration = 5
      )
    }
  })

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
    pick <- req(derive_phen_file_pick())
    phenomapr_busy_show("Loading phenotype table...", pick$name)
    on.exit(phenomapr_busy_hide(), add = TRUE)
    df <- tryCatch(
      parse_metadata_upload(pick$datapath, pick$name),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8); NULL
      }
    )
    if (is.null(df)) return()
    state$derive_phen <- df
    cols <- colnames(df)
    updateSelectInput(session, "derive_id_col", choices = cols, selected = cols[1L])
    updateSelectInput(session, "derive_pheno_col", choices = cols,
                      selected = if (length(cols) > 1L) cols[2L] else cols[1L])
    updateSelectInput(session, "derive_time_col", choices = cols)
    updateSelectInput(session, "derive_event_col", choices = cols)
  })

  observeEvent(input$derive_run, {
    req(state$derive_bulk, state$derive_phen)
    phenomapr_busy_show(
      "Deriving phenotype signature...",
      sprintf("Bulk + phenotype | %s outcome", input$derive_type %||% "binary")
    )
    on.exit(phenomapr_busy_hide(), add = TRUE)
    bin_pos <- if (input$derive_binary_positive %in% c("first", "second")) {
      input$derive_binary_positive
    } else "second"
    ref <- tryCatch(
      PhenoMapR::derive_reference_from_bulk(
        bulk_expression  = state$derive_bulk,
        phenotype        = state$derive_phen,
        sample_id_column = input$derive_id_col,
        phenotype_column = if (input$derive_type == "survival") NULL else input$derive_pheno_col,
        phenotype_type   = input$derive_type,
        survival_time    = if (input$derive_type == "survival") input$derive_time_col else NULL,
        survival_event   = if (input$derive_type == "survival") input$derive_event_col else NULL,
        normalize        = isTRUE(input$derive_normalize),
        hugo_species     = input$derive_hugo_species %||% "human",
        binary_positive_reference = bin_pos,
        verbose          = TRUE
      ),
      error = function(e) {
        showNotification(paste0("derive_reference_from_bulk failed: ",
                                conditionMessage(e)),
                         type = "error", duration = 10); NULL
      }
    )
    if (!is.null(ref)) {
      state$reference <- ref
      state$reference_label <- sprintf(
        "Custom (derived from bulk + phenotype, type = %s)",
        attr(ref, "phenotype_type") %||% input$derive_type
      )
      showNotification(
        sprintf("Reference derived (%s genes).", .fmt_int(nrow(ref))),
        type = "message", duration = 5
      )
    }
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
    datatable(df, rownames = FALSE,
              options = list(pageLength = 10, dom = "tip", scrollX = TRUE))
  })

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
    datatable(df, rownames = FALSE,
              options = list(pageLength = 10, dom = "tip", scrollX = TRUE))
  })

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
    if (is.character(ref)) {
      tagList(
        tags$p(tags$strong("Built-in: "), ref),
        tags$p(tags$strong("Cancer type: "), input$cancer_type %||% "(none)"),
        tags$p(tags$strong("|z| cutoff: "), input$z_score_cutoff),
        tags$p(tags$strong("Direction sign: "), input$reference_sign)
      )
    } else {
      tagList(
        tags$p(tags$strong("Source: "), state$reference_label),
        tags$p(tags$strong("Genes: "), .fmt_int(nrow(ref))),
        tags$p(tags$strong("z column: "), colnames(ref)[1L])
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
    datatable(cov, rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })

  output$reference_signature_plot <- renderPlot({
    req(state$reference)
    if (is.data.frame(state$reference) || is.matrix(state$reference)) {
      if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        showNotification("Install ComplexHeatmap for the signature heatmap.",
                         type = "warning", duration = 5)
        return(NULL)
      }
      PhenoMapR::plot_reference_signature(state$reference)
    } else {
      # Built-in: gate on `cancer_type` actually being a valid cancer type
      # for the currently selected reference. When the user switches
      # reference (e.g. PRECOG -> TCGA) Shiny fires `input$reference_choice`
      # immediately but `input$cancer_type` only updates asynchronously
      # via `updateSelectInput()`. Without this guard, the renderer
      # transiently sees a stale (reference, cancer_type) pair like
      # ("tcga", "Adrenocortical") and `get_top_prognostic_genes()` throws.
      valid_cts <- get_cancer_types(state$reference)
      req(input$cancer_type, input$cancer_type %in% valid_cts)
      tg <- PhenoMapR::get_top_prognostic_genes(
        reference = state$reference, cancer_type = input$cancer_type,
        n = 500L, direction = "both"
      )
      ggplot(tg, aes(x = z_score)) +
        .geom_rounded_histogram(bins = 40, fill = "#2A9D8F", color = "white") +
        labs(x = "Reference z-score", y = "Count",
             title = sprintf("Top %d genes by |z| (%s · %s)",
                             nrow(tg), state$reference, input$cancer_type %||% "")) +
        theme_minimal(base_size = 13)
    }
  })

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
    on.exit(phenomapr_busy_hide(), add = TRUE)
    sign <- tryCatch(as.integer(input$reference_sign), error = function(e) 1L)
    scores <- tryCatch(
      run_phenomap_with_progress(
        expression = state$expression,
        reference = state$reference,
        cancer_type = input$cancer_type,
        z_score_cutoff = input$z_score_cutoff,
        pseudobulk = isTRUE(input$pseudobulk),
        group_by = input$pseudobulk_group_by,
        assay = input$score_assay,
        slot = input$score_slot,
        reference_sign = sign
      ),
      error = function(e) {
        showNotification(paste0("PhenoMap failed: ", conditionMessage(e)),
                         type = "error", duration = 10); NULL
      }
    )
    if (is.null(scores)) return()
    state$scores <- scores
    state$groups <- NULL  # invalidate downstream
    state$markers <- NULL
    nm <- colnames(scores)
    updateSelectInput(session, "groups_score_col", choices = nm, selected = nm[1L])
    showNotification(
      sprintf("Scored %s (%.1fs).",
              .fmt_n_units(nrow(scores), "row"),
              attr(scores, "elapsed_s") %||% NA_real_),
      type = "message", duration = 5
    )
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
      df <- data.frame(score = as.numeric(s[[cn]]))
      ttl <- sprintf("PhenoMapR score distribution (%s)", cn)
    }
    PhenoMapR::plot_score_distribution(df, score_column = "score", main = ttl)
  })

  # Score table. Augments the raw PhenoMapR score column(s) with a
  # matching `scaled_<col>` column (mean 0, sd 1) for each numeric score
  # column. The scaled value makes it easy to compare cells across
  # references that produce different raw magnitudes.
  output$score_table <- renderDT({
    req(state$scores)
    s <- state$scores
    s$cell_id <- rownames(s)
    score_cols <- setdiff(colnames(s), "cell_id")
    numeric_cols <- score_cols[vapply(s[score_cols], is.numeric, logical(1))]
    for (col in numeric_cols) {
      v <- as.numeric(scale(s[[col]]))
      s[[paste0("scaled_", col)]] <- v
    }
    ordered_cols <- c("cell_id", unlist(lapply(score_cols, function(col) {
      if (col %in% numeric_cols) c(col, paste0("scaled_", col)) else col
    }), use.names = FALSE))
    s <- s[, ordered_cols, drop = FALSE]
    datatable(
      s, rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE)
    ) |>
      DT::formatRound(intersect(ordered_cols,
                                c(numeric_cols, paste0("scaled_", numeric_cols))),
                      digits = 3)
  })

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
    if (!nrow(d)) return(NULL)
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
    base <- ggplot(d, aes(x = rank, y = score)) +
      labs(x = "Rank by PhenoMapR score", y = "PhenoMapR score") +
      theme_minimal(base_size = 13)
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
      return(p)
    }
    if (color_by == "source" && "source" %in% colnames(d)) {
      # Reuse the PhenoMapR brand palette so the "Source" coloring in this
      # plot matches the "Cells per source / group" and "Cell type ×
      # source composition" plots on the Data tab. Sort levels by their
      # appearance order in the input metadata so the color mapping is
      # stable across refreshes.
      d$source <- factor(as.character(d$source),
                         levels = unique(as.character(d$source)))
      return(
        base +
          geom_point(aes(color = source), size = 0.7, alpha = 0.85) +
          scale_color_phenomapr_d() +
          labs(color = "Source") +
          legend_pt_override
      )
    }
    base +
      geom_point(aes(color = score), size = 0.7, alpha = 0.85) +
      scale_color_gradient2(
        low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
        midpoint = 0,
        name = "Score"
      )
  })

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
    if (!nrow(d)) return(NULL)

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

    # Compute Wilcoxon p-values:
    #   - Source mapped: per-cell-type contrast between the two sources;
    #     only significant brackets are returned (significant_only = TRUE).
    #   - Source not mapped: pairwise Wilcoxon between cell types,
    #     brackets stacked above the boxes for the significant pairs only.
    # Both helpers honour the median-based `cell_levels` ordering so the
    # bracket xmin/xmax positions match the plot's x-axis.
    if (has_source) {
      pval_df <- celltype_source_pvalues(
        dl, "Score", "Cell_type", "Source",
        cell_levels = cell_levels, significant_only = TRUE
      )
    } else {
      pval_df <- celltype_pairwise_pvalues(
        dl, "Score", "Cell_type", cell_levels = cell_levels
      )
    }

    dl$Cell_type <- factor(dl$Cell_type, levels = cell_levels)

    pal_cells <- tryCatch(
      PhenoMapR::get_celltype_palette(cell_levels),
      error = function(e) NULL
    )

    if (has_source) {
      # Source levels are colored with the PhenoMapR brand palette so the
      # Source coloring on this plot lines up 1:1 with the source-keyed
      # plots elsewhere ("Cells per source / group", "Cell type × source
      # composition", score rank plot, UMAP). We freeze the level order
      # via factor() so the palette → level mapping is stable.
      src_levels <- sort(unique(dl$Source))
      dl$Source <- factor(dl$Source, levels = src_levels)
      src_pal <- setNames(pm_brand_palette(length(src_levels)), src_levels)
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
          fill  = guide_legend(order = 1),
          color = guide_legend(order = 2)
        )
      title <- "PhenoMapR Score distribution by Source and Cell Type"
    } else {
      p <- ggplot(dl, aes(x = Cell_type, y = scale(Score), fill = Cell_type)) +
        .geom_rounded_boxplot(
          outlier.alpha    = 0.5,
          median.linewidth = 0.5,
          outlier.shape    = 21
        )
      title <- "PhenoMapR Score distribution by Cell Type"
    }

    p <- p +
      labs(x = NULL, y = "PhenoMapR score (scaled)", title = title) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        plot.title      = element_text(hjust = 0.5)
      )
    if (!is.null(pal_cells)) {
      p <- p + scale_fill_manual(values = pal_cells, name = "Cell type")
    } else {
      p <- p + labs(fill = "Cell type")
    }

    if (!is.null(pval_df) && nrow(pval_df) > 0) {
      use_signif <- requireNamespace("ggsignif", quietly = TRUE)
      if (use_signif) {
        # Suppress benign "Ignoring unknown aesthetics" warnings from
        # geom_signif manual mode — cosmetic, clutters the Shiny console.
        suppressWarnings(
          p <- p + ggsignif::geom_signif(
            data        = pval_df,
            aes(xmin = xmin, xmax = xmax,
                annotations = label, y_position = y_pos),
            manual      = TRUE,
            tip_length  = 0.02,
            textsize    = 4,
            color       = "black",
            inherit.aes = FALSE
          )
        )
      } else {
        # ggsignif fallback: hand-drawn brackets.
        seg_df <- data.frame(
          x = pval_df$xmin, xend = pval_df$xmax,
          y = pval_df$y_pos, yend = pval_df$y_pos
        )
        tick_df <- data.frame(
          x = c(pval_df$xmin, pval_df$xmax),
          xend = c(pval_df$xmin, pval_df$xmax),
          y    = c(pval_df$y_pos, pval_df$y_pos),
          yend = c(pval_df$y_pos - 0.1, pval_df$y_pos - 0.1)
        )
        text_df <- data.frame(
          x = (pval_df$xmin + pval_df$xmax) / 2,
          y = pval_df$y_pos + 0.15,
          label = pval_df$label
        )
        p <- p +
          geom_segment(data = seg_df,
                       aes(x = x, xend = xend, y = y, yend = yend),
                       inherit.aes = FALSE, color = "black") +
          geom_segment(data = tick_df,
                       aes(x = x, xend = xend, y = y, yend = yend),
                       inherit.aes = FALSE, color = "black") +
          geom_text(data = text_df,
                    aes(x = x, y = y, label = label),
                    inherit.aes = FALSE, size = 4)
      }
    }

    p
  })

  # ------------------------------------------------------------------------
  # UMAP / embedding tab
  # ------------------------------------------------------------------------
  # Custom embedding upload (overrides object-resident reductions).
  uploaded_embedding <- reactiveVal(NULL)

  observeEvent(umap_upload_pick(), {
    pick <- req(umap_upload_pick())
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
  })

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

  output$umap_plot <- renderPlot({
    emb <- current_embedding()
    req(emb)
    color_by <- input$umap_color_by %||% "score"
    ct <- cell_table()
    df <- emb
    df <- dplyr::left_join(df, ct %||% data.frame(cell_id = character(0)), by = "cell_id")
    pt_size <- input$umap_point_size %||% 0.8
    pt_alpha <- input$umap_point_alpha %||% 0.75
    # Detect spatial-frame embeddings (set by `extract_embedding()` when
    # the user picks the synthetic "spatial" reduction). Spatial plots
    # need (a) equal aspect so tissue isn't squashed and (b) a reversed
    # y-axis since image-space coordinates have origin at top-left.
    is_spatial <- isTRUE(any(df$is_spatial))

    base <- ggplot(df, aes(x = dim1, y = dim2)) +
      labs(x = unique(df$dim1_name)[1L] %||% "dim1",
           y = unique(df$dim2_name)[1L] %||% "dim2") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
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
      # Pass the *updated* df (with `score_to_plot`) explicitly to
      # geom_point so the layer sees the new column — ggplot2 4.0.0
      # deprecated `%+%` for swapping plot-level data after-the-fact, so
      # supplying the data at the layer level keeps things forward-
      # compatible without rebuilding `base`.
      p <- base +
        geom_point(data = df, aes(color = score_to_plot),
                   size = pt_size, alpha = pt_alpha) +
        scale_color_gradient2(
          low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
          midpoint = 0, limits = c(-lim, lim),
          oob = scales::squish, name = legend_name
        )
    } else if (color_by == "cell_type" && "cell_type" %in% colnames(df)) {
      pal <- tryCatch(
        PhenoMapR::get_celltype_palette(as.character(unique(df$cell_type))),
        error = function(e) NULL
      )
      p <- base +
        geom_point(aes(color = cell_type), size = pt_size, alpha = pt_alpha) +
        labs(color = "Cell type")
      if (!is.null(pal)) p <- p + scale_color_manual(values = pal)
    } else if (color_by == "source" && "source" %in% colnames(df)) {
      # Apply the PhenoMapR brand palette so the Source coloring here
      # matches the Source coloring on the Data tab and Score tab plots.
      df$source <- factor(as.character(df$source),
                          levels = unique(as.character(df$source)))
      p <- base +
        geom_point(aes(color = source), size = pt_size, alpha = pt_alpha) +
        scale_color_phenomapr_d() +
        labs(color = "Source") +
        ggplot2::guides(color = ggplot2::guide_legend(
          override.aes = list(size = 4, alpha = 1)
        ))
    } else if (color_by == "group" && "group" %in% colnames(df)) {
      p <- base +
        geom_point(aes(color = group), size = pt_size, alpha = pt_alpha) +
        scale_color_manual(values = c(
          "Most Adverse"   = "#B2182B",
          "Most Favorable" = "#2166AC",
          "Other"          = "#BBBBBB"
        ), na.value = "#E0E0E0", name = "Group")
    } else {
      p <- base + geom_point(size = pt_size, alpha = pt_alpha, color = "#555555")
    }

    if (isTRUE(input$umap_facet_source) && "source" %in% colnames(df)) {
      p <- p + facet_wrap(~ source)
    }
    p
  })

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
        percentile = input$percentile,
        score_columns = if (nzchar(input$groups_score_col %||% "")) input$groups_score_col else NULL
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
    grp_col <- grep("^phenotype_group_", colnames(g), value = TRUE)[1L]
    if (is.na(grp_col)) return(tags$p(tags$em("No phenotype group column found.")))
    tbl <- as.data.frame(table(g[[grp_col]], useNA = "ifany"))
    colnames(tbl) <- c("Phenotype group", "Cells")
    tbl$Cells <- .fmt_int(tbl$Cells)
    tagList(
      tags$p(tags$strong("Group column: "), grp_col),
      tag_table(tbl)
    )
  })

  output$group_by_celltype_plot <- renderPlot({
    req(state$groups)
    g <- state$groups
    grp_col <- grep("^phenotype_group_", colnames(g), value = TRUE)[1L]
    if (is.na(grp_col)) return(NULL)
    ct_col <- input$meta_cell_type_col
    if (is.null(ct_col) || ct_col == "(none)" || !ct_col %in% colnames(g)) {
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
    ggplot(df_count, aes(x = cell_type, y = frac, fill = group)) +
      .geom_rounded_stack(color = NA) +
      scale_y_continuous(labels = scales::percent_format()) +
      scale_fill_manual(values = c(
        "Most Adverse"    = "#B2182B",
        "Most Favorable" = "#2166AC",
        "Other"           = "#BBBBBB"
      )) +
      labs(x = "Cell type", y = "Fraction of cells", fill = "Group") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
  })

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
    on.exit(phenomapr_busy_hide(), add = TRUE)

    grp_col <- grep("^phenotype_group_", colnames(state$groups), value = TRUE)[1L]
    if (is.na(grp_col)) {
      showNotification("Could not find a phenotype_group_* column.",
                       type = "error", duration = 6); return()
    }
    ct_col_use <- input$meta_cell_type_col
    if (ct_col_use == "(none)" || !ct_col_use %in% colnames(state$groups)) {
      ct_col_use <- NULL
    }

    markers <- tryCatch(
      PhenoMapR::find_phenotype_markers(
        expression = state$expression,
        group_labels = state$groups,
        group_column = grp_col,
        cell_id_column = "cell_id",
        marker_scope = input$marker_scope,
        cell_type_column = ct_col_use,
        min.pct = input$marker_min_pct,
        logfc.threshold = input$marker_logfc,
        pval_threshold = input$marker_pval,
        max_cells_per_ident = input$marker_maxcells,
        verbose = TRUE
      ),
      error = function(e) {
        showNotification(paste0("find_phenotype_markers failed: ",
                                conditionMessage(e)),
                         type = "error", duration = 10); NULL
      }
    )
    if (is.null(markers)) return()
    state$markers <- markers
    showNotification(sprintf(
      "Found %s adverse + %s favorable markers.",
      .fmt_int(nrow(markers$adverse_markers)),
      .fmt_int(nrow(markers$favorable_markers))
    ), type = "message", duration = 5)
  })

  output$adverse_markers_tbl <- renderDT({
    req(state$markers)
    df <- state$markers$adverse_markers
    if (!nrow(df)) return(datatable(df, options = list(dom = "t")))
    datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })
  output$favorable_markers_tbl <- renderDT({
    req(state$markers)
    df <- state$markers$favorable_markers
    if (!nrow(df)) return(datatable(df, options = list(dom = "t")))
    datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })

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
    on.exit(phenomapr_busy_hide(), add = TRUE)

    # For AnnData inputs, the marker heatmap only needs the marker genes —
    # we pass that list straight to extract_expression_matrix() so that
    # subsetting happens on the Python side and we don't materialise the
    # full expression matrix in R.
    marker_genes <- unique(c(
      if (!is.null(state$markers$adverse_markers))   state$markers$adverse_markers$gene,
      if (!is.null(state$markers$favorable_markers)) state$markers$favorable_markers$gene
    ))
    marker_genes <- marker_genes[nzchar(marker_genes)]

    expr_mat <- tryCatch(
      extract_expression_matrix(
        state$expression,
        assay = if (nzchar(input$score_assay %||% "")) input$score_assay else NULL,
        slot = input$score_slot %||% "data",
        gene_subset = if (length(marker_genes)) marker_genes else NULL
      ),
      error = function(e) {
        showNotification(paste0("Could not extract expression matrix: ",
                                conditionMessage(e)),
                         type = "error", duration = 8); NULL
      }
    )
    if (is.null(expr_mat)) return()

    # Build per-cell metadata: cell_id, group, score, cell_type
    grp_col <- grep("^phenotype_group_", colnames(state$groups), value = TRUE)[1L]
    if (is.na(grp_col)) {
      showNotification("No phenotype_group_* column on the groups table — re-tag groups first.",
                       type = "error", duration = 6); return()
    }
    score_name <- colnames(state$scores)[1L]
    scores_df <- data.frame(
      cell_id = rownames(state$scores),
      .score = as.numeric(state$scores[[1L]]),
      stringsAsFactors = FALSE
    )
    colnames(scores_df)[2L] <- score_name
    meta_df <- data.frame(
      cell_id = as.character(state$groups$cell_id),
      stringsAsFactors = FALSE
    )
    meta_df$.group <- state$groups[[grp_col]]
    colnames(meta_df)[2L] <- grp_col
    meta_df <- dplyr::left_join(meta_df, scores_df, by = "cell_id")

    ct_col_use <- input$meta_cell_type_col
    has_ct <- !is.null(ct_col_use) && nzchar(ct_col_use) && ct_col_use != "(none)" &&
      !is.null(state$metadata) && ct_col_use %in% colnames(state$metadata)
    if (has_ct) {
      md <- state$metadata
      id_col <- input$meta_cell_id_col
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

    # Decide which heatmap to draw: cell-type-specific when both the markers
    # were computed that way AND a cell-type column is mapped.
    is_ct <- input$marker_scope == "cell_type_specific" && has_ct
    heatmap_type <- if (is_ct) "cell_type_specific" else "global"

    marker_heatmap_args(list(
      markers = state$markers,
      expr_mat = expr_mat,
      meta = meta_df,
      cell_id_col = "cell_id",
      group_col = grp_col,
      score_col = score_name,
      celltype_col = if (has_ct) ct_col_use else NULL,
      heatmap_type = heatmap_type,
      top_n_markers = input$hm_top_n %||% 20L,
      n_mark_labels = input$hm_n_labels %||% 5L,
      use_raster = FALSE
    ))
    showNotification("Heatmap ready -- drawing...", type = "message", duration = 3)
  })

  output$marker_heatmap <- renderImage(
    {
      args <- marker_heatmap_args(); req(args)

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
      return(datatable(
        data.frame(
          note = paste0(
            "Reference diagnostics are only available for the built-in ",
            "signatures (PRECOG, TCGA, Pediatric PRECOG, ICI PRECOG). ",
            "Switch the 'Phenotype signature source' above to view top ",
            "prognostic genes."
          )
        ),
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
        n           = input$top_genes_n,
        direction   = input$top_genes_dir
      ),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8)
        NULL
      }
    )
    req(!is.null(tg))
    datatable(tg, rownames = FALSE,
              options = list(pageLength = 15, scrollX = TRUE))
  })

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
