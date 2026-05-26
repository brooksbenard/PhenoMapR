# Launch the PhenoMapR Shiny app

Starts a local Shiny session that walks users through the full PhenoMapR
workflow: upload data, choose a reference, score, tag phenotype groups,
and find marker genes. The app is a thin UI over the same exported
functions you would call from a script (`PhenoMap`,
`define_phenotype_groups`, `find_phenotype_markers`,
`derive_reference_from_bulk`, etc.).

## Usage

``` r
run_app(
  host = "127.0.0.1",
  port = NULL,
  launch.browser = interactive(),
  max_upload_size_mb = Inf,
  ...
)
```

## Arguments

- host:

  Character. Network interface to bind to (default `"127.0.0.1"`; use
  `"0.0.0.0"` to expose on all interfaces for remote use).

- port:

  Integer or `NULL` (random free port).

- launch.browser:

  Logical. Open the app in the default browser (`TRUE`); on headless
  servers set `FALSE`.

- max_upload_size_mb:

  Numeric. Per-file upload limit in megabytes. Defaults to `Inf` (no
  limit), so users on a workstation or HPC node with plenty of memory
  can upload very large `.rds` / `.h5ad` objects without hitting Shiny's
  default 5 MB cap. Set to a finite number only if you want to cap
  uploads (e.g. on a shared web server).

- ...:

  Additional arguments passed to
  [`runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Value

Invisibly, the return value of
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Details

On a remote workstation or HPC node, run the app with
`run_app(host = "0.0.0.0", port = 3838, launch.browser = FALSE)` and
point a browser at `http://<server>:<port>`. No data leaves the machine
the app runs on.

## Examples

``` r
if (FALSE) { # \dontrun{
# Local laptop use:
PhenoMapR::run_app()

# Remote workstation, listen on all interfaces:
PhenoMapR::run_app(host = "0.0.0.0", port = 3838, launch.browser = FALSE)
} # }
```
