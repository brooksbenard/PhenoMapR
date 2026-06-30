# Shared helpers for PhenoMapR real-data integration stress tests.

.integration_state <- new.env(parent = emptyenv())
.integration_state$results <- list()
.integration_state$config <- NULL
.integration_state$options <- list(
  update_fixtures = FALSE,
  stress = FALSE,
  full = FALSE,
  download = TRUE,
  skip_smokes = FALSE
)

`%||%` <- function(x, y) if (is.null(x)) y else x

integration_init <- function(options = list(), root = NULL) {
  if (is.null(root)) root <- integration_repo_root()
  .integration_state$config <- integration_load_config(root)
  .integration_state$results <- list()
  opts <- modifyList(.integration_state$options, options)
  .integration_state$options <- opts
  dir.create(.integration_state$config$fixtures, recursive = TRUE, showWarnings = FALSE)
  dir.create(.integration_state$config$reports, recursive = TRUE, showWarnings = FALSE)
  invisible(.integration_state$config)
}

integration_get_option <- function(name) {
  .integration_state$options[[name]]
}

integration_fixture_path <- function(name) {
  file.path(.integration_state$config$fixtures, name)
}

integration_record <- function(id, label, status, details = list(), elapsed_sec = NA_real_) {
  entry <- list(
    id = id,
    label = label,
    status = status,
    details = details,
    elapsed_sec = elapsed_sec,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  .integration_state$results[[length(.integration_state$results) + 1L]] <- entry
  prefix <- switch(status,
    pass = "PASS",
    fail = "FAIL",
    skip = "SKIP",
    status
  )
  msg <- sprintf("[%s] %s: %s", prefix, id, label)
  if (status == "skip" && !is.null(details$reason)) {
    msg <- paste0(msg, " (", details$reason, ")")
  }
  if (status == "fail" && !is.null(details$error)) {
    msg <- paste0(msg, " — ", details$error)
  }
  message(msg)
  invisible(entry)
}

integration_assert <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

integration_run_test <- function(id, label, expr, skip_if = NULL) {
  if (!is.null(skip_if)) {
    reason <- skip_if
    if (is.character(skip_if) && length(skip_if) == 1L) reason <- skip_if
    if (isTRUE(skip_if)) reason <- "precondition not met"
    if (!isFALSE(skip_if)) {
      return(integration_record(id, label, "skip", list(reason = as.character(reason))))
    }
  }
  t0 <- proc.time()
  res <- tryCatch(
    {
      out <- force(expr)
      list(ok = TRUE, details = if (is.list(out)) out else list())
    },
    error = function(e) list(ok = FALSE, details = list(error = conditionMessage(e)))
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]
  if (isTRUE(res$ok)) {
    integration_record(id, label, "pass", res$details, elapsed_sec = elapsed)
  } else {
    integration_record(id, label, "fail", res$details, elapsed_sec = elapsed)
  }
}

integration_jaccard <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  if (!length(a) && !length(b)) return(1)
  if (!length(a) || !length(b)) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

integration_top_genes <- function(df, n = 20L) {
  if (is.null(df) || !nrow(df)) return(character())
  if ("p_val_adj" %in% names(df)) {
    df <- df[order(df$p_val_adj, df$p_val), , drop = FALSE]
  } else if ("p_val" %in% names(df)) {
    df <- df[order(df$p_val), , drop = FALSE]
  }
  head(as.character(df$gene), n)
}

integration_check_golden <- function(fixture_name, observed, threshold = 0.8) {
  path <- integration_fixture_path(fixture_name)
  if (isTRUE(integration_get_option("update_fixtures")) || !file.exists(path)) {
    saveRDS(observed, path)
    return(list(updated = TRUE, overlap = 1, path = path))
  }
  expected <- readRDS(path)
  overlap <- integration_jaccard(observed, expected)
  integration_assert(
    overlap >= threshold,
    sprintf(
      "Golden fixture %s overlap %.3f < %.3f",
      fixture_name, overlap, threshold
    )
  )
  list(updated = FALSE, overlap = overlap, path = path)
}

integration_load_package <- function(root) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(root, quiet = TRUE, helpers = FALSE, export_all = FALSE)
  } else {
    library(PhenoMapR)
  }
}

integration_source_shiny_helpers <- function() {
  shiny_dir <- .integration_state$config$shiny_dir
  helpers <- file.path(shiny_dir, "helpers.R")
  integration_assert(file.exists(helpers), paste("Missing", helpers))
  env <- new.env(parent = globalenv())
  sys.source(helpers, envir = env)
  if (exists("clear_shiny_demo_pool_cache", envir = env, inherits = FALSE)) {
    env$clear_shiny_demo_pool_cache()
  }
  env
}

integration_run_smokes <- function() {
  root <- .integration_state$config$root
  smoke_dir <- file.path(root, "inst", "benchmarks")
  smokes <- sort(list.files(smoke_dir, pattern = "_smoke\\.R$", full.names = TRUE))
  smokes <- smokes[!grepl("integration", smokes)]
  old <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = NA_character_)
  Sys.setenv(PHENOMAPR_SHINY_DIR = .integration_state$config$shiny_dir)
  on.exit({
    if (is.na(old)) Sys.unsetenv("PHENOMAPR_SHINY_DIR") else Sys.setenv(PHENOMAPR_SHINY_DIR = old)
  }, add = TRUE)
  for (s in smokes) {
    id <- paste0("smoke:", basename(s))
    integration_run_test(
      id,
      paste("UI smoke:", basename(s)),
      expr = {
        out <- system2("Rscript", c(s), stdout = TRUE, stderr = TRUE)
        code <- attr(out, "status") %||% 0L
        integration_assert(
          code == 0L,
          paste("Smoke failed:", paste(tail(out, 8L), collapse = "\n"))
        )
        list(lines = length(out))
      }
    )
  }
}

integration_summarize <- function() {
  res <- .integration_state$results
  status <- vapply(res, `[[`, character(1), "status")
  list(
    total = length(res),
    pass = sum(status == "pass"),
    fail = sum(status == "fail"),
    skip = sum(status == "skip"),
    results = res
  )
}

integration_write_json_report <- function(path = NULL) {
  if (is.null(path)) {
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    path <- file.path(.integration_state$config$reports, paste0("stress_test_", stamp, ".json"))
  }
  summary <- integration_summarize()
  json <- jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE)
  writeLines(as.character(json), path)
  message("Wrote JSON report: ", path)
  invisible(path)
}

integration_write_html_report <- function(path = NULL) {
  if (is.null(path)) {
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    path <- file.path(.integration_state$config$reports, paste0("stress_test_", stamp, ".html"))
  }
  summary <- integration_summarize()
  rows <- summary$results
  trs <- vapply(rows, function(r) {
    details <- r$details
    detail_txt <- if (length(details)) {
      paste(vapply(names(details), function(nm) {
        val <- details[[nm]]
        if (length(val) > 1L) val <- paste(val, collapse = ", ")
        sprintf("<li><strong>%s:</strong> %s</li>", nm, htmltools::htmlEscape(as.character(val)))
      }, character(1)), collapse = "")
    } else ""
    sprintf(
      "<tr class=\"%s\"><td>%s</td><td>%s</td><td>%s</td><td>%.2f</td><td><ul>%s</ul></td></tr>",
      r$status,
      htmltools::htmlEscape(r$id),
      htmltools::htmlEscape(r$label),
      htmltools::htmlEscape(r$status),
      r$elapsed_sec %||% NA_real_,
      detail_txt
    )
  }, character(1))
  html <- paste0(
    "<!DOCTYPE html><html><head><meta charset=\"utf-8\">",
    "<title>PhenoMapR integration stress test</title>",
    "<style>body{font-family:system-ui,sans-serif;margin:2rem}",
    "table{border-collapse:collapse;width:100%}",
    "td,th{border:1px solid #ccc;padding:.4rem .6rem;vertical-align:top}",
    "tr.pass td:nth-child(3){color:green}",
    "tr.fail td:nth-child(3){color:red}",
    "tr.skip td:nth-child(3){color:#888}</style></head><body>",
    "<h1>PhenoMapR integration stress test</h1>",
    "<p>Generated: ", format(Sys.time()), "</p>",
    "<p>Pass: ", summary$pass, " | Fail: ", summary$fail,
    " | Skip: ", summary$skip, " | Total: ", summary$total, "</p>",
    "<table><thead><tr><th>ID</th><th>Label</th><th>Status</th>",
    "<th>Seconds</th><th>Details</th></tr></thead><tbody>",
    paste(trs, collapse = "\n"),
    "</tbody></table></body></html>"
  )
  writeLines(html, path)
  message("Wrote HTML report: ", path)
  invisible(path)
}

integration_write_reports <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("jsonlite not installed; skipping JSON report")
  } else {
    integration_write_json_report()
  }
  if (!requireNamespace("htmltools", quietly = TRUE)) {
    warning("htmltools not installed; skipping HTML report")
  } else {
    integration_write_html_report()
  }
}
