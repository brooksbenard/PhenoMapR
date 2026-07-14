#!/usr/bin/env Rscript
## Smoke tests for the welcome-page cleanup:
##
##   The "Server / remote use" collapsible card -- the section that
##   advertised `PhenoMapR::run_app(host = "0.0.0.0", port = 3838)`
##   to remote-server users -- was removed from the Home tab at the
##   user's request. The deployment recipe still lives in the
##   package README and the ?run_app help page; the welcome screen
##   is now focused on the in-app workflow.
##
##   This script asserts that:
##     1. inst/shiny/app.R no longer contains the visible "Server /
##        remote use" tags$summary header (only the explanatory
##        comment about why the section was removed).
##     2. The companion `.welcome-server` / `.welcome-code` CSS
##        rules are gone from inst/shiny/www/styles.css.
##     3. The remote-launch command itself
##        (`PhenoMapR::run_app(host = "0.0.0.0", port = 3838)`) is no
##        longer rendered into the welcome card.

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

app_R      <- system.file("shiny", "app.R",            package = "PhenoMapR")
styles_css <- system.file("shiny", "www", "styles.css", package = "PhenoMapR")
dev_dir <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  app_R      <- file.path(dev_dir, "app.R")
  styles_css <- file.path(dev_dir, "www", "styles.css")
}
stopifnot(file.exists(app_R), file.exists(styles_css))
app_src    <- paste(readLines(app_R,      warn = FALSE), collapse = "\n")
styles_src <- paste(readLines(styles_css, warn = FALSE), collapse = "\n")

## 1. The user-visible "Server / remote use" header tag is gone.
stopifnot_msg(
  !grepl('tags\\$span\\("Server / remote use"\\)', app_src, fixed = FALSE),
  "tags$span(\"Server / remote use\") header was removed from app.R"
)

## 2. The collapsible card chrome (welcome-server class + the
##    chevron caret + the welcome-code <pre><code> block) are all
##    gone from app.R.
stopifnot_msg(
  !grepl("welcome-server", app_src, fixed = TRUE),
  "no .welcome-server class is referenced from app.R"
)
stopifnot_msg(
  !grepl("welcome-code", app_src, fixed = TRUE),
  "no .welcome-code class is referenced from app.R"
)

## 3. The remote-launch command is no longer *rendered* into the
##    Home tab as a <pre><code> block. The string itself can still
##    appear in source-level comments (the file-header docblock
##    keeps it as a hint for developers), so we look only for the
##    actual rendered HTML pattern, not for the substring.
stopifnot_msg(
  !grepl('<pre class="welcome-code"><code>', app_src, fixed = TRUE),
  "no <pre class=\"welcome-code\"><code> block is rendered into the Home tab"
)
stopifnot_msg(
  !grepl('PhenoMapR::run_app\\(host = "0\\.0\\.0\\.0", port = 3838\\)</code></pre>',
         app_src),
  "the run_app(host=...) <code></code></pre> tail was removed from app.R"
)

## 4. The companion CSS rules are gone from styles.css. We look
##    specifically for actual CSS selectors -- `.welcome-X`
##    followed (after optional whitespace) by a CSS rule-block
##    opener `{`, a child combinator `>`, sibling combinators
##    `+` / `~`, an attribute bracket `[`, a pseudo-class colon
##    `:`, or a comma starting another selector. This skips the
##    explanatory breadcrumb comment that mentions the names in
##    backticks.
selector_pattern <- function(class_token) {
  paste0("\\.", class_token, "[^\\w-]*?[{>+~\\[:,]")
}
stopifnot_msg(
  !grepl(selector_pattern("welcome-server"), styles_src, perl = TRUE),
  ".welcome-server CSS selectors were removed from styles.css"
)
stopifnot_msg(
  !grepl(selector_pattern("welcome-code"), styles_src, perl = TRUE),
  ".welcome-code CSS selectors were removed from styles.css"
)

## 5. The (lighter-touch) `.welcome-section-body p code` rule that
##    is shared with the rest of the welcome panel survived the
##    cleanup -- removing it would regress the inline `<code>`
##    formatting on the rest of the Home tab.
stopifnot_msg(
  grepl("\\.welcome-section-body p code", styles_src),
  ".welcome-section-body p code rule still exists for inline <code>"
)

cat("\nAll welcome-page cleanup smoke tests passed.\n")
