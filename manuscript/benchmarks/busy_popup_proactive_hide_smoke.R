#!/usr/bin/env Rscript
## Smoke tests for the busy-popup visibility model.
##
## Architecture
## ============
## R-driven custom-message visibility, with client-side debounce:
##
##   * R calls phenomapr_busy_show(...) -> sends `phenomapr-busy-show`
##     custom message. JS handleShow() stores the text and starts a
##     SHOW_DELAY_MS timer. Popup renders only if no busy_hide arrives
##     before the timer fires.
##
##   * R calls phenomapr_busy_hide() -> sends `phenomapr-busy-hide`
##     custom message. JS handleHide() cancels any pending show timer
##     and hides the popup if visible.
##
##   * Critically, R-side observers must defer their heavy synchronous
##     work into later::later(fn, delay = 0) AFTER calling
##     phenomapr_busy_show() so the websocket queue can flush the show
##     message to the browser BEFORE libuv gets blocked. Without that
##     deferral the show + hide messages queue together server-side
##     and arrive back-to-back on the client (popup never appears),
##     because R is single-threaded and a long synchronous observer
##     blocks libuv's I/O loop. See helpers.R's busy-overlay
##     file-comment for the full rationale.
##
## Defence in depth: shiny:idle triggers a (delayed) advisory hide,
## shiny:disconnected forces hide, and a 3-min absolute watchdog
## force-hides any popup that survives all of the above.
##
## This smoke test asserts the structural shape of the JS so future
## regressions get caught before reaching users.

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

helpers_R <- system.file("shiny", "helpers.R", package = "PhenoMapR")
dev_dir <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  helpers_R <- file.path(dev_dir, "helpers.R")
}
stopifnot(file.exists(helpers_R))
helpers_src <- paste(readLines(helpers_R, warn = FALSE), collapse = "\n")

## ---- Tunable timing constants ---------------------------------------------
stopifnot_msg(
  grepl("var SHOW_DELAY_MS\\s*=\\s*[0-9]+;", helpers_src),
  "helpers.R declares the SHOW_DELAY_MS debounce constant"
)
stopifnot_msg(
  grepl("var IDLE_SETTLE_MS\\s*=\\s*[0-9]+;", helpers_src),
  "helpers.R declares the IDLE_SETTLE_MS settle constant"
)
stopifnot_msg(
  grepl("var ABSOLUTE_TIMEOUT\\s*=\\s*[0-9]+;", helpers_src),
  "helpers.R declares the ABSOLUTE_TIMEOUT hard-cap constant"
)

## ---- R-side message handlers DRIVE visibility ------------------------------
## handleShow() schedules the SHOW_DELAY_MS debounce timer; if no
## hide arrives before it fires, the popup is rendered. handleHide()
## cancels any pending timer and hides the popup if visible. This is
## the visibility gate; previous busy-class-poll model could never
## see "busy" during a synchronous observer body, so it could not
## drive the popup in the very scenario it was needed for.
hs_anchor <- regexpr("function handleShow\\(", helpers_src)
stopifnot_msg(hs_anchor > 0L, "located handleShow() body")
hs_window <- substr(helpers_src, hs_anchor,
                    min(nchar(helpers_src), hs_anchor + 1500L))
stopifnot_msg(
  grepl("SHOW_DELAY_MS", hs_window),
  "handleShow() applies the SHOW_DELAY_MS debounce"
)
stopifnot_msg(
  grepl("setTimeout", hs_window),
  "handleShow() uses setTimeout for the debounce timer"
)
stopifnot_msg(
  grepl("renderShow\\(\\)", hs_window),
  "handleShow() ultimately renders the popup once the debounce elapses"
)
stopifnot_msg(
  grepl("pendingShow", hs_window),
  "handleShow() tracks a pendingShow flag so handleHide() can cancel"
)

hh_anchor <- regexpr("function handleHide\\(", helpers_src)
stopifnot_msg(hh_anchor > 0L, "located handleHide() body")
hh_window <- substr(helpers_src, hh_anchor,
                    min(nchar(helpers_src), hh_anchor + 800L))
stopifnot_msg(
  grepl("clearPendingShow|clearTimeout", hh_window),
  "handleHide() cancels any pending show timer"
)
## handleHide() goes through forceHide(), the unified hide path that
## probes the DOM directly so a drifted JS `isVisible` flag cannot
## leave the popup stuck on screen. forceHide() is itself defined in
## helpers.R; renderHide() is kept around as a thin forwarder.
stopifnot_msg(
  grepl("forceHide\\(", hh_window) || grepl("renderHide\\(\\)", hh_window),
  "handleHide() hides the popup if it was visible (forceHide / renderHide)"
)
stopifnot_msg(
  grepl("resetMessage\\(\\)", hh_window),
  "handleHide() resets the stored message text"
)
stopifnot_msg(
  grepl("lastHideAt", hh_window),
  "handleHide() records lastHideAt so the watchdog can detect stale-hide"
)

## ---- shiny:idle is an ADVISORY safety hide --------------------------------
## If R signals busy_show and then forgets to call busy_hide, the
## popup still goes away when shiny becomes idle for IDLE_SETTLE_MS.
idle_anchor <- regexpr("function onShinyIdle\\(\\)", helpers_src)
stopifnot_msg(idle_anchor > 0L, "located onShinyIdle() body")
idle_window <- substr(helpers_src, idle_anchor,
                      min(nchar(helpers_src), idle_anchor + 2000L))
stopifnot_msg(
  grepl("IDLE_SETTLE_MS", idle_window),
  "onShinyIdle() schedules an advisory hide after IDLE_SETTLE_MS"
)
stopifnot_msg(
  grepl("forceHide\\(", idle_window) || grepl("renderHide\\(\\)", idle_window),
  "onShinyIdle() ultimately hides the popup if it stays idle"
)
## The advisory hide must probe the DOM in addition to the JS flag,
## because the JS `isVisible` flag can drift out of sync with the
## actual overlay element under hot reload / Shiny reconnect /
## partial-render paths. Without the DOM check the popup can get
## stuck on screen even though shiny is fully idle.
stopifnot_msg(
  grepl("overlayDomVisible\\(\\)", idle_window),
  "onShinyIdle() probes the DOM (overlayDomVisible) not just the isVisible flag"
)
stopifnot_msg(
  grepl("clearAllFileInputLoading\\(\\)", idle_window),
  "onShinyIdle() also clears the file-input loading bar"
)
stopifnot_msg(
  grepl("dismissedThisRun = false", idle_window),
  "onShinyIdle() clears dismissedThisRun so the next op can re-show"
)

## ---- shiny:disconnected force-hides ---------------------------------------
disc_anchor <- regexpr("function onShinyDisconnected\\(\\)", helpers_src)
stopifnot_msg(disc_anchor > 0L, "located onShinyDisconnected() body")
disc_window <- substr(helpers_src, disc_anchor,
                      min(nchar(helpers_src), disc_anchor + 600L))
stopifnot_msg(
  grepl("forceHide\\(", disc_window) || grepl("renderHide\\(\\)", disc_window),
  "onShinyDisconnected() force-hides the popup"
)

## ---- Watchdog: 3-minute absolute hard cap ---------------------------------
wd_anchor <- regexpr("function watchdogTick\\(\\)", helpers_src)
stopifnot_msg(wd_anchor > 0L, "located watchdogTick() body")
wd_window <- substr(helpers_src, wd_anchor,
                    min(nchar(helpers_src), wd_anchor + 1500L))
stopifnot_msg(
  grepl("ABSOLUTE_TIMEOUT", wd_window),
  "watchdogTick() enforces the ABSOLUTE_TIMEOUT hard cap"
)
stopifnot_msg(
  grepl("forceHide\\(", wd_window) || grepl("renderHide\\(\\)", wd_window),
  "watchdogTick() can force-hide a stuck popup"
)
## The watchdog must also catch the "stale-hide" case: R has signalled
## busy_hide more recently than busy_show but the popup is still up
## on screen (e.g. handleHide ran before renderShow could flip
## isVisible, or some race re-added the .is-visible class). Without
## this case the popup can stay stuck on screen indefinitely after
## the operation completes -- exactly the regression this guards.
stopifnot_msg(
  grepl("lastHideAt", wd_window) && grepl("lastShowAt", wd_window),
  "watchdogTick() detects stale-hide via lastHideAt / lastShowAt"
)
stopifnot_msg(
  grepl("requestAnimationFrame\\(watchdogTick\\)", helpers_src),
  "watchdogTick() runs as a requestAnimationFrame loop"
)

## ---- Custom-message handlers are registered -------------------------------
stopifnot_msg(
  grepl('Shiny\\.addCustomMessageHandler\\("phenomapr-busy-show",\\s*handleShow\\)',
        helpers_src),
  "phenomapr-busy-show is wired to handleShow()"
)
stopifnot_msg(
  grepl('Shiny\\.addCustomMessageHandler\\("phenomapr-busy-hide",\\s*handleHide\\)',
        helpers_src),
  "phenomapr-busy-hide is wired to handleHide()"
)

## CRITICAL: Shiny.addCustomMessageHandler enforces handler.length === 1
## via `if (handler.length !== 1) console.error(...)`. A zero-arg or
## multi-arg handler silently fails to register with an "Uncaught
## handler must be a function that takes one argument." error in
## the browser console. handleShow already takes its payload arg
## (`function handleShow(p)`); handleHide ignores its payload but
## MUST still declare exactly one parameter, hence the `_msg`
## placeholder. Removing it breaks the entire hide path AND pins
## dismissedThisRun to true after a manual dismiss (which then
## suppresses every subsequent show). Guard against that regression.
stopifnot_msg(
  grepl("function handleShow\\([A-Za-z_][A-Za-z0-9_]*\\)", helpers_src),
  "handleShow declares exactly one parameter (Shiny enforces handler.length === 1)"
)
stopifnot_msg(
  grepl("function handleHide\\([A-Za-z_][A-Za-z0-9_]*\\)", helpers_src),
  "handleHide declares exactly one parameter (Shiny enforces handler.length === 1)"
)

## ---- Event listeners must register shiny:idle / shiny:disconnected --------
stopifnot_msg(
  grepl('document\\.addEventListener\\("shiny:idle",\\s*onShinyIdle',
        helpers_src),
  "shiny:idle is wired via native addEventListener"
)
stopifnot_msg(
  grepl('document\\.addEventListener\\("shiny:disconnected",\\s*onShinyDisconnected',
        helpers_src),
  "shiny:disconnected is wired via native addEventListener"
)
stopifnot_msg(
  grepl('\\$j\\(document\\)\\.on\\("shiny:idle",\\s*onShinyIdle\\)',
        helpers_src),
  "shiny:idle is also wired via jQuery .on()"
)

## ---- R-side observers MUST defer heavy work via later::later --------------
## This is the OTHER half of the architecture. Without it, the popup
## still wouldn't appear during a slow synchronous observer because
## the show + hide custom messages would queue together inside the
## observer's flushReact and arrive at the client back-to-back. We
## audit the canonical heavy observers in app.R to confirm each
## phenomapr_busy_show is followed by a later::later block.
app_R <- system.file("shiny", "app.R", package = "PhenoMapR")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  app_R <- file.path(dev_dir, "app.R")
}
stopifnot(file.exists(app_R))
app_src <- paste(readLines(app_R, warn = FALSE), collapse = "\n")

heavy_show_calls <- c(
  "Reading expression file",
  "Loading cell metadata",
  "Cleaning matrix",
  "Computing PhenoMap scores",
  "Finding marker genes",
  "Building marker heatmap"
)
for (label in heavy_show_calls) {
  pat <- paste0(
    "phenomapr_busy_show\\([^)]*",
    gsub("\\.", "\\\\.", label, fixed = TRUE),
    "[\\s\\S]{0,1500}?later::later\\(function\\(\\)"
  )
  stopifnot_msg(
    grepl(pat, app_src, perl = TRUE),
    sprintf(
      "app.R: '%s' busy_show is followed by a later::later() deferral",
      label
    )
  )
}

cat("\nAll busy-popup smoke tests passed.\n")
