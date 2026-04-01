#' Package imports for base and stats functions
#'
#' @name PhenoMapR-imports
#' @importFrom stats na.omit sd median cor setNames
#' @importFrom utils data head
#' @importFrom dplyr %>%
#' @keywords internal
NULL

# Tidy eval / ggplot2 computed vars; avoid R CMD check false positives
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".data", ".env", "x"))
}
