#' Package imports for base and stats functions
#'
#' @name PhenoMapR-imports
#' @importFrom stats na.omit sd median cor setNames
#' @importFrom graphics hist abline legend
#' @importFrom utils data head
#' @importFrom dplyr %>%
#' @keywords internal
NULL

# Tidy eval pronouns used in dplyr::filter(); avoid R CMD check false positive
if (getRversion() >= "2.15.1") utils::globalVariables(c(".data", ".env"))
