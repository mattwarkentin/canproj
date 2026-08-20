#' Extract Model Projections
#'
#' Extract projections from a canproj model object.
#'
#' @param object Model object.
#' @param ... Additional arguments passed to methods.
#'
#' @export
get_projections <- function(object, ...) {
  UseMethod("get_projections")
}
