#' rcea: package to perform cumulative effects assessments
#'
#' @docType package
#' @name rcea
#' 
#' @description 
#' Package to perform cumulative effects assessments.
#' 
#' @importFrom cli symbol
#' @importFrom glue glue
#' @importFrom grDevices dev.off png colorRampPalette
#' @importFrom graphics par box layout lines mtext
#' @importFrom graphics polygon text
#' @importFrom rlang sym
#' @importFrom stats setNames
# needed to get the namespace attached see https://tshafer.com/blog/2020/08/r-packages-s3-methods
#' @importFrom stars st_as_stars
#' @importFrom yaml yaml.load_file write_yaml read_yaml
NULL

globalVariables(c(
  "cumulative_effects", "cumulative_footprint", "direct", "effect",
  "effect_i", "effect_j", "effect_k", "effect.x", "effect.y", "i",
  "id_cell", "j", "k", "km2", "M", "pid_i", "pid_j", "pid_k", "pj",
  "pk", "pos", "pos_i", "pos_j", "pos_k", "presence", "psum", "Sensitivity",
  "sensitivity_1", "setNames", "sp", "Species", "type", "uid",
  "value", "vc_id", "weight", "x", "y"
))




# ------------------------------------------------------------------------------
# Gracieuset\u00e9 de Kevin Cazelles: https://github.com/KevCaz
# my simple(r) version of use template
use_template <- function(template, save_as = stdout(), pkg = "rcea", ...) {
  template <- readLines(
    fs::path_package(package = pkg, template)
  )
  # NB by default whisker forward the parent envi and I used this
  writeLines(whisker::whisker.render(template, ...), save_as)
}

#' Check if folder exists and create if not
#'
#' @param path path of folder to use as output, create if it does not already exist
#'
#' @export
chk_create <- function(path) {
  if (!file.exists(path)) dir.create(path, recursive = TRUE)
}