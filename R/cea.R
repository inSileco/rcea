#' Cumulative effects assessments
#'
#' Assessment of cumulative effects using the Halpern et al. 2008 method.
#'
#' @eval arguments(c("drivers", "vc", "sensitivity"))
#' @param exportAs string, "SpatRaster" or "matrix".
#' @param align alignment policy, one of "error", "reproject", "template".
#' @param template optional SpatRaster to align to when `align = "template"`.
#' @param cores number of threads for terra::app.
#' @param filename optional path to write the output raster stack (useful for large jobs).
#' @param engine calculation engine, one of "matrix" (fast, in-memory) or "terra" (streaming/chunked).
#'
#' @examples
#' drv_paths <- system.file(
#'   "extdata/rasters",
#'   c("pressure_shipping.tif", "pressure_climate.tif"),
#'   package = "rcea"
#' )
#' vc_paths <- system.file(
#'   "extdata/rasters",
#'   c("vc_cod.tif", "vc_salmon.tif"),
#'   package = "rcea"
#' )
#' drivers <- terra::rast(drv_paths)
#' names(drivers) <- c("shipping", "climate")
#' vc <- terra::rast(vc_paths)
#' names(vc) <- c("cod", "salmon")
#' sens <- matrix(
#'   c(
#'     0.8, 0.5,
#'     0.2, 0.7
#'   ),
#'   nrow = 2,
#'   dimnames = list(c("cod", "salmon"), c("shipping", "climate"))
#' )
#' ce <- cea(drivers, vc, sens, exportAs = "SpatRaster")
#' ce
#'
#' @export
cea <- function(drivers,
                vc,
                sensitivity,
                exportAs = "SpatRaster",
                align = "error",
                template = NULL,
                cores = NULL,
                filename = NULL,
                engine = "matrix") {
  engine <- match.arg(engine, c("matrix", "terra"))
  exportAs <- match.arg(exportAs, c("SpatRaster", "matrix"))

  # Align rasters if needed
  out_align <- align_pair(drivers, vc, align = align, template = template)
  drivers <- out_align$drivers
  vc <- out_align$vc

  nmDr <- names(drivers)
  nmVC <- names(vc)
  sensitivity <- sensitivity[nmVC, nmDr, drop = FALSE]

  if (engine == "matrix") {
    # In-memory matrix path: build all vc_driver layers at once
    dr_mat <- terra::values(drivers, mat = TRUE)
    vc_mat <- terra::values(vc, mat = TRUE)
    out_mat <- NULL
    for (j in seq_along(nmDr)) {
      exp_j <- vc_mat * dr_mat[, j]
      eff_j <- sweep(exp_j, MARGIN = 2, sensitivity[, j], `*`)
      out_mat <- cbind(out_mat, eff_j)
    }
    layer_names <- as.vector(outer(nmVC, nmDr, paste, sep = "_"))
    if (exportAs == "matrix") {
      colnames(out_mat) <- layer_names
      out_mat <- with_template(out_mat, drivers[[1]])
      out_mat
    } else {
      matrix_to_raster(out_mat, drivers[[1]], layer_names)
    }
  } else {
    # Combined stack to avoid building exposure separately
    stk <- c(vc, drivers)
    nvc <- terra::nlyr(vc)
    ndr <- terra::nlyr(drivers)
    layer_names <- as.vector(outer(nmVC, nmDr, paste, sep = "_"))

    fun_ce <- function(x) {
      v <- x[seq_len(nvc)]
      d <- x[(nvc + 1):(nvc + ndr)]
      eff <- (v %o% d) * sensitivity
      as.vector(eff)
    }

    args <- list(x = stk, fun = fun_ce)
    if (!is.null(cores)) args$cores <- cores
    if (!is.null(filename)) args$filename <- filename
    res <- do.call(terra::app, args)
    names(res) <- layer_names
    if (exportAs == "matrix") {
      raster_to_matrix(res)
    } else {
      res
    }
  }
}

#' Run CEA from a cube
#'
#' Convenience wrapper around `cea()` that uses stacks and sensitivity stored in
#' an `rcea_cube`.
#'
#' @param cube rcea_cube with `stack$drivers`, `stack$vc`, and optionally `sensitivity`.
#' @param sensitivity Optional vulnerability matrix; overrides `cube$sensitivity` when provided.
#' @param ... Additional arguments passed to `cea()`.
#'
#' @examples
#' layers_path <- system.file("extdata/catalog/layers.csv", package = "rcea")
#' groups_path <- system.file("extdata/catalog/groups.yaml", package = "rcea")
#' catalog <- load_catalog(layers_path, groups_path)
#' sens <- matrix(
#'   1,
#'   nrow = 2,
#'   ncol = 2,
#'   dimnames = list(
#'     c("vc_cod", "vc_salmon"),
#'     c("pressure_shipping", "pressure_climate")
#'   )
#' )
#' cube <- make_cube(catalog, sensitivity = sens)
#' cube <- stack_layers(cube)
#' ce <- cea_cube(cube, engine = "matrix")
#' ce
#'
#' @export
cea_cube <- function(cube, sensitivity = NULL, ...) {
  if (!inherits(cube, "rcea_cube")) stop("cube must be an rcea_cube.", call. = FALSE)
  if (is.null(cube$stack)) stop("cube$stack is empty; call stack_layers() first.", call. = FALSE)

  if (!is.null(cube$stack) && is.list(cube$stack)) {
    drivers <- cube$stack$drivers
    vc <- cube$stack$vc
  } else if (inherits(cube$stack, "SpatRaster")) {
    stop("cube$stack must include drivers and vc stacks; rebuild with stack_layers().", call. = FALSE)
  } else {
    stop("cube$stack must be a list with drivers and vc stacks.", call. = FALSE)
  }

  if (is.null(drivers) || is.null(vc)) {
    stop("cube$stack must include both drivers and vc stacks.", call. = FALSE)
  }

  if (is.null(sensitivity)) sensitivity <- cube$sensitivity
  if (is.null(sensitivity)) {
    stop("sensitivity not provided and cube$sensitivity is NULL.", call. = FALSE)
  }

  cea(drivers, vc, sensitivity, ...)
}
