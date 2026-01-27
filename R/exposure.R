#' Cumulative exposure assessment
#'
#' Assessment of the exposure (i.e. overlap) between valued components and environmental drivers in the context of cumulative effects assessments.
#'
#' @eval arguments(c("drivers","vc"))
#' @param exportAs string, "SpatRaster" (stack of layers).
#' @param align alignment policy, one of "error", "reproject", "template".
#' @param template optional SpatRaster to align to when `align = "template"`.
#' @param cores number of threads for terra::app (default uses terraOptions).
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
#' expo <- exposure(drivers, vc, exportAs = "SpatRaster")
#' expo
#' names(expo)
#'
#' @export
exposure <- function(drivers,
                     vc,
                     exportAs = c("SpatRaster"),
                     align = "error",
                     template = NULL,
                     cores = NULL) {

  exportAs <- match.arg(exportAs)
  out <- align_pair(drivers, vc, align = align, template = template)
  drivers <- out$drivers
  vc <- out$vc

  dr_names <- names(drivers)
  vc_names <- names(vc)

  # Compute overlap for each driver with all vc
  exp_list <- lapply(seq_along(dr_names), function(i) {
    r <- vc * drivers[[i]]
    names(r) <- vc_names
    r
  })
  names(exp_list) <- dr_names

  # Stack all exposures with layer names vc_driver
  layers <- mapply(function(r, dr) {
    names(r) <- paste0(names(r), "_", dr)
    r
  }, exp_list, names(exp_list), SIMPLIFY = FALSE)
  exp_r <- terra::rast(layers)
  # ensure names are carried through
  names(exp_r) <- unlist(lapply(layers, names), use.names = FALSE)
  exp_r
}

#' Run exposure from a cube
#'
#' Convenience wrapper around `exposure()` that uses stacks stored in
#' an `rcea_cube`.
#'
#' @param cube rcea_cube with `stack$drivers` and `stack$vc`.
#' @param ... Additional arguments passed to `exposure()`.
#'
#' @examples
#' layers_path <- system.file("extdata/catalog/layers.csv", package = "rcea")
#' groups_path <- system.file("extdata/catalog/groups.yaml", package = "rcea")
#' catalog <- load_catalog(layers_path, groups_path)
#' cube <- make_cube(catalog)
#' cube <- stack_layers(cube)
#' expo <- exposure_cube(cube)
#' expo
#'
#' @export
exposure_cube <- function(cube, ...) {
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

  exposure(drivers, vc, ...)
}

# Internal alignment helper
align_pair <- function(drivers, vc, align = "error", template = NULL) {
  same_crs <- terra::compareGeom(drivers, vc, stopOnError = FALSE, crs = TRUE, ext = FALSE, rowcol = FALSE)
  same_res <- all.equal(terra::res(drivers), terra::res(vc))
  same_ext <- terra::ext(drivers) == terra::ext(vc)
  if (isTRUE(same_crs) && isTRUE(same_res) && isTRUE(same_ext)) {
    return(list(drivers = drivers, vc = vc))
  }

  if (align == "error") {
    stop("Drivers and vc differ in CRS/resolution/extent; set align = 'reproject' or 'template'.", call. = FALSE)
  }

  if (align == "template") {
    if (is.null(template)) stop("Template required when align = 'template'.", call. = FALSE)
    tmpl <- template
  } else {
    tmpl <- drivers[[1]]
  }

  drivers <- terra::project(drivers, tmpl)
  vc <- terra::project(vc, tmpl)
  list(drivers = drivers, vc = vc)
}
