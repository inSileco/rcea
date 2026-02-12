#' Construct an rcea result object
#'
#' @param data SpatRaster or matrix with template attribute.
#' @param layer_map data.frame with columns: layer, vc, driver.
#' @param meta list of metadata.
#' @return rcea_result object.
#'
#' @keywords internal
make_result <- function(data, layer_map, meta = list()) {
  structure(
    list(
      data = data,
      layer_map = layer_map,
      meta = meta
    ),
    class = "rcea_result"
  )
}

# internal helper to build layer names/map
#' @keywords internal
build_layer_map <- function(drivers, vc) {
  dr_names <- names(drivers)
  vc_names <- names(vc)
  layer_names <- as.vector(outer(vc_names, dr_names, paste, sep = "_"))
  layer_map <- data.frame(
    layer = layer_names,
    vc = rep(vc_names, times = length(dr_names)),
    driver = rep(dr_names, each = length(vc_names)),
    stringsAsFactors = FALSE
  )
  list(layer_names = layer_names, layer_map = layer_map)
}

#' @export
print.rcea_result <- function(x, ...) {
  cat("<rcea_result>\n")
  cat("  data:", class(x$data)[1], "\n")
  if (inherits(x$data, "SpatRaster")) {
    cat("  layers:", terra::nlyr(x$data), "\n")
  } else if (is.matrix(x$data)) {
    cat("  layers:", ncol(x$data), "\n")
  }
  if (!is.null(x$layer_map)) {
    cat("  layer_map:", nrow(x$layer_map), "rows\n")
  }
  invisible(x)
}

#' Export rcea_result to raster files and metadata
#'
#' @param result rcea_result
#' @param dir output directory
#' @param prefix filename prefix
#' @param overwrite logical, overwrite existing files
#' @param metadata name of metadata file (csv)
#' @return invisible list with paths
#'
#' @export
export_result <- function(result, dir, prefix = "result", overwrite = FALSE, metadata = "metadata.csv") {
  if (!inherits(result, "rcea_result")) stop("result must be an rcea_result.", call. = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  data <- result$data
  if (is.matrix(data)) {
    template <- attr(data, "template")
    if (is.null(template)) stop("Matrix result lacks template attribute.", call. = FALSE)
    data <- matrix_to_raster(data, template, colnames(data))
  }

  tif_path <- file.path(dir, paste0(prefix, ".tif"))
  terra::writeRaster(data, tif_path, overwrite = overwrite)

  meta_path <- file.path(dir, metadata)
  if (!is.null(result$layer_map)) {
    utils::write.csv(result$layer_map, meta_path, row.names = FALSE)
  }

  invisible(list(raster = tif_path, metadata = meta_path))
}

# internal helpers
is_result <- function(x) inherits(x, "rcea_result")

extract_result <- function(x) {
  if (is_result(x)) {
    return(list(data = x$data, map = x$layer_map))
  }
  if (inherits(x, "SpatRaster") || is.matrix(x)) {
    return(list(data = x, map = attr(x, "layer_map")))
  }
  list(data = x, map = NULL)
}
