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
build_layer_map <- function(drivers, vc, driver_meta = NULL, vc_meta = NULL) {
  dr_names <- names(drivers)
  vc_names <- names(vc)
  layer_names <- as.vector(outer(vc_names, dr_names, paste, sep = "_"))

  normalize_meta <- function(meta, layer_names, id_col) {
    if (is.null(meta)) {
      out <- data.frame(layer = layer_names, stringsAsFactors = FALSE)
    } else {
      out <- as.data.frame(meta, stringsAsFactors = FALSE)
      if (!"layer" %in% names(out)) {
        if ("layer_id" %in% names(out)) {
          out$layer <- out$layer_id
        } else {
          stop("layer metadata must include `layer` or `layer_id` column.", call. = FALSE)
        }
      }
      out <- out[out$layer %in% layer_names, , drop = FALSE]
      out <- out[match(layer_names, out$layer), , drop = FALSE]
    }

    if (!id_col %in% names(out)) {
      out[[id_col]] <- out$layer
    }

    if ("time" %in% names(out)) {
      time_chr <- as.character(out$time)
      time_chr[!nzchar(time_chr)] <- NA_character_
      if (!"year" %in% names(out)) {
        out$year <- ifelse(grepl("^\\d{4}", time_chr), substr(time_chr, 1, 4), NA_character_)
      }
      if (!"month" %in% names(out)) {
        out$month <- ifelse(grepl("^\\d{4}-\\d{2}", time_chr), substr(time_chr, 1, 7), NA_character_)
      }
    }

    out
  }

  dr_meta <- normalize_meta(driver_meta, dr_names, "driver_id")
  vc_meta <- normalize_meta(vc_meta, vc_names, "vc_id")

  vc_idx <- rep(seq_along(vc_names), times = length(dr_names))
  dr_idx <- rep(seq_along(dr_names), each = length(vc_names))

  layer_map <- data.frame(
    layer = layer_names,
    vc = vc_names[vc_idx],
    driver = dr_names[dr_idx],
    stringsAsFactors = FALSE
  )

  vc_cols <- setdiff(names(vc_meta), "layer")
  for (nm in vc_cols) {
    layer_map[[paste0("vc_", nm)]] <- vc_meta[[nm]][vc_idx]
  }

  dr_cols <- setdiff(names(dr_meta), "layer")
  for (nm in dr_cols) {
    layer_map[[paste0("driver_", nm)]] <- dr_meta[[nm]][dr_idx]
  }

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
