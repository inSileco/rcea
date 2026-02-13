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
#' @param driver_meta Optional data.frame with metadata for driver layers (matched by `layer` or `layer_id`).
#' @param vc_meta Optional data.frame with metadata for valued component layers (matched by `layer` or `layer_id`).
#' @param pair_by Optional vector of metadata keys that must match between VC and drivers (e.g. "month", "year", "time").
#'   When set, only matched VC-driver pairs are computed.
#' @param pair_missing Pairing behavior when one side is missing a `pair_by` key:
#'   `"broadcast"` (default) allows the pair; `"strict"` requires both sides to be present and equal.
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
                engine = "matrix",
                driver_meta = NULL,
                vc_meta = NULL,
                pair_by = NULL,
                pair_missing = c("broadcast", "strict")) {
  engine <- match.arg(engine, c("matrix", "terra"))
  exportAs <- match.arg(exportAs, c("SpatRaster", "matrix"))
  pair_missing <- match.arg(pair_missing)

  # Align rasters if needed
  out_align <- align_pair(drivers, vc, align = align, template = template)
  drivers <- out_align$drivers
  vc <- out_align$vc

  if (is.null(driver_meta)) driver_meta <- attr(drivers, "layer_meta")
  if (is.null(vc_meta)) vc_meta <- attr(vc, "layer_meta")

  nmDr <- names(drivers)
  nmVC <- names(vc)
  nvc <- length(nmVC)
  ndr <- length(nmDr)

  layer_info <- build_layer_map(drivers, vc, driver_meta = driver_meta, vc_meta = vc_meta)
  full_layer_names <- layer_info$layer_names
  full_layer_map <- layer_info$layer_map

  vc_keys_primary <- if ("vc_vc_id" %in% names(full_layer_map)) {
    as.character(full_layer_map$vc_vc_id[seq_len(nvc)])
  } else {
    as.character(full_layer_map$vc[seq_len(nvc)])
  }
  dr_keys_primary <- if ("driver_driver_id" %in% names(full_layer_map)) {
    as.character(full_layer_map$driver_driver_id[seq.int(1, nrow(full_layer_map), by = nvc)])
  } else {
    as.character(full_layer_map$driver[seq.int(1, nrow(full_layer_map), by = nvc)])
  }
  vc_keys_fallback <- as.character(full_layer_map$vc[seq_len(nvc)])
  dr_keys_fallback <- as.character(full_layer_map$driver[seq.int(1, nrow(full_layer_map), by = nvc)])

  if (is.null(rownames(sensitivity)) || is.null(colnames(sensitivity))) {
    stop("sensitivity must have row and column names for VC and driver IDs.", call. = FALSE)
  }

  vc_keys <- if (all(unique(vc_keys_primary) %in% rownames(sensitivity))) {
    vc_keys_primary
  } else {
    vc_keys_fallback
  }
  dr_keys <- if (all(unique(dr_keys_primary) %in% colnames(sensitivity))) {
    dr_keys_primary
  } else {
    dr_keys_fallback
  }

  miss_vc <- setdiff(unique(vc_keys), rownames(sensitivity))
  miss_dr <- setdiff(unique(dr_keys), colnames(sensitivity))
  if (length(miss_vc)) {
    stop("Missing VC IDs in sensitivity rownames: ", paste(miss_vc, collapse = ", "), call. = FALSE)
  }
  if (length(miss_dr)) {
    stop("Missing driver IDs in sensitivity colnames: ", paste(miss_dr, collapse = ", "), call. = FALSE)
  }
  sensitivity <- sensitivity[vc_keys, dr_keys, drop = FALSE]

  pair_mask <- matrix(TRUE, nrow = nvc, ncol = ndr)
  if (!is.null(pair_by)) {
    pair_keys <- unique(sub("^(driver_|vc_)", "", as.character(pair_by)))
    for (k in pair_keys) {
      vc_col <- paste0("vc_", k)
      dr_col <- paste0("driver_", k)
      if (!vc_col %in% names(full_layer_map) || !dr_col %in% names(full_layer_map)) {
        stop("`pair_by` key '", k, "' requires both `", vc_col, "` and `", dr_col, "` in layer metadata.", call. = FALSE)
      }
      vc_vals <- as.character(full_layer_map[[vc_col]][seq_len(nvc)])
      dr_vals <- as.character(full_layer_map[[dr_col]][seq.int(1, nrow(full_layer_map), by = nvc)])
      pair_mask <- pair_mask & outer(vc_vals, dr_vals, function(v, d) {
        v <- trimws(v)
        d <- trimws(d)
        v_missing <- is.na(v) | !nzchar(v)
        d_missing <- is.na(d) | !nzchar(d)
        if (pair_missing == "strict") {
          !v_missing & !d_missing & (v == d)
        } else {
          v_missing | d_missing | (v == d)
        }
      })
    }
    if (!any(pair_mask)) {
      stop("No valid VC-driver pairs matched `pair_by` keys.", call. = FALSE)
    }
  }

  keep_idx <- which(as.vector(pair_mask))
  layer_names <- full_layer_names[keep_idx]
  layer_map <- full_layer_map[keep_idx, , drop = FALSE]

  if (engine == "matrix") {
    # In-memory matrix path: build vc_driver layers in allowed pair mask
    dr_mat <- terra::values(drivers, mat = TRUE)
    vc_mat <- terra::values(vc, mat = TRUE)
    out_list <- list()
    out_names <- character(0)
    for (j in seq_along(nmDr)) {
      idx_v <- which(pair_mask[, j])
      if (!length(idx_v)) next
      exp_j <- vc_mat[, idx_v, drop = FALSE] * dr_mat[, j]
      eff_j <- sweep(exp_j, MARGIN = 2, sensitivity[idx_v, j], `*`)
      out_list[[length(out_list) + 1]] <- eff_j
      out_names <- c(out_names, paste(nmVC[idx_v], nmDr[j], sep = "_"))
    }
    out_mat <- do.call(cbind, out_list)
    colnames(out_mat) <- out_names

    if (exportAs == "matrix") {
      out_mat <- with_template(out_mat, drivers[[1]])
      attr(out_mat, "layer_map") <- layer_map
      make_result(out_mat, layer_map, meta = list(engine = engine, exportAs = exportAs))
    } else {
      out_r <- matrix_to_raster(out_mat, drivers[[1]], out_names)
      attr(out_r, "layer_map") <- layer_map
      make_result(out_r, layer_map, meta = list(engine = engine, exportAs = exportAs))
    }
  } else {
    # Combined stack to avoid building exposure separately
    stk <- c(vc, drivers)

    fun_ce <- function(x) {
      v <- x[seq_len(nvc)]
      d <- x[(nvc + 1):(nvc + ndr)]
      eff <- (v %o% d) * sensitivity
      as.vector(eff)[keep_idx]
    }

    args <- list(x = stk, fun = fun_ce)
    if (!is.null(cores)) args$cores <- cores
    if (!is.null(filename)) args$filename <- filename
    res <- do.call(terra::app, args)
    names(res) <- layer_names
    if (exportAs == "matrix") {
      out_mat <- raster_to_matrix(res)
      attr(out_mat, "layer_map") <- layer_map
      make_result(out_mat, layer_map, meta = list(engine = engine, exportAs = exportAs))
    } else {
      attr(res, "layer_map") <- layer_map
      make_result(res, layer_map, meta = list(engine = engine, exportAs = exportAs))
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
#'     c("cod", "salmon"),
#'     c("shipping", "climate")
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

  cea(
    drivers,
    vc,
    sensitivity,
    driver_meta = attr(drivers, "layer_meta"),
    vc_meta = attr(vc, "layer_meta"),
    ...
  )
}
