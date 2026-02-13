#' Extract layers from a stack or matrix
#'
#' @param dat SpatRaster, matrix (with `template` attribute), or data.frame.
#' @param layer_ids Optional exact layer names to keep.
#' @param drivers Optional driver names to keep (matches suffix after `_`).
#' @param vcs Optional valued component names to keep (matches prefix before `_`).
#' @return Same type as input, with selected layers.
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
#' # Keep only the shipping driver layers
#' layers_extract(ce, drivers = "shipping")
#'
#' @export
layers_extract <- function(dat, layer_ids = NULL, drivers = NULL, vcs = NULL) {
  res <- extract_result(dat)
  dat <- res$data
  layer_map <- res$map

  if (is.matrix(dat)) {
    template <- attr(dat, "template")
    nms <- colnames(dat)
    keep <- rep(TRUE, length(nms))
    if (!is.null(layer_ids)) keep <- keep & nms %in% layer_ids
    if (!is.null(layer_map)) {
      idx <- rep(TRUE, nrow(layer_map))
      if (!is.null(vcs)) idx <- idx & layer_map$vc %in% vcs
      if (!is.null(drivers)) idx <- idx & layer_map$driver %in% drivers
      keep <- keep & nms %in% layer_map$layer[idx]
    } else {
      if (!is.null(vcs)) keep <- keep & sub("_.*", "", nms) %in% vcs
      if (!is.null(drivers)) keep <- keep & sub(".*_", "", nms) %in% drivers
    }
    dat <- dat[, keep, drop = FALSE]
    dat <- with_template(dat, template)
    attr(dat, "layer_map") <- if (is.null(layer_map)) NULL else layer_map[layer_map$layer %in% colnames(dat), , drop = FALSE]
    return(dat)
  }

  if (inherits(dat, "SpatRaster")) {
    nms <- names(dat)
    keep <- rep(TRUE, length(nms))
    if (!is.null(layer_ids)) keep <- keep & nms %in% layer_ids
    if (!is.null(layer_map)) {
      idx <- rep(TRUE, nrow(layer_map))
      if (!is.null(vcs)) idx <- idx & layer_map$vc %in% vcs
      if (!is.null(drivers)) idx <- idx & layer_map$driver %in% drivers
      keep <- keep & nms %in% layer_map$layer[idx]
    } else {
      if (!is.null(vcs)) keep <- keep & sub("_.*", "", nms) %in% vcs
      if (!is.null(drivers)) keep <- keep & sub(".*_", "", nms) %in% drivers
    }
    out <- dat[[which(keep)]]
    attr(out, "layer_map") <- if (is.null(layer_map)) NULL else layer_map[layer_map$layer %in% names(out), , drop = FALSE]
    return(out)
  }

  # data.frame fallback (legacy)
  if (is.data.frame(dat)) {
    if (!is.null(drivers) && "drivers" %in% colnames(dat)) {
      dat <- dat[dat$drivers %in% drivers, , drop = FALSE]
    }
    if (!is.null(layer_ids)) {
      dat <- dplyr::select(dat, dplyr::any_of(c("x", "y", "drivers")), dplyr::all_of(layer_ids))
    }
    if (!is.null(vcs)) {
      dat <- dplyr::select(dat, dplyr::any_of(c("x", "y", "drivers")), dplyr::all_of(vcs))
    }
    return(dat)
  }

  dat
}

#' Aggregate layers in a stack or matrix
#'
#' @param dat SpatRaster, matrix (with `template` attribute), or data.frame.
#' @param layer_ids Optional exact layer names to keep before aggregating.
#' @param drivers Optional driver names to keep (matches suffix after `_`).
#' @param vcs Optional valued component names to keep (matches prefix before `_`).
#' @param by Aggregation mode: "both" (single layer), "drivers" (aggregate over drivers per VC),
#'   "vcs" (aggregate over VCs per driver), or "none" (no aggregation).
#' @param group_by Optional metadata columns in `layer_map` used to group layers before aggregation.
#'   When set, this takes precedence over `by`.
#' @param fun Aggregation function: "sum", "mean", "median", "min", "max", "sd".
#' @param exportAs string, "SpatRaster" or "matrix".
#' @return Aggregated object (SpatRaster by default).
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
#' # Aggregate all selected layers to a single layer
#' layers_aggregate(ce, by = "both")
#' # Aggregate over drivers (one layer per VC)
#' layers_aggregate(ce, by = "drivers")
#' # Use a different aggregation function
#' layers_aggregate(ce, by = "drivers", fun = "mean")
#'
#' @export
layers_aggregate <- function(dat,
                             layer_ids = NULL,
                             drivers = NULL,
                             vcs = NULL,
                             by = c("both", "drivers", "vcs", "none"),
                             group_by = NULL,
                             fun = c("sum", "mean", "median", "min", "max", "sd"),
                             exportAs = c("SpatRaster", "matrix")) {
  by <- match.arg(by)
  fun <- match.arg(fun)
  exportAs <- match.arg(exportAs)

  res <- extract_result(dat)
  dat <- res$data
  layer_map <- res$map

  dat <- layers_extract(dat, layer_ids = layer_ids, drivers = drivers, vcs = vcs)
  if (is_result(dat)) {
    layer_map <- dat$layer_map
    dat <- dat$data
  } else if (is.null(layer_map)) {
    layer_map <- attr(dat, "layer_map")
  }

  agg_fun <- function(x) {
    switch(fun,
      "sum" = rowSums(x, na.rm = TRUE),
      "mean" = rowMeans(x, na.rm = TRUE),
      "median" = apply(x, 1, stats::median, na.rm = TRUE),
      "min" = apply(x, 1, min, na.rm = TRUE),
      "max" = apply(x, 1, max, na.rm = TRUE),
      "sd" = apply(x, 1, stats::sd, na.rm = TRUE)
    )
  }

  resolve_grouping <- function(layer_names) {
    if (is.null(group_by)) return(NULL)
    if (is.null(layer_map)) {
      stop("`group_by` requires layer metadata (`layer_map`) but none was found.", call. = FALSE)
    }

    map <- layer_map[match(layer_names, layer_map$layer), , drop = FALSE]
    if (anyNA(map$layer)) {
      stop("`group_by` could not align selected layers with `layer_map`.", call. = FALSE)
    }

    missing_group_cols <- setdiff(group_by, names(map))
    if (length(missing_group_cols)) {
      stop("Unknown `group_by` columns: ", paste(missing_group_cols, collapse = ", "), call. = FALSE)
    }

    grp_df <- map[, group_by, drop = FALSE]
    labels <- if (length(group_by) == 1) {
      as.character(grp_df[[1]])
    } else {
      apply(grp_df, 1, function(x) {
        x_chr <- as.character(x)
        x_chr[is.na(x_chr) | !nzchar(x_chr)] <- "NA"
        paste(paste0(group_by, "=", x_chr), collapse = "|")
      })
    }
    labels[is.na(labels) | !nzchar(labels)] <- "NA"

    first <- !duplicated(labels)
    grp_labels <- labels[first]
    grp_index <- lapply(grp_labels, function(lbl) which(labels == lbl))
    names(grp_index) <- grp_labels

    grp_map <- cbind(
      data.frame(layer = grp_labels, vc = NA_character_, driver = NA_character_, stringsAsFactors = FALSE),
      grp_df[first, , drop = FALSE]
    )
    rownames(grp_map) <- NULL

    list(index = grp_index, map = grp_map)
  }

  if (is.matrix(dat)) {
    template <- attr(dat, "template")
    nms <- colnames(dat)
    grp <- resolve_grouping(nms)
    if (!is.null(layer_map)) {
      vc_ids <- layer_map$vc[match(nms, layer_map$layer)]
      dr_ids <- layer_map$driver[match(nms, layer_map$layer)]
    } else {
      vc_ids <- sub("_.*", "", nms)
      dr_ids <- sub(".*_", "", nms)
    }

    out <- if (!is.null(grp)) {
      tmp <- lapply(grp$index, function(idx) agg_fun(dat[, idx, drop = FALSE]))
      out <- do.call(cbind, tmp)
      colnames(out) <- names(grp$index)
      out
    } else {
      switch(by,
        "none" = dat,
        "drivers" = {
          vc_levels <- unique(vc_ids)
          tmp <- lapply(vc_levels, function(v) agg_fun(dat[, vc_ids == v, drop = FALSE]))
          out <- do.call(cbind, tmp)
          colnames(out) <- vc_levels
          out
        },
        "vcs" = {
          dr_levels <- unique(dr_ids)
          tmp <- lapply(dr_levels, function(d) agg_fun(dat[, dr_ids == d, drop = FALSE]))
          out <- do.call(cbind, tmp)
          colnames(out) <- dr_levels
          out
        },
        "both" = {
          out <- matrix(agg_fun(dat), ncol = 1)
          colnames(out) <- "cumulative"
          out
        }
      )
    }

    out <- with_template(out, template)
    if (exportAs == "SpatRaster") {
      if (is.null(template)) stop("Matrix input lacks template attribute; cannot reconstruct raster.", call. = FALSE)
      out_r <- matrix_to_raster(out, template, colnames(out))
      if (!is.null(grp)) {
        attr(out_r, "layer_map") <- grp$map
      } else if (!is.null(layer_map)) {
        if (by == "drivers") {
          attr(out_r, "layer_map") <- data.frame(layer = names(out_r), vc = names(out_r), driver = NA, stringsAsFactors = FALSE)
        } else if (by == "vcs") {
          attr(out_r, "layer_map") <- data.frame(layer = names(out_r), vc = NA, driver = names(out_r), stringsAsFactors = FALSE)
        } else if (by == "both") {
          attr(out_r, "layer_map") <- data.frame(layer = "cumulative", vc = NA, driver = NA, stringsAsFactors = FALSE)
        } else {
          attr(out_r, "layer_map") <- layer_map[layer_map$layer %in% names(out_r), , drop = FALSE]
        }
      }
      return(out_r)
    }
    if (!is.null(grp)) {
      attr(out, "layer_map") <- grp$map
    } else if (!is.null(layer_map)) {
      if (by == "drivers") {
        attr(out, "layer_map") <- data.frame(layer = colnames(out), vc = colnames(out), driver = NA, stringsAsFactors = FALSE)
      } else if (by == "vcs") {
        attr(out, "layer_map") <- data.frame(layer = colnames(out), vc = NA, driver = colnames(out), stringsAsFactors = FALSE)
      } else if (by == "both") {
        attr(out, "layer_map") <- data.frame(layer = "cumulative", vc = NA, driver = NA, stringsAsFactors = FALSE)
      } else {
        attr(out, "layer_map") <- layer_map[layer_map$layer %in% colnames(out), , drop = FALSE]
      }
    }
    return(out)
  }

  if (inherits(dat, "SpatRaster")) {
    nms <- names(dat)
    grp <- resolve_grouping(nms)
    if (!is.null(layer_map)) {
      vc_ids <- layer_map$vc[match(nms, layer_map$layer)]
      dr_ids <- layer_map$driver[match(nms, layer_map$layer)]
    } else {
      vc_ids <- sub("_.*", "", nms)
      dr_ids <- sub(".*_", "", nms)
    }

    out <- if (!is.null(grp)) {
      tmp <- lapply(grp$index, function(idx) {
        terra::app(dat[[idx]], get(fun, mode = "function"), na.rm = TRUE)
      })
      names(tmp) <- names(grp$index)
      terra::rast(tmp)
    } else {
      switch(by,
        "none" = dat,
        "drivers" = {
          vc_levels <- unique(vc_ids)
          tmp <- lapply(vc_levels, function(v) {
            terra::app(dat[[vc_ids == v]], get(fun, mode = "function"), na.rm = TRUE)
          })
          names(tmp) <- vc_levels
          terra::rast(tmp)
        },
        "vcs" = {
          dr_levels <- unique(dr_ids)
          tmp <- lapply(dr_levels, function(d) {
            terra::app(dat[[dr_ids == d]], get(fun, mode = "function"), na.rm = TRUE)
          })
          names(tmp) <- dr_levels
          terra::rast(tmp)
        },
        "both" = {
          terra::app(dat, get(fun, mode = "function"), na.rm = TRUE)
        }
      )
    }

    if (exportAs == "matrix") {
      out_mat <- raster_to_matrix(out)
      if (!is.null(grp)) {
        attr(out_mat, "layer_map") <- grp$map
      } else if (!is.null(layer_map)) {
        if (by == "drivers") {
          attr(out_mat, "layer_map") <- data.frame(layer = colnames(out_mat), vc = colnames(out_mat), driver = NA, stringsAsFactors = FALSE)
        } else if (by == "vcs") {
          attr(out_mat, "layer_map") <- data.frame(layer = colnames(out_mat), vc = NA, driver = colnames(out_mat), stringsAsFactors = FALSE)
        } else if (by == "both") {
          attr(out_mat, "layer_map") <- data.frame(layer = "cumulative", vc = NA, driver = NA, stringsAsFactors = FALSE)
        } else {
          attr(out_mat, "layer_map") <- layer_map[layer_map$layer %in% colnames(out_mat), , drop = FALSE]
        }
      }
      return(out_mat)
    }
    if (!is.null(grp)) {
      attr(out, "layer_map") <- grp$map
    } else if (!is.null(layer_map)) {
      if (by == "drivers") {
        attr(out, "layer_map") <- data.frame(layer = names(out), vc = names(out), driver = NA, stringsAsFactors = FALSE)
      } else if (by == "vcs") {
        attr(out, "layer_map") <- data.frame(layer = names(out), vc = NA, driver = names(out), stringsAsFactors = FALSE)
      } else if (by == "both") {
        attr(out, "layer_map") <- data.frame(layer = "cumulative", vc = NA, driver = NA, stringsAsFactors = FALSE)
      } else {
        attr(out, "layer_map") <- layer_map[layer_map$layer %in% names(out), , drop = FALSE]
      }
    }
    return(out)
  }

  # data.frame fallback (legacy)
  if (is.data.frame(dat)) {
    if (!is.null(group_by)) {
      stop("`group_by` is not supported for legacy data.frame inputs.", call. = FALSE)
    }
    if (by == "none") {
      return(dat)
    }
    if (by == "drivers") {
      return(dat |>
        dplyr::select(-drivers) |>
        dplyr::group_by(x, y) |>
        dplyr::summarise(dplyr::across(dplyr::everything(), \(x) sum(x, na.rm = TRUE))) |>
        dplyr::ungroup())
    }
    if (by == "vcs") {
      return(dat |>
        dplyr::mutate(value = rowSums(dplyr::pick(-x, -y, -drivers), na.rm = TRUE)) |>
        dplyr::select(x, y, drivers, value) |>
        tidyr::pivot_wider(names_from = "drivers", values_from = "value"))
    }
    if (by == "both") {
      return(dat |>
        dplyr::mutate(value = rowSums(dplyr::pick(-x, -y, -drivers), na.rm = TRUE)) |>
        dplyr::select(x, y, drivers, value) |>
        tidyr::pivot_wider(names_from = "drivers", values_from = "value") |>
        dplyr::mutate(cumulative = rowSums(dplyr::pick(-x, -y), na.rm = TRUE)) |>
        dplyr::select(x, y, cumulative))
    }
  }

  dat
}

#' Compute per-area effects for layered outputs
#'
#' Calculates per-unit-area values for each driver/VC layer using an area mask
#' derived from `vc`.
#'
#' @param dat SpatRaster or matrix (with `template` attribute) containing layers
#'   named as `vc_driver`.
#' @param vc SpatRaster or matrix (with `template` attribute) of valued components.
#' @param drivers Optional driver names to keep (matches suffix after `_`).
#' @param vcs Optional valued component names to keep (matches prefix before `_`).
#' @param layer_ids Optional exact layer names to keep before aggregation.
#' @param area_unit Unit for area calculation passed to `terra::cellSize()` (e.g., "km").
#' @return data.frame with one row per VC and one column per driver.
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
#' layers_per_area(ce, vc)
#'
#' @export
layers_per_area <- function(dat,
                            vc,
                            drivers = NULL,
                            vcs = NULL,
                            layer_ids = NULL,
                            area_unit = "km") {
  res <- extract_result(dat)
  dat <- res$data
  layer_map <- res$map

  dat <- layers_extract(dat, layer_ids = layer_ids, drivers = drivers, vcs = vcs)
  if (is_result(dat)) {
    layer_map <- dat$layer_map
    dat <- dat$data
  } else if (is.null(layer_map)) {
    layer_map <- attr(dat, "layer_map")
  }

  # Subset vc if requested
  if (!is.null(vcs)) {
    if (inherits(vc, "SpatRaster")) {
      vc <- vc[[names(vc) %in% vcs]]
    } else if (is.matrix(vc)) {
      template <- attr(vc, "template")
      vc <- vc[, colnames(vc) %in% vcs, drop = FALSE]
      vc <- with_template(vc, template)
    }
  }

  # Compute area per VC
  if (inherits(vc, "SpatRaster")) {
    area_r <- terra::cellSize(vc[[1]], unit = area_unit)
    vc_area <- terra::global(vc * area_r, "sum", na.rm = TRUE)[["sum"]]
    names(vc_area) <- names(vc)
  } else if (is.matrix(vc)) {
    template <- attr(vc, "template")
    if (is.null(template)) stop("Matrix vc lacks template attribute; cannot compute area.", call. = FALSE)
    area_r <- terra::cellSize(template, unit = area_unit)
    area_vec <- terra::values(area_r, mat = FALSE)
    vc_area <- colSums(vc * area_vec, na.rm = TRUE)
    names(vc_area) <- colnames(vc)
  } else {
    stop("vc must be a SpatRaster or matrix with template attribute.", call. = FALSE)
  }

  # Sum effects per layer
  if (inherits(dat, "SpatRaster")) {
    layer_names <- names(dat)
    layer_sums <- terra::global(dat, "sum", na.rm = TRUE)[["sum"]]
  } else if (is.matrix(dat)) {
    layer_names <- colnames(dat)
    layer_sums <- colSums(dat, na.rm = TRUE)
  } else {
    stop("dat must be a SpatRaster or matrix with template attribute.", call. = FALSE)
  }

  if (!is.null(layer_map)) {
    vc_ids <- layer_map$vc[match(layer_names, layer_map$layer)]
    dr_ids <- layer_map$driver[match(layer_names, layer_map$layer)]
  } else {
    vc_ids <- sub("_.*", "", layer_names)
    dr_ids <- sub(".*_", "", layer_names)
  }
  df <- data.frame(vc = vc_ids, drivers = dr_ids, total = layer_sums)

  df <- tidyr::pivot_wider(df, names_from = drivers, values_from = total)
  df <- df[match(names(vc_area), df$vc), , drop = FALSE]

  driver_cols <- setdiff(names(df), "vc")
  df[driver_cols] <- df[driver_cols] / vc_area[df$vc]
  unit_obj <- switch(area_unit,
    "km" = units::as_units(1, "km")^-2,
    "m" = units::as_units(1, "m")^-2,
    "ha" = units::as_units(1, "ha")^-1,
    units::as_units(1, area_unit)^-2
  )
  df[driver_cols] <- lapply(df[driver_cols], function(x) x * unit_obj)
  df
}
