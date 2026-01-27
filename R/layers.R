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
  if (is.matrix(dat)) {
    template <- attr(dat, "template")
    nms <- colnames(dat)
    keep <- rep(TRUE, length(nms))
    if (!is.null(layer_ids)) keep <- keep & nms %in% layer_ids
    if (!is.null(vcs)) keep <- keep & sub("_.*", "", nms) %in% vcs
    if (!is.null(drivers)) keep <- keep & sub(".*_", "", nms) %in% drivers
    dat <- dat[, keep, drop = FALSE]
    return(with_template(dat, template))
  }

  if (inherits(dat, "SpatRaster")) {
    nms <- names(dat)
    keep <- rep(TRUE, length(nms))
    if (!is.null(layer_ids)) keep <- keep & nms %in% layer_ids
    if (!is.null(vcs)) keep <- keep & sub("_.*", "", nms) %in% vcs
    if (!is.null(drivers)) keep <- keep & sub(".*_", "", nms) %in% drivers
    return(dat[[which(keep)]])
  }

  # data.frame fallback (legacy)
  if (is.data.frame(dat)) {
    if (!is.null(drivers) && "drivers" %in% colnames(dat)) {
      dat <- dplyr::filter(dat, drivers %in% .env$drivers)
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
                             fun = c("sum", "mean", "median", "min", "max", "sd"),
                             exportAs = c("SpatRaster", "matrix")) {
  by <- match.arg(by)
  fun <- match.arg(fun)
  exportAs <- match.arg(exportAs)

  dat <- layers_extract(dat, layer_ids = layer_ids, drivers = drivers, vcs = vcs)

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

  if (is.matrix(dat)) {
    template <- attr(dat, "template")
    nms <- colnames(dat)
    vc_ids <- sub("_.*", "", nms)
    dr_ids <- sub(".*_", "", nms)

    out <- switch(by,
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

    out <- with_template(out, template)
    if (exportAs == "SpatRaster") {
      if (is.null(template)) stop("Matrix input lacks template attribute; cannot reconstruct raster.", call. = FALSE)
      return(matrix_to_raster(out, template, colnames(out)))
    }
    return(out)
  }

  if (inherits(dat, "SpatRaster")) {
    nms <- names(dat)
    vc_ids <- sub("_.*", "", nms)
    dr_ids <- sub(".*_", "", nms)

    out <- switch(by,
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

    if (exportAs == "matrix") {
      return(raster_to_matrix(out))
    }
    return(out)
  }

  # data.frame fallback (legacy)
  if (is.data.frame(dat)) {
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
  dat <- layers_extract(dat, layer_ids = layer_ids, drivers = drivers, vcs = vcs)

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

  vc_ids <- sub("_.*", "", layer_names)
  dr_ids <- sub(".*_", "", layer_names)
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
