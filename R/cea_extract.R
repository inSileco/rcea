#' Extract attributes and aggregate cumulative effects assessment results
#'
#' Function to extract specific drivers and vc of interest from the drivers, valued component, exposure and cumulative effects assessments results, and to aggregate data over drivers, valued components, or drivers and valued components.
#'
#' @param dat stars object, either the drivers, valued components, exposure, cumulative effects assessment or network-sacale cumulative effects assessment results
#' @param dr_sel string, name of drivers to extract 
#' @param vc_sel string, name of valued components to extract
#' @param cumul_fun function to apply on the stars object, one of "drivers" for the cumulative effects of each drivers on all value components, "vc" for the cumulative effects of all drivers on each valued component, "full" for the cumulative effects of all drivers on all valued components, "footprint" for the cumulative footprint of drivers or valued components, and "none" to keep data as is, i.e. to extract specific drivers and valued components from the assessment results.
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
#' drivers <- terra::rast(drv_paths); names(drivers) <- c("shipping", "climate")
#' vc <- terra::rast(vc_paths); names(vc) <- c("cod", "salmon")
#' sens <- matrix(
#'   c(0.8, 0.5,
#'     0.2, 0.7),
#'   nrow = 2,
#'   dimnames = list(c("cod", "salmon"), c("shipping", "climate"))
#' )
#' ce <- cea(drivers, vc, sens, exportAs = "SpatRaster")
#' # Footprint across all layers
#' cea_extract(ce, cumul_fun = "footprint")
#' # Sum per VC
#' cea_extract(ce, cumul_fun = "vc")
#'
#' @export
cea_extract <- function(dat, dr_sel = NULL, vc_sel = NULL, cumul_fun = "none") {
  # Select drivers and vc to extract
  dat <- select_attr(dat, dr_sel, vc_sel)

  # Apply relevant aggregations, if applicable
  dat <- aggr(dat, cumul_fun)

  # Return
  dat
}


select_attr <- function(dat, dr_sel = NULL, vc_sel = NULL) {
  # Select drivers, if applicable
  if (!is.null(dr_sel)) {
    if (inherits(dat, "SpatRaster")) {
      nms <- names(dat)
      # direct match
      idx <- nms %in% dr_sel
      if (!any(idx)) {
        idx <- vapply(nms, function(x) sub(".*_", "", x) %in% dr_sel, logical(1))
      }
      dat <- dat[[which(idx)]]
    } else {
      if ("drivers" %in% colnames(dat)) {
        dat <- dplyr::filter(dat, drivers %in% dr_sel)
      } else {
        dat <- dplyr::select(dat, x, y, dplyr::all_of(dr_sel))
      }
    }
  }
  
  # Select valued components, if applicable
  if (!is.null(vc_sel)) {
    if (inherits(dat, "SpatRaster")) {
      nms <- names(dat)
      idx <- nms %in% vc_sel
      if (!any(idx)) {
        idx <- vapply(nms, function(x) sub("_.*", "", x) %in% vc_sel, logical(1))
      }
      dat <- dat[[which(idx)]]
    } else {
      dat <- dplyr::select(dat, x, y, dplyr::any_of("drivers"), dplyr::all_of(vc_sel))
    }
  }
  
  # Return
  dat
}

cumul <- function(dat) {
    stars::st_apply(dat, c("x","y"), sum, na.rm = TRUE)    
}

cumul_vc <- function(dat) {
  if (inherits(dat, "SpatRaster")) {
    vc_ids <- unique(sub("_.*", "", names(dat)))
    out <- lapply(vc_ids, function(v) {
      terra::app(dat[[grep(paste0("^", v, "_"), names(dat))]], sum, na.rm = TRUE)
    })
    names(out) <- vc_ids
    terra::rast(out)
  } else {
    dat |>
    dplyr::select(-drivers) |>
    dplyr::group_by(x,y) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        \(x) sum(x, na.rm = TRUE)
      )
    ) |>
    dplyr::ungroup()
  }
}

cumul_drivers <- function(dat) {
  if (inherits(dat, "SpatRaster")) {
    dr_ids <- unique(sub(".*_", "", names(dat)))
    out <- lapply(dr_ids, function(d) {
      terra::app(dat[[grep(paste0("_", d, "$"), names(dat))]], sum, na.rm = TRUE)
    })
    names(out) <- dr_ids
    terra::rast(out)
  } else {
    dat |>
    dplyr::mutate(value = rowSums(dplyr::pick(-x,-y,-drivers), na.rm = TRUE)) |>
    dplyr::select(x,y,drivers,value) |>
    tidyr::pivot_wider(names_from = "drivers", values_from = "value")
  }
}

cumul_full <- function(dat) {
  if (inherits(dat, "SpatRaster")) {
    terra::app(dat, sum, na.rm = TRUE)
  } else {
    dat |>
    dplyr::mutate(value = rowSums(dplyr::pick(-x,-y,-drivers), na.rm = TRUE)) |>
    dplyr::select(x,y,drivers,value) |>
    tidyr::pivot_wider(names_from = "drivers", values_from = "value") |>
    dplyr::mutate(cumulative_effects = rowSums(dplyr::pick(-x,-y), na.rm = TRUE)) |>
    dplyr::select(x,y,cumulative_effects)
  }
}

cumul_footprint <- function(dat) {
  if (inherits(dat, "SpatRaster")) {
    terra::app(dat, sum, na.rm = TRUE)
  } else {
    dat |>
    dplyr::mutate(cumulative_footprint = rowSums(dplyr::pick(-x,-y), na.rm = TRUE)) |>
    dplyr::select(x,y,cumulative_footprint) 
  }
}

aggr <- function(dat, cumul_fun) {
  switch(
   cumul_fun, 
   "drivers" = cumul_drivers(dat),
   "vc" = cumul_vc(dat),
   "full" = cumul_full(dat),
   "footprint" = cumul_footprint(dat),
   "none" = dat
  )
}
