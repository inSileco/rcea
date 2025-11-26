#' Cumulative effects assessments
#'
#' Assessment of cumulative effects using the Halpern et al. 2008 method.
#'
#' @eval arguments(c("drivers", "vc", "sensitivity"))
#' @param exportAs string, the type of object that should be created, either a "list" or a "stars" object.
#'
#' @export
#'
#' @examples
#' # Data
#' drivers <- rcea:::drivers
#' vc <- rcea:::vc
#' sensitivity <- rcea:::sensitivity
#'
#' # Species-scale effects
#' (halpern <- cea(drivers, vc, sensitivity, "stars"))
#' plot(halpern)
#' halpern <- merge(halpern, name = "vc") |>
#'   split("drivers")
#' plot(halpern)
#' # do not work
#' # get_cekm_cea(halpern, vc)
cea <- function(drivers, vc, sensitivity, exportAs = "list") {
  # needed as 
  # requireNamespace("stars", quietly = TRUE) 
  # Exposure
  dat <- exposure(drivers, vc)

  # Sensitivity
  nmDr <- names(dat)
  nmVC <- names(dat[[1]])
  sensitivity <- sensitivity[nmVC, nmDr]

  # Effect of drivers on valued components (D * VC * u)
  for (i in seq_len(length(dat))) {
    dat[[i]] <- sweep(dat[[i]], MARGIN = 2, sensitivity[, i], `*`)
  }

  # Return
  if (exportAs == "list") {
    dat
  } else if (exportAs == "stars") {
    xy <- sf::st_coordinates(vc)
    drNames <- names(drivers)
    make_stars(dat, drivers, vc)
  }
}

#' @describeIn cea get effects per km2
#' @param dat TODO
#' @export
get_cekm_cea <- function(dat, vc) {
  dat2 <- dat
  # CEA as data.frame
  dat <- as.data.frame(dat)

  # vc as data.frame
  vc_df <- as.data.frame(vc) |>
    dplyr::select(-x, -y)

  # Index of vc
  vc_index <- data.frame(
    vc = colnames(vc_df),
    vc_id = seq_len(ncol(vc_df))
  )

  # Calculate area, i.e. number of cells (assuming 1km2 grid cells)
  vc_df <- vc_df |>
    dplyr::mutate(id_cell = 1:dplyr::n()) |>
    tidyr::pivot_longer(cols = -c(id_cell), names_to = "vc", values_to = "presence") |>
    dplyr::group_by(vc) |>
    dplyr::summarise(km2 = sum(presence, na.rm = TRUE)) |>
    dplyr::left_join(vc_index, by = "vc") |>
    dplyr::ungroup()

  # Total effects per vc
  dat <- dplyr::select(dat, -x, -y) |>
    dplyr::group_by(drivers) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        \(x) sum(x, na.rm = TRUE)
      )
    ) |>
    dplyr::ungroup() |>
    tidyr::pivot_longer(cols = -c(drivers), names_to = "vc", values_to = "cea") |>
    dplyr::left_join(vc_df, by = "vc") |>
    dplyr::mutate(cea = cea / km2) |>
    dplyr::select(-vc, -km2) |>
    tidyr::pivot_wider(names_from = drivers, values_from = cea) |>
    dplyr::arrange(vc_id)

  # Return
  dat
}
