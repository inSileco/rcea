#' Cumulative exposure assessment
#'
#' Assessment of the exposure (i.e. overlap) between valued components and environmental drivers in the context of cumulative effects assessments
#'
#' @eval arguments(c("drivers","vc"))
#' @param exportAs string, the type of object that should be created, either a "list" or a "stars" object
#'
#' @export
#'
#' @examples
#' # Data
#' drivers <- rcea:::drivers
#' vc <- rcea:::vc
#'
#' # Exposure
#' (expo <- exposure(drivers, vc, "stars"))
#' plot(expo)
#' expo <- merge(expo, name = "vc") |>
#'   split("drivers")
#' plot(expo)
exposure <- function(drivers, vc, exportAs = c("list", "stars")) {

  #
  exportAs <- match.arg(exportAs)
  # Drivers
  dr_df <- as.data.frame(drivers) |>
    dplyr::select(-x, -y)

  # Valued components
  vc_df <- as.data.frame(vc) |>
    dplyr::select(-x, -y)

  # Exposure of valued components to drivers (Dj * VCi)
  dat <- apply(
    dr_df,
    MARGIN = 2,
    function(x) {
      sweep(vc_df, MARGIN = 1, x, `*`)
    }
  )

  # Return
  if (exportAs == "list") {
    dat
  } else if (exportAs == "stars") {
    make_stars(dat, drivers, vc)
  }
}
