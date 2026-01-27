#' Normalize AOI input
#'
#' Accepts sf object, bbox list, or WKT string and returns a terra SpatVector.
#'
#' @param aoi sf object, bbox list with xmin/xmax/ymin/ymax/crs, or WKT string
#' @return SpatVector
#' @examples
#' bbox_path <- system.file("extdata/aoi/bbox.yml", package = "rcea")
#' bbox <- yaml::read_yaml(bbox_path)
#' aoi1 <- normalize_aoi(bbox)
#' aoi1
#'
#' geojson_path <- system.file("extdata/aoi/aoi.geojson", package = "rcea")
#' aoi2 <- normalize_aoi(geojson_path)
#' aoi2
#' @export
normalize_aoi <- function(aoi) {
  if (inherits(aoi, "SpatVector")) {
    return(aoi)
  }
  if (inherits(aoi, "sf")) {
    return(terra::vect(aoi))
  }

  if (is.list(aoi) && all(c("xmin", "xmax", "ymin", "ymax") %in% names(aoi))) {
    v <- terra::vect(
      matrix(
        c(
          aoi$xmin, aoi$ymin,
          aoi$xmax, aoi$ymin,
          aoi$xmax, aoi$ymax,
          aoi$xmin, aoi$ymax,
          aoi$xmin, aoi$ymin
        ),
        ncol = 2, byrow = TRUE
      ),
      type = "polygons",
      crs = aoi$crs %||% NA
    )
    return(v)
  }

  if (is.character(aoi) && length(aoi) == 1) {
    v <- try(terra::vect(aoi), silent = TRUE)
    if (!inherits(v, "try-error")) {
      return(v)
    }
  }

  stop("Unsupported AOI format. Provide sf, bbox list, WKT/GeoJSON string, or SpatVector.", call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
