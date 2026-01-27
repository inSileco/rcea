make_cekm_rasters <- function() {
  tmpl <- terra::rast(ncols = 2, nrows = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:3857")

  d1 <- tmpl
  d2 <- tmpl
  terra::values(d1) <- c(1, 2, 3, 4)
  terra::values(d2) <- c(5, 6, 7, 8)
  drivers <- c(d1, d2)
  names(drivers) <- c("shipping", "climate")

  v1 <- tmpl
  v2 <- tmpl
  terra::values(v1) <- c(1, 0, 1, 0)
  terra::values(v2) <- c(0.5, 0.5, 0.5, 0.5)
  vc <- c(v1, v2)
  names(vc) <- c("cod", "salmon")

  sensitivity <- matrix(
    c(0.8, 0.5,
      0.2, 0.7),
    nrow = 2,
    dimnames = list(c("cod", "salmon"), c("shipping", "climate"))
  )

  list(drivers = drivers, vc = vc, sensitivity = sensitivity)
}

test_that("layers_per_area computes per-area effects", {
  dat <- make_cekm_rasters()
  ce_mat <- cea(dat$drivers, dat$vc, dat$sensitivity, exportAs = "matrix", engine = "matrix")

  res <- layers_per_area(ce_mat, dat$vc)

  layer_sums <- colSums(ce_mat, na.rm = TRUE)
  area_r <- terra::cellSize(dat$vc[[1]], unit = "km")
  area_vec <- terra::values(area_r, mat = FALSE)
  vc_area <- colSums(terra::values(dat$vc, mat = TRUE) * area_vec, na.rm = TRUE)

  expected <- data.frame(
    vc = c("cod", "salmon"),
    shipping = c(unname(layer_sums["cod_shipping"]), unname(layer_sums["salmon_shipping"])) / vc_area,
    climate = c(unname(layer_sums["cod_climate"]), unname(layer_sums["salmon_climate"])) / vc_area,
    stringsAsFactors = FALSE
  )

  rownames(expected) <- NULL
  res_df <- as.data.frame(res)
  driver_cols <- setdiff(names(res_df), "vc")
  res_df[driver_cols] <- lapply(res_df[driver_cols], units::drop_units)
  expect_equal(res_df, expected)
})
