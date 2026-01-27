make_layers_rasters <- function() {
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

test_that("layers_aggregate defaults to full aggregation (both)", {
  dat <- make_layers_rasters()
  ce_mat <- cea(dat$drivers, dat$vc, dat$sensitivity, exportAs = "matrix", engine = "matrix")

  res <- layers_aggregate(ce_mat)
  expect_true(inherits(res, "SpatRaster"))

  expected <- rowSums(ce_mat, na.rm = TRUE)
  out <- terra::values(res, mat = TRUE)
  expect_equal(out[, 1], expected)
})

test_that("layers_aggregate supports VC aggregation on matrix", {
  dat <- make_layers_rasters()
  ce_mat <- cea(dat$drivers, dat$vc, dat$sensitivity, exportAs = "matrix", engine = "matrix")

  res <- layers_aggregate(ce_mat, by = "drivers", exportAs = "matrix")

  expected <- cbind(
    cod = rowSums(ce_mat[, c("cod_shipping", "cod_climate")], na.rm = TRUE),
    salmon = rowSums(ce_mat[, c("salmon_shipping", "salmon_climate")], na.rm = TRUE)
  )
  attr(res, "template") <- NULL
  dimnames(res) <- NULL
  dimnames(expected) <- NULL
  expect_equal(res, expected)
})

test_that("layers_extract selects drivers on matrix input", {
  dat <- make_layers_rasters()
  ce_mat <- cea(dat$drivers, dat$vc, dat$sensitivity, exportAs = "matrix", engine = "matrix")

  res <- layers_extract(ce_mat, drivers = "shipping")
  expect_equal(colnames(res), c("cod_shipping", "salmon_shipping"))
})

test_that("layers_aggregate supports mean on matrix", {
  dat <- make_layers_rasters()
  ce_mat <- cea(dat$drivers, dat$vc, dat$sensitivity, exportAs = "matrix", engine = "matrix")

  res <- layers_aggregate(ce_mat, by = "both", fun = "mean", exportAs = "matrix")
  expected <- matrix(rowMeans(ce_mat, na.rm = TRUE), ncol = 1)
  attr(res, "template") <- NULL
  dimnames(res) <- NULL
  dimnames(expected) <- NULL
  expect_equal(res, expected)
})
