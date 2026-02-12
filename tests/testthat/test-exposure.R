make_exposure_rasters <- function() {
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

  list(drivers = drivers, vc = vc)
}

test_that("exposure returns expected layers and values", {
  dat <- make_exposure_rasters()
  expo <- exposure(dat$drivers, dat$vc, exportAs = "SpatRaster")
  expo_r <- expo$data

  expect_true(inherits(expo, "rcea_result"))
  expect_true(inherits(expo_r, "SpatRaster"))
  expect_equal(names(expo_r), c("cod_shipping", "salmon_shipping", "cod_climate", "salmon_climate"))

  dr_mat <- terra::values(dat$drivers, mat = TRUE)
  vc_mat <- terra::values(dat$vc, mat = TRUE)

  exp1 <- vc_mat * dr_mat[, 1]
  exp2 <- vc_mat * dr_mat[, 2]
  expected <- cbind(exp1, exp2)
  colnames(expected) <- names(expo_r)

  out <- terra::values(expo_r, mat = TRUE)
  expect_equal(out, expected)
})

test_that("exposure supports matrix export", {
  dat <- make_exposure_rasters()
  expo_mat <- exposure(dat$drivers, dat$vc, exportAs = "matrix")
  expo_data <- expo_mat$data

  expect_true(inherits(expo_mat, "rcea_result"))
  expect_true(is.matrix(expo_data))
  expect_equal(colnames(expo_data), c("cod_shipping", "salmon_shipping", "cod_climate", "salmon_climate"))
  expect_true(inherits(attr(expo_data, "template"), "SpatRaster"))
})

test_that("exposure errors on misaligned rasters when align = 'error'", {
  dat <- make_exposure_rasters()
  vc2 <- terra::rast(ncols = 2, nrows = 2, xmin = 0, xmax = 4, ymin = 0, ymax = 4, crs = "EPSG:3857")
  terra::values(vc2) <- c(1, 1, 1, 1)
  vc2 <- c(vc2, vc2)
  names(vc2) <- names(dat$vc)

  expect_error(
    exposure(dat$drivers, vc2, align = "error"),
    "Drivers and vc differ"
  )
})

test_that("exposure aligns with template", {
  dat <- make_exposure_rasters()
  vc2 <- terra::rast(ncols = 4, nrows = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4, crs = "EPSG:3857")
  terra::values(vc2) <- rep(1, terra::ncell(vc2))
  vc2 <- c(vc2, vc2)
  names(vc2) <- names(dat$vc)

  expo <- exposure(dat$drivers, vc2, align = "template", template = dat$drivers[[1]])
  expo_r <- expo$data

  expect_equal(terra::ncell(expo_r), terra::ncell(dat$drivers[[1]]))
  expect_equal(terra::res(expo_r), terra::res(dat$drivers[[1]]))
  expect_equal(as.vector(terra::ext(expo_r)), as.vector(terra::ext(dat$drivers[[1]])))
})

test_that("exposure_cube matches exposure", {
  layers_path <- system.file("extdata/catalog/layers.csv", package = "rcea")
  groups_path <- system.file("extdata/catalog/groups.yaml", package = "rcea")
  if (layers_path == "" || groups_path == "") {
    skip("extdata catalog not available")
  }

  catalog <- load_catalog(layers_path, groups_path)
  cube <- make_cube(catalog)
  cube <- stack_layers(cube)

  expo1 <- exposure(cube$stack$drivers, cube$stack$vc)
  expo2 <- exposure_cube(cube)

  expect_equal(terra::values(expo1$data, mat = TRUE), terra::values(expo2$data, mat = TRUE))
})
