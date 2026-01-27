make_test_rasters <- function() {
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

test_that("cea matrix engine returns expected layers and values", {
  dat <- make_test_rasters()
  ce <- cea(dat$drivers, dat$vc, dat$sensitivity, engine = "matrix")

  expect_true(inherits(ce, "SpatRaster"))
  expect_equal(terra::nlyr(ce), 4)
  expect_equal(names(ce), c("cod_shipping", "salmon_shipping", "cod_climate", "salmon_climate"))

  dr_mat <- terra::values(dat$drivers, mat = TRUE)
  vc_mat <- terra::values(dat$vc, mat = TRUE)

  exp1 <- vc_mat * dr_mat[, 1]
  exp2 <- vc_mat * dr_mat[, 2]
  eff1 <- sweep(exp1, 2, dat$sensitivity[, 1], `*`)
  eff2 <- sweep(exp2, 2, dat$sensitivity[, 2], `*`)
  expected <- cbind(eff1, eff2)
  colnames(expected) <- names(ce)

  out <- terra::values(ce, mat = TRUE)
  expect_equal(out, expected)
})

test_that("cea reorders sensitivity by names", {
  dat <- make_test_rasters()
  sens_rev <- dat$sensitivity[c("salmon", "cod"), c("climate", "shipping")]

  ce1 <- cea(dat$drivers, dat$vc, dat$sensitivity, engine = "matrix")
  ce2 <- cea(dat$drivers, dat$vc, sens_rev, engine = "matrix")

  expect_equal(terra::values(ce1, mat = TRUE), terra::values(ce2, mat = TRUE))
})

test_that("cea terra engine matches matrix engine on small raster", {
  dat <- make_test_rasters()
  ce_mat <- cea(dat$drivers, dat$vc, dat$sensitivity, engine = "matrix")
  ce_terra <- cea(dat$drivers, dat$vc, dat$sensitivity, engine = "terra")

  expect_equal(terra::values(ce_mat, mat = TRUE), terra::values(ce_terra, mat = TRUE))
})

test_that("cea errors on misaligned rasters when align = 'error'", {
  dat <- make_test_rasters()
  vc2 <- terra::rast(ncols = 2, nrows = 2, xmin = 0, xmax = 4, ymin = 0, ymax = 4, crs = "EPSG:3857")
  terra::values(vc2) <- c(1, 1, 1, 1)
  vc2 <- c(vc2, vc2)
  names(vc2) <- names(dat$vc)

  expect_error(
    cea(dat$drivers, vc2, dat$sensitivity, align = "error", engine = "matrix"),
    "Drivers and vc differ"
  )
})

test_that("cea_cube uses cube stacks and sensitivity", {
  layers_path <- system.file("extdata/catalog/layers.csv", package = "rcea")
  groups_path <- system.file("extdata/catalog/groups.yaml", package = "rcea")
  if (layers_path == "" || groups_path == "") {
    skip("extdata catalog not available")
  }

  catalog <- load_catalog(layers_path, groups_path)
  sens <- matrix(
    1,
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("vc_cod", "vc_salmon"),
      c("pressure_shipping", "pressure_climate")
    )
  )
  cube <- make_cube(catalog, sensitivity = sens)
  cube <- stack_layers(cube)

  ce1 <- cea(cube$stack$drivers, cube$stack$vc, cube$sensitivity, engine = "matrix")
  ce2 <- cea_cube(cube, engine = "matrix")

  expect_equal(terra::values(ce1, mat = TRUE), terra::values(ce2, mat = TRUE))
})
