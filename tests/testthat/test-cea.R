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
  ce_r <- ce$data

  expect_true(inherits(ce, "rcea_result"))
  expect_true(inherits(ce_r, "SpatRaster"))
  expect_equal(terra::nlyr(ce_r), 4)
  expect_equal(names(ce_r), c("cod_shipping", "salmon_shipping", "cod_climate", "salmon_climate"))

  dr_mat <- terra::values(dat$drivers, mat = TRUE)
  vc_mat <- terra::values(dat$vc, mat = TRUE)

  exp1 <- vc_mat * dr_mat[, 1]
  exp2 <- vc_mat * dr_mat[, 2]
  eff1 <- sweep(exp1, 2, dat$sensitivity[, 1], `*`)
  eff2 <- sweep(exp2, 2, dat$sensitivity[, 2], `*`)
  expected <- cbind(eff1, eff2)
  colnames(expected) <- names(ce_r)

  out <- terra::values(ce_r, mat = TRUE)
  expect_equal(out, expected)
})

test_that("cea supports matrix export", {
  dat <- make_test_rasters()
  ce_mat <- cea(dat$drivers, dat$vc, dat$sensitivity, exportAs = "matrix", engine = "matrix")
  ce_data <- ce_mat$data

  expect_true(inherits(ce_mat, "rcea_result"))
  expect_true(is.matrix(ce_data))
  expect_equal(colnames(ce_data), c("cod_shipping", "salmon_shipping", "cod_climate", "salmon_climate"))
  expect_true(inherits(attr(ce_data, "template"), "SpatRaster"))

  ce_r <- rcea:::matrix_to_raster(ce_data, attr(ce_data, "template"), colnames(ce_data))
  attr(ce_data, "template") <- NULL
  attr(ce_data, "layer_map") <- NULL
  out <- unname(terra::values(ce_r, mat = TRUE))
  attr(out, "layer_map") <- NULL
  expect_equal(out, unname(ce_data))
})

test_that("cea reorders sensitivity by names", {
  dat <- make_test_rasters()
  sens_rev <- dat$sensitivity[c("salmon", "cod"), c("climate", "shipping")]

  ce1 <- cea(dat$drivers, dat$vc, dat$sensitivity, engine = "matrix")
  ce2 <- cea(dat$drivers, dat$vc, sens_rev, engine = "matrix")

  expect_equal(terra::values(ce1$data, mat = TRUE), terra::values(ce2$data, mat = TRUE))
})

test_that("cea terra engine matches matrix engine on small raster", {
  dat <- make_test_rasters()
  ce_mat <- cea(dat$drivers, dat$vc, dat$sensitivity, engine = "matrix")
  ce_terra <- cea(dat$drivers, dat$vc, dat$sensitivity, engine = "terra")

  expect_equal(terra::values(ce_mat$data, mat = TRUE), terra::values(ce_terra$data, mat = TRUE))
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
      c("cod", "salmon"),
      c("shipping", "climate")
    )
  )
  cube <- make_cube(catalog, sensitivity = sens)
  cube <- stack_layers(cube)

  ce1 <- cea(cube$stack$drivers, cube$stack$vc, cube$sensitivity, engine = "matrix")
  ce2 <- cea_cube(cube, engine = "matrix")

  expect_equal(terra::values(ce1$data, mat = TRUE), terra::values(ce2$data, mat = TRUE))
})

test_that("cea pair_by filters to matched temporal pairs", {
  tmpl <- terra::rast(ncols = 2, nrows = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:3857")

  d1 <- tmpl
  d2 <- tmpl
  v1 <- tmpl
  v2 <- tmpl
  terra::values(d1) <- c(1, 2, 3, 4)
  terra::values(d2) <- c(10, 20, 30, 40)
  terra::values(v1) <- c(1, 1, 1, 1)
  terra::values(v2) <- c(2, 2, 2, 2)

  drivers <- c(d1, d2)
  vc <- c(v1, v2)
  names(drivers) <- c("shipping_2020-01", "shipping_2020-02")
  names(vc) <- c("cod_2020-01", "cod_2020-02")

  driver_meta <- data.frame(
    layer = names(drivers),
    driver_id = "shipping",
    time = c("2020-01-15", "2020-02-15"),
    stringsAsFactors = FALSE
  )
  vc_meta <- data.frame(
    layer = names(vc),
    vc_id = "cod",
    time = c("2020-01-01", "2020-02-01"),
    stringsAsFactors = FALSE
  )

  sens <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list("cod", "shipping"))

  ce <- cea(
    drivers, vc, sens,
    exportAs = "matrix",
    engine = "matrix",
    driver_meta = driver_meta,
    vc_meta = vc_meta,
    pair_by = "month"
  )

  expect_equal(colnames(ce$data), c("cod_2020-01_shipping_2020-01", "cod_2020-02_shipping_2020-02"))
  expect_equal(ncol(ce$data), 2)

  expected <- cbind(
    terra::values(v1, mat = TRUE)[, 1] * terra::values(d1, mat = TRUE)[, 1] * 0.5,
    terra::values(v2, mat = TRUE)[, 1] * terra::values(d2, mat = TRUE)[, 1] * 0.5
  )
  out <- ce$data
  attr(out, "layer_map") <- NULL
  attr(out, "template") <- NULL
  expect_equal(unname(out), unname(expected))
})

test_that("cea pair_by errors when requested keys are missing", {
  dat <- make_test_rasters()
  expect_error(
    cea(dat$drivers, dat$vc, dat$sensitivity, engine = "matrix", pair_by = "month"),
    "requires both"
  )
})

test_that("cea pair_by broadcasts when one side has missing key values", {
  tmpl <- terra::rast(ncols = 2, nrows = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:3857")
  d1 <- tmpl
  d2 <- tmpl
  v1 <- tmpl
  v2 <- tmpl
  terra::values(d1) <- 1
  terra::values(d2) <- 2
  terra::values(v1) <- 3
  terra::values(v2) <- 4

  drivers <- c(d1, d2)
  vc <- c(v1, v2)
  names(drivers) <- c("shipping_jan", "coastal_static")
  names(vc) <- c("cod_jan", "cod_feb")

  driver_meta <- data.frame(
    layer = names(drivers),
    driver_id = c("shipping", "coastal"),
    month = c("2020-01", NA),
    stringsAsFactors = FALSE
  )
  vc_meta <- data.frame(
    layer = names(vc),
    vc_id = "cod",
    month = c("2020-01", "2020-02"),
    stringsAsFactors = FALSE
  )

  sens <- matrix(1, nrow = 1, ncol = 2, dimnames = list("cod", c("shipping", "coastal")))

  ce <- cea(
    drivers, vc, sens,
    exportAs = "matrix",
    engine = "matrix",
    driver_meta = driver_meta,
    vc_meta = vc_meta,
    pair_by = "month"
  )

  expect_setequal(
    colnames(ce$data),
    c("cod_jan_shipping_jan", "cod_jan_coastal_static", "cod_feb_coastal_static")
  )
  expect_equal(ncol(ce$data), 3)
})

test_that("cea pair_by strict drops pairs with missing key values", {
  tmpl <- terra::rast(ncols = 1, nrows = 1, xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "EPSG:3857")
  drivers <- c(tmpl, tmpl)
  vc <- c(tmpl, tmpl)
  names(drivers) <- c("shipping_jan", "coastal_static")
  names(vc) <- c("cod_jan", "cod_feb")
  terra::values(drivers) <- 1
  terra::values(vc) <- 1

  driver_meta <- data.frame(
    layer = names(drivers),
    driver_id = c("shipping", "coastal"),
    month = c("2020-01", NA),
    stringsAsFactors = FALSE
  )
  vc_meta <- data.frame(
    layer = names(vc),
    vc_id = "cod",
    month = c("2020-01", "2020-02"),
    stringsAsFactors = FALSE
  )

  sens <- matrix(1, nrow = 1, ncol = 2, dimnames = list("cod", c("shipping", "coastal")))
  ce <- cea(
    drivers, vc, sens,
    exportAs = "matrix",
    engine = "matrix",
    driver_meta = driver_meta,
    vc_meta = vc_meta,
    pair_by = "month",
    pair_missing = "strict"
  )

  expect_equal(colnames(ce$data), "cod_jan_shipping_jan")
  expect_equal(ncol(ce$data), 1)
})
