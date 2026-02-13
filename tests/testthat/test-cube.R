catalog_fixture <- function() {
  layers_path <- system.file("extdata/catalog/layers.csv", package = "rcea")
  groups_path <- system.file("extdata/catalog/groups.yaml", package = "rcea")
  if (layers_path == "" || groups_path == "") {
    testthat::skip("extdata catalog not available")
  }
  load_catalog(layers_path, groups_path)
}

test_that("make_cube stores catalog, aoi, and sensitivity", {
  catalog <- catalog_fixture()
  aoi_path <- system.file("extdata/aoi/aoi.geojson", package = "rcea")
  if (aoi_path == "") {
    testthat::skip("extdata AOI not available")
  }

  sens <- matrix(
    1,
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("cod", "salmon"),
      c("shipping", "climate")
    )
  )

  cube <- make_cube(catalog, aoi = aoi_path, sensitivity = sens)

  expect_true(inherits(cube, "rcea_cube"))
  expect_true(inherits(cube$aoi, "SpatVector"))
  expect_equal(cube$sensitivity, sens)
})

test_that("stack_layers builds drivers and vc stacks", {
  catalog <- catalog_fixture()
  cube <- make_cube(catalog)
  cube <- stack_layers(cube)

  expect_true(is.list(cube$stack))
  expect_true(inherits(cube$stack$drivers, "SpatRaster"))
  expect_true(inherits(cube$stack$vc, "SpatRaster"))

  expected_dr <- catalog$layers$layer_id[catalog$layers$type %in% "pressure"]
  expected_vc <- catalog$layers$layer_id[catalog$layers$type %in% "vc"]

  expect_setequal(names(cube$stack$drivers), expected_dr)
  expect_setequal(names(cube$stack$vc), expected_vc)
})

test_that("stack_layers supports explicit ids", {
  catalog <- catalog_fixture()
  expected_dr <- catalog$layers$layer_id[catalog$layers$type %in% "pressure"]
  expected_vc <- catalog$layers$layer_id[catalog$layers$type %in% "vc"]

  cube <- make_cube(catalog)
  cube <- stack_layers(cube, drivers_id = expected_dr[1], vc_id = expected_vc[1])

  expect_equal(terra::nlyr(cube$stack$drivers), 1)
  expect_equal(terra::nlyr(cube$stack$vc), 1)
  expect_equal(names(cube$stack$drivers), expected_dr[1])
  expect_equal(names(cube$stack$vc), expected_vc[1])
})

test_that("collect_layers returns stacks with or without collect", {
  catalog <- catalog_fixture()
  cube <- make_cube(catalog)
  cube <- stack_layers(cube)

  stacks <- collect_layers(cube, collect = FALSE)
  expect_true(is.list(stacks))
  expect_true(inherits(stacks$drivers, "SpatRaster"))
  expect_true(inherits(stacks$vc, "SpatRaster"))

  stacks2 <- collect_layers(cube, collect = TRUE)
  expect_true(is.list(stacks2))
  expect_true(inherits(stacks2$drivers, "SpatRaster"))
  expect_true(inherits(stacks2$vc, "SpatRaster"))
})

test_that("stack_layers supports collect argument", {
  catalog <- catalog_fixture()
  cube <- make_cube(catalog)
  cube <- stack_layers(cube, collect = TRUE)

  expect_true(is.list(cube$stack))
  expect_true(inherits(cube$stack$drivers, "SpatRaster"))
  expect_true(inherits(cube$stack$vc, "SpatRaster"))
})

test_that("stack_layers aligns to template when requested", {
  td <- tempdir()
  r1 <- terra::rast(ncols = 4, nrows = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4, crs = "EPSG:3857")
  r2 <- terra::rast(ncols = 8, nrows = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 8, crs = "EPSG:3857")
  r3 <- terra::rast(ncols = 6, nrows = 6, xmin = 1, xmax = 7, ymin = 1, ymax = 7, crs = "EPSG:3857")
  terra::values(r1) <- 1
  terra::values(r2) <- 2
  terra::values(r3) <- 3

  p1 <- file.path(td, "dr1.tif")
  p2 <- file.path(td, "dr2.tif")
  p3 <- file.path(td, "vc1.tif")
  terra::writeRaster(r1, p1, overwrite = TRUE)
  terra::writeRaster(r2, p2, overwrite = TRUE)
  terra::writeRaster(r3, p3, overwrite = TRUE)

  catalog <- list(
    layers = data.frame(
      layer_id = c("shipping", "climate", "cod"),
      path = c(p1, p2, p3),
      type = c("pressure", "pressure", "vc"),
      stringsAsFactors = FALSE
    ),
    groups = list()
  )
  cube <- make_cube(catalog)
  tmpl <- terra::rast(ncols = 5, nrows = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "EPSG:3857")
  cube <- stack_layers(cube, align = "template", template = tmpl)

  expect_true(isTRUE(terra::compareGeom(cube$stack$drivers, tmpl, stopOnError = FALSE)))
  expect_true(isTRUE(terra::compareGeom(cube$stack$vc, tmpl, stopOnError = FALSE)))
})

test_that("stack_layers errors on misaligned rasters when align = 'error'", {
  td <- tempdir()
  r1 <- terra::rast(ncols = 4, nrows = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4, crs = "EPSG:3857")
  r2 <- terra::rast(ncols = 8, nrows = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 8, crs = "EPSG:3857")
  r3 <- terra::rast(ncols = 6, nrows = 6, xmin = 1, xmax = 7, ymin = 1, ymax = 7, crs = "EPSG:3857")
  terra::values(r1) <- 1
  terra::values(r2) <- 2
  terra::values(r3) <- 3

  p1 <- file.path(td, "dr1_error.tif")
  p2 <- file.path(td, "dr2_error.tif")
  p3 <- file.path(td, "vc1_error.tif")
  terra::writeRaster(r1, p1, overwrite = TRUE)
  terra::writeRaster(r2, p2, overwrite = TRUE)
  terra::writeRaster(r3, p3, overwrite = TRUE)

  catalog <- list(
    layers = data.frame(
      layer_id = c("shipping", "climate", "cod"),
      path = c(p1, p2, p3),
      type = c("pressure", "pressure", "vc"),
      stringsAsFactors = FALSE
    ),
    groups = list()
  )
  cube <- make_cube(catalog)

  expect_error(
    stack_layers(cube, align = "error"),
    "extents do not match"
  )
})

test_that("cea_cube carries catalog time metadata into layer_map grouping", {
  td <- tempdir()
  r1 <- terra::rast(ncols = 2, nrows = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:3857")
  r2 <- r1
  r3 <- r1
  terra::values(r1) <- 1
  terra::values(r2) <- 2
  terra::values(r3) <- 1

  p1 <- file.path(td, "shipping_2020_01.tif")
  p2 <- file.path(td, "shipping_2020_02.tif")
  p3 <- file.path(td, "cod.tif")
  terra::writeRaster(r1, p1, overwrite = TRUE)
  terra::writeRaster(r2, p2, overwrite = TRUE)
  terra::writeRaster(r3, p3, overwrite = TRUE)

  catalog <- list(
    layers = data.frame(
      layer_id = c("shipping_2020_01", "shipping_2020_02", "cod"),
      path = c(p1, p2, p3),
      type = c("pressure", "pressure", "vc"),
      group = c("shipping", "shipping", "fish"),
      driver_id = c("shipping", "shipping", NA),
      vc_id = c(NA, NA, "cod"),
      time = c("2020-01-15", "2020-02-15", NA),
      stringsAsFactors = FALSE
    ),
    groups = list()
  )

  cube <- make_cube(catalog)
  cube <- stack_layers(cube)

  sens <- matrix(1, nrow = 1, ncol = 2, dimnames = list("cod", c("shipping_2020_01", "shipping_2020_02")))
  ce <- cea_cube(cube, sensitivity = sens, exportAs = "matrix", engine = "matrix")

  expect_true("driver_month" %in% names(ce$layer_map))
  grouped <- layers_aggregate(ce, group_by = "driver_month", exportAs = "matrix")
  expect_equal(ncol(grouped), 2)
  expect_setequal(colnames(grouped), c("2020-01", "2020-02"))
})

test_that("preprocess_layers applies global transforms and normalization", {
  td <- tempdir()
  dr <- terra::rast(ncols = 2, nrows = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:3857")
  vc <- dr
  terra::values(dr) <- c(0, 1, 3, 7)
  terra::values(vc) <- c(1, 1, 1, 1)

  p_dr <- file.path(td, "pre_dr.tif")
  p_vc <- file.path(td, "pre_vc.tif")
  terra::writeRaster(dr, p_dr, overwrite = TRUE)
  terra::writeRaster(vc, p_vc, overwrite = TRUE)

  catalog <- list(
    layers = data.frame(
      layer_id = c("shipping", "cod"),
      path = c(p_dr, p_vc),
      type = c("pressure", "vc"),
      group = c("shipping", "fish"),
      stringsAsFactors = FALSE
    ),
    groups = list()
  )

  cube <- make_cube(catalog)
  cube <- stack_layers(cube)
  cube <- preprocess_layers(cube, drivers_transform = "log1p", drivers_normalize = "minmax")

  out_dr <- as.numeric(terra::values(cube$stack$drivers[[1]], mat = TRUE)[, 1])
  expected <- log1p(c(0, 1, 3, 7))
  expected <- (expected - min(expected)) / (max(expected) - min(expected))
  expect_equal(out_dr, expected)

  out_vc <- as.numeric(terra::values(cube$stack$vc[[1]], mat = TRUE)[, 1])
  expect_equal(out_vc, rep(1, 4))
})

test_that("preprocess_layers allows catalog overrides", {
  td <- tempdir()
  dr <- terra::rast(ncols = 1, nrows = 2, xmin = 0, xmax = 1, ymin = 0, ymax = 2, crs = "EPSG:3857")
  vc <- dr
  terra::values(dr) <- c(9, 16)
  terra::values(vc) <- c(4, 9)

  p_dr <- file.path(td, "pre2_dr.tif")
  p_vc <- file.path(td, "pre2_vc.tif")
  terra::writeRaster(dr, p_dr, overwrite = TRUE)
  terra::writeRaster(vc, p_vc, overwrite = TRUE)

  catalog <- list(
    layers = data.frame(
      layer_id = c("shipping", "cod"),
      path = c(p_dr, p_vc),
      type = c("pressure", "vc"),
      group = c("shipping", "fish"),
      transform = c("none", "sqrt"),
      normalize = c("max", "none"),
      stringsAsFactors = FALSE
    ),
    groups = list()
  )

  cube <- make_cube(catalog)
  cube <- stack_layers(cube)
  cube <- preprocess_layers(
    cube,
    drivers_transform = "log1p",
    drivers_normalize = "none",
    vc_transform = "none",
    vc_normalize = "none",
    use_catalog_overrides = TRUE
  )

  out_dr <- as.numeric(terra::values(cube$stack$drivers[[1]], mat = TRUE)[, 1])
  expect_equal(out_dr, c(9, 16) / 16)

  out_vc <- as.numeric(terra::values(cube$stack$vc[[1]], mat = TRUE)[, 1])
  expect_equal(out_vc, sqrt(c(4, 9)))

  dr_meta <- attr(cube$stack$drivers, "layer_meta")
  vc_meta <- attr(cube$stack$vc, "layer_meta")
  expect_equal(dr_meta$applied_transform[dr_meta$layer == "shipping"], "none")
  expect_equal(dr_meta$applied_normalize[dr_meta$layer == "shipping"], "max")
  expect_equal(vc_meta$applied_transform[vc_meta$layer == "cod"], "sqrt")
})
