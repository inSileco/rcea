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
