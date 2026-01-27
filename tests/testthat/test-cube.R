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
      c("vc_cod", "vc_salmon"),
      c("pressure_shipping", "pressure_climate")
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
