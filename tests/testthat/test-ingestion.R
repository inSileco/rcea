test_that("load_catalog validates structure and paths", {
  tmp <- tempfile("catalog_")
  dir.create(tmp)

  tif <- file.path(tmp, "dummy.tif")
  r <- terra::rast(ncols = 2, nrows = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:3857")
  r[] <- 1
  terra::units(r) <- "idx"
  terra::writeRaster(r, tif, overwrite = TRUE)

  layers_path <- file.path(tmp, "layers.csv")
  groups_path <- file.path(tmp, "groups.yaml")

  layers <- data.frame(
    layer_id = "shipping",
    path = tif,
    type = "pressure",
    group = "shipping",
    stringsAsFactors = FALSE
  )
  write.csv(layers, layers_path, row.names = FALSE)

  groups <- list(pressure = list(shipping = list(members = "shipping")))
  yaml::write_yaml(groups, groups_path)

  out <- load_catalog(layers_path, groups_path)
  expect_true(is.list(out))
  expect_true(all(c("layers", "groups") %in% names(out)))
  expect_equal(out$layers$layer_id, "shipping")
  expect_equal(out$layers$crs, "EPSG:3857")
  expect_equal(out$layers$res_x, 1)
  expect_equal(out$layers$res_y, 1)
  expect_equal(out$layers$units, "idx")
})

test_that("validate_catalog flags bad inputs", {
  tmp <- tempfile("catalog_")
  dir.create(tmp)

  tif <- file.path(tmp, "dummy.tif")
  r <- terra::rast(ncols = 2, nrows = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:3857")
  r[] <- 1
  terra::writeRaster(r, tif, overwrite = TRUE)

  base_layers <- data.frame(
    layer_id = "l1",
    path = tif,
    type = "pressure",
    group = "shipping",
    stringsAsFactors = FALSE
  )

  groups <- list(pressure = list(shipping = list(members = "l1")))

  bad_type <- base_layers
  bad_type$type <- "bad"
  expect_error(validate_catalog(bad_type, groups), "Invalid type values")

  bad_group <- base_layers
  bad_group$group <- "nope"
  expect_error(validate_catalog(bad_group, groups), "Groups not defined")

  bad_path <- base_layers
  bad_path$path <- file.path(tmp, "missing.tif")
  expect_error(validate_catalog(bad_path, groups), "Layer paths not found")

  dup <- rbind(base_layers, base_layers)
  expect_error(validate_catalog(dup, groups), "Duplicate layer_id")

  bad_raster <- base_layers
  bad_raster$path <- file.path(tmp, "not_a_raster.tif")
  writeLines("not a raster", bad_raster$path)
  expect_error(validate_catalog(bad_raster, groups), "Unable to read raster")
})

test_that("normalize_aoi supports bbox, WKT, and SpatVector", {
  bbox <- list(xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "EPSG:4326")
  v1 <- normalize_aoi(bbox)
  expect_true(inherits(v1, "SpatVector"))

  wkt <- "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"
  v2 <- normalize_aoi(wkt)
  expect_true(inherits(v2, "SpatVector"))

  v3 <- normalize_aoi(v1)
  expect_true(inherits(v3, "SpatVector"))
})

test_that("normalize_aoi errors on unsupported input", {
  expect_error(normalize_aoi(123), "Unsupported AOI format")
})
