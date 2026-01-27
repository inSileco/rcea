test_that("load_catalog validates structure and paths", {
  tmp <- tempfile("catalog_")
  dir.create(tmp)

  tif <- file.path(tmp, "dummy.tif")
  file.create(tif)

  layers_path <- file.path(tmp, "layers.csv")
  groups_path <- file.path(tmp, "groups.yaml")

  layers <- data.frame(
    layer_id = "pressure_shipping",
    path = tif,
    type = "pressure",
    group = "shipping",
    units = "idx",
    crs = "EPSG:3857",
    res_x = 1,
    res_y = 1,
    stringsAsFactors = FALSE
  )
  write.csv(layers, layers_path, row.names = FALSE)

  groups <- list(pressure = list(shipping = list(members = "pressure_shipping")))
  yaml::write_yaml(groups, groups_path)

  out <- load_catalog(layers_path, groups_path)
  expect_true(is.list(out))
  expect_true(all(c("layers", "groups") %in% names(out)))
  expect_equal(out$layers$layer_id, "pressure_shipping")
})

test_that("validate_catalog flags bad inputs", {
  tmp <- tempfile("catalog_")
  dir.create(tmp)

  tif <- file.path(tmp, "dummy.tif")
  file.create(tif)

  base_layers <- data.frame(
    layer_id = "l1",
    path = tif,
    type = "pressure",
    group = "shipping",
    units = "idx",
    crs = "EPSG:3857",
    res_x = 1,
    res_y = 1,
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

  bad_res <- base_layers
  bad_res$res_x <- NA
  expect_error(validate_catalog(bad_res, groups), "res_x/res_y must be finite")
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
