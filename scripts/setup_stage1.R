#!/usr/bin/env Rscript

# Script to scaffold Stage 1 assets:
# - catalog/layers.csv
# - catalog/groups.yaml
# - config.yml
# - synthetic rasters in inst/extdata/rasters
# - example AOI polygon and bbox in inst/extdata/aoi
#
# Run from the package root. Edit the parameters below to tweak the scaffold.

suppressPackageStartupMessages({
  library(terra)
  library(yaml)
})

# ------------------------------------------------------------------------------
# Parameters
proj_crs <- "EPSG:3857"
rast_ncol <- 20
rast_nrow <- 20
rast_extent <- c(xmin = 0, xmax = 20000, ymin = 0, ymax = 20000) # metres
set.seed(42)

# ------------------------------------------------------------------------------
root <- normalizePath(".", winslash = "/")
paths <- list(
  catalog_dir   = file.path(root, "catalog"),
  layers_csv    = file.path(root, "catalog", "layers.csv"),
  groups_yaml   = file.path(root, "catalog", "groups.yaml"),
  config_yml    = file.path(root, "config.yml"),
  raster_dir    = file.path(root, "inst", "extdata", "rasters"),
  aoi_dir       = file.path(root, "inst", "extdata", "aoi")
)
dir.create(paths$catalog_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$raster_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$aoi_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Template raster
r_template <- rast(
  ncols = rast_ncol,
  nrows = rast_nrow,
  xmin = rast_extent["xmin"],
  xmax = rast_extent["xmax"],
  ymin = rast_extent["ymin"],
  ymax = rast_extent["ymax"],
  crs  = proj_crs
)

mk_grad <- function(r, slope = 1, noise = 0.1) {
  vals <- matrix(rep(seq(0, slope, length.out = nrow(r)), each = ncol(r)), nrow = nrow(r))
  vals <- vals + matrix(runif(ncell(r), -noise, noise), nrow = nrow(r))
  r[] <- vals
  r
}

mk_binary <- function(r, p = 0.5) {
  r[] <- as.numeric(runif(ncell(r)) < p)
  r
}

# Synthetic layers
pressure_shipping <- mk_grad(r_template, slope = 1, noise = 0.2)
pressure_climate  <- mk_grad(r_template, slope = 0.5, noise = 0.15)
vc_cod            <- mk_binary(r_template, p = 0.35)
vc_salmon         <- mk_binary(r_template, p = 0.25)

mask_center <- r_template
mask_center[] <- NA
center_xy <- cbind(
  mean(c(rast_extent["xmin"], rast_extent["xmax"])),
  mean(c(rast_extent["ymin"], rast_extent["ymax"])))
center_cell <- terra::cellFromXY(mask_center, center_xy)
mask_center[center_cell] <- 1
dist_r <- terra::distance(mask_center)
mask_vals <- as.numeric(dist_r <= 12000)
mask_vals[is.na(mask_vals)] <- 0
mask_study_area <- r_template
mask_study_area[] <- mask_vals

# Write rasters
wopt <- list(filetype = "GTiff", gdal = c("COMPRESS=LZW"), datatype = "FLT4S", overwrite = TRUE)
terra::writeRaster(pressure_shipping, file.path(paths$raster_dir, "pressure_shipping.tif"), wopt = wopt)
terra::writeRaster(pressure_climate,  file.path(paths$raster_dir, "pressure_climate.tif"),  wopt = wopt)
terra::writeRaster(vc_cod,            file.path(paths$raster_dir, "vc_cod.tif"),            wopt = wopt)
terra::writeRaster(vc_salmon,         file.path(paths$raster_dir, "vc_salmon.tif"),         wopt = wopt)
terra::writeRaster(mask_study_area,   file.path(paths$raster_dir, "mask_study_area.tif"),   wopt = wopt)

# ------------------------------------------------------------------------------
# Catalog layers.csv
res_vals <- terra::res(r_template)
layers <- data.frame(
  layer_id       = c("shipping", "climate", "cod", "salmon", "mask_study_area"),
  path           = file.path("inst", "extdata", "rasters",
                             c("pressure_shipping.tif", "pressure_climate.tif", "vc_cod.tif", "vc_salmon.tif", "mask_study_area.tif")),
  type           = c("pressure", "pressure", "vc", "vc", "aux"),
  group          = c("shipping", "climate", "fish", "fish", "mask"),
  units          = c("index", "index", "presence", "presence", "mask"),
  crs            = proj_crs,
  res_x          = res_vals[1],
  res_y          = res_vals[2],
  time           = NA_character_,
  tags           = c("pressure;shipping", "pressure;climate", "vc;marine", "vc;marine", "aux;mask"),
  requires_auth  = FALSE,
  version        = "0.1",
  checksum       = NA_character_,
  nodata         = NA_real_,
  datatype       = "FLOAT32",
  source         = "synthetic"
)
write.csv(layers, paths$layers_csv, row.names = FALSE)

# groups.yaml
groups <- list(
  pressure = list(
    shipping = list(members = c("shipping")),
    climate  = list(members = c("climate"))
  ),
  vc = list(
    fish = list(members = c("cod", "salmon"))
  ),
  aux = list(
    mask = list(members = c("mask_study_area"))
  )
)
write_yaml(groups, file = paths$groups_yaml)

# config.yml (defaults aligned with implementation plan)
config <- list(
  paths = list(
    cache_dir    = "cache",
    temp_dir     = NULL,
    catalog_dir  = "catalog",
    vrt_cache_dir = "cache/vrt"
  ),
  parallel = list(
    default_cores = "auto",
    memfrac       = 0.6,
    terra_parallel = TRUE
  ),
  alignment = list(
    policy        = "error",
    template_path = NULL,
    warn_on_align = TRUE,
    nodata_default = NULL
  ),
  crs = list(
    enforce_projected = TRUE,
    target_crs        = NULL,
    area_units        = "km2"
  ),
  remote = list(
    allowed_schemes = c("file", "http", "https", "s3", "gs"),
    require_cog     = TRUE,
    requires_auth_default = FALSE
  ),
  cache = list(
    use_cache       = TRUE,
    evict_after_days = 14,
    hash_salt       = NULL
  ),
  logging = list(
    verbosity = "info",
    progress  = TRUE
  ),
  shiny = list(
    preview_aggregate = 2,
    session_cache     = TRUE
  ),
  eo_time = list(
    time_enabled = TRUE
  )
)
write_yaml(config, file = paths$config_yml)

# ------------------------------------------------------------------------------
# AOI examples
poly_coords <- matrix(
  c(2000, 2000,
    18000, 2000,
    18000, 18000,
    2000, 18000,
    2000, 2000),
  ncol = 2,
  byrow = TRUE
)
aoi <- terra::vect(list(poly_coords), type = "polygons", crs = proj_crs)
terra::writeVector(aoi, file.path(paths$aoi_dir, "aoi.geojson"), filetype = "GeoJSON", overwrite = TRUE)

bbox <- list(
  crs  = proj_crs,
  xmin = 5000,
  ymin = 5000,
  xmax = 15000,
  ymax = 15000
)
write_yaml(bbox, file.path(paths$aoi_dir, "bbox.yml"))

message("Stage 1 scaffold created under ", root)
