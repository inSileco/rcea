#' Cube object constructor
#'
#' @param catalog List with layers/groups from `load_catalog()`
#' @param config List from `load_config()`
#' @param aoi Optional AOI (sf/bbox/WKT/SpatVector)
#' @return An object of class `rcea_cube`
make_cube <- function(catalog, config, aoi = NULL) {
  if (!is.null(aoi)) {
    aoi <- normalize_aoi(aoi)
  }
  structure(
    list(
      catalog = catalog,
      config = config,
      aoi = aoi,
      stack = NULL
    ),
    class = "rcea_cube"
  )
}

#' Stack layers into a SpatRaster
#'
#' @param cube rcea_cube
#' @param layer_ids Optional vector of layer_ids to stack (default all)
#' @param cache Use cache (logical)
#' @return SpatRaster
stack_layers <- function(cube, layer_ids = NULL, cache = cube$config$cache$use_cache) {
  lay <- cube$catalog$layers
  if (!is.null(layer_ids)) {
    lay <- lay[lay$layer_id %in% layer_ids, , drop = FALSE]
  }
  if (nrow(lay) == 0) stop("No layers to stack", call. = FALSE)

  # Build absolute paths for local files
  files <- lay$path
  is_remote <- grepl("^(http|https|s3|gs)://", files)
  files[!is_remote] <- file.path(cube$config$paths$catalog_dir %||% ".", files[!is_remote])

  # Optional VRT cache (deterministic key)
  if (cache) {
    key <- digest::digest(list(
      layer_ids = lay$layer_id,
      aoi = if (is.null(cube$aoi)) NULL else terra::ext(cube$aoi),
      template = cube$config$alignment$template_path,
      align = cube$config$alignment$policy,
      salt = cube$config$cache$hash_salt
    ))
    vrt_dir <- cube$config$paths$vrt_cache_dir %||% file.path(cube$config$paths$cache_dir, "vrt")
    dir.create(vrt_dir, showWarnings = FALSE, recursive = TRUE)
    vrt_path <- file.path(vrt_dir, paste0("stack_", key, ".vrt"))
  } else {
    vrt_path <- NULL
  }

  if (!is.null(vrt_path) && file.exists(vrt_path)) {
    r <- terra::vrt(vrt_path)
  } else {
    r <- terra::rast(files)
    # Apply AOI if present
    if (!is.null(cube$aoi)) {
      r <- terra::crop(r, cube$aoi)
      r <- terra::mask(r, cube$aoi)
    }
    if (!is.null(vrt_path)) {
      terra::vrt(r, filename = vrt_path, overwrite = TRUE)
      r <- terra::vrt(vrt_path)
    }
  }
  names(r) <- lay$layer_id

  cube$stack <- r
  cube
}

#' Materialize stack to memory or return SpatRaster
#'
#' @param cube rcea_cube
#' @param collect Logical, if TRUE returns in-memory copy
collect_layers <- function(cube, collect = FALSE) {
  if (is.null(cube$stack)) stop("Stack not built; call stack_layers() first.", call. = FALSE)
  if (collect) {
    return(terra::deepcopy(cube$stack))
  }
  cube$stack
}

#' @export
print.rcea_cube <- function(x, ...) {
  cat("<rcea_cube>\n")
  cat("  layers:", nrow(x$catalog$layers), "\n")
  cat("  groups:", paste(names(x$catalog$groups), collapse = ", "), "\n")
  cat("  aoi:", if (is.null(x$aoi)) "none" else "set", "\n")
  invisible(x)
}
