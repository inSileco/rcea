#' Cube object constructor
#'
#' @param catalog List with layers/groups from `load_catalog()`
#' @param aoi Optional AOI (sf/bbox/WKT/SpatVector)
#' @param sensitivity Optional vulnerability matrix for CEA
#' @return An object of class `rcea_cube`
#' @examples
#' layers_path <- system.file("extdata/catalog/layers.csv", package = "rcea")
#' groups_path <- system.file("extdata/catalog/groups.yaml", package = "rcea")
#' catalog <- load_catalog(layers_path, groups_path)
#' sens <- matrix(
#'   1,
#'   nrow = 2,
#'   ncol = 2,
#'   dimnames = list(
#'     c("cod", "salmon"),
#'     c("shipping", "climate")
#'   )
#' )
#' cube <- make_cube(catalog, sensitivity = sens)
#' cube
#' @export
make_cube <- function(catalog, aoi = NULL, sensitivity = NULL) {
  if (!is.null(aoi)) {
    aoi <- normalize_aoi(aoi)
  }
  if (is.null(sensitivity) && !is.null(catalog$sensitivity)) {
    sensitivity <- catalog$sensitivity
  }
  structure(
    list(
      catalog = catalog,
      aoi = aoi,
      sensitivity = sensitivity,
      stack = NULL
    ),
    class = "rcea_cube"
  )
}

#' Stack layers into SpatRaster objects
#'
#' @param cube rcea_cube
#' @param drivers_id Optional vector of driver layer_ids to stack
#' @param vc_id Optional vector of valued component layer_ids to stack
#' @param catalog_dir Base directory for relative paths in the catalog
#' @param cache Use cache (logical)
#' @param cache_dir Base cache directory
#' @param vrt_cache_dir Optional VRT cache directory (defaults under cache_dir)
#' @param cache_salt Optional salt for cache key
#' @param collect Logical, if TRUE returns in-memory copies of stacks
#' @param align alignment policy, one of "error", "reproject", "template"
#' @param template optional SpatRaster used when `align = "template"`
#' @return rcea_cube with `stack` as a list of SpatRaster objects
#' @examples
#' layers_path <- system.file("extdata/catalog/layers.csv", package = "rcea")
#' groups_path <- system.file("extdata/catalog/groups.yaml", package = "rcea")
#' catalog <- load_catalog(layers_path, groups_path)
#' cube <- make_cube(catalog)
#' cube <- stack_layers(cube)
#' names(cube$stack$drivers)
#' names(cube$stack$vc)
#' @export
stack_layers <- function(cube,
                         drivers_id = NULL,
                         vc_id = NULL,
                         catalog_dir = "catalog",
                         cache = FALSE,
                         cache_dir = "cache",
                         vrt_cache_dir = NULL,
                         cache_salt = NULL,
                         collect = FALSE,
                         align = "error",
                         template = NULL) {
  align <- match.arg(align, c("error", "reproject", "template"))

  lay <- cube$catalog$layers
  if (is.null(drivers_id) && is.null(vc_id)) {
    drivers_id <- lay$layer_id[lay$type %in% "pressure"]
    vc_id <- lay$layer_id[lay$type %in% "vc"]
  }

  lay_dr <- if (is.null(drivers_id)) lay[0, , drop = FALSE] else lay[lay$layer_id %in% drivers_id, , drop = FALSE]
  lay_vc <- if (is.null(vc_id)) lay[0, , drop = FALSE] else lay[lay$layer_id %in% vc_id, , drop = FALSE]

  if (!is.null(drivers_id) && nrow(lay_dr) == 0) stop("No driver layers to stack", call. = FALSE)
  if (!is.null(vc_id) && nrow(lay_vc) == 0) stop("No valued component layers to stack", call. = FALSE)
  if (nrow(lay_dr) == 0 && nrow(lay_vc) == 0) stop("No layers to stack", call. = FALSE)

  build_stack <- function(lay_sub, tag) {
    if (nrow(lay_sub) == 0) {
      return(NULL)
    }

    files <- lay_sub$path
    is_remote <- grepl("^(http|https|s3|gs)://", files)
    missing_local <- !is_remote & !file.exists(files)
    if (any(missing_local)) {
      files[missing_local] <- file.path(catalog_dir %||% ".", files[missing_local])
    }

    if (cache) {
      key <- digest::digest(list(
        tag = tag,
        layer_ids = lay_sub$layer_id,
        aoi = if (is.null(cube$aoi)) NULL else terra::ext(cube$aoi),
        salt = cache_salt
      ))
      if (is.null(vrt_cache_dir)) vrt_cache_dir <- file.path(cache_dir, "vrt")
      vrt_dir <- vrt_cache_dir
      dir.create(vrt_dir, showWarnings = FALSE, recursive = TRUE)
      vrt_path <- file.path(vrt_dir, paste0("stack_", key, ".vrt"))
    } else {
      vrt_path <- NULL
    }

    if (!is.null(vrt_path) && file.exists(vrt_path) && align == "error") {
      r <- terra::vrt(vrt_path)
    } else {
      if (align == "error") {
        r <- terra::rast(files)
      } else {
        target <- if (align == "template") template else terra::rast(files[1])
        if (is.null(target)) {
          stop("Template required when align = 'template'.", call. = FALSE)
        }

        r_list <- lapply(files, function(f) {
          ri <- terra::rast(f)
          if (!isTRUE(terra::compareGeom(ri, target, stopOnError = FALSE))) {
            ri <- terra::project(ri, target)
          }
          ri
        })
        r <- terra::rast(r_list)
      }

      if (!is.null(cube$aoi)) {
        aoi <- cube$aoi
        if (!isTRUE(terra::same.crs(r, aoi))) {
          aoi <- terra::project(aoi, r)
        }
        r <- terra::crop(r, aoi)
        r <- terra::mask(r, aoi)
      }
      if (!is.null(vrt_path) && align == "error") {
        terra::vrt(r, filename = vrt_path, overwrite = TRUE)
        r <- terra::vrt(vrt_path)
      }
    }
    names(r) <- lay_sub$layer_id
    layer_meta <- lay_sub
    layer_meta$layer <- lay_sub$layer_id
    attr(r, "layer_meta") <- layer_meta
    r
  }

  drivers_stack <- build_stack(lay_dr, "drivers")
  vc_stack <- build_stack(lay_vc, "vc")

  if (!is.null(drivers_stack) && !is.null(vc_stack)) {
    aligned <- align_pair(drivers_stack, vc_stack, align = align, template = template)
    drivers_stack <- aligned$drivers
    vc_stack <- aligned$vc
  }

  if (!is.null(drivers_stack)) {
    dr_meta <- lay_dr
    dr_meta$layer <- lay_dr$layer_id
    attr(drivers_stack, "layer_meta") <- dr_meta
  }
  if (!is.null(vc_stack)) {
    vc_meta <- lay_vc
    vc_meta$layer <- lay_vc$layer_id
    attr(vc_stack, "layer_meta") <- vc_meta
  }

  cube$stack <- list(drivers = drivers_stack, vc = vc_stack)
  if (collect) {
    cube$stack <- lapply(cube$stack, function(x) {
      if (inherits(x, "SpatRaster")) return(terra::deepcopy(x))
      x
    })
  }
  cube
}

#' Materialize stack to memory or return SpatRaster
#'
#' @param cube rcea_cube
#' @param collect Logical, if TRUE returns in-memory copy
#' @examples
#' layers_path <- system.file("extdata/catalog/layers.csv", package = "rcea")
#' groups_path <- system.file("extdata/catalog/groups.yaml", package = "rcea")
#' catalog <- load_catalog(layers_path, groups_path)
#' cube <- make_cube(catalog)
#' cube <- stack_layers(cube)
#' stacks <- collect_layers(cube, collect = FALSE)
#' inherits(stacks$drivers, "SpatRaster")
#' @export
collect_layers <- function(cube, collect = FALSE) {
  if (is.null(cube$stack)) stop("Stack not built; call stack_layers() first.", call. = FALSE)
  if (collect) {
    if (inherits(cube$stack, "SpatRaster")) {
      return(terra::deepcopy(cube$stack))
    }
    if (is.list(cube$stack)) {
      return(lapply(cube$stack, function(x) {
        if (inherits(x, "SpatRaster")) {
          return(terra::deepcopy(x))
        }
        x
      }))
    }
    return(cube$stack)
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
