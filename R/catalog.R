#' Load and validate catalog
#'
#' @param layers Path to layers.csv
#' @param groups Path to groups.yaml
#' @return A list with `layers` tibble and `groups` list.
#' @examples
#' layers_path <- system.file("extdata/catalog/layers.csv", package = "rcea")
#' groups_path <- system.file("extdata/catalog/groups.yaml", package = "rcea")
#' catalog <- load_catalog(layers_path, groups_path)
#' str(catalog, max.level = 1)
#' catalog$layers$layer_id
#'
#' @export
load_catalog <- function(layers = file.path("catalog", "layers.csv"),
                         groups = file.path("catalog", "groups.yaml")) {
  if (!file.exists(layers)) stop("layers.csv not found at ", layers, call. = FALSE)
  if (!file.exists(groups)) stop("groups.yaml not found at ", groups, call. = FALSE)

  lay <- utils::read.csv(layers, stringsAsFactors = FALSE)
  if ("path" %in% names(lay)) {
    is_remote <- grepl("^(http|https|s3|gs)://", lay$path)
    is_abs <- grepl("^(/|[A-Za-z]:[\\\\/]|\\\\\\\\)", lay$path)
    rel <- !is_remote & !is_abs
    if (any(rel)) {
      base_dir <- dirname(layers)
      lay$path[rel] <- file.path(base_dir, lay$path[rel])
    }
  }
  grp <- yaml::read_yaml(groups)

  validate_catalog(lay, grp)
}

validate_catalog <- function(layers, groups) {
  required <- c("layer_id", "path", "type", "group")
  missing <- setdiff(required, names(layers))
  if (length(missing)) stop("Missing required columns in layers.csv: ", paste(missing, collapse = ", "), call. = FALSE)

  if (!"units" %in% names(layers)) layers$units <- NA_character_

  # type enum
  valid_types <- c("pressure", "vc", "aux")
  if (any(!layers$type %in% valid_types)) {
    bad <- unique(layers$type[!layers$type %in% valid_types])
    stop("Invalid type values: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  # group exists in groups.yaml namespace
  grp_names <- unlist(groups, recursive = TRUE, use.names = TRUE)
  # groups.yaml is nested; collect top-level names
  allowed_groups <- unique(unlist(lapply(groups, names)))
  if (any(!layers$group %in% allowed_groups)) {
    bad <- unique(layers$group[!layers$group %in% allowed_groups])
    stop("Groups not defined in groups.yaml: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  # paths exist (relative or absolute), skip existence check for remote
  missing_paths <- !file.exists(layers$path) & !grepl("^(http|https|s3|gs)://", layers$path)
  if (any(missing_paths)) {
    stop("Layer paths not found: ", paste(layers$path[missing_paths], collapse = ", "), call. = FALSE)
  }

  # unique IDs
  if (any(duplicated(layers$layer_id))) stop("Duplicate layer_id in layers.csv", call. = FALSE)

  meta <- lapply(seq_len(nrow(layers)), function(i) {
    layer_id <- layers$layer_id[[i]]
    path <- layers$path[[i]]

    r <- tryCatch(
      suppressWarnings(terra::rast(path)),
      error = function(e) {
        stop(
          "Unable to read raster for layer '", layer_id, "' at '", path, "': ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    r_crs_desc <- tryCatch(terra::crs(r, describe = TRUE), error = function(e) NULL)
    if (is.data.frame(r_crs_desc) && nrow(r_crs_desc) > 0 &&
        !is.na(r_crs_desc$authority[1]) && !is.na(r_crs_desc$code[1])) {
      crs <- paste0(r_crs_desc$authority[1], ":", r_crs_desc$code[1])
    } else {
      crs <- terra::crs(r)
    }

    r_res <- terra::res(r)
    if (length(r_res) < 2 || any(!is.finite(r_res[1:2]))) {
      stop("Raster resolution is invalid for layer '", layer_id, "' at '", path, "'", call. = FALSE)
    }

    raster_units <- tryCatch(terra::units(r)[1], error = function(e) NA_character_)
    if (!is.character(raster_units) || !nzchar(raster_units)) raster_units <- NA_character_

    list(
      crs = crs,
      res_x = as.numeric(r_res[1]),
      res_y = as.numeric(r_res[2]),
      units = raster_units
    )
  })

  layers$crs <- vapply(meta, `[[`, character(1), "crs")
  layers$res_x <- vapply(meta, `[[`, numeric(1), "res_x")
  layers$res_y <- vapply(meta, `[[`, numeric(1), "res_y")

  provided_units <- as.character(layers$units)
  provided_units[!nzchar(provided_units)] <- NA_character_
  derived_units <- vapply(meta, `[[`, character(1), "units")
  layers$units <- ifelse(is.na(provided_units), derived_units, provided_units)

  if (any(!is.finite(layers$res_x)) || any(!is.finite(layers$res_y))) {
    stop("Derived res_x/res_y must be finite numerics", call. = FALSE)
  }

  list(layers = layers, groups = groups)
}
