#' Load and validate catalog
#'
#' @param layers Path to layers.csv
#' @param groups Path to groups.yaml
#' @param sensitivity Optional path to sensitivity.csv. If `NULL`, a file named
#'   `sensitivity.csv` in the same directory as `layers` is used when present.
#' @return A list with `layers` tibble, `groups` list, and optional `sensitivity` matrix.
#' @examples
#' layers_path <- system.file("extdata/catalog/layers.csv", package = "rcea")
#' groups_path <- system.file("extdata/catalog/groups.yaml", package = "rcea")
#' catalog <- load_catalog(layers_path, groups_path)
#' str(catalog, max.level = 1)
#' catalog$layers$layer_id
#'
#' @export
load_catalog <- function(layers = file.path("catalog", "layers.csv"),
                         groups = file.path("catalog", "groups.yaml"),
                         sensitivity = NULL) {
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

  sens <- resolve_sensitivity(sensitivity, layers)
  out <- validate_catalog(lay, grp)
  out$sensitivity <- sens
  out
}

resolve_sensitivity <- function(sensitivity, layers_path) {
  base_dir <- dirname(layers_path)
  sens_path <- sensitivity
  if (is.null(sens_path)) {
    candidate <- file.path(base_dir, "sensitivity.csv")
    if (!file.exists(candidate)) return(NULL)
    sens_path <- candidate
  } else {
    is_abs <- grepl("^(/|[A-Za-z]:[\\\\/]|\\\\\\\\)", sens_path)
    if (!is_abs) sens_path <- file.path(base_dir, sens_path)
    if (!file.exists(sens_path)) {
      stop("sensitivity.csv not found at ", sens_path, call. = FALSE)
    }
  }

  read_sensitivity_csv(sens_path)
}

read_sensitivity_csv <- function(path) {
  dat <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(dat) < 2) stop("sensitivity.csv must have one id column and at least one driver column.", call. = FALSE)

  id_col <- if ("vc_id" %in% names(dat)) "vc_id" else names(dat)[1]
  ids <- as.character(dat[[id_col]])
  if (any(!nzchar(ids) | is.na(ids))) stop("sensitivity.csv contains empty VC ids.", call. = FALSE)
  if (any(duplicated(ids))) stop("Duplicate VC ids in sensitivity.csv.", call. = FALSE)

  val_cols <- setdiff(names(dat), id_col)
  if (!length(val_cols)) stop("sensitivity.csv has no driver columns.", call. = FALSE)
  if (any(duplicated(val_cols))) stop("Duplicate driver ids in sensitivity.csv header.", call. = FALSE)

  vals <- dat[, val_cols, drop = FALSE]
  vals[] <- lapply(vals, as.numeric)
  if (any(is.na(as.matrix(vals)))) {
    stop("sensitivity.csv must contain numeric values only.", call. = FALSE)
  }

  m <- as.matrix(vals)
  rownames(m) <- ids
  colnames(m) <- val_cols
  m
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
