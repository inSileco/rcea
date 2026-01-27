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
  required <- c("layer_id", "path", "type", "group", "units", "crs", "res_x", "res_y")
  missing <- setdiff(required, names(layers))
  if (length(missing)) stop("Missing required columns in layers.csv: ", paste(missing, collapse = ", "), call. = FALSE)

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

  # numeric res
  if (any(!is.finite(layers$res_x)) || any(!is.finite(layers$res_y))) {
    stop("res_x/res_y must be finite numerics", call. = FALSE)
  }

  list(layers = layers, groups = groups)
}
