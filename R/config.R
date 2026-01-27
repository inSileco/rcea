#' Load package configuration
#'
#' Reads a YAML config file and applies defaults for missing entries.
#'
#' @param path Path to config.yml
#' @return A list with config entries.
load_config <- function(path = "config.yml") {
  if (!file.exists(path)) {
    stop("Config file not found: ", path, call. = FALSE)
  }
  cfg <- yaml::read_yaml(path)

  defaults <- list(
    paths = list(
      cache_dir = "cache",
      temp_dir = NULL,
      catalog_dir = "catalog",
      vrt_cache_dir = "cache/vrt"
    ),
    parallel = list(
      default_cores = "auto",
      memfrac = 0.6,
      terra_parallel = TRUE
    ),
    alignment = list(
      policy = "error",
      template_path = NULL,
      warn_on_align = TRUE,
      nodata_default = NULL
    ),
    crs = list(
      enforce_projected = TRUE,
      target_crs = NULL,
      area_units = "km2"
    ),
    remote = list(
      allowed_schemes = c("file", "http", "https", "s3", "gs"),
      require_cog = TRUE,
      requires_auth_default = FALSE
    ),
    cache = list(
      use_cache = TRUE,
      evict_after_days = 14,
      hash_salt = NULL
    ),
    logging = list(
      verbosity = "info",
      progress = TRUE
    ),
    shiny = list(
      preview_aggregate = 2,
      session_cache = TRUE
    ),
    eo_time = list(
      time_enabled = TRUE
    )
  )

  cfg <- utils::modifyList(defaults, cfg, keep.null = TRUE)
  cfg
}

#' Apply configuration to environment and terra options
#'
#' @param cfg A config list from `load_config()`
apply_config <- function(cfg) {
  if (!is.null(cfg$paths$temp_dir)) {
    terra::terraOptions(tempdir = cfg$paths$temp_dir)
  }
  if (!is.null(cfg$parallel$memfrac)) {
    terra::terraOptions(memfrac = cfg$parallel$memfrac)
  }
  if (!is.null(cfg$parallel$terra_parallel)) {
    terra::terraOptions(parallel = cfg$parallel$terra_parallel)
  }
  invisible(cfg)
}
