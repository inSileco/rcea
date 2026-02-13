#' Preprocess stacked layers before CEA
#'
#' Apply optional per-layer transforms and normalization to cube stacks. Global
#' defaults can be overridden by catalog columns (`transform`, `normalize`) when
#' available.
#'
#' @param cube rcea_cube with `stack_layers()` already run.
#' @param use_catalog_overrides logical; if TRUE, use `transform` and
#'   `normalize` values from layer metadata when present.
#' @param drivers_transform Default transform for driver layers. One of
#'   `"none"`, `"log1p"`, `"log10"`, `"sqrt"`.
#' @param drivers_normalize Default normalization for driver layers. One of
#'   `"none"`, `"minmax"`, `"max"`.
#' @param vc_transform Default transform for valued component layers.
#' @param vc_normalize Default normalization for valued component layers.
#'
#' @return rcea_cube with transformed stacks.
#'
#' @export
preprocess_layers <- function(cube,
                              use_catalog_overrides = TRUE,
                              drivers_transform = "none",
                              drivers_normalize = "none",
                              vc_transform = "none",
                              vc_normalize = "none") {
  if (!inherits(cube, "rcea_cube")) stop("cube must be an rcea_cube.", call. = FALSE)
  if (is.null(cube$stack)) stop("cube$stack is empty; call stack_layers() first.", call. = FALSE)
  if (!is.list(cube$stack)) stop("cube$stack must be a list with drivers and vc stacks.", call. = FALSE)

  valid_transform <- c("none", "log1p", "log10", "sqrt")
  valid_normalize <- c("none", "minmax", "max")
  if (!drivers_transform %in% valid_transform) stop("Invalid drivers_transform.", call. = FALSE)
  if (!vc_transform %in% valid_transform) stop("Invalid vc_transform.", call. = FALSE)
  if (!drivers_normalize %in% valid_normalize) stop("Invalid drivers_normalize.", call. = FALSE)
  if (!vc_normalize %in% valid_normalize) stop("Invalid vc_normalize.", call. = FALSE)

  apply_transform <- function(r, method) {
    switch(method,
      "none" = r,
      "log1p" = terra::app(r, function(x) {
        x[x < -1] <- NA_real_
        log1p(x)
      }),
      "log10" = terra::app(r, function(x) {
        x[x <= 0] <- NA_real_
        log10(x)
      }),
      "sqrt" = terra::app(r, function(x) {
        x[x < 0] <- NA_real_
        sqrt(x)
      })
    )
  }

  apply_normalize <- function(r, method) {
    if (method == "none") return(r)

    if (method == "minmax") {
      stats <- terra::global(r, c("min", "max"), na.rm = TRUE)
      mn <- as.numeric(stats[1, "min"])
      mx <- as.numeric(stats[1, "max"])
      if (!is.finite(mn) || !is.finite(mx)) return(r)
      if (mx == mn) return(r * 0)
      return((r - mn) / (mx - mn))
    }

    if (method == "max") {
      mx <- as.numeric(terra::global(r, "max", na.rm = TRUE)[1, 1])
      if (!is.finite(mx) || mx == 0) return(r * 0)
      return(r / mx)
    }

    r
  }

  resolve_methods <- function(stack, default_transform, default_normalize) {
    n <- terra::nlyr(stack)
    transform <- rep(default_transform, n)
    normalize <- rep(default_normalize, n)
    names(transform) <- names(stack)
    names(normalize) <- names(stack)

    if (isTRUE(use_catalog_overrides)) {
      meta <- attr(stack, "layer_meta")
      if (is.data.frame(meta) && nrow(meta)) {
        meta_layer <- if ("layer" %in% names(meta)) meta$layer else meta$layer_id
        idx <- match(names(stack), meta_layer)

        if ("transform" %in% names(meta)) {
          override <- as.character(meta$transform[idx])
          override[is.na(override) | !nzchar(override)] <- NA_character_
          bad <- unique(override[!is.na(override) & !override %in% valid_transform])
          if (length(bad)) stop("Invalid transform values in metadata: ", paste(bad, collapse = ", "), call. = FALSE)
          transform[!is.na(override)] <- override[!is.na(override)]
        }
        if ("normalize" %in% names(meta)) {
          override <- as.character(meta$normalize[idx])
          override[is.na(override) | !nzchar(override)] <- NA_character_
          bad <- unique(override[!is.na(override) & !override %in% valid_normalize])
          if (length(bad)) stop("Invalid normalize values in metadata: ", paste(bad, collapse = ", "), call. = FALSE)
          normalize[!is.na(override)] <- override[!is.na(override)]
        }
      }
    }

    list(transform = transform, normalize = normalize)
  }

  preprocess_stack <- function(stack, default_transform, default_normalize) {
    if (is.null(stack)) return(NULL)
    if (!inherits(stack, "SpatRaster")) stop("stack entries must be SpatRaster.", call. = FALSE)

    methods <- resolve_methods(stack, default_transform, default_normalize)
    layers <- lapply(seq_len(terra::nlyr(stack)), function(i) {
      ri <- stack[[i]]
      ri <- apply_transform(ri, methods$transform[[i]])
      ri <- apply_normalize(ri, methods$normalize[[i]])
      names(ri) <- names(stack)[[i]]
      ri
    })
    out <- terra::rast(layers)
    names(out) <- names(stack)

    meta <- attr(stack, "layer_meta")
    if (is.data.frame(meta) && nrow(meta)) {
      meta_layer <- if ("layer" %in% names(meta)) meta$layer else meta$layer_id
      idx <- match(names(out), meta_layer)
      meta$applied_transform <- methods$transform[idx]
      meta$applied_normalize <- methods$normalize[idx]
      attr(out, "layer_meta") <- meta
    }

    out
  }

  cube$stack$drivers <- preprocess_stack(cube$stack$drivers, drivers_transform, drivers_normalize)
  cube$stack$vc <- preprocess_stack(cube$stack$vc, vc_transform, vc_normalize)
  cube
}
