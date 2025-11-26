#' Network-scale cumulative effects assessments for a single focal species
#'
#' Assessment of cumulative effects and related metrics using the Beauchesne et al. 2021 method.
#'
#' @eval arguments(c("focus","drivers", "vc", "sensitivity", "metaweb", "trophic_sensitivity", "weights", "output", "output_format"))
#'
#' @examples
#' # Data
#' drivers <- rcea:::drivers
#' vc <- rcea:::vc
#' sensitivity <- rcea:::sensitivity
#' metaweb <- rcea:::metaweb
#' trophic_sensitivity <- rcea::trophic_sensitivity
#'
#' \dontrun{
#' # Network-scale effects for individual species
#' ncea_species(focus = "vc1", drivers, vc, sensitivity, metaweb, trophic_sensitivity)
#' ncea_species(focus = "vc1", drivers, vc, sensitivity, metaweb, trophic_sensitivity, output_format = "COG")
#' }
#' @export
ncea_species <- function(focus, drivers, vc, sensitivity, metaweb, trophic_sensitivity, w_d = 0.5, w_i = 0.25, output = "output/ncea", output_format = "geotiff") {
  # ----------------------------------------------------------------------------
  # NOTE: Function to extract results
  ncea_species_res <- function(dat, metric) {
    # Unlist
    tmp <- data.table::rbindlist(dat)

    if (metric == "species_contribution") {
      out <- here::here(output, metric)
      chk_create(out)
      tmp |>
        dplyr::group_by(vc_id) |>
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(drnm),
            \(x) sum(x, na.rm = TRUE) / area
          )
        ) |>
        dplyr::ungroup() |>
        dplyr::left_join(vc_index, by = "vc_id") |>
        dplyr::rename(vc_from = vc) |>
        dplyr::select(-vc_id) |>
        dplyr::mutate(vc = focus) |>
        vroom::vroom_write(
          file = here::here(out, glue::glue("{focus}.csv")),
          delim = ","
        )
    } else {
      # Effects / km2
      out <- here::here(output, "cekm")
      chk_create(out)
      tmp |>
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(drnm),
            \(x) sum(x, na.rm = TRUE) / area
          )
        ) |>
        dplyr::mutate(vc = focus) |>
        vroom::vroom_write(
          file = here::here(out, glue::glue("{focus}_{metric}.csv")),
          delim = ","
        )

      # Spatial distribution of effects
      dat <- dplyr::left_join(xy, tmp, by = "id_cell") |>
        dplyr::select(-id_cell, -vc_id) |>
        stars::st_as_stars() |>
        merge()
      sf::st_crs(dat) <- prj

      out <- here::here(output, metric)
      chk_create(out)
      if (output_format == "geotiff") {
        stars::write_stars(
          dat,
          dsn = here::here(out, glue::glue("{focus}.tif")),
        )
      } else if (output_format == "COG") {
        as(dat, "Raster") |>
          terra::rast() |>
          terra::writeRaster(
            filename = here::here(out, glue::glue("{focus}.tif")),
            filetype = "COG",
            gdal = c("COMPRESS=LZW", "TILED=YES", "OVERVIEW_RESAMPLING=AVERAGE"),
            overwrite = TRUE
          )
      }
    }
  }
  # ----------------------------------------------------------------------------------------

  # w_d + 2*w_i = 1
  stopifnot(w_d + 2 * w_i == 1)

  # Create output folder
  chk_create(output)

  # Projection
  prj <- sf::st_crs(drivers)

  # Area of focal species
  n <- sum(vc[focus][[1]], na.rm = TRUE)
  area <- vc |>
    stars::st_res() |>
    prod() * n

  # xy coordinates
  xy <- sf::st_coordinates(drivers) |>
    dplyr::mutate(id_cell = 1:dplyr::n())

  # Drivers and vcs as data.frames
  drivers <- as.data.frame(drivers) |>
    dplyr::mutate(id_cell = 1:dplyr::n()) |>
    dplyr::select(-x, -y)
  vc <- as.data.frame(vc) |>
    dplyr::select(-x, -y) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.logical))

  # Index of vc
  vc_index <- data.frame(
    vc = colnames(vc),
    vc_id = 1:ncol(vc)
  )

  # Add cell_id to vc
  vc <- dplyr::mutate(vc, id_cell = 1:dplyr::n())

  # Only select cells in which focal species is found
  vc <- dplyr::filter(vc, get(focus))
  drivers <- dplyr::filter(drivers, id_cell %in% vc$id_cell)

  # Create object with cell id and remove cell_id from objects
  id_cell <- vc$id_cell
  vc <- dplyr::select(vc, -id_cell) |> as.matrix()
  drivers <- dplyr::select(drivers, -id_cell) |> as.matrix()

  # Names
  vcnm <- colnames(vc)
  drnm <- colnames(drivers)

  # ID of focal species
  focusID <- which(colnames(vc) == focus)

  # Lists to store results
  ncea_res <- list()
  ncea_res$species_contribution <- list()
  ncea_res$direct <- list()
  ncea_res$indirect <- list()

  # 3-species motifs for full metaweb
  motifs <- triads(metaweb, trophic_sensitivity)

  for (i in 1:nrow(vc)) {
    v <- as.matrix(vc[i, ])
    d <- drivers[i, ]

    # id of species in cell i
    idsp <- which(v)

    # Direct effects
    direct_effects <- sweep(sensitivity[idsp, ], MARGIN = 2, d, `*`)

    # Direct pathways
    direct_pathways <- data.frame(
      vc_id = idsp,
      effect = as.logical(rowSums(direct_effects, na.rm = TRUE))
    )

    # Add species id to direct_effects
    direct_effects <- as.data.frame(direct_effects) |>
      dplyr::mutate(vc_id = idsp)

    # Select triads involving species found in cell j & focal species
    uid <- (motifs$vc_i %in% idsp & motifs$vc_j %in% idsp & motifs$vc_k %in% idsp) &
      (motifs$vc_i == focusID | motifs$vc_j == focusID | motifs$vc_k == focusID)
    dat <- motifs[uid, ]

    if (nrow(dat) > 0) {
      dat <- dplyr::left_join(dat, direct_pathways, by = c("vc_i" = "vc_id")) |>
        dplyr::left_join(direct_pathways, by = c("vc_j" = "vc_id")) |>
        dplyr::left_join(direct_pathways, by = c("vc_k" = "vc_id")) |>
        dplyr::filter(effect_i == effect.x & effect_j == effect.y & effect_k == effect)

      tmp <- dat[, c("vc_i", "vc_j", "vc_k")]
      sens <- dat[, c("TS_i", "TS_j", "TS_k")]
      tmp2 <- tmp[rep(1:nrow(dat), each = 3), ]
      dat <- data.frame(
        vc_id = rep(c(t(tmp)), each = 3),
        interaction = c(t(tmp2)),
        Sensitivity = rep(c(t(sens)), each = 3)
      ) |>
        dplyr::filter(vc_id == focusID)
    } else {
      dat <- data.frame(
        vc_id = focusID,
        interaction = focusID,
        Sensitivity = 1
      )
    }

    # Add direct effects
    dat <- dplyr::left_join(dat, direct_effects, by = c("interaction" = "vc_id")) |>
      dplyr::mutate(
        direct = vc_id == interaction,
        M = sum(direct)
      )

    # Weights
    dat$weight <- dat$direct
    dat$weight[dat$weight] <- w_d
    dat$weight[!dat$weight] <- w_i

    # Effects, i.e. multiply driver columns by trophic sensitivity and weights
    dat <- dat |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(drnm),
          \(x) (x * weight * Sensitivity) / M
        )
      ) |>
      dplyr::select(
        vc_id, interaction, direct, dplyr::all_of(drnm),
        -weight, -Sensitivity, -M
      )

    # Species contribution to indirect effects
    ncea_res$species_contribution[[i]] <- get_species_contribution_sp(dat) |>
      dplyr::rename(vc_id = interaction) |>
      dplyr::mutate(id_cell = id_cell[i])

    # Direct & indirect effects
    direct_indirect <- get_direct_indirect_sp(dat)

    # Direct
    ncea_res$direct[[i]] <- dplyr::filter(direct_indirect, direct) |>
      dplyr::select(-direct) |>
      dplyr::mutate(id_cell = id_cell[i])

    # Indirect
    ncea_res$indirect[[i]] <- dplyr::filter(direct_indirect, !direct) |>
      dplyr::select(-direct) |>
      dplyr::mutate(id_cell = id_cell[i])

    # Net effects
    ncea_res$net[[i]] <- get_net_sp(dat) |>
      dplyr::mutate(id_cell = id_cell[i])
  }




  # Direct, indirect and net effects
  ncea_species_res(ncea_res$direct, "direct")
  ncea_species_res(ncea_res$indirect, "indirect")
  ncea_species_res(ncea_res$net, "net")
  ncea_res$direct <- ncea_res$indirect <- ncea_res$net <- NULL

  # Species contribution summarized by species and area
  out <- here::here(output, "species_contribution")
  chk_create(out)
  ncea_species_res(ncea_res$species_contribution, "species_contribution")
}


# ==============================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------------------------------------------------------
get_direct_indirect_sp <- function(motif_effects) {
  notDr <- c("id_cell", "vc_id", "interaction", "Sensitivity", "direct", "M", "weight")
  drNames <- colnames(motif_effects)
  drNames <- drNames[!drNames %in% notDr]
  dplyr::group_by(motif_effects, vc_id, direct) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(drNames),
        \(x) sum(x, na.rm = TRUE)
      )
    ) |>
    dplyr::ungroup()
}


# ==============================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------------------------------------------------------
get_net_sp <- function(motif_effects) {
  notDr <- c("id_cell", "vc_id", "interaction", "Sensitivity", "direct", "M", "weight")
  drNames <- colnames(motif_effects)
  drNames <- drNames[!drNames %in% notDr]
  dplyr::group_by(motif_effects, vc_id) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(drNames),
        \(x) sum(x, na.rm = TRUE)
      )
    ) |>
    dplyr::ungroup()
}

# ==============================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------------------------------------------------------
get_species_contribution_sp <- function(motif_effects) {
  notDr <- c("id_cell", "vc_id", "interaction", "Sensitivity", "direct", "M", "weight")
  drNames <- colnames(motif_effects)
  drNames <- drNames[!drNames %in% notDr]
  dplyr::filter(motif_effects, !direct) |>
    dplyr::group_by(interaction) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(drNames),
        \(x) sum(x, na.rm = TRUE)
      )
    ) |>
    dplyr::ungroup()
}
