#' Network-scale cumulative effects assessments
#'
#' Assessment of cumulative effects and related metrics using the Beauchesne et al. 2021 method.
#'
#' @eval arguments(c("drivers", "vc", "sensitivity", "metaweb", "trophic_sensitivity", "weights", "motif_effects"))
#' @param exportAs string, the type of object that should be created, either multiple "data.frame" or "stars" objects
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
#' # Network-scale effects
#' beauchesne <- ncea(drivers, vc, sensitivity, metaweb, trophic_sensitivity)
#' plot(beauchesne$net)
#' plot(beauchesne$direct)
#' plot(beauchesne$indirect)
#' }
#' @export
ncea <- function(drivers, vc, sensitivity, metaweb, trophic_sensitivity, w_d = 0.5, w_i = 0.25, exportAs = "stars") {
  # 3-species motifs for full metaweb
  motifs <- triads(metaweb, trophic_sensitivity)

  # Direct effects, i.e. Halpern approach
  direct_effect <- cea(drivers, vc, sensitivity) |>
    make_array()

  # Pathways of direct effect
  direct_pathways <- cea_pathways(direct_effect, vc)

  if (nrow(direct_pathways) > 0) {
    # Pathways of indirect effect
    indirect_pathways <- ncea_pathways_(direct_pathways, motifs)

    # Effects for motifs in each cell
    motif_summary <- ncea_motifs(direct_effect, indirect_pathways)

    # Measure effects on each motif
    motif_effects <- ncea_effects(motif_summary, w_d, w_i)

    # params for stars object creation
    xy <- sf::st_coordinates(drivers)
    drNames <- names(drivers)

    # Species contribution to indirect effects
    species_contribution <- get_species_contribution(motif_effects) |>
      dplyr::rename(vc_id = interaction)

    # Direct & indirect effects
    direct_indirect <- get_direct_indirect(motif_effects)
    direct <- dplyr::filter(direct_indirect, direct) |>
      dplyr::select(-direct)
    indirect <- dplyr::filter(direct_indirect, !direct) |>
      dplyr::select(-direct)

    # Net effects
    net <- get_net(motif_effects)

    # Effects / km2
    cekm <- get_cekm_ncea(motif_effects, vc)

    if (exportAs == "stars") {
      species_contribution <- make_stars(species_contribution, drivers, vc)
      direct <- make_stars(direct, drivers, vc)
      indirect <- make_stars(indirect, drivers, vc)
      net <- make_stars(net, drivers, vc)
    }

    # Return
    list(
      # motif_effects = motif_effects,
      xy = xy,
      net = net,
      direct = direct,
      indirect = indirect,
      species_contribution = species_contribution,
      cekm = cekm
    )
  } else {
    NULL
  }
}


#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn ncea transform effects assessment into binary 2D matrix to assess the presence of an effect to a valued component in a specific grid cell
#' @export
#' @param effect TODO
cea_binary <- function(effect) {
  (effect / effect) |>
    apply(c(1, 2), sum, na.rm = TRUE) |>
    as.data.frame() |>
    dplyr::mutate(id_cell = 1:dplyr::n()) |>
    tidyr::pivot_longer(cols = -c(id_cell), names_to = "vc", values_to = "effect") |>
    dplyr::mutate(effect = as.logical(effect))
}


#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn ncea assess all triads of interest from metaweb and attach trophic sensitivities
#' @export
triads <- function(metaweb, trophic_sensitivity) {
  motifs <- as.matrix(metaweb) |>
    # Motif census of metaweb of interest
    motifcensus::motif_census_triplet() |>
    # Reduce number of of columns
    dplyr::select(psum, i, j, k, pid_i, pid_j, pid_k) |>
    # Select motifs of interest
    # 5:exploitative competition; 12:linear chain; 28:apparent competition; 36:omnivory
    dplyr::filter(psum %in% c(5, 12, 28, 36)) |>
    # Give species proper id, motifcensus puts them back to 0
    dplyr::mutate(i = i + 1, j = j + 1, k = k + 1) |>
    # Rename columns and pivot wider
    ## The next steps are all to reposition species in proper order of i,j,k
    ## These will then be used to identify all possible pathways of effect
    dplyr::rename(i_vc = i, j_vc = j, k_vc = k, i_pos = pid_i, j_pos = pid_j, k_pos = pid_k) |>
    dplyr::mutate(uid = 1:dplyr::n()) |>
    tidyr::pivot_longer(
      cols = c("i_vc", "j_vc", "k_vc", "i_pos", "j_pos", "k_pos"),
      names_to = c("sp", ".value"),
      names_sep = "_"
    ) |>
    dplyr::group_by(uid, psum) |>
    dplyr::arrange(uid, psum, pos) |> # Reposition species in proper order
    dplyr::ungroup() |>
    dplyr::mutate(sp = rep(c("i", "j", "k"), dplyr::n() / 3)) |>
    tidyr::pivot_wider(
      id_cols = c(uid, psum),
      names_from = sp,
      values_from = c(vc, pos)
    ) |>
    dplyr::select(-uid, -pos_i, -pos_j, -pos_k)

  # Now we can join with trophic sensitivity data
  ## Start by pivoting sensitivities wider to obtain a single line per pathway of effect
  sensitivity <- tidyr::pivot_wider(
    trophic_sensitivity,
    id_cols = c("motifID", "pathID", "pi", "pj", "pk"),
    names_from = Species,
    values_from = sensitivity_1
  ) |>
    dplyr::rename(
      effect_i = pi, effect_j = pj, effect_k = pk,
      TS_i = i, TS_j = j, TS_k = k
    )

  ## Join with censusmotif data
  dat <- dplyr::left_join(
    motifs,
    sensitivity,
    by = c("psum" = "motifID"),
    relationship = "many-to-many"
  )

  # Return
  dat
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn ncea pathways of direct effect
#' @export
cea_pathways <- function(effect, vc) {
  # Binary effects (for pathways)
  bin <- cea_binary(effect)

  # vc as data.frame
  vc_df <- as.data.frame(vc) |>
    dplyr::select(-x, -y)

  # Index of vc
  vc_index <- data.frame(
    vc = colnames(vc_df),
    vc_id = 1:ncol(vc_df)
  )

  # Species in each cell and presence of direct effect
  vc_df <- vc_df |>
    dplyr::mutate(id_cell = 1:dplyr::n()) |>
    tidyr::pivot_longer(cols = -c(id_cell), names_to = "vc", values_to = "presence") |>
    tidyr::drop_na() |>
    dplyr::left_join(vc_index, by = c("vc")) |>
    dplyr::left_join(bin, by = c("id_cell", "vc")) |>
    dplyr::select(-vc, -presence)
}

# ==============================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------------------------------------------------------
#' @describeIn ncea assess all triads of interest from metaweb
#' @export
#' @param vc_id TODO
# Difference starts here
ncea_pathways <- function(vc_id, motifs) {
  uid <- motifs$vc_i %in% vc_id$vc_id &
    motifs$vc_j %in% vc_id$vc_id &
    motifs$vc_k %in% vc_id$vc_id
  dat <- motifs[uid, ]
  id_cell <- vc_id$id_cell[1]
  vc_id <- vc_id[, c("vc_id", "effect")]

  if (nrow(dat) > 0) {
    dat <- dplyr::left_join(dat, vc_id, by = c("vc_i" = "vc_id")) |>
      dplyr::left_join(vc_id, by = c("vc_j" = "vc_id")) |>
      dplyr::left_join(vc_id, by = c("vc_k" = "vc_id")) |>
      dplyr::filter(effect_i == effect.x & effect_j == effect.y & effect_k == effect)

    tmp <- dat[, c("vc_i", "vc_j", "vc_k")]
    sens <- dat[, c("TS_i", "TS_j", "TS_k")]
    tmp2 <- tmp[rep(1:nrow(dat), each = 3), ]
    dat <- data.frame(
      vc_id = rep(c(t(tmp)), each = 3),
      interaction = c(t(tmp2)),
      Sensitivity = rep(c(t(sens)), each = 3)
    )

    # Add missing species
    uid <- !vc_id$vc_id %in% dat$vc_id
    if (any(uid)) {
      add <- data.frame(
        vc_id = vc_id$vc_id[uid],
        interaction = vc_id$vc_id[uid],
        Sensitivity = 1
      )
      dat <- rbind(dat, add)
    }
  } else {
    dat <- data.frame(
      vc_id = vc_id$vc_id,
      interaction = vc_id$vc_id,
      Sensitivity = 1
    )
  }
  dat$id_cell <- id_cell

  # Return
  dat
}

# ----------------------------------------------------------------------------------------
#' @describeIn ncea apply ncea_pathways to get pathways of indirect effect and trophic sensitivity for all cells
#' @export
#' @param direct_pathways TODO
#' @param motifs TODO
ncea_pathways_ <- function(direct_pathways, motifs) {
  dat <- dplyr::group_by(direct_pathways, id_cell) |>
    dplyr::group_split() |>
    lapply(function(x) ncea_pathways(x, motifs))

  dat
}

# ========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------------------------------------------------------------------------
#' @describeIn ncea get effects of drivers for species in all motifs in each cell
#' @param direct_effect TODO
#' @param indirect_pathways TODO
#' @export
ncea_motifs <- function(direct_effect, indirect_pathways) {
  drNames <- dimnames(direct_effect)[3][[1]]
  lapply(
    indirect_pathways,
    function(x) {
      cbind(x, direct_effect[x$id_cell[1], x$interaction, ])
    }
  ) |>
    dplyr::bind_rows() |>
    dplyr::select(
      id_cell,
      vc_id,
      interaction,
      Sensitivity,
      dplyr::all_of(drNames)
    ) |>
    dplyr::mutate(direct = vc_id == interaction) |>
    dplyr::group_by(id_cell, vc_id) |>
    dplyr::mutate(M = sum(direct)) |>
    dplyr::ungroup()
}


#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn ncea evaluate effects for all motifs using trophic sensitivity and effect weights
#' @param motif_summary TODO
#' @export
ncea_effects <- function(motif_summary, w_d = 0.5, w_i = 0.25) {
  # w_d + 2*w_i = 1
  stopifnot(w_d + 2 * w_i == 1)
  notDr <- c("id_cell", "vc_id", "interaction", "Sensitivity", "direct", "M", "weight")
  drNames <- colnames(motif_summary)
  drNames <- drNames[!drNames %in% notDr]

  # Direct & indirect weights
  w <- data.frame(
    direct = c(TRUE, FALSE),
    weight = c(w_d, w_i)
  )
  motif_summary <- dplyr::left_join(motif_summary, w, by = "direct")

  # Effects, i.e. multiply driver columns by trophic sensitivity and weights
  motif_effects <- motif_summary |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(drNames),
        \(x) (x * weight * Sensitivity) / M
      )
    ) |>
    dplyr::select(
      id_cell, M, vc_id, interaction, direct, dplyr::all_of(drNames),
      -weight, -Sensitivity, -M
    )

  # Return
  motif_effects
}


#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn ncea get contribution of species to indirect effects
#' @export
get_species_contribution <- function(motif_effects) {
  notDr <- c("id_cell", "vc_id", "interaction", "Sensitivity", "direct", "M", "weight")
  drNames <- colnames(motif_effects)
  drNames <- drNames[!drNames %in% notDr]
  dplyr::filter(motif_effects, !direct) |>
    dplyr::group_by(id_cell, interaction) |>
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
#' @describeIn ncea get direct and indirect effects of drivers
#' @export
get_direct_indirect <- function(motif_effects) {
  notDr <- c("id_cell", "vc_id", "interaction", "Sensitivity", "direct", "M", "weight")
  drNames <- colnames(motif_effects)
  drNames <- drNames[!drNames %in% notDr]
  dplyr::group_by(motif_effects, id_cell, vc_id, direct) |>
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
#' @describeIn ncea get net effects of drivers
#' @export
get_net <- function(motif_effects) {
  notDr <- c("id_cell", "vc_id", "interaction", "Sensitivity", "direct", "M", "weight")
  drNames <- colnames(motif_effects)
  drNames <- drNames[!drNames %in% notDr]
  dplyr::group_by(motif_effects, id_cell, vc_id) |>
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
#' @describeIn ncea get effects per km2
#' @export
get_cekm_ncea <- function(motif_effects, vc) {
  # Driver names
  notDr <- c("id_cell", "vc_id", "interaction", "Sensitivity", "direct", "M", "weight")
  drNames <- colnames(motif_effects)
  drNames <- drNames[!drNames %in% notDr]

  # vc as data.frame
  vc_df <- as.data.frame(vc) |>
    dplyr::select(-x, -y)

  # Index of vc
  vc_index <- data.frame(
    vc = colnames(vc_df),
    vc_id = 1:ncol(vc_df)
  )

  # Calculate area, i.e. number of cells (assuming 1km2 grid cells)
  vc_df <- vc_df |>
    dplyr::mutate(id_cell = 1:dplyr::n()) |>
    tidyr::pivot_longer(cols = -c(id_cell), names_to = "vc", values_to = "presence") |>
    dplyr::group_by(vc) |>
    dplyr::summarise(km2 = sum(presence, na.rm = TRUE)) |>
    dplyr::left_join(vc_index, by = "vc") |>
    dplyr::select(-vc) |>
    dplyr::ungroup()

  # Direct & indirect effects
  direct_indirect <- get_direct_indirect(motif_effects) |>
    dplyr::left_join(vc_df, by = "vc_id") |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(drNames),
        \(x) x / km2
      )
    ) |>
    dplyr::group_by(vc_id, direct) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(drNames),
        \(x) sum(x, na.rm = TRUE)
      )
    ) |>
    dplyr::mutate(
      type =
        dplyr::case_when(
          direct ~ "direct",
          !direct ~ "indirect"
        )
    ) |>
    dplyr::select(-direct) |>
    dplyr::ungroup()

  # Net effects
  net <- direct_indirect |>
    dplyr::select(-type) |>
    dplyr::group_by(vc_id) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(drNames),
        \(x) sum(x, na.rm = TRUE)
      )
    ) |>
    dplyr::mutate(type = "net") |>
    dplyr::ungroup()


  # Summary table
  cekm <- dplyr::bind_rows(direct_indirect, net) |>
    dplyr::arrange(vc_id) |>
    dplyr::select(vc_id, type, dplyr::all_of(drNames))

  # Return
  cekm
}




# ==============================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------------------------------------------------------
#' @describeIn ncea split the assessment in smaller parts for larger analyses that run into memory issues and need to be run in parallel
#' @export
#' @param output TODO
#' @param niter TODO
#' @param run TODO
ncea_split <- function(drivers, vc, sensitivity, metaweb, trophic_sensitivity, w_d = 0.5, w_i = 0.25, output = "output/ncea", niter = NULL, run = NULL) {
  # Output folder
  chk_create(output)

  # 3-species motifs for full metaweb
  out <- here::here(output, "motifs.csv")
  if (!file.exists(out)) {
    motifs <- triads(metaweb, trophic_sensitivity)
    vroom::vroom_write(motifs, out, delim = ",")
  } else {
    motifs <- vroom::vroom(out)
  }

  # Direct effects, i.e. Halpern approach
  out <- here::here(output, "direct_effect.RData")
  if (!file.exists(out)) {
    direct_effect <- cea(drivers, vc, sensitivity) |>
      make_array()
    save(direct_effect, file = out)
  } else {
    load(out)
  }

  # Pathways of direct effect
  out <- here::here(output, "direct_pathways.csv")
  if (!file.exists(out)) {
    direct_pathways <- cea_pathways(direct_effect, vc)
    vroom::vroom_write(direct_pathways, out, delim = ",")
  } else {
    direct_pathways <- vroom::vroom(out)
  }

  if (nrow(direct_pathways) > 0) {
    out <- list()
    out$motif <- here::here(output, "motif_summary")
    out$species_contribution <- here::here(output, "partial", "species_contribution")
    out$direct <- here::here(output, "partial", "direct")
    out$indirect <- here::here(output, "partial", "indirect")
    out$net <- here::here(output, "partial", "net")
    lapply(out, chk_create)

    # Motif summaries
    files <- dir(out$motif, full.names = TRUE)
    if (length(files) != niter) {
      # Pathways of indirect effect and effects for motifs in each cell
      temp <- dplyr::group_by(direct_pathways, id_cell) |>
        dplyr::group_split()
      rm(direct_pathways) # save memory

      # Iterate over cells
      if (is.null(niter)) stop("You must provide the number of iterations into which you wish to split your assessment, i.e. the number of cells per iteration that will be computed.")

      # Iteration data.frame
      iter <- seq(1, length(temp), length.out = niter + 1) |> ceiling()
      iter <- data.frame(
        from = iter[1:(niter)],
        to = c(iter[2:(niter)] - 1, length(temp))
      )
    }

    # Iterate and export
    files <- c(files, "a") # To avoid error in if statement (not pretty, but efficient and it works)
    for (i in run) {
      if (!any(stringr::str_detect(files, sprintf("%04d", i)))) {
        # Range
        beg <- iter$from[i]
        end <- iter$to[i]

        # Pathways of indirect effect
        dat <- lapply(
          temp[beg:end],
          function(x) ncea_pathways(x, motifs)
        )

        # Effects for motifs in each cell
        dat <- ncea_motifs(direct_effect, dat)

        # Export (to test multiple w_d or w_i)
        vroom::vroom_write(
          dat,
          sprintf(paste0(out$motif, "/motif_summary.%04d.csv"), i),
          delim = ","
        )
      }
    }
    rm(temp, direct_effect) # save memory

    for (i in run) {
      # Load motif_summary
      dat <- vroom::vroom(sprintf(paste0(out$motif, "/motif_summary.%04d.csv"), i))

      # Measure effects on each motif
      dat <- ncea_effects(dat, w_d, w_i)

      # Species contribution to indirect effects
      get_species_contribution(dat) |>
        dplyr::rename(vc_id = interaction) |>
        vroom::vroom_write(
          sprintf(paste0(out$species_contribution, "/species_contribution.%04d.csv"), i),
          delim = ","
        )

      # Direct & indirect effects
      direct_indirect <- get_direct_indirect(dat)

      ## Direct effects
      dplyr::filter(direct_indirect, direct) |>
        dplyr::select(-direct) |>
        vroom::vroom_write(
          sprintf(paste0(out$direct, "/direct.%04d.csv"), i),
          delim = ","
        )

      ## Indirect effects
      dplyr::filter(direct_indirect, !direct) |>
        dplyr::select(-direct) |>
        vroom::vroom_write(
          sprintf(paste0(out$indirect, "/indirect.%04d.csv"), i),
          delim = ","
        )
      rm(direct_indirect) # save memory

      # Net effects
      get_net(dat) |>
        vroom::vroom_write(
          sprintf(paste0(out$net, "/net.%04d.csv"), i),
          delim = ","
        )
    }

    # Combine and export rasters if all iterations are present
    # NOTE: In theory, if parallelized, the last run to complete should complete the assessment
    files <- dir(out$net, full.names = TRUE)
    if (length(files) == niter) {
      # Species contribution to indirect effects
      dir(out$species_contribution, full.names = TRUE) |>
        purrr::map(vroom::vroom) |>
        purrr::list_rbind() |>
        make_stars(drivers, vc) |>
        export_stars(output, "species_contribution", length(vc))

      ## Direct effects
      dir(out$direct, full.names = TRUE) |>
        purrr::map(vroom::vroom) |>
        purrr::list_rbind() |>
        make_stars(drivers, vc) |>
        export_stars(output, "direct", length(vc))

      ## Indirect effects
      dir(out$indirect, full.names = TRUE) |>
        purrr::map(vroom::vroom) |>
        purrr::list_rbind() |>
        make_stars(drivers, vc) |>
        export_stars(output, "indirect", length(vc))

      # Net effects
      dir(out$net, full.names = TRUE) |>
        purrr::map(vroom::vroom) |>
        purrr::list_rbind() |>
        make_stars(drivers, vc) |>
        export_stars(output, "net", length(vc))

      # # Effects / km2
      # out <- paste0(output,"/cekm/")
      # chk_create(out)
      # get_cekm_ncea(motif_effects, vc) |>
      # vroom::vroom_write(paste0(out,"cekm.csv"), delim = ",")
    } else {
      NULL
    }
  }
}

# ==============================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------------------------------------------------------
export_stars <- function(dat, out, metric, n) {
  # Create output
  out <- paste0(out, "/", metric, "/")
  chk_create(out)
  nm <- names(dat)

  # Export
  for (i in seq_len(n)) {
    stars::write_stars(dat[i], file.path(out, paste0(nm[i], ".tif")))
  }
}
