#' Trophic sensitivity
#'
#' A dataset containing the simplified values for trophic sensitivity as measured by Beauchesne et al. 2021 (DOI: 10.1111/ele.13841) and used by Beauchesne et al. 2020 for network-scale cumulative effects assessments
#'
#' @format ## `trophic_sensitivity`
#' A data frame with 124 rows and 10 columns:
#' \describe{
#'   \item{Motif}{Name of motifs: apparent competition (ap); disconnected (di); exploitative competition (ex); omnivory (om); partially connected (pa); tri-trophic interaction (tt)}
#'   \item{Species}{Position of species in motif (x,y,z)}
#'   \item{px, py, pz}{Pathways of effects, whether species x, y or z are affected by disturbances}
#'   \item{Sensitivity}{Trophic sensitivity scaled between 0 and 1}
#'   \item{sensitivity_original}{Original value of trophic sensitivity}
#'   \item{pathID}{Unique identifier of pathway of effect}
#'   \item{speciesID}{Numeric ID for species position in motifs}
#'   \item{motifID}{Numeric ID for motifs}
#' }
#' @source <https://github.com/david-beauchesne/FoodWeb-MultiStressors/blob/master/Data/vulnerability.RData>
#' @source <https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.13841>
#' @source <https://semaphore.uqar.ca/id/eprint/1922/1/David_Beauchesne_decembre2020.pdf>
"trophic_sensitivity"