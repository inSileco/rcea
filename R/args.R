# ------------------------------------------------------------------------------
# List of arguments used throughout the package
# NOTE: Documented here in order to avoid unnecessary repetition in documentation and allow for easy modifications in the future if necessary.
arguments <- function(arg) {
  dat <- list()
  dat$drivers <- "@param drivers distribution and intensity of environmental drivers as stars object"
  dat$vc <- "@param vc distribution of valued components as stars object"
  dat$sensitivity <- "@param sensitivity matrix of environmental drivers and valued component, with same name as those used in `drivers` and `vc`"
  dat$metaweb <- "@param metaweb matrix of valued component by valued component describing the binary interations structuring the network of valued components"
  dat$trophic_sensitivity <- "@param trophic_sensitivity data.frame of trophic sensitivities, default from Beauchesne. Available as data package with `data(trophic_sensitivity)`"
  dat$weights <- "@param w_d,w_i weight for the direct (`w_d`) and indirect (`w_i`) modules when calculating network-scale cea scores; w_d + 2*w_i should be equal to 1."
  dat$output <- "@param output relative path to export results of assessment."
  dat$motif_effects <- "@param motif_effects TODO"
  dat$output_format <- "@param output_format output format, one of `geotiff` or `COG`."
  dat$focus <- "@param focus named argument, string with name of valued component on which to perform the assessment."

  dat[names(dat) %in% arg] |>
    unlist()
}
