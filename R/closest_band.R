#' This function identifies the closest spectral band among a selection of bands
#' @param target_wl  numeric. central wavelength of spectral band of interest
#' @param wl numeric. vector of spectral bands
#'
#' @return selected_band rank of selected band
#' @export
#'
closest_band <- function(target_wl, wl){
  selected_band <- which(abs(target_wl-wl)==min(abs(target_wl-wl)))
  return(selected_band)
}
