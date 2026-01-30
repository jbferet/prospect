#' This function defines a regression model to estimate N from R only or T only
#' @param lambda  numeric. spectral bands corresponding to experimental data
#' @param spec_prospect list. Includes optical constants refractive index,
#' specific absorption coefficients and corresponding spectral bands
#' @param refl numeric. Measured reflectance data
#' @param tran numeric. Measured Transmittance data
#' @param OptWL_R list. optimal wavelengths used to estimate N from R only
#' @param OptWL_T list. optimal wavelengths used to estimate N from T only
#'
#' @return Nprior prior estimation of N with R or T only
#' @importFrom stats lm runif
#' @export
#'
get_N_prior <- function(lambda, spec_prospect = NULL, refl = NULL, tran = NULL,
                        OptWL_R = list(NIR = 800,SWIR = 1131),
                        OptWL_T = list(NIR = 753,SWIR = 1121)) {

  if (is.null(spec_prospect))
    spec_prospect <- prospect::spec_prospect_full_range
  # if prior information based on reflectance
  if (is.null(tran)) {
    # if required spectral bands in the original data
    if (OptWL_R$SWIR %in% spec_prospect$lambda) {
      OptWL <- OptWL_R$SWIR
      # else if close to spectral bands of original data
    } else if (!OptWL_R$SWIR %in% spec_prospect$lambda &
               min(abs(OptWL_R$SWIR-spec_prospect$lambda))<10) {
      OptWL <- spec_prospect$lambda[closest_band(target_wl = OptWL_R$SWIR,
                                                 wl = spec_prospect$lambda)]
      # else if NIR band available
    } else if (!OptWL_R$SWIR %in% spec_prospect$lambda &
               OptWL_R$NIR %in% spec_prospect$lambda) {
      print_msg(cause = 'Nprior_noRefl')
      OptWL <- OptWL_R$NIR
      # else if close to NIR band available
    } else if (!OptWL_R$NIR %in% spec_prospect$lambda &
               min(abs(OptWL_R$NIR-spec_prospect$lambda))<10) {
      print_msg(cause = 'Nprior_noRefl')
      OptWL <- spec_prospect$lambda[closest_band(target_wl = OptWL_R$NIR,
                                                 wl = spec_prospect$lambda)]
      # else if NIR band available
    } else if (!OptWL_R$SWIR %in% spec_prospect$lambda &
               !OptWL_R$NIR %in% spec_prospect$lambda) {
      print_msg(cause = 'Nprior_noRefl_atAll')
    }
  } else if (is.null(refl)) {
    if (OptWL_T$SWIR %in% spec_prospect$lambda) {
      OptWL <- OptWL_T$SWIR
    } else if (!OptWL_T$SWIR %in% spec_prospect$lambda &
               min(abs(OptWL_T$SWIR-spec_prospect$lambda))<10) {
      OptWL <- spec_prospect$lambda[closest_band(target_wl = OptWL_T$SWIR,
                                                 wl = spec_prospect$lambda)]
    } else if (!OptWL_T$SWIR %in% spec_prospect$lambda &
               OptWL_T$NIR %in% spec_prospect$lambda) {
      print_msg(cause = 'Nprior_noTran')
      OptWL <- OptWL_T$NIR
    } else if (!OptWL_T$NIR %in% spec_prospect$lambda &
               min(abs(OptWL_T$NIR-spec_prospect$lambda))<10) {
      print_msg(cause = 'Nprior_noTran')
      OptWL <- spec_prospect$lambda[closest_band(target_wl = OptWL_T$NIR,
                                                 wl = spec_prospect$lambda)]
    } else if (!OptWL_T$SWIR %in% spec_prospect$lambda &
               !OptWL_T$NIR %in% spec_prospect$lambda) {
      print_msg(cause = 'Nprior_noTran_atAll')
    }
  }

  # get the subdomain corresponding to OptWL
  subset_data <- fit_spectral_data(spec_prospect = spec_prospect, lambda = lambda,
                                   refl = refl, tran = tran,
                                   user_domain = c(OptWL, OptWL), ul_bounds = TRUE)

  # create a LUT using only the spectral band of interest
  nb_simulations <- 1000
  input_prospect <- data.frame('chl' = 0.5 + 100 * runif(nb_simulations),
                               'car' = 0.5 + 20 * runif(nb_simulations),
                               'ewt' = 0.001 + 0.02 * runif(nb_simulations),
                               'lma' = 0.001 + 0.01 * runif(nb_simulations),
                               'n_struct' = 1 + 1.5 * runif(nb_simulations))
  LUT <- prospect_lut(input_prospect = input_prospect,
                      spec_prospect = subset_data$spec_prospect)
  # fit a linear model between
  if (is.null(tran)) {
    Ratio <- data.frame(c(unlist(LUT$reflectance))) / (1 - data.frame(c(unlist(LUT$reflectance))))
    Ratio_Meas <- data.frame(c(unlist(subset_data$refl))) / (1 - data.frame(c(unlist(subset_data$refl))))
  } else if (is.null(refl)) {
    Ratio <- (1 - data.frame(c(unlist(LUT$transmittance)))) / data.frame(c(unlist(LUT$transmittance)))
    Ratio_Meas <- (1 - data.frame(c(unlist(subset_data$tran)))) / data.frame(c(unlist(subset_data$tran)))
  }
  N_Model <- lm(c(LUT$input_prospect$n_struct) ~ c(Ratio[[1]]))
  NpriorMOD <- N_Model$coefficients[2] * Ratio + N_Model$coefficients[1]
  Nprior <- N_Model$coefficients[2] * Ratio_Meas + N_Model$coefficients[1]
  rownames(Nprior) <- NULL
  colnames(Nprior) <- 'n_struct'
  return(Nprior)
}

#' @rdname Get_Nprior-deprecated
#' @export
#'
Get_Nprior <- function(lambda, SpecPROSPECT = NULL, Refl = NULL, Tran = NULL,
                       OptWL_R = list(NIR = 800,SWIR = 1131),
                       OptWL_T = list(NIR = 753,SWIR = 1121)) {

  .Deprecated("Get_Nprior")
  get_N_prior(lambda = lambda, spec_prospect = SpecPROSPECT,
              refl = Refl, tran = Tran, OptWL_R = OptWL_R, OptWL_T = OptWL_T)
}
