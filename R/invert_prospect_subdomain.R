#' Function performing PROSPECT inversion on a spectral subset combining an
#' initial feature subset (initial_features) and an additional feature subset (additional_features)
#'
#' @param additional_features spectral bands (in nm) to be added to the initial_features
#' @param refl numeric. matrix of reflectances (n spectral bands x p samples)
#' @param tran numeric. matrix of transmittances (n spectral bands x p samples)
#' @param lambda numeric. spectral bands corresponding to reflectance and transmittance measurements
#' @param measured_bioch numeric. value of biophysical/biochemical parameter to estimate for each sample
#' @param target_parm character. name of the parameter.
#' Should be picked between "chl", "car", "ant", "brown", "ewt", "prot", "cbc", "n_struct"
#' @param parms_to_estimate  character vector. Parameters to estimate (can be 'ALL')
#' @param initial_features Initial feature set (in nm)
#' @param prospect_version  character. Version of prospect model used for the inversion: 'D' or 'PRO'
#' @param options list.
#' - init_values data.frame. Default values of PROSPECT parameters. During optimization,
#' they are used either as initialization values for parameters to estimate,
#' or as fix values for other parameters.
#' Parameters not compatible with PROSPECT_version are not taken into account.
#' - merit_function character. name of the function to be used as merit function
#' with given criterion to minimize (default = RMSE)
#' - xlub data.frame. Boundaries of the parameters to estimate.
#' The data.frame must have columns corresponding to \code{Parms2Estimate}
#' first line being the lower boundaries and second line the upper boundaries.
#' - estimate_brown_pigments boolean. should brown pigments be accounted for during inversion?
#' - estimate_alpha boolean. should alpha be accounted for during inversion?
#' - verbose boolean. set true to get info about adjustment of tolerance or initialization
#' - progressBar boolean. show progressbar?
#'
#' @return list containing estimated value of the target and RMSE compared to measured value
#' @importFrom pracma rmserr
#' @export

invert_prospect_subdomain <- function(additional_features, refl, tran,
                                      lambda, measured_bioch,
                                      target_parm, parms_to_estimate, initial_features,
                                      prospect_version = "D", options = NULL){

  options <- set_options_prospect(fun = 'invert_prospect_subdomain',
                                  options = options)
  spectral_domain <- c(initial_features, additional_features)
  # Fit spectral data to match PROSPECT with user optical properties
  subset_lop <- fit_spectral_data(spec_prospect = options$spec_prospect,
                                  lambda = lambda, refl = refl, tran = tran,
                                  user_domain = spectral_domain)
  # Invert PROSPECT with optimal spectral information
  estimated_val <- c()
  for (i in seq_len(length.out = subset_lop$nbSamples)){
    if (options$verbose)
      print(i)
    options$verbose <- FALSE
    res <- invert_prospect(refl = subset_lop$refl[,i],
                           tran = subset_lop$tran[,i],
                           prospect_version = prospect_version,
                           parms_to_estimate = parms_to_estimate,
                           options = options)
    # only save results for variable of interest
    estimated_val <- c(estimated_val,res[[target_parm]])
  }
  Perf <- pracma::rmserr(measured_bioch,estimated_val,summary = FALSE)$rmse
  rm(refl,tran,subset_lop)
  gc()
  return(list('Estimated' =  estimated_val, 'Perf' =  Perf))
}
