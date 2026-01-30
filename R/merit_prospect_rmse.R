#' Merit function for PROSPECT inversion
#'
#' @param x numeric. Vector of input variables to estimate
#' @param spec_prospect list. Includes optical constants refractive index,
#' specific absorption coefficients and corresponding spectral bands
#' @param refl  numeric. measured reflectance data
#' @param tran  numeric. measured Transmittance data
#' @param input_prospect dataframe. set of PROSPECT input variables
#' @param parms_to_estimate  numeric. location of variables from input_prospect
#' to be estimated through inversion
#'
#' @return fc estimates of the parameters
#' @export

merit_prospect_rmse <- function(x, spec_prospect, refl, tran, input_prospect,
                                parms_to_estimate){
  x[x < 0] <- 0
  input_prospect[parms_to_estimate] <- x
  RT <- do.call(what = 'prospect',
                args = c(list('spec_prospect' = spec_prospect),
                         input_prospect))
  # RT <- do.call("prospect", args = list(spec_prospect = spec_prospect,
  #                                       input_prospect = input_prospect,
  #                                       check = F))
  fcr <- fct <- 0
  if (!is.null(refl))
    fcr <- sqrt(sum((refl - RT$reflectance)**2) / length(RT$reflectance))
  if (!is.null(tran))
    fct <- sqrt(sum((tran - RT$transmittance)**2) / length(RT$transmittance))
  fc <- fcr + fct
  return(fc)
}
