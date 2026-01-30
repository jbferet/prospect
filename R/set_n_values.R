#' Function set N values if reflectance or transmittance is missing
#'
#' @param refl dataframe. reflectance data
#' @param tran dataframe. transmittance data
#' @param spec_prospect dataframe. spectral properties used in PROSPECT
#'
#' @return Nprior dataframe. N values corresponding to refl or tran
#' @export
#'
set_n_values <- function(refl, tran, spec_prospect){

  if (!is.null(refl)){
    Nprior <- data.frame('N' = rep(1.5,ncol(refl)))
  } else if (!is.null(tran)){
    Nprior <- data.frame('N' = rep(1.5,ncol(tran)))
  } else {
    message('no leaf optical properties provided to set_n_values')
  }
  if (is.null(refl) | is.null(tran)){
    # compute prior estimate of N
    message('computing prior estimation of N as both R & T are not provided')
    Nprior <- get_N_prior(spec_prospect = spec_prospect,
                          lambda = spec_prospect$lambda,
                          refl = refl, tran = tran)
  }
  return(Nprior)
}

#' @rdname SetNValues-deprecated
#' @export
#'
SetNValues <- function(Refl, Tran, SpecPROSPECT){

  .Deprecated("SetNValues")
  set_n_values(refl = Refl, tran = Tran, spec_prospect = SpecPROSPECT)
}
