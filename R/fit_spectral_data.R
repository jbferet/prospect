#' This function adapts spec_prospect accordingly to experimental data
#' or to a spectral domain defined by user_domain
#' @param lambda numeric. Spectral bands corresponding to experimental data
#' @param spec_prospect list. Includes optical constants: refractive index,
#' specific absorption coefficients and corresponding spectral bands
#' @param refl numeric. Measured reflectance data
#' @param tran numeric. Measured transmittance data
#' @param user_domain numeric. either Lower and upper bounds for domain of
#' interest (optional) or list of spectral bands of interest
#' @param ul_bounds boolean. set to TRUE if user_domain only includes lower and
#' upper band, set to FALSE if user_domain is a list of spectral bands (in nm)
#'
#' @return list including spectral properties at the new resolution
#' @import dplyr
#' @export

fit_spectral_data <- function(lambda, spec_prospect = NULL,
                              refl = NULL, tran = NULL,
                              user_domain = NULL, ul_bounds = FALSE) {
  # default: simulates leaf optics using full spectral range
  if (is.null(spec_prospect))
    spec_prospect <- prospect::spec_prospect_full_range
  # convert refl and tran into dataframe if needed
  if (inherits(refl, what = c('numeric', 'matrix')))
    refl <- data.frame(refl)
  if (inherits(tran, what = c('numeric', 'matrix')))
    tran <- data.frame(tran)
  # if (class(refl)[1]%in%c('numeric', 'matrix')) refl <- data.frame(refl)
  # if (class(tran)[1]%in%c('numeric', 'matrix')) tran <- data.frame(tran)
  # Adjust LOP: check common spectral domain between PROSPECT and leaf optics
  if (!is.null(refl))
    refl <- refl %>% filter(lambda%in%spec_prospect$lambda)
  if (!is.null(tran))
    tran <- tran %>% filter(lambda%in%spec_prospect$lambda)
  lambda <- lambda[lambda%in%spec_prospect$lambda]
  # Adjust PROSPECT
  lb <- lambda
  spec_prospect  <- spec_prospect %>% filter(spec_prospect$lambda%in%lb)
  # if user_domain is defined
  if (is.null(user_domain))
    user_domain <- lambda
  if (ul_bounds==TRUE)
    user_domain <- seq(min(user_domain), max(user_domain))
  if (!is.null(refl))
    refl <- refl %>% filter(lambda%in%user_domain)
  if (!is.null(tran))
    tran <- tran %>% filter(lambda%in%user_domain)
  lambda <- lambda[lambda%in%user_domain]
  # Adjust PROSPECT
  spec_prospect  <- spec_prospect %>% filter(spec_prospect$lambda%in%user_domain)
  if (any(!user_domain%in%lambda))
    message('leaf optics out of range defined by user_domain')
  RT <- reshape_lop_for_inversion(refl = refl, tran = tran,
                                  spec_prospect = spec_prospect)
  return(list("spec_prospect" = spec_prospect, "lambda" = lambda,
              "refl" = RT$refl, "tran" = RT$tran, "nb_samples" = RT$nb_samples))
}

#' @rdname FitSpectralData-deprecated
#' @export
#'
FitSpectralData <- function(lambda, SpecPROSPECT = NULL,
                            Refl = NULL, Tran = NULL,
                            UserDomain = NULL, UL_Bounds = FALSE) {

  .Deprecated("FitSpectralData")
  fit_spectral_data(lambda = lambda, spec_prospect = SpecPROSPECT,
                    refl = Refl, tran = Tran, user_domain = UserDomain,
                    ul_bounds = UL_Bounds)
}
