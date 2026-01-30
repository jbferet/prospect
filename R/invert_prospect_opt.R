#' This function uses optimal configuration identified by Spafford et al. (2020) to estimate leaf chemistry
#' prior information on N is provided if only reflectance or only Transmittance is available
#' @param lambda  numeric. spectral bands corresponding to experimental data
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param refl  numeric. Measured reflectance data
#' @param tran  numeric. Measured transmittance data
#' @param prospect_version  character. Version of prospect model used for the
#' inversion: 'D', 'PRO'.
#' @param parms_to_estimate  character vector. Parameters to estimate (can be 'ALL')
#' @param options list.
#' - init_values data.frame. Default values of PROSPECT parameters. During optimization,
#' they are used either as initialization values for parameters to estimate,
#' or as fix values for other parameters.
#' Parameters not compatible with PROSPECT_version are not taken into account.
#' - merit_function character. name of the function to be used as merit function
#' with given criterion to minimize (default = RMSE)
#' - xlub data.frame. Boundaries of the parameters to estimate.
#' The data.frame must have columns corresponding to \code{parms_to_estimate}
#' first line being the lower boundaries and second line the upper boundaries.
#' - estimate_brown_pigments boolean. should brown pigments be accounted for during inversion?
#' - estimate_alpha boolean. should alpha be accounted for during inversion?
#' - verbose boolean. set true to get info about adjustment of tolerance or initialization
#' - progressBar boolean. show progressbar?
#'
#' @return Nprior vector corresponding to teh prior estimation of N based on R only or T only
#' @importFrom stats lm runif
#' @importFrom progress progress_bar
#' @export
#'
invert_prospect_opt <- function(lambda, refl = NULL, tran = NULL,
                                parms_to_estimate = 'ALL',
                                prospect_version = 'D', options = NULL){

  options <- set_options_prospect(fun = 'invert_prospect_opt', options = options)
  spec_prospect <- options$spec_prospect
  merit_function <- options$merit_function
  xlub <- options$xlub
  init_values <- options$init_values
  progressBar <- options$progressBar
  verbose <- options$verbose
  estimate_alpha <- options$estimate_alpha
  estimate_brown_pigments <- options$estimate_brown_pigments

  if (is.null(spec_prospect))
    spec_prospect <- prospect::spec_prospect_full_range
  # versions implemented in the package
  if (!prospect_version %in% c('D','PRO'))
    print_msg(cause = 'WrongVersion')
  # if brown pigments required
  if ('brown' %in% parms_to_estimate) {
    print_msg(cause = 'NoBrown_OPT')
    parms_to_estimate <- parms_to_estimate[-match('brown',parms_to_estimate)]
  }
  # define optimal domain for the different constituents
  optimal_domain <- optimal_domain_R <- optimal_domain_T <- prospect::optimal_domains_rt
  # check if list of parameters applicable to PROSPECT version
  parms_checked <- check_prospect_parms(prospect_version = prospect_version,
                                        parms_to_estimate = parms_to_estimate,
                                        estimate_brown_pigments = FALSE,
                                        estimate_alpha = FALSE, xlub = xlub,
                                        init_values = init_values)
  parms_to_estimate <- parms_checked$parms_to_estimate
  ant_init <- parms_checked$init_values$ant
  if (prospect_version == 'PRO' & 'lma' %in% parms_to_estimate)
    print_msg(cause = 'version_PROSPECT_Invert_PROSPECT_OPT')
  # adjust leaf optics and convert to data frame if needed
  adjusted_lop <- fit_spectral_data(spec_prospect = spec_prospect,
                                    lambda = lambda, refl = refl, tran = tran)
  # prior assessment of N if R or T missing
  Nprior <- set_n_values(refl = adjusted_lop$refl,
                         tran = adjusted_lop$tran,
                         spec_prospect = adjusted_lop$spec_prospect)
  parm_estimated <- list()
  if (is.null(adjusted_lop$refl) | is.null(adjusted_lop$tran))
    parm_estimated$n_struct <- Nprior
  if ('n_struct' %in% parms_to_estimate)
    parms_to_estimate <- parms_to_estimate[-which(parms_to_estimate=='n_struct')]
  # initialize parameters for each sample (useful when using prior n_struct)
  list_init <- set_init_parm(parms_to_estimate = parms_to_estimate,
                             parm_estimated = parm_estimated,
                             prospect_version = prospect_version,
                             Nprior = Nprior, ant_init = ant_init,
                             optimal_domain = optimal_domain,
                             refl = adjusted_lop$refl, tran = adjusted_lop$tran)
  parm_estimated <- list_init$parm_estimated

  for (parm in list_init$parms_to_estimate){
    if (length(parm_estimated[[parm]])==0){
      # Fit spectral data to match PROSPECT with user optical properties
      OptLOP <- fit_spectral_data(spec_prospect = adjusted_lop$spec_prospect,
                                  lambda = adjusted_lop$lambda,
                                  refl = adjusted_lop$refl,
                                  tran = adjusted_lop$tran,
                                  user_domain = optimal_domain[[parm]],
                                  ul_bounds = list_init$ul_bounds[[parm]])
      parm_display <- parm
      if (parm=='ewt' & 'lma'%in%parms_to_estimate)
        parm_display <- 'ewt & lma'
      if (progressBar==TRUE){
        pb <- progress::progress_bar$new(
          format = paste("Estimation of ",parm_display," [:bar] :percent in :elapsedfull , estimated time remaining :eta"),
          total = adjusted_lop$nb_samples, clear = FALSE, width= 100)
      }
      for (i in seq_len(adjusted_lop$nb_samples)){
        options <- set_options_prospect(fun = 'invert_prospect', options = options)
        options$init_values <- list_init$init_values[[parm]][i,]
        options$spec_prospect <- OptLOP$spec_prospect
        res <- invert_prospect(refl = OptLOP$refl[[i]],
                               tran = OptLOP$tran[[i]],
                               prospect_version = list_init$prospect_version[[parm]],
                               parms_to_estimate = list_init$parms_to_estimate_tmp[[parm]],
                               options = options)
        parm_estimated[[parm]][i] <- res[[parm]]
        if (parm=='ewt' & !is.na(match('lma',parms_to_estimate)))
          parm_estimated$lma[i] <- res$lma
        if (progressBar==TRUE) pb$tick()
      }
    }
  }
  parm_estimated <- data.frame(parm_estimated)
  return(parm_estimated)
}



#' @rdname Invert_PROSPECT_OPT-deprecated
#' @export
#'
Invert_PROSPECT_OPT <- function(SpecPROSPECT = NULL,
                                Refl = NULL, Tran = NULL,
                                InitValues = data.frame(chl = 40,
                                                        car = 10,
                                                        ant = 0.1,
                                                        brown = 0.0,
                                                        ewt = 0.01,
                                                        lma = 0.01,
                                                        prot = 0.001,
                                                        cbc = 0.009,
                                                        n_struct = 1.5,
                                                        alpha = 40),
                                Parms2Estimate = "ALL",
                                PROSPECT_version = "D",
                                MeritFunction = "merit_prospect_rmse",
                                xlub = data.frame(chl = c(1e-4, 150),
                                                  car = c(1e-4, 25),
                                                  ant = c(0, 50),
                                                  brown = c(0, 4),
                                                  ewt = c(1e-8, 0.1),
                                                  lma = c(1e-8, 0.06),
                                                  prot = c(1e-7, 0.006),
                                                  cbc = c(1e-6, 0.054),
                                                  n_struct = c(0.5, 4),
                                                  alpha = c(10, 90)),
                                Est_Brown_Pigments = FALSE, Est_alpha = FALSE,
                                verbose = FALSE, progressBar = TRUE) {

  .Deprecated("Invert_PROSPECT_OPT")
  options <- set_options_prospect(fun = 'invert_prospect_opt')
  options$spec_prospect <- SpecPROSPECT
  options$merit_function <- MeritFunction
  options$xlub <- xlub
  options$init_values <- InitValues
  options$progressBar <- progressBar
  options$verbose <- verbose
  options$estimate_alpha <- Est_alpha
  options$estimate_brown_pigments <- Est_Brown_Pigments
  invert_prospect_opt(refl = Refl, tran = Tran, parms_to_estimate = Parms2Estimate,
                      prospect_version = PROSPECT_version, options = options)
}
