#' Function to set initial parameterization when performing PROSPECT inversion
#' using optimal configuration
#'
#' @param parms_to_estimate  character. Parameters to estimate
#' @param parm_estimated character. PROSPECT parameters to be estimated
#' @param prospect_version  character. Version of prospect model used for the
#' inversion: 'D' or 'PRO'
#' @param Nprior numeric prior estimation of N
#' @param ant_init numeric prior estimation of ANT
#' @param optimal_domain vector. optimal spectral domain
#' @param refl dataframe reflectance data
#' @param tran dataframe transmittance data
#'
#' @return list containing parms_to_estimate, parm_estimated, ul_bounds, prospect_version.
#' parms_to_estimate_tmp, init_values
#' @export
#'
set_init_parm <- function(parms_to_estimate, parm_estimated, prospect_version,
                          Nprior, ant_init, optimal_domain, refl, tran){

  ul_bounds <- prospect_version_tmp <-
    parms_to_estimate_tmp <- init_values <- list()
  for (parm in parms_to_estimate){
    parm_estimated[[parm]] <- c()
    init_values[[parm]] <- list()
    if (length(optimal_domain[[parm]])==2){
      ul_bounds[[parm]] <- TRUE
    } else {
      ul_bounds[[parm]] <- FALSE
    }
    if (parm == "ant"){
      print_msg(cause = 'NoOpt_ANT')
      prospect_version_tmp[[parm]] <- prospect_version
      parms_to_estimate_tmp[[parm]] <- c('chl', 'car', 'ant')
      init_values[[parm]] <- data.frame(chl = 40, car = 10, ant = ant_init,
                                       brown = 0, ewt = 0.01, lma = 0.01,
                                       n_struct = Nprior)
    }
    if (parm == "chl" | parm == "car"){
      prospect_version_tmp[[parm]] <- prospect_version
      parms_to_estimate_tmp[[parm]] <- c('chl', 'car', 'ant')
      init_values[[parm]] <- data.frame(chl = 40, car = 10, ant = ant_init,
                                       brown = 0, ewt = 0.01, lma = 0.01,
                                       n_struct = Nprior)
    }
    if (parm == "ewt"){
      prospect_version_tmp[[parm]] <- prospect_version
      if (!prospect_version =='PRO'){
        parms_to_estimate_tmp[[parm]] <- c('ewt', 'lma')
        init_values[[parm]] <- data.frame(chl = 0, car = 0, ant = 0, brown = 0,
                                         ewt = 0.01, lma = 0.01,
                                         n_struct = Nprior)
      } else if (prospect_version =='PRO'){
        parms_to_estimate_tmp[[parm]] <- c('ewt', 'prot', 'cbc')
        init_values[[parm]] <- data.frame(chl = 0, car = 0, ant = 0, brown = 0,
                                         ewt = 0.01, lma = 0.00, prot = 0.001,
                                         cbc = 0.009, n_struct = Nprior)
      }
    }
    if (parm == "lma" & !prospect_version == 'PRO'){
      prospect_version_tmp[[parm]] <- prospect_version
      parms_to_estimate_tmp[[parm]] <- c('ewt', 'lma')
      init_values[[parm]] <- data.frame(chl = 0, car = 0, ant = 0, brown = 0,
                                       ewt = 0.01, lma = 0.01, n_struct = Nprior)
    }
    if (parm == "prot" | parm == "cbc"){
      prospect_version_tmp[[parm]] <- 'PRO'
      parms_to_estimate_tmp[[parm]] <- c('ewt', 'prot', 'cbc')
      init_values[[parm]] <- data.frame(chl = 0, car = 0, ant = 0, brown = 0,
                                       ewt = 0.01, lma = 0.00, prot = 0.001,
                                       cbc = 0.009, n_struct = Nprior)
    }
    if (!is.null(refl) & !is.null(tran))
      parms_to_estimate_tmp[[parm]] <- c(parms_to_estimate_tmp[[parm]],'n_struct')
  }
  return(list('parms_to_estimate' = parms_to_estimate,
              'parm_estimated' = parm_estimated,
              'ul_bounds' = ul_bounds,
              'prospect_version' = prospect_version_tmp,
              'parms_to_estimate_tmp' = parms_to_estimate_tmp,
              'init_values'=init_values))
}



#' @rdname SetInitParm-deprecated
#' @export
#'
SetInitParm <- function(Parms2Estimate, ParmEst, PROSPECT_version, Nprior,
                        ANTinit, OptDomain, Refl, Tran){

  .Deprecated("SetInitParm")
  set_init_parm(parms_to_estimate = Parms2Estimate, parm_estimated = ParmEst,
                prospect_version = PROSPECT_version, Nprior = Nprior,
                ant_init = ANTinit, optimal_domain = OptDomain,
                refl = Refl, tran = Tran)
}
