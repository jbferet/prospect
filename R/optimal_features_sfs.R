#' Function performing optimal feature selection based on sequential forward feature slection
#'
#' @param refl numeric. matrix of reflectances (n spectral bands x p samples)
#' @param tran numeric. matrix of transmittances (n spectral bands x p samples)
#' @param lambda numeric. spectral bands corresponding to reflectance and transmittance measurements
#' @param measured_bioch numeric. value of biophysical/biochemical parameter to estimate for each sample
#' @param target_parm character. name of the parameter.
#' Should be picked between "chl", "car", "ant", "brown", "ewt", "prot", "cbc", "n_struct"
#' @param parms_to_estimate  character vector. Parameters to estimate (can be 'ALL')
#' @param spectral_domain vector. defines minimum and maximum wavelengths of the spectral domain (in nm).
#' Assumes 1nm spectral sampling
#' @param spectral_width vector. width of individual spectral features (in nm)
#' @param number_features vector. number of features to be identified
#' @param prospect_version character.
#' @param nbCPU numeric. defines number of CPU for multithread processing
#' @param continue boolean. set to TRUE if the function has already been run for a lower value of number_features
#' @param already_done list. contains output of function optimal_features_SFS previouly run with lower value of number_features
#' @param options list.
#'
#' @return list containing :
#' - spectral_features = the Spectral Features identified as most relevant for estimation of target_parm (vector)
#' - spectral_features_List = the Spectral Features identified as most relevant for estimation of target_parm (list)
#' - estimated_parms = the estimated target_parm parameter for each additional spectral feature
#' - RMSE = the minimum RMSE corresponding to the estimation of target_parm for each asditional spectral feature
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @export

optimal_features_sfs <- function(refl = NULL, tran = NULL, lambda,
                                 measured_bioch, target_parm,
                                 parms_to_estimate = "ALL",
                                 spectral_domain, spectral_width,
                                 number_features, prospect_version = 'D',
                                 nbCPU = 1, continue = FALSE,
                                 already_done = NULL, options = NULL){

  # split spectral_domain into subdomains based on spectral_width
  full_domain <- seq(spectral_domain[1],spectral_domain[2])
  x <- seq_along(full_domain)
  subdomains <- split(full_domain,ceiling(x/spectral_width))
  nb_subdomains <- length(subdomains)
  # spectral bands identified as optimal
  if (continue){
    # number of features to start from
    featStart <- length(already_done$RMSE)+1
    estimated_all <- already_done$estimated_parms
    initial_features_list <- already_done$spectral_features_List
    initial_features <- already_done$spectral_features
    minRMSE <- already_done$RMSE

    # number of features to start from
    for (i in seq_len(length(already_done$RMSE))){
      for (j in rev(seq_len(length.out = length(subdomains)))){
        if (already_done$spectral_features_List[[i]][1] %in% subdomains[[j]])
          subdomains[[j]] <- NULL
      }
    }
  } else {
    featStart <- 1
    # for incremental number of spectral features from 1 to number_features
    estimated_all <- initial_features_list <- list()
    initial_features <- minRMSE <- c()
  }

  if (number_features>=featStart){
    for (feat in featStart:number_features){
      message(paste('Identify Feature #',feat))
      Perf <- c()
      plan(multisession, workers = nbCPU) ## Parallelize using four cores
      schedule <- ceiling(length(subdomains)/nbCPU)

      estimated_parm <- future_lapply(X = subdomains,
                                      FUN = invert_prospect_subdomain,
                                      refl = refl, tran = tran,
                                      lambda = lambda, measured_bioch = measured_bioch,
                                      target_parm = target_parm,
                                      parms_to_estimate  = parms_to_estimate,
                                      initial_features = initial_features,
                                      prospect_version = prospect_version,
                                      options = options,
                                      future.scheduling = schedule)
      plan(sequential)
      Perf <- c()
      for (i in seq_len(length(estimated_parm)))
        Perf <- c(Perf,estimated_parm[[i]]$Perf)
      BestPerf <- which(Perf==min(Perf, na.rm = TRUE))
      estimated_all[[feat]] <- estimated_parm[[BestPerf[1]]]$Estimated
      minRMSE <- c(minRMSE,min(Perf,na.rm = TRUE))
      print(min(Perf,na.rm = TRUE))
      message('Spectral domain identified:')
      print(subdomains[[BestPerf[1]]])
      initial_features <- c(initial_features,subdomains[[BestPerf[1]]])
      initial_features_list[[feat]] <- subdomains[[BestPerf[1]]]
      subdomains[[BestPerf[1]]] <- NULL
    }
  }
  return(list('spectral_features' = initial_features,
              'spectral_features_List' = initial_features_list,
              'estimated_parms' = estimated_all, 'RMSE' = minRMSE))
}
