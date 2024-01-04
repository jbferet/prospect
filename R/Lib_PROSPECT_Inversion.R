# ==============================================================================
# prospect
# Lib_PROSPECT_Inversion.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de Boissieu <fdeboiss@gmail.com>
# Copyright 2020/04 Jean-Baptiste FERET
# ==============================================================================
# This Library includes functions dedicated to PROSPECT inversion
# ==============================================================================

#' Performs PROSPECT inversion using iterative optimization
#' in order to define estimate a set of user defined parameters
#' @param SpecPROSPECT list. Includes optical constants refractive index,
#' specific absorption coefficients and corresponding spectral bands
#' @param Refl  numeric. Measured reflectance data
#' @param Tran  numeric. Measured Transmittance data
#' @param Parms2Estimate  character. Parameters to estimate (can be 'ALL')
#' @param InitValues  data.frame. Default values of PROSPECT parameters. During optimization,
#' they are used either as initialization values for parameters to estimate,
#' or as fix values for other parameters.
#' Parameters not compatible with PROSPECT_version are not taken into account.
#' @param PROSPECT_version  character. prospect version used for inversion: 'D' or 'PRO'
#' See details.
#' @param MeritFunction  character. name of the function to be used as merit function
#' with given criterion to minimize (default = RMSE)
#' @param xlub data.frame. Boundaries of the parameters to estimate.
#' The data.frame must have columns corresponding to \code{Parms2Estimate}
#' first line being the lower boundaries and second line the upper boundaries.
#' @param Est_Brown_Pigments boolean. should brown pigments be accounted for during inversion?
#' @param Est_alpha boolean. should alpha be accounted for during inversion?
#' @param verbose boolean. set true to get info about adjustment of tolerance or initialization
#' @param progressBar boolean. show progressbar?
#'
#'
#' @return OutPROSPECT estimated values corresponding to Parms2Estimate
#' @importFrom progress progress_bar
#' @details
#' Two versions of prospect are available for inversion.
#' The version is depending on the parameters taken into account:
#'
#' | Version  | D                                      | PRO
#' | :------: |:--------------------------------------:|:--------------------------------------:|
#' | CHL      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | CAR      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | ANT      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | BROWN    |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | EWT      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | LMA      |`r emojifont::emoji('white_check_mark')`|                                        |
#' | PROT     |                                        |`r emojifont::emoji('white_check_mark')`|
#' | CBC      |                                        |`r emojifont::emoji('white_check_mark')`|
#' | N        |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#'
#' Argument `InitValues` is expecting a default value for each of the parameters as well as an `alpha` value.
#' @importFrom pracma fmincon
#' @export
#' @md
Invert_PROSPECT <- function(SpecPROSPECT = NULL,
                            Refl = NULL, Tran = NULL,
                            InitValues = data.frame(
                              CHL = 40, CAR = 10,
                              ANT = 0.1, BROWN = 0.0,
                              EWT = 0.01, LMA = 0.01,
                              PROT = 0.001, CBC = 0.009,
                              N = 1.5, alpha = 40),
                            Parms2Estimate = "ALL",
                            PROSPECT_version = "D",
                            MeritFunction = "Merit_PROSPECT_RMSE",
                            xlub = data.frame(
                              CHL = c(1e-4, 150), CAR = c(1e-4, 25),
                              ANT = c(0, 50), BROWN = c(0, 1),
                              EWT = c(1e-8, 0.1), LMA = c(1e-8, 0.06),
                              PROT = c(1e-7, 0.006), CBC = c(1e-6, 0.054),
                              N = c(0.5, 4), alpha = c(10, 90)),
                            Est_Brown_Pigments = FALSE, Est_alpha = FALSE,
                            verbose = FALSE, progressBar = TRUE) {

  if (is.null(SpecPROSPECT)) SpecPROSPECT <- prospect::SpecPROSPECT_FullRange
  # check if list of parameters applicable to PROSPECT version
  parms_checked <- check_prospect_parms(PROSPECT_version = PROSPECT_version,
                                        Parms2Estimate = Parms2Estimate,
                                        Est_Brown_Pigments = Est_Brown_Pigments,
                                        Est_alpha = Est_alpha,
                                        xlub = xlub,
                                        InitValues = InitValues)
  # check if data class is compatible and convert into data.frame
  RT <- reshape_lop4inversion(Refl = Refl,
                              Tran = Tran,
                              SpecPROSPECT = SpecPROSPECT)
  OutPROSPECT <- list()
  if (progressBar==TRUE){
    pb <- progress::progress_bar$new(
      format = "Inverting PROSPECT [:bar] :percent in :elapsedfull, estimated time remaining :eta",
      total = RT$nbSamples, clear = FALSE, width= 100)
  }
  for (idsample in seq_len(RT$nbSamples)){
    res <- tryInversion(x0 = parms_checked$InitValues,
                        MeritFunction = MeritFunction,
                        SpecPROSPECT = SpecPROSPECT,
                        Refl = RT$Refl[[idsample]], Tran = RT$Tran[[idsample]],
                        Parms2Estimate = parms_checked$Parms2Estimate,
                        lb = parms_checked$lb, ub = parms_checked$ub,
                        verbose = verbose)
    if (NA %in% res$par){
      ModifyInit <- match(parms_checked$Parms2Estimate,
                          names(parms_checked$InitValues))
      updateInitValues <- parms_checked$InitValues
      updateInitValues[ModifyInit] <- 1.1*updateInitValues[ModifyInit]
      res <- tryInversion(x0 = updateInitValues,
                          MeritFunction = MeritFunction,
                          SpecPROSPECT = SpecPROSPECT,
                          Refl = RT$Refl[[idsample]], Tran =RT$Tran[[idsample]],
                          Parms2Estimate = parms_checked$Parms2Estimate,
                          lb = parms_checked$lb, ub = parms_checked$ub,
                          verbose = verbose)
    }
    names(res$par) <- parms_checked$Parms2Estimate
    OutPROSPECT[[idsample]] <- parms_checked$InitValues
    OutPROSPECT[[idsample]][names(res$par)] <- res$par
    if (progressBar==TRUE) pb$tick()
  }
  OutPROSPECT <- do.call(rbind,OutPROSPECT)
  return(OutPROSPECT)
}

#' Function handling error during inversion
#'
#' @param x0 numeric. Vector of input variables to estimate
#' @param MeritFunction  character. name of the merit function
#' @param SpecPROSPECT list. Includes optical constants refractive index,
#' specific absorption coefficients and corresponding spectral bands
#' @param Refl  numeric. measured reflectance data
#' @param Tran  numeric. measured Transmittance data
#' @param Parms2Estimate  character vector. Parameters to estimate
#' to be estimated through inversion
#' @param lb numeric. Lower bound
#' @param ub numeric. Upper bound
#' @param verbose boolean. Set TRUE for information on adjustment of tolerance
#' during inversion.
#'
#' @return res estimates of the parameters
#' @details
#' This function is based on \code{\link[pracma]{fmincon}}.
#' @importFrom pracma fmincon
#' @export

tryInversion <- function(x0, MeritFunction, SpecPROSPECT, Refl, Tran,
                         Parms2Estimate, lb, ub, verbose = FALSE) {

  res <-list('par' = NA*vector(length = length(Parms2Estimate)))
  TolRange <- seq(-14,-2,1)
  for (i in TolRange){
    Tolerance <- 10**(i)
    if (is.na(res$par[1])){
      res <- tryCatch(
        {
          res <- fmincon(
            x0 = as.numeric(x0[Parms2Estimate]),
            fn = MeritFunction, gr = NULL,
            SpecPROSPECT = SpecPROSPECT, Refl = Refl, Tran = Tran,
            Input_PROSPECT = x0, Parms2Estimate = Parms2Estimate,
            method = "SQP", A = NULL, b = NULL, Aeq = NULL, beq = NULL,
            lb = as.numeric(lb), ub = as.numeric(ub), hin = NULL, heq = NULL,
            tol = Tolerance, maxfeval = 2000, maxiter = 1000
          )
        },
        error = function(cond) {
          if (verbose) print_msg(cause = 'AdjustTol',
                                 args = list('Tolerance' = Tolerance))
          return(list('par' = NA*vector(length = length(Parms2Estimate))))
        },
        finally = {}
      )
    }
  }

  # test if one of the parameters to be estimated reached lower or upper bound
  attempt <- 0
  names(res$par) <- Parms2Estimate
  reinit <- FALSE
  for (parm in Parms2Estimate){
    if (isTRUE(all.equal(res$par[[parm]], lb[[parm]]))){
      reinit <- TRUE
      x0[parm] <- 0.5*(x0[parm]+lb[[parm]])
    }
    if (isTRUE(all.equal(res$par[[parm]], ub[[parm]]))){
      reinit <- TRUE
      x0[parm] <- 0.5*(x0[parm]+ub[[parm]])
    }
  }
  if (is.na(as.numeric(res$par[[1]]))){
    for (parm in Parms2Estimate) x0[parm] <- lb[[parm]]+runif(1)*(ub[[parm]]-lb[[parm]])
    reinit <- TRUE
  }
  # perform inversion with readjusted initial values if lower/upper band reached
  while (reinit==TRUE & attempt<2){
    attempt <- attempt+1
    if (verbose) print_msg(cause = 'ULBounds')
    res <- list('par' = NA*vector(length = length(Parms2Estimate)))
    for (i in TolRange){
      Tolerance <- 10**(i)
      if (is.na(res$par[1])){
        res <- tryCatch(
          {
            res <- fmincon(
              x0 = as.numeric(x0[Parms2Estimate]),
              fn = MeritFunction,
              gr =NULL,
              SpecPROSPECT = SpecPROSPECT,
              Refl = Refl, Tran = Tran,
              Input_PROSPECT = x0,
              Parms2Estimate = Parms2Estimate,
              method = "SQP", A = NULL, b = NULL, Aeq = NULL, beq = NULL,
              lb = as.numeric(lb), ub = as.numeric(ub), hin = NULL, heq = NULL,
              tol = Tolerance, maxfeval = 2000, maxiter = 1000)
          },
          error = function(cond) {
            if (verbose) print_msg(cause = 'AdjustTol',
                                   args = list('Tolerance' = Tolerance))
            return(list('par' = NA*vector(length = length(Parms2Estimate))))
          },
          finally = {}
        )
      }
    }
    names(res$par) <- Parms2Estimate
    reinit <- FALSE
    for (parm in Parms2Estimate){
      if (isTRUE(all.equal(res$par[[parm]], lb[[parm]]))){
        reinit <- TRUE
        x0[parm] <- 0.5*(x0[parm]+lb[[parm]])
      }
      if (isTRUE(all.equal(res$par[[parm]], ub[[parm]]))){
        reinit <- TRUE
        x0[parm] <- 0.5*(x0[parm]+ub[[parm]])
      }
    }
  }
  return(res)
}

#' Merit function for PROSPECT inversion
#'
#' @param x numeric. Vector of input variables to estimate
#' @param SpecPROSPECT list. Includes optical constants refractive index,
#' specific absorption coefficients and corresponding spectral bands
#' @param Refl  numeric. measured reflectance data
#' @param Tran  numeric. measured Transmittance data
#' @param Input_PROSPECT dataframe. set of PROSPECT input variables
#' @param Parms2Estimate  numeric. location of variables from Input_PROSPECT
#' to be estimated through inversion
#'
#' @return fc estimates of the parameters
#' @export
Merit_PROSPECT_RMSE <- function(x, SpecPROSPECT,
                                Refl, Tran,
                                Input_PROSPECT,
                                Parms2Estimate) {
  x[x < 0] <- 0
  Input_PROSPECT[Parms2Estimate] <- x
  RT <- do.call("PROSPECT", c(list(SpecPROSPECT = SpecPROSPECT), Input_PROSPECT))
  fcr <- fct <- 0
  if (!is.null(Refl)) fcr <- sqrt(sum((Refl - RT$Reflectance)**2) / length(RT$Reflectance))
  if (!is.null(Tran)) fct <- sqrt(sum((Tran - RT$Transmittance)**2) / length(RT$Transmittance))
  fc <- fcr + fct
  return(fc)
}

#' This function defines a regression model to estimate N from R only or T only
#' @param lambda  numeric. spectral bands corresponding to experimental data
#' @param SpecPROSPECT list. Includes optical constants refractive index,
#' specific absorption coefficients and corresponding spectral bands
#' @param Refl  numeric. Measured reflectance data
#' @param Tran  numeric. Measured Transmittance data
#' @param OptWL_R  list. optimal wavelengths used to estimate N from R only
#' @param OptWL_T  list. optimal wavelengths used to estimate N from T only
#'
#' @return Nprior prior estimation of N with R or T only
#' @importFrom stats lm runif
#' @export
Get_Nprior <- function(lambda, SpecPROSPECT = NULL, Refl = NULL, Tran = NULL,
                       OptWL_R = list(NIR = 800,SWIR = 1131),
                       OptWL_T = list(NIR = 753,SWIR = 1121)) {

  if (is.null(SpecPROSPECT)) SpecPROSPECT <- prospect::SpecPROSPECT_FullRange
  # if prior information based on Reflectance
  if (is.null(Tran)) {
    # if required spectral bands in the original data
    if (OptWL_R$SWIR %in% SpecPROSPECT$lambda) {
      OptWL <- OptWL_R$SWIR
      # else if close to spectral bands of original data
    } else if (!OptWL_R$SWIR %in% SpecPROSPECT$lambda & min(abs(OptWL_R$SWIR-SpecPROSPECT$lambda))<10) {
      OptWL <- SpecPROSPECT$lambda[which(abs(OptWL_R$SWIR-SpecPROSPECT$lambda)==min(abs(OptWL_R$SWIR-SpecPROSPECT$lambda)))]
      # else if NIR band available
    } else if (!OptWL_R$SWIR %in% SpecPROSPECT$lambda & OptWL_R$NIR %in% SpecPROSPECT$lambda) {
      print_msg(cause = 'Nprior_noRefl')
      OptWL <- OptWL_R$NIR
      # else if close to NIR band available
    } else if (!OptWL_R$NIR %in% SpecPROSPECT$lambda & min(abs(OptWL_R$NIR-SpecPROSPECT$lambda))<10) {
      print_msg(cause = 'Nprior_noRefl')
      OptWL <- SpecPROSPECT$lambda[which(abs(OptWL_R$NIR-SpecPROSPECT$lambda)==min(abs(OptWL_R$NIR-SpecPROSPECT$lambda)))]
      # else if NIR band available
    } else if (!OptWL_R$SWIR %in% SpecPROSPECT$lambda & !OptWL_R$NIR %in% SpecPROSPECT$lambda) {
      print_msg(cause = 'Nprior_noRefl_atAll')
    }
  } else if (is.null(Refl)) {
    if (OptWL_T$SWIR %in% SpecPROSPECT$lambda) {
      OptWL <- OptWL_T$SWIR
    } else if (!OptWL_T$SWIR %in% SpecPROSPECT$lambda & min(abs(OptWL_T$SWIR-SpecPROSPECT$lambda))<10) {
      OptWL <- SpecPROSPECT$lambda[which(abs(OptWL_T$SWIR-SpecPROSPECT$lambda)==min(abs(OptWL_T$SWIR-SpecPROSPECT$lambda)))]
    } else if (!OptWL_T$SWIR %in% SpecPROSPECT$lambda & OptWL_T$NIR %in% SpecPROSPECT$lambda) {
      print_msg(cause = 'Nprior_noTran')
      OptWL <- OptWL_T$NIR
    } else if (!OptWL_T$NIR %in% SpecPROSPECT$lambda & min(abs(OptWL_T$NIR-SpecPROSPECT$lambda))<10) {
      print_msg(cause = 'Nprior_noTran')
      OptWL <- SpecPROSPECT$lambda[which(abs(OptWL_T$NIR-SpecPROSPECT$lambda)==min(abs(OptWL_T$NIR-SpecPROSPECT$lambda)))]
    } else if (!OptWL_T$SWIR %in% SpecPROSPECT$lambda & !OptWL_T$NIR %in% SpecPROSPECT$lambda) {
      print_msg(cause = 'Nprior_noTran_atAll')
    }
  }

  # get the subdomain corresponding to OptWL
  SubData <- FitSpectralData(SpecPROSPECT = SpecPROSPECT, lambda = lambda,
                             Refl = Refl, Tran = Tran,
                             UserDomain = c(OptWL, OptWL), UL_Bounds = TRUE)

  # create a LUT using only the spectral band of interest
  nbSim <- 1000
  Input_PROSPECT <- data.frame('CHL' = 0.5 + 100 * runif(nbSim),
                               'CAR' = 0.5 + 20 * runif(nbSim),
                               'EWT' = 0.001 + 0.02 * runif(nbSim),
                               'LMA' = 0.001 + 0.01 * runif(nbSim),
                               'N' = 1 + 1.5 * runif(nbSim))
  LUT <- PROSPECT_LUT(Input_PROSPECT = Input_PROSPECT,
                      SpecPROSPECT = SubData$SpecPROSPECT)
  # fit a linear model between
  if (is.null(Tran)) {
    Ratio <- data.frame(c(unlist(LUT$Reflectance))) / (1 - data.frame(c(unlist(LUT$Reflectance))))
    Ratio_Meas <- data.frame(c(unlist(SubData$Refl))) / (1 - data.frame(c(unlist(SubData$Refl))))
  } else if (is.null(Refl)) {
    Ratio <- (1 - data.frame(c(unlist(LUT$Transmittance)))) / data.frame(c(unlist(LUT$Transmittance)))
    Ratio_Meas <- (1 - data.frame(c(unlist(SubData$Tran)))) / data.frame(c(unlist(SubData$Tran)))
  }
  N_Model <- lm(c(LUT$Input_PROSPECT$N) ~ c(Ratio[[1]]))
  NpriorMOD <- N_Model$coefficients[2] * Ratio + N_Model$coefficients[1]
  Nprior <- N_Model$coefficients[2] * Ratio_Meas + N_Model$coefficients[1]
  rownames(Nprior) <- NULL
  colnames(Nprior) <- 'N'
  return(Nprior)
}

#' This function uses optimal configuration identified by Spafford et al. (2020) to estimate leaf chemistry
#' prior information on N is provided if only Reflectance or only Transmittance is available
#' @param lambda  numeric. spectral bands corresponding to experimental data
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param Refl  numeric. Measured reflectance data
#' @param Tran  numeric. Measured Transmittance data
#' @param PROSPECT_version  character. Version of prospect model used for the inversion: '5', '5B', 'D', 'DB', 'PRO', 'PROB',
#' See details.
#' @param Parms2Estimate  character vector. Parameters to estimate (can be 'ALL')
#' @param InitValues  data.frame. Default values of PROSPECT parameters. During optimization,
#' @param xlub data.frame. Boundaries of the parameters to estimate.
#' The data.frame must have columns corresponding to \code{Parms2Estimate} first line being
#' the lower boundaries and second line the upper boundaries.
#' @param verbose  boolean. set to TRUE to display sample number to be inverted
#' they are used either as initialization values for parameters to estimate,
#' or as fix values for other parameters.
#' Parameters not compatible with PROSPECT_version are not taken into account.
#' @param progressBar boolean. show progressbar?
#'
#' @return Nprior vector corresponding to teh prior estimation of N based on R only or T only
#' @importFrom stats lm runif
#' @importFrom progress progress_bar
#' @export
Invert_PROSPECT_OPT <- function(lambda, SpecPROSPECT = NULL, Refl = NULL, Tran = NULL,
                                PROSPECT_version = 'D', Parms2Estimate = 'ALL',
                                InitValues = data.frame(
                                  CHL = 40, CAR = 10, ANT = 0.1, BROWN = 0.01, EWT = 0.01,
                                  LMA = 0.01, PROT = 0.001, CBC = 0.009, N = 1.5, alpha = 40),
                                xlub = data.frame(
                                  CHL = c(1e-4, 150), CAR = c(1e-4, 25), ANT = c(0, 50),
                                  BROWN = c(0, 1), EWT = c(1e-8, 0.1), LMA = c(1e-6, .06),
                                  PROT = c(1e-7, .006), CBC = c(1e-6, .054), N = c(.5, 4),
                                  alpha = c(10, 90)),
                                verbose = FALSE,
                                progressBar = TRUE) {

  if (is.null(SpecPROSPECT)) SpecPROSPECT <- prospect::SpecPROSPECT_FullRange
  # versions implemented in the package
  if (!PROSPECT_version%in%c('D','PRO')) print_msg(cause = 'WrongVersion')
  # if brown pigments required
  if ('BROWN'%in%Parms2Estimate) {
    print_msg(cause = 'NoBrown_OPT')
    Parms2Estimate <- Parms2Estimate[-match('BROWN',Parms2Estimate)]
  }
  # define optimal domain for the different constituents
  OptDomain <-OptDomain_R <- OptDomain_T <- prospect::OptDomain_RT
  # check if list of parameters applicable to PROSPECT version
  parms_checked <- check_prospect_parms(PROSPECT_version = PROSPECT_version,
                                        Parms2Estimate = Parms2Estimate,
                                        Est_Brown_Pigments = FALSE,
                                        Est_alpha = FALSE,
                                        xlub = xlub,
                                        InitValues = InitValues)
  Parms2Estimate <- parms_checked$Parms2Estimate
  ANTinit <- parms_checked$InitValues$ANT
  if (PROSPECT_version == 'PRO' & 'LMA'%in%Parms2Estimate){
    print_msg(cause = 'version_PROSPECT_Invert_PROSPECT_OPT')
  }
  # adjust leaf optics and convert to data frame if needed
  AdjustedLOP <- FitSpectralData(SpecPROSPECT = SpecPROSPECT,
                                 lambda = lambda,
                                 Refl = Refl, Tran = Tran)
  # prior assessment of N if R or T missing
  Nprior <- SetNValues(Refl = AdjustedLOP$Refl,
                       Tran = AdjustedLOP$Tran,
                       SpecPROSPECT = AdjustedLOP$SpecPROSPECT)
  ParmEst <- list()
  if (is.null(AdjustedLOP$Refl) | is.null(AdjustedLOP$Tran)) ParmEst$N <- Nprior
  if ('N'%in%Parms2Estimate) Parms2Estimate <- Parms2Estimate[-which(Parms2Estimate=='N')]
  # initialize parameters for each sample (useful when using prior N)
  List_Init <- SetInitParm(Parms2Estimate = Parms2Estimate,
                           ParmEst = ParmEst,
                           PROSPECT_version = PROSPECT_version,
                           Nprior = Nprior,
                           ANTinit = ANTinit,
                           OptDomain = OptDomain,
                           Refl = AdjustedLOP$Refl,
                           Tran = AdjustedLOP$Tran)
  ParmEst <- List_Init$ParmEst

  for (parm in List_Init$Parms2Estimate){
    if (length(ParmEst[[parm]])==0){
      # Fit spectral data to match PROSPECT with user optical properties
      OptLOP <- FitSpectralData(SpecPROSPECT = AdjustedLOP$SpecPROSPECT,
                                lambda = AdjustedLOP$lambda,
                                Refl = AdjustedLOP$Refl,
                                Tran = AdjustedLOP$Tran,
                                UserDomain = OptDomain[[parm]],
                                UL_Bounds = List_Init$UL_Bounds[[parm]])
      ParmDisplay <- parm
      if (parm=='EWT' & 'LMA'%in%Parms2Estimate) ParmDisplay <- 'EWT & LMA'
      if (progressBar==TRUE){
        pb <- progress::progress_bar$new(
          format = paste("Estimation of ",ParmDisplay," [:bar] :percent in :elapsedfull , estimated time remaining :eta"),
          total = AdjustedLOP$nbSamples, clear = FALSE, width= 100)
      }
      for (i in seq_len(AdjustedLOP$nbSamples)){
        res <- Invert_PROSPECT(SpecPROSPECT = OptLOP$SpecPROSPECT,
                               Refl = OptLOP$Refl[[i]],
                               Tran = OptLOP$Tran[[i]],
                               PROSPECT_version = List_Init$PROSPECT_version[[parm]],
                               Parms2Estimate = List_Init$Parms2Estimate_tmp[[parm]],
                               InitValues = List_Init$InitValues[[parm]][i,],
                               xlub = xlub,
                               verbose = verbose)
        ParmEst[[parm]][i] <- res[[parm]]
        if (parm=='EWT' & !is.na(match('LMA',Parms2Estimate))) ParmEst$LMA[i] <- res$LMA
        if (progressBar==TRUE) pb$tick()
      }
    }
  }
  ParmEst <- data.frame(ParmEst)
  return(ParmEst)
}


#' Function performing optimal feature selection based on sequential forward feature slection
#'
#' @param Refl numeric. matrix of reflectances (n spectral bands x p samples)
#' @param Tran numeric. matrix of transmittances (n spectral bands x p samples)
#' @param lambda numeric. spectral bands corresponding to reflectance and transmittance measurements
#' @param BiochTruth numeric. value of biophysical/biochemical parameter to estimate for each sample
#' @param Target character. name of the parameter.
#' Should be picked between "CHL", "CAR", "ANT", "BROWN", "EWT", "PROT", "CBC", "N"
#' @param Parms2Estimate  character vector. Parameters to estimate (can be 'ALL')
#' @param InitValues  data.frame. Default values of PROSPECT parameters. During optimization,
#' they are used either as initialization values for parameters to estimate,
#' or as fix values for other parameters.
#' Parameters not compatible with PROSPECT_version are not taken into account.
#' @param xlub data.frame. Boundaries of the parameters to estimate.
#' The data.frame must have columns corresponding to \code{Parms2Estimate} first line being
#' the lower boundaries and second line the upper boundaries.
#' @param SpecPROSPECT list. Optical constants for a spectral domain broader than spectral_domain
#' @param spectral_domain vector. defines minimum and maximum wavelengths of the spectral domain (in nm).
#' Assumes 1nm spectral sampling
#' @param spectral_width vector. width of individual spectral features (in nm)
#' @param number_features vector. number of features to be identified
#' @param PROSPECT_version  character. Version of prospect model used for the inversion: '5', '5B', 'D', 'DB', 'PRO', 'PROB',
#' @param MeritFunction  character. name of the function to be used as merit function
#' with given criterion to minimize (default = RMSE)
#' @param Est_alpha boolean. should alpha be estimated or not?
#' @param verbose boolean. set to TRUE to display sample number to be inverted
#' @param nbCPU numeric. defines number of CPU for multithread processing
#'
#' @param Continue boolean. set to TRUE if the function has already been run for a lower value of number_features
#' @param AlreadyDone list. contains output of function optimal_features_SFS previouly run with lower value of number_features
#' @return list containing :
#' - SpectralFeatures = the Spectral Features identified as most relevant for estimation of Target (vector)
#' - SpectralFeatures_List = the Spectral Features identified as most relevant for estimation of Target (list)
#' - EstimatedParm = the estimated Target parameter for each additional spectral feature
#' - RMSE = the minimum RMSE corresponding to the estimation of Target for each asditional spectral feature
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @export

optimal_features_SFS <- function(Refl = NULL, Tran = NULL, lambda, BiochTruth,
                                 Target, Parms2Estimate = "ALL",
                                 InitValues = data.frame(
                                   CHL = 40, CAR = 10, ANT = 0.1, BROWN = 0.01,
                                   EWT = 0.01, LMA = 0.01, PROT = 0.001,
                                   CBC = 0.009, N = 1.5, alpha = 40),
                                 xlub = data.frame(
                                   CHL = c(1e-4, 150), CAR = c(1e-4, 25),
                                   ANT = c(0, 50), BROWN = c(0, 1),
                                   EWT = c(1e-8, 0.1), LMA = c(1e-6, .06),
                                   PROT = c(1e-7, .006), CBC = c(1e-6, .054),
                                   N = c(.5, 4), alpha = c(10, 90)),
                                 SpecPROSPECT, spectral_domain, spectral_width,
                                 number_features,PROSPECT_version = 'D',
                                 MeritFunction = "Merit_PROSPECT_RMSE",
                                 Est_alpha = FALSE, verbose = FALSE,
                                 nbCPU = 1,Continue = FALSE, AlreadyDone = NULL){

  # split spectral_domain into subdomains based on spectral_width
  FullDomain <- seq(spectral_domain[1],spectral_domain[2])
  x <- seq_along(FullDomain)
  subdomains <- split(FullDomain,ceiling(x/spectral_width))
  nb_subdomains <- length(subdomains)
  # spectral bands identified as optimal
  if (Continue){
    # number of features to start from
    featStart <- length(AlreadyDone$RMSE)+1
    Estimated_All <- AlreadyDone$EstimatedParm
    Initial_Features_list <- AlreadyDone$SpectralFeatures_List
    Initial_Features <- AlreadyDone$SpectralFeatures
    minRMSE <- AlreadyDone$RMSE

    # number of features to start from
    for (i in seq_len(length(AlreadyDone$RMSE))){
      for (j in rev(seq_len(length.out = length(subdomains)))){
        if (AlreadyDone$SpectralFeatures_List[[i]][1] %in% subdomains[[j]]){
          subdomains[[j]] <- NULL
        }
      }
    }
  } else {
    featStart <- 1
    # for incremental number of spectral features from 1 to number_features
    Estimated_All <- list()
    Initial_Features_list <- list()
    Initial_Features <- minRMSE <- c()
  }

  if (number_features>=featStart){
    for (feat in featStart:number_features){
      message(paste('Identify Feature #',feat))
      Perf <- c()
      plan(multisession, workers = nbCPU) ## Parallelize using four cores
      schedule <- ceiling(length(subdomains)/nbCPU)
      Estimated_Parm <- future_lapply(subdomains,FUN = Invert_PROSPECT_subdomain,
                                      Refl = Refl, Tran = Tran,
                                      SpecPROSPECT = SpecPROSPECT,
                                      lambda = lambda, BiochTruth = BiochTruth,
                                      Target = Target,
                                      Parms2Estimate  = Parms2Estimate,
                                      Initial_Features = Initial_Features,
                                      PROSPECT_version = PROSPECT_version,
                                      MeritFunction = MeritFunction,
                                      Est_alpha = Est_alpha,
                                      InitValues = InitValues,
                                      xlub = xlub, verbose = FALSE,
                                      future.scheduling = schedule)
      plan(sequential)
      Perf <- c()
      for (i in seq_len(length(Estimated_Parm))){
        Perf <- c(Perf,Estimated_Parm[[i]]$Perf)
      }
      BestPerf <- which(Perf==min(Perf,na.rm = TRUE))
      Estimated_All[[feat]] <- Estimated_Parm[[BestPerf[1]]]$Estimated
      minRMSE <- c(minRMSE,min(Perf,na.rm = TRUE))
      print(min(Perf,na.rm = TRUE))
      message('Spectral domain identified:')
      print(subdomains[[BestPerf[1]]])
      Initial_Features <- c(Initial_Features,subdomains[[BestPerf[1]]])
      Initial_Features_list[[feat]] <- subdomains[[BestPerf[1]]]
      subdomains[[BestPerf[1]]] <- NULL
    }
  }
  return(list('SpectralFeatures' = Initial_Features,
              'SpectralFeatures_List' = Initial_Features_list,
              'EstimatedParm' = Estimated_All,
              'RMSE' = minRMSE))
}


#' Function performing PROSPECT inversion on a spectral subset combining an
#' initial feature subset (Initial_Features) and an additional feature subset (New_Features)
#'
#' @param New_Features spectral bands (in nm) to be added to the Initial_Features
#' @param Refl numeric. matrix of reflectances (n spectral bands x p samples)
#' @param Tran numeric. matrix of transmittances (n spectral bands x p samples)
#' @param SpecPROSPECT list. Optical constants for a spectral domain broader than spectral_domain
#' @param lambda numeric. spectral bands corresponding to reflectance and transmittance measurements
#' @param BiochTruth numeric. value of biophysical/biochemical parameter to estimate for each sample
#' @param Target character. name of the parameter.
#' Should be picked between "CHL", "CAR", "ANT", "BROWN", "EWT", "PROT", "CBC", "N"
#' @param Parms2Estimate  character vector. Parameters to estimate (can be 'ALL')
#' @param Initial_Features Initial feature set (in nm)
#' @param InitValues  data.frame. Default values of PROSPECT parameters. During optimization,
#' they are used either as initialization values for parameters to estimate,
#' or as fix values for other parameters.
#' Parameters not compatible with PROSPECT_version are not taken into account.
#' @param PROSPECT_version  character. Version of prospect model used for the inversion: 'D' or 'PRO'
#' @param Est_Brown_Pigments boolean. should brown pigments be accounted for during inversion?
#' @param Est_alpha boolean. should alpha be accounted for during inversion?
#' @param MeritFunction  character. name of the function to be used as merit function
#' with given criterion to minimize (default = RMSE)
#' @param xlub data.frame. Boundaries of the parameters to estimate.
#' The data.frame must have columns corresponding to \code{Parms2Estimate} first line being
#' the lower boundaries and second line the upper boundaries.
#' @param verbose boolean. set to TRUE to display sample number to be inverted
#' @return list containing estimated value of the target and RMSE compared to measured value
#' @importFrom pracma rmserr
#' @export

Invert_PROSPECT_subdomain <- function(New_Features, Refl, Tran, SpecPROSPECT,
                                      lambda, BiochTruth,
                                      Target, Parms2Estimate, Initial_Features,
                                      InitValues = data.frame(
                                        CHL = 40, CAR = 10, ANT = 0.1,
                                        BROWN = 0.01, EWT = 0.01, LMA = 0.01,
                                        PROT = 0.001, CBC = 0.009, N = 1.5,
                                        alpha = 40),
                                      PROSPECT_version = "D",
                                      MeritFunction = "Merit_PROSPECT_RMSE",
                                      Est_Brown_Pigments = FALSE,
                                      Est_alpha = FALSE,
                                      xlub = data.frame(
                                        CHL = c(1e-4, 150), CAR = c(1e-4, 25),
                                        ANT = c(0, 50), BROWN = c(0, 1),
                                        EWT = c(1e-8, 0.1), LMA = c(1e-6, 0.06),
                                        PROT = c(1e-7, .006), CBC = c(1e-6, 0.054),
                                        N = c(.5, 4), alpha = c(10, 90)),
                                      verbose = FALSE){

  SpectralDomain <- c(Initial_Features,New_Features)
  # Fit spectral data to match PROSPECT with user optical properties
  SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT, lambda=lambda,
                             Refl=Refl, Tran=Tran, UserDomain = SpectralDomain)
  # Invert PROSPECT with optimal spectral information
  Est <- c()
  for (i in seq_len(length.out = SubData$nbSamples)){
    if (verbose) print(i)
    res <- Invert_PROSPECT(SpecPROSPECT = SubData$SpecPROSPECT,
                           Refl = SubData$Refl[,i],Tran = SubData$Tran[,i],
                           PROSPECT_version = PROSPECT_version,
                           Parms2Estimate = Parms2Estimate,
                           InitValues = InitValues,
                           MeritFunction = MeritFunction,
                           xlub = xlub,
                           Est_Brown_Pigments = Est_Brown_Pigments,
                           Est_alpha = Est_alpha)
    # only save results for variable of interest
    Est <- c(Est,res[[Target]])
  }
  Perf <- pracma::rmserr(BiochTruth,Est,summary = FALSE)$rmse
  rm(Refl,Tran,SubData)
  gc()
  return(list("Estimated" =  Est,"Perf" =  Perf))
}

#' Function printing message
#' - if leaf optics do not match with spectral domain
#' - if leaf optics do not match with expected data class
#' - if confusion is identified for prospect version to use
#' - if missing input variables during PROSPECT_LUT simulations
#' @param cause character. cause raising message print
#' @param args list. optional arguments used to print message
#'
#' @return none
#' @export

print_msg <- function(cause, args = NULL){
  if (cause == 'lengthLOP'){
    message('MISMATCH between leaf optical properties and spectral domain')
    message('Make sure leaf optical properties match with the spectral domain')
    message('defined in SpecPROSPECT$lambda')
    message('if not, use the function "FitSpectralData"')
    stop()
  } else if (cause == 'classLOP'){
    message('the class of leaf optical properties is not identified')
    message('Make sure leaf optical properties are either data.frame, or matrix, or numeric vector')
    stop()
  } else if (cause == 'version_PROSPECT'){
    message("PROT and/or CBC are not set to 0")
    message("LMA is not set to 0 neither, which is physically incorrect")
    message("(LMA = PROT + CBC)")
    message("We assume that PROSPECT-PRO was called and set LMA to 0")
    message("Please correct input parameters LMA, PROT and/or CBC if needed")
  } else if (cause == 'version_PROSPECT_Invert_PROSPECT_OPT'){
    message('LMA is not estimated directly using PROSPECT-PRO')
    message('Please run inversion for PROSPECT-PRO and sum PROT and CBC')
    message('if you want to get estimated LMA from PROSPECT-PRO')
  } else if (cause == 'Missing_Input'){
    message('_________________ WARNING _________________')
    message('  The list of parameters provided for LUT  ')
    message('      in the Input_PROSPECT variable       ')
    message('does not include the full set of parameters')
    message('')
    message(' The following parameters will be set to  ')
    message('   their default value as given here      ')
    WhichMissing <- args$Expected[which(!names(args$Expected) %in% args$Input)]
    print(WhichMissing)
  } else if (cause == 'diff_RTshape'){
    message('Reflectance and Transmittance do not have the same number of samples')
    message('please fix that')
    stop()
  } else if (cause == 'WrongVersion'){
    stop('PROSPECT_version not available. Choice is limited to "5", "5B", "D", "DB", "PRO", "PROB".')
  } else if (cause == 'NoBrown_OPT'){
    message('brown pigments are not accounted for when performing optimal estimation of leaf chemistry')
    message('model version excluding brown pigments will be used instead')
  } else if (cause == 'AdjustTol'){
    message('Adjusting Tolerance value for iterative optimization:  ',args$Tolerance)
  } else if (cause == 'ULBounds'){
    message('lower or upper bound reached for one or several parameters to estimate')
    message('re-run inversion with updated initial values to attempt proper convergence')
  } else if (cause == 'Nprior_noRefl'){
    warning("________________________ WARNING _______________________")
    warning("The optimal prior estimation of N using Reflectance only")
    warning("requires information at 1131nm.")
    warning("The reflectance does not include this spectral band")
    warning("Using reflectance at 800 nm instead")
  } else if (cause == 'Nprior_noRefl_atAll'){
    warning("________________________ WARNING _______________________")
    warning("The spectral information of the reflectance provided here")
    warning("does not contain the spectral bands required to estimate")
    warning("prior information about N.")
    warning("The process will stop")
    stop()
  } else if (cause == 'Nprior_noTran'){
    warning("________________________ WARNING _______________________")
    warning("The optimal prior estimation of N using Transmittance only")
    warning("requires information at 1121nm.")
    warning("The Transmittance does not include this spectral band")
    warning("Using Transmittance at 753 nm instead")
  } else if (cause == 'Nprior_noTran_atAll'){
    warning("________________________ WARNING _______________________")
    warning("The spectral information of the Transmittance provided here")
    warning("does not contain the spectral bands required to estimate")
    warning("prior information about N.")
    warning("The proces will stop")
    stop()
  } else if (cause == 'NoOpt_ANT'){
    message('Currently no optimal estimation for anthocyanins')
    message('PROSPECT inversion will be performed using full spectral information')
  } else if (cause == 'missing_InitValues'){
    stop('missing prospect parameters in "InitValues"')
  } else if (cause == 'missing_xlub'){
    stop('missing prospect parameters in "xlub"')
  } else if (cause == 'unknown_parm'){
    message('_________________ WARNING _________________')
    message(' unknown parameter requested for inversion ')
    print(args$Parms2Estimate[which(!args$Parms2Estimate %in% args$allParms)])
    message('will be ignored')
  }
  return(invisible())
}


#' Function to reshape reflectance and transmittance, and convert into dataframes
#'
#' @param Refl  numeric. measured reflectance data
#' @param Tran  numeric. measured Transmittance data
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#'
#' @return RT list of leaf optics converted into dataframe and coresponding number of samples
#' @export

reshape_lop4inversion <- function(Refl, Tran, SpecPROSPECT){
  # check if Refl/Tran is data frame, matrix or vector
  # convert into data frame if not
  RT <- list('Refl' = Refl, 'Tran' = Tran)
  for (lop in names(RT)){
    if (!is.null(RT[[lop]])){
      classLOP <- class(RT[[lop]])[1]
      if (!classLOP %in% c('numeric', 'matrix', 'data.frame', 'data.table')) print_msg(cause = 'classLOP')
      if (classLOP %in% c('numeric', 'matrix')) RT[[lop]] <- data.frame(RT[[lop]])
      nbwl <- nrow(RT[[lop]])
      # check if the number of spectral bands is compatible with SpecPROSPECT
      if (nbwl==length(SpecPROSPECT$lambda)) rownames(RT[[lop]]) <- SpecPROSPECT$lambda
      if (!nbwl==length(SpecPROSPECT$lambda)) print_msg(cause = 'lengthLOP')
    }
  }
  # check if Refl & Tran have the same number of samples
  if (!is.null(RT$Refl) & !is.null(RT$Tran)){
    if (!ncol(RT$Refl)==ncol(RT$Tran)) print_msg(cause = 'diff_RTshape')
  }
  RT$nbSamples <- unique(c(ncol(RT$Refl),ncol(RT$Tran)))
  return(RT)
}


#' Function to check good agreement between prospect version and parameter list
#'
#' @param PROSPECT_version character. version used for inversion: 'D' or 'PRO'
#' @param Est_Brown_Pigments boolean. should brown pigments be accounted for during inversion?
#' @param Est_alpha boolean. should alpha be accounted for during inversion?
#' @param Parms2Estimate  character vector. Parameters to estimate (can be 'ALL')
#' @param xlub data.frame. Boundaries of the parameters to estimate.
#' The data.frame must have columns corresponding to \code{Parms2Estimate} first
#' line being the lower boundaries and second line the upper boundaries.
#' @param InitValues  dataframe. Default values of PROSPECT parameters. During
#' optimization, they are used either as initialization values for parameters to
#' estimate, or as fix values for other parameters.
#' Parameters not compatible with PROSPECT_version are not taken into account.
#'
#' @return list of parameters to estimate & corresponding lower/upper boundaries
#' @export

check_prospect_parms <- function(PROSPECT_version,
                                 Parms2Estimate,
                                 Est_Brown_Pigments,
                                 Est_alpha,
                                 xlub,
                                 InitValues){

  # check if version required is available
  if (!PROSPECT_version %in% c('D', 'PRO')) print_msg(cause = 'WrongVersion')
  # define parameters if PROSPECT-D or PROSPECT-PRO
  if (PROSPECT_version == 'D') allParms <- c('N', 'CHL', 'CAR', 'ANT', 'EWT', 'LMA')
  if (PROSPECT_version == 'PRO') allParms <- c('N', 'CHL', 'CAR', 'ANT', 'EWT', 'PROT', 'CBC')
  # add brown pigments and alpha angle if required
  if (Est_Brown_Pigments==TRUE) allParms <- c(allParms, "BROWN")
  if (Est_alpha==TRUE) allParms <- c(allParms, "alpha")

  # add default values to xlub in case they were not defined
  xlub_default <- data.frame(CHL = c(1e-4, 150), CAR = c(1e-4, 25),
                             ANT = c(0, 50), BROWN = c(0, 1),
                             EWT = c(1e-8, 0.1), LMA = c(1e-6, .06),
                             PROT = c(1e-7, .006), CBC = c(1e-6, .054),
                             N = c(.5, 4), alpha = c(10, 90))
  AddedParm_LB <- setdiff(names(xlub_default),names(xlub))
  for (ad in AddedParm_LB) xlub[[ad]] <- xlub_default[[ad]]

  InitValues_default = data.frame(CHL = 40, CAR = 10,
                                  ANT = 0.1, BROWN = 0.0,
                                  EWT = 0.01, LMA = 0.01,
                                  PROT = 0.001, CBC = 0.009,
                                  N = 1.5, alpha = 40)
  AddedParm_init <- setdiff(names(InitValues_default),names(InitValues))
  for (ad in AddedParm_init) InitValues[[ad]] <- InitValues_default[[ad]]

  # if 'ALL' is provided, then assess all parameters available
  if ("ALL" %in% Parms2Estimate) Parms2Estimate <- allParms
  # if unknown parameter is provided, then warn
  if (!all(Parms2Estimate %in% allParms)) print_msg(cause = 'unknown_parm',
                                                    args = list('Parms2Estimate' =Parms2Estimate,
                                                                'allParms'= allParms))
  Parms2Estimate <- allParms[allParms %in% Parms2Estimate]
  # if missing value in xlub defining lower and upper boudaries
  if (!all(allParms %in% names(xlub))) print_msg(cause = 'missing_xlub')
  # if missing value in InitValues defining initial values for inversion
  if (!all(allParms %in% names(InitValues))) print_msg(cause = 'missing_InitValues')
  InitValues <- InitValues[allParms[allParms %in% names(InitValues)]]
  if (PROSPECT_version == "PRO"){
    InitValues$LMA <- 0
    if (is.null(InitValues$PROT)) InitValues$PROT <- 0.001
    if (is.null(InitValues$CBC)) InitValues$CBC <- 0.009
  }
  if (PROSPECT_version == "D"){
    InitValues$PROT <- InitValues$CBC <- 0
    if (is.null(InitValues$LMA)) InitValues$LMA <- 0.01
  }
  if (Est_Brown_Pigments==FALSE) InitValues$BROWN <- 0
  xlub <- data.frame(xlub[Parms2Estimate])
  lb <- xlub[1,]
  ub <- xlub[2,]
  return(list('lb' = lb, 'ub' = ub, 'Parms2Estimate' = Parms2Estimate,
              'InitValues' = InitValues))
}

#' Function set N values if reflectance or transmittance is missing
#'
#' @param Refl dataframe. reflectance data
#' @param Tran dataframe. transmittance data
#' @param SpecPROSPECT dataframe. spectral properties used in PROSPECT
#'
#' @return Nprior dataframe. N values corresponding to Refl or Tran
#' @export

SetNValues <- function(Refl, Tran, SpecPROSPECT){

  if (is.null(Refl) | is.null(Tran)){
    # compute prior estimate of N
    message('computing prior estimation of N as both R & T are not provided')
    Nprior <- Get_Nprior(SpecPROSPECT = SpecPROSPECT,
                         lambda = SpecPROSPECT$lambda,
                         Refl = Refl, Tran = Tran)
  } else {
    Nprior <- data.frame('N' = rep(1.5,ncol(Refl)))
  }
  return(Nprior)
}

#' Function to set initial parameterization when performing PROSPECT inversion
#' using optimal configuration
#'
#' @param Parms2Estimate  character. Parameters to estimate
#' @param ParmEst character. PROSPECT parameters to be estimated
#' @param PROSPECT_version  character. Version of prospect model used for the
#' inversion: '5', '5B', 'D', 'DB', 'PRO', 'PROB',
#' @param Nprior numeric prior estimation of N
#' @param ANTinit numeric prior estimation of ANT
#' @param OptDomain vector. optimal spectral domain
#' @param Refl dataframe reflectance data
#' @param Tran dataframe transmittance data
#'
#' @return list containing Parms2Estimate, ParmEst, UL_Bounds, PROSPECT_version.
#' Parms2Estimate_tmp, InitValues
#' @export

SetInitParm <- function(Parms2Estimate, ParmEst, PROSPECT_version, Nprior,
                        ANTinit, OptDomain, Refl, Tran){

  UL_Bounds <- PROSPECT_version_tmp <-
    Parms2Estimate_tmp <- InitValues <- list()
  for (parm in Parms2Estimate){
    ParmEst[[parm]] <- c()
    InitValues[[parm]] <- list()
    if (length(OptDomain[[parm]])==2){
      UL_Bounds[[parm]] <- TRUE
    } else {
      UL_Bounds[[parm]] <- FALSE
    }
    if (parm == "ANT"){
      if (!is.na(match(PROSPECT_version,'5'))){
        message ('Cannot estimate anthocyanins using PROSPECT-5')
        Parms2Estimate <- Parms2Estimate[-which(Parms2Estimate=='ANT')]
      } else {
        print_msg(cause = 'NoOpt_ANT')
        PROSPECT_version_tmp[[parm]] <- PROSPECT_version
        Parms2Estimate_tmp[[parm]] <- c('CHL', 'CAR', 'ANT')
      }
      InitValues[[parm]] <- data.frame(CHL=40, CAR=10, ANT=ANTinit,
                                       BROWN=0, EWT=0.01, LMA=0.01, N=Nprior)
    }
    if (parm == "CHL" | parm == "CAR"){
      PROSPECT_version_tmp[[parm]] <- PROSPECT_version
      if (!PROSPECT_version =='5'){
        Parms2Estimate_tmp[[parm]] <- c('CHL', 'CAR', 'ANT')
      } else if (PROSPECT_version =='5'){
        Parms2Estimate_tmp[[parm]] <- c('CHL', 'CAR')
      }
      InitValues[[parm]] <- data.frame(CHL=40, CAR=10, ANT=ANTinit, BROWN=0,
                                       EWT=0.01, LMA=0.01, N=Nprior)
    }
    if (parm == "EWT"){
      PROSPECT_version_tmp[[parm]] <- PROSPECT_version
      if (!PROSPECT_version =='PRO'){
        Parms2Estimate_tmp[[parm]] <- c('EWT', 'LMA')
        InitValues[[parm]] <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01,
                                         LMA=0.01, N=Nprior)
      } else if (PROSPECT_version =='PRO'){
        Parms2Estimate_tmp[[parm]] <- c('EWT', 'PROT', 'CBC')
        InitValues[[parm]] <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01,
                                         LMA=0.00, PROT=0.001, CBC=0.009,
                                         N=Nprior)
      }
    }
    if (parm == "LMA" & !PROSPECT_version == 'PRO'){
      PROSPECT_version_tmp[[parm]] <- PROSPECT_version
      Parms2Estimate_tmp[[parm]] <- c('EWT', 'LMA')
      InitValues[[parm]] <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01,
                                       LMA=0.01, N=Nprior)
    }

    if (parm == "PROT" | parm == "CBC"){
      PROSPECT_version_tmp[[parm]] <- 'PRO'
      Parms2Estimate_tmp[[parm]] <- c('EWT', 'PROT', 'CBC')
      InitValues[[parm]] <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01,
                                       LMA=0.00, PROT=0.001, CBC=0.009,
                                       N=Nprior)
    }
    if (!is.null(Refl) & !is.null(Tran)){
      Parms2Estimate_tmp[[parm]] <- c(Parms2Estimate_tmp[[parm]],'N')
    }
  }
  return(list('Parms2Estimate' = Parms2Estimate,
              'ParmEst' = ParmEst,
              'UL_Bounds' = UL_Bounds,
              'PROSPECT_version' = PROSPECT_version_tmp,
              'Parms2Estimate_tmp' = Parms2Estimate_tmp,
              'InitValues'=InitValues))
}


#' #' Function plotting results of PROSPECT inversion for a given biophysical property
#' #'
#' #' @param BP_df dataframe. should include fields 'estimated', 'measured', and 'config'
#' #' @param Labs character. labels for X and Y axes
#' #' @param MinMax numeric. min and max axis values
#' #' @param Colors character. colors corresponding to the different configurations in BP_df
#' #' @param stats boolean. should statistics be displayed on figure
#' #' @param filename character. path for the file to be saved
#' #' @return plotxy ggplot object
#' #' @import ggplot2
#' #' @import dplyr
#' #' @importFrom Metrics rmse
#' #' @export
#' plotinv <- function(BP_df, Labs = NULL, MinMax = NULL, Colors = NULL,
#'                     stats = TRUE, filename = NULL){
#'
#'   if (is.null(BP_df$config)) BP_df$config <- NA
#'   if (is.null(Labs)) Labs <- c('Estimated', 'Measured')
#'   if (is.null(Colors)) Colors <- rainbow(unique(BP_df$config))
#'
#'
#'   measest_vals <- c(BP_df$measured, BP_df$estimated)
#'   if (is.null(MinMax)) {
#'     MinMax0 <- c(min(measest_vals, na.rm = T),
#'                  max(measest_vals, na.rm = T))
#'     MinMax <- c(min(0,MinMax0[1]-0.1*diff(MinMax0)),
#'                 MinMax0[2]+0.1*diff(MinMax0))
#'
#'   }
#'   plotxy <- ggplot2::ggplot(BP_df, aes(x = estimated,
#'                                        y = measured,
#'                                        group = config)) +
#'     geom_point(aes(pch = config, color = config, stroke = 1)) +
#'     theme(aspect.ratio = 1) +
#'     xlim(MinMax[1], MinMax[2]) + ylim(MinMax[1], MinMax[2]) +
#'     scale_color_manual(values = Colors) +
#'     scale_shape_manual(values = c(20,20)) +
#'     labs(x = Labs$x, y = Labs$y) +
#'     theme(legend.position = "bottom",
#'           legend.title = element_text(color = "white"),
#'           legend.text = element_text(size = 14),
#'           axis.text = element_text(size=15),
#'           axis.title.x = element_text(size=16, face="bold"),
#'           axis.title.y = element_text(size=16, face="bold")) +
#'     guides(fill = guide_legend(nrow = 2),size = 'none')
#'
#'   # Add 1:1 line
#'   plotxy <- plotxy + geom_abline(slope = 1, intercept = 0,linetype='dashed',size=1.25)
#'   # add statistics if needed
#'   BP_conf <- group_split(BP_df %>% group_by(config))
#'   if (stats ==TRUE){
#'     # add R2 value
#'     for (conf in 1:length(unique(BP_df$config))){
#'       chl.lm = lm(estimated ~ measured, data = BP_conf[[conf]])
#'       R2 <- summary(chl.lm)$r.squared
#'       R2.expr <- paste("bolditalic(R ^ 2) == ",format(round(R2, 2),nsmall = 2))
#'       plotxy <- plotxy + annotate("text", vjust = "bottom",
#'                                   x = MinMax[1] +0.00*diff(MinMax),
#'                                   y = MinMax[1] +(1.07-0.1*conf)*diff(MinMax),
#'                                   label = R2.expr, parse = TRUE, hjust = 0, size=5,
#'                                   color = Colors[conf])
#'     }
#'     # add RMSE value
#'     for (conf in 1:length(unique(BP_df$config))){
#'       RMSE <- Metrics::rmse(actual = BP_conf[[conf]]$measured,
#'                             predicted = BP_conf[[conf]]$estimated)
#'       RMSE.expr <- paste("bolditalic(RMSE) ==",format(round(RMSE, 2),nsmall = 2))
#'       plotxy <- plotxy + annotate("text", vjust = "bottom",
#'                                   x = MinMax[1] +0.3*diff(MinMax),
#'                                   y = MinMax[1] +(1.07-0.1*conf)*diff(MinMax),
#'                                   label = RMSE.expr, parse = TRUE, hjust = 0, size=5,
#'                                   color = Colors[conf])
#'     }
#'   }
#'   if (!is.null(filename)){
#'     ggsave(filename = filename, plot = plotxy, device = "png",
#'            scale = 1, width = 12, height = 12, units = "cm", dpi = 300)
#'   }
#'   return(plotxy)
#' }
