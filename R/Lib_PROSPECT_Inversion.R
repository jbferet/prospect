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
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param Refl  numeric. Measured reflectance data
#' @param Tran  numeric. Measured Transmittance data
#' @param Parms2Estimate  character vector. Parameters to estimate (can be 'ALL')
#' @param InitValues  data.frame. Default values of PROSPECT parameters. During optimization,
#' they are used either as initialization values for parameters to estimate,
#' or as fix values for other parameters.
#' Parameters not compatible with PROSPECT_version are not taken into account.
#' @param PROSPECT_version  character. Version of prospect model used for the inversion: '5', '5B', 'D', 'DB', 'PRO', 'PROB',
#' See details.
#' @param MeritFunction  character. name of the function to be used as merit function
#' with given criterion to minimize (default = RMSE)
#' @param xlub data.frame. Boundaries of the parameters to estimate.
#' The data.frame must have columns corresponding to \code{Parms2Estimate} first line being
#' the lower boundaries and second line the upper boundaries.
#' @param alphaEst boolean. should alpha be estimated or not?
#'
#' @return OutPROSPECT estimated values corresponding to Parms2Estimate
#' @details
#' Six versions of prospect are available for inversion.
#' The version is depending on the parameters taken into account:
#'
#' | Version  | 5                                       | 5B                                    | D                                      | DB                                     | PRO                                    | PROB
#' | :------: |:--------------------------------------:|:--------------------------------------:|:--------------------------------------:|:--------------------------------------:|:--------------------------------------:|:---------------------------------------:|
#' | CHL      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | CAR      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | ANT      |                                        |                                        |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | BROWN    |                                        |`r emojifont::emoji('white_check_mark')`|                                        |`r emojifont::emoji('white_check_mark')`|                                        |`r emojifont::emoji('white_check_mark')`|
#' | EWT      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | LMA      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|                                        |
#' | PROT     |                                        |                                        |                                        |                                        |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | CBC      |                                        |                                        |                                        |                                        |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | N        |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#'
#' Argument `InitValues` is expecting a default value for each of the parameters as well as an `alpha` value.
#' @importFrom pracma fmincon
#' @export
#' @md
Invert_PROSPECT <- function(SpecPROSPECT, Refl = NULL, Tran = NULL,
                            InitValues = data.frame(
                              CHL = 40, CAR = 10, ANT = 0.1, BROWN = 0.01, EWT = 0.01,
                              LMA = 0.01, PROT = 0.001, CBC = 0.009, N = 1.5, alpha = 40),
                            Parms2Estimate = "ALL",
                            PROSPECT_version = "D",
                            MeritFunction = "Merit_RMSE_PROSPECT",
                            xlub = data.frame(
                              CHL = c(1e-4, 150), CAR = c(1e-4, 25), ANT = c(0, 50),
                              BROWN = c(0, 1), EWT = c(1e-7, .08), LMA = c(1e-6, .04),
                              PROT = c(1e-7, .005), CBC = c(1e-6, .04), N = c(.5, 4),
                              alpha = c(10, 90)),
                            alphaEst = FALSE) {

  # define PROSPECT input parameters
  if (PROSPECT_version == "5") {
    allParms <- c("CHL", "CAR", "EWT", "LMA", "N")
  } else if (PROSPECT_version == "5B") {
    allParms <- c("CHL", "CAR", "BROWN", "EWT", "LMA", "N")
  } else if (PROSPECT_version == "D") {
    allParms <- c("CHL", "CAR", "ANT", "EWT", "LMA", "N")
  } else if (PROSPECT_version == "DB") {
    allParms <- c("CHL", "CAR", "ANT", "BROWN", "EWT", "LMA", "N")
  } else if (PROSPECT_version == "PRO") {
    allParms <- c("CHL", "CAR", "ANT", "EWT", "PROT", "CBC", "N")
  } else if (PROSPECT_version == "PROB") {
    allParms <- c("CHL", "CAR", "ANT", "BROWN", "EWT", "PROT", "CBC", "N")
  } else {
    stop('PROSPECT_version not available. Choice is limited to "5", "5B", "D", "DB", "PRO", "PROB".')
  }
  if (alphaEst==TRUE){
    allParms <- c(allParms, "alpha")
  }
  if ("ALL" %in% Parms2Estimate) {
    Parms2Estimate <- allParms
  }

  Parms2Estimate <- allParms[allParms %in% Parms2Estimate]
  if (!all(allParms %in% names(xlub))) {
    stop('Some prospect parameters are missing in argument "InitValues".')
  }
  InitValues <- InitValues[allParms[allParms %in% names(InitValues)]]
  if (PROSPECT_version == "PRO" | PROSPECT_version == "PROB"){
    InitValues$LMA = 0
  }
  if (!all(Parms2Estimate %in% names(xlub))) {
    stop('Boundaries are missing for some parameters. Please make sure all parameters to estimate have a boundary defined in argument "xlub".')
  }
  xlub <- xlub[, Parms2Estimate]

  # update init value and lower/upper boundaries for inversion based on Vars2Estimate
  lb <- xlub[1, ]
  ub <- xlub[2, ]
  # run inversion procedure with standard parameterization
  res <- tryInversion(InitValues, MeritFunction, SpecPROSPECT, Refl, Tran, Parms2Estimate, lb, ub)

  names(res$par) = Parms2Estimate
  OutPROSPECT <- InitValues
  OutPROSPECT[names(res$par)] <- res$par
  return(OutPROSPECT)
}

#' Function handling error during inversion
#'
#' @param x0 numeric. Vector of input variables to estimate
#' @param MeritFunction  character. name of the function to be used as merit function
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param Refl  numeric. measured reflectance data
#' @param Tran  numeric. measured Transmittance data
#' @param Parms2Estimate  character vector. Parameters to estimate
#' to be estimated through inversion
#' @param lb numeric. Lower bound
#' @param ub numeric. Upper bound
#' @param verbose boolean. Set to TRUE if you want information about adjustment of tolerance during inversion.
#'
#' @return fc estimates of the parameters
#' @details
#' This function is based on \code{\link[pracma]{fmincon}}.
#' @importFrom pracma fmincon
#' @export

tryInversion <- function(x0, MeritFunction, SpecPROSPECT, Refl, Tran, Parms2Estimate, lb, ub, verbose = FALSE) {

  res <-list()
  res$par <- NA * c(1:length(Parms2Estimate))
  for (i in seq(-14,-2,1)){
    Tolerance = 10**(i)
    if (is.na(res$par[1])){
      res <- tryCatch(
        {
          res <- fmincon(
            x0 = as.numeric(x0[Parms2Estimate]), fn = MeritFunction, gr = NULL,
            SpecPROSPECT = SpecPROSPECT, Refl = Refl, Tran = Tran,
            Input_PROSPECT = x0, Parms2Estimate = Parms2Estimate,
            method = "SQP", A = NULL, b = NULL, Aeq = NULL, beq = NULL,
            lb = as.numeric(lb), ub = as.numeric(ub), hin = NULL, heq = NULL, tol = Tolerance,
            maxfeval = 2000, maxiter = 1000
          )
        },
        error = function(cond) {
          if (verbose){
            message('Adjusting Tolerance value for iterative optimization:  ',Tolerance)
          }
          res <- list()
          res$par <- NA * c(1:length(Parms2Estimate))
          return(res)
        },
        finally = {
        }
      )
    }
  }
  return(res)
}

#' Merit function for PROSPECT inversion
#'
#' @param x numeric. Vector of input variables to estimate
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param Refl  numeric. measured reflectance data
#' @param Tran  numeric. measured Transmittance data
#' @param Input_PROSPECT dataframe. set of PROSPECT input variables
#' @param Parms2Estimate  numeric. location of variables from Input_PROSPECT
#' to be estimated through inversion
#'
#' @return fc estimates of the parameters
#' @export
Merit_RMSE_PROSPECT <- function(x, SpecPROSPECT, Refl, Tran, Input_PROSPECT, Parms2Estimate) {
  x[x < 0] <- 0
  Input_PROSPECT[Parms2Estimate] <- x
  RT <- do.call("PROSPECT", c(list(SpecPROSPECT = SpecPROSPECT), Input_PROSPECT))
  fc <- CostVal_RMSE(RT, Refl, Tran)
  return(fc)
}

#' Value of the cost criterion to minimize during PROSPECT inversion
#' @param RT  list. Simulated reflectance and transmittance
#' @param Refl  numeric. Reflectance on which PROSPECT ins inverted
#' @param Tran  numeric. Transmittance on which PROSPECT ins inverted
#'
#' @return fc sum of squared difference between simulated and measured leaf optical properties
#' @export
CostVal_RMSE <- function(RT, Refl, Tran) {
  if (is.null(Tran)) {
    fc <- sqrt(sum((Refl - RT$Reflectance)**2) / length(RT$Reflectance))
  } else if (is.null(Refl)) {
    fc <- sqrt(sum((Tran - RT$Transmittance)**2) / length(RT$Transmittance))
  } else {
    fc <- sqrt(sum((Refl - RT$Reflectance)**2) / length(RT$Reflectance) + sum((Tran - RT$Transmittance)**2) / length(RT$Transmittance))
  }
  return(fc)
}

#' This function adapts SpecPROSPECT accordingly to experimental data
#' or to a spectral domain defined by UserDomain
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param lambda  numeric. Spectral bands corresponding to experimental data
#' @param Refl  numeric. Measured reflectance data
#' @param Tran  numeric. Measured Transmittance data
#' @param UserDomain  numeric. either Lower and upper bounds for domain of interest (optional)
#' or list of spectral bands of interest
#' @param UL_Bounds boolean. set to TRUE if UserDomain only includes lower and upper band,
#' set to FALSE if UserDomain is a list of spectral bands (in nm)
#'
#' @return res list including spectral properties at the new resolution
#' @importFrom utils tail
#' @export
FitSpectralData <- function(SpecPROSPECT, lambda, Refl = NULL, Tran = NULL, UserDomain = NULL, UL_Bounds = TRUE) {
  LowerPROSPECT <- SpecPROSPECT$lambda[1]
  UpperPROSPECT <- tail(SpecPROSPECT$lambda, n = 1)
  LowerLOP <- lambda[1]
  UpperLOP <- tail(lambda, n = 1)
  # if only need to fit PROSPECT data with spctral dat aprovided by user
  if (is.null(UserDomain)) {
    if (LowerPROSPECT > LowerLOP) {
      warning("________________________ WARNING _______________________")
      warning("User spectal data will be shrinked to start at the same ")
      warning("         Spectral band as SpecPROSPECT, which is        ")
      print(LowerPROSPECT)
      LowerBand_LOP <- which(abs(lambda - LowerPROSPECT) == min(abs(lambda - LowerPROSPECT)))
      LowerBand_Spec <- 1
    } else if (LowerPROSPECT < LowerLOP) {
      LowerBand_Spec <- which(abs(SpecPROSPECT$lambda - LowerLOP) == min(abs(SpecPROSPECT$lambda - LowerLOP)))
      LowerBand_LOP <- 1
    } else if (LowerPROSPECT == LowerLOP) {
      LowerBand_Spec <- 1
      LowerBand_LOP <- 1
    }
    if (UpperPROSPECT < UpperLOP) {
      warning("________________________ WARNING _______________________")
      warning("  User spectal data will be shrinked to end at the same ")
      warning("         Spectral band as SpecPROSPECT, which is        ")
      print(UpperPROSPECT)
      UpperBand_LOP <- which(abs(lambda - UpperPROSPECT) == min(abs(lambda - UpperPROSPECT)))
      UpperBand_Spec <- length(SpecPROSPECT$lambda)
    } else if (UpperPROSPECT > UpperLOP) {
      UpperBand_Spec <- which(abs(SpecPROSPECT$lambda - UpperLOP) == min(abs(SpecPROSPECT$lambda - UpperLOP)))
      UpperBand_LOP <- length(lambda)
    } else if (UpperPROSPECT == UpperLOP) {
      UpperBand_Spec <- length(SpecPROSPECT$lambda)
      UpperBand_LOP <- length(lambda)
    }
    # if user specifies a spectral domain which is different from PROSPECT and user data
  } else if (!is.null(UserDomain)) {
    LowerUser <- min(UserDomain)
    UpperUser <- max(UserDomain)
    if (LowerLOP > LowerUser | UpperLOP < UpperUser | LowerPROSPECT > LowerUser | UpperPROSPECT < UpperUser) {
      if (LowerPROSPECT > LowerUser | UpperPROSPECT < UpperUser) {
        warning("________________________ WARNING _______________________")
        warning("  The spectral domain defined in UserDomain provided as ")
        warning(" input in function FitSpectralData does not match with  ")
        warning("       the spectral domain covered by PROSPECT          ")
        warning("                                                        ")
        warning("                 PLEASE ADJUST UserDomain               ")
        warning("                                                        ")
        stop()
      }
      if (LowerLOP > LowerUser | UpperLOP < UpperUser) {
        warning("________________________ WARNING _______________________")
        warning("  The spectral domain defined in UserDomain provided as ")
        warning(" input in function FitSpectralData does not match with  ")
        warning("       the spectral domain covered by user data         ")
        warning("                                                        ")
        warning("                 PLEASE ADJUST UserDomain               ")
        warning("                                                        ")
        stop()
      }
    } else {
      if (LowerLOP <= LowerUser) {
        LowerBand_LOP <- which(abs(lambda - LowerUser) == min(abs(lambda - LowerUser)))
      }
      if (UpperLOP >= UpperUser) {
        UpperBand_LOP <- which(abs(lambda - UpperUser) == min(abs(lambda - UpperUser)))
      }
      if (LowerPROSPECT <= LowerUser) {
        LowerBand_Spec <- which(abs(SpecPROSPECT$lambda - LowerUser) == min(abs(SpecPROSPECT$lambda - LowerUser)))
      }
      if (UpperPROSPECT >= UpperUser) {
        UpperBand_Spec <- which(abs(SpecPROSPECT$lambda - UpperUser) == min(abs(SpecPROSPECT$lambda - UpperUser)))
      }
    }
  }
  SubSpecPROSPECT <- SpecPROSPECT[LowerBand_Spec:UpperBand_Spec, ]
  Sublambda <- lambda[LowerBand_LOP:UpperBand_LOP]
  if (!length(Sublambda) == length(SubSpecPROSPECT$lambda)) {
    warning("______________________ WARNING _____________________")
    warning("       PROSPECT expects 1nm spectal sampling        ")
    warning("The data provided as input shows unexpected sampling")
    warning("    Please prepare your data accordingly before     ")
    warning("             running PROSPECT inversion             ")
    warning("                The process will stop               ")
    stop()
  }
  SubRefl <- SubTran <- NULL
  if (!is.null(Refl)) {
    SubRefl <- Refl[LowerBand_LOP:UpperBand_LOP, ]
    if (ncol(Refl)==1){
      SubRefl <- matrix(SubRefl,ncol = 1)
    }
  }
  if (!is.null(Tran)) {
    SubTran <- Tran[LowerBand_LOP:UpperBand_LOP, ]
    if (ncol(Tran)==1){
      SubTran <- matrix(SubTran,ncol = 1)
    }
  }
  # in case a list of spectral bands has been provided, not only boundaries
  if (!UL_Bounds & !is.null(UserDomain)){
    # select spectral bands defined in UL_Bounds
    spectralBands <- unique(as.integer(UserDomain))
    SpectralLocation <- match(spectralBands,lambda)
    SubRefl <- SubRefl[SpectralLocation,]
    SubTran <- SubTran[SpectralLocation,]
    Sublambda <- Sublambda[SpectralLocation]
    SubSpecPROSPECT <- SubSpecPROSPECT[SpectralLocation,]
  }
  res <- list("SpecPROSPECT" = SubSpecPROSPECT, "lambda" = Sublambda, "Refl" = SubRefl, "Tran" = SubTran)
  return(res)
}

#' This function defines a regression model to estimate N from R only or T only
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param lambda  numeric. spectral bands corresponding to experimental data
#' @param Refl  numeric. Measured reflectance data
#' @param Tran  numeric. Measured Transmittance data
#'
#' @return Nprior vector corresponding to teh prior estimation of N based on R only or T only
#' @importFrom stats lm runif
#' @export
Get_Nprior <- function(SpecPROSPECT, lambda, Refl = NULL, Tran = NULL) {

  # definition of the optimal spectral band based on data available
  OptWL_R <- OptWL_T <- list()
  OptWL_R$NIR <- 800
  OptWL_R$SWIR <- 1131
  OptWL_T$NIR <- 753
  OptWL_T$SWIR <- 1121
  # if prior information based on Reflectance
  if (is.null(Tran)) {
    if (OptWL_R$SWIR %in% SpecPROSPECT$lambda) {
      OptWL <- OptWL_R$SWIR
    } else if (!OptWL_R$SWIR %in% SpecPROSPECT$lambda & OptWL_R$NIR %in% SpecPROSPECT$lambda) {
      warning("________________________ WARNING _______________________")
      warning("The optimal prior estimation of N using Reflectance only")
      warning("requires information at 1131nm.")
      warning("The reflectance does not include this spectral band")
      warning("Using reflectance at 800 nm instead")
      OptWL <- OptWL_R$NIR
    } else if (!OptWL_R$SWIR %in% SpecPROSPECT$lambda & !OptWL_R$NIR %in% SpecPROSPECT$lambda) {
      warning("________________________ WARNING _______________________")
      warning("The spectral information of the reflectance provided here")
      warning("does not contain the spectral bands required to estimate")
      warning("prior information about N.")
      warning("The proces will stop")
      stop()
    }
  } else if (is.null(Refl)) {
    if (OptWL_T$SWIR %in% SpecPROSPECT$lambda) {
      OptWL <- OptWL_T$SWIR
    } else if (!OptWL_T$SWIR %in% SpecPROSPECT$lambda & OptWL_T$NIR %in% SpecPROSPECT$lambda) {
      warning("________________________ WARNING _______________________")
      warning("The optimal prior estimation of N using Transmittance only")
      warning("requires information at 1121nm.")
      warning("The Transmittance does not include this spectral band")
      warning("Using Transmittance at 753 nm instead")
      OptWL <- OptWL_T$NIR
    } else if (!OptWL_T$SWIR %in% SpecPROSPECT$lambda & !OptWL_T$NIR %in% SpecPROSPECT$lambda) {
      warning("________________________ WARNING _______________________")
      warning("The spectral information of the Transmittance provided here")
      warning("does not contain the spectral bands required to estimate")
      warning("prior information about N.")
      warning("The proces will stop")
      stop()
    }
  }

  # get the subdomain corresponding to OptWL
  SubData <- FitSpectralData(SpecPROSPECT = SpecPROSPECT, lambda = lambda, Refl = Refl, Tran = Tran, UserDomain = c(OptWL, OptWL))
  SubSpecPROSPECT <- SubData$SpecPROSPECT
  Sublambda <- SubData$lambda
  SubRefl <- SubData$Refl
  subTran <- SubData$Tran

  # create a LUT using only the spectral band of interest
  nbSim <- 1000
  CHL <- 0.5 + 100 * runif(nbSim)
  CAR <- 0.5 + 20 * runif(nbSim)
  EWT <- 0.001 + 0.02 * runif(nbSim)
  LMA <- 0.001 + 0.01 * runif(nbSim)
  N <- 1 + 1.5 * runif(nbSim)
  Input_PROSPECT <- data.frame(CHL, CAR, EWT, LMA, N)
  LUT <- PROSPECT_LUT(SubSpecPROSPECT, Input_PROSPECT)
  # fit a linear model between
  if (is.null(Tran)) {
    Ratio <- LUT$Reflectance / (1 - LUT$Reflectance)
    Ratio_Meas <- SubRefl / (1 - SubRefl)
  } else if (is.null(Refl)) {
    Ratio <- (1 - LUT$Transmittance) / LUT$Transmittance
    Ratio_Meas <- (1 - subTran) / subTran
  }
  N_Model <- lm(matrix(LUT$Input_PROSPECT$N) ~ matrix(Ratio))
  NpriorMOD <- N_Model$coefficients[2] * Ratio + N_Model$coefficients[1]
  Nprior <- N_Model$coefficients[2] * Ratio_Meas + N_Model$coefficients[1]
  return(Nprior)
}


#' This function uses optimal configuration identified by Spafford et al. (2020) to estimate leaf chemistry
#' prior information on N is provided if only Reflectance or only Transmittance is available
#' @param SpecPROSPECT list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#' @param lambda  numeric. spectral bands corresponding to experimental data
#' @param Refl  numeric. Measured reflectance data
#' @param Tran  numeric. Measured Transmittance data
#' @param PROSPECT_version  character. Version of prospect model used for the inversion: '5', '5B', 'D', 'DB', 'PRO', 'PROB',
#' See details.
#' @param Parms2Estimate  character vector. Parameters to estimate (can be 'ALL')
#' @param InitValues  data.frame. Default values of PROSPECT parameters. During optimization,
#' they are used either as initialization values for parameters to estimate,
#' or as fix values for other parameters.
#' Parameters not compatible with PROSPECT_version are not taken into account.
#'
#' @return Nprior vector corresponding to teh prior estimation of N based on R only or T only
#' @importFrom stats lm runif
#' @export
Invert_PROSPECT_OPT <- function(SpecPROSPECT, lambda, Refl = NULL, Tran = NULL,
                                PROSPECT_version = 'D',Parms2Estimate = 'ALL',
                                InitValues = data.frame(
                                  CHL = 40, CAR = 10, ANT = 0.1, BROWN = 0.01, EWT = 0.01,
                                  LMA = 0.01, PROT = 0.001, CBC = 0.009, N = 1.5, alpha = 40
                                )) {

  # define optimal domain for the different constituents
  OptDomain_RT <- list('CHL' = c(700,720), 'CAR' = c(520,560), 'ANT' = c(400,800),
                      'EWT' = c(1700,2400), 'LMA' = c(1700,2400), 'PROT' = c(2125,2175),
                      'CBC' = c(2025,2350))
  OptDomain_R <- list('CHL' = c(700,720), 'CAR' = c(400,800), 'ANT' = c(400,800),
                      'EWT' = c(1200,2400), 'LMA' = c(1700,2400), 'PROT' = c(2125,2175),
                      'CBC' = c(2025,2350))
  OptDomain_T <- list('CHL' = c(700,720), 'CAR' = c(520,560), 'ANT' = c(400,800),
                      'EWT' = c(1700,2400), 'LMA' = c(1700,2400), 'PROT' = c(2125,2175),
                      'CBC' = c(2025,2350))

  # adjust depending on avalibility of R and T
  if (is.null(Refl)){
    OptDomain <- OptDomain_T
  } else if (is.null(Tran)){
      OptDomain <- OptDomain_R
  }  else {
    OptDomain <- OptDomain_RT
  }

  # define PROSPECT input parameters
  # if all parameters Parms2Estimate = "ALL"
  # Which parameters for PROSPECT_version?
  if (PROSPECT_version == "5") {
    allParms <- c("CHL", "CAR", "EWT", "LMA")
    ANTinit <- 0.0
  } else if (PROSPECT_version == "D") {
    allParms <- c("CHL", "CAR", "ANT", "EWT", "LMA")
    ANTinit <- 0.1
  } else if (PROSPECT_version == "PRO") {
    allParms <- c("CHL", "CAR", "ANT", "EWT", "PROT", "CBC")
    ANTinit <- 0.1
  } else if (PROSPECT_version == "5B" | PROSPECT_version == "DB" | PROSPECT_version == "PROB") {
    message('brown pigments are not accounted for when performing optimal estimation of leaf chemistry')
    message('model version excluding brown pigments will be used instead')
    if (PROSPECT_version == "5B") {
      allParms <- c("CHL", "CAR", "EWT", "LMA")
      PROSPECT_version <- "5"
      ANTinit <- 0.0
    } else if (PROSPECT_version == "DB") {
      allParms <- c("CHL", "CAR", "ANT", "EWT", "LMA")
      PROSPECT_version <- "D"
      ANTinit <- 0.1
    } else if (PROSPECT_version == "PROB") {
      allParms <- c("CHL", "CAR", "ANT", "EWT", "PROT", "CBC")
      PROSPECT_version <- "PRO"
      ANTinit <- 0.1
    }
  } else {
    stop('PROSPECT_version not available. Choice is limited to "5", "5B", "D", "DB", "PRO", "PROB".')
  }
  if ("ALL" %in% Parms2Estimate) {
    Parms2Estimate <- allParms
  }
  # get the number of samples to be processed
  if (!is.null(Refl)){
    nbSamples <- ncol(Refl)
  } else if (!is.null(Tran)){
      nbSamples <- ncol(Tran)
  }

  # if reflectance or transmittance is not available
  if (is.null(Refl) | is.null(Tran)){
    # compute prior estimate of N
    message('computing prior estimation of N as both R & T are not provided')
    Nprior <- Get_Nprior(SpecPROSPECT, lambda, Refl = Refl, Tran = Tran)
    # estimation of chlorophylls: Prior N and optimal spectral domain
    ParmEst = list()
    ParmEst$N = Nprior
    for (parm in Parms2Estimate){
      ParmEst[[parm]] = c()
      if (parm == "ANT"){
        message('Currently no optimal estimation for anthocyanins')
        message('PROSPECT inversion will be performed using full spectral information')
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        if (!PROSPECT_version =='5'){Parms2EstimateBis <- c('CHL', 'CAR', 'ANT')}
        else if (PROSPECT_version =='5'){Parms2EstimateBis <- c('CHL', 'CAR')}
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          InitValues <- data.frame(CHL=40, CAR=10, ANT=ANTinit, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior[i])
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_version,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$ANT
        }
      }
      if (parm == "CHL"){
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        if (!PROSPECT_version =='5'){Parms2EstimateBis <- c('CHL', 'CAR', 'ANT')}
        else if (PROSPECT_version =='5'){Parms2EstimateBis <- c('CHL', 'CAR')}
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          InitValues <- data.frame(CHL=40, CAR=10, ANT=ANTinit, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior[i])
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_version,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$CHL
        }
      }
      if (parm == "CAR"){
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        if (!PROSPECT_version =='5'){Parms2EstimateBis <- c('CHL', 'CAR', 'ANT')}
        else if (PROSPECT_version =='5'){Parms2EstimateBis <- c('CHL', 'CAR')}
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          InitValues <- data.frame(CHL=40, CAR=10, ANT=ANTinit, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior[i])
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_version,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$CAR
        }
      }
      if (parm == "EWT"){
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        if (!PROSPECT_version =='PRO'){Parms2EstimateBis <- c('EWT', 'LMA')}
        else if (PROSPECT_version =='PRO'){Parms2EstimateBis <- c('EWT', 'PROT', 'CBC')}
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          if (!PROSPECT_version =='PRO'){
            InitValues <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior[i])
          } else if (!PROSPECT_version =='PRO'){
            InitValues <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01, LMA=0.00, PROT=0.001, PROT=0.009, N=Nprior[i])
          }
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_version,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$EWT
        }
      }
      if (parm == "LMA"){
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        if (PROSPECT_version == 'PRO'){
          message('LMA is not estimated directly using PROSPECT-PRO')
          message('Please run inversion for PROSPECT-PRO and sum PROT and CBC')
          message('if you want to get estimated LMA from PROSPECT-PRO')
          message('')
        } else {
          # Fit spectral data to match PROSPECT with user optical properties
          SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
          SubSpecPROSPECT = SubData$SpecPROSPECT
          Sublambda       = SubData$lambda
          SubRefl         = SubData$Refl
          SubTran         = SubData$Tran
          # Parms2EstimateBis <- c('EWT', 'PROT', 'CBC')
          Parms2EstimateBis <- c('EWT', 'LMA')
          # Invert PROSPECT with optimal spectral information
          for (i in 1:nbSamples){
            print(i)
            InitValues <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior[i])
            res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_version,
                                   Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
            ParmEst[[parm]][i] = res$LMA
          }
        }
      }
      if (parm == "PROT"){
        PROSPECT_versionPRO = 'PRO'
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        Parms2EstimateBis <- c('EWT', 'PROT', 'CBC')
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          InitValues <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01, LMA=0.00, PROT=0.001, CBC=0.009, N=Nprior[i])
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_versionPRO,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$PROT
        }
      }
      if (parm == "CBC"){
        PROSPECT_versionPRO = 'PRO'
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        Parms2EstimateBis <- c('EWT', 'PROT', 'CBC')
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          InitValues <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01, LMA=0.00, PROT=0.001, CBC=0.009, N=Nprior[i])
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_versionPRO,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$PROT
        }
      }
    }
  # if reflectance and transmittance are available
  } else {
    Nprior <- Get_Nprior(SpecPROSPECT, lambda, Refl = NULL, Tran = Tran)
    ParmEst = list()
    ParmEst$N = Nprior
    for (parm in Parms2Estimate){
      ParmEst[[parm]] = c()
      if (parm == "ANT"){
        message('Currently no optimal estimation for anthocyanins')
        message('PROSPECT inversion will be performed using full spectral information')
        message('Optimal estimation of ',parm)
        message('Using optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        if (!PROSPECT_version =='5'){Parms2EstimateBis <- c('CHL', 'CAR', 'ANT','N')}
        else if (PROSPECT_version =='5'){Parms2EstimateBis <- c('CHL', 'CAR','N')}
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          InitValues <- data.frame(CHL=40, CAR=10, ANT=ANTinit, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior[i])
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_version,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$ANT
        }
      }
      if (parm == "CHL"){
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        if (!PROSPECT_version =='5'){
          Parms2EstimateBis <- c('CHL', 'CAR', 'ANT','N')
        } else if (PROSPECT_version =='5'){
            Parms2EstimateBis <- c('CHL', 'CAR','N')
        }
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          InitValues <- data.frame(CHL=40, CAR=10, ANT=ANTinit, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior[i])
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_version,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$CHL
        }
      }
      if (parm == "CAR"){
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        if (!PROSPECT_version =='5'){Parms2EstimateBis <- c('CHL', 'CAR', 'ANT','N')}
        else if (PROSPECT_version =='5'){Parms2EstimateBis <- c('CHL', 'CAR','N')}
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          InitValues <- data.frame(CHL=40, CAR=10, ANT=ANTinit, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior[i])
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_version,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$CAR
        }
      }
      if (parm == "EWT"){
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        if (!PROSPECT_version =='PRO'){Parms2EstimateBis <- c('EWT', 'LMA','N')}
        else if (PROSPECT_version =='PRO'){Parms2EstimateBis <- c('EWT', 'PROT', 'CBC','N')}
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          if (!PROSPECT_version =='PRO'){
            InitValues <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior[i])
          } else if (!PROSPECT_version =='PRO'){
            InitValues <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01, LMA=0.00, PROT=0.001, PROT=0.009, N=Nprior[i])
          }
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_version,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$EWT
        }
      }
      if (parm == "LMA"){
        message('Optimal estimation of ',parm)
        message('Using Reflectance only, no prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        if (PROSPECT_version == 'PRO'){
          message('LMA is not estimated directly using PROSPECT-PRO')
          message('Please run inversion for PROSPECT-PRO and sum PROT and CBC')
          message('if you want to get estimated LMA from PROSPECT-PRO')
          message('')
        } else {
          # Fit spectral data to match PROSPECT with user optical properties
          SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
          SubSpecPROSPECT = SubData$SpecPROSPECT
          Sublambda       = SubData$lambda
          SubRefl         = SubData$Refl
          SubTran         = SubData$Tran
          Parms2EstimateBis <- c('EWT', 'LMA', 'N')
          # Invert PROSPECT with optimal spectral information
          for (i in 1:nbSamples){
            print(i)
            InitValues <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior[i])
            res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_version,
                                   Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
            ParmEst[[parm]][i] = res$LMA
          }
        }
      }
      if (parm == "PROT"){
        PROSPECT_versionPRO = 'PRO'
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        Parms2EstimateBis <- c('EWT', 'PROT', 'CBC','N')
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          InitValues <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01, LMA=0.00, PROT=0.001, CBC=0.009, N=Nprior[i])
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_versionPRO,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$PROT
        }
      }
      if (parm == "CBC"){
        PROSPECT_versionPRO = 'PRO'
        message('Optimal estimation of ',parm)
        message('Using Prior N value and optimal spectral subdomain (nm) :')
        message('from ', OptDomain[[parm]][1],' nm to ', OptDomain[[parm]][2],' nm')
        message('')
        # Fit spectral data to match PROSPECT with user optical properties
        SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = OptDomain[[parm]])
        SubSpecPROSPECT = SubData$SpecPROSPECT
        Sublambda       = SubData$lambda
        SubRefl         = SubData$Refl
        SubTran         = SubData$Tran
        Parms2EstimateBis <- c('EWT', 'PROT', 'CBC','N')
        # Invert PROSPECT with optimal spectral information
        for (i in 1:nbSamples){
          print(i)
          InitValues <- data.frame(CHL=0, CAR=0, ANT=0, BROWN=0, EWT=0.01, LMA=0.00, PROT=0.001, CBC=0.009, N=Nprior[i])
          res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], PROSPECT_version = PROSPECT_versionPRO,
                                 Parms2Estimate = Parms2EstimateBis, InitValues = InitValues)
          ParmEst[[parm]][i] = res$PROT
        }
      }
    }
  }
  return(ParmEst)
}
