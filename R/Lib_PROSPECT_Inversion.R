# ==============================================================================
# prospect
# Lib_PROSPECT_Inversion.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Florian de Boissieu <fdeboiss@gmail.com>
# Copyright 2019/11 Jean-Baptiste FERET
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
#' @param x0  data.frame. Default values of PROSPECT parameters. During optimization,
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
#'
#' @importFrom pracma fmincon
#' @export
#' @md
Invert_PROSPECT  <- function(SpecPROSPECT, Refl = NULL, Tran = NULL,
                             x0 = data.frame(CHL = 40, CAR = 10, ANT = 0.1, BROWN = 0.01, EWT = 0.01,
                                             LMA = 0.008, PROT = 0.001, CBC = 0.008, N = 1.5, alpha = 40),
                             Parms2Estimate = 'ALL',
                             PROSPECT_version = 'D',
                             MeritFunction = 'Merit_RMSE_PROSPECT',
                             xlub = data.frame(CHL = c(1e-4, 150), CAR = c(1e-4, 25), ANT = c(0, 20),
                                               BROWN = c(0, 1), EWT = c(1e-7, .08), LMA = c(1e-7, .4),
                                               PROT = c(0, .005), CBC = c(0, .04), N = c(.5, 3),
                                               alpha = c(10, 90))){

  # define PROSPECT input parameters
  if (PROSPECT_version == '5'){
    allParms = c('CHL', 'CAR', 'EWT', 'LMA', 'N')
  } else if (PROSPECT_version == '5B'){
    allParms = c('CHL', 'CAR', 'BROWN', 'EWT', 'LMA', 'N')
  } else if (PROSPECT_version == 'D'){
    allParms = c('CHL', 'CAR', 'ANT', 'EWT', 'LMA', 'N')
  } else if (PROSPECT_version == 'DB'){
    allParms = c('CHL', 'CAR', 'ANT', 'BROWN', 'EWT', 'LMA', 'N')
  } else if (PROSPECT_version == 'PRO'){
    allParms = c('CHL', 'CAR', 'ANT', 'EWT', 'PROT', 'CBC', 'N')
  } else if (PROSPECT_version == 'PROB'){
    allParms = c('CHL', 'CAR', 'ANT', 'BROWN', 'EWT', 'PROT', 'CBC', 'N')
  }else{
    stop('PROSPECT_version not available. Choice is limited to "5", "5B", "D", "DB", "PRO", "PROB".')
  }


  if ('ALL'%in%Parms2Estimate)
    Parms2Estimate = allParms

  Parms2Estimate = allParms[ allParms %in% Parms2Estimate]
  if(!all(allParms %in% names(xlub)))
    stop('Some prospect parameters are missing in argument "x0".')
  x0 = x0[allParms[allParms %in% names(x0)]]
  if(!all(Parms2Estimate %in% names(xlub)))
    stop('Boundaries are missing for some parameters. Please make sure all parameters to estimate have a boundary defined in argument "xlub".')
  xlub = xlub[, Parms2Estimate]

  # update init value and lower/upper boundaries for inversion based on Vars2Estimate
  lb    = xlub[1,]
  ub    = xlub[2,]
  # run inversion procedure with standard parameterization
  res <- tryInversion(x0,MeritFunction,SpecPROSPECT,Refl,Tran,Parms2Estimate,lb,ub)

  OutPROSPECT  = x0
  OutPROSPECT[names(res$par)]= res$par
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
#'
#' @return fc estimates of the parameters
#' @details
#' This function is based on \code{\link[pracma]{fmincon}}.
#' @importFrom pracma fmincon
#' @export

tryInversion <- function(x0, MeritFunction, SpecPROSPECT, Refl, Tran, Parms2Estimate, lb, ub) {
  res <- tryCatch(
    {
      res   = fmincon(x0 = as.numeric(x0[Parms2Estimate]), fn = MeritFunction, gr = NULL,
                      SpecPROSPECT = SpecPROSPECT, Refl = Refl, Tran = Tran,
                      Input_PROSPECT = x0, Parms2Estimate = Parms2Estimate,
                      method = "SQP", A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                      lb = as.numeric(lb), ub = as.numeric(ub), hin = NULL, heq = NULL,tol = 1e-15,
                      maxfeval = 2000, maxiter = 1000)
    },
    error=function(cond) {
      message("Inversion failed on a sample")
      message("Error message obtained:")
      message(cond)
      message("")
      message("NA values will be set for this sample")
      # Choose a return value in case of error
      res <- list()
      res$par <- NA*c(1:length(Parms2Estimate))
      return(res)
    },
    finally={
    }
  )
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

  x[x<0] = 0
  Input_PROSPECT[Parms2Estimate] = x
  RT = do.call('PROSPECT', c(list(SpecPROSPECT = SpecPROSPECT), Input_PROSPECT))
  fc = CostVal_RMSE(RT, Refl, Tran)
  return(fc)
}


#' Value of the cost criterion to minimize during PROSPECT inversion
#' @param RT  list. Simulated reflectance and transmittance
#' @param Refl  numeric. Reflectance on which PROSPECT ins inverted
#' @param Tran  numeric. Transmittance on which PROSPECT ins inverted
#'
#' @return fc sum of squared difference between simulated and measured leaf optical properties
#' @export
CostVal_RMSE  <- function(RT, Refl, Tran) {

  if (is.null(Tran)){
    fc = sqrt(sum((Refl-RT$Reflectance)**2)/length(RT$Reflectance))
  } else if (is.null(Refl)){
    fc = sqrt(sum((Tran-RT$Transmittance)**2)/length(RT$Transmittance))
  } else {
    fc=sqrt(sum((Refl-RT$Reflectance)**2)/length(RT$Reflectance)+sum((Tran-RT$Transmittance)**2)/length(RT$Transmittance))
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
#' @param UserDomain  numeric. Lower and upper bounds for domain of interest (optional)
#'
#' @return res list including spectral properties at the new resolution
#' @importFrom utils tail
#' @export
FitSpectralData <- function(SpecPROSPECT,lambda,Refl=NULL,Tran=NULL,UserDomain=NULL) {

  LowerPROSPECT = SpecPROSPECT$lambda[1]
  UpperPROSPECT = tail(SpecPROSPECT$lambda, n=1)
  LowerLOP      = lambda[1]
  UpperLOP      = tail(lambda, n=1)
  # if only need to fit PROSPECT data with spctral dat aprovided by user
  if (is.null(UserDomain)){
    if (LowerPROSPECT>LowerLOP){
      warning('________________________ WARNING _______________________')
      warning('User spectal data will be shrinked to start at the same ')
      warning('         Spectral band as SpecPROSPECT, which is        ')
      print(LowerPROSPECT)
      LowerBand_LOP = which(abs(lambda-LowerPROSPECT)==min(abs(lambda-LowerPROSPECT)))
      LowerBand_Spec= 1
    } else  if (LowerPROSPECT<LowerLOP){
      LowerBand_Spec  = which(abs(SpecPROSPECT$lambda-LowerLOP)==min(abs(SpecPROSPECT$lambda-LowerLOP)))
      LowerBand_LOP = 1
    } else  if (LowerPROSPECT==LowerLOP){
      LowerBand_Spec= 1
      LowerBand_LOP = 1
    }
    if (UpperPROSPECT<UpperLOP){
      warning('________________________ WARNING _______________________')
      warning('  User spectal data will be shrinked to end at the same ')
      warning('         Spectral band as SpecPROSPECT, which is        ')
      print(UpperPROSPECT)
      UpperBand_LOP = which(abs(lambda-UpperPROSPECT)==min(abs(lambda-UpperPROSPECT)))
      UpperBand_Spec= length(SpecPROSPECT$lambda)
    } else  if (UpperPROSPECT>UpperLOP){
      UpperBand_Spec  = which(abs(SpecPROSPECT$lambda-UpperLOP)==min(abs(SpecPROSPECT$lambda-UpperLOP)))
      UpperBand_LOP = length(lambda)
    } else  if (UpperPROSPECT==UpperLOP){
      UpperBand_Spec= length(SpecPROSPECT$lambda)
      UpperBand_LOP = length(lambda)
    }
    # if user specifies a spectral domain which is different from PROSPECT and user data
  } else if (!is.null(UserDomain)){
    LowerUser     = UserDomain[1]
    UpperUser     = UserDomain[2]
    if (LowerLOP>LowerUser | UpperLOP<UpperUser | LowerPROSPECT>LowerUser | UpperPROSPECT<UpperUser){
      if (LowerPROSPECT>LowerUser | UpperPROSPECT<UpperUser){
        warning('________________________ WARNING _______________________')
        warning('  The spectral domain defined in UserDomain provided as ')
        warning(' input in function FitSpectralData does not match with  ')
        warning('       the spectral domain covered by PROSPECT          ')
        warning('                                                        ')
        warning('                 PLEASE ADJUST UserDomain               ')
        warning('                                                        ')
        stop()
      }
      if (LowerLOP>LowerUser | UpperLOP<UpperUser){
        warning('________________________ WARNING _______________________')
        warning('  The spectral domain defined in UserDomain provided as ')
        warning(' input in function FitSpectralData does not match with  ')
        warning('       the spectral domain covered by user data         ')
        warning('                                                        ')
        warning('                 PLEASE ADJUST UserDomain               ')
        warning('                                                        ')
        stop()
      }
    } else {
      if (LowerLOP<=LowerUser){
        LowerBand_LOP  = which(abs(lambda-LowerUser)==min(abs(lambda-LowerUser)))
      }
      if (UpperLOP>=UpperUser){
        UpperBand_LOP  = which(abs(lambda-UpperUser)==min(abs(lambda-UpperUser)))
      }
      if (LowerPROSPECT<=LowerUser){
        LowerBand_Spec  = which(abs(SpecPROSPECT$lambda-LowerUser)==min(abs(SpecPROSPECT$lambda-LowerUser)))
      }
      if (UpperPROSPECT>=UpperUser){
        UpperBand_Spec  = which(abs(SpecPROSPECT$lambda-UpperUser)==min(abs(SpecPROSPECT$lambda-UpperUser)))
      }
    }
  }
  SubSpecPROSPECT = SpecPROSPECT[LowerBand_Spec:UpperBand_Spec,]
  Sublambda       = lambda[LowerBand_LOP:UpperBand_LOP]
  if (!length(Sublambda)==length(SubSpecPROSPECT$lambda)){
    warning('______________________ WARNING _____________________')
    warning('       PROSPECT expects 1nm spectal sampling        ')
    warning('The data provided as input shows unexpected sampling')
    warning('    Please prepare your data accordingly before     ')
    warning('             running PROSPECT inversion             ')
    warning('                The process will stop               ')
    stop()
  }
  SubRefl = SubTran = NULL
  if (!is.null(Refl)){
    SubRefl = Refl[LowerBand_LOP:UpperBand_LOP,]
  }
  if (!is.null(Tran)){
    SubTran = Tran[LowerBand_LOP:UpperBand_LOP,]
  }
  res       = list('SpecPROSPECT'=SubSpecPROSPECT,'lambda'=Sublambda,'Refl'=SubRefl,'Tran'=SubTran)
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
Get_Nprior <- function(SpecPROSPECT,lambda,Refl=NULL,Tran=NULL) {

  # definition of the optimal spectral band based on data available
  OptWL_R = OptWL_T = list()
  OptWL_R$NIR     = 800
  OptWL_R$SWIR    = 1131
  OptWL_T$NIR     = 753
  OptWL_T$SWIR    = 1121
  # if prior information based on Reflectance
  if (is.null(Tran)){
    if (OptWL_R$SWIR%in%SpecPROSPECT$lambda){
      OptWL = OptWL_R$SWIR
    } else if (!OptWL_R$SWIR%in%SpecPROSPECT$lambda & OptWL_R$NIR%in%SpecPROSPECT$lambda ){
      warning('________________________ WARNING _______________________')
      warning('The optimal prior estimation of N using Reflectance only')
      warning('requires information at 1131nm.')
      warning('The reflectance does not include this spectral band')
      warning('Using reflectance at 800 nm instead')
      OptWL = OptWL_R$NIR
    } else if (!OptWL_R$SWIR%in%SpecPROSPECT$lambda & !OptWL_R$NIR%in%SpecPROSPECT$lambda ){
      warning('________________________ WARNING _______________________')
      warning('The spectral information of the reflectance provided here')
      warning('does not contain the spectral bands required to estimate')
      warning('prior information about N.')
      warning('The proces will stop')
      stop()
    }
  } else if (is.null(Refl)){
    if (OptWL_T$SWIR%in%SpecPROSPECT$lambda){
      OptWL = OptWL_T$SWIR
    } else if (!OptWL_T$SWIR%in%SpecPROSPECT$lambda & OptWL_T$NIR%in%SpecPROSPECT$lambda ){
      warning('________________________ WARNING _______________________')
      warning('The optimal prior estimation of N using Transmittance only')
      warning('requires information at 1121nm.')
      warning('The Transmittance does not include this spectral band')
      warning('Using Transmittance at 753 nm instead')
      OptWL = OptWL_T$NIR
    } else if (!OptWL_T$SWIR%in%SpecPROSPECT$lambda & !OptWL_T$NIR%in%SpecPROSPECT$lambda ){
      warning('________________________ WARNING _______________________')
      warning('The spectral information of the Transmittance provided here')
      warning('does not contain the spectral bands required to estimate')
      warning('prior information about N.')
      warning('The proces will stop')
      stop()
    }
  }

  # get the subdomain corresponding to OptWL
  SubData = FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl,Tran=Tran,UserDomain = c(OptWL,OptWL))
  SubSpecPROSPECT = SubData$SpecPROSPECT
  Sublambda       = SubData$lambda
  SubRefl         = SubData$Refl
  subTran         = SubData$Tran

  # create a LUT using only the spectral band of interest
  nbSim   = 1000
  CHL     = 0.5+100*runif(nbSim)
  CAR     = 0.5+20*runif(nbSim)
  EWT     = 0.001+0.02*runif(nbSim)
  LMA     = 0.001+0.01*runif(nbSim)
  N       = 1+1.5*runif(nbSim)
  Input_PROSPECT = data.frame('CHL'=CHL,'CAR'=CAR,'EWT'=EWT,'LMA'=LMA,'N'=N)
  LUT   = PROSPECT_LUT(SubSpecPROSPECT,Input_PROSPECT)
  # fit a linear model between
  if (is.null(Tran)){
    Ratio       = LUT$Reflectance/(1-LUT$Reflectance)
    Ratio_Meas  = SubRefl/(1-SubRefl)
  } else if (is.null(Refl)){
    Ratio       = (1-LUT$Transmittance)/LUT$Transmittance
    Ratio_Meas  = (1-subTran)/subTran
  }
  N_Model   = lm(matrix(LUT$Input_PROSPECT$N) ~ matrix(Ratio))
  NpriorMOD = N_Model$coefficients[2]*Ratio+N_Model$coefficients[1]
  Nprior    = N_Model$coefficients[2]*Ratio_Meas+N_Model$coefficients[1]
  return(Nprior)
}
