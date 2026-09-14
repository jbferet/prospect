#' Deprecated functions in prospect
#'
#' These functions are provided for compatibility with older versions of prospect
#' but are deprecated.
#'
#' @param x0,MeritFunction,SpecPROSPECT,Refl,Tran,lb,ub,verbose See the new functions that replace the old ones
#' @param Input_PROSPECT,CHL,CAR,ANT,BROWN,EWT,LMA,PROT,CBC,N,alpha See the new functions that replace the old ones
#' @param lambda,UserDomain,UL_Bounds See the new functions that replace the old ones
#' @param InitValues,Parms2Estimate,PROSPECT_version,xlub,Est_Brown_Pigments,Est_alpha,progressBar See the new functions that replace the old ones
#' @param OptWL_R,OptWL_T,check,ParmEst,Nprior,ANTinit,OptDomain See the new functions that replace the old ones
#'
#'
#' @rdname deprecated
#' @name deprecated
#'
NULL




#' @rdname deprecated
#' @export
#'
define_Input_PROSPECT <- function(Input_PROSPECT, CHL = NULL, CAR = NULL,
                                  ANT = NULL, BROWN = NULL, EWT = NULL,
                                  LMA = NULL, PROT = NULL, CBC = NULL,
                                  N = NULL, alpha = NULL){
  .Deprecated(old = "define_Input_PROSPECT",
              new = "define_input_prospect")
  define_input_prospect(input_prospect = Input_PROSPECT, chl = CHL, car = CAR,
                        ant = ANT, brown = BROWN, ewt = EWT, lma = LMA,
                        prot = PROT, cbc = CBC, n_struct = N, alpha = alpha)
}

#' @rdname deprecated
#' @export
#'
FitSpectralData <- function(lambda, SpecPROSPECT = NULL,
                            Refl = NULL, Tran = NULL,
                            UserDomain = NULL, UL_Bounds = FALSE) {

  .Deprecated(old = "FitSpectralData",
              new = "fit_spectral_data")
  fit_spectral_data(lambda = lambda, spec_prospect = SpecPROSPECT,
                    refl = Refl, tran = Tran, user_domain = UserDomain,
                    ul_bounds = UL_Bounds)
}


#' @rdname deprecated
#' @export
#'
Get_Nprior <- function(lambda, SpecPROSPECT = NULL, Refl = NULL, Tran = NULL,
                       OptWL_R = list(NIR = 800,SWIR = 1131),
                       OptWL_T = list(NIR = 753,SWIR = 1121)) {

  .Deprecated(old = "Get_Nprior",
              new = "get_N_prior")
  get_N_prior(lambda = lambda, spec_prospect = SpecPROSPECT,
              refl = Refl, tran = Tran, OptWL_R = OptWL_R, OptWL_T = OptWL_T)
}


#' @rdname deprecated
#' @export
#'
Invert_PROSPECT <- function(
    SpecPROSPECT = NULL, Refl = NULL, Tran = NULL,
    InitValues = data.frame(chl = 40, car = 10, ant = 0.1, brown = 0.0,
                            ewt = 0.01, lma = 0.01, prot = 0.001, cbc = 0.009,
                            n_struct = 1.5, alpha = 40),
    Parms2Estimate = "ALL", PROSPECT_version = "D",
    MeritFunction = "merit_prospect_rmse",
    xlub = data.frame(chl = c(1e-4, 150), car = c(1e-4, 25), ant = c(0, 50),
                      brown = c(0, 4), ewt = c(1e-8, 0.1), lma = c(1e-8, 0.06),
                      prot = c(1e-7, 0.006), cbc = c(1e-6, 0.054),
                      n_struct = c(0.5, 4), alpha = c(10, 90)),
    Est_Brown_Pigments = FALSE, Est_alpha = FALSE, verbose = FALSE,
    progressBar = TRUE) {

  .Deprecated(old = "Invert_PROSPECT",
              new = "invert_prospect")
  options <- set_options_prospect(fun = 'invert_prospect')
  options$spec_prospect <- SpecPROSPECT
  options$init_values <- InitValues
  options$merit_function <- MeritFunction
  options$xlub <- xlub
  options$init_values <- InitValues
  options$progressBar <- progressBar
  options$verbose <- verbose
  options$estimate_alpha <- Est_alpha
  options$estimate_brown_pigments <- Est_Brown_Pigments
  invert_prospect(refl = Refl, tran = Tran, parms_to_estimate = Parms2Estimate,
                  prospect_version = PROSPECT_version, options = options)
}



#' @rdname deprecated
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

  .Deprecated(old = "Invert_PROSPECT_OPT",
              new = "invert_prospect_opt")
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


#' @rdname deprecated
#' @export
#'
PROSPECT <- function(SpecPROSPECT = NULL, Input_PROSPECT = NULL,
                     N = 1.5, CHL = 40.0, CAR = 8.0, ANT = 0.0, BROWN = 0.0,
                     EWT = 0.01, LMA = NULL, PROT = 0, CBC = 0, alpha = 40.0,
                     check = TRUE) {

  .Deprecated(old = "PROSPECT",
              new = "prospect")
  prospect(spec_prospect = SpecPROSPECT, input_prospect = Input_PROSPECT,
           n_struct = N, chl = CHL, car = CAR, ant = ANT, brown = BROWN,
           ewt = EWT, lma = LMA, prot = PROT, cbc = CBC, alpha = alpha,
           check = check)
}


#' @rdname deprecated
#' @export
#'
PROSPECT_LUT <- function(Input_PROSPECT, SpecPROSPECT = NULL) {

  .Deprecated(old = "PROSPECT_LUT",
              new = "prospect_lut")
  prospect_lut(input_prospect = Input_PROSPECT,
               spec_prospect = SpecPROSPECT)
}


#' @rdname deprecated
#' @export
#'
reshape_lop4inversion <- function(Refl, Tran, SpecPROSPECT){

  .Deprecated(old = "reshape_lop4inversion",
              new = "reshape_lop_for_inversion")
  reshape_lop_for_inversion(refl = Refl, tran = Tran,
                            spec_prospect = SpecPROSPECT)
}


#' @rdname deprecated
#' @export
#'
SetInitParm <- function(Parms2Estimate, ParmEst, PROSPECT_version, Nprior,
                        ANTinit, OptDomain, Refl, Tran){

  .Deprecated(old = "SetInitParm", new = "set_init_parm")
  set_init_parm(parms_to_estimate = Parms2Estimate, parm_estimated = ParmEst,
                prospect_version = PROSPECT_version, Nprior = Nprior,
                ant_init = ANTinit, optimal_domain = OptDomain,
                refl = Refl, tran = Tran)
}


#' @rdname deprecated
#' @export
#'
SetNValues <- function(Refl, Tran, SpecPROSPECT){

  .Deprecated(old = "SetNValues", new = "set_n_values")
  set_n_values(refl = Refl, tran = Tran, spec_prospect = SpecPROSPECT)
}


#' @rdname deprecated
#' @export
#'
tryInversion <- function(x0, MeritFunction, SpecPROSPECT, Refl, Tran,
                         Parms2Estimate, lb, ub, verbose = FALSE) {

  .Deprecated(old = "tryInversion",
              new = "try_inversion")
  try_inversion(x0 = x0, merit_function = MeritFunction,
                spec_prospect = SpecPROSPECT, refl = Refl, tran = Tran,
                parms_to_estimate = Parms2Estimate, lb = lb, ub = ub,
                verbose = verbose)
}
