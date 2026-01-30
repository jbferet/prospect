#' Function handling error during inversion
#'
#' @param x0 numeric. Vector of input variables to estimate
#' @param merit_function  character. name of the merit function
#' @param spec_prospect list. Includes optical constants refractive index,
#' specific absorption coefficients and corresponding spectral bands
#' @param refl  numeric. measured reflectance data
#' @param tran  numeric. measured transmittance data
#' @param parms_to_estimate  character vector. Parameters to estimate
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

try_inversion <- function(x0, merit_function, spec_prospect, refl, tran,
                          parms_to_estimate, lb, ub, verbose = FALSE) {

  res <-list('par' = NA*vector(length = length(parms_to_estimate)))
  TolRange <- seq(-14,-2,1)
  for (i in TolRange){
    Tolerance <- 10**(i)
    if (is.na(res$par[1])){
      res <- tryCatch({
        res <- fmincon(
          x0 = as.numeric(x0[parms_to_estimate]),
          fn = merit_function, gr = NULL,
          spec_prospect = spec_prospect, refl = refl, tran = tran,
          input_prospect = x0, parms_to_estimate = parms_to_estimate,
          method = "SQP", A = NULL, b = NULL, Aeq = NULL, beq = NULL,
          lb = as.numeric(lb), ub = as.numeric(ub), hin = NULL, heq = NULL,
          tol = Tolerance, maxfeval = 2000, maxiter = 1000)
      }, error = function(cond) {
        if (verbose)
          print_msg(cause = 'AdjustTol',
                    args = list('Tolerance' = Tolerance))
        return(list('par' = NA*vector(length = length(parms_to_estimate))))
      },
      finally = {}
      )
    }
  }

  # test if one of the parameters to be estimated reached lower or upper bound
  attempt <- 0
  names(res$par) <- parms_to_estimate
  reinit <- FALSE
  for (parm in parms_to_estimate){
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
    for (parm in parms_to_estimate)
      x0[parm] <- lb[[parm]]+runif(1)*(ub[[parm]]-lb[[parm]])
    reinit <- TRUE
  }
  # perform inversion with readjusted initial values if lower/upper band reached
  while (reinit==TRUE & attempt<2){
    attempt <- attempt+1
    if (verbose)
      print_msg(cause = 'ULBounds')
    res <- list('par' = NA*vector(length = length(parms_to_estimate)))
    for (i in TolRange){
      Tolerance <- 10**(i)
      if (is.na(res$par[1])){
        res <- tryCatch({
            res <- fmincon(
              x0 = as.numeric(x0[parms_to_estimate]),
              fn = merit_function,
              gr =NULL,
              spec_prospect = spec_prospect,
              refl = refl, tran = tran,
              input_prospect = x0,
              parms_to_estimate = parms_to_estimate,
              method = "SQP", A = NULL, b = NULL, Aeq = NULL, beq = NULL,
              lb = as.numeric(lb), ub = as.numeric(ub), hin = NULL, heq = NULL,
              tol = Tolerance, maxfeval = 2000, maxiter = 1000)
          }, error = function(cond) {
            if (verbose)
              print_msg(cause = 'AdjustTol',
                        args = list('Tolerance' = Tolerance))
            return(list('par' = NA*vector(length = length(parms_to_estimate))))
          },
          finally = {}
        )
      }
    }
    names(res$par) <- parms_to_estimate
    reinit <- FALSE
    for (parm in parms_to_estimate){
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



#' @rdname tryInversion-deprecated
#' @export
#'
tryInversion <- function(x0, MeritFunction, SpecPROSPECT, Refl, Tran,
                         Parms2Estimate, lb, ub, verbose = FALSE) {

  .Deprecated("tryInversion")
  try_inversion(x0 = x0, merit_function = MeritFunction,
                spec_prospect = SpecPROSPECT, refl = Refl, tran = Tran,
                parms_to_estimate = Parms2Estimate, lb = lb, ub = ub,
                verbose = verbose)
}
