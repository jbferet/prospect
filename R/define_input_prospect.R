#' This function produces a data frame from all prospect input variables if not
#' defined already
#'
#' @param input_prospect numeric. prospect input parameters
#' @param chl numeric. Chlorophyll content (microg.cm-2)
#' @param car numeric. carotenoid content (microg.cm-2)
#' @param ant numeric. anthocyanin content (microg.cm-2)
#' @param brown numeric. brown pigment content (Arbitrary units)
#' @param ewt numeric. Equivalent Water Thickness (g.cm-2)
#' @param lma numeric. Leaf Mass per Area (g.cm-2)
#' @param prot numeric. protein content  (g.cm-2)
#' @param cbc numeric. Nonprot carbon-based constituent content (g.cm-2)
#' @param n_struct numeric. Leaf structure parameter
#' @param alpha numeric. Solid angle for incident light at surface of leaf
#'
#' @return list. updated lma, prot and cbc
#' @export

define_input_prospect <- function(input_prospect, chl = NULL, car = NULL,
                                  ant = NULL, brown = NULL, ewt = NULL,
                                  lma = NULL, prot = NULL, cbc = NULL,
                                  n_struct = NULL, alpha = NULL){

  default_prospect <- data.frame('chl' = 40.0, 'car' = 8.0, 'ant' = 0.0,
                                 'brown' = 0.0, 'ewt' = 0.01, 'lma' = 0.0,
                                 'prot'= 0.0, 'cbc' = 0.0, 'n_struct' = 1.5,
                                 'alpha' = 40.0)
  if (is.null(input_prospect)){
    dm_val <- check_version_prospect(lma = lma, prot = prot, cbc = cbc)
    input_prospect <- data.frame('n_struct' = n_struct, 'chl' = chl, 'car' = car,
                                 'ant' = ant, 'brown' = brown, 'ewt' = ewt,
                                 'lma' = dm_val$lma, 'prot'= dm_val$prot,
                                 'cbc' = dm_val$cbc, 'alpha' = alpha)
  } else if (!is.null(input_prospect)){
    names(input_prospect) <- tolower(names(input_prospect))
    if (!is.null(input_prospect$n)){
      input_prospect$n_struct <- input_prospect$n
      input_prospect$n <- NULL
    }
    missing <- which(!names(default_prospect)%in%names(input_prospect))
    if (length(missing)>0)
      input_prospect <- cbind(input_prospect, default_prospect[missing])
    dm_val <- check_version_prospect(lma = input_prospect$lma,
                                     prot = input_prospect$prot,
                                     cbc = input_prospect$cbc)
    input_prospect$lma <- dm_val$lma
    input_prospect$prot <- dm_val$prot
    input_prospect$cbc <- dm_val$cbc
  }
  return(input_prospect)
}

#' @rdname define_Input_PROSPECT-deprecated
#' @export
#'
define_Input_PROSPECT <- function(Input_PROSPECT, CHL = NULL, CAR = NULL,
                                  ANT = NULL, BROWN = NULL, EWT = NULL,
                                  LMA = NULL, PROT = NULL, CBC = NULL,
                                  N = NULL, alpha = NULL){

  .Deprecated("define_Input_PROSPECT")
  define_input_prospect(input_prospect = Input_PROSPECT, chl = CHL, car = CAR,
                        ant = ANT, brown = BROWN, ewt = EWT, lma = LMA,
                        prot = PROT, cbc = CBC, n_struct = N, alpha = alpha)
}

