#' Complete the list of PROSPECT parameters with default values
#'
#' @param input_prospect input parameters sent to PROSPECT by user
#' @param parm_to_add Parameters to be added to input parameters
#' @param expected_parms full set of parameters expected to run PROSPECT
#'
#' @return input_prospect
#' @export
#'
complete_input_prospect <- function(input_prospect, parm_to_add, expected_parms) {
  ii <- 0
  nb_samples <- length(input_prospect[[1]])
  nb_inputs <- length(input_prospect)
  for (i in parm_to_add) {
    ii <- ii + 1
    nb_inputs <- nb_inputs + 1
    input_prospect[[nb_inputs]] <- matrix(expected_parms[[i]],
                                         ncol = 1, nrow = nb_samples)
    names(input_prospect)[[nb_inputs]] <- names(expected_parms)[[i]]
  }
  return(data.frame('chl' = matrix(input_prospect$chl, ncol = 1),
                    'car' = matrix(input_prospect$car, ncol = 1),
                    'ant' = matrix(input_prospect$ant, ncol = 1),
                    'brown' = matrix(input_prospect$brown, ncol = 1),
                    'ewt' = matrix(input_prospect$ewt, ncol = 1),
                    'lma' = matrix(input_prospect$lma, ncol = 1),
                    'prot' = matrix(input_prospect$prot, ncol = 1),
                    'cbc' = matrix(input_prospect$cbc, ncol = 1),
                    'n_struct' = matrix(input_prospect$n_struct, ncol = 1),
                    'alpha' = matrix(input_prospect$alpha, ncol = 1)))
}
