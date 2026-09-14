#' Performs PROSPECT inversion using iterative optimization
#' in order to define estimate a set of user defined parameters
#' @param refl  numeric. Measured reflectance data
#' @param tran  numeric. Measured Transmittance data
#' @param parms_to_estimate  character. Parameters to estimate (can be 'ALL')
#' @param prospect_version  character. prospect version used for inversion: 'D' or 'PRO'
#' See details.
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
#' @return output_prospect estimated values corresponding to parms_to_estimate
#' @importFrom progress progress_bar
#' @details
#' Two versions of prospect are available for inversion.
#' The version is depending on the parameters taken into account:
#'
#' | Version  | D                                      | PRO
#' | :------: |:--------------------------------------:|:--------------------------------------:|
#' | chl      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | car      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | ant      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | brown    |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | ewt      |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#' | lma      |`r emojifont::emoji('white_check_mark')`|                                        |
#' | prot     |                                        |`r emojifont::emoji('white_check_mark')`|
#' | cbc      |                                        |`r emojifont::emoji('white_check_mark')`|
#' | n_struct |`r emojifont::emoji('white_check_mark')`|`r emojifont::emoji('white_check_mark')`|
#'
#' Argument `init_values` is expecting a default value for each of the parameters as well as an `alpha` value.
#' @importFrom pracma fmincon
#' @export
#' @md
#'
invert_prospect <- function(refl = NULL, tran = NULL,
                            parms_to_estimate = 'ALL',
                            prospect_version = 'D', options = NULL) {

  options <- set_options_prospect(fun = 'invert_prospect', options = options)
  spec_prospect <- options$spec_prospect
  merit_function <- options$merit_function
  xlub <- options$xlub
  init_values <- options$init_values
  progressBar <- options$progressBar
  verbose <- options$verbose
  estimate_alpha <- options$estimate_alpha
  estimate_brown_pigments <- options$estimate_brown_pigments
  # check if list of parameters applicable to PROSPECT version
  parms_checked <- check_prospect_parms(prospect_version = prospect_version,
                                        parms_to_estimate = parms_to_estimate,
                                        estimate_brown_pigments = estimate_brown_pigments,
                                        estimate_alpha = estimate_alpha,
                                        xlub = xlub,
                                        init_values = init_values)
  # complement initial values
  parms_checked$init_values <- define_input_prospect(input_prospect = parms_checked$init_values)
  # check if data class is compatible and convert into data.frame
  RT <- reshape_lop_for_inversion(refl = refl, tran = tran,
                                  spec_prospect = spec_prospect)
  output_prospect <- list()
  if (progressBar==TRUE){
    pb <- progress::progress_bar$new(
      format = "inverting prospect [:bar] :percent in :elapsedfull, estimated time remaining :eta",
      total = RT$nb_samples, clear = FALSE, width= 100)
  }
  for (idsample in seq_len(RT$nb_samples)){
    res <- try_inversion(x0 = parms_checked$init_values,
                         merit_function = merit_function,
                         spec_prospect = spec_prospect,
                         refl = RT$refl[[idsample]], tran = RT$tran[[idsample]],
                         parms_to_estimate = parms_checked$parms_to_estimate,
                         lb = parms_checked$lb, ub = parms_checked$ub,
                         verbose = verbose)
    if (NA %in% res$par){
      modify_init <- match(parms_checked$parms_to_estimate,
                          names(parms_checked$init_values))
      update_init_values <- parms_checked$init_values
      update_init_values[modify_init] <- 1.1*update_init_values[modify_init]
      res <- try_inversion(x0 = update_init_values,
                          merit_function = merit_function,
                          spec_prospect = spec_prospect,
                          refl = RT$refl[[idsample]], tran =RT$Tran[[idsample]],
                          parms_to_estimate = parms_checked$parms_to_estimate,
                          lb = parms_checked$lb, ub = parms_checked$ub,
                          verbose = verbose)
    }
    names(res$par) <- parms_checked$parms_to_estimate
    output_prospect[[idsample]] <- parms_checked$init_values
    output_prospect[[idsample]][names(res$par)] <- res$par
    if (progressBar==TRUE)
      pb$tick()
  }
  output_prospect <- do.call(rbind,output_prospect)
  return(output_prospect)
}
