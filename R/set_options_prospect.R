#' set options
#'
#' @param fun character. name of the function which has optional parameters
#' @param options list. including
#' - nb_clusters numeric. number of clusters
#'
#' @return options with default values when missing
#' @export

set_options_prospect <- function(fun, options = NULL){

  if (fun %in% c('invert_prospect',
                 'invert_prospect_opt',
                 'invert_prospect_subdomain',
                 'optimal_features_sfs')){
    if (is.null(options$init_values))
      options$init_values <- data.frame(chl = 40, car = 10, ant = 0.1,
                                        brown = 0.0, ewt = 0.01, lma = 0.01,
                                        prot = 0.001, cbc = 0.009,
                                        n_struct = 1.5, alpha = 40)
    if (is.null(options$spec_prospect))
      options$spec_prospect <- prospect::spec_prospect_full_range
    if (is.null(options$merit_function))
      options$merit_function <- "merit_prospect_rmse"
    if (is.null(options$xlub))
      options$xlub <- data.frame(chl = c(1e-4, 150), car = c(1e-4, 25),
                                 ant = c(0, 50), brown = c(0, 4),
                                 ewt = c(1e-8, 0.1), lma = c(1e-8, 0.06),
                                 prot = c(1e-7, 0.006), cbc = c(1e-6, 0.054),
                                 n_struct = c(0.5, 4), alpha = c(10, 90))
    if (is.null(options$estimate_brown_pigments))
      options$estimate_brown_pigments <- FALSE
    if (is.null(options$estimate_alpha))
      options$estimate_alpha <- FALSE
    if (is.null(options$verbose))
      options$verbose <- FALSE
    if (is.null(options$progressBar))
      options$progressBar <- TRUE
  }
  return(options)
}
