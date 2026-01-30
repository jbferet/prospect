#' Function to check good agreement between prospect version and parameter list
#'
#' @param prospect_version character. version used for inversion: 'D' or 'PRO'
#' @param estimate_brown_pigments boolean. should brown pigments be accounted for during inversion?
#' @param estimate_alpha boolean. should alpha be accounted for during inversion?
#' @param parms_to_estimate  character vector. Parameters to estimate (can be 'ALL')
#' @param xlub data.frame. Boundaries of the parameters to estimate.
#' The data.frame must have columns corresponding to \code{parms_to_estimate} first
#' line being the lower boundaries and second line the upper boundaries.
#' @param init_values  dataframe. Default values of PROSPECT parameters. During
#' optimization, they are used either as initialization values for parameters to
#' estimate, or as fix values for other parameters.
#' Parameters not compatible with prospect_version are not taken into account.
#'
#' @return list of parameters to estimate & corresponding lower/upper boundaries
#' @export

check_prospect_parms <- function(prospect_version, parms_to_estimate,
                                 estimate_brown_pigments, estimate_alpha, xlub,
                                 init_values){

  # check if version required is available
  if (!prospect_version %in% c('D', 'PRO'))
    print_msg(cause = 'WrongVersion')
  # define parameters if PROSPECT-D or PROSPECT-PRO
  if (prospect_version == 'D')
    all_parms <- c('n_struct', 'chl', 'car', 'ant', 'ewt', 'lma')
  if (prospect_version == 'PRO')
    all_parms <- c('n_struct', 'chl', 'car', 'ant', 'ewt', 'prot', 'cbc')
  # add brown pigments and alpha angle if required
  if (estimate_brown_pigments==TRUE)
    all_parms <- c(all_parms, "brown")
  if (estimate_alpha==TRUE)
    all_parms <- c(all_parms, "alpha")

  # add default values to xlub in case they were not defined
  xlub_default <- data.frame(chl = c(1e-4, 150), car = c(1e-4, 25),
                             ant = c(0, 50), brown = c(0, 4),
                             ewt = c(1e-8, 0.1), lma = c(1e-6, .06),
                             prot = c(1e-7, .006), cbc = c(1e-6, .054),
                             n_struct = c(.5, 4), alpha = c(10, 90))
  added_parm_lb <- setdiff(names(xlub_default),names(xlub))
  for (ad in added_parm_lb)
    xlub[[ad]] <- xlub_default[[ad]]

  init_values_default = data.frame(chl = 40, car = 10, ant = 0.1, brown = 0.0,
                                   ewt = 0.01, lma = 0.01, prot = 0.001,
                                   cbc = 0.009, n_struct = 1.5, alpha = 40)
  added_param_init <- setdiff(names(init_values_default),names(init_values))
  for (ad in added_param_init)
    init_values[[ad]] <- init_values_default[[ad]]

  # if 'ALL' is provided, then assess all parameters available
  if ("ALL" %in% parms_to_estimate)
    parms_to_estimate <- all_parms
  # if unknown parameter is provided, then warn
  if (!all(parms_to_estimate %in% all_parms))
    print_msg(cause = 'unknown_parm',
              args = list('parms_to_estimate' = parms_to_estimate,
                          'all_parms'= all_parms))
  parms_to_estimate <- all_parms[all_parms %in% parms_to_estimate]
  # if missing value in xlub defining lower and upper boudaries
  if (!all(all_parms %in% names(xlub)))
    print_msg(cause = 'missing_xlub')
  # if missing value in init_values defining initial values for inversion
  if (!all(all_parms %in% names(init_values)))
    print_msg(cause = 'missing_InitValues')
  init_values <- init_values[all_parms[all_parms %in% names(init_values)]]
  if (prospect_version == "PRO"){
    init_values$lma <- 0
    if (is.null(init_values$prot))
      init_values$prot <- 0.001
    if (is.null(init_values$cbc))
      init_values$cbc <- 0.009
  }
  if (prospect_version == "D"){
    init_values$prot <- init_values$cbc <- 0
    if (is.null(init_values$lma))
      init_values$lma <- 0.01
  }
  if (estimate_brown_pigments==FALSE)
    init_values$brown <- 0
  xlub <- data.frame(xlub[parms_to_estimate])
  lb <- xlub[1,]
  ub <- xlub[2,]
  return(list('lb' = lb, 'ub' = ub, 'parms_to_estimate' = parms_to_estimate,
              'init_values' = init_values))
}

