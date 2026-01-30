#' computation of a LUT of leaf optical properties using a set of
#' leaf chemical & structural parameters
#'
#' @param input_prospect dataframe. list of PROSPECT input parameters.
#' @param spec_prospect  list. spectral constants
#' refractive index, specific absorption coefficients & spectral bands
#'
#' @return list. LUT including leaf reflectance and transmittance
#' @export
#'
prospect_lut <- function(input_prospect, spec_prospect = NULL) {

  if (is.null(spec_prospect))
    spec_prospect <- prospect::spec_prospect_full_range
  # expected PROSPECT input parameters
  expected_parms <- data.frame('chl' = 0, 'car' = 0, 'ant' = 0, 'brown' = 0,
                               'ewt' = 0, 'lma' = 0, 'prot' = 0, 'cbc' = 0,
                               'n_prior' = 1.5, 'alpha' = 40)
  # parameters provided
  inOK <- names(input_prospect)
  # identify missing elements
  parm_to_add <- which(!names(expected_parms) %in% inOK)
  # check if all parameters are included.
  if (length(parm_to_add) > 0)
    print_msg(cause = 'Missing_Input',
              args = list('Input' = inOK, 'Expected' = expected_parms))

  # re-order missing elements the end of the list using default value
  input_prospect <- complete_input_prospect(input_prospect = input_prospect,
                                            parm_to_add = parm_to_add,
                                            expected_parms = expected_parms)
  # print number of samples to be simulated
  nb_samples <- nrow(input_prospect)
  messageLUT <- paste('A LUT with', nb_samples, 'samples will be produced')
  cat(colour_to_ansi('green'), messageLUT, "\033[0m\n")

  # run PROSPECT for nb_samples
  run_list_PROSPECT <- function(input_prospect, spec_prospect){
    lut_tmp <- do.call(prospect,
                       c(list(spec_prospect = spec_prospect), input_prospect))
    return(lut_tmp)
  }
  indiv_leaves <- split(input_prospect, factor(seq_len(nb_samples)))
  lut_tmp <- lapply(X = indiv_leaves,
                    FUN = run_list_PROSPECT,
                    spec_prospect = spec_prospect)
  lut_Refl <- as.data.frame(lapply(lut_tmp,'[[', 'reflectance'))
  lut_Tran <- as.data.frame(lapply(lut_tmp,'[[', 'transmittance'))
  names(lut_Refl) <- names(lut_Tran) <- paste0('sample_',seq_len(nb_samples))
  return(list('reflectance' = lut_Refl,
              'transmittance' = lut_Tran,
              'input_prospect' = input_prospect))
}


#' @rdname PROSPECT_LUT-deprecated
#' @export
#'
PROSPECT_LUT <- function(Input_PROSPECT, SpecPROSPECT = NULL) {

  .Deprecated("PROSPECT_LUT")
  prospect_lut(input_prospect = Input_PROSPECT,
               spec_prospect = SpecPROSPECT)
}
