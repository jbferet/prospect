#' Function printing message
#' - if leaf optics do not match with spectral domain
#' - if leaf optics do not match with expected data class
#' - if confusion is identified for prospect version to use
#' - if missing input variables during PROSPECT_LUT simulations
#' @param cause character. cause raising message print
#' @param args list. optional arguments used to print message
#'
#' @return none
#' @export

print_msg <- function(cause, args = NULL){
  if (cause == 'lengthLOP'){
    message('MISMATCH between leaf optical properties and spectral domain')
    message('Make sure leaf optical properties match with the spectral domain')
    message('defined in spec_prospect$lambda')
    message('if not, use the function "fit_spectral_data"')
    stop()
  }

  if (cause == 'classLOP'){
    message('the class of leaf optical properties is not identified')
    message('Make sure leaf optical properties are either data.frame, or matrix, or numeric vector')
    stop()
  }

  if (cause == 'version_PROSPECT'){
    message("prot and/or cbc are not set to 0")
    message("lma is not set to 0 neither, which is physically incorrect")
    message("(lma = prot + cbc)")
    message("We assume that PROSPECT-PRO was called and set lma to 0")
    message("Please correct input parameters lma, prot and/or cbc if needed")
  }

  if (cause == 'version_PROSPECT_Invert_PROSPECT_OPT'){
    message('lma is not estimated directly using PROSPECT-PRO')
    message('Please run inversion for PROSPECT-PRO and sum prot and cbc')
    message('if you want to get estimated lma from PROSPECT-PRO')
  }

  if (cause == 'Missing_Input'){
    message('_________________ WARNING _________________')
    message('  The list of parameters provided for LUT  ')
    message('      in the input_prospect variable       ')
    message('does not include the full set of parameters')
    message('')
    message(' The following parameters will be set to  ')
    message('   their default value as given here      ')
    WhichMissing <- args$Expected[which(!names(args$Expected) %in% args$Input)]
    print(WhichMissing)
  }

  if (cause == 'diff_RTshape'){
    message('Reflectance and Transmittance do not have the same number of samples')
    message('please fix that')
    stop()
  }

  if (cause == 'WrongVersion')
    stop('PROSPECT_version not available. Choice is limited to "D" and "PRO".')

  if (cause == 'NoBrown_OPT'){
    message('brown pigments are not accounted for when performing optimal estimation of leaf chemistry')
    message('model version excluding brown pigments will be used instead')
  }
  if (cause == 'AdjustTol')
    message('Adjusting Tolerance value for iterative optimization:  ',args$Tolerance)

  if (cause == 'ULBounds'){
    message('lower or upper bound reached for one or several parameters to estimate')
    message('re-run inversion with updated initial values to attempt proper convergence')
  }

  if (cause == 'Nprior_noRefl'){
    warning("________________________ WARNING _______________________")
    warning("The optimal prior estimation of N using Reflectance only")
    warning("requires information at 1131nm.")
    warning("The reflectance does not include this spectral band")
    warning("Using reflectance at 800 nm instead")
  }

  if (cause == 'Nprior_noRefl_atAll'){
    warning("________________________ WARNING _______________________")
    warning("The spectral information of the reflectance provided here")
    warning("does not contain the spectral bands required to estimate")
    warning("prior information about n_struct.")
    warning("The process will stop")
    stop()
  }

  if (cause == 'Nprior_noTran'){
    warning("________________________ WARNING _______________________")
    warning("The optimal prior estimation of N using Transmittance only")
    warning("requires information at 1121nm.")
    warning("The Transmittance does not include this spectral band")
    warning("Using Transmittance at 753 nm instead")
  }

  if (cause == 'Nprior_noTran_atAll'){
    warning("________________________ WARNING _______________________")
    warning("The spectral information of the Transmittance provided here")
    warning("does not contain the spectral bands required to estimate")
    warning("prior information about N.")
    warning("The proces will stop")
    stop()
  }

  if (cause == 'NoOpt_ANT'){
    message('Currently no optimal estimation for anthocyanins')
    message('PROSPECT inversion will be performed using full spectral information')
  }

  if (cause == 'missing_InitValues')
    stop('missing prospect parameters in "init_values"')

  if (cause == 'missing_xlub')
    stop('missing prospect parameters in "xlub"')

  if (cause == 'unknown_parm'){
    message('_________________ WARNING _________________')
    message(' unknown parameter requested for inversion ')
    print(args$parms_to_estimate[which(!args$parms_to_estimate %in% args$all_parms)])
    message('will be ignored')
  }
  return(invisible())
}
