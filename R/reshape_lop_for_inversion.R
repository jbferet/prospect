#' Function to reshape reflectance and transmittance, and convert into dataframes
#'
#' @param refl  numeric. measured reflectance data
#' @param tran  numeric. measured transmittance data
#' @param spec_prospect list. Includes optical constants
#' refractive index, specific absorption coefficients and corresponding spectral bands
#'
#' @return RT list of leaf optics converted into dataframe and coresponding number of samples
#' @export
#'
reshape_lop_for_inversion <- function(refl, tran, spec_prospect){
  # check if refl/tran is data frame, matrix or vector
  # convert into data frame if not
  RT <- list('refl' = refl, 'tran' = tran)
  for (lop in names(RT)){
    if (!is.null(RT[[lop]])){
      classLOP <- class(RT[[lop]])[1]
      if (!classLOP %in% c('numeric', 'matrix', 'data.frame', 'data.table'))
        print_msg(cause = 'classLOP')
      if (classLOP %in% c('numeric', 'matrix'))
        RT[[lop]] <- data.frame(RT[[lop]])
      nbwl <- nrow(RT[[lop]])
      # check if the number of spectral bands is compatible with spec_prospect
      if (nbwl==length(spec_prospect$lambda))
        rownames(RT[[lop]]) <- spec_prospect$lambda
      if (!nbwl==length(spec_prospect$lambda))
        print_msg(cause = 'lengthLOP')
    }
  }
  # check if refl & tran have the same number of samples
  if (!is.null(RT$refl) & !is.null(RT$tran)){
    if (!ncol(RT$refl)==ncol(RT$tran))
      print_msg(cause = 'diff_RTshape')
  }
  RT$nb_samples <- unique(c(ncol(RT$refl), ncol(RT$tran)))
  return(RT)
}



#' @rdname reshape_lop4inversion-deprecated
#' @export
#'
reshape_lop4inversion <- function(Refl, Tran, SpecPROSPECT){

  .Deprecated("reshape_lop4inversion")
  reshape_lop_for_inversion(refl = Refl, tran = Tran,
                            spec_prospect = SpecPROSPECT)
}

