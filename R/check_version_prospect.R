#' This function checks if the input parameters are defined as expected
#' to run either PROSPECT-D or PROSPECT-PRO
#' @param lma numeric. content corresponding to lma
#' @param prot numeric. content corresponding to protein content
#' @param cbc numeric. content corresponding to carbon based constituents
#'
#' @return list. updated lma, prot and cbc
#' @export

check_version_prospect <- function(lma, prot, cbc){
  # PROSPECT-D as default value
  if (is.null(lma) & prot == 0 & cbc == 0)
    lma <- 0.008
  # PROSPECT-PRO if prot or cbc are not NULL
  if (is.null(lma) & (prot > 0 | cbc > 0))
    lma <- 0
  # if calling PROSPECT-PRO (protein content or cbc defined by user)
  # then set lma to 0 in any case
  if (!lma==0 & (prot > 0 | cbc > 0)) {
    print_msg('version_PROSPECT')
    lma <- 0
  }
  return(list('lma' = lma, 'prot' = prot, 'cbc' = cbc))
}
