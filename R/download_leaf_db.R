#' Complete the list of PROSPECT parameters with default values
#'
#' @param url_db character. URL for online repository where to download data
#' @param db_name character. name of the database available online
#'
#' @return list. Includes leaf chemistry, refl, tran & number of samples
#' @export
#'
download_leaf_db <- function(url_db = NULL, db_name = 'ANGERS'){
  # repository where data are stored
  if (is.null(url_db))
    url_db <- 'https://gitlab.com/jbferet/myshareddata/raw/master/LOP/'
  # download leaf chemistry and optical properties
  leaf_chemistry <- data.table::fread(file.path(url_db,db_name,'DataBioch.txt'))
  refl <- data.table::fread(file.path(url_db,db_name,'ReflectanceData.txt'))
  tran <- data.table::fread(file.path(url_db,db_name,'TransmittanceData.txt'))
  # Get wavelengths corresponding to the reflectance & transmittance measurements
  lambda <- refl$wavelength
  refl$wavelength <- tran$wavelength <- NULL
  # Get the number of samples
  nb_samples <- ncol(refl)
  return(list('leaf_chemistry' = leaf_chemistry, 'lambda' = lambda,
              'refl' = refl, 'tran' = tran, 'nb_samples' = nb_samples))
}
