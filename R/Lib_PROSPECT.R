# ==============================================================================
# prospect
# Lib_PROSPECT.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de Boissieu <fdeboiss@gmail.com>
# Copyright 2020/04 Jean-Baptiste FERET
# ==============================================================================
# This Library includes functions dedicated to PROSPECT simulation
# ==============================================================================

#' core function running PROSPECT
#' This function allows simulations using PROSPECT-D or PROSPECT-PRO depending
#' on the parameterization.
#  This code includes numerical optimizations proosed in the FLUSPECT code
#  Authors: Wout Verhoef, Christiaan van der Tol (tol@itc.nl), Joris Timmermans,
#  Date: 2007
#  Update from PROSPECT to FLUSPECT: January 2011 (CvdT)
#'
#' @param SpecPROSPECT list. Includes spectral constants derived from
#' SpecPROSPECT_FullRange: refractive index, specific absorption coefficients
#' and corresponding spectral bands
#' @param N numeric. Leaf structure parameter
#' @param CHL numeric. Chlorophyll content (microg.cm-2)
#' @param CAR numeric. Carotenoid content (microg.cm-2)
#' @param ANT numeric. Anthocyanin content (microg.cm-2)
#' @param BROWN numeric. Brown pigment content (Arbitrary units)
#' @param EWT numeric. Equivalent Water Thickness (g.cm-2)
#' @param LMA numeric. Leaf Mass per Area (g.cm-2)
#' @param PROT numeric. protein content  (g.cm-2)
#' @param CBC numeric. NonProt Carbon-based constituent content (g.cm-2)
#' @param alpha numeric. Solid angle for incident light at surface of leaf
#'
#' @return leaf directional-hemispherical reflectance and transmittance
#' @importFrom expint expint
#' @export
PROSPECT <- function(SpecPROSPECT = NULL, N = 1.5, CHL = 40.0,
                     CAR = 8.0, ANT = 0.0, BROWN = 0.0, EWT = 0.01,
                     LMA = NULL, PROT = 0, CBC = 0, alpha = 40.0) {

  # default: simulates leaf optics using full spectral range
  if (is.null(SpecPROSPECT)) SpecPROSPECT <- prospect::SpecPROSPECT_FullRange
  # check if LMA, PROT and CBC are correctly parameterized
  dm_val <- check_version_prospect(LMA, PROT, CBC)
  # compute total absorption corresponding to each homogeneous layer
  Kall <- (CHL * SpecPROSPECT$SAC_CHL +
             CAR * SpecPROSPECT$SAC_CAR +
             ANT * SpecPROSPECT$SAC_ANT +
             BROWN * SpecPROSPECT$SAC_BROWN +
             EWT * SpecPROSPECT$SAC_EWT +
             dm_val$LMA * SpecPROSPECT$SAC_LMA +
             dm_val$PROT * SpecPROSPECT$SAC_PROT +
             dm_val$CBC * SpecPROSPECT$SAC_CBC) / N

  # Non-conservative scattering (normal case) when Kall > 0
  j <- which(Kall <= 0)
  t1 <- (1 - Kall) * exp(-Kall)
  t2 <- (Kall * Kall) * expint(Kall)
  tau <- t1 + t2
  if (length(j)>0) tau[j] <- 1
  # tau <- matrix(1, ncol = 1, nrow = length(t1))
  # tau[j] <- t1[j] + t2[j]

  # ***********************************************************************
  # reflectance and transmittance of one layer
  # ***********************************************************************
  # Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969),
  # Interaction of isotropic ligth with a compact plant leaf, J. Opt.
  # Soc. Am., 59(10):1376-1379.
  # ***********************************************************************
  # reflectivity and transmissivity at the interface
  # ***********************************************************************
  if (alpha == 40) {
    talf <- SpecPROSPECT$calctav_40
  } else {
    talf <- calctav(alpha, SpecPROSPECT$nrefrac)
  }
  ralf <- 1 - talf
  t12 <- SpecPROSPECT$calctav_90
  r12 <- 1 - t12
  t21 <- t12 / (SpecPROSPECT$nrefrac**2)
  r21 <- 1 - t21

  # top surface side
  denom <- 1 - (r21 * r21 * (tau**2))
  Ta <- (talf * tau * t21) / denom
  Ra <- ralf + (r21 * tau * Ta)
  # bottom surface side
  t <- t12 * tau * t21 / denom
  r <- r12 + (r21 * tau * t)

  # ***********************************************************************
  # reflectance and transmittance of N layers
  # Stokes equations to compute properties of next N-1 layers (N real)
  # Normal case
  # ***********************************************************************
  # Stokes G.G. (1862), On the intensity of the light reflected from
  # or transmitted through a pile of plates, Proc. Roy. Soc. Lond.,
  # 11:545-556.
  # ***********************************************************************
  D <- sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
  rq <- r**2
  tq <- t**2
  a <- (1 + rq - tq + D) / (2 * r)
  b <- (1 - rq + tq + D) / (2 * t)

  bNm1 <- b**(N - 1)
  bN2 <- bNm1**2
  a2 <- a**2
  denom <- a2 * bN2 - 1
  Rsub <- a * (bN2 - 1) / denom
  Tsub <- bNm1 * (a2 - 1) / denom

  # 	Case of zero absorption
  j <- which(r + t >= 1)
  Tsub[j] <- t[j] / (t[j] + (1 - t[j]) * (N - 1))
  Rsub[j] <- 1 - Tsub[j]

  # leaf reflectance and transmittance : combine top layer with next N-1 layers
  denom <- 1 - Rsub * r
  tran <- Ta * Tsub / denom
  refl <- Ra + (Ta * Rsub * t) / denom
  return(data.frame('wvl' = SpecPROSPECT$lambda,
                    'Reflectance' = refl,
                    'Transmittance' = tran))
}

#' computation of transmissivity of a dielectric plane surface,
#' averaged over all directions of incidence and over all polarizations.
#'
#' @param alpha numeric. max incidence angle of solid angle of incident light
#' @param nr numeric. refractive index
#'
#' @return numeric. Transmissivity of a dielectric plane surface
#' @export
calctav <- function(alpha, nr) {
  # Stern F. (1964), Transmission of isotropic radiation across an
  # interface between two dielectrics, Appl. Opt., 3(1):111-113.
  # Allen W.A. (1973), Transmission of isotropic light across a
  # dielectric surface in two and three dimensions, J. Opt. Soc. Am.,
  # 63(6):664-666.
  # ***********************************************************************

  rd <- pi / 180
  n2 <- nr**2
  np <- n2 + 1
  nm <- n2 - 1
  a <- (nr + 1) * (nr + 1) / 2
  k <- -(n2 - 1) * (n2 - 1) / 4
  sa <- sin(alpha * rd)

  b2 <- (sa**2) - (np / 2)
  if (alpha == 90) {
    b1 <- 0 * b2
  } else {
    b1 <- sqrt((b2**2) + k)
  }
  b <- b1 - b2
  b3 <- b**3
  a3 <- a**3
  ts <- ((k**2) / (6 * b3) + (k / b) - b / 2) - ((k**2) / (6 * a3) + (k / a) - (a / 2))

  tp1 <- -2 * n2 * (b - a) / (np**2)
  tp2 <- -2 * n2 * np * log(b / a) / (nm**2)
  tp3 <- n2 * ((1 / b) - (1 / a)) / 2
  tp4 <- 16 * n2**2 * ((n2**2) + 1) * log(((2 * np * b) - (nm**2)) / ((2 * np * a) - (nm**2))) / ((np**3) * (nm**2))
  tp5 <- 16 * (n2**3) * (1 / ((2 * np * b) - (nm**2)) - (1 / (2 * np * a - (nm**2)))) / (np**3)
  tp <- tp1 + tp2 + tp3 + tp4 + tp5
  tav <- (ts + tp) / (2 * (sa**2))
  return(tav)
}

#' This function checks if the input parameters are defined as expected
#' to run either PROSPECT-D or PROSPECT-PRO
#' @param LMA numeric. content corresponding to LMA
#' @param PROT numeric. content corresponding to protein content
#' @param CBC numeric. content corresponding to carbon based constituents
#'
#' @return list. updated LMA, PROT and CBC
#' @export

check_version_prospect <- function(LMA, PROT, CBC){
  # PROSPECT-D as default value
  if (is.null(LMA) & PROT == 0 & CBC == 0) LMA <- 0.008
  # PROSPECT-PRO if PROT or CBC are not NULL
  if (is.null(LMA) & (PROT > 0 | CBC > 0)) LMA <- 0
  # if calling PROSPECT-PRO (protein content or CBC defined by user)
  # then set LMA to 0 in any case
  if (!LMA==0 & (PROT > 0 | CBC > 0)) {
    print_msg('version_PROSPECT')
    LMA <- 0
  }
  return(list('LMA' = LMA, 'PROT' = PROT, 'CBC' = CBC))
}

#' This function adapts SpecPROSPECT accordingly to experimental data
#' or to a spectral domain defined by UserDomain
#' @param lambda numeric. Spectral bands corresponding to experimental data
#' @param SpecPROSPECT list. Includes optical constants: refractive index,
#' specific absorption coefficients and corresponding spectral bands
#' @param Refl numeric. Measured reflectance data
#' @param Tran numeric. Measured Transmittance data
#' @param UserDomain numeric. either Lower and upper bounds for domain of
#' interest (optional) or list of spectral bands of interest
#' @param UL_Bounds boolean. set to TRUE if UserDomain only includes lower and
#' upper band, set to FALSE if UserDomain is a list of spectral bands (in nm)
#'
#' @return list including spectral properties at the new resolution
#' @import dplyr
#' @importFrom utils tail
#' @export

FitSpectralData <- function(lambda, SpecPROSPECT = NULL,
                            Refl = NULL, Tran = NULL,
                            UserDomain = NULL, UL_Bounds = FALSE) {
  # default: simulates leaf optics using full spectral range
  if (is.null(SpecPROSPECT)) SpecPROSPECT <- prospect::SpecPROSPECT_FullRange
  # convert Refl and Tran into dataframe if needed
  if (class(Refl)[1]%in%c('numeric', 'matrix')) Refl <- data.frame(Refl)
  if (class(Tran)[1]%in%c('numeric', 'matrix')) Tran <- data.frame(Tran)
  # Adjust LOP: check common spectral domain between PROSPECT and leaf optics
  if (!is.null(Refl)) Refl <- Refl %>% filter(lambda%in%SpecPROSPECT$lambda)
  if (!is.null(Tran)) Tran <- Tran %>% filter(lambda%in%SpecPROSPECT$lambda)
  lambda <- lambda[lambda%in%SpecPROSPECT$lambda]
  # Adjust PROSPECT
  lb <- lambda
  SpecPROSPECT  <- SpecPROSPECT %>% filter(SpecPROSPECT$lambda%in%lb)
  # if UserDomain is defined
  if (is.null(UserDomain)) UserDomain <- lambda
  if (UL_Bounds==TRUE) UserDomain <- seq(min(UserDomain), max(UserDomain))
  if (!is.null(Refl)) Refl <- Refl %>% filter(lambda%in%UserDomain)
  if (!is.null(Tran)) Tran <- Tran %>% filter(lambda%in%UserDomain)
  lambda <- lambda[lambda%in%UserDomain]
  # Adjust PROSPECT
  SpecPROSPECT  <- SpecPROSPECT %>% filter(SpecPROSPECT$lambda%in%UserDomain)
  if (any(!UserDomain%in%lambda)){
    message('leaf optics out of range defined by UserDomain')
  }
  RT <- reshape_lop4inversion(Refl = Refl,
                              Tran = Tran,
                              SpecPROSPECT = SpecPROSPECT)
  return(list("SpecPROSPECT" = SpecPROSPECT, "lambda" = lambda,
              "Refl" = RT$Refl, "Tran" = RT$Tran, "nbSamples" = RT$nbSamples))
}

#' computation of a LUT of leaf optical properties using a set of
#' leaf chemical & structural parameters
#'
#' @param Input_PROSPECT dataframe. list of PROSPECT input parameters.
#' @param SpecPROSPECT list. spectral constants
#' refractive index, specific absorption coefficients & spectral bands
#'
#' @return list. LUT including leaf reflectance and transmittance
#' @export
PROSPECT_LUT <- function(Input_PROSPECT, SpecPROSPECT = NULL) {

  if (is.null(SpecPROSPECT)) SpecPROSPECT <- prospect::SpecPROSPECT_FullRange
  # expected PROSPECT input parameters
  ExpectedParms <- data.frame('CHL' = 0, 'CAR' = 0, 'ANT' = 0, 'BROWN' = 0,
                              'EWT' = 0, 'LMA' = 0, 'PROT' = 0, 'CBC' = 0,
                              'N' = 1.5, 'alpha' = 40)
  # parameters provided
  inOK <- names(Input_PROSPECT)
  # identify missing elements
  Parm2Add <- which(!names(ExpectedParms) %in% inOK)
  # check if all parameters are included.
  if (length(Parm2Add) > 0) print_msg(cause = 'Missing_Input',
                                      args = list('Input' = inOK,
                                                  'Expected' = ExpectedParms))

  # re-order missing elements the end of the list using default value
  Input_PROSPECT <- Complete_Input_PROSPECT(Input_PROSPECT = Input_PROSPECT,
                                            Parm2Add = Parm2Add,
                                            ExpectedParms = ExpectedParms)
  # print number of samples to be simulated
  nbSamples <- nrow(Input_PROSPECT)
  messageLUT <- paste('A LUT with', nbSamples, 'samples will be produced')
  cat(colour_to_ansi('green'), messageLUT, "\033[0m\n")

  # run PROSPECT for nbSamples
  run_list_PROSPECT <- function(Input_PROSPECT, SpecPROSPECT){
    LUT_tmp <- do.call(PROSPECT,
                       c(list(SpecPROSPECT = SpecPROSPECT),
                         Input_PROSPECT))
    return(LUT_tmp)
  }
  indiv_leaves <- split(Input_PROSPECT,
                        factor(seq(1:nbSamples)))
  LUT_tmp <- lapply(X = indiv_leaves,
                    FUN = run_list_PROSPECT,
                    SpecPROSPECT = SpecPROSPECT)
  LUT_Refl <- as.data.frame(lapply(LUT_tmp,'[[', 'Reflectance'))
  LUT_Tran <- as.data.frame(lapply(LUT_tmp,'[[', 'Transmittance'))
  names(LUT_Refl) <- names(LUT_Tran) <- paste0('sample_',seq(1,nbSamples))
  return(list('Reflectance' = LUT_Refl,
              'Transmittance' = LUT_Tran,
              'Input_PROSPECT' = Input_PROSPECT))
}

#' Complete the list of PROSPECT parameters with default values
#'
#' @param urldb character. URL for online repository where to download data
#' @param dbName character. name of the database available online
#'
#' @return list. Includes leaf chemistry, refl, tran & number of samples
#' @importFrom data.table fread
#' @export

download_LeafDB <- function(urldb = NULL,
                            dbName = 'ANGERS'){
  # repository where data are stored
  if (is.null(urldb)) urldb <- 'https://gitlab.com/jbferet/myshareddata/raw/master/LOP/'
  # download leaf chemistry and optical properties
  DataBioch <- data.table::fread(file.path(urldb,dbName,'DataBioch.txt'))
  Refl <- data.table::fread(file.path(urldb,dbName,'ReflectanceData.txt'))
  Tran <- data.table::fread(file.path(urldb,dbName,'TransmittanceData.txt'))
  # Get wavelengths corresponding to the reflectance & transmittance measurements
  lambda <- Refl$wavelength
  Refl$wavelength <- Tran$wavelength <- NULL
  # Get the number of samples
  nbSamples <- ncol(Refl)
  return(list('DataBioch' = DataBioch,
              'lambda' = lambda,
              'Refl' = Refl,
              'Tran' = Tran,
              'nbSamples' = nbSamples))
}


#' Complete the list of PROSPECT parameters with default values
#'
#' @param Input_PROSPECT input parameters sent to PROSPECT by user
#' @param Parm2Add Parameters to be added to input parameters
#' @param ExpectedParms full set of parameters expected to run PROSPECT
#'
#' @return Input_PROSPECT
#' @export

Complete_Input_PROSPECT <- function(Input_PROSPECT, Parm2Add, ExpectedParms) {
  ii <- 0
  nbSamples <- length(Input_PROSPECT[[1]])
  nbInputs <- length(Input_PROSPECT)
  for (i in Parm2Add) {
    ii <- ii + 1
    nbInputs <- nbInputs + 1
    Input_PROSPECT[[nbInputs]] <- matrix(ExpectedParms[[i]],
                                         ncol = 1, nrow = nbSamples)
    names(Input_PROSPECT)[[nbInputs]] <- names(ExpectedParms)[[i]]
  }
  return(data.frame('CHL' = matrix(Input_PROSPECT$CHL, ncol = 1),
                    'CAR' = matrix(Input_PROSPECT$CAR, ncol = 1),
                    'ANT' = matrix(Input_PROSPECT$ANT, ncol = 1),
                    'BROWN' = matrix(Input_PROSPECT$BROWN, ncol = 1),
                    'EWT' = matrix(Input_PROSPECT$EWT, ncol = 1),
                    'LMA' = matrix(Input_PROSPECT$LMA, ncol = 1),
                    'PROT' = matrix(Input_PROSPECT$PROT, ncol = 1),
                    'CBC' = matrix(Input_PROSPECT$CBC, ncol = 1),
                    'N' = matrix(Input_PROSPECT$N, ncol = 1),
                    'alpha' = matrix(Input_PROSPECT$alpha, ncol = 1)))
}


#' Convert plain text colour to ANSI code
#'
#' @param colour colour in plain text ("red", "green", etc.) to convert to ANSI
#'
#' @return string representing provided colour as ANSI encoding
#'
#' @examples
#' colour_to_ansi("red") # gives: "\033[31m"
#' @export

colour_to_ansi <- function(colour) {
  # Note ANSI colour codes
  colour_codes <- list("black" = 30,
                       "red" = 31,
                       "green" = 32,
                       "yellow" = 33,
                       "blue" = 34,
                       "magenta" = 35,
                       "cyan" = 36,
                       "white" = 37)

  # Create ANSI version of colour
  ansi_colour <- paste0("\033[", colour_codes[[colour]], "m")
  return(ansi_colour)
}

