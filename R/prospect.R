#' core function running PROSPECT
#' This function allows simulations using PROSPECT-D or PROSPECT-PRO depending
#' on the parameterization.
#  This code includes numerical optimizations proosed in the FLUSPECT code
#  Authors: Wout Verhoef, Christiaan van der Tol (tol@itc.nl), Joris Timmermans,
#  Date: 2007
#  Update from PROSPECT to FLUSPECT: January 2011 (CvdT)
#'
#' @param spec_prospect list. Includes spectral constants derived from
#' spec_prospect_full_range: refractive index, specific absorption coefficients
#' and corresponding spectral bands
#' @param input_prospect list. Includes all prospect input parameters
#' @param n_struct numeric. Leaf structure parameter
#' @param chl numeric. Chlorophyll content (microg.cm-2)
#' @param car numeric. Carotenoid content (microg.cm-2)
#' @param ant numeric. Anthocyanin content (microg.cm-2)
#' @param brown numeric. Brown pigment content (Arbitrary units)
#' @param ewt numeric. Equivalent Water Thickness (g.cm-2)
#' @param lma numeric. Leaf Mass per Area (g.cm-2)
#' @param prot numeric. protein content  (g.cm-2)
#' @param cbc numeric. NonProt Carbon-based constituent content (g.cm-2)
#' @param alpha numeric. Solid angle for incident light at surface of leaf
#' @param check boolean. set to TRUE to check input data format
#'
#' @return leaf directional-hemispherical reflectance and transmittance
#' @importFrom expint expint
#' @export
#'
prospect <- function(spec_prospect = NULL, input_prospect = NULL,
                     n_struct = 1.5, chl = 40.0, car = 8.0, ant = 0.0, brown = 0.0,
                     ewt = 0.01, lma = NULL, prot = 0, cbc = 0, alpha = 40.0,
                     check = TRUE) {

  # define PROSPECT input in a dataframe
  input_prospect <- define_input_prospect(input_prospect = input_prospect,
                                          chl = chl, car = car, ant = ant,
                                          brown = brown, ewt = ewt, lma = lma,
                                          prot = prot, cbc = cbc,
                                          n_struct = n_struct, alpha)
  # if (check) input_prospect <- define_input_prospect(input_prospect, chl, car,
  #                                                    ant, brown, ewt, lma, prot,
  #                                                    cbc, n_struct, alpha)
  # default: simulates leaf optics using full spectral range
  if (is.null(spec_prospect))
    spec_prospect <- prospect::spec_prospect_full_range
  # compute total absorption corresponding to each homogeneous layer
  Kall <- (input_prospect$chl * spec_prospect$sac_chl +
             input_prospect$car * spec_prospect$sac_car +
             input_prospect$ant * spec_prospect$sac_ant +
             input_prospect$brown * spec_prospect$sac_brown +
             input_prospect$ewt * spec_prospect$sac_ewt +
             input_prospect$lma * spec_prospect$sac_lma +
             input_prospect$prot * spec_prospect$sac_prot +
             input_prospect$cbc * spec_prospect$sac_cbc) / input_prospect$n_struct

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
  if (input_prospect$alpha == 40) {
    talf <- spec_prospect$calctav_40
  } else {
    talf <- calctav(alpha = input_prospect$alpha, nr = spec_prospect$nrefrac)
  }
  ralf <- 1 - talf
  t12 <- spec_prospect$calctav_90
  r12 <- 1 - t12
  t21 <- t12 / (spec_prospect$nrefrac**2)
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

  bNm1 <- b**(input_prospect$n_struct - 1)
  bN2 <- bNm1**2
  a2 <- a**2
  denom <- a2 * bN2 - 1
  Rsub <- a * (bN2 - 1) / denom
  Tsub <- bNm1 * (a2 - 1) / denom

  # 	Case of zero absorption
  j <- which(r + t >= 1)
  Tsub[j] <- t[j] / (t[j] + (1 - t[j]) * (input_prospect$n_struct - 1))
  Rsub[j] <- 1 - Tsub[j]

  # leaf reflectance and transmittance : combine top layer with next N-1 layers
  denom <- 1 - Rsub * r
  tran <- Ta * Tsub / denom
  refl <- Ra + (Ta * Rsub * t) / denom
  return(data.frame('wvl' = spec_prospect$lambda,
                    'reflectance' = refl,
                    'transmittance' = tran))
}

#' @rdname prospect-deprecated
#' @export
#'
PROSPECT <- function(SpecPROSPECT = NULL, Input_PROSPECT = NULL,
                     N = 1.5, CHL = 40.0, CAR = 8.0, ANT = 0.0, BROWN = 0.0,
                     EWT = 0.01, LMA = NULL, PROT = 0, CBC = 0, alpha = 40.0,
                     check = TRUE) {

  .Deprecated("prospect")
  prospect(spec_prospect = SpecPROSPECT, input_prospect = Input_PROSPECT,
           n_struct = N, chl = CHL, car = CAR, ant = ANT, brown = BROWN,
           ewt = EWT, lma = LMA, prot = PROT, cbc = CBC, alpha = alpha,
           check = check)
}
