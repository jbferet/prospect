#' SpecPROSPECT_FullRange: optical constants defined for PROSPECT
#'
#' Corresponds to spectral bands, refractive index and specific absorption
#' coefficient for each chemica constituent, defined over the spectral domain
#' ranging from 400 nm to 2500 nm
#'
#' The specific absorption coefficients were calibrated using experimental data.
#' The details of the calibration should be found in the following publications:
#'  http://dx.doi.org/10.1016/j.rse.2017.03.004
#'  https://doi.org/10.1016/j.rse.2020.112173
#'
SpecPROSPECT_FullRange <- read.table(file = 'data-raw/dataSpec_PRO.txt',
                           header = TRUE,
                           sep = '\t')

#' calctav_90 & calctav_40: transmissivity of a dielectric plane surface,
#' averaged over all directions of incidence and over all polarizations for
#' solid angle of 90 and 40 degrees
#'
calctav_90 <- prospect::calctav(90, prospect::SpecPROSPECT_FullRange$nrefrac)
calctav_40 <- prospect::calctav(40, prospect::SpecPROSPECT_FullRange$nrefrac)
SpecPROSPECT_FullRange$calctav_90 <- calctav_90
SpecPROSPECT_FullRange$calctav_40 <- calctav_40

#' OptDomain_RT: optimal spectral domains defined to assess the different leaf
#' chemical constituents.
#' These optimal spectral domains were defined in
#' Féret et al. (2019) https://doi.org/10.1016/j.rse.2018.11.002
#' Féret et al. (2021) https://doi.org/10.1016/j.rse.2020.112173
#' Spafford et al. (2021) https://doi.org/10.1016/j.rse.2020.112176
#'
OptDomain_RT <- list('CHL' = seq(700,720),
                     'CAR' = seq(520,560),
                     'ANT' = seq(400,800),
                     'EWT' = seq(1700,2400),
                     'LMA' = seq(1700,2400),
                     'PROT' = c(seq(2100,2139), seq(2160,2179)),
                     'CBC' = c(seq(1480,1499), seq(1560,1579), seq(1760,1799),
                               seq(2040,2059), seq(2120,2139), seq(2160,2239),
                               seq(2260,2279), seq(2340,2359), seq(2380,2399)))



usethis::use_data(SpecPROSPECT_FullRange, OptDomain_RT,
                  internal = FALSE,
                  overwrite = TRUE)
