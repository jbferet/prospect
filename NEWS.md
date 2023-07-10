# prospect v1.2.7

## fixes
- removed reference to multiprocess, multisession used exclusively from future

# prospect v1.2.6

## fixes
fixed issue raised by morelju on 2023/05/31
--> bug when running prospect with default parameterization because NULL value included in dataframe
LRT_default <- PROSPECT(SpecPROSPECT)
Now runs

# prospect v1.2.5

## addition
- added function `plotinv` to produce scatteerplots between measured and estimated values

# prospect v1.2.4

## Fixes
- Fixed bug when setting Parms2Estimate = "ALL" with Invert_PROSPECT_OPT

# prospect v1.2.3

## Fixes
- Fixed bug initial inversion outputs NA

# prospect v1.2.2

## Fixes
- added default xlub which are not specified by user, and result in an error


# prospect v1.2.1

## Changes
- optional progress bar when inverting PROSPECT
- possibility to modify lower and upper bounds for parameters to estimate from inversion

# prospect v1.2.0

## Changes
- now accepts dataframes, matrices and vectors as input leaf optical properties
- optimal inversion codes rewritten and simplified
- added progress bar when inverting PROSPECT on several samples
- updated vignettes

# prospect v1.1.0
- added verbose as input in invert_PROSPECT_OPT
- modified info displayed when performing optimal inversion (Invert_PROSPECT_OPT) 
- changed lambda to Sublambda in Invert_PROSPECT_OPT
 
## Fixes
- Added library NlcOptim in function `tryInversion`
- fixed bug confusing CBC with proteins when performing optimal spectral domain inversion
- fixed input parameters for FitSpectralData: added `UL_Bounds=TRUE` to specify that UserDomain corresponds to upper and lower boundaries, not explicit definition of spectral bands

# prospect v1.0.1

## Changes
- Added NEWS.md

## Fixes
- Fixed bug occuring when only reflectance or only transmittance defined
- Fixed case when input data for inversion is not 1 nm sampling

# prospect v1.0.0
First public release in Gitlab

The package prospect includes the PROSPECT leaf model in its two recent versions: 
PROSPECT-D, accounting for chlorophyll, carotenoids, anthocyanins, equivalent water thickness and leaf mass per area
PROSPECT-PRO, accounting for chlorophyll, carotenoids, anthocyanins, equivalent water thickness, proteins and carbon-based constituents

Brown pigments car be included additionally to the other chemical constituents for both versions of the model.

The package also includes methods for the inversion of PROSPECT using leaf optical properties (directional-hemispherical reflectance, directional-hemispherical transmittance, or both) in order to estimate leaf chemistry and structure parameter.
These methods are based on iteerative optimization of a merit function minimizing RMSE between measured and estimated leaf optical properties, but custom functions can be implemented by users.

Optimal subdomains and prior N estimation as introduced in recent publications have also been included in the package.

References: 
[PROSPECT-D](https://doi.org/10.1016/j.rse.2017.03.004)
[PROSPECT-PRO](https://doi.org/10.1016/j.rse.2020.112173)
[Prior N estimation for inversion using R or T only](https://doi.org/10.1016/j.rse.2020.112176)
