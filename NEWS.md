# prospect v1.6.2
## fix
- fixed wrong merge for correction of JOSS paper

## addition
- added merit function Merit_PROSPECT_dRMSE to perform inversion based on 1st derivative

## change
- moved data.table from Imports to Suggestions

# prospect v1.6.1
## addition
- added input dataframe InputPROSPECT to function PROSPECT

# prospect v1.6.0
Released version of the software with the latest changes from the review

# prospect v1.5.1
version following reviews for publication in JOSS

## addition
- added a function download_LeafDB to get leaf databases from gitlab repository

## Changes
- modified defaulty merit function as one function instead of two
- updated documentation
- updated function FitSpectralData

# prospect v1.5.0
version following reviews for publication in JOSS

## Changes
- change license to MIT (Issue#3:  Consider a different license)
- include plain text version of specific absorption coefficients with variable prospect::SpecPROSPECT_FullRange (Issue#4:  Consider including plain-text versions of PROSPECT coefficients)
- remove function plotinv and corresponding required packages (Issue#5: Move randomcoloR to Suggests)
- pre-calculate calctav for angles 40 degrees (default) and 90 degrees (always needed) and add them to SpecPROSPECT_FullRange (Issue#6: Pre-calculate calctav)
- added tests for PROSPECT, Invert_PROSPECT and Invert_PROSPECT_OPT (Issue#7: Add automated unit tests)
- add CONTRIBUTION.md and CODE_OF_CONDUCT.md (Issue#8: Add contribution guide)
- update documentation and harmonize with corresponding code for JOSS manuscript (Issue#10: harmonize web tutorial with JOSS manuscript examples)
- overall update of the code to work with dataframes more systematically
- only two versions now defined in the package (D and PRO), and the use of brown pigments is accessible during inversion by setting Est_Brown_Pigments == TRUE

## fixes
- fixed error occuring when PROSPECT-PRO inversion performed to assess EWT

# prospect v1.4.0

## fixes
deal with default values when seting proteins, CBC and LMA to non constant value

# prospect v1.3.0

new release including JOSS draft

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
