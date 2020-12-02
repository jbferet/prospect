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
