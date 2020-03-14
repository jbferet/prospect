# v0.0.0.9001

## Changes

- `Invert_PROSPECT`: `x0` has now the role of initial value or fixed value of prospect parameter depending on parameters specified in `Parms2Estimate`.
- `tryInversion`: `xinit` becomes `x0` to fit with names of fmincon, remove `InPROSPECT`  (useless as all can be deduced from x0). Also, within the function, the tolerance argument of `fmincon` is set to 1e-15 (instead of 1e-6) as it shows results more coherent with matlab version of `fmincon`.
- `Merit_RMSE_PROSPECT`: `xinit` becomes `x` as it is used at each step of optimisation and not only at initialization.
- `PROSPECT`: remove arg `Input_PROSPECT` as it is useless. To call the function with a list/data.frame of arguments use do.call().

