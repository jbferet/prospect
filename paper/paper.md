---
title: '`prospect`: an R package to link leaf optical properties with their chemical and structural properties with the leaf model PROSPECT'
tags:
  - R
  - leaf spectroscopy
  - vegetation biophysical properties
  - remote sensing
  - physical modeling
  - iterative optimization
authors:
  - name: Jean-Baptiste Féret
    orcid: 0000-0002-0151-1334
    corresponding: true
    affiliation: 1
  - name: Florian de Boissieu
    orcid: 0000-0002-2185-9952
    affiliation: 1
affiliations:
 - name: TETIS, INRAE, AgroParisTech, CIRAD, CNRS, Université Montpellier, Montpellier, France
   index: 1
date: 17 July 2023
bibliography: paper.bib
---

# Summary

PROSPECT simulates leaf optical properties based on a limited set of light 
absorbing chemical constituents, and a unique leaf structure parameter to 
account for light scattering.
We present `prospect`, an R package which includes the most recent versions 
of the model, complemented with multiple inversion routines to assess leaf 
chemistry from leaf optical properties.

# Statement of need

The capacity to measure, map and monitor vegetation traits corresponding to 
biophysical and chemical properties is crucial to better understand ecosystem 
and agrosystem functions, as well as carbon, water and energy budgets. 
At leaf scale, these vegetation traits are linked to their 
optical properties through absorption and scattering mechanisms. 
Techniques based on spectroscopy have developed in the past decades to provide 
rapid, accurate and non-destructive assessment of leaf chemical composition. 
Physical models aim at simulating optical properties of leaves from their 
chemical composition, and they are also used to assess chemical composition 
from leaf optical properties based on inversion techniques. 
The model PROSPECT (leaf optical PROperties SPECtra) is currently the most popular 
physical model used for the simulation of leaf optical properties. 
PROSPECT simulates leaf directional-hemispherical reflectance and transmittance 
from a combination of chemical constituents and their corresponding specific 
absorption coefficients. 
It uses a simplified representation of leaf structure through the generalised 
plate model [@allen1970] to account for scattering. 
PROSPECT is also coupled with canopy reflectance models to analyze Earth 
observation data, such as satellite imagery. Hence, it is a key component for 
remote sensing applications dedicated to vegetation monitoring. 

Multiple versions aiming at expanding the pool of chemical constituents accounted 
for by PROSPECT have been released since its first version [@jacquemoud1990].
Successive versions introduced carotenoids [@feret2008] and anthocyanins [@feret2017], 
to simulate leaf optical properties from juvenile to mature and senescent 
development stages. 
The latest version, PROSPECT-PRO, separates dry matter constituents into 
proteins and carbon based constituents [@feret2021]. 
In parallel with the updated versions of the model, model inversion strategies 
are also developed to improve the assessment of leaf chemical constituents 
[@feret2019; @spafford2021].

An overview of different implementations of PROSPECT since [@feret2008] is 
currently maintained at [this webpage](http://teledetection.ipgp.jussieu.fr/prosail/).
This includes distributions in matlab, R and fortran programming languages. 
It also provides links to distribution corresponding to the coupling of PROSPECT 
with vegetation models such as COSINE dedicated to close-range imaging 
spectroscopy [@jay_physically-based_2016] and PROSAIL for canopy reflectance modeling 
[verhoef_coupled_2007, @verhoef_unified_2007, @jacquemoud_prospect+_2009]. 
Note that PROSPECT is also available in packages written in [python](https://github.com/jgomezdans/prosail), 
[Julia](https://github.com/RemoteSensingTools/CanopyOptics.jl) and [R](https://github.com/ashiklom/rrtm).

PROSPECT versions 4 and 5 developed by [@feret2008] are deprecated. 
Therefore PROSPECT-D and PROSPECT-PRO are recommended versions. 
The R package `prospect` includes recent versions of PROSPECT (D and PRO), 
and aims at providing up-to-date developments in terms of future model version 
and inversion strategies. 
This includes parameterizable inversion routines, to ease physical model inversion 
for beginners, and to help advanced users in designing and testing their own 
inversion strategy. 


# Overview

## PROSPECT simulation in forward mode

The application of PROSPECT in forward mode results in the simulation of leaf 
directional-hemispherical reflectance and transmittance from a set of 
chemical constituents and a unique leaf structure parameter usually identified as `N`
and corresponding to the number of homogeneous layers introduced in the 
generalized plate model. 

Two versions of the model PROSPECT are implemented in `prospect`: PROSPECT-D 
[@feret2017] and PROSPECT-PRO [@feret2021]. 
Table \ref{table:1} lists the leaf chemical constituents included in the 
different versions which can be specified when calling the function `PROSPECT`.

| Version 	|  `D` 	| `PRO` |
|---------	|:----:	|:----:	|
| CHL     	| **X** | **X** |
| CAR   	| **X** | **X** |
| ANT   	| **X** | **X** |
| BROWN 	| **X** | **X** |
| EWT   	| **X** | **X** |
| LMA   	| **X** |     	|
| PROT  	|      	| **X** |
| CBC   	|    	| **X** |

: Versions of the model PROSPECT available in the `prospect` package and 
corresponding chemical constituents (CHL: chlorophylls; CAR: carotenoids; 
ANT: anthocyanins; BROWN: brown pigments; EWT: equivalent water thickness; 
LMA: leaf mass per area; PROT: proteins; CBC: carbon based constituents).\label{table:1}

As PROSPECT is a relatively simple and computationally efficient model, 
iterative optimization is the most widespread method to invert PROSPECT 
and assess leaf chemistry and structure from their optical properties. 
It usually takes less than 1 second to perform PROSPECT inversion, which is acceptable
when processing experimental datasets of leaf optical properties including 100-1000 samples.
Iterative optimization aims at minimizing a cost function comparing 
measured and simulated leaf optical properties. 
This procedure is based on the function `fmincon` included in the package `pracma`.

Various inversion strategies using iterative optimization are described in the literature. 
These inversion strategies differ either by the cost function, or by the 
selection of specific spectral domains used to retrieve one or several leaf 
biophysical properties, or by the introduction of prior information. 
The default cost function, `CostVal_RMSE`, corresponds to the root mean square 
of the mean quadratic difference between measured and simulated leaf optical 
properties (reflectance and/or transmittance).
Users can define their own cost function. 
Table \ref{table:2} provides information on the optimal spectral range used to 
assess leaf chemical constituents from their optical properties, as identified 
by [@feret2019; @spafford2021].

| Constituent|    Optimal spectral domain     	|   Versions   |
|----------- |---------------------------------	|:-----------: |
| CHL        |           700 -- 720           	| `D`, `PRO`   |
| CAR        |           520 -- 560            	| `D`, `PRO`   |
| ANT        |           400 -- 800           	| `D`, `PRO`   |
| BROWN      |               NA               	| NA           |
| EWT        |          1700 -- 2400           	| `D`, `PRO`   |
| LMA        |          1700 -- 2400           	| `D`          |
| PROT       | 2100 -- 2139; 2160 -- 2179     	| `PRO`        |
| CBC        | 1480 -- 1499;	1560 -- 1579;	1760 -- 1799;	2040 -- 2059;	2120 -- 2139;	2160 -- 2239;	2260 -- 2279;	2340 -- 2359;	2380 -- 2399 | `PRO`        |

: Optimal spectral domains selected to assess vegetation chemical constituents
 from leaf optical properties (CHL: chlorophylls; CAR: carotenoids; 
ANT: anthocyanins; BROWN: brown pigments; EWT: equivalent water thickness; 
LMA: leaf mass per area; PROT: proteins; CBC: carbon based constituents).\label{table:2}

# Example 1: running PROSPECT in forward mode

## Individual simulations

PROSPECT simulates leaf directional-hemispherical reflectance and transmittance 
from the leaf structure parameter and a combination of chemical constituents 
when running in forward mode.
PROSPECT relies on specific absorption coefficients corresponding to each 
chemical constituent accounted for by the model, and on 
the leaf refractive index. 
These optical constants are identical for all leaves, and are defined over the 
spectral domain ranging from 400 nm to 2500 nm.
They are accessible through the variable `SpecPROSPECT_FullRange`, a data frame 
which is automatically loaded with the package. 
Multiple functions in the package `prospect` expect a variable `SpecPROSPECT` as input. 
`SpecPROSPECT_FullRange` is used as default value when this input variable is not defined. 
Users can adjust `SpecPROSPECT_FullRange` in order to adjust the spectral characteristics 
for simulation, including spectral domain, sampling or sensor-specicif spectral response. 

These two examples illustrate how to run PROSPECT. 
The function `PROSPECT` identifies the version to be used: PROSPECT-D is used 
if LMA is defined, while PROSPECT-PRO is used if proteins (PROT) and carbon 
based constituents (CBC) are defined. 
If LMA, PROT and CBC are defined simultaneously, PROSPECT-PRO is used and LMA 
is set to 0.
Figure \ref{fig:LOP} compares simulated leaf optical properties. 
Here, the differences between PROSPECT-D and PROSPECT-PRO are mainly driven by the 
difference set for the `N` structure parameter.

```r
# Load prospect package
library(prospect)
# Run PROSPECT-D
LRT_D <- PROSPECT(CHL = 45, CAR = 10, ANT = 0.2, 
                  EWT = 0.012, LMA = 0.010, N = 1.3)
# Run PROSPECT-PRO
LRT_PRO <- PROSPECT(CHL = 45, CAR = 10, ANT = 0.2, 
                    EWT = 0.012, PROT = 0.001, CBC = 0.009, N = 1.7)
```

![Leaf optical properties simulated with PROSPECT-D and PROSPECT-PRO. Different values of N were defined to highlight differences in simulated leaf optics \label{fig:LOP}](compare_RT_PROSPECT_PRO_D.png){ width=85% }

The spectral domain covered with PROSPECT simulations can be adjusted with the function 
`FitSpectralData`, which adjust information from `SpecPROSPECT_FullRange` to a user-defined 
spectral domain. 

```r
# define spectral bands for Visible / Near InfraRed (VNIR) simulation
wvlRange_VNIR <- seq(400,1000)
# adjust spectral properties used in PROSPECT to VNIR domain
VNIR <- FitSpectralData(lambda = wvlRange_VNIR)
# Run PROSPECT-D in VNIR domain
LRT_VNIR <- PROSPECT(SpecPROSPECT = VNIR$SpecPROSPECT,
                     N = 1.4, CHL = 30, CAR = 6, EWT = 0.02, LMA = 0.01)
```

## Simulation of a Look-Up-Table

Look-Up-Tables (LUTs) are widely used in order to infer leaf characteristics from PROSPECT. 
The function `PROSPECT_LUT` computes a LUT based on a list of input parameters.
Undefined parameters are set to their default value. Vectors of values are expected 
to be the same length.
The output of `PROSPECT_LUT` is a list containing a data frame including the 
input parameters, a reflectance data frame and a transmittance data frame.

```r
# define input parametrs for PROSPECT
Input_PROSPECT <- data.frame('CHL' = 100*runif(1000), 
                             'CAR' = 25*runif(1000), 
                             'ANT' = 2*runif(1000), 
                             'EWT' = 0.04*runif(1000), 
                             'LMA' = 0.02*runif(1000), 
                             'N' = 1+2*runif(1000))
# produce a LUT defined over the VSWIR domain covered by PROSPECT
LUT <- PROSPECT_LUT(Input_PROSPECT = Input_PROSPECT)
# produce a LUT defined over the VNIR domain
LUT_VNIR <- PROSPECT_LUT(SpecPROSPECT = VNIR$SpecPROSPECT,
                         Input_PROSPECT = Input_PROSPECT)
```

# Example 2: PROSPECT inversion using iterative optimization

The package `prospect` offers possibilities to adjust these parameters for inversion. 
Here, we will illustrate different types of inversion with an experimental database 
named __ANGERS__.
This database was used to calibrate the model PROSPECT, and is among the most 
popular public data sets in the domain of leaf spectroscopy.

## Downloading the ANGERS experimental dataset

A version of the ANGERS data set can be downloaded from a gitlab repository.

```r
# use data.table library
library(data.table)
# download ANGERS leaf optics database from gitlab repository
gitlab_Rep <- 'https://gitlab.com/jbferet/myshareddata/raw/master/LOP'
# download leaf biochemical constituents and leaf optical properties
DataBioch <- data.table::fread(file.path(gitlab_Rep,'ANGERS/DataBioch.txt'))
Refl<- data.table::fread(file.path(gitlab_Rep,'ANGERS/ReflectanceData.txt'))
Tran <- data.table::fread(file.path(gitlab_Rep,'ANGERS/TransmittanceData.txt'))
# Get the wavelengths corresponding to reflectance and transmittance measurements
lambda <- Refl$wavelength
Refl$wavelength <- Tran$wavelength <- NULL
```

## PROSPECT inversion using the full spectral information

The spectral domains covered by `SpecPROSPECT` and the leaf optical properties 
are expected to match when performing inversion on the full spectral domain 
covered by the data.
The function `FitSpectralData` automatically harmonizes the spectral domain 
for `SpecPROSPECT` and for the leaf optical properties. 
`UserDomain` corresponding to the wavelengths of interest can also be provided in addition. 
Here the following R code aims at adjusting `SpecPROSPECT` based on the spectral 
domain covered by the leaf optical properties.


```r
# Adjust spectral domain for SpecPROSPECT to fit leaf optical properties 
SubData <- FitSpectralData(SpecPROSPECT = SpecPROSPECT_FullRange,
                           Refl = Refl, Tran = Tran,
                           lambda = lambda, UserDomain = lambda)
```

The main inversion procedure is called with the function `Invert_PROSPECT`, 
which minimizes a cost function. 
They can also select the biophysical 
properties to assess: 

The full set or a selection of chemical constituents can be assessed from 
PROSPECT inversion. 
This list of parameters is defined with the variable `Parms2Estimate`. 
The default parameterization of PROSPECT inversion assesses all parmeters 
listed in Table \ref{table:1} except BROWN, which can be assessed by setting 
`Est_Brown_Pigments = TRUE` as input for `Invert_PROSPECT`.
The value set for parameters which are not assessed is then defined with 
the `InitValues` input variable. 

```r
# Assess all parameters using PROSPECT inversion applied to full spectral data
res_all_WL <- Invert_PROSPECT(SpecPROSPECT = SubData$SpecPROSPECT, 
                              Refl = SubData$Refl,  
                              Tran = SubData$Tran)
```

## PROSPECT inversion using optimal spectral domains for each constituent

The function `Invert_PROSPECT_OPT` performs PROSPECT inversion using optimal 
spectral domains defined in [@feret2019; @spafford2021]. 
The optimal spectral domains are specific to each constituent, and the function 
automatically adjusts the spectral domain of the leaf optical properties provided 
by user. 

```r
# Assess a set of parameters using PROSPECT inversion with optimal spectral domains
Parms2Estimate <- c('CHL', 'CAR', 'EWT', 'LMA')
res_opt_WL <- Invert_PROSPECT_OPT(lambda = lambda, 
                                  Refl = Refl, Tran = Tran
                                  Parms2Estimate = Parms2Estimate)
```

## Performances of the two types of inversion: Comparison with ANGERS data

Figure \ref{fig:scatter} displays the outputs of the inversion obtained either from 
`Invert_PROSPECT` or from `Invert_PROSPECT_OPT` for each 
leaf chemical constituent.

![Estimation of chlorophyll content, carotenoid content, EWT and LMA from PROSPECT inversion applied on the ANGERS data set. `full WL` corresponds to the inversion performed with the full spectral information (function `Invert_PROSPECT`); `opt WL` corresponds to the inversion performed with the optimal spectral information (function `Invert_PROSPECT_OPT`). \label{fig:scatter}](PROSPECT_Inversions.png){ width=90% }

# Conclusion

We have described `prospect`, an R package dedicated to the PROSPECT leaf model. 
`prospect` can run different versions of the model in direct mode to simulate 
directional-hemispherical reflectance and transmittance.
`prospect` also includes inversion routines to assess leaf structure and chemical 
constituent content either from directional-hemispherical reflectance and 
transmittance, or from reflectance or transmittance only. 
`prospect` provides latest advances in terms of model version and inversion 
procedures to the leaf spectroscopy community. 
`prospect` is coupled with the canopy model SAIL through the R package [`prosail`](https://jbferet.gitlab.io/prosail/index.html).
`prosail` is dedicated to applications focusing on Earth observation imagery analysis 
and allows simulation of canopy reflectance for multispectral and hyperspectral sensors.
Hybrid inversions based on physical modeling and machine learning are also implemented 
in `prosail` to assess vegetation traits from imagery data.


# Availability

prospect is an open-source software package made available under the MIT license. 
It can be installed through GitHub repository: remotes::install_github("jbferet/prospect"). 
Tutorials are available at [https://jbferet.gitlab.io/prospect/](https://jbferet.gitlab.io/prospect/).

# Acknowledgements

The authors acknowledge financial support from Agence Nationale de la Recherche 
(BioCop project — ANR-17-CE32-0001)
We are grateful to Stéphane Jacquemoud and Frédéric Baret for the development of 
the initial version of the PROSPECT model. 
We also warmly thank Luc Bidel, Christophe François and Gabriel Pavan who 
collected the ANGERS data set.

# References
