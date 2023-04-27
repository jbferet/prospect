---
title: 'PROSPECT: an R package to link leaf optical properties with their chemical and structural properties with the leaf model PROSPECT'
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
date: 21 March 2023
bibliography: paper.bib
---

# Summary

PROSPECT is the most widely used leaf model, allowing simulation of leaf optical properties 
based on a limited set of chemical constituents absorbing light, and a unique leaf 
structure parameter to account for light scattering.

We present `prospect`, an R package which includes PROSPECT-D [@feret2017] and PROSPECT-PRO [@feret2021], 
its two most recent versions.
The package also provides inversion routines based on iterative optimization to estimate leaf chemistry from LOP.
Various inversion strategies are implemented in the package, including strategies recently developed 
[@feret2019, @spafford2021] aiming at estimating leaf chemical constituents based on optimal spectral domains.

# Statement of need

Leaf optical properties are linked to leaf biophysical properties through absorption 
and scattering mechanisms. Physical models describe the interactions between these properties.
PROSPECT takes advantage of specific absorption coefficients from chemical constituents 
such as leaf pigments, water and proteins, and uses a simplified representation of leaf structure 
through the generalised plate model [@allen1970] to simulate leaf directional-hemispherical 
reflectance and transmittance.
Leaf optical properties contribute significantly to vegetation reflectance acquired from Earth observation systems.
The capacity to map and monitor leaf chemical properties is key to better understand ecosystems 
and agrosystems functions, as well as carbon, water and energy budgets, and models and methods 
giving access to these core leaf constituents based on observations at leaf and canopy scale are 
important components.

Techniques based on leaf spectroscopy have developed in the past decades to provide a 
rapid, accurate and non-destructive estimation of leaf constituents. Physical modeling has an
important role in these techniques, as it can be used to directly estimate leaf chemical content
based on model inversion. It can also generate simulated leaf optical properties, in order to adjust 
statistical models or to identify spectral indices leading to accurate estimation of leaf chemistry
based on parcimonious spectral information [@feret2011]. 

The model PROSPECT is widely used in the remote sensing community, 
alone or combined with canopy reflectance models such as SAIL [@verhoef2007] and DART [@gastellu2015]. 
Multiple versions have been released since the first version of PROSPECT [@jacquemoud1990], aiming at expanding 
its capacity to account for more chemical constituents. Recent versions introduced carotenoids [@feret2008],  
and anthocyanins [@feret2017], allowing simulating and analyzing leaf optical properties corresponding to different 
leaf development stages, from juvenile to mature and senescent. The latest version, PROSPECT-PRO, 
separated dry matter constituents into proteins and carbon based constituents. In parallel with the 
development of new versions of the model, advances in model inversion were also performed 
to improve the estimation of leaf chemical constituents [@feret2019; @spafford2021].

The application of an appropriate model inversion strategy is critical to ensure optimal 
estimation of leaf constituents, as it was evidenced for the estimation of leaf mass per area (LMA) [@feret2019], 
which was reported as poorly retrieved from PROSPECT inversion in earlier studies. 

To ensure access to latest advances in terms of model version and inversion strategies, 
we present `prospect`, an R package including the most recent versions of PROSPECT, 
and parameterizable inversion routines, allowing users to design their own inversion strategy.


# Overview

Two versions of the model PROSPECT are implemented in `prospect`: PROSPECT-D [@feret2017] and 
PROSPECT-PRO [@feret2021]. For each version, it is possible to also include the influence brown pigments (BROWN), 
which appear during senescence, by ending the name of the version with `B`. 

Table \ref{table:1} provides an overview of the leaf chemical constituents included in the different versions which can be specified when calling the function `PROSPECT`.

| Version 	| `D` 	| `DB` 	| `PRO` | `PROB`|
|---------	|------	|:----:	|:----:	|:----:	|
| CHL     	| **X** | **X** | **X** | **X** |
| CAR   	| **X** | **X** | **X** | **X** |
| ANT   	| **X** | **X** | **X** | **X** |
| EWT   	| **X** | **X** | **X** | **X** |
| LMA   	| **X** | **X** |     	|     	|
| PROT  	|     	|      	| **X** | **X** |
| CBC   	|      	|    	| **X** | **X** |
| BROWN 	|   	| **X** |    	| **X** |

: Versions of the model PROSPECT available in the `prospect` package and 
corresponding chemical constituents (CHL: chlorophylls; CAR: carotenoids; 
ANT: anthocyanins; EWT: equivalent water thickness; LMA: leaf mass per area; 
PROT: proteins; CBC: carbon based constituents; 
BROWN: brown pigments).\label{table:1}

The inversion procedure is based on the function `fmincon` implemented in the package `pracma`. 
`fmincon` is an optimization algorithm aiming at finding the minimum of multivariable 
functions with nonlinear constraints. 


Table \ref{table:2} provides information on the optimal spectral range used to estimate leaf 
chemical constituents from their optical properties, as identified by [@feret2019; @spafford2021].

| Constituent 	|    Optimal spectral domain 	|            Versions     	|
|------------	|-----------------------------	|:-----------------------:	|
| CHL       	|           700 -- 720       	| `D`, `DB`, `PRO`, `PROB`	|
| CAR       	|           520 -- 560        	| `D`, `DB`, `PRO`, `PROB`	|
| ANT       	|           400 -- 800       	| `D`, `DB`, `PRO`, `PROB` 	|
| BROWN     	|               NA           	| NA                       	|
| EWT       	|          1700 -- 2400       	| `D`, `DB`, `PRO`, `PROB`	|
| LMA       	|          1700 -- 2400       	| `D`, `DB`               	|
| PROT      	| 2100 -- 2139 ; 2160 -- 2179 	| `PRO`, `PROB`            	|
| CBC       	| 1480 -- 1499;	1560 -- 1579;	 1760 -- 1799; 2040 -- 2059;	 2120 -- 2139; 2160 -- 2239;	 2260 -- 2279; 2340 -- 2359;	 2380 -- 2399 	| `PRO`, `PROB`                      	|

: Optimal spectral domains selected to estimate vegetation chemical constituents from leaf 
optical properties (CHL: chlorophylls; CAR: carotenoids; 
ANT: anthocyanins; EWT: equivalent water thickness; LMA: leaf mass per area; 
PROT: proteins; CBC: carbon based constituents; 
BROWN: brown pigments).\label{table:2}

# Example 1: running PROSPECT in forward mode

## Individual simulation with PROSPECT-D and PROSPECT-PRO

This section illustrates how to run PROSPECT in forward mode in order to simulate 
leaf directional-hemispherical reflectance and transmittance from leaf structure parameter and a combination 
of chemical constituents.

The variable `SpecPROSPECT` is automatically available as a dataframe when loading `prospect`. 
`SpecPROSPECT` includes the leaf refractive index and all specific absorption coefficents defined on the spectral 
range from 400 nm to 2500 nm.

The first example uses PROSPECT-D. The function autmatically identifies the version to be used (PROSPECT-D), as 
the input parameters correspond to this version: anthocyanins are set to a value, brown pigments, proteins, and 
carbon based constituents (CBC) are not defined.

```r
# Load prospect package
library(prospect)
# Run PROSPECT-D
LRT_D <- PROSPECT(SpecPROSPECT, CHL = 45, CAR = 10, ANT = 0.2, 
                  EWT = 0.012, LMA = 0.010, N = 1.3)
```

In this second example, PROSPECT-PRO is used as proteins and CBC are defined.  

```r
# Run PROSPECT-PRO
LRT_PRO <- PROSPECT(SpecPROSPECT, CHL = 45, CAR = 10, ANT = 0.2, 
                    EWT = 0.012, PROT = 0.001,  CBC = 0.009, N = 1.7)
```

The leaf optical properties can then be compared. Here, the differences between PROSPECT-D and PROSPECT-PRO 
are mainly driven by the difference set for the `N` parameter.

```r
R_D <- data.frame('wvl' = LRT_D$wvl, 'RT' = 100*LRT_D$Reflectance,
                  'LOP' = 'R (PROSPECT-D)')
T_D <- data.frame('wvl' = LRT_D$wvl, 'RT' = 100*(1-LRT_D$Transmittance),
                  'LOP' = 'T (PROSPECT-D)')
R_PRO <- data.frame('wvl' = LRT_PRO$wvl, 'RT' = 100*LRT_PRO$Reflectance,
                    'LOP' = 'R (PROSPECT-PRO)')
T_PRO <- data.frame('wvl' = LRT_PRO$wvl, 'RT' = 100*(1-LRT_PRO$Transmittance),
                    'LOP' = 'T (PROSPECT-PRO)')
LRT_df <- rbind(R_D, T_D, R_PRO, T_PRO)

# plot reflectance and transmittance for both models
library(ggplot2)
RTplot <- ggplot2::ggplot(LRT_df, aes(x=wvl, y=RT, group=LOP)) +
  geom_line(aes(linetype=LOP, color=LOP),linewidth=1.00)+
  scale_color_manual(values=c('#FF9999','red1','#9999FF','blue4'))+
  scale_size_manual(values=c(5, 5))+
  labs(x='Wavelength (nm)',y='Reflectance     (%)     100-Transmittance') +
  theme(legend.position = "bottom",
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) +
  scale_linetype_manual(values=c("solid","solid","solid","solid"))

filename = 'compare_RT_PROSPECT_PRO_D.png'
ggsave(filename, plot = RTplot, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm", dpi = 300)
```
![Leaf optical properties simulated with PROSPECT-D and PROSPECT-PRO. Different values of N were defined to highlight differences in simulated leaf optics \label{fig:AIM}](compare_RT_PROSPECT_PRO_D.png){ width=85% }


## Simulation of a look up table with PROSPECT-D

Look-Up-Tables (LUTs) are widely used in order to infer leaf charactristics from PROSPECT, based on minimization techniques. 
The function `PROSPECT_LUT` allows computation of a LUT directly based on a list of input parameters.
Undefined parameters are set to their default value. Vectors of values are expected to be the same length.
The output of `PROSPECT_LUT` is a list containing a dataframe including the input parameters, as well as two matrices 
corresponding to leaf directional-hemispherical reflectance and transmittance.

```r
CHL <- 100*runif(1000)
CAR <- 25*runif(1000)
ANT <- 2*runif(1000)
EWT <- 0.04*runif(1000)
LMA <- 0.02*runif(1000)
N   <- 1+2*runif(1000)
Input_PROSPECT <- data.frame('CHL' = CHL, 'CAR' = CAR, 'ANT' = ANT, 
                             'EWT' = EWT, 'LMA' = LMA, 'N' = N)
LUT <- PROSPECT_LUT(SpecPROSPECT,Input_PROSPECT)
```

# Example 2: PROSPECT inversion using iterative optimization

Several approaches can be used to perform PROSPECT inversion. As PROSPECT is a relatively simple and 
computationally efficient model, inversion based on iterative optimization is one of the most popular 
method to invert PROSPECT and estimate leaf chemistry and structure from their optical properties. 
Here, the iterative optimization is based on the minimization of a multivariable function with nonlinear constraints. 
This procedure is based on the function `fmincon` included in the package `pracma`, which uses 
Sequential Quadratic Programming.

Various inversion strategies have been proposed in the literature since the first version of PROSPECT. 
These inversion strategies may differ in various ways: either by the cost function, or by the 
selection of specific spectral domains aiming at optimizing the retrieval of one or several 
leaf biophysical properties, or by the introduction of prior information. 

The package `prospect` offers possibilities to adjust these parameters for inversion. 
Here, we will illustrate different types of inversion with an experimental database named __ANGERS__.
This database was used to calibrate the model PROSPECT, and is among the most popular public 
datasets in the domain of leaf spectroscopy.

## Downloading the ANGERS dataset and adjusting leaf optics

A version of the ANGERS dataset is hosted on gitlab, and can be directly downloaded from R.

```r
# use data.table library
library(data.table)

# download ANGERS leaf optics database from gitlab repository
gitlab_Rep <- 'https://gitlab.com/jbferet/myshareddata/raw/master/LOP'
dbName <- 'ANGERS'

# download leaf biochemical constituents and leaf optical properties
fileName <- list('DataBioch.txt','ReflectanceData.txt','TransmittanceData.txt')
DataBioch <- data.table::fread(file.path(gitlab_Rep,dbName,fileName[[1]]))
Refl<- data.table::fread(file.path(gitlab_Rep,dbName,fileName[[2]]))
Tran <- data.table::fread(file.path(gitlab_Rep,dbName,fileName[[3]]))

# Get the wavelengths corresponding to reflectance and transmittance measurements  
lambda <- unlist(Refl$V1, use.names=FALSE)
Refl$V1 <- Tran$V1 <- NULL
```

## PROSPECT inversion using the full spectral range available

The spectral domain covered by `SpecPROSPECT` and the leaf optical properties are expected to match 
when performing inversion on the full spectral domain covered by the data.
This can be done automatically using `FitSpectralData`

```r
# Adjust spectral domain for SpecPROSPECT to fit leaf optical properties 
SubData <- FitSpectralData(SpecPROSPECT = SpecPROSPECT, # spectral properties
                           lambda = lambda, # spectral bands
                           Refl = Refl, # Reflectance data (matrix or dataframe)
                           Tran = Tran, # Transmittance data (matrix or dataframe)
                           UserDomain = lambda) # spectral domain of interest
						   
SubSpecPROSPECT <- SubData$SpecPROSPECT # spectral properties adjusted to current data
Sublambda <- SubData$lambda # updated spectral domain defined by UserDomain
SubRefl <- SubData$Refl # updated reflectance domain defined by UserDomain
SubTran <- SubData$Tran # updated transmittance defined by UserDomain
```

The main inversion procedure is called with the function `Invert_PROSPECT`.
It is based on the minimization of a cost function. The default cost function, 
`CostVal_RMSE`, corresponds to the root mean square of the mean quadratic difference between 
measured and simulated leaf optical properties (either reflectance and transmittance, or one of them only).
Users can define their own cost function. They can also select the biophyscal properties to estimate: 
The full set (`Parms2Estimate` set to `ALL`) or a reduced set of variables to estimate can 
be defined in `Parms2Estimate`. The value set for parameters which are not estimated is then defined 
with the `InitValues` input variable. 

```r
# Estimate all parameters using PROSPECT inversion applied to full spectral data
res_all_WL <- Invert_PROSPECT(SpecPROSPECT = SubSpecPROSPECT, 
                              Refl = SubRefl, 
                              Tran = SubTran)
```

## PROSPECT inversion using the optimal spectral domains specific to each constituent

The function `Invert_PROSPECT_OPT` performs PROSPECT inversion using optimal spctral domains
defined in [@feret2019; @spafford2021]. The optimal spectral domains are specific to each constituent, 
and the function automatically adjusts the spectral domain of the leaf optical properties provided 
by user. 

```r
# Estimate a set of parameters using PROSPECT inversion with optimal spectral domains
Parms2Estimate  <- c('CHL', 'CAR', 'EWT', 'LMA')
res_opt_WL <- Invert_PROSPECT_OPT(SpecPROSPECT = SpecPROSPECT, 
                                  lambda = lambda, 
                                  Refl = Refl, 
                                  Tran = Tran,
                                  Parms2Estimate = Parms2Estimate)
```

## Comparing performances of the two types of inversion with experimental data

The outputs of the inversion can then be plotted and compared with scatterplots.

```r
# define axis labels for each leaf chemical constituent
Labs <- list(CHL = list(x='Estimated CHL (µg/cm-2)', y= 'Measured CHL (µg/cm-2)'), 
             CAR = list(x='Estimated CAR (µg/cm-2)', y= 'Measured CAR (µg/cm-2)'), 
             EWT = list(x='Estimated EWT (mg/cm-2)', y= 'Measured EWT (mg/cm-2)'), 
             LMA = list(x='Estimated LMA (mg/cm-2)', y= 'Measured LMA (mg/cm-2)'))
# define colors
Colors <- list(CHL = c('gray69','green4'), CAR = c('gray69','orange2'),
               EWT = c('gray69','blue2'), LMA = c('gray69','red3'))
# multiplying factor to convert unit
fact <- list('CHL' = 1, 'CAR' = 1, 'EWT' = 1000, 'LMA' = 1000)
plotBP <- list()
for (parm in Parms2Estimate){
  Est_all_WL <- data.frame('measured'= fact[[parm]]*DataBioch[[parm]],
                           'estimated'= fact[[parm]]*res_all_WL[[parm]],
                           'config' = 'full WL')
  EST_opt_WL <- data.frame('measured'= fact[[parm]]*DataBioch[[parm]],
                           'estimated'= fact[[parm]]*res_opt_WL[[parm]],
                           'config' = 'opt WL')
  Meas_Est <- rbind(Est_all_WL, EST_opt_WL)
  plotBP[[parm]] <- plotinv(Meas_Est, stats = TRUE, 
                            Labs = Labs[[parm]], Colors = Colors[[parm]])
}

plotALL <- gridExtra::grid.arrange(plotBP$CHL, plotBP$CAR, 
                                   plotBP$EWT, plotBP$LMA, 
                                   ncol = 2, nrow = 2)
ggsave(filename = 'PROSPECT_Inversions.png', plot = plotALL, device = "png", 
       scale = 1, width = 24, height = 24, units = "cm", dpi = 300)
```


![Estimation of chlorophyll content, carotenoid content, EWT and LMA from PROSPECT inversion applied on the ANGERS dataset. full WL corresponds to the inversion performed with the full spectral information; opt WL corresponds to the inversion performed with the optimal spectral information](PROSPECT_Inversions.png){ width=80% }

# Conclusion

We have described `prospect`, an R package dedicated to the PROSPECT leaf model. 
`prospect` can run different versions of the model in direct mode, in order to simulate 
directional-hemispherical reflectance and transmittance for individual leaves or LUTs.
`prospect` also include a selection of inversion procedures based on iterative optimization, 
in order to estimate leaf structure and chemical constituent content either from 
directional-hemispherical reflectance and transmittance, or from reflectance or transmittance only. 

We described these procedures in examples, and used experimental data to illustrate inversion procedures. 
The `prospect` package aims at providing latest advances in terms of model version and inversion procedures
to the leaf spectroscopy community. 


# Acknowledgements

The authors acknowledge financial support from Agence Nationale de la Recherche (BioCop project — ANR-17-CE32-0001)
We are grateful to Stephane Jacquemoud and Frederic Baret for the development of the initial version of 
the PROSPECT model. We also warmly thank Luc Bidel, Christophe François and Gabriel Pavan who collected the ANGERS dataset.

# References
