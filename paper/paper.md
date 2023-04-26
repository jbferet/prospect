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

Leaf optical properties are directly linked to their biochemical content and structural 
properties through light absorption and scattering. They are also an important driver 
of vegetation reflectance acquired from Earth observation systems.
Hence, physical modeling of leaf optical properties is key to better understand 
and monitor ecosystems and agrosystems.
PROSPECT is the most widely used leaf model, allowing simulation of leaf optical properties 
(directional-hemispherical reflectance and transmittance) based on a limited set of chemical 
constituents with specific absorption coefficients and a unique leaf structure parameter 
to account for scattering.

We present `prospect`, an R package which contains
the most recent versions of the model, including PROSPECT-D [@feret2017] and PROSPECT-PRO [@feret2021].
The package also includes inversion routines aiming at estimating leaf chemistry 
directly from leaf optical properties, either leaf reflectance or leaf transmittance, or both.
These inversion routines, and are based on iterative optimization. Various inversion strategies
are implemented in the package, including strategies recently developed [@feret2019, @spafford2021] 
aiming at estimating leaf chemical constituents based on optimal spectral domains, and 
estimating leaf structure when only reflectance or transmittance is used.

# Statement of need

The estimation of vegetation biophysical and chemical properties is 
becoming increasingly important for the monitoring of ecosystem and agrosystem functions. 
A set of key biochemical traits, including leaf pigment content, water content, nitrogen content
and leaf mass per area (LMA) absorb light in specific domains. 
These constituents also influence canopy reflectance measured from Earth observation systems. 
In situ data collection at leaf scale is critical to provide ground truth when upscaling 
and mapping vegetation traits from satellite, airborne or UAV imagery. 

Techniques based on leaf spectroscopy have developed in the past decades in order to provide a 
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
leaf development stages, from juvenile to mature and senescent. The latest version named 
PROSPECT-PRO included proteins and carbon based constituents as complementary fractions of the 
constituent identified as dry matter (corresponding to LMA). In parallel with the 
development of new versions of the model, advances in model inversion were also performed 
to improve the estimation of leaf chemical constituents [@feret2019; @spafford2021].

The application of an appropriate model inversion strategy is critical to ensure optimal 
estimation of leaf constituents, as it was evidenced for the estimation of LMA [@feret2019], 
which was reported as poorly retrieved from PROSPECT inversion in earlier studies. 

To ensure access to latest advances in terms of model version and inversion strategies, 
we present `prospect`, an R package including the most recent versions of PROSPECT, 
and parameterizable inversion routines, allowing users to design their own inversion strategy.


# Overview

Two versions of the model PROSPECT are implemented in `prospect`: PROSPECT-D [@feret2017] and 
PROSPECT-PRO [@feret2021]. [@feret2017] highly recommended using PROSPECT-D instead of the previous 
version PROSPECT-5, even for vegetation without anthocyanins. Nevertheless, the function `PROSPECT`
allows specifying version `5`. This version is actually a version of PROSPECT-D with anthocyanin 
content set to 0.
For each version, it is possible to also include the influence brown pigments (BROWN), which appear during senescence, 
by ending the name of the version with `B`. Otherwise, `prospect` assumes that leaves contain no brown pigments.

(Table \ref{table:1} provides an overview of the leaf chemical constituents included in the different versions which can be specified when calling the function `PROSPECT`.

| Version 	|  `5`  	|  `5B`  	|  `D`  	|  `DB`  	|  `PRO`  	|  `PROB`  	|
|---------	|--------	|:------:	|:------:	|:------:	|:------:	|:------:	|
| CHL     	|  **X**  	|  **X**  	|  **X**  	|  **X**  	|  **X**  	|  **X**  	|
| CAR   	|  **X**  	|  **X**  	|  **X**  	|  **X**  	|  **X**  	|  **X**  	|
| ANT   	|        	|       	|  **X**  	|  **X**  	|  **X**  	|  **X**  	|
| BROWN 	|       	|  **X**  	|        	|  **X**  	|          	|  **X**  	|
| EWT   	|  **X**  	|  **X**  	|  **X**  	|  **X**  	|  **X**  	|  **X**  	|
| LMA   	|  **X**  	|  **X**  	|  **X**  	|  **X**  	|         	|          	|
| PROT  	|         	|          	|         	|         	|  **X**  	|  **X**  	|
| CBC   	|          	|        	|         	|         	|  **X**  	|  **X**  	|

: Versions of the model PROSPECT  accepted as input of the `PROSPECT` function and 
corresponding chemical constituents.\label{table:1}

The inversion procedure is based on the function `fmincon` implemented in the package `pracma`. 
`fmincon` is an optimization algorithm aiming at finding the minimum of multivariable 
functions with nonlinear constraints. 


(Table \ref{table:2} provides information on the spectral range used to estimate leaf chemical constituents from 
their optical properties, for the two main inversion strategies.

| Constituent 	|    Optimal spectral domain    	|                 Versions              	|
|------------	|---------------------------------	|:-------------------------------------:	|
| CHL       	|           700 -- 720           	| `5` | `5B` | `D` | `DB` | `PRO` | `PROB`	|
| CAR       	|           520 -- 560           	| `5` | `5B` | `D` | `DB` | `PRO` | `PROB`	|
| ANT       	|           400 -- 800           	|  |  | `D` | `DB` | `PRO` | `PROB`	|
| BROWN     	|               NA              	| NA                                	|
| EWT       	|          1700 -- 2400           	| `5` | `5B` | `D` | `DB` | `PRO` | `PROB`	|
| LMA       	|          1700 -- 2400           	| `5` | `5B` | `D` | `DB` |  | 	|
| PROT      	|   2100 -- 2139 ; 2160 -- 2179   	|  |  |  |  | `PRO` | `PROB`	|
| CBC       	|   2100 -- 2139 ; 2160 -- 2179   	|  |  |  |  | `PRO` | `PROB`	|


: Versions of the model PROSPECT  accepted as input of the `PROSPECT` function and 
corresponding chemical constituents.\label{table:3}


# Example 1: running PROSPECT in forward mode

## Individual simulation with PROSPECT-D and PROSPECT-PRO

This section illustrates how to run two different versions of PROSPECT in forward mode in order to simulate 
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

In this second example, as proteins and CBC are defined, PROSPECT-PRO is called. If user also defines a non-null 
value for LMA in addition to PROT and CBC, the value of LMA is automatically set to 0, as PROSPECT-PRO does not 
use LMA as input value. A warning is then displayed. 

```r
# Load prospect package
library(prospect)

# Run PROSPECT-PRO
LRT_PRO <- PROSPECT(SpecPROSPECT, CHL = 45, CAR = 10, ANT = 0.2, 
 EWT = 0.012, PROT = 0.001,  CBC = 0.009, N = 1.7)
```

The leaf optical properties can then be compared.

```r
R_D_df <- data.frame('wvl' = LRT_D$wvl,
  'RT' = 100*LRT_D$Reflectance,
  'LOP' = matrix('R PROSPECT-D',nrow = length(LRT_D$Reflectance),ncol = 1))

T_D_df <- data.frame('wvl' = LRT_D$wvl,
  'RT' = 100*(1-LRT_D$Transmittance),
  'LOP' = matrix('T PROSPECT-D',nrow = length(LRT_D$Transmittance),ncol = 1))

R_PRO_df <- data.frame('wvl' = LRT_PRO$wvl,
  'RT' = 100*LRT_PRO$Reflectance,
  'LOP' = matrix('R PROSPECT-PRO',nrow = length(LRT_PRO$Reflectance),ncol = 1))

T_PRO_df <- data.frame('wvl' = LRT_PRO$wvl,
  'RT' = 100*(1-LRT_PRO$Transmittance),
  'LOP' = matrix('T PROSPECT-PRO',nrow = length(LRT_PRO$Transmittance),ncol = 1))

LRT_df <- rbind(R_D_df, T_D_df, R_PRO_df, T_PRO_df)

library(ggplot2)
plot <- ggplot2::ggplot(LRT_df, aes(x=wvl, y=RT, group=LOP)) +
  geom_line(aes(linetype=LOP, color=LOP),linewidth=1.00)+
  scale_color_manual(values=c('#FF9999','red1','#9999FF','blue4'))+
  scale_size_manual(values=c(5, 5))+
  labs(x='Wavelength (nm)',y='Reflectance        (%)        1-Transmittance') +
  theme(legend.position="bottom",
        axis.text=element_text(size=12),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) +
  scale_linetype_manual(values=c("solid","solid","solid","solid"))

filename = 'compare_RT_PROSPECT_PRO_D.png'
ggsave(filename, plot = last_plot(), device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)
```

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

Here, the iterative optimization is based on the minimization of a nonlinear constrained multivariable function. 
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
# download ANGERS leaf optics database from gitlab repository
dbName <- 'ANGERS'
gitlab_Rep <- 'https://gitlab.com/jbferet/myshareddata/raw/master/LOP'

# download leaf biochemical constituents and leaf optical properties
fileName <- list('DataBioch.txt','ReflectanceData.txt','TransmittanceData.txt')
DataBioch <- fread(file.path(gitlab_Rep,dbName,fileName[[1]]))
Refl<- fread(file.path(gitlab_Rep,dbName,fileName[[2]]))
Tran <- fread(file.path(gitlab_Rep,dbName,fileName[[3]]))

# Get the wavelengths corresponding to reflectance and transmittance measurements  
lambda <- unlist(Refl$V1, use.names=FALSE)
Refl$V1 <- Tran$V1 <- NULL
```

## PROSPECT inversion using the full spectral range available

The spectral domain covered by `SpecPROSPECT` and the leaf optical properties are expected to match 
when performing inversion on the full spectral domain covered by the data.
This can be done automatically using `FitSpectralData`

```r
SubData <- FitSpectralData(SpecPROSPECT = SpecPROSPECT,      # spectral properties over 400-2500 nm
        lambda = lambda,                  # spectral bands corresponding to refl and Tran
        Refl = Refl,   # Reflectance data (matrix or dataframe)
        Tran = Tran,   # Transmittance data (matrix or dataframe)
        UserDomain = lambda)              # spectral domain of interest
```

The main inversion procedure is called with the function `Invert_PROSPECT`.
It is based on the minimization of a merit function. The default merit function is based on the 
root mean squared difference between measured and simulated leaf optical properties, 




```r
```

## PROSPECT inversion using the optimal spectral domains specific to each constituent

```r
```



# Conclusion

We have described `prospect`, an R package dedicated to the PROSPECT leaf model. 
`prospect` can run different versions of the model in direct mode, in order to simulate 
directional-hemispherical reflectance and transmittance or LUTs. 
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
