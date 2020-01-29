# __prospect__ <img src="man/figures/logo.png" align="right" alt="" width="200" />

# An R package for the simulation of leaf optical properties based on their biochemical and biophysical properties using the PROSPECT leaf model. 

[![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg)](https://www.r-project.org/Licenses/GPL-3)
[![Build Status](https://gitlab.com/jbferet/prospect/badges/master/pipeline.svg)](https://gitlab.com/jbferet/prospect/pipelines/latest)

# 1 Install

After installing package `devtools`, the package `prospect` can be installed with the following command line in R session:
```
devtools::install_gitlab('jbferet/prospect')
```

... if you are already on this webpage, but prospect is still not publicly available... Lucky you!!!
then install the `getPass` package, and run this command line:

```
devtools::install_git('https://gitlab.com/jbferet/prospect',credentials = git2r::cred_user_pass('Your_Gitlab_UserName',getPass::getPass())) 
```

# 2 Tutorial

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- ```{r include = FALSE} -->
<!-- knitr::opts_chunk$set( -->
<!--   collapse = TRUE, -->
<!--   comment = "#>", -->
<!--   fig.path = "man/figures/README-", -->
<!--   out.width = "100%" -->
<!-- ) -->
<!-- ``` -->

A tutorial vignette is available [here](https://jbferet.gitlab.io/prospect/articles/prospect.html).

# 3 Citation

If you use **prospect**, please cite the following references:

Féret, J.-B., Berger, K., de Boissieu, F. & Malenovský 2020. Estimation of leaf protein and carbon-based constituent content from optical properties with the PROSPECT-PRO model. Remote Sensing of Environment.

Féret, J.-B., Gitelson, A.A., Noble, S.D. & Jacquemoud, S. 2017. PROSPECT-D: Towards modeling leaf optical properties through a complete lifecycle. Remote Sensing of Environment. 193, 204–215. http://dx.doi.org/10.1016/j.rse.2017.03.004
