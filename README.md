Using Spatial Capture-Recapture Models to Estimate Abundances, Densities, and Interspecific Associations in a Carnivore Community
================

Notes
-----

This repository contains data and code from:

Forsyth, D.M., Ramsey, D.S.L., Woodford, L. (2018). "Using Spatial Capture-Recapture Models to Estimate Abundances, Densities, and Interspecific Associations in a Carnivore Community" *Journal of Wildlife Management*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1459045.svg)](https://doi.org/10.5281/zenodo.1459045)

Getting started
---------------

**File descriptions**:

*SPA nimble.R* – R code used to conduct density estimation for foxes using the spatial unmarked model of Ramsey et al 2015. Uses functions defined in *Nimble functions.R*

*SMR nimble.R* – R functions used to conduct density estimation for dingos and feral cats using spatial mark-resight models.Uses functions defined in *Nimble functions.R* and *SMR functions.r*.

*Gudgenby interaction analysis.R* – R code used to estimate the spatial association among the estimated activity centres for foxes, cats and dogs. Uses estimates from the posterior distributions of activity centres provided in Results/.

Prerequisites
-------------

These scripts require packages *nimble*, *spatstat*, *survey*, *readxl*. *MCMCvis*, and *ggmcmc*.
