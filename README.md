# README

This repository provides code to reproduce the results of the following paper:

> Nikolaos Ignatiadis, Sujayam Saha, Dennis L. Sun & Omkar Muralidharan (2021) **Empirical Bayes Mean Estimation With Nonparametric Errors Via Order Statistic Regression on Replicated Data**, Journal of the American Statistical Association, DOI: 10.1080/01621459.2021.1967164

The paper is also available on arXiv:
> Nikolaos Ignatiadis, Sujayam Saha, Dennis L. Sun & Omkar Muralidharan (2021).  **Empirical Bayes mean estimation with nonparametric errors via order statistic regression on replicated data.** [[arXiv]](https://arxiv.org/abs/1911.05970)


The main method has been implemented in Julia and is available as the [Aurora.jl](https://github.com/nignatiadis/Aurora.jl) Julia package.

## AuroraR
[![R build status](https://github.com/nignatiadis/AuroraPaper/workflows/R-CMD-check/badge.svg)](https://github.com/nignatiadis/AuroraPaper/actions) 
[![codecov](https://codecov.io/gh/nignatiadis/AuroraPaper/branch/main/graph/badge.svg?token=j8JBrFe6Ks)](https://codecov.io/gh/nignatiadis/AuroraPaper)


In the subdirectory `AuroraR` of this repository we provide a R package that wraps the Julia package and includes code for the different methods and simulations in the paper.

The R package may be installed as follows:

```r
devtools::install_github("nignatiadis/AuroraPaper", subdir="AuroraR")
```
We note that the R package also wraps the nonparametric maximum likelihood (NPMLE) functionality from the [REBayes](https://cran.r-project.org/web/packages/REBayes/index.html) package. In turn, REBayes requires a working installation of the [Mosek](https://www.mosek.com/) convex optimization solver; we used Version 9.2.

## Reproduction code

The folder `simulation_scripts` contains the code for the simulation studies (that calls functions from the `AuroraR` package above). Concretely:

* `simulation_scripts/homoskedastic_simulations.R` runs the simulations of Section 6.1. It can be called in the terminal via `R homoskedastic_simulations.R arg1 dir`, where `arg1` can take on values 1,2,.., 6 (corresponding to different combinations of prior/likelihood) and `dir` is the directory in which Julia is installed.
* `simulation_scripts/heteroskedastic_simulations.R` runs the simulations of Section 6.2. It can be called in the terminal via `R heteroskedastic_simulations.R arg1 dir` where `arg1` can take values 1,2 or 3 (for the three simulation settings considered) and `dir` is the Julia directory.
* `simulation_scripts/pareto_simulations.R` runs the simulations of Section 6.3. (Warning: the case with 100,000 units and 100 replicates is slow and may take >20 hours per Monte Carlo replicate.)

The `vignettes` folder contains R Markdown files that reproduce the figures from the paper, some of which require files that are generated from the three previous scripts. 

* `vignettes/motivation.Rmd` reproduces Figure 1.
* `vignettes/homoskedastic_simulations_plots.Rmd` reproduces Figure 2. 
* `vignettes/location_family_auroral_coefficients.Rmd` reproduces Figure 3.
* `vignettes/heteroskedastic_simulations_plots.Rmd` reproduces Figure 4.
* `vignettes/pareto_simulations_plots.Rmd` reproduces Figure 5.
