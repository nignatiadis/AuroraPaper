# README

This repository provides code to reproduce the results of the following paper:


> Ignatiadis, N., Saha, S., Sun D. L., & Muralidharan, O. (2019).  **Empirical Bayes mean estimation with nonparametric errors via order statistic regression.** [[arXiv]](https://arxiv.org/abs/1911.05970)


The main method has been implemented in Julia and is available as the [Aurora.jl](https://github.com/nignatiadis/Aurora.jl) Julia package.
In the subdirectory `AuroraR` of this repository we provide a R package that wraps the Julia package and includes code for the different methods and simulations in the paper.

The R package may be installed as follows:

```r
devtools::install_github("nignatiadis/AuroraPaper", subdir="AuroraR")
```

The folder `simulation_scripts` contains the code for the simulation studies (that calls functions from the `AuroraR` package above). Concretely:
* `simulation_scripts/homoskedastic_simulations.R` runs the simulations of Section 6.1.
* `simulation_scripts/heteroskedastic_simulations.R` runs the simulations of Section 6.2.
* `simulation_scripts/pareto_simulations.R` runs the simulations of Section 6.3.

The `vignettes` folder contains R Markdown files that reproduce the figures from the paper, some of which require files that are generated from the three previous scripts. 

* `vignettes/motivation.Rmd` reproduces Figure 1.
* `vignettes/homoskedastic_simulations_plots.Rmd` reproduces Figure 2. 
* `vignettes/location_family_auroral_coefficients.Rmd` reproduces Figure 3.
* `vignettes/heteroskedastic_simulations_plots.Rmd` reproduces Figure 4.
* `vignettes/pareto_simulations_plots.Rmd` reproduces Figure 5.
