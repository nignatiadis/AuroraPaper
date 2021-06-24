#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

arg1 <- as.integer(args[1])
julia_dir <- args[2]

print(arg1)
print(julia_dir)

library(AuroraR)
library(tidyverse)


#julia_dir = "/Applications/Julia-1.5.app/Contents/Resources/julia/bin"

aurora_setup(JULIA_HOME=julia_dir)


methods_list_homosk <- list( auroral = function(Zs, sigmas) auroral(Zs),
                             auroraKNN = function(Zs, sigmas) aurora_knn(Zs),
                             location_mean = function(Zs, sigmas) location_mean(Zs),
                             location_median = function(Zs, sigmas) location_median(Zs),
                             location_midrange = function(Zs, sigmas) location_midrange(Zs),
                             ccl = function(Zs, sigmas) ebcc(Zs),
                             james_stein = function(Zs, sigmas) james_stein(Zs, sigma=mean(sigmas)),
                             gmleb_oracle = function(Zs, sigmas) gmleb(Zs, sigma=sigmas),
                             gmleb = function(Zs, sigmas) gmleb(Zs, sigma=apply(Zs, 1 , sd)/sqrt(ncol(Zs)))
                            )


evaluate_single_homoskedastic_sim <- function(A, Likelihood, Prior="Normal", n=10000, B=10, sigma=sqrt(0.4)){
  tmp_sim <- switch(Prior,
                "Normal"=switch(Likelihood,
                    "Normal"=normal_normal_sim(n, B, sigma, A, prior_mu=0.5),
                    "Laplace"=normal_laplace_sim(n, B, sigma, A, prior_mu=0.5),
                    "Rectangular"=normal_uniform_sim(n, B, sigma, A, prior_mu=0.5)
                    ),
                "ThreeComponent"=switch(Likelihood,
                                "Normal"=three_component_normal_sim(n, B, sigma, A),
                                "Laplace"=three_component_laplace_sim(n, B, sigma, A),
                                "Rectangular"=three_component_uniform_sim(n, B, sigma, A)
                    )
                )
  tmp_error <- sapply(methods_list_homosk, function(f) mean( (f(tmp_sim$Zs, tmp_sim$sigma) - tmp_sim$true_mus)^2))
  tmp_df <- tibble( method=names(methods_list_homosk),
                    error = tmp_error,
                    A=A,
                    Prior=Prior,
                    Likelihood = Likelihood)
  tmp_df
}

#---------------------------------------------
# Simulation Options
#---------------------------------------------
set.seed(1)
B <- 10
n <- 10000
sigma <- sqrt(0.4)
nreps <- 100
As <- seq(0.4, to=2.5, by=0.2)^2
#---------------------------------------------

homoskedastic_setting_names <- c("Normal", "Laplace", "Rectangular")
prior_setting_names <- c("Normal","ThreeComponent")

homoskedastic_param_df <- expand.grid(Likelihood=homoskedastic_setting_names, Prior=prior_setting_names)

Likelihood <- as.character(homoskedastic_param_df$Likelihood[arg1])
Prior <- as.character(homoskedastic_param_df$Prior[arg1])
setting_df <- data.frame(Likelihood = Likelihood,
                          Prior = Prior,
                          A = rep(As, nreps))

homoskedastic_res <- bind_rows(lapply(1:nrow(setting_df),
                              function(i) evaluate_single_homoskedastic_sim(setting_df$A[i],
                                                                            setting_df$Likelihood[i],
                                                                            Prior=setting_df$Prior[i],
                                                                            n=n, B=B, sigma=sigma)
                              )
                            )

saveRDS(homoskedastic_res,
       file.path("files", paste0(paste("homoskedastic", Prior, Likelihood, sep="_"),".rds")))
