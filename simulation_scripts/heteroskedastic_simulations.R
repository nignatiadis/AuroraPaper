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



methods_list_heterosk <- list( auroral = function(Zs, sigmas) auroral(Zs),
                               auroraKNN = function(Zs, sigmas) aurora_knn(Zs),
                               location_mean = function(Zs, sigmas) location_mean(Zs),
                               ccl = function(Zs, sigmas) ebcc(Zs),
                               gmleb_oracle = function(Zs, sigmas) gmleb(Zs, sigma=sigmas),
                               gmleb = function(Zs, sigmas) gmleb(Zs, sigma=apply(Zs, 1 , sd)/sqrt(ncol(Zs))),
                               grouplinear_oracle = function(Zs, sigmas) grouplinear(Zs, sigma=sigmas),
                               grouplinear = function(Zs, sigmas) grouplinear(Zs, sigma=apply(Zs, 1 , sd)/sqrt(ncol(Zs))),
                               sure_oracle = function(Zs, sigmas) sure_grandmean(Zs, sigma=sigmas),
                               sure = function(Zs, sigmas) sure_grandmean(Zs, sigma=apply(Zs, 1 , sd)/sqrt(ncol(Zs)))
                              )


evaluate_single_heterosk_sim <- function(Likelihood="WMBZ_C", n=10000, B=10,
                                                sigma_squared_lower=0.1, sigma_squared_upper=1){
  tmp_sim <- switch(Likelihood,
                    "WMBZ_A" = wmbz_example_a_sim(n, B, sigma_squared_lower, sigma_squared_upper, A=0.5),
                    "WMBZ_C" = wmbz_example_c_sim(n, B, sigma_squared_lower, sigma_squared_upper),
                    "WMBZ_F" = wmbz_example_f_sim(n, B, sigma_squared_lower, sigma_squared_upper))

  tmp_error <- sapply(methods_list_heterosk, function(f) mean( (f(tmp_sim$Zs, tmp_sim$sigmas) - tmp_sim$true_mus)^2))
  tmp_df <- tibble( method=names(methods_list_heterosk),
                    error = tmp_error,
                    sigma_squared_lower = sigma_squared_lower,
                    n = n, 
                    B = B,
                    sigma_squared_upper = sigma_squared_upper,
                    Likelihood = Likelihood)
  tmp_df
}

heteroskedastic_setting_names <- c("WMBZ_A", "WMBZ_C", "WMBZ_F")
setting_name <- heteroskedastic_setting_names[arg1]

sigma_squared_lower <- 0.1
sigma_squared_upper_list <- seq(1, to=3, length=10)
nreps <- 100

heteroskedastic_param_df <-tibble( sigma_squared_upper=rep(sigma_squared_upper_list, each=nreps), Likelihood=setting_name)

heteroskedastic_res <- bind_rows(lapply(1:nrow(heteroskedastic_param_df),
                                          function(i) {lik <- heteroskedastic_param_df$Likelihood[i];
                                          sigma_squared_upper <- heteroskedastic_param_df$sigma_squared_upper[i];
                                          evaluate_single_heterosk_sim(n=10000,
                                                                              B=10,
                                                                              Likelihood=setting_name,
                                                                              sigma_squared_lower=0.1,
                                                                              sigma_squared_upper=sigma_squared_upper)}

                                        )
                                 )

saveRDS(heteroskedastic_res,
        file.path("files", paste0(paste("heteroskedastic", setting_name, sep="_"),".rds")))


