#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

arg1 <- as.integer(args[1])
arg2 <- 1 + (as.integer(args[2]) -1 ) %% 10


print(arg1)
print(arg2)

library(AuroraR)
library(tidyverse)


julia_dir <- "/home/users/ignat/julia/julia-1.6.1/bin"

aurora_setup(JULIA_HOME=julia_dir)

methods_list_pareto <- list( auroral = function(Zs, sigmas) auroral(Zs),
                             auroraKNN = function(Zs, sigmas) aurora_knn(Zs, kKNN=100),
                             ccl = function(Zs, sigmas) ebcc(Zs),
                             location_mean = function(Zs, sigmas) location_mean(Zs),
                             location_median = function(Zs, sigmas) location_median(Zs),
                             pareto_mle = function(Zs, sigmas) pareto_mle(Zs)
                            )




evaluate_single_pareto_sim <- function(n=10000, B=10, pareto_alpha=3, unif_lb=0.5, unif_ub=10){
  tmp_sim <- uniform_pareto_sim(n, B, pareto_alpha=pareto_alpha, unif_lb=unif_lb, unif_ub=unif_ub)
  tmp_error <- sapply(methods_list_pareto, function(f) mean( (f(tmp_sim$Zs, tmp_sim$sigmas) - tmp_sim$true_mus)^2))
  tmp_df <- tibble( method=names(methods_list_pareto),
                    error=tmp_error,
                    B=B,
                    n=n,
                    pareto_alpha=pareto_alpha,
                    unif_lb=unif_lb,
                    unif_ub=unif_ub,
                    N_B_string=paste0("K=",B,", N=",n))
  tmp_df
}


pareto_alpha <- 3
unif_ubs <- seq(2.5, to=4, length=10)
nreps <- 100

n <- switch(arg1, 10000, 10000, 100000)
B <- switch(arg1, 20, 100, 100)

unif_ub <- unif_ubs[arg2]

for (i in 1:nreps){
  pareto_res <- evaluate_single_pareto_sim(n=n,
                                           B=B,
                                           unif_lb=2,
                                           unif_ub=unif_ub,
                                           pareto_alpha=pareto_alpha)

  rand_string = paste0(sample(LETTERS, 5, TRUE), collapse="")

  saveRDS(pareto_res,
        file.path("files", paste0(paste("pareto", n, B, rand_string, sep="_"),".rds")))

}