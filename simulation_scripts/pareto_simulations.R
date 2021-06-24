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

methods_list_pareto <- list( auroral = function(Zs, sigmas) auroral(Zs),
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

bla <- evaluate_single_pareto_sim()


pareto_alphas <- 3
unif_ubs <- seq(2.5, to=4, length=10)
nreps <- 100

pareto_prelim_param <- expand.grid(n=c(10000,100000),
                                   B=c(20,100))

n = pareto_prelim_param$n[arg1]
B = pareto_prelim_param$B[arg1]

pareto_param_df <- expand.grid( pareto_alpha=rep(pareto_alphas, each=nreps), #force correct number of Monte Carlo Replicates
                                unif_ub=unif_ubs,
                                n = n,
                                B = B
                                )
pareto_res <- bind_rows(lapply(1:nrow(pareto_param_df),
                                 function(i) {pareto_alpha <- pareto_param_df$pareto_alpha[i];
                                 unif_ub <- pareto_param_df$unif_ub[i];
                                 n <- pareto_param_df$n[i];
                                 B <- pareto_param_df$B[i];
                                 evaluate_single_pareto_sim(n=n,
                                                            B=B,
                                                            unif_lb=2,
                                                            unif_ub=unif_ub,
                                                            pareto_alpha=pareto_alpha)}))

saveRDS(pareto_res,
        file.path("files", paste0(paste("pareto", n, B, sep="_"),".rds")))
