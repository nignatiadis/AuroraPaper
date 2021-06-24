
#-------------------------
# Simulation functions
#--------------------------
# Below we implement some basic simulation functions for homoskedastic location families. T
# They are named as {prior distribution}_{sampling distribution}_sim

#' Homoskedastic location family simulations with Normal prior
#'
#' @param m Integer, number of units.
#' @param B Integer, number of replicates.
#' @param sigma Numeric, noise standard deviation (of replicate average, i.e.,
#'     each noisy measurements has a standard deviation equal to sqrt(B)*sigma).
#' @param A Numeric, prior variance
#' @param prior_mu Numeric, mean of the prior (defaults to 0).
#' @template simulation_output
#' @name homoskedasticloc
NULL

#' @describeIn homoskedasticloc Normal prior, Normal likelihood
#' @export
normal_normal_sim <- function(m, B, sigma, A, prior_mu=0){
  mus <- stats::rnorm(m, prior_mu, sqrt(A))
  Zs <- matrix( stats::rnorm(m*B, 0, sigma*sqrt(B)),  nrow=m) + mus
  res <- list(true_mus = mus, Zs = Zs, m=m, sigma=sigma, B=B, A=A, prior_mu=prior_mu)
  return(res)
}

#' @describeIn homoskedasticloc Normal prior, Laplace likelihood
#' @export
normal_laplace_sim <- function(m, B, sigma, A, prior_mu=0){
  mus <- stats::rnorm(m, prior_mu, sqrt(A))
  Zs <- matrix( stats::rexp(m*B, 1/sigma/sqrt(B/2)) *(2*stats::rbinom(m*B, 1, 0.5)-1),  nrow=m) + mus
  res <- list(true_mus = mus, Zs = Zs, m=m, sigma=sigma, B=B, A=A, prior_mu=prior_mu)
  res
}

#' @describeIn homoskedasticloc Normal prior, Rectangular (Uniform) likelihood
#' @export
normal_uniform_sim <- function(m, B, sigma, A, prior_mu=0){
  mus <- stats::rnorm(m, prior_mu, sqrt(A))
  # For U[-M,M], sigma^2 = M^2/3
  Zs <- matrix( stats::runif(m*B, -sigma*sqrt(B*3), + sigma*sqrt(B*3)),  nrow=m) + mus
  res <- list(true_mus = mus, Zs = Zs, m=m, sigma=sigma, B=B, A=A, prior_mu=prior_mu)
  res
}



#' Homoskedastic location family simulations with Three component prior
#'
#' Here the prior assigns equal probability to -c, 0, +c, where c is specified
#' through the parameter A.
#'
#' @param m Integer, number of units.
#' @param B Integer, number of replicates.
#' @param sigma Numeric, noise standard deviation (of replicate average, i.e.,
#'     each noisy measurements has a standard deviation equal to sqrt(B)*sigma).
#' @param A Numeric, prior variance (defaults to 9). Note that the c on which
#'    G puts its mass is equal to sqrt(1.5*A).
#' @template simulation_output
#' @name threecomponentloc
NULL

#' @describeIn threecomponentloc Three component prior, Normal likelihood
#' @export
three_component_normal_sim <- function(m, B, sigma, A = 9){
  effect_sizes <- sqrt(1.5*A)*(-1:1)
  mus <- sample(effect_sizes, size=m, replace=TRUE)
  Zs <- matrix( stats::rnorm(m*B, 0, sigma*sqrt(B)),  nrow=m) + mus
  res <- list(true_mus = mus, Zs = Zs, m=m, sigma=sigma, B=B, A=A)
  return(res)
}

#' @describeIn threecomponentloc Three component prior, Laplace likelihood
#' @export
three_component_laplace_sim <- function(m, B, sigma, A = 9){
  effect_sizes <- sqrt(1.5*A)*(-1:1)
  mus <- sample(effect_sizes, size=m, replace=TRUE)
  Zs <- matrix( stats::rexp(m*B, 1/sigma/sqrt(B/2)) *(2*stats::rbinom(m*B, 1, 0.5)-1),  nrow=m) + mus
  res <- list(true_mus = mus, Zs = Zs, m=m, sigma=sigma, B=B, A=A)
  res
}

#' @describeIn threecomponentloc Three component prior, Rectangular (Uniform) likelihood
#' @export
three_component_uniform_sim <- function(m, B, sigma, A = 9){
  effect_sizes <- sqrt(1.5*A)*(-1:1)
  mus <- sample(effect_sizes, size=m, replace=TRUE)
  # For U[-M,M], sigma^2 = M^2/3
  Zs <- matrix( stats::runif(m*B, -sigma*sqrt(B*3), + sigma*sqrt(B*3)),  nrow=m) + mus
  res <- list(true_mus = mus, Zs = Zs, m=m, sigma=sigma, B=B, A=A)
  res
}





#' Heteroskedastic location family simulations
#'
#' The variance of the replicates of unit i vary across units (randomly), as follows:
#'
#'  sigma_i^2 ~ U[sigma_squared_lower, sigma_squared_upper]/B,
#'
#' where B is the number of replicates per unit and
#' sigma_squared_lower, sigma_squared_upper are simulation parameters.
#'
#' Note that sigma_i^2 is the variance of the replicate average, i.e.,
#' i.e.,  each noisy measurements has a standard deviation equal to sqrt(B)*sigma_i)
#'
#' @param m Integer, number of units.
#' @param B Integer, number of replicates.
#' @param sigma_squared_lower Numeric, lower limit of noise variance.
#' @param sigma_squared_upper Numeric, upper limit of noise variance.
#' @template simulation_output
#' @references The simulations here in are adapted from the paper
#'       "Group-linear empirical Bayes estimates for a heteroscedastic normal mean"
#'        by A Weinstein, Z Ma, LD Brown, CH Zhang (JASA, 2018)
#' @name heteroskedasticloc

NULL

#' @describeIn heteroskedasticloc Here mu_i ~ N(0, A). Then Z_ij ~ N(mu_i, B*sigma_i^2).
#' @param A Numeric, Prior noise variance.
#' @export
wmbz_example_a_sim <- function(m, B, sigma_squared_lower, sigma_squared_upper, A=1){
  mus <- stats::rnorm(m, 0, sqrt(A))
  sigmas_squared <- sort(stats::runif(m, min=sigma_squared_lower, max=sigma_squared_upper))
  sigmas <- sqrt(sigmas_squared)
  Zs <- matrix( stats::rnorm(m*B, 0, sqrt(B)),  nrow=m)*sigmas + mus
  res <- list(true_mus = mus, Zs = Zs,
              sigmas = sigmas,
              m = m, B = B, A=A,
              sigma_squared_lower = sigma_squared_lower,
              sigma_squared_upper = sigma_squared_upper)
  return(res)
}



#' @describeIn heteroskedasticloc Here mu_i = sigma_i^2 deterministically. Then Z_ij ~ N(mu_i, B*sigma_i^2).
#' @export
wmbz_example_c_sim <- function(m, B, sigma_squared_lower, sigma_squared_upper){
  sigmas_squared <- sort(stats::runif(m, min=sigma_squared_lower, max=sigma_squared_upper))
  mus <- sigmas_squared
  sigmas <- sqrt(sigmas_squared)
  Zs <- matrix( stats::rnorm(m*B, 0, sqrt(B)),  nrow=m)*sigmas + mus
  res <- list(true_mus = mus, Zs = Zs, sigmas = sigmas,
              m = m, B = B,
              sigma_squared_lower = sigma_squared_lower,
              sigma_squared_upper = sigma_squared_upper)
  return(res)
}

#' @describeIn heteroskedasticloc Here mu_i = sigma_i^2 deterministically. Then Z_ij = mu_i + sqrt(3B)*sigma_i*U[-1, 1].
#' @export
wmbz_example_f_sim <- function(m, B, sigma_squared_lower, sigma_squared_upper){
  sigmas_squared <- sort(stats::runif(m, min=sigma_squared_lower, max=sigma_squared_upper))
  mus <- sigmas_squared
  sigmas <- sqrt(sigmas_squared)
  Zs <-  matrix( stats::runif(m*B, -sqrt(B*3), sqrt(B*3)),  nrow=m)*sigmas + mus
  res <- list(true_mus = mus, Zs = Zs, sigmas = sigmas,
              m = m, B = B,
              sigma_squared_lower = sigma_squared_lower,
              sigma_squared_upper = sigma_squared_upper)
  return(res)
}

#' Simulation: Pareto likelihood, Gamma Prior
#'
#' mu_i ~ Gamma(gamma_shape, gamma_scale), Z_ij ~ Pareto( mean = mu_i, alpha)
#'
#' @param m Integer, number of units.
#' @param B Integer, number of replicates.
#' @param pareto_alpha Numeric, tail index of Pareto distribution (defaults to 3, note
#'           that it needs to be >2 for the variance to exist).
#' @param gamma_shape Numeric, shape of the prior Gamma distribution (defaults to 1).
#' @param gamma_scale Numeric, scale of the prior Gamma distribution (defaults to 1).
#' @template simulation_output
#' @export
gamma_pareto_sim <- function(m, B, pareto_alpha=3, gamma_shape=1, gamma_scale=1){
  mus <- stats::rgamma(m, shape=gamma_shape, scale=gamma_scale)
  xs <- mus*(pareto_alpha - 1)/pareto_alpha
  exp_rvs <- matrix( stats::rexp(m*B, rate=pareto_alpha),  nrow=m)
  Zs <- exp(exp_rvs)*xs
  if (pareto_alpha > 2){
    sigmas_squared <- xs^2*pareto_alpha/((pareto_alpha-1)^2*(pareto_alpha-2))
    sigmas <- sqrt(sigmas_squared)
  } else {
    sigmas <- rep(m, Inf)
  }
  res <- list(true_mus = mus, Zs = Zs, sigmas=sigmas, B=B,
              m=m, pareto_alpha = pareto_alpha, shape=gamma_shape, scale=gamma_scale)
  res
}

#' Simulation: Pareto likelihood, Uniform Prior
#'
#' mu_i ~ Uniform[unif_lb, unif_ub], Z_ij ~ Pareto( mean = mu_i, alpha)
#'
#' @param m Integer, number of units.
#' @param B Integer, number of replicates.
#' @param pareto_alpha Numeric, tail index of Pareto distribution (defaults to 3, note
#'           that it needs to be >2 for the variance to exist).
#' @param unif_lb Numeric, lower support point of the prior Uniform distribution (defaults to 0.1).
#' @param unif_ub Numeric, upper support point of the prior Uniform distribution (defaults to 1).
#' @template simulation_output
#' @export
uniform_pareto_sim <- function(m, B, pareto_alpha=3, unif_lb=0.1, unif_ub=1){
  mus <- stats::runif(m, unif_lb, unif_ub)
  xs <- mus*(pareto_alpha - 1)/pareto_alpha
  exp_rvs <- matrix( stats::rexp(m*B, rate=pareto_alpha),  nrow=m)
  Zs <- exp(exp_rvs)*xs
  if (pareto_alpha > 2){
    sigmas_squared <- xs^2*pareto_alpha/((pareto_alpha-1)^2*(pareto_alpha-2))
    sigmas <- sqrt(sigmas_squared)
  } else {
    sigmas <- rep(m, Inf)
  }
  res <- list(true_mus = mus, Zs = Zs, sigmas=sigmas,
              m=m, B=B, pareto_alpha = pareto_alpha,
              unif_lb=unif_lb, unif_ub=unif_ub)
  res
}
