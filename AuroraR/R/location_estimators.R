#' Simple location estimators
#'
#' @template z_matrix
#' @template method_output
#' @name simpleloc
NULL

#' @describeIn simpleloc row-wise sample mean (average)
#' @export
location_mean <- function(Zs) {rowMeans(Zs)}

#' @describeIn simpleloc row-wise median
#' @export
location_median <- function(Zs) {apply(Zs, 1, stats::median)}

#' @describeIn simpleloc row-wise midrange
#' @export
location_midrange <- function(Zs) {apply(Zs, 1, function(x) mean(range(x)))}


pareto_mle_one_unit <- function(xs){
  x <- min(xs)
  alpha <- 1/mean(log(xs/x))
  alpha*x/(alpha-1)
}

#' @describeIn simpleloc Maximum likelihood estimator of mean,
#'  assuming each row consists of i.i.d. samples from the Pareto distribution
#'  with unknown parameters.
#' @export
pareto_mle <- function(Zs) {apply(Zs, 1, pareto_mle_one_unit)}





#' James-Stein estimator
#'
#' @template z_matrix
#' @param sigma Numeric, noise standard deviation (of replicate average, i.e.,
#'     each noisy measurements has a standard deviation equal to sqrt(K)*sigma).
#' @param positive_part Bool, whether to use the positive part James-Stein estimator
#'    which uniformly dominates the James-Stein estimator (defaults to TRUE).
#' @template method_output
#' @export
james_stein <- function(Zs, sigma=1, positive_part=TRUE){
  zs <- rowMeans(Zs) # reduction to sufficient statistics
  mu_bar <- mean(zs)
  zs_norm_squared <- sum( (zs - mu_bar)^2 )
  m <- length(zs)
  shrink_factor <- 1 - (m-3)*sigma^2/zs_norm_squared
  if (positive_part){
    shrink_factor <- max(0, shrink_factor)
  }
  z_pred <- (zs - mu_bar)*shrink_factor + mu_bar
  z_pred
}

#' GMLEB assuming Gaussian average measurements (CLT)
#'
#' @template z_matrix
#' @param sigma Numeric, noise standard deviation (of replicate average, i.e.,
#'     each noisy measurements has a standard deviation equal to sqrt(K)*sigma).
#' @template method_output
#' @references This is just a wrapper around the implementation of GLmix in the REBayes package
#'      (cite REBayes if you use this function).
#' @export
gmleb <- function(Zs, sigma=1 ){
  zs <- rowMeans(Zs) # reduction to sufficient statistics
  npmle_fit <- REBayes::GLmix(zs, sigma=sigma)
  z_pred <- stats::predict(npmle_fit, zs)
  z_pred
}

# SURE and GroupLinear estimators


# bayes rule for fixed lambda,mu
normal_posterior_mean <- function(X,sigma_squared,lambda,mu){
  lambda/(lambda+sigma_squared) * X + sigma_squared/(lambda+sigma_squared) * mu
}

sure_grandmean_objective <- function(lambda, X, sigma_squared){
  ss <- sigma_squared
  n <- length(X)
  sum(  ( ss/(ss+lambda) )^2 * (X - mean(X))^2 + ss/(ss+lambda) * (lambda - ss + 2/n * ss)  )
}



#' SURE shrinkage in heteroskedastic Gaussian problems towards grand mean
#'
#' This code is directly adapted from https://github.com/MaZhuang/grouplinear,
#' with minor changes.
#'
#' @template z_matrix
#' @param sigma Numeric, noise standard deviation (of replicate average, i.e.,
#'     each noisy measurements has a standard deviation equal to sqrt(K)*sigma).
#' @template method_output
#' @references Xie, X., Kou, S. C., & Brown, L. D. (2012). SURE estimates for a heteroscedastic hierarchical model.
#'          Journal of the American Statistical Association, 107(500), 1465-1479.
#' @export
sure_grandmean <- function(Zs, sigma=1){
  zs <- rowMeans(Zs) # reduction to sufficient statistics

  n <- length(zs)
  if (length(sigma) == 1){
    sigma <-  rep(sigma,n)
  }
  sigma_squared <- sigma^2
  lambda <- stats::optimize(sure_grandmean_objective,lower=0,upper=1000,
                     X=zs,sigma_squared=sigma_squared)$minimum
  normal_posterior_mean(zs,sigma_squared,lambda,mean(zs))
}






## spherically symmetric estimator with c_n = c^*_n
spher_shrink <- function(x.,v.){
  n. <- length(x.)
  if ( (n.==1) | (stats::var(x.)==0) ) x. else {
    cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0)
    bhat <- min( cstar*mean(v.)/stats::var(x.), 1 )
    x. - bhat*(x. - mean(x.))
  }
}

#' Group-linear estimator
#'
#' This code is directly adapted from https://github.com/MaZhuang/grouplinear,
#' with only minor changes.
#'
#' @template z_matrix
#' @param sigma Numeric, noise standard deviation (of replicate average, i.e.,
#'     each noisy measurements has a standard deviation equal to sqrt(K)*sigma).
#' @param nbreak Integer, number of breaks for group-linear estimators. Defaults to floor(N^(1/3)).
#' @template method_output
#' @references "Group-linear empirical Bayes estimates for a heteroscedastic normal mean"
#'        by A Weinstein, Z Ma, LD Brown, CH Zhang (JASA, 2018)
#' @export
grouplinear <- function( Zs, sigma=1, nbreak=floor(nrow(Zs)^(1/3)) ){  # default: bin log(v) into same NUMBER (=n^(1/3) of intervals
  zs <- rowMeans(Zs)
  n <- length(zs)

  if (length(sigma) == 1){
    sigma <-  rep(sigma,n)
  }

  sigma_squared <- sigma^2

  splitby=cut(log(sigma_squared), nbreak, labels=FALSE)
  xsub <- split(zs,splitby)
  vsub <- split(sigma_squared,splitby)
  indexsub <- split(1:n,splitby)
  thetahatsub <- mapply(spher_shrink, xsub, vsub)
  indexsub_unlist <- as.vector( unlist(indexsub) )
  thetahatsub_unlist <- as.vector( unlist(thetahatsub) )
  thetahat <- thetahatsub_unlist[order(indexsub_unlist)]
  return(thetahat)
}


