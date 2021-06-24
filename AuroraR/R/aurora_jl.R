.aurora <- new.env(parent = emptyenv())

#' Setup AuroraR
#'
#' This function initializes Julia and the Aurora.jl package.
#' The first time will be long since it includes precompilation.
#' Additionally, this will install Julia and the required packages
#' if they are missing.
#'
#' @param installJulia logical, whether to install Julia automatically when Julia is not found,
#'           whose default value is TRUE.
#' @param pkg_check logical, check for Aurora.jl package and install if necessary.
#' @param ... arguments passed to \code{JuliaCall::julia_setup}.
#'
#' @references Similar functions are available from other R packages that call JuliaCall,
#' for example: \cr
#' https://github.com/SciML/diffeqr/blob/master/R/diffeqr.R \cr
#' https://github.com/Non-Contradiction/ipoptjlr/blob/master/R/IPOPT.R
#' @export
aurora_setup <- function(installJulia=TRUE, pkg_check=TRUE, ...) {
  .aurora$julia <- JuliaCall::julia_setup(installJulia=installJulia,...)
  if(pkg_check) JuliaCall::julia_install_package_if_needed("Aurora")
  .aurora$julia$library("Aurora")
  JuliaCall::autowrap("Auroral")
  JuliaCall::autowrap("CoeyCunningham")
  JuliaCall::autowrap("AuroraKNN")
  JuliaCall::autowrap("Aurora.FittedAuroral")
  JuliaCall::autowrap("Aurora.FittedAuroraKNN")
  .aurora$julia$source(system.file("julia/aurora_wrappers.jl", package = "AuroraR"))
  NULL
}

#' Auroral estimators
#'
#' @template z_matrix
#' @param predictions_only Boolean, whether to return only vector of predictions, or
#'    to also return the objects containing more informations about the fits. Defaults to TRUE.
#' @template method_output
#' @name auroralmethods
NULL

#' @describeIn auroralmethods Auroral: Aurora with linear regression
#' @export
auroral <- function(Zs, predictions_only=TRUE){
  auroral_fit <- .aurora$julia$call("auroral", Zs)
  if (predictions_only) {
    return(auroral_fit$mus)
  } else {
    return(auroral_fit)
  }
}

#' @describeIn auroralmethods Coey-Cunningham: Similar to Aurora with linear regression,
#'  but K-1 replicates are first summarized to their sample mean
#' @export
ebcc <- function(Zs, predictions_only=TRUE){
  ebcc_fit <- .aurora$julia$call("ebcc", Zs)
  if (predictions_only) {
    return(ebcc_fit$mus)
  } else {
    return(ebcc_fit)
  }
}


#' Aurora with K-Nearest Neighbors
#'
#' @template z_matrix
#' @param loocv Boolean, whether to choose K by leave-one-out-cross-validation. Defaults to TRUE.
#' @param kKNN Integer. If loocv = FALSE, this specifies the number of neighbors to use.
#'        If loocv = TRUE, then K is chosen by LOOCV from the set {1, ..., kKNN}. Defaults to 1000.
#' @param predictions_only Boolean, whether to return only vector of predictions, or
#'    to also return the objects containing more informations about the fits. Defaults to TRUE.
#' @template method_output
#' @export
aurora_knn <- function(Zs, loocv=TRUE, kKNN = 1000,  predictions_only=TRUE){
  aurora_knn_fit <- .aurora$julia$call("aurora_knn", Zs, kKNN, loocv)
  if (predictions_only) {
    return(aurora_knn_fit$mus)
  } else {
    return(aurora_knn_fit)
  }
}

