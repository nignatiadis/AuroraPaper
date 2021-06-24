#' @return A list with the following entries:
#'      \item{m}{The number of units in the simulation.}
#'      \item{B}{The number of replicates per unit in the simulation.}
#'      \item{Zs}{A m*B matrix of noisy measurements. Each row corresponds to a different units.
#'     The K measurements in the i-th row are the replicate measurements for the i-th unit.}
#'      \item{sigmas}{The noise standard deviation of the measurements of each unit. Can be either a single number (if
#'         the simulation is homoskedastic) or a vector of length m.}
#'      \item{...}{Other entries that contains simulation specific parameters.}

