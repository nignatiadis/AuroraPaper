% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulations.R
\name{homoskedasticloc}
\alias{homoskedasticloc}
\alias{normal_normal_sim}
\alias{normal_laplace_sim}
\alias{normal_uniform_sim}
\title{Homoskedastic location family simulations with Normal prior}
\usage{
normal_normal_sim(m, B, sigma, A, prior_mu = 0)

normal_laplace_sim(m, B, sigma, A, prior_mu = 0)

normal_uniform_sim(m, B, sigma, A, prior_mu = 0)
}
\arguments{
\item{m}{Integer, number of units.}

\item{B}{Integer, number of replicates.}

\item{sigma}{Numeric, noise standard deviation (of replicate average, i.e.,
each noisy measurements has a standard deviation equal to sqrt(B)*sigma).}

\item{A}{Numeric, prior variance}

\item{prior_mu}{Numeric, mean of the prior (defaults to 0).}
}
\value{
A list with the following entries:
     \item{m}{The number of units in the simulation.}
     \item{B}{The number of replicates per unit in the simulation.}
     \item{Zs}{A m*B matrix of noisy measurements. Each row corresponds to a different units.
    The K measurements in the i-th row are the replicate measurements for the i-th unit.}
     \item{sigmas}{The noise standard deviation of the measurements of each unit. Can be either a single number (if
        the simulation is homoskedastic) or a vector of length m.}
     \item{...}{Other entries that contains simulation specific parameters.}
}
\description{
Homoskedastic location family simulations with Normal prior
}
\section{Functions}{
\itemize{
\item \code{normal_normal_sim}: Normal prior, Normal likelihood

\item \code{normal_laplace_sim}: Normal prior, Laplace likelihood

\item \code{normal_uniform_sim}: Normal prior, Rectangular (Uniform) likelihood
}}

