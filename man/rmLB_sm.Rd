% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rmLB_sm}
\alias{rmLB_sm}
\title{Sample from the matrix Langevin-Bingham on the Stiefel manifold.}
\usage{
rmLB_sm(num, mJ, mH, mC, mX, ir)
}
\arguments{
\item{num}{number of observations or sample size.}

\item{mJ}{symmetric ip*ip matrix}

\item{mH}{symmetric ir*ir matrix}

\item{mC}{ip*ir matrix}

\item{mX}{ip*ir matrix, the initial value}

\item{ir}{ir}
}
\value{
an array containing a sample of draws from the matrix Langevin-Bingham on the Stiefel manifold.
}
\description{
This function draws a sample from the matrix Langevin-Bingham on the Stiefel manifold.
}
\details{
The matrix Langevin-Bingham distribution on the Stiefel manifold has the density kernel:
\deqn{f(X) \propto \mathrm{etr}\{ H X' J X + C' X \}}
where \eqn{X} satisfies \eqn{X'X = I_r}, and \eqn{H} and \eqn{J} are symmetric matrices.
}
\author{
Yukai Yang, \email{yukai.yang@statistik.uu.se}

#' @section References:
Hoff, P. D. (2009) "Simulation of the Matrix Bingham—von Mises—Fisher Distribution, With Applications to Multivariate and Relational Data", Journal of Computational and Graphical Statistics, Vol. 18, pp. 438-456.
}
\keyword{utils}
