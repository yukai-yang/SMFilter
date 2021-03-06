% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{FDist2}
\alias{FDist2}
\title{Compute the squared Frobenius distance between two matrices.}
\usage{
FDist2(mX, mY)
}
\arguments{
\item{mX}{a \eqn{p \times r} matrix where \eqn{p \geq r}.}

\item{mY}{another \eqn{p \times r} matrix where \eqn{p \geq r}.}
}
\value{
the Frobenius distance.
}
\description{
This function Compute the squared Frobenius distance between two matrices.
}
\details{
The Frobenius distance between two matrices is defined to be
\deqn{d(X, Y) = \sqrt{ \mathrm{tr} \{ A' A \} }}
where \eqn{A = X - Y}.

The Frobenius distance is a possible measure of the distance between two points on the Stiefel manifold.
}
\examples{

FDist2(runif_sm(1,4,2)[1,,], runif_sm(1,4,2)[1,,])

}
\author{
Yukai Yang, \email{yukai.yang@statistik.uu.se}
}
\keyword{utils}
