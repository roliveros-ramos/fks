% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fks-methods.R
\name{residuals.freekt}
\alias{residuals.freekt}
\title{Compute Residuals For Free-Knot Spline}
\usage{
\method{residuals}{freekt}(object, ...)
}
\arguments{
\item{object}{An object of class "\code{freekt}" obtained by using one of
the fitting algorithms.}

\item{\dots}{Additional arguments to be passed to residuals.freekt.
Currently ignored.}
}
\value{
A vector containing the residuals.
}
\description{
This function computes the residuals, given the optimal values for lambda
and the locations of the knots.
}
\examples{

x = 0:30/30
truey = x*sin(10*x)
set.seed(10556)
y = truey + rnorm(31, 0, 0.2)
xy.freekt = freelsgen(x, y, degree = 2, numknot = 2, 555)
plot(x, residuals(xy.freekt))

}
\seealso{
\code{\link{fitted.freekt}} to compute the fitted values.
}
\author{
Steven Spiriti
}
