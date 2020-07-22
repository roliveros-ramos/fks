#' Fit Free-Knot Splines To Data
#' 
#' These functions fit free-knot splines to data with one independent variable
#' and one dependent variable.  It is assumed that the number of knots is known
#' in advance.  \code{freelsgen} and \code{freelsgold} fit least-squares
#' splines with no penalty, while \code{freepsgen} and \code{freepsgold} fit
#' penalized splines.  \code{freelsgen} and \code{freepsgen} use a genetic
#' algorithm, while \code{freelsgold} and \code{freepsgold} use a blind search
#' augmented with a golden section algorithm.
#' 
#' 
#' @name freeknotfit
#' @aliases freelsgen freelsgold freepsgen freepsgold
#' @param x A vector containing the values of the independent variable.
#' @param y A vector containing the values of the dependent variable.
#' @param degree The degree of the spline fit.
#' @param numknot The number of knots.
#' @param seed The value of the initial seed.  Defaults to 5.
#' @param stream The value of the initial stream to be used for parallel
#' programming.  Defaults to 0.
#' @return An object of class "\code{freekt}" containing the following
#' components:
#' 
#' \item{x}{A vector containing the x values.} \item{y}{A vector containing the
#' y values.} \item{degree}{The degree of the spline fit.} \item{seed}{The
#' value of the initial seed.} \item{stream}{The value of the stream.}
#' \item{lambda}{The optimum amount of penalty.  This is automatically equal to
#' 0 for \code{freelsgen} and \code{freelsgold}.} \item{optknot}{A vector
#' containing the optimal knots.} \item{tracehat}{The trace of the hat matrix
#' for the optimal fit.} \item{GCV}{The value of generalized cross validation
#' (GCV) for the optimal fit.} \item{GSJS}{The GSJS estimator, an estimator of
#' the variance of the data.} \item{call}{The function call.}
#' @author Steven Spiriti, Philip Smith, and Pierre Lecuyer
#' @seealso \code{\link{fit.search.numknots}} for the case where the number of
#' knots is not specified in advance.
#' @references Eubank, R. (1999), \emph{Nonparametric Regression and Spline
#' Smoothing}, New York: Marcel Dekker, Inc., Second ed.\cr \cr Spiriti, S.,
#' Eubank, R., Smith, P., Young, D., "Knot Selection for Least-Squares and
#' Penalized Splines," \emph{Journal of Statistical Computation and
#' Simulation}, in press.
#' @keywords
#' @examples
#' 
#' x = 0:30/30
#' truey = x*sin(10*x)
#' set.seed(10556)
#' y = truey + rnorm(31, 0, 0.2)
#' xy.freekt = freelsgen(x, y, degree = 2, numknot = 2, 555)
#' plot.freekt(xy.freekt, xfit = 0:1000/1000)
#'
NULL

freelsgen = function(x, y, degree, numknot, seed = 5, stream = 0) {
  n = length(x)
  ord = degree + 1
  optknot = rep(0, times = numknot)
  tracehat = 0
  GCV = 0
  GSJS = 0
  result = .C("freelsgen", as.integer(n), as.double(x), as.double(y), 
               as.integer(ord), as.integer(numknot), as.integer(seed), 
               as.integer(stream), as.double(optknot), 
               as.double(tracehat), as.double(GCV), as.double(GSJS), PACKAGE = "fks")
  answer = list(x = x, y = y, degree = as.integer(degree),
                 seed = as.integer(seed), stream = as.integer(stream),
                 lambda = 0, optknot = result[[8]], 
                 tracehat = result[[9]], GCV = result[[10]], 
                 GSJS = result[[11]], call = match.call())
  class(answer) = "freekt"
  return(answer)
}

freelsgold = function(x, y, degree, numknot, seed = 5, stream = 0) {
  n = length(x)
  ord = degree + 1
  optknot = rep(0, times = numknot)
  tracehat = 0
  GCV = 0
  GSJS = 0
  result = .C("freelsgold", as.integer(n), as.double(x), as.double(y), 
               as.integer(ord), as.integer(numknot), as.integer(seed), 
               as.integer(stream), as.double(optknot), 
               as.double(tracehat), as.double(GCV), as.double(GSJS), PACKAGE = "fks")
  answer = list(x = x, y = y, degree = as.integer(degree),
                 seed = as.integer(seed), stream = as.integer(stream),
                 lambda = 0, optknot = result[[8]], 
                 tracehat = result[[9]], GCV = result[[10]], 
                 GSJS = result[[11]], call = match.call())
  class(answer) = "freekt"
  return(answer)
}

freepsgen = function(x, y, degree, numknot, seed = 5, stream = 0) {
  n = length(x)
  ord = degree + 1
  optknot = rep(0, times = numknot)
  lambda = 0
  tracehat = 0
  GCV = 0
  GSJS = 0
  result = .C("freepsgen", as.integer(n), as.double(x), as.double(y), 
               as.integer(ord), as.integer(numknot), as.integer(seed), 
               as.integer(stream), as.double(lambda), 
               as.double(optknot), as.double(tracehat), 
               as.double(GCV), as.double(GSJS), PACKAGE = "fks")
  answer = list(x = x, y = y, degree = as.integer(degree),
                 seed = as.integer(seed), stream = as.integer(stream),
                 lambda = result[[8]], optknot = result[[9]], 
                 tracehat = result[[10]], GCV = result[[11]], 
                 GSJS = result[[12]], call = match.call())
  class(answer) = "freekt"
  return(answer)
}

freepsgold = function(x, y, degree, numknot, seed = 5, stream = 0) {
  n = length(x)
  ord = degree + 1
  optknot = rep(0, times = numknot)
  lambda = 0
  tracehat = 0
  GCV = 0
  GSJS = 0
  result = .C("freepsgold", as.integer(n), as.double(x), as.double(y), 
               as.integer(ord), as.integer(numknot), as.integer(seed), 
               as.integer(stream), as.double(lambda), 
               as.double(optknot), as.double(tracehat), 
               as.double(GCV), as.double(GSJS), PACKAGE = "fks")
  answer = list(x = x, y = y,  degree = as.integer(degree),
                 seed = as.integer(seed), stream = as.integer(stream),
                 lambda = result[[8]], optknot = result[[9]], 
                 tracehat = result[[10]], GCV = result[[11]], 
                 GSJS = result[[12]], call = match.call())
  class(answer) = "freekt"
  return(answer)
}
