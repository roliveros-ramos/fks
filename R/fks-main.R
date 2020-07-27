
# Main exported functions -------------------------------------------------

#' Perform a Search on the Number of Knots and Fit Free-Knot Splines To Data
#' Using the Optimal Number of Knots
#' 
#' This function fits free-knot splines to data using every value for the
#' number of knots between \code{minknot} and \code{maxknot}.  The number of
#' knots is then chosen to optimize a fit criterion.  The free-knot spline with
#' the optimum number of knots is returned.
#' 
#' 
#' @param x A vector containing the values of the independent variable.
#' @param y A vector containing the values of the dependent variable.
#' @param degree The degree of the spline fit.
#' @param minknot The minimum number of knots to search.  Defaults to 1.
#' @param maxknot The maximum number of knots to search.  Defaults to 5.
#' @param alg The spline-fitting algorithm.  Choices are "LS" for least-squares
#' and "PS" for P-splines.  Defaults to "LS."
#' @param search The random search algorithm.  Choices are "genetic" for a
#' genetic algorithm and "golden" for a blind random search with golden section
#' adjustment.  Defaults to "genetic."
#' @param knotnumcrit The criterion to be used for determining the number of
#' knots.  Choices are "GCV" for generalized cross-validation, "AIC" for the
#' Akaike information criterion, "AICc" for corrected Akaike information
#' criterion, "BIC" for Bayesian information criterion, "adjAIC" for an
#' adjusted version of the Akaiki information criterion, and "adjGCV" for an
#' adjusted version of generalized cross-validation.  Defaults to "adjGCV."
#' @param k The amount of penalty when AIC is used.  Has no effect with
#' criteria other than AIC.  Defaults to 2.
#' @param d The amount of penalty when adjGCV is used.  Has no effect with
#' criteria other than adjGCV.  Defaults to 3.
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
#' (GCV) for the optimal fit.)} \item{GSJS}{The GSJS estimator, an estimator of
#' the variance of the data.} \item{call}{The function call.}
#' @author Steven Spiriti
#' @seealso \code{\link{fitcriteria}} for the fit criteria,
#' \code{\link{freeknotfit}} for the free-knot spline algorithms.
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
#' xy.freekt = fit.search.numknots(x, y, degree = 2, minknot = 1, maxknot = 3, seed = 555)
#' plot.freekt(xy.freekt, xfit = 0:1000/1000)
#' 
#' @export
fks = function(x, y, k=10, degree=3, prec=NULL, ...) {
  
  if(length(k)==1) k = c(1, k)
  if(length(k)==2) k = sort(k)
  if(length(k)>2) stop("Incorrect k specification.")
  
  dat = data.frame(x=x, y=y)
  
  xknots = fit.search.numknots(x, y, minknot=k[1], maxknot=k[2], degree=degree)
  knots = c(rep(min(x), degree+1), xknots$optknot, rep(max(x), degree+1))
  
  if(!is.null(prec)) knots = round(knots, prec)
  
  m = c(degree, 2)
  kopt = length(knots) - m[1] - 1 
  mod = gam(y ~ s(x, bs="bs", m=m, k=kopt), knots=list(x=knots),
            family = gaussian(link="identity"), dat=dat)
  
  h = 1e-7
  newdat = data.frame(x=unique(knots))
  values = predict(mod, newdata = newdat, type = "response")
  newdat$x = newdat$x + h
  hvalues = predict(mod, newdata = newdat, type = "response") 
  deriv = (hvalues-values)/h
  
  mod$knots = data.frame(knots=unique(knots), value=values, deriv=deriv)
  
  return(mod)
  
}


#' Perform a Search on the Number of Knots and Fit Free-Knot Splines To Data
#' Using the Optimal Number of Knots (freeknotsplines legacy version)
#' 
#' This function fits free-knot splines to data using every value for the
#' number of knots between \code{minknot} and \code{maxknot}.  The number of
#' knots is then chosen to optimize a fit criterion.  The free-knot spline with
#' the optimum number of knots is returned.
#' 
#' 
#' @param x A vector containing the values of the independent variable.
#' @param y A vector containing the values of the dependent variable.
#' @param degree The degree of the spline fit.
#' @param minknot The minimum number of knots to search.  Defaults to 1.
#' @param maxknot The maximum number of knots to search.  Defaults to 5.
#' @param alg The spline-fitting algorithm.  Choices are "LS" for least-squares
#' and "PS" for P-splines.  Defaults to "LS."
#' @param search The random search algorithm.  Choices are "genetic" for a
#' genetic algorithm and "golden" for a blind random search with golden section
#' adjustment.  Defaults to "genetic."
#' @param knotnumcrit The criterion to be used for determining the number of
#' knots.  Choices are "GCV" for generalized cross-validation, "AIC" for the
#' Akaike information criterion, "AICc" for corrected Akaike information
#' criterion, "BIC" for Bayesian information criterion, "adjAIC" for an
#' adjusted version of the Akaiki information criterion, and "adjGCV" for an
#' adjusted version of generalized cross-validation.  Defaults to "adjGCV."
#' @param k The amount of penalty when AIC is used.  Has no effect with
#' criteria other than AIC.  Defaults to 2.
#' @param d The amount of penalty when adjGCV is used.  Has no effect with
#' criteria other than adjGCV.  Defaults to 3.
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
#' (GCV) for the optimal fit.)} \item{GSJS}{The GSJS estimator, an estimator of
#' the variance of the data.} \item{call}{The function call.}
#' @author Steven Spiriti
#' @seealso \code{\link{fitcriteria}} for the fit criteria,
#' \code{\link{freeknotfit}} for the free-knot spline algorithms.
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
#' xy.freekt = fit.search.numknots(x, y, degree = 2, minknot = 1, maxknot = 3, seed = 555)
#' plot.freekt(xy.freekt, xfit = 0:1000/1000)
#' 
#' @export
fit.search.numknots = function(x, y, degree=3, minknot = 1, maxknot = 5, 
                                alg = "LS", search = "genetic",
                                knotnumcrit = "adjGCV", k = 2,
                                d = 3, seed = 5, stream = 0) {
  bestcrit = Inf
  funcname = ""
  answer = NULL
  if ((alg == "LS") && (search == "genetic"))
    funcname = "freelsgen"
  if ((alg == "LS") && (search == "golden"))
    funcname = "freelsgold"
  if ((alg == "PS") && (search == "genetic"))
    funcname = "freepsgen"
  if ((alg == "PS") && (search == "golden"))
    funcname = "freepsgold"
  for (numknot in seq(from = minknot, to = maxknot)) {
      currcall = call(funcname, x, y, degree, numknot, seed, stream)
      currfit = eval(currcall)
      currcrit = switch(knotnumcrit, GCV = currfit$GCV, 
          AIC = AIC(currfit, k = k), adjAIC = adjAIC.freekt(currfit),
          AICc = AICc.freekt(currfit), BIC = BIC(currfit), 
          adjGCV = adjGCV.freekt(currfit, d))
      msg = sprintf("Number of knots = %d, %s = %0.15g", numknot, knotnumcrit, currcrit)
      message(msg)
      if (currcrit < bestcrit)
      {
          bestcrit = currcrit
          answer = currfit
      }     
  }
  return(answer)
}
