
#' Plot Fitted Values for Free-Knot Spline
#' 
#' This function plots the fit obtained using a free-knot spline.
#' 
#' 
#' @param x An object of class "\code{freekt}" obtained by using one of the
#' fitting algorithms.
#' @param xfit A vector of x values at which to plot the fitted values.
#' Defaults to the x values of the data.
#' @param linecolor The color of the line.  Defaults to blue.
#' @param lwd The line width.  It is passed to the \code{lines} function.
#' Defaults to 1.
#' @param lty The line type.  It is passed to the \code{lines} function.
#' Defaults to 1.
#' @param \dots Additional arguments to be passed to the \code{plot} function.
#' @return A plot of the data, together with the spline estimator.
#' @author Steven Spiriti
#' @seealso \code{\link{fitted.freekt}} to compute the fitted values.
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
#' @export 
plot.freekt = function(x, xfit = x$x, linecolor="blue", lwd = 1, lty = 1, ...) {
  xfit = as.vector(xfit)
  yfit = fitted.freekt(x, xfit)
  plot(x$x, x$y, ...)
  lines(xfit[order(xfit)], yfit[order(xfit)], col=linecolor, 
        lwd = lwd, lty = lty)
}


#' Summarize Free-Knot Spline Fit
#' 
#' This function displays a summary of the fit obtained using a free-knot
#' spline.
#' 
#' 
#' @param object An object of class "\code{freekt}" obtained by using one of
#' the fitting algorithms.
#' @param \dots Additional arguments to be passed to summary.freekt.  Currently
#' ignored.
#' @return A table containing the values of the optimal amount of penalty (when
#' applicable), the optimal knots, sum of squared errors (SSE), and generalized
#' cross-validation (GCV).
#' @author Steven Spiriti
#' @seealso \code{\link{freeknotfit}} for the fitting algorithms.
#' @keywords
#' @examples
#' 
#' x = 0:30/30
#' truey = x*sin(10*x)
#' set.seed(10556)
#' y = truey + rnorm(31, 0, 0.2)
#' xy.freekt = freelsgen(x, y, degree = 2, numknot = 2, 555)
#' summary.freekt(xy.freekt)
#' 
#' @export
summary.freekt = function(object, ...) {
  answer = NULL
  if (object$lambda != 0)
  {
    currline = c(object$lambda, 
                 rep(NA, times = length(object$optknot)-1))
    answer = rbind(answer, currline)
  }
  currline = object$optknot
  answer = rbind(answer, currline)
  currline = c(object$GCV, 
               rep(NA, times = length(object$optknot)-1))
  answer = rbind(answer, currline)
  RSS = sum((residuals(object))^2)
  currline = c(RSS, rep(NA, times = length(object$optknot)-1))
  answer = rbind(answer, currline)
  if (object$lambda != 0)
    rownames(answer) = 
    c("Optimal lambda", "Optimal knots", "GCV", "RSS") 
  else
    rownames(answer) = c("Optimal knots", "GCV", "RSS") 
  colnames(answer) = rep("", times = length(object$optknot))
  class(answer) = "table"
  print(answer)
}

#' Compute Residuals For Free-Knot Spline
#' 
#' This function computes the residuals, given the optimal values for lambda
#' and the locations of the knots.
#' 
#' 
#' @param object An object of class "\code{freekt}" obtained by using one of
#' the fitting algorithms.
#' @param \dots Additional arguments to be passed to residuals.freekt.
#' Currently ignored.
#' @return A vector containing the residuals.
#' @author Steven Spiriti
#' @seealso \code{\link{fitted.freekt}} to compute the fitted values.
#' @keywords
#' @examples
#' 
#' x = 0:30/30
#' truey = x*sin(10*x)
#' set.seed(10556)
#' y = truey + rnorm(31, 0, 0.2)
#' xy.freekt = freelsgen(x, y, degree = 2, numknot = 2, 555)
#' plot(x, residuals(xy.freekt))
#' 
#' @export
residuals.freekt = function(object, ...) {
  fit = fitted.freekt(object)
  return(object$y - fit)
}

#' Compute Coefficients of B-Splines For Free-Knot Splines
#' 
#' This function computes the coefficients of the B-splines for free-knot
#' splines, given the amount of the penalty (if applicable) and the locations
#' of the knots.
#' 
#' 
#' @param object An object of class "\code{freekt}" obtained by using one of
#' the fitting algorithms.
#' @param \dots Additional arguments to be passed to coef.freekt.  Currently
#' ignored.
#' @return A vector containing the coefficients of the B-splines.
#' @author Steven Spiriti
#' @seealso \code{\link{fitted.freekt}} to compute the fitted values and
#' \code{\link{residuals.freekt}} to compute the residuals.
#' @keywords
#' @examples
#' 
#' x = 0:30/30
#' truey = x*sin(10*x)
#' set.seed(10556)
#' y = truey + rnorm(31, 0, 0.2)
#' xy.freekt = freelsgen(x, y, degree = 2, numknot = 2, 555)
#' coef.freekt(xy.freekt)
#' 
#' @export
coef.freekt = function(object, ...) {
  xdat = object$x
  ydat = object$y
  optknot = object$optknot
  ord = object$degree + 1
  lambda = object$lambda
  fulloptknot = c(rep(min(xdat), ord), optknot, rep(max(xdat), ord))  # includes endpoints
  Xmat = splineDesign(fulloptknot, xdat, ord)
  if ((lambda == 0) | (length(optknot) == 0))
    coef = solve(t(Xmat)%*%Xmat, t(Xmat)%*%ydat) 
  else
  {
    numknots = length(optknot)
    Amat = chgbasismat(fulloptknot, ord)
    Istar = diag(c(rep(0, times = ord), rep(1, times = numknots))) 
    coef = solve(t(Xmat)%*%Xmat + lambda*t(Amat)%*%Istar%*%Amat, 
                 t(Xmat)%*%ydat) 
  }   
  return(coef)
}

#' Compute Fitted Values For Free-Knot Spline
#' 
#' This function computes the fitted values, given the amount of the penalty
#' (if applicable) and the locations of the knots.
#' 
#' 
#' @param object An object of class "\code{freekt}" obtained by using one of
#' the fitting algorithms.
#' @param xfit A vector of x values at which to compute the fitted values.
#' Defaults to the x values of the data.
#' @param \dots Additional arguments to be passed to fitted.freekt.  Currently
#' ignored.
#' @return A vector containing the fitted values.
#' @author Steven Spiriti
#' @seealso \code{\link{residuals.freekt}} for the residuals.
#' @keywords
#' @examples
#' 
#' x = 0:30/30
#' truey = x*sin(10*x)
#' set.seed(10556)
#' y = truey + rnorm(31, 0, 0.2)
#' xy.freekt = freelsgen(x, y, degree = 2, numknot = 2, 555)
#' fitted.freekt(xy.freekt)
#' 
#' @export
fitted.freekt = function(object, xfit = object$x, ...) {        
  xdat = object$x
  ydat = object$y
  optknot = object$optknot
  ord = object$degree + 1
  fulloptknot = c(rep(min(xdat), ord), optknot, rep(max(xdat), ord))  # includes endpoints
  coef = coef.freekt(object)
  yfit = splineDesign(fulloptknot, xfit, ord) %*%coef
  return(yfit)
}


# Fit measures ------------------------------------------------------------


#' Fit Criteria for Free-Knot Splines
#' 
#' These functions compute various criteria for determining the fit of a
#' free-knot spline.  \code{AIC.freekt} computes the Akaike Information
#' Criterion, with \code{k} determining the amount of the penalty.
#' \code{AICc.freekt} computes the corrected Akaike Information Criterion.
#' \code{BIC.freekt} computes the Bayesian Information Criterion, also known as
#' Schwarz Information Criterion.  \code{adjAIC.freekt} computes an adjusted
#' Akaike Information Criterion with the penalty increased to account for the
#' greater flexibility of free knots.  \code{adjGCV.freekt} computes an
#' adjusted GCV with the degrees of freedom increased to account for the
#' greater flexibility of free knots.
#' 
#' 
#' @aliases fitcriteria AIC.freekt AICc.freekt BIC.freekt adjAIC.freekt
#' adjGCV.freekt
#' @param object An object of class "\code{freekt}" obtained by using one of
#' the fitting algorithms.
#' @param k The amount of the penalty.  Used only for \code{AIC.freekt}.
#' @param d The amount of the penalty.  Used only for \code{adjGCV.freekt}.
#' @param \dots Additional arguments to be passed to the \code{AIC.freekt} and
#' \code{BIC.freekt} functions.
#' @return Returns the value of the specified fit criterion.
#' @author Steven Spiriti
#' @seealso \code{\link{fit.search.numknots}}, which uses these fit criteria to
#' determine the number of knots.
#' @references Spiriti, S., Eubank, R., Smith, P., Young, D., "Knot Selection
#' for Least-Squares and Penalized Splines," \emph{Journal of Statistical
#' Computation and Simulation}, in press.
#' @keywords
NULL


#' @export
AIC.freekt = function(object, ..., k = 2) {
  answer = 0
  RSS = sum((residuals(object))^2)
  n = length(object$x)
  npar = object$tracehat
  answer = n * log(RSS/n) + k * npar
  return(answer)
}

#' @export
AICc.freekt = function(object) {
  answer = 0
  RSS = sum((residuals(object))^2)
  n = length(object$x)
  npar = object$tracehat
  answer = n * log(RSS/n) + 2 * npar + 
    2 * npar * (npar + 1) /(n - npar - 1)
  return(answer)
}

#' @export
BIC.freekt = function(object, ...) {
  answer = 0
  RSS = sum((residuals(object))^2)
  n = length(object$x)
  npar = object$tracehat
  answer = n * log(RSS/n) + log(n) * npar
  return(answer)
}

#' @export
adjGCV.freekt = function(object, d = 3) {
  RSS = sum((residuals(object))^2)
  n = length(object$x)
  adjtrace = object$tracehat + d * length(object$optknot)
  answer = (RSS/n) / (1 - adjtrace / n)^2
  return(answer)     
}

#' @export
adjAIC.freekt = function(object) {
  answer = 0
  RSS = sum((residuals(object))^2)
  n = length(object$x)
  npar = object$tracehat
  effdim = 2 * npar - object$degree - 1
  answer = n * log(RSS/n) +  2 * effdim
  return(answer)
}


