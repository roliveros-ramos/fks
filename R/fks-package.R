#' Free-Knot Splines
#' 
#' This package is for fitting free-knot splines for data with one independent
#' variable and one dependent variable.
#' 
#' Four fitting methods are included for the case where the number of knots is
#' known in advance; for details, see \code{\link{freeknotfit}}.  In addition,
#' methods are available to compute the fitted values, the residuals, and the
#' coefficients of the splines, and to plot the results, along with a method to
#' summarize the results.  Finally, a function (see
#' \code{\link{fit.search.numknots}}) is provided to optimize the number of
#' knots using a specified fit criterion.  Several fit criteria are provided
#' (see \code{\link{fitcriteria}}).
#' 
#' @name fks-package
#' @docType package
#' @author Steven Spiriti, Philip Smith, and Pierre Lecuyer
#' @references Spiriti, S., Eubank, R., Smith, P., Young, D. 2013. "Knot Selection
#' for Least-Squares and Penalized Splines," \emph{Journal of Statistical
#' Computation and Simulation}, 83(6):1020-1036.
#' @keywords internal
"_PACKAGE"
