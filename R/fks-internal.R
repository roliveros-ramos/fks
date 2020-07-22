
#' Compute Change-of-Basis Matrix to Convert Truncated Power Basis to B-Spline
#' Basis
#' 
#' This function computes a change-of-basis matrix that converts from the
#' truncated power basis to the B-spline basis.  It is used by
#' \code{coef.freekt} and should not be called directly by the user.
#' 
#' 
#' @param knot The vector of knots.  The first and last knots are repeated
#' \code{ord} times to make the length equal to the dimension of the spline
#' space.
#' @param ord The order of the spline, which is one more than the degree.
#' @return The change-of-basis matrix that converts from the truncated power
#' basis to the B-spline basis.
#' @author Steven Spiriti
#' @references Smith, P. (1982), "Hypothesis Testing B-Spline Regression,"
#' \emph{Communications in Statistics, Series B}, 11, 143-157.\cr \cr Spiriti,
#' S., Eubank, R., Smith, P., Young, D., "Knot Selection for Least-Squares and
#' Penalized Splines," \emph{Journal of Statistical Computation and
#' Simulation}, in press.
chgbasismat = function(knot, ord) {
  dimmat = length(knot) - ord
  answer = matrix(0, nrow = dimmat, ncol = dimmat)
  for (j in 0:(ord-1))
  {
    brow = splineDesign(knot, knot[1], ord, j)
    brow = as.vector(brow/factorial(j))
    answer[j + 1, ] = brow
  }
  nknot = sort(-1*knot)
  for (j in 1:(dimmat - ord))
  {
    brow = splineDesign(knot, knot[ord + j], ord, ord - 1)
    brow2 = splineDesign(nknot, nknot[length(knot) - ord - (j - 1)],
                         ord, ord - 1)
    brow2 = brow2[dimmat:1]
    brow = brow + (-1)^ord * brow2
    brow = as.vector(brow/factorial(ord - 1))
    answer[ord + j, ] = brow
  }
  return(answer)
}


