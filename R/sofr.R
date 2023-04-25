#' A wrapper for several fitting algorithms
#'
#' @param X a number of observations times nodes matrix.
#' @param Y a number of observations times reponses matrix.
#' @param argvals a set of argument values. \code{NULL} for FEM basis
#' @param nodes a 2-column matrix with the nodes. Needed for "fpls_tps".
#' @param nbasis number of TPS basis to use in method "fpls_tps".
#' @param ncomp number of components, integer.
#' @param center logical, indicating if data should be to centered.
#' @param basisobj a basis object from package \code{fda}.
#' @param method one of: "r1fpls_fem", "fpls_tps", "r1fpls_bs", or "fpls_bs".
#' @param penalty a vector of size ncomp indicating the penalty for each component.
#' @param tol convergence tolerance.
#' @param verbose logical, indicating if messages should be printed.
#' @param stripped logical.  If \code{TRUE} the calculations are stripped as
#' much as possible for speed; this is meant for use with cross-validation or
#' simulations when only the coefficients are needed.  Defaults to
#' \code{FALSE}. Inspired by package \code{pls}.
#' @param ... currently not used.
#'
#' @return  a sofr model
#' @export
#'
#' @examples
#' # 1D example:
#' # library(pls)
#' # library(fda)
#'
#' # Octane number:
#' Y <- as.matrix( pls::gasoline$octane )
#'
#' # Gasoline NIR spectra:
#' X <- as.matrix( as.data.frame( pls::gasoline$NIR ) )
#'
#' # Wavenumber:
#' argvals <- seq(900, 1700, by = 2)
#'
#' # Ruppert's law:  nbasis = nbreaks + norder - 2  and norder = degree + 1
#' n_breaks <- min(round(length(argvals)/4), 40)
#' n_basis <- n_breaks + (3+1) - 2
#'
#' # B-spline basis:
#' bs_basis <- fda::create.bspline.basis(rangeval = range(argvals),
#'                                  nbasis = n_basis)
#'
#'
#'
#' res_ps <-  sofr(X = X, Y = Y, argvals = argvals,
#'                 ncomp = 3, center = TRUE, penalty = 100,
#'                 basisobj = bs_basis, method = "r1fpls_bs", stripped = FALSE)
#'
#' predict(res_ps, newdata = X)
#'
#' # 2D example:
#' # Generate data (50 samples, 100 nodes):
#' x <- seq(0, 1, length.out = 10)
#' y <- seq(0, 1, length.out = 10)
#'
#' L <- generate_2d_data(x, y, 50, 3, 0.95)
#'
#' X <- L[["X"]]
#' Y <- L[["Y"]]
#' FEM_basis <- L[["basisobj"]]
#'
#' results_fpls_fem <- sofr(X, Y, ncomp = 3, center = TRUE,
#' basisobj = FEM_basis, method = "r1fpls_fem",  penalty = 0,
#' verbose = TRUE )
sofr <- function(X,
                 Y,
                 argvals = NULL,
                 nodes = NULL,
                 nbasis = NULL,
                 ncomp = 3,
                 center = TRUE,
                 basisobj,
                 method = NULL,
                 penalty = 0,
                 tol = .Machine$double.eps^0.5,
                 verbose = TRUE,
                 stripped = FALSE,
                 ...
) {


  if (!is.matrix(Y)) {
    Y <- tryCatch(
      {
        # Just to highlight: if you want to use more than one
        # R expression in the "try" part then you'll have to
        # use curly brackets.
        # 'tryCatch()' will return the last evaluated expression
        # in case the "try" part was completed successfully

        as.matrix(Y)
      },
      error=function(cond) {
        message("Converting Y to a matrix does not work. Check the error:")
        message(cond)
        # Choose a return value in case of error
        return(Y)
      },
      warning=function(cond) {
        message("Converting Y to a matrix caused a warning:")
        message(cond)
        # Choose a return value in case of warning
        return(Y)
      } )
  }


  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("Y", 1:dim(Y)[2])
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", 1:dim(X)[2])
  }

  # Select fit function:

  method <- match.arg(method, c("r1fpls_fem", "r1fpls_bs", "fpls_bs", "fpls_tps"))

  fitFunc <- switch(method,
                    r1fpls_fem = r1fpls_fem,
                    r1fpls_bs = r1fpls_bs,
                    fpls_bs = fpls_bs,
                    fpls_tps = fpls_tps
  )

  res <- fitFunc(X = X, Y = Y, argvals = argvals, ncomp = ncomp, center = center,
                 basisobj = basisobj, penalty = penalty, nodes = nodes,
                 nbasis = nbasis, verbose = verbose, ... )


  return(res)

}
