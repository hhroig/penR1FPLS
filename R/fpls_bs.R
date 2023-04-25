#' Penalized functional PLS based on the (1D) basis representation of the data
#' Method by Aguilera et al. 2016. Default are B-spline basis.
#'
#' @param X a number of observations times nodes matrix.
#' @param Y a number of observations times reponses matrix.
#' @param argvals a set of argument values.
#' @param ncomp number of components, integer.
#' @param center logical, indicating if data should be to centered.
#' @param basisobj a basis object from package \code{fda}.
#' @param penalty penalty value, numeric.
#' @param tol convergence tolerance.
#' @param verbose logical, indicating if messages should be printed.
#' @param stripped logical.  If \code{TRUE} the calculations are stripped as
#' much as possible for speed; this is meant for use with cross-validation or
#' simulations when only the coefficients are needed.  Defaults to
#' \code{FALSE}. Inspired by package \code{pls}.
#' @param ... further arguments.  Currently not used
#'
#' @return an fpls_bs model.
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
#' res_ps <-  fpls_bs(X = X, Y = Y, argvals = argvals,
#'                 ncomp = 3, center = TRUE, penalty = 100,
#'                 basisobj = bs_basis, stripped = FALSE)
#'
#' predict(res_ps, newdata = X)
fpls_bs <- function(X,
                    Y,
                    argvals,
                    ncomp = 3,
                    center = TRUE,
                    basisobj,
                    penalty = 0,
                    tol = .Machine$double.eps^0.5,
                    verbose = TRUE,
                    stripped = FALSE,
                    ...
) {

  tictoc::tic("FPLS-FDA")

  if (center) {

    # center:
    Xc <- scale(X, scale = F)
    Yc <- scale(Y, scale = F)

    # get mean of Y and X to add in the end:
    Y_mean <- attr(Yc, "scaled:center")
    X_mean <- attr(Xc, "scaled:center")


  }else {

    Xc <- X
    Yc <- Y
    X_mean <- rep(0, ncol(X))
    Y_mean <- 0

  } # center data

  # number of nodes and samples:
  n_samp <- nrow(X)
  n_nodes <- ncol(X)
  n_colY <- ncol(Y)

  # Check if there are as many penalties as components ncomp:
  if ( length(penalty) > 1 ) {
    cat("penalty should be a single numeric value.\n")
    penalty <- penalty[1]
  }


  # Evaluate basis functions:
  # (rows corresponding to argument values and columns to basis functions)
  Psi <- fda::eval.basis(evalarg = argvals,
                         basisobj = basisobj,
                         Lfdobj=0, returnMatrix=TRUE)
  tPsi <- Matrix::t(Psi)

  # Matrix of inner products (mass):
  R0 <- fda::bsplinepen(basisobj = basisobj, Lfdobj = 0)

  # Penalty matrix:
  if (basisobj$type == "bspline") { # Psplines penalty matrix

    dorder <- 2
    delta <- diff(diag(basisobj$nbasis),
                  differences = dorder) #d-order differences
    P <- Matrix::t(delta) %*% delta

  }

  # Matrix L with penalty:
  LLprim <- R0 + penalty*P
  L <- expm::sqrtm(LLprim)
  inv_t_L <- Matrix::solve(L)

  # Represent X using a B-spline basis:
  bsplineX <- fda::Data2fd(argvals = argvals,
                           y = Matrix::t(Xc),
                           basisobj = basisobj  )

  # Coeff. matrix of size obs. times num. basis:
  A <- Matrix::t(bsplineX[["coefs"]])

  data_pls <- A %*% R0 %*% inv_t_L

  # PLS model:
  mvpls_model <- pls::plsr(Yc ~ data_pls,
                           ncomp =  ncomp,
                           method = "oscorespls",
                           center = FALSE,
                           scale = FALSE)


  # Rertuns:
  if (stripped) { # fast return (for CV)

    ret <- list(argvals = argvals,
                basisobj = basisobj,
                R0 = R0,
                inv_t_L = inv_t_L,
                X_mean = X_mean,
                Y_mean = Y_mean,
                mvpls_model = mvpls_model,
                ncomp = ncomp,
                elapsed = tictoc::toc(quiet = !verbose) )

    class(ret) <- "fpls_bs"

  }else {         # full computations

    # Get MV-model components:
    TT <- as.matrix(mvpls_model[["scores"]])

    V <-  as.matrix(mvpls_model[["Yscores"]]) # scores Y

    C <- as.matrix(mvpls_model[["loadings"]] )# regress. loadings A_Phi_L_train

    D <-  as.matrix(mvpls_model[["Yloadings"]]) # regress. loadings A_Phi_L_train

    W <- as.matrix( mvpls_model[["loading.weights"]] )

    # B <- as.matrix(mvpls_model[["coefficients"]] )


    # Coefficient function:

    Beta_hat <- array(NA, dim = c(n_nodes, n_colY, ncomp))

    for (h in 1:ncomp) {

      # # Coeffs. of the basis representation of Beta(p):
      W_star <- W[, 1:h, drop = F] %*% Matrix::solve( Matrix::t(C[, 1:h, drop = F]) %*% W[, 1:h, drop = F] )
      beta_coeff_h <- inv_t_L %*% W_star %*% Matrix::t(D[, 1:h, drop = F])

      # Beta_hat observed in p_1, ..., p_m
      Beta_hat[ , , h] <- Psi %*% beta_coeff_h
    }



    ret <- list(argvals = argvals,
                basisobj = basisobj,
                R0 = R0,
                inv_t_L = inv_t_L,
                mvpls_model = mvpls_model,
                ncomp = ncomp,
                V = V,
                TT = TT,
                C = C,
                D = D,
                W = W,
                X_mean = X_mean,
                Y_mean = Y_mean,
                fitted.values = mvpls_model$fitted.values + Y_mean,
                coefficient_function = Beta_hat,
                elapsed = tictoc::toc(quiet = !verbose)
    )

    class(ret) <- "fpls_bs"

  }

  return(ret)


} # end: FPLS function
