#' Penalized functional PLS based on the (3D) TPS representation of the data.
#' Inspired by Aguilera et al. 2016.
#'
#' @param X a number of observations times nodes matrix.
#' @param Y a number of observations times reponses matrix.
#' @param nodes a set of x, y, z nodes.
#' @param ncomp number of components, integer.
#' @param center logical, indicating if data should be to centered.
#' @param nbasis number of TPS basis to use.
#' @param penalty penalty value, numeric.
#' @param tol convergence tolerance.
#' @param verbose logical, indicating if messages should be printed.
#' @param stripped logical.  If \code{TRUE} the calculations are stripped as
#' much as possible for speed; this is meant for use with cross-validation or
#' simulations when only the coefficients are needed.  Defaults to
#' \code{FALSE}. Inspired by package \code{pls}.
#' @param R0 (mass) matrix of inner products between basis functions.
#' @param P penalty matrix.
#' @param ... further arguments.  Currently not used
#'
#' @return an fpls-tps model.
#' @export
#'
#' @examples
#' #
fpls_tps3 <- function(X,
                     Y,
                     nodes,
                     ncomp = 3,
                     center = TRUE,
                     nbasis,
                     penalty = 0,
                     tol = .Machine$double.eps^0.5,
                     verbose = TRUE,
                     stripped = FALSE,
                     R0 = NULL,
                     P = NULL,
                     ...
) {

  tictoc::tic("FPLS-TPS")

  if (center) {

    # center:
    Xc <- scale(X, scale = FALSE)
    Yc <- scale(Y, scale = FALSE)

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

  # Coeff. matrix of size obs. times num. basis:
  A <- matrix(NA, nrow = nrow(X), ncol = nbasis)

  for (ind_surface in 1:nrow(X)) {

    gam_fit <- mgcv::gam(Xc[ind_surface, ] ~ s(nodes[ , 1],
                                               nodes[ , 2],
                                               nodes[ , 3],
                                               bs = "tp",
                                               k = nbasis))

    A[ind_surface, ] <- gam_fit$coefficients

  }
  colnames(A) <- names(gam_fit$coefficients)

  # Evaluate basis functions:
  # (rows corresponding to argument values and columns to basis functions)
  # X is approximated by A %*% t(Psi):
  Psi <- stats::model.matrix(gam_fit)
  tPsi <- Matrix::t(Psi)

  # Matrix of inner products (mass):
  if (is.null(R0)) {
    stop("Please provide the R0 matrix")
  }


  # Penalty matrix:
  if ( is.null(P) ) { # Psplines penalty matrix

    dorder <- 2
    delta <- diff(diag(nbasis),
                  differences = dorder) #d-order differences
    P <- Matrix::t(delta) %*% delta

  }

  # Matrix L with penalty:
  LLprim <- R0 + penalty*P
  L <- expm::sqrtm(LLprim)
  inv_t_L <- Matrix::solve(L)

  data_pls <- A %*% R0 %*% inv_t_L

  # PLS model:
  mvpls_model <- pls::plsr(Yc ~ data_pls,
                           ncomp =  ncomp,
                           method = "oscorespls",
                           center = FALSE,
                           scale = FALSE )


  # Rertuns:
  if (stripped) { # fast return (for CV)

    ret <- list(nodes = nodes,
                nbasis = nbasis,
                R0 = R0,
                X_mean = X_mean,
                Y_mean = Y_mean,
                inv_t_L = inv_t_L,
                mvpls_model = mvpls_model,
                ncomp = ncomp,
                elapsed = tictoc::toc(quiet = !verbose) )

    class(ret) <- "fpls_tps"

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


    ret <- list(nodes = nodes,
                nbasis = nbasis,
                R0 = R0,
                inv_t_L = inv_t_L,
                mvpls_model = mvpls_model,
                gam_model = gam_fit,
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

    class(ret) <- "fpls_tps"

  }

  return(ret)


} # end: FPLS function
