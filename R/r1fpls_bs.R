#' Penalized rank-one functional PLS using fda-package basis (B-spline basis as default)
#'
#' @param X a number of observations times nodes matrix.
#' @param Y a number of observations times reponses matrix.
#' @param argvals a set of argument values.
#' @param ncomp number of components, integer.
#' @param center logical, indicating if data should be to centered.
#' @param basisobj a basis object from package \code{fda}.
#' @param penalty a vector of size ncomp indicating the penalty for each component.
#' @param tol convergence tolerance.
#' @param verbose logical, indicating if messages should be printed.
#' @param stripped logical.  If \code{TRUE} the calculations are stripped as
#' much as possible for speed; this is meant for use with cross-validation or
#' simulations when only the coefficients are needed.  Defaults to
#' \code{FALSE}. Inspired by package \code{pls}.
#' @param ... further arguments.  Currently not used
#'
#' @return an r1fpls_bs model.
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
#' res_ps <-  r1fpls_bs(X = X, Y = Y, argvals = argvals,
#'                 ncomp = 3, center = TRUE, penalty = 100,
#'                 basisobj = bs_basis, stripped = FALSE)
#'
#' predict(res_ps, newdata = X)
r1fpls_bs <- function(X,
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

  # number of nodes and samples:
  n_samp <- nrow(X)
  n_nodes <- ncol(X)
  n_colY <- ncol(Y)

  # Check if there are as many penalties as components ncomp:
  if ( length(penalty) != ncomp ) {
    penalty <- rep(penalty, times = ncomp)
  }


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

  # Evaluate basis functions:
  # (rows corresponding to argument values and columns to basis functions)
  Psi <- fda::eval.basis(evalarg = argvals,
                         basisobj = basisobj,
                         Lfdobj=0, returnMatrix=TRUE)
  tPsi <- Matrix::t(Psi)

  # Penalty matrix:
  if (basisobj$type == "bspline") { # Psplines penalty matrix

    dorder <- 2
    delta <- diff(diag(basisobj$nbasis),
                  differences = dorder) #d-order differences
    P <- Matrix::t(delta) %*% delta

  }


  # Initialize:
  W <- Matrix::Matrix(data = 0, nrow = n_nodes, ncol = ncomp, sparse = TRUE)
  V <- Matrix::Matrix(data = 0, nrow = n_colY, ncol = ncomp, sparse = TRUE)
  C <- Matrix::Matrix(data = 0, nrow = n_nodes, ncol = ncomp, sparse = TRUE)
  D <- Matrix::Matrix(data = 0, nrow = n_colY, ncol = ncomp, sparse = TRUE)
  TT <- Matrix::Matrix(data = 0, nrow = n_samp, ncol = ncomp, sparse = TRUE)


  Beta_hat <- array(0, dim = c(n_nodes, n_colY, ncomp))

  Y_hat <- array(0, dim = c(n_samp, n_colY, ncomp))

  # Loop FPLS components:

  Xc_old <- Xc # as(Xc, "CsparseMatrix")
  Yc_old <- Yc # as(Yc, "CsparseMatrix")

  for (h in 1:ncomp) {

    if (verbose) {
      cat("\nComputing component", h, "\n")
      cat("  Using penalty =", penalty[h], "\n")
    }

    # Compute M:
    M = as.matrix(Matrix::t(Yc_old) %*% Xc_old)

    # SVD:
    Msvd <- svd(M)

    # Initialize weights:
    v <- as.matrix( Msvd[["u"]][, 1] )
    w <- as.matrix( Msvd[["v"]][, 1] )

    # Loop on v, w:
    w_change <- tol + 100
    v_change <- tol + 100

    count_conv <- 0

    while (w_change > tol & v_change > tol) {

      # Update v:

      if (ncol(Yc_old) > 1) {

        v <-  M %*% w
        v  <-  as.matrix( v/norm(v, type = "2") )

        # Update z:
        z <- Matrix::t(M) %*% v

      }else {

        v <- Matrix::Matrix(data = 1, nrow = 1, ncol = 1, sparse = TRUE)

        # Update z:
        z <-  Matrix::t(M)

      }

      if (penalty[h] == 0 || basisobj$type != "bspline") { # no penalties

        w_K <- Matrix::solve( tPsi %*% Psi) %*% tPsi %*% z

      }else {                  # penalized

        w_K <- Matrix::solve( tPsi %*% Psi + (penalty[h]*P) ) %*% tPsi %*% z

      }

      # Evaluate:
      w <- Psi %*% w_K

      # save w and v old:
      v_old <- v
      w_old <- w

      #control convergence:
      w_change <- sum(abs(w_old - w))
      v_change <- sum(abs(v_old - v))

      # count convergence steps:
      count_conv <- count_conv + 1

    } # end: v, w loop

    if (verbose) {
      cat("     Convergence in", count_conv, "step(s)\n\n")
    }

    # Normalize w according to functional norm:
    w_norm <- num_int_1d(argvals = argvals, f_obs = w^2)
    w <- w/ sqrt( w_norm )

    # Save w and v:
    W[ , h] <- w
    V[ , h] <- v # not needed

    # Scores, also called PLS components:
    for (tt_row in 1:n_samp) {

      TT[tt_row , h] <- num_int_1d(argvals = argvals,
                                   f_obs = Xc_old[tt_row, ] * w)

    }

    tt  <- TT[, h]

    # Deflate coefficients:
    c <- (Matrix::t(Xc_old) %*% tt) / (  sum(tt^2) )
    d <- (Matrix::t(Yc_old) %*% tt) / ( sum(tt^2) )

    # Save coefficients:
    D[ , h] <- d
    C[ , h] <- c

    # Deflate:
    Xcnew <- Xc_old - tt %*% Matrix::t(c)
    Ycnew <- Yc_old - tt %*% Matrix::t(d)

    # Update roles:
    Xc_old <- Xcnew
    Yc_old <- Ycnew



    # Parameter/coefficient function:
    Phi <- W[, 1:h, drop = FALSE]
    if (ncol(Phi) > 1) {

      for (k_Phi in 2:ncol(Phi)) {

        sub <- 0

        for (i_sub in 1:(k_Phi - 1)) {

          sub <- sub +
            num_int_1d(argvals = argvals,
                       f_obs = (C[, i_sub]*W[, k_Phi])) * Phi[, i_sub]

        } # loop: i_sub

        Phi[, k_Phi] <- Phi[, k_Phi] - sub

      } # loop: k_Phi

    } # id: ncol(Phi)

    Beta_hat[ , , h] <- as.matrix(Phi %*% Matrix::t(D[, 1:h, drop = FALSE]))

    if (!stripped) {
      # Fitted values:
      Y_hat[ , , h] <- as.matrix(TT[, 1:h, drop = F] %*% Matrix::t(D[, 1:h, drop = F])) + Y_mean
    }


  } # end: ncomp loop




  # Rertuns:
  if (stripped) { # fast return (for CV)

    ret <- list(argvals = argvals,
                W = W,
                V = V,
                coefficient_function = Beta_hat,
                X_mean = X_mean,
                Y_mean = Y_mean,
                elapsed = tictoc::toc(quiet = !verbose) )

    class(ret) <- "r1fpls_bs"

  }else {         # full computations


    ret <- list(argvals = argvals,
                W = W,
                V = V,
                TT = TT,
                C = C,
                D = D,
                fitted.values = Y_hat,
                coefficient_function = Beta_hat,
                Y_mean = Y_mean,
                X_mean = X_mean,
                elapsed = tictoc::toc(quiet = !verbose)
    )

    class(ret) <- "r1fpls_bs"

  }

  return(ret)


} # end: FPLS function
