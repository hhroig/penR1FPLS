#' Penalized rank-one functional PLS using FE basis functions
#'
#' @param X a number of observations times nodes matrix.
#' @param Y a number of observations times reponses matrix.
#' @param ncomp number of components, integer.
#' @param center logical, indicating if data should be to centered.
#' @param basisobj a Finite Elements basis as in the fdaPDE package.
#' @param penalty a vector of size ncomp indicating the penalty for each component.
#' @param tol convergence tolerance.
#' @param verbose logical, indicating if messages should be printed.
#' @param stripped logical.  If \code{TRUE} the calculations are stripped as
#' much as possible for speed; this is meant for use with cross-validation or
#' simulations when only the coefficients are needed.  Defaults to
#' \code{FALSE}. Inspired by package \code{pls}.
#' @param ... further arguments.  Currently not used
#'
#' @return an r1fpls_fem model.
#' @export
#'
#' @examples
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
#' results_fpls_fem <- r1fpls_fem(X, Y, ncomp = 3, center = TRUE,
#'                              basisobj = FEM_basis, penalty = 10,
#'                              verbose = TRUE )
#'
#' true_coeff <- L[["coefficient_function"]]
#' image(matrix(true_coeff, 10, 10), main = "True coefficient function")
#'
#' image(matrix(results_fpls_fem[["coefficient_function"]][ , , 3], 10, 10),
#'       main = "Estimated coefficient function")
r1fpls_fem <- function(X,
                       Y,
                       ncomp = 3,
                       center = TRUE,
                       basisobj,
                       penalty = 0,
                       tol = .Machine$double.eps^0.5,
                       verbose = TRUE,
                       stripped = FALSE,
                       ...
) {

  tictoc::tic("FPLS-FEM")

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

  # compute mass matrix:
  R0 <- mass_mat_fun(FEMbasis = basisobj)

  # Initialize:
  W <- Matrix::Matrix(data = 0, nrow = n_nodes, ncol = ncomp, sparse = TRUE)
  V <- Matrix::Matrix(data = 0, nrow = n_colY, ncol = ncomp, sparse = TRUE)
  C <- Matrix::Matrix(data = 0, nrow = n_nodes, ncol = ncomp, sparse = TRUE)
  D <- Matrix::Matrix(data = 0, nrow = n_colY, ncol = ncomp, sparse = TRUE)
  TT <- Matrix::Matrix(data = 0, nrow = n_samp, ncol = ncomp, sparse = TRUE)

  W_star_tD <- array(0, dim = c(n_nodes, n_colY, ncomp))
  Beta_hat <- W_star_tD

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
    M <- as.matrix(Matrix::t(Yc_old) %*% Xc_old)

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



      # Solution of w via linear system:
      Psi <- Matrix::.sparseDiagonal(n = n_nodes)
      # Psi <- Matrix::Diagonal(n = n_nodes, x = 1)

      if (penalty[h] == 0) { # no penalties

        # a_solve <- Psi  # same as: t(Psi)%*%Psi
        # b_solve <- z    # same as: t(Psi) %*% z

        w_K <- Matrix::solve( a = Psi, b = z )

      }else {                  # penalized

        # compute stiffness matrix:
        R1 <- stiff_mat_fun(FEMbasis = basisobj)

        a_solve <- cbind(

          rbind(Psi,
                penalty[h]*R1 ),

          rbind(penalty[h]*R1,
                -penalty[h]*R0)
        )

        # a_solve <- Matrix::Matrix(data = NA,
        #                           nrow = 2*n_nodes ,
        #                           ncol = 2*n_nodes ,
        #                           sparse = TRUE )
        #
        # a_solve[1:n_nodes, 1:n_nodes] <- Psi
        # a_solve[(n_nodes+1):(2*n_nodes), 1:n_nodes] <- penalty[h]*R1
        # a_solve[1:n_nodes, (n_nodes+1):(2*n_nodes)] <- penalty[h]*R1
        # a_solve[(n_nodes+1):(2*n_nodes), (n_nodes+1):(2*n_nodes)] <- -penalty[h]*R0

        b_solve <- rbind(z,  # same as: t(Psi) %*% z
                         Matrix::Matrix(nrow = nrow(z),
                                        ncol = ncol(z), data = 0,
                                        sparse = TRUE)
        )

        # Compute w_K
        # but w = Psi %*% w_K = Identity w_K (no need to "evaluate"):
        w <- Matrix::solve( a = a_solve,
                            b = b_solve,
                            sparse = TRUE )[1:n_nodes]

      }

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
    w_norm <- Matrix::t(w) %*% R0 %*% w
    w <- w/sqrt(as.numeric(w_norm))

    # Save w and v:
    W[ , h] <- w
    V[ , h] <- v

    # Scores, also called PLS components:
    tt  <- Xc_old %*% R0 %*% w # as.matrix(Xc_old %*% R0 %*% w)

    TT[ , h] <- tt

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


    # Parameter/coefficient function aux:

    AUX_inv <- Matrix::solve( Matrix::t(C[, 1:h, drop = F]) %*% R0  %*% W[, 1:h, drop = F] )

    W_star_tD[ , , h] <- as.matrix(R0 %*% W[, 1:h, drop = F] %*% AUX_inv %*% Matrix::t(D[, 1:h, drop = F]) )

    if (!stripped) {

      # Parameter/coefficient function:
      Beta_hat[ , , h] <- as.matrix(W[, 1:h, drop = F] %*% AUX_inv %*% Matrix::t(D[, 1:h, drop = F]))

      # Fitted values:
      Y_hat[ , , h] <- as.matrix(TT[, 1:h, drop = F] %*% Matrix::t(D[, 1:h, drop = F])) + Y_mean

    }


  } # end: ncomp loop



  # Rertuns:

  if (stripped) { # fast return (for CV)

    ret <- list(W = W,
                V = V,
                coefficients = W_star_tD,
                X_mean = X_mean,
                Y_mean = Y_mean,
                elapsed = tictoc::toc(quiet = !verbose) )

    class(ret) <- "r1fpls_fem"

  }else {         # full computations


    ret <- list(W = W,
                V = V,
                TT = TT,
                C = C,
                D = D,
                coefficients = W_star_tD,
                fitted.values = Y_hat,
                coefficient_function = Beta_hat,
                Y_mean = Y_mean,
                X_mean = X_mean,
                elapsed = tictoc::toc(quiet = !verbose)
    )

    class(ret) <- "r1fpls_fem"

  }

  return(ret)


} # end: FPLS function
