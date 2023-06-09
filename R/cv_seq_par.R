#' Sequential crossvalidation in parallel (different penalty for each component).
#'
#' @param X a number-of-observations times nodes matrix.
#' @param Y a number-of-observations times reponses matrix.
#' @param center logical, indicating if data should be centered.
#' @param argvals a set of argument values. Not needed for FEM basis.
#' @param nodes a 2-column matrix with the nodes. Needed for "fpls_tps".
#' @param nbasis number of TPS basis to use in method "fpls_tps".
#' @param penalty_vec a vector of possible penalties.
#' @param ncomp number of components, integer.
#' @param folds a user defined list of folds (as generated by caret::createFolds())
#' or an integer indicating the number of folds.
#' @param basisobj a Finite Elements basis as in the fdaPDE package.
#' @param R0 (mass) matrix of inner products between TPS basis functions.
#' @param P penalty matrix (optional for method "fpls_tps").
#' @param method only supported by: "r1fpls_fem".
#' @param tol convergence tolerance.
#' @param verbose logical, indicating if messages should be printed.
#' @param stripped logical.  If \code{TRUE} the calculations are stripped as
#' much as possible for speed. Particularly, if \code{FALSE} (default) it computes
#' the final models using the best combination of penalties.
#' Inspired by package \code{pls}.
#'
#' @return A list of crossvalidates erros (CVEs) and penalties giving the minimum
#' CVEs per number of components.
#' @export
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#'
#' @examples
#' # 2D example:
#'
#' # Generate data (30 samples, 100 nodes):
#' x <- seq(0, 1, length.out = 10)
#' y <- seq(0, 1, length.out = 10)
#'
#' L <- generate_2d_data(x, y, 30, 3, 0.95)
#'
#' X <- L[["X"]]
#' Y <- L[["Y"]]
#' FEM_basis <- L[["basisobj"]]
#'
#' cv_fem <- cv_seq_par(X = X, Y = Y, penalty_vec = c(0.01, 10, 100),
#'                   ncomp = 3, folds = 5, basisobj = FEM_basis,
#'                   method = "r1fpls_fem",
#'                   verbose = TRUE, stripped = FALSE)
cv_seq_par <- function(X,
                       Y,
                       center = TRUE,
                       argvals = NULL,
                       nodes = NULL,
                       nbasis = NULL,
                       penalty_vec,
                       ncomp = min(10, ncol(X)),
                       folds = 5,
                       basisobj = NULL,
                       R0 = NULL,
                       P = NULL,
                       method = NULL,
                       tol = .Machine$double.eps^0.5,
                       verbose = TRUE,
                       stripped = TRUE ) {

  tictoc::tic("Crossvalidation")

  # Initialize grid for the first component
  penalty_grid <- matrix(penalty_vec)
  colnames(penalty_grid) <- "ncomp_1"

  if (is.numeric(folds)) {

    num_folds <- folds
    folds <- caret::createFolds(Y, k = num_folds)

  }else if (is.list(folds)) {

    num_folds <- length(folds)

  }

  # Initialize CVEs:
  CVEs_ncomp <- array(data = NA, dim = ncomp) # averaged
  names(CVEs_ncomp) <- paste0("ncomp_", 1:ncomp)

  MSE_ncomp_fold <- matrix(data = NA,
                           nrow = ncomp,
                           ncol = num_folds) # MSE per component per fold
  colnames(MSE_ncomp_fold) <- paste0("fold_", 1:num_folds)
  rownames(MSE_ncomp_fold) <- paste0("ncomp_", 1:ncomp)

  for (ncomp_i in 1:ncomp) {

    if (verbose) {
      cat("Component ", ncomp_i, "/", ncomp, "\n")
    }

    i <- row_lambda <- NULL

    MSE_lambda_fold <-
      foreach::foreach(i = 1:num_folds,
                       .packages = c("penR1FPLS"),
                       .combine = "cbind") %:%
      foreach::foreach(row_lambda = 1:nrow(penalty_grid),
                       .packages = c("penR1FPLS"),
                       .combine = 'c' ) %dopar%
      {

        # build train
        Y_fold_train <- Y[-folds[[i]], , drop = F]
        X_fold_train <- X[-folds[[i]], , drop = F]

        # build test:
        Y_fold_test <- Y[folds[[i]], , drop = F]
        X_fold_test <- X[folds[[i]], , drop = F]

        res_fpls <- sofr(X = X_fold_train,
                         Y = Y_fold_train,
                         argvals = argvals,
                         nodes = nodes,
                         nbasis = nbasis,
                         ncomp = ncomp_i,
                         center = center,
                         basisobj = basisobj,
                         R0 = R0,
                         P = P,
                         method = method,
                         penalty = as.numeric(penalty_grid[row_lambda, ]),
                         tol = tol,
                         verbose = FALSE,
                         stripped = stripped )

        mean(as.numeric(Y_fold_test - stats::predict(object = res_fpls,
                                                     newdata = X_fold_test)[, , ncomp_i] )^2)


      } # nested loop


    # Averaged MSE_fold:
    CVEs_ncomp_lambda <- rowMeans(MSE_lambda_fold)

    # Best penalties per component:
    sel_lambda <- which.min(CVEs_ncomp_lambda)

    best_penalties <- penalty_grid[sel_lambda, , drop = F]

    if (ncomp_i < ncomp) {

      grid_names <- c(colnames(best_penalties), paste0("ncomp_", ncomp_i+1))

      penalty_grid <- data.frame(best_penalties,
                                 penalty_vec, row.names = NULL)

      colnames(penalty_grid) <- grid_names
    }


    # Save the folds-averaged CV error:
    CVEs_ncomp[ncomp_i] <- CVEs_ncomp_lambda[sel_lambda]

    # Save MSEs per fold, for the best lambda:
    MSE_ncomp_fold[ncomp_i, ] <- MSE_lambda_fold[sel_lambda, ]

  } # loop in ncomp: number of components


  names(CVEs_ncomp) <- paste0("ncomp_", 1:ncomp)
  names(best_penalties) <- paste0("ncomp_", 1:ncomp)
  colnames(MSE_ncomp_fold) <- paste0("fold_", 1:num_folds)
  rownames(MSE_ncomp_fold) <- paste0("ncomp_", 1:ncomp)


  if (stripped) {
    ret <- list(
      CVEs_ncomp = CVEs_ncomp,
      MSE_ncomp_fold = MSE_ncomp_fold,
      best_penalties = best_penalties,
      elapsed = tictoc::toc(quiet = !verbose)
    )

  }else {

    if (verbose) {
      cat("Fitting final model\n")
    }

    final_model <- sofr(X = X,
                        Y = Y,
                        argvals = argvals,
                        nodes = nodes,
                        nbasis = nbasis,
                        ncomp = ncomp,
                        center = center,
                        basisobj = basisobj,
                        R0 = R0,
                        P = P,
                        method = method,
                        penalty = as.numeric(best_penalties),
                        tol = tol,
                        verbose = verbose,
                        stripped = stripped )

    ret <- list(
      CVEs_ncomp = CVEs_ncomp,
      MSE_ncomp_fold = MSE_ncomp_fold,
      best_penalties = best_penalties,
      final_model = final_model,
      elapsed = tictoc::toc(quiet = !verbose)
    )

  }



  return(ret)


}
