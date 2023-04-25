#' Predict new data for FPLS-FEM model (R1FPLS algorithm using FE basis)
#'
#' @param object an \code{fpls_fem} model.
#' @param newdata new data matrix.
#' @param ... further arguments.  Currently not used
#'
#' @return predicted values
#' @export
#'
#' @examples
#' #' # Generate data (50 samples, 100 nodes):
#' x <- seq(0, 1, length.out = 10)
#' y <- seq(0, 1, length.out = 10)
#'
#' L <- generate_2d_data(x, y, 50, 3, 0.95)
#'
#' X <- L[["X"]]
#' Y <- L[["Y"]]
#' FEM_basis <- L[["basisobj"]]
#'
#' res_fpls_fem <- r1fpls_fem(X, Y, ncomp = 3, center = TRUE,
#'                              basisobj = FEM_basis, penalty = 10,
#'                              verbose = FALSE, stripped = FALSE )
#'
#' Y_pred <- predict(res_fpls_fem, X)
#' Y_hat <- fitted.values(res_fpls_fem)
predict.r1fpls_fem <-  function(object, newdata, ...){

  Xc <- scale(newdata, center = object$X_mean, scale = FALSE)

  pred <- array(NA, dim = c(nrow(Xc),  dim(object[["coefficients"]])[2:3] ) )

  for (h in 1:(dim(pred)[3])) {

    pred[, , h] <- Xc %*% (object$coefficients[ , , h, drop = F]) + (object$Y_mean)

  }


  return(pred)

}

#' Predict new data for FPLS-FDA model (R1FPLS algorithm using an fda package basis)
#'
#' @param object an \code{fpls_fda} model.
#' @param newdata new data matrix.
#' @param ... further arguments.  Currently not used
#'
#' @return predicted values
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
predict.r1fpls_bs <-  function(object, newdata, ...){

  Xc <- scale(newdata, center = object$X_mean, scale = FALSE)

  pred <- array(NA, dim = c(nrow(Xc),  dim(object[["coefficient_function"]])[2:3] ) )

  for (h in 1:(dim(pred)[3])) {

    # pred[, , h] <- Xc %*% (object$coefficients[ , , h, drop = F]) + (object$Y_mean)

    for (j in 1:ncol(pred)) {
      for (i in 1:nrow(pred)) {

        pred[i , j, h] <- num_int_1d(argvals = object$argvals,
                                     f_obs = Xc[i, ] * object$coefficient_function[ , j, h]) +
          object$Y_mean

      } # i
    } # j
  } # h





  return(pred)

}

#' Predict new data for FPLS-FDA model (penalized FPLS algorithm using an fda basis,
#' as in Aguilera et al. 2016)
#'
#' @param object an \code{fpls_fda} model.
#' @param newdata new data matrix.
#' @param ... further arguments.  Currently not used
#'
#' @return predicted values
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
#'                 ncomp = 3, penalty = 100,
#'                 basisobj = bs_basis, stripped = FALSE)
#'
#' predict(res_ps, newdata = X)
predict.fpls_bs <-  function(object, newdata, ...){

  newXc <- scale(newdata, center = object$X_mean, scale = FALSE)

  # Represent newdata using a B-spline basis:
  bsplineXtest <- fda::Data2fd(argvals = object$argvals,
                               y = Matrix::t(newXc),
                               basisobj = object$basisobj  )

  # Coeff. matrix of size obs. times num. basis:
  A_test <- Matrix::t( bsplineXtest[["coefs"]] )

  # Input for mv-pls:
  new_data_pls <- A_test %*% object$R0 %*% object$inv_t_L

  pred <- stats::predict(object$mvpls_model, new_data_pls) + object$Y_mean

  return(pred)

}



#' Predict new data for FPLS-TPS model (penalized FPLS algorithm using a TPS basis,
#' inspired by Aguilera et al. 2016 methodology)
#'
#' @param object an \code{fpls_tps} model.
#' @param newdata new data matrix.
#' @param ... further arguments.  Currently not used
#'
#' @return predicted values
#' @export
#'
#' @examples
#' #' # Generate data (50 samples, 100 nodes):
#' x <- seq(0, 1, length.out = 10)
#' y <- seq(0, 1, length.out = 10)
#'
#' L <- generate_2d_data(x, y, 50, 3, 0.95)
#'
#' X <- L[["X"]]
#' Y <- L[["Y"]]
#' nodes <- L[["mesh"]][["nodes"]]
#'
#' res_fpls_tps <- fpls_tps(X, Y, ncomp = 3, nodes = nodes, nbasis = 5,
#'                          penalty = 0, verbose = FALSE, stripped = FALSE )
#'
#' Y_pred <- predict(res_fpls_tps, X)
#' Y_hat <- fitted.values(res_fpls_tps)
#' plot(Y_pred, Y_hat)
#' abline(0, 1)
predict.fpls_tps <-  function(object, newdata, ...){

  newdata <- scale(newdata, center = object$X_mean, scale = FALSE)

  # Coeff. matrix of size obs. times num. basis:
  A <- matrix(NA, nrow = nrow(newdata), ncol = object$nbasis)

  for (ind_surface in 1:nrow(newdata)) {

    gam_fit <- mgcv::gam(newdata[ind_surface, ] ~ s(object$nodes[ , 1],
                                                    object$nodes[ , 2],
                                                    bs = "tp",
                                                    k = object$nbasis ))

    A[ind_surface, ] <- gam_fit$coefficients

  }
  colnames(A) <- names(gam_fit$coefficients)

  # Input for mv-pls:
  new_data_pls <- A %*% object$R0 %*% object$inv_t_L

  pred <- stats::predict(object$mvpls_model, new_data_pls) + object$Y_mean

  return(pred)

}

