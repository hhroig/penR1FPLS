# Mass matrix from fdaPDE (R0):
#
mass_mat_fun <- fdaPDE:::CPP_get.FEM.Mass.Matrix

# Stiffness matrix from fdaPDE (R1):
#
stiff_mat_fun <- fdaPDE:::CPP_get.FEM.Stiff.Matrix


# 1D numerical integration
#
num_int_1d <- function(argvals, f_obs) {

  f <- stats::splinefun(x = argvals, y = f_obs)
  # f <- stats::approxfun(x = argvals, y = f_obs)

  subds <- 100

  res <- "no_convergence"

  count_iter <- 0

  while (!inherits(res, "numeric") && count_iter < 20) {

    res <- try( stats::integrate(f,
                                 min(argvals),
                                 max(argvals),
                                 subdivisions = 10*subds,
                                 stop.on.error = FALSE)$value, silent = TRUE)

    count_iter <- count_iter + 1

  }



  if (!inherits(res, "numeric")) {

    message("Integral did not converge with 10^", count_iter+2, "subdivisions.\n")
    message(res)

  }

  # stopifnot( inherits(res, "numeric") )

  return(res)


}

# 2D Integration, taken from:
# https://stackoverflow.com/questions/64375216/method-for-calculating-volume-under-a-surface-in-r

getVolume <- function(df) {
  #find triangular tesselation of (x,y) grid
  res  <-  geometry::delaunayn(as.matrix(df[ , -3]), options="Qz", output.options=TRUE)
  #calulates sum of truncated prism volumes
  sum(mapply(function(triPoints,A) A/3*sum(df[triPoints,"z"]),
             split.data.frame(res$tri,seq_along(res$areas)),
             res$areas))
}


# Integrated mean square error for the evaluation of the estimated beta (2D domain)
#
imse_beta_2d <- function(beta_hat,
                         beta_true,
                         R0,
                         domain_area){

  imse <- as.numeric((t(beta_hat)%*%R0%*%beta_hat -
                        2*t(beta_hat)%*%R0%*%beta_true +
                        t(beta_true)%*%R0%*%beta_true)/domain_area)

  return(imse)

}


# Integrated mean square error for the evaluation of the estimated beta (1D domain)
#
imse_beta_1d <- function(argvals, beta_true, beta_hat) {

  imse <- num_int_1d(argvals = argvals,
                     f_obs = (beta_true - beta_hat)^2 ) / abs(diff(range(argvals)))

  return(imse)

}
