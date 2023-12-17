#' 3D data generation
#'
#' @param nodes 3D nodes (expanded grid).
#' @param tetrahedrons tetrahedrons as in fdaPDE.
#' @param beta_num an integer 1:3 of possible betas.
#' @param Rsq coefficient of determinaiton between a clean Y and a noisy Y.
#' @param num_samples number of samples.
#'
#' @return a list with the simulated data.
#' @export
#'
#' @examples
#'
#' library(fdaPDE)
#' data(sphere3Ddata)
#'
#' nodes=sphere3Ddata$nodes
#' tetrahedrons=sphere3Ddata$tetrahedrons
#'
#' new3d_data <- generate_3d_data(nodes = nodes, tetrahedrons = tetrahedrons,
#' num_samples = 30, beta_num = 3, Rsq = 0.9)
#'
#' # library(plotly)
#' #
#' # gen_beta <- as.numeric(new3d_data[["coefficient_function"]] )
#' #
#' # fig_coeff <- plot_ly(x = ~nodes[, 1], y = ~nodes[, 2], z = ~nodes[, 3],
#' #                      marker = list(color = ~ gen_beta,
#' #                      colorscale = c('#FFE1A1', '#683531'),
#' #                      showscale = TRUE))
#' # fig_coeff <- fig_coeff %>% add_markers( alpha = 0.2)
#' # fig_coeff
#' #
#' #
#' # X1 <- as.numeric(new3d_data[["X"]][1, ] ) # first observation
#' #
#' # fig_X1 <- plot_ly(x = ~nodes[, 1], y = ~nodes[, 2], z = ~nodes[, 3],
#' #                      marker = list(color = ~ X1, colorscale = c('#FFE1A1',
#' #                       '#683531'), showscale = TRUE))
#' # fig_X1 <- fig_X1 %>% add_markers( alpha = 0.2)
#' # fig_X1
generate_3d_data <- function(nodes,
                             tetrahedrons,
                             num_samples = 100,
                             beta_num = 3,
                             Rsq = 0.95) {

  # create 3D mesh:
  mesh <- fdaPDE::create.mesh.3D(nodes = nodes, tetrahedrons = tetrahedrons)

  # FEM basis:
  FEM_basis  <-  fdaPDE::create.FEM.basis(mesh)

  # Mass matrix:
  R0 <- mass_mat_fun(FEMbasis = FEM_basis)

  # Generate X:
  X = NULL

  for(ii in 1:num_samples){
    a1 = stats::rnorm(1, mean = 1, sd = 0.2)
    a2 = stats::rnorm(1, mean = 1, sd = 0.2)
    a3 = stats::rnorm(1, mean = 1, sd = 0.2)

    func_evaluation = numeric(nrow(mesh$nodes))

    for (i in 1:nrow(mesh$nodes)){

      func_evaluation[i] = a1* cos(2*pi*mesh$nodes[i,1]) +
        a2* cos(2*pi*mesh$nodes[i,2]) +
        a3* sin(2*pi*mesh$nodes[i,3])^2
      1

    }
    data = func_evaluation + stats::rnorm(nrow(mesh$nodes), mean = 0, sd = 0.2)
    X = rbind(X, data)
  }


  # Generate beta(x, y):
  if (beta_num == 1) {

    # centered at (0, 0, 0):
    r  <-  0.4 # set the r parameter
    beta_z  <-  5*exp(-((nodes[, 1] )^2 + (nodes[, 2] )^2 + (nodes[, 3] )^2)/( 2*r^2 ))

  }else if (beta_num == 2) {

    # top right corner
    beta_z <- 5*exp(-((nodes[, 1] - 0.25)^2 + (nodes[, 2] - 0.25)^2 + (nodes[, 3] - 0.25)^2 ) /( 2*0.4^2 ))

  }else if (beta_num == 3) {

    # bottom left + top right corner:
    beta_z <- 5*exp(-((nodes[, 1] - 0.25)^2 + (nodes[, 2] - 0.25)^2 + (nodes[, 3] - 0.25)^2 )/( 2*0.4^2 )) +
      5*exp(-((nodes[, 1] + 0.3)^2 + (nodes[, 2] + 0.3)^2 + (nodes[, 3] + 0.3)^2 )/( 2*0.2^2 ))

  }


  beta_true <- as.matrix(beta_z)

  # Center X:
  Xc <- scale(X, scale = F)

  # No-noise Y:
  Y_clean <- as.matrix(Xc %*% R0 %*% beta_true)

  # Variance of errors:
  var_e <- (1/Rsq - 1)*stats::var(Y_clean)

  # Noisy Y:
  Y <- Y_clean +
    as.matrix(stats::rnorm(length(Y_clean), mean = 0, sd = sqrt(var_e)))

  return(list(X = X,
              Y = Y,
              Y_clean = Y_clean,
              basisobj = FEM_basis,
              mesh = mesh,
              coefficient_function = beta_true ))

}
