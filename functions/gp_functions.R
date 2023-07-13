# Functions for GP predictions

calc_point_distances = function(coords, coords2){
  if(missing(coords2)){
    distances <- as.matrix(dist(coords, upper=TRUE, diag=TRUE))
  } else {
    distances <- raster::pointDistance(coords, coords2, lonlat = FALSE)
  }

  return(distances)
}

predict_gp <- function(iter, alpha, rho, obs_gp_function, obs_distances, pred_distances,
                       obs_to_pred_distances, beta, Xpred, OS, rho_alpha_median = TRUE,
                       phi = NULL, dist = "negbin", return_gp = FALSE,
                       kernel = c("quad","matern32","matern52","exp")){
  # Estimated covariance kernel for the fitted points
  kern_quad<- function(x, alpha, rho) {
    return(alpha^2 * exp(-x^2/ (2*rho^2)))
  }
  kern_matern32<- function(x, alpha, rho) {
    return(alpha^2 *(1 + (sqrt(3)*x)/rho) * exp(-(sqrt(3)*x)/rho))
  }
  kern_matern52<- function(x, alpha, rho) {
    return(alpha^2 *(1 + (sqrt(5)*x)/(3*rho^2)) * exp(-(sqrt(5)*x)/rho))
  }
  kern_exp<- function(x, alpha, rho) {
    return(alpha^2 * exp(-x/rho))
  }

  kernel <- match.arg(kernel, c("quad","matern32","matern52","exp"))
  if(identical(kernel, "quad")) kern<- kern_quad
  if(identical(kernel, "matern32")) kern<- kern_matern32
  if(identical(kernel, "matern52")) kern<- kern_matern52
  if(identical(kernel, "exp")) kern<- kern_exp

  npred<- dim(pred_distances)[1]
  obs_gp_function<- as.vector(as.matrix(obs_gp_function[iter,]))
  if(rho_alpha_median) {
    alpha<- median(alpha)
    rho<- median(rho)
  } else {
  alpha<- alpha[iter]
  rho<- rho[iter]
  }
  beta<- beta[iter,]
  if(dist == "negbin") {
  phi <- phi[iter]
  }
  K_obs <- kern(obs_distances, alpha, rho) + diag(1e-9, dim(obs_distances)[1])


  # covariance between prediction points and fitted points
  K_new <- kern(obs_to_pred_distances, alpha, rho)

  # Estimated covariance kernel for prediction points
  K_star <- kern(pred_distances, alpha, rho) + diag(1e-9, npred)

  gc(verbose = F)

  # Estimated mean for prediction points
  gp_pred<- t(K_new) %*% solve(K_obs, obs_gp_function) +
              MASS::mvrnorm(1, mu = rep(0, npred), Sigma = K_star - t(K_new) %*% solve(K_obs, K_new))

  if(return_gp) {
    return(gp_pred)
  }

  lambda <- exp(Xpred %*% beta + log(OS) + gp_pred)
  gc(verbose = F)
  if(dist == "poisson") {
  return(rpois(npred, lambda))
    } else if(dist == "negbin") {
      return(rnbinom(n = npred, mu = lambda, size = phi))
    }

}

