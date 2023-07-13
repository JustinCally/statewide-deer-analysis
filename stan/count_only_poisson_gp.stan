data {
  int<lower=0> N;                      // number of observations
  real delta;                          // bin width
  int<lower=1> n_site;                 // site
  int<lower=1> n_distance_bins;        // distance bins
  int<lower=1> n_gs;                   // number of group sizes
  vector[n_gs] gs;                     // group sizes
  vector[n_distance_bins] midpts;      // midpt of each interval
  real<lower=1> max_distance;         // truncation distance
  int<lower=1> max_int_dist;          // max distance as an integer
  real<lower=0> theta_frac;           // fraction of camera view
  array[n_site] int effort;           // effort
  array[n_site, n_gs] int n_obs;      //number of observations
  array[n_site, n_distance_bins, n_gs] int y; // observations matrix

  // availability prior information
  real<lower=0> bshape;               // availability shape
  real<lower=0> bscale;               // availability scale

  // camera trap detection parameters
  int<lower=0> det_ncb;                 // number of covariates for detection model
  matrix[n_site, det_ncb] det_model_matrix; // detection model matrix
  array[n_gs, n_distance_bins] real pa ; // averaged availability for multiple individuals

  // Abundance/occupancy model matrix
  int<lower = 1> m_psi;
  matrix[n_site, m_psi] X_psi;
  // GP param
  array[n_site] vector[2] coords;
}

transformed data {
  // vector[n_distance_bins] pi; // availability for point
  vector[n_site] survey_area;
    for(i in 1:n_site) {
      survey_area[i] = theta_frac * effort[i] * pi() * (max_distance/1000)^2;
    }

//  for (i in 1:n_distance_bins) {
//    //pi[i] = delta/max_distance; // availability function (line transect)
//    pi[i] = (2 * delta * midpts[i])/max_distance^2; // point transect
//  }

}

parameters {
  // Site RE
  real<lower=0> site_sd;
  vector[n_site] site_raw;
 // abundance parameters
  simplex[n_gs] eps_ngs; // random variation in group size
  vector[m_psi] beta_psi;
  // detection parameters
  vector[det_ncb] beta_det;
  //real<lower=0> theta;
  // temporal availability parameters
  real<lower=0, upper=1> activ;
  // GP params
  real<lower=0> rho;    // GP length scale
  real<lower=0> alpha;  // GP marginal sd
  vector[n_site] eta;
}

transformed parameters {
  // re
  vector[n_site] eps_site;
  // distance parameters
  array[n_site] real log_sigma;
  array[n_site] real sigma;
  array[n_site, n_distance_bins] real p_raw; // detection probability
  array[n_site, n_distance_bins, n_gs] real p_raw_scale; // detection probability for point independence model and accounting for availability
  array[n_site, n_distance_bins, n_gs] real<upper=0> log_p_raw; // log detection probability accounting for availability
  array[n_site, n_gs] real log_p; // log probability with summed exponential for the multinomial logit model
  array[n_site, n_gs] real<lower=0,upper=1> p; // log probability with summed exponential for the multinomial logit model
  // abundance parameters
  array[n_site, n_gs] real<lower=0> lambda;
  // activity parameters
  real log_activ = log(activ);
  vector[n_site] log_lambda_psi;
  // GP params
  vector[n_site] gp_predict = cholesky_decompose(gp_exp_quad_cov(coords, alpha, rho) +
                                diag_matrix(rep_vector(1e-9, n_site))) * eta;

for(n in 1:n_site) {
  log_sigma[n] = det_model_matrix[n,] * beta_det;
  sigma[n] = exp(log_sigma[n]);
  for (i in 1:n_distance_bins) {
      // assuming a half-normal detection fn from line
       p_raw[n,i] = exp(-(midpts[i])^2/(2*sigma[n]^2)); // half-normal function (pg 419 of AHM)
      // hazard function
      // p_raw[n,i] = 1 - exp(-(midpts[i]/theta)^(-sigma[n])); // hazard function (pg 507 of AHM)
    for(j in 1:n_gs) {
      p_raw_scale[n,i,j] = p_raw[n,i]*pa[j,i]; //  pr(animal occurs and is detected in bin i)
      log_p_raw[n,i,j] = log(p_raw_scale[n,i,j]);
      }
  }

  eps_site[n] = site_sd * site_raw[n]; // site random effect sd centreing
// define log lambda
  log_lambda_psi[n] = X_psi[n,] * beta_psi + gp_predict[n] + eps_site[n];
// convert to occupancy psi
  // logit_psi[n] = inv_cloglog(log_lambda_psi[n]);
  // log_psi[n] = log_inv_logit(logit_psi[n]);
  // log1m_psi[n] = log1m_inv_logit(logit_psi[n]);

  for(j in 1:n_gs) {
    log_p[n,j] = log_sum_exp(log_p_raw[n,,j]);
    p[n,j] = exp(log_p[n,j]);
    // model site abundance
    lambda[n,j] = exp(log_lambda_psi[n] + log_p[n,j] + log_activ + log(eps_ngs[j])) .* survey_area[n];
  }
  }
    }

model {
  beta_det ~ normal(0, 4); // prior for sigma
  eps_ngs ~ uniform(0, 1); // prior for group size effect
  beta_psi ~ normal(0, 2); // prior for poisson model
  activ ~ beta(bshape, bscale);  //informative prior
  // GP priors
  eta ~ std_normal();
  rho ~ inv_gamma(16.152, 5.38782);
  alpha ~ normal(0, 1);
  site_sd ~ normal(0, 1);
  site_raw ~ std_normal();

  for(n in 1:n_site) {
  for(j in 1:n_gs) {
  n_obs[n,j] ~ poisson(lambda[n,j]);
  y[n,,j] ~ multinomial_logit(to_vector(log_p_raw[n,,j]));
  }
}
}

generated quantities {
  array[n_site,n_gs] real n_obs_pred;
  array[n_site, n_gs] real n_obs_true;
  array[n_site] real N_site;
  array[n_site] real N_site_pred;
  array[n_site, max_int_dist+1] real DetCurve;
  array[n_site, n_gs] real log_lik1;
  array[n_site, n_gs] real log_lik2;
  array[n_site] real log_lik;
  vector[n_site] Site_lambda;
  vector[n_site] psi;
  // array[npc] real pred;
  // real Nhat;


for(n in 1:n_site) {
  for(j in 1:n_gs) {
  log_lik1[n,j] = multinomial_logit_lpmf(y[n,,j] | to_vector(log_p_raw[n,,j])); //for loo
  log_lik2[n,j] =  poisson_lpmf(n_obs[n,j] | lambda[n,j]); //for loo
  n_obs_true[n, j] = gs[j] * (poisson_rng(exp(log_lambda_psi[n] + log(eps_ngs[j]))));
  n_obs_pred[n,j] = gs[j] * (poisson_rng(exp(log_lambda_psi[n] + log_p[n,j] + log_activ + log(eps_ngs[j])) .* survey_area[n]));
    }
    log_lik[n] = log_sum_exp(log_sum_exp(log_lik1[n,]), log_sum_exp(log_lik2[n,]));
    Site_lambda[n] = exp(log_lambda_psi[n]);
    N_site[n] = sum(n_obs_true[n,]);
    N_site_pred[n] = sum(n_obs_pred[n,]);
    // Occupancy probability transformation
    psi[n] = inv_cloglog(log_lambda_psi[n]);

  // loop over distance bins
  for(j in 0:max_int_dist) { // get DS predictions for distance 0 to max bin distance
    DetCurve[n, j+1] = exp(-(j+0.5)^2/(2*sigma[n]^2)); // half normal
    // DetCurve[n, j+1] =  1 - exp(-(j+0.5/theta)^(-sigma[n])); //hazard rate
    }
}

}
