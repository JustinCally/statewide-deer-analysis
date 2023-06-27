functions {

  /* Efficient computation of the horseshoe prior
   * Args:
   *   zb: standardized population-level coefficients
   *   global: global horseshoe parameters
   *   local: local horseshoe parameters
   *   scale_global: global scale of the horseshoe prior
   *   c2: positive real number for regularization
   * Returns:
   *   population-level coefficients following the horseshoe prior
   */
  vector horseshoe(vector zb, vector[] local, real[] global,
                   real scale_global, real c2) {
    int K = rows(zb);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    vector[K] lambda2 = square(lambda);
    real tau = global[1] * sqrt(global[2]) * scale_global;
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return zb .* lambda_tilde * tau;
  }
}

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

  // summary of whether species is known to be present at each site
  int<lower = 0, upper = 1> any_seen[n_site];

  // number of surveys at each site
  int<lower = 0> n_survey[n_site];

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
  // Prediction for abundance model
  int<lower = 1> npc;
  matrix[npc, m_psi] X_pred_psi;
  vector[npc] prop_pred;
  // negbinom scale
  real reciprocal_phi_scale;
  // hs params
  real<lower=0> hs_df; // horseshoe prior param
  real<lower=0> hs_df_global;
  real<lower=0> hs_df_slab;
  real<lower=0> hs_scale_global;
  real<lower=0> hs_scale_slab;
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
  vector[1] beta_intercept;
  // detection parameters
  vector[det_ncb] beta_det;
  //real<lower=0> theta;
  // temporal availability parameters
  real<lower=0, upper=1> activ;
  real reciprocal_phi;
  // local parameters for horseshoe prior
  vector[m_psi-1] zb;
  vector<lower=0>[m_psi-1] hs_local[2];
  // horseshoe shrinkage parameters
  real<lower=0> hs_global[2];
  real<lower=0> hs_c2;
}

transformed parameters {
  // distance parameters
  vector[n_site] eps_site;
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
  // real log_theta = log(theta);
  vector[n_site] log_lambda_psi;
  // negbin dispersion
  real phi;
  phi = 1. / reciprocal_phi;
  vector[m_psi-1] beta_covs;
  vector[m_psi] beta_psi;
  // vector[n_site] logit_psi;
  // vector[n_site] log_psi;
  // vector[n_site] log1m_psi;// site random effects

  beta_covs = horseshoe(zb, hs_local, hs_global, hs_scale_global, hs_scale_slab^2 * hs_c2);
  beta_psi = append_row(beta_intercept, beta_covs);

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
  log_lambda_psi[n] = X_psi[n,] * beta_psi + eps_site[n];
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
  eps_site ~ student_t(4, 0, 1);
  eps_ngs ~ uniform(0, 1); // prior for group size effect
  beta_intercept ~ normal(0, 2); // prior for poisson model
  activ ~ beta(bshape, bscale);  //informative prior
  site_sd ~ normal(0, 1);
  site_raw ~ std_normal();
  reciprocal_phi ~ cauchy(0., reciprocal_phi_scale);
  //log_theta ~ normal(2,2);

  for(n in 1:n_site) {
  for(j in 1:n_gs) {
  n_obs[n,j] ~ neg_binomial_2(lambda[n,j], phi);
  y[n,,j] ~ multinomial_logit(to_vector(log_p_raw[n,,j]));
  }

}

   // priors including all constants: poisson
  target += normal_lpdf(zb | 0, 1);
  target += normal_lpdf(hs_local[1] | 0, 1)
    - 101 * log(0.5);
  target += inv_gamma_lpdf(hs_local[2] | 0.5 * hs_df, 0.5 * hs_df);
  target += normal_lpdf(hs_global[1] | 0, 1)
    - 1 * log(0.5);
  target += inv_gamma_lpdf(hs_global[2] | 0.5 * hs_df_global, 0.5 * hs_df_global);
  target += inv_gamma_lpdf(hs_c2 | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
}

generated quantities {
  array[n_site,n_gs] real n_obs_pred;
  array[n_site, n_gs] real n_obs_true;
  array[n_site] real N_site;
  array[n_site] real N_site_pred;
  array[n_site, max_int_dist+1] real DetCurve;
  array[n_site, n_gs] real log_lik;
  vector[n_site] Site_lambda;
  vector[n_site] psi;
  array[npc] real pred;
  real Nhat;


for(n in 1:n_site) {
  for(j in 1:n_gs) {
  log_lik[n,j] = multinomial_logit_lpmf(y[n,,j] | to_vector(log_p_raw[n,,j])); //for loo
  n_obs_true[n, j] = gs[j] * (neg_binomial_2_rng(exp(log_lambda_psi[n] + log(eps_ngs[j])), phi));
  n_obs_pred[n,j] = gs[j] * (neg_binomial_2_rng(exp(log_lambda_psi[n] + log_p[n,j] + log_activ + log(eps_ngs[j])), phi));
    }
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

for(i in 1:npc) {
  pred[i] = neg_binomial_2_rng(exp(X_pred_psi[i,] * beta_psi), phi) * prop_pred[i];
}
Nhat = sum(pred);
}
