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

  //transect level information
  int<lower=1> trans;                  // total number of transects across all sites for all methods
  int<lower = 0, upper = 1> y2[trans]; // transect based binomial detections
  int<lower = 0, upper = trans> start_idx[n_site];
  int<lower = 0, upper = trans> end_idx[n_site];
  int<lower=1> trans_det_ncb;           // number of covariates for transect detection model
  matrix[trans, trans_det_ncb] trans_pred_matrix; // transect detection model matrix
  int<lower=1>  n_max[n_site]; // max for poisson RN

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
 // abundance parameters
  simplex[n_gs] eps_ngs; // random variation in group size
  vector[1] beta_intercept;
  // detection parameters
  vector[det_ncb] beta_det;
  //real<lower=0> theta;
  // transect detection parameters
  vector[trans_det_ncb] beta_trans_det;
  // temporal availability parameters
  real<lower=0, upper=1> activ;
  real reciprocal_phi;
  // local parameters for horseshoe prior
  vector[m_psi-1] zb;
  vector<lower=0>[m_psi-1] hs_local[2];
  // horseshoe shrinkage parameters
  real<lower=0> hs_global[2];
  real<lower=0> hs_c2;
  // GP params
  real<lower=0> rho;    // GP length scale
  real<lower=0> alpha;  // GP marginal sd
  vector[n_site] eta;
}

transformed parameters {
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
  // real log_theta = log(theta);
  vector[trans] logit_trans_p = trans_pred_matrix * beta_trans_det; // observation process model
  vector[trans] r = inv_logit(logit_trans_p);
  vector[n_site] log_lambda_psi;
  // negbin dispersion
  real phi;
  phi = 1. / reciprocal_phi;
  vector[m_psi-1] beta_covs;
  vector[m_psi] beta_psi;
  // vector[n_site] logit_psi;
  // vector[n_site] log_psi;
  // vector[n_site] log1m_psi;// site random effects
  // GP params
  matrix[n_site, n_site] cov = gp_exp_quad_cov(coords, alpha, rho) +
                                diag_matrix(rep_vector(1e-9, n_site));
  matrix[n_site, n_site] L_cov = cholesky_decompose(cov);
  vector[n_site] gp_predict = L_cov * eta;

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

// define log lambda
  log_lambda_psi[n] = X_psi[n,] * beta_psi + gp_predict[n];
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
  beta_intercept ~ normal(0, 2); // prior for poisson model
  beta_trans_det ~ normal(0, 1); // beta for transect detection
  activ ~ beta(bshape, bscale);  //informative prior
  reciprocal_phi ~ cauchy(0., reciprocal_phi_scale);
  //log_theta ~ normal(2,2);
  // GP priors
  eta ~ std_normal();
  rho ~ inv_gamma(72.2999, 20.2811);
  alpha ~ normal(3, 1);

  for(n in 1:n_site) {
  for(j in 1:n_gs) {
  n_obs[n,j] ~ neg_binomial_2(lambda[n,j], phi);
  y[n,,j] ~ multinomial_logit(to_vector(log_p_raw[n,,j]));
  }

  // Royle-Nichols implementation in STAN (looping over possible discrete values of N)
// https://discourse.mc-stan.org/t/royle-and-nichols/14150
// https://discourse.mc-stan.org/t/identifiability-across-levels-in-occupancy-model/5340/2
if (n_survey[n] > 0) {
  vector[n_max[n] - any_seen[n] + 1] lp;
// seen
    if(any_seen[n] == 0){ // not seen
      lp[1] = neg_binomial_2_lpmf(0 | exp(log_lambda_psi[n]), phi);
    }
// not seen
// lp 1 simplification (not necessary)
    else lp[1] = neg_binomial_2_lpmf(1 | exp(log_lambda_psi[n]), phi) +
    bernoulli_lpmf(y2[start_idx[n]:end_idx[n]] | r[start_idx[n]:end_idx[n]]);
     // loop through possible values for maximum count (km2)
    for (j in 2:(n_max[n] - any_seen[n] + 1)){
      lp[j] = neg_binomial_2_lpmf(any_seen[n] + j - 1 | exp(log_lambda_psi[n]), phi)
      + bernoulli_lpmf(y2[start_idx[n]:end_idx[n]] | 1 - (1 - r[start_idx[n]:end_idx[n]])^(any_seen[n] + j - 1));
    }
    target += log_sum_exp(lp);
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
  // array[npc] real pred;
  // real Nhat;


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

// for(i in 1:npc) {
//   pred[i] = neg_binomial_2_rng(exp(X_pred_psi[i,] * beta_psi), phi) * prop_pred[i];
// }
// Nhat = sum(pred);
}
