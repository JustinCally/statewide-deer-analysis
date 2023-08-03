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
  // negbinom scale
  // real reciprocal_phi_scale;

  //transect level information
  int<lower=1> trans;                  // total number of transects across all sites for all methods
  int<lower = 0, upper = 1> y2[trans]; // transect based binomial detections
  int<lower = 0, upper = trans> start_idx[n_site];
  int<lower = 0, upper = trans> end_idx[n_site];
  int<lower=1> trans_det_ncb;           // number of covariates for transect detection model
  matrix[trans, trans_det_ncb] trans_pred_matrix; // transect detection model matrix
  int<lower=1>  n_max[n_site]; // max for poisson RN

  // Prediction data
  int<lower=1> npc;                 // number of prediction grid cells
  matrix[npc, m_psi] X_pred_psi; // pred matrix
  vector[npc] prop_pred; //offset
  // bioregion RE
  int<lower=1> np_bioreg;
  int<lower=1> site_bioreg[n_site];
  int<lower=1> pred_bioreg[npc];
    // region data
  int<lower=1> np_reg;
  int<lower=1> site_reg[n_site];
  int<lower=1> pred_reg[npc];
}

transformed data {
  // vector[n_distance_bins] pi; // availability for point
  vector[n_site] survey_area;
  vector[n_site] cam_seen;
    for(i in 1:n_site) {
      survey_area[i] = theta_frac * effort[i] * pi() * (max_distance/1000)^2;
      cam_seen[i] = sum(n_obs[i,]);
    }

}

parameters {
 // abundance parameters
  array[n_site] simplex[n_gs] eps_ngs; // random variation in group size
  vector[m_psi] beta_psi;
  vector[det_ncb] beta_det;
  // transect detection parameters
  vector[trans_det_ncb] beta_trans_det;
  // temporal availability parameters
  real<lower=0, upper=1> activ;
  // bioregion RE
  real<lower=0> bioregion_sd;
  vector[np_bioreg] bioregion_raw;
}

transformed parameters {
  // random effects
  vector[np_bioreg] eps_bioregion; // bioregion random effect
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
  // lp_site for RN model
  vector[n_site] log_lambda_psi;

  for(b in 1:np_bioreg) {
    eps_bioregion[b] = bioregion_sd * bioregion_raw[b];
  }

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

  log_lambda_psi[n] = X_psi[n,] * beta_psi + eps_bioregion[site_bioreg[n]];

  for(j in 1:n_gs) {
    log_p[n,j] = log_sum_exp(log_p_raw[n,,j]);
    p[n,j] = exp(log_p[n,j]);
    // model site abundance
    lambda[n,j] = exp(log_lambda_psi[n] + log_p[n,j] + log_activ + log(eps_ngs[n,j])) .* survey_area[n];
  }

    }
}

model {
  beta_det ~ normal(0, 4); // prior for sigma
  beta_psi ~ normal(0, 3); // prior for intercept in poisson model
  activ ~ beta(bshape, bscale);  //informative prior
  bioregion_sd ~ normal(0, 2);
  bioregion_raw ~ normal(0,1);

  for(n in 1:n_site) {
  eps_ngs[n] ~ uniform(0, 1); // prior for group size effect
  for(j in 1:n_gs) {
  target += poisson_lpmf(n_obs[n,j] | lambda[n,j]);
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
  vector[n_site] av_gs_site;
  array[npc] real pred;
  array[npc] real pred_trunc;
  array[np_reg] real Nhat_reg;
  real av_gs;
  // array[np_reg] real Nhat_reg_design;
  real Nhat;
  real Nhat_trunc;
  int trunc_counter;
  trunc_counter = 0;


for(n in 1:n_site) {
  for(j in 1:n_gs) {
  log_lik1[n,j] = multinomial_logit_lpmf(y[n,,j] | to_vector(log_p_raw[n,,j])); //for loo
  log_lik2[n,j] = poisson_lpmf(n_obs[n,j] | lambda[n,j]); //for loo
  n_obs_true[n, j] = gs[j] * (poisson_log_rng(log_lambda_psi[n] + log(eps_ngs[n,j])));
  n_obs_pred[n,j] = gs[j] * poisson_rng(lambda[n,j]);
    }
    // get loglik on a site level
    log_lik[n] = log_sum_exp(log_sum_exp(log_lik1[n,]), log_sum_exp(log_lik2[n,]));
    Site_lambda[n] = exp(log_lambda_psi[n]);
    N_site[n] = sum(n_obs_true[n,]);
    N_site_pred[n] = sum(n_obs_pred[n,]);
    // Occupancy probability transformation
    psi[n] = inv_cloglog(log_lambda_psi[n]);
    av_gs_site[n] = sum(gs .* eps_ngs[n]);

  // loop over distance bins
  for(j in 0:max_int_dist) { // get DS predictions for distance 0 to max bin distance
    DetCurve[n, j+1] = exp(-(j+0.5)^2/(2*sigma[n]^2)); // half normal
    // DetCurve[n, j+1] =  1 - exp(-(j+0.5/theta)^(-sigma[n])); //hazard rate
    }
}
  av_gs = mean(av_gs_site);

for(i in 1:np_reg) Nhat_reg[i] = 0;

for(i in 1:npc) {
  pred[i] = poisson_log_rng(X_pred_psi[i,] * beta_psi + eps_bioregion[pred_bioreg[i]]) * prop_pred[i] * av_gs; //offset
  if(pred[i] > max(N_site)) {
    pred_trunc[i] = max(N_site);
    trunc_counter += 1;
  } else {
    pred_trunc[i] = pred[i];
  } // upper limit placed at highest site estimate
  Nhat_reg[pred_reg[i]] += pred[i];
}
Nhat = sum(pred);
Nhat_trunc = sum(pred_trunc);
}
