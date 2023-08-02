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
  array[n_site, n_gs] int n_obs;      //number of observations at each site and for varying group sizes
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
  int<lower = 1> m_psi; // number of covariates
  matrix[n_site, m_psi] X_psi; // model matrix

  // hs params
  real<lower=0> hs_df; // horseshoe prior param: see brms documentation of hs() for explantations of data
  real<lower=0> hs_df_global;
  real<lower=0> hs_df_slab;
  real<lower=0> hs_scale_global;
  real<lower=0> hs_scale_slab;

  //transect level information
  int<lower=1> trans;                  // total number of transects across all sites for all methods
  int<lower = 0, upper = 1> y2[trans]; // transect based binomial detections
  int<lower = 0, upper = trans> start_idx[n_site]; // first "transect" of the site, which is the camera detection
  int<lower = 0, upper = trans> end_idx[n_site]; // last transect of the site
  int<lower=1> trans_det_ncb;           // number of covariates for transect detection model
  matrix[trans, trans_det_ncb] trans_pred_matrix; // transect detection model matrix
  int<lower=1>  n_max[n_site]; // max for poisson RN

  // Prediction data
  int<lower=1> npc;                 // number of prediction grid cells
  matrix[npc, m_psi] X_pred_psi;    // pred matrix
  vector[npc] prop_pred;            // offset for prediction
  // bioregion RE
  int<lower=1> np_bioreg;           // number of bioregions
  int<lower=1> site_bioreg[n_site]; // bioregion of each site
  int<lower=1> pred_bioreg[npc];    // bioregion of each prediction grid
    // DEECA region data
  int<lower=1> np_reg;              // number of regions
  int<lower=1> site_reg[n_site];    // region of each site
  int<lower=1> pred_reg[npc];       // region of each prediction grid
}

transformed data {
  vector[n_site] survey_area;       // total survey area, calculated from snapshot moments and camera field of view
  vector[n_site] cam_seen;          // whether any have been seen at a site on camera
    for(i in 1:n_site) {
      survey_area[i] = theta_frac * effort[i] * pi() * (max_distance/1000)^2;
      cam_seen[i] = sum(n_obs[i,]);
    }

//  for (i in 1:n_distance_bins) {
//    //pi[i] = delta/max_distance; // availability function (line transect)
//    pi[i] = (2 * delta * midpts[i])/max_distance^2; // point transect
//  }

}

parameters {
 // abundance parameters
  simplex[n_gs] eps_ngs; // random variation in group size
  vector[1] beta_intercept; // intercept parameter (not subject to horseshoe)
  // detection parameters
  vector[det_ncb] beta_det; // coefficients for detection submodel
  // transect detection parameters
  vector[trans_det_ncb] beta_trans_det; // coefficients for transect submodel
  // temporal availability parameters
  real<lower=0, upper=1> activ; // average activity of deer, equates to the time they are available for detection
  // bioregion RE
  real<lower=0> bioregion_sd; // standard deviation for bioregion random effect
  vector[np_bioreg] bioregion_raw; // bioregion-level effect
  // local parameters for horseshoe prior
  vector[m_psi-1] zb;
  vector<lower=0>[m_psi-1] hs_local[2];
  // horseshoe shrinkage parameters
  real<lower=0> hs_global[2];
  real<lower=0> hs_c2;
}

transformed parameters {
  // random effects
  vector[np_bioreg] eps_bioregion; // bioregion random effect
  // distance parameters
  array[n_site] real log_sigma; // log of distance sampling detection function sigma
  array[n_site] real sigma; // distance sampling detection function
  array[n_site, n_distance_bins] real p_raw; // detection probability
  array[n_site, n_distance_bins, n_gs] real p_raw_scale; // detection probability for point independence model and accounting for availability
  array[n_site, n_distance_bins, n_gs] real<upper=0> log_p_raw; // log detection probability accounting for availability
  array[n_site, n_gs] real log_p; // log probability with summed exponential for the multinomial logit model
  array[n_site, n_gs] real<lower=0,upper=1> p; // log probability with summed exponential for the multinomial logit model
  // abundance parameters
  array[n_site, n_gs] real<lower=0> lambda; // expected number of groups for the camera surveys
  // activity parameters
  real log_activ = log(activ);
  vector[trans] logit_trans_p = trans_pred_matrix * beta_trans_det; // observation process model
  // lp_site for RN model
  vector[n_site] lp_site; // log-probability for each site based on RN loops
  vector[trans] r = inv_logit(logit_trans_p); // probability of detection for each transect observation/period, the first at each site is the camera-level detection rate
  vector[n_site] log_lambda_psi; //log lambda for estimating number of groups per site
  // hs betas
  vector[m_psi-1] beta_covs; // covariates with horseshow prior
  vector[m_psi] beta_psi; // covariates combined with intercept

  beta_covs = horseshoe(zb, hs_local, hs_global, hs_scale_global, hs_scale_slab^2 * hs_c2); // horseshoe bior placed on beta covariates
  beta_psi = append_row(beta_intercept, beta_covs); // combine intercept and covariates

  for(b in 1:np_bioreg) {
    eps_bioregion[b] = bioregion_sd * bioregion_raw[b]; // generate bioregion random effect
  }

for(n in 1:n_site) {
  log_sigma[n] = det_model_matrix[n,] * beta_det; // calculate sigma from detection function
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

  log_lambda_psi[n] = X_psi[n,] * beta_psi + eps_bioregion[site_bioreg[n]]; //site-level estimate of abundance of groups

// for various group sizes break down the estimated detection rates and the estimated camera counts
  for(j in 1:n_gs) {
    log_p[n,j] = log_sum_exp(log_p_raw[n,,j]);
    p[n,j] = exp(log_p[n,j]);
    // model site abundance taking into account true abundance, distance sampling (log_p),
    // temporal availability (log_activ), proportional group size (epgs_ngs), and spatial
    // availability/probability of camera-level detection (r[start_idx[n]])
    // offset is the survey area, which includes the snapshot moments and area in camera field of view
    lambda[n,j] = exp(log_lambda_psi[n] + log_p[n,j] + log_activ + log(eps_ngs[j]) + log(r[start_idx[n]])) .* survey_area[n];
  }

// Royle-Nichols implementation in STAN (looping over possible discrete values of N)
// https://discourse.mc-stan.org/t/royle-and-nichols/14150
// https://discourse.mc-stan.org/t/identifiability-across-levels-in-occupancy-model/5340/2
if (n_survey[n] > 0) {
  vector[n_max[n] - any_seen[n] + 1] lp;
// seen
    if(any_seen[n] == 0){ // not seen
      lp[1] = poisson_lpmf(0 | exp(log_lambda_psi[n]));
    }
// not seen
// lp 1 simplification (not necessary)
    else lp[1] = poisson_lpmf(1 | exp(log_lambda_psi[n])) +
    bernoulli_lpmf(y2[start_idx[n]:end_idx[n]] | r[start_idx[n]:end_idx[n]]);
     // loop through possible values for maximum count (km2)
    for (j in 2:(n_max[n] - any_seen[n] + 1)){
      lp[j] = poisson_lpmf(any_seen[n] + j - 1 | exp(log_lambda_psi[n]))
      + bernoulli_lpmf(y2[start_idx[n]:end_idx[n]] | 1 - (1 - r[start_idx[n]:end_idx[n]])^(any_seen[n] + j - 1));
    }
    lp_site[n] = log_sum_exp(lp);
  } else{
    lp_site[n] = 0;
  }

  }
    }

model {
  beta_det ~ normal(0, 4); // prior for sigma
  eps_ngs ~ uniform(0, 1); // prior for group size effect
  beta_intercept ~ normal(-3, 3); // prior for intercept in poisson model
  beta_trans_det ~ normal(0, 3); // beta for transect detection
  activ ~ beta(bshape, bscale);  //informative prior
  bioregion_sd ~ normal(0, 2); // prior for bioregion RE SD
  bioregion_raw ~ normal(0,1); // prior for bioregion RE effect

  for(n in 1:n_site) {
  for(j in 1:n_gs) {
  target += poisson_lpmf(n_obs[n,j] | lambda[n,j]);
  y[n,,j] ~ multinomial_logit(to_vector(log_p_raw[n,,j]));
  }
  target += lp_site[n];
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
  array[n_site, n_gs] real log_lik1;
  array[n_site, n_gs] real log_lik2;
  array[n_site] real log_lik;
  vector[n_site] Site_lambda;
  vector[n_site] psi;
  array[npc] real pred;
  array[npc] real pred_trunc;
  array[np_reg] real Nhat_reg;
  real av_gs;
  real Nhat;
  real Nhat_trunc;
  int trunc_counter;
  trunc_counter = 0;


for(n in 1:n_site) {
  for(j in 1:n_gs) {
  log_lik1[n,j] = multinomial_logit_lpmf(y[n,,j] | to_vector(log_p_raw[n,,j])); //for loo
  log_lik2[n,j] = poisson_lpmf(n_obs[n,j] | lambda[n,j]); //for loo
  n_obs_true[n, j] = gs[j] * (poisson_log_rng(log_lambda_psi[n] + log(eps_ngs[j])));
  n_obs_pred[n,j] = gs[j] * poisson_rng(lambda[n,j]);
    }
    // get loglik on a site level
    log_lik[n] = log_sum_exp(log_sum_exp(log_sum_exp(log_lik1[n,]), log_sum_exp(log_lik2[n,])), lp_site[n]);
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
  av_gs = sum(gs .* eps_ngs);

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
