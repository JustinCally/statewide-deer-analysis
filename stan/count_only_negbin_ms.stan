functions {

    /* Half-normal function
   * Args:
   *   sigma: sigma
   *   midpoints: midpoints
   * Returns:
   *   detection probability
   */
    vector halfnorm(real sigma, vector midpoints) {
      int bins = rows(midpoints);
      vector[bins] p_raw; // detection probability

      p_raw = exp(-(midpoints)^2/(2*sigma^2));
      return p_raw;
    }

      /* Hazard function
   * Args:
   *   sigma: sigma
   *   theta: theta
   *   midpoints: midpoints
   * Returns:
   *   detection probability
   */
    vector hazard(real sigma, real theta, vector midpoints) {
      int bins = rows(midpoints);
      vector[bins] p_raw; // detection probability

      p_raw = 1 - exp(-(midpoints/sigma)^(-theta));
      return p_raw;
    }

  vector prob_dist(real sigma, real theta, int keyfun, vector midpoints){
  int bins = rows(midpoints);
  vector[bins] out; // detection probability

  if(keyfun == 0){
    out = halfnorm(sigma, midpoints);
  } else if(keyfun == 1){
    out = hazard(sigma, theta, midpoints);
  }
  return out;
  }
}
data {
  int<lower=0> N;                      // number of observations
  int<lower=0> S;                      // number of species
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
  array[n_site, n_gs, S] int n_obs;      //number of observations
  array[n_site, n_distance_bins, n_gs] int y; // observations matrix

  // summary of whether species is known to be present at each site
  int<lower = 0, upper = 1> any_seen[S, n_site];

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

  // key function, 0 = halfnorm
  int keyfun;
}

transformed data {
  // vector[n_distance_bins] pi; // availability for point
  vector[n_site] survey_area;
  array[S, n_site] real cam_seen;
    for(i in 1:n_site) {
      survey_area[i] = theta_frac * effort[i] * pi() * (max_distance/1000)^2;
      for(s in 1:S) {
      cam_seen[s,i] = sum(n_obs[i,,s]);
      }
    }

}

parameters {
 // abundance parameters
  //array[n_site] simplex[n_gs] eps_ngs; // random variation in group size
  array[S] vector[m_psi] beta_psi;
  vector[det_ncb] beta_det;
  real log_theta;
  // temporal availability parameters
  real<lower=0, upper=1> activ;
  // bioregion RE
  array[S] real<lower=0> bioregion_sd;
  array[S] vector[np_bioreg] bioregion_raw;
  // site_sd
  array[S] real<lower=0> site_sd;
  array[S] vector[n_site] site_raw;
  // eps group size params
  array[S] vector[n_gs] zeta;
  array[S] matrix[n_site, n_gs] eps_raw;
  array[S] real<lower=0> grp_sd;
  // od
  array[S] real od_mu;
}

transformed parameters {
  // random effects
  array[S] vector[np_bioreg] eps_bioregion; // bioregion random effect
  array[S] vector[n_site] site_re; // bioregion random effect
  // distance parameters
  array[n_site] real log_sigma;
  array[n_site] real sigma;
  array[n_site] vector[n_distance_bins] p_raw; // detection probability
  array[n_site, n_distance_bins, n_gs] real p_raw_scale; // detection probability for point independence model and accounting for availability
  array[n_site, n_distance_bins, n_gs] real<upper=0> log_p_raw; // log detection probability accounting for availability
  array[n_site, n_gs] real log_p; // log probability with summed exponential for the multinomial logit model
  array[n_site, n_gs] real<lower=0,upper=1> p; // log probability with summed exponential for the multinomial logit model
  // abundance parameters
  array[S, n_site, n_gs] real<lower=0> lambda;
  // activity parameters
  real log_activ = log(activ);
  array[S] vector[n_site] log_lambda_psi;
  // eps group size
  array[S, n_site] simplex[n_gs] eps_ngs; // random variation in group size
  array[S, n_site] vector[n_gs] epsi;
  array[S] matrix[n_site, n_gs] eps_site;
  real<lower=0> theta = exp(log_theta);
  array[S] real<lower=0> od; // bioregion random effect

  for(s in 1:S) {
    for(b in 1:np_bioreg) {
    eps_bioregion[s,b] = bioregion_sd[s] * bioregion_raw[s,b];
  }
  od[s] = 1/od_mu[s];
  }



for(n in 1:n_site) {
  log_sigma[n] = det_model_matrix[n,] * beta_det;
  sigma[n] = exp(log_sigma[n]);
  p_raw[n] = prob_dist(sigma[n], theta, keyfun, midpts);
  for (i in 1:n_distance_bins) {
      // assuming a half-normal detection fn from line
      // p_raw[n,i] = exp(-(midpts[i])^2/(2*sigma[n]^2)); // half-normal function (pg 419 of AHM)
      // hazard function
      // p_raw[n,i] = 1 - exp(-(midpts[i]/theta)^(-sigma[n])); // hazard function (pg 507 of AHM)
    for(j in 1:n_gs) {
      p_raw_scale[n,i,j] = p_raw[n,i]*pa[j,i]; //  pr(animal occurs and is detected in bin i)
      log_p_raw[n,i,j] = log(p_raw_scale[n,i,j]);
      }
  }

  for(s in 1:S) {
  site_re[s,n] = site_sd[s] * site_raw[s,n];
  log_lambda_psi[s,n] = X_psi[n,] * beta_psi[s] + eps_bioregion[s,site_bioreg[n]] + site_re[s,n];

  for(j in 1:n_gs) {
    eps_site[s, n,j] = grp_sd[s] * eps_raw[s,n,j];
    epsi[s,n,j] = exp(zeta[s,j] + eps_site[s,n,j]);
  }
  eps_ngs[s,n,] = epsi[s,n,]/sum(epsi[s,n,]);

  }

  for(j in 1:n_gs) {
    log_p[n,j] = log_sum_exp(log_p_raw[n,,j]);
    p[n,j] = exp(log_p[n,j]);
    // model site abundance
    for(s in 1:S) {
    lambda[s,n,j] = exp(log_lambda_psi[s,n] + log_p[n,j] + log_activ + log(eps_ngs[s,n,j])) .* survey_area[n];
    }
  }
    }
}

model {

  for(s in 1:S) {
    beta_psi[s] ~ normal(0, 3); // prior for intercept in poisson model
    bioregion_sd[s] ~ normal(0,2);
    bioregion_raw[s,] ~ normal(0,1);
    to_vector(eps_raw[s,,]) ~ std_normal();
    grp_sd[s] ~ normal(0, 1);
    zeta[s,] ~ normal(0, 2);
    // od
    od_mu[s] ~ normal(0,0.1);
    site_sd[s] ~ normal(0,2);
    site_raw[s,] ~ normal(0,1);
  }
  beta_det ~ normal(0, 4); // prior for sigma
  activ ~ beta(bshape, bscale);  //informative prior
  log_theta ~ normal(0,2); // prior for theta

  for(n in 1:n_site) {
  for(j in 1:n_gs) {
  for(s in 1:S) {
  target += neg_binomial_2_lpmf(n_obs[n,j,s] | lambda[s,n,j], od[s]);
        }
  y[n,,j] ~ multinomial_logit(to_vector(log_p_raw[n,,j]));
  }
}
}

generated quantities {
  array[S, n_site,n_gs] real n_obs_pred;
  array[S, n_site, n_gs] real n_obs_true;
  array[S, n_site] real N_site;
  array[S, n_site] real N_site_pred;
  array[n_site, max_int_dist+1] real DetCurve;
  array[n_site, n_gs] real log_lik1;
  array[S, n_site, n_gs] real log_lik2;
  array[n_site, n_gs] real log_lik2_site;
  array[S, n_site] real log_lik2_species;
  array[n_site] real log_lik;
  array[n_site] real log_lik_det;
  //array[S, n_site] real Site_lambda;
  real av_gs[S];
  array[S] simplex[n_gs] eps_gs_ave;
  array[S, npc] real pred;
  //array[S, npc] real pred_trunc;
  array[S, np_reg] real Nhat_reg;
  array[S, np_bioreg] real Nhat_bioreg;
  // array[np_reg] real Nhat_reg_design;
  real Nhat[S];
  //real Nhat_trunc[S];
  real Nhat_sum;
  int trunc_counter[S];
  trunc_counter[S] = 0;

for(s in 1:S) {
  eps_gs_ave[s] = exp(zeta[s])/sum(exp(zeta[s]));
  av_gs[s] = sum(gs .* eps_gs_ave[s]); //average group size
}


for(n in 1:n_site) {
  for(j in 1:n_gs) {
  log_lik1[n,j] = multinomial_logit_lpmf(y[n,,j] | to_vector(log_p_raw[n,,j])); //for loo
  for(s in 1:S) {
  log_lik2[s,n,j] = neg_binomial_2_lpmf(n_obs[n,j,s] | lambda[s,n,j], od[s]); //for loo
  n_obs_true[s, n, j] = gs[j] * (neg_binomial_2_log_rng(log_lambda_psi[s,n] + log(eps_ngs[s,n,j]), od[s]));
  n_obs_pred[s, n,j] = gs[j] * neg_binomial_2_rng(lambda[s,n,j], od[s]);
  }
  log_lik2_site[n, j] = log_sum_exp(log_lik2[,n,j]);
    }
    // get loglik on a site level
    log_lik_det[n] = log_sum_exp(log_lik1[n,]);
    log_lik[n] = log_sum_exp(log_lik_det[n],
    log_sum_exp(log_lik2_site[n,]));
      for(s in 1:S) {
    //Site_lambda[s,n] = exp(log_lambda_psi[s,n]);
    log_lik2_species[s, n] = log_sum_exp(log_lik2[s,n,]);
    N_site[s,n] = sum(n_obs_true[s,n,]);
    N_site_pred[s,n] = sum(n_obs_pred[s,n,]);
      }

  // loop over distance bins
  if(keyfun == 0) {
  for(j in 0:max_int_dist) { // get DS predictions for distance 0 to max bin distance
    DetCurve[n, j+1] = exp(-(j+0.5)^2/(2*sigma[n]^2)); // half normal
    }
  } else if(keyfun == 1) {
    for(j in 0:max_int_dist) { // get DS predictions for distance 0 to max bin distance
    DetCurve[n, j+1] =  1 - exp(-((j+0.5)/sigma[n])^(-theta)); //hazard rate
    }
  }
}

for(s in 1:S) {
for(i in 1:np_reg) Nhat_reg[s,i] = 0;
for(i in 1:np_bioreg) Nhat_bioreg[s,i] = 0;

for(i in 1:npc) {
  pred[s,i] = neg_binomial_2_log_rng(X_pred_psi[i,] * beta_psi[s] + eps_bioregion[s, pred_bioreg[i]], od[s]) * prop_pred[i] * av_gs[s]; //offset
  if(pred[s,i] > max(N_site[s,])) {
    trunc_counter[s] += 1;
  }
  // for Hog Deer limit to Gippsland region
  if(pred_reg[i] != 2 && s == 4) {
    pred[s,i] = 0;
  }
  Nhat_reg[s,pred_reg[i]] += pred[s,i];
  Nhat_bioreg[s,pred_bioreg[i]] += pred[s,i];
}
Nhat[s] = sum(pred[s,]);
//Nhat_trunc[s] = sum(pred_trunc[s,]);
}
Nhat_sum = sum(Nhat);
}
