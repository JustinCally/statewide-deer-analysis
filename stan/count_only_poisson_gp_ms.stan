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

}

parameters {
 // abundance parameters
  array[S] simplex[n_gs] eps_ngs; // random variation in group size
  array[S] vector[m_psi] beta_psi;
  // detection parameters
  vector[det_ncb] beta_det;
  //real<lower=0> theta;
  // temporal availability parameters
  real<lower=0, upper=1> activ;
  // GP params
  real<lower=0> rho;    // GP length scale
  array[S] real<lower=0> alpha;  // GP marginal sd
  array[S] vector[n_site] eta;

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
  array[S, n_site, n_gs] real<lower=0> lambda;
  // activity parameters
  real log_activ = log(activ);
  array[S] vector[n_site] log_lambda_psi;
  // GP params
  array[S] vector[n_site] gp_predict;

  for(k in 1:S) {
    gp_predict[k] = cholesky_decompose(gp_exp_quad_cov(coords, alpha[k], rho) +
                                diag_matrix(rep_vector(1e-9, n_site))) * eta[k];
    // define log lambda
   log_lambda_psi[k] = X_psi * beta_psi[k] + gp_predict[k];
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


  for(j in 1:n_gs) {
    log_p[n,j] = log_sum_exp(log_p_raw[n,,j]);
    p[n,j] = exp(log_p[n,j]);
    // model site abundance
    for(k in 1:S) {
    lambda[k,n,j] = exp(log_lambda_psi[k,n] + log_p[n,j] + log_activ + log(eps_ngs[k,j])) .* survey_area[n];
    }
  }

  }
    }

model {
  beta_det ~ normal(0, 4); // prior for sigma
  for(k in 1:S) {
  eps_ngs[k] ~ uniform(0, 1); // prior for group size effect
  beta_psi[k] ~ normal(0, 2); // prior for poisson model
  eta[k] ~ std_normal();
  alpha[k] ~ normal(0, 2);
  }
  activ ~ beta(bshape, bscale);  //informative prior
  //log_theta ~ normal(2,2);
  // GP priors
  rho ~ inv_gamma(8.91924, 1.72902);

  for(n in 1:n_site) {
  for(j in 1:n_gs) {
    for(k in 1:S) {
  n_obs[n,j,k] ~ poisson(lambda[k,n,j]);
    }
  y[n,,j] ~ multinomial_logit(to_vector(log_p_raw[n,,j]));
  }
}
}

generated quantities {
  array[S, n_site,n_gs] real n_obs_pred;
  array[S, n_site, n_gs] real n_obs_true;
  array[S, n_site] real N_site;
  array[n_site] real N_site_all;
  array[S, n_site] real N_site_pred;
  array[n_site] real N_site_pred_all;
  array[n_site, max_int_dist+1] real DetCurve;
  array[n_site, n_gs] real log_lik1;
  array[n_site, n_gs] real log_lik2;
  array[n_site] real log_lik;
  array[S, n_site] real Site_lambda;
  array[S, n_site] real psi;
  // array[npc] real pred;
  // real Nhat;


for(n in 1:n_site) {
  for(j in 1:n_gs) {
  log_lik1[n,j] = multinomial_logit_lpmf(y[n,,j] | to_vector(log_p_raw[n,,j])); //for loo
  log_lik2[n,j] =  poisson_lpmf(n_obs[n,j,] | lambda[,n,j]); //for loo
  for(k in 1:S) {
  n_obs_true[k,n,j] = gs[j] * (poisson_rng(exp(log_lambda_psi[k,n] + log(eps_ngs[k,j]))));
  n_obs_pred[k,n,j] = gs[j] * (poisson_rng(exp(log_lambda_psi[k,n] + log_p[n,j] + log_activ + log(eps_ngs[k,j])) .* survey_area[n]));
  }
    }
    // get loglik on a site level
    log_lik[n] = log_sum_exp(log_sum_exp(log_lik1[n,]), log_sum_exp(log_lik2[n,]));
    for(k in 1:S) {
    Site_lambda[k, n] = exp(log_lambda_psi[k,n]);
    N_site[k,n] = sum(n_obs_true[k,n,]);
    N_site_pred[k,n] = sum(n_obs_pred[k,n,]);
    // Occupancy probability transformation
    psi[k,n] = inv_cloglog(log_lambda_psi[k,n]);
    }
    N_site_all[n] = sum(N_site[,n]);
    N_site_pred_all[n] = sum(N_site_pred[,n]);

  // loop over distance bins
  for(j in 0:max_int_dist) { // get DS predictions for distance 0 to max bin distance
    DetCurve[n, j+1] = exp(-(j+0.5)^2/(2*sigma[n]^2)); // half normal
    // DetCurve[n, j+1] =  1 - exp(-(j+0.5/theta)^(-sigma[n])); //hazard rate
    }
}

}
