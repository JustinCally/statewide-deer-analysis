functions {
  /* Half-normal function */
  vector halfnorm(real sigma, vector midpoints) {
    int bins = rows(midpoints);
    vector[bins] p_raw;
    p_raw = exp( - square(midpoints) / (2 * square(sigma)) );
    return p_raw;
  }

  /* Hazard function */
  vector hazard(real sigma, real theta, vector midpoints) {
    int bins = rows(midpoints);
    vector[bins] p_raw;
    p_raw = 1 - exp( - pow( midpoints / sigma, -theta ) );
    return p_raw;
  }

  vector prob_dist(real sigma, real theta, int keyfun, vector midpoints){
    int bins = rows(midpoints);
    vector[bins] out;
    if (keyfun == 0) {
      out = halfnorm(sigma, midpoints);
    } else if (keyfun == 1) {
      out = hazard(sigma, theta, midpoints);
    }
    return out;
  }

  /**
   * Weighted mean group size with detection adjustment
   */
  real weighted_prob_mean(int n_gs, vector gs, vector freqs, vector det_p) {
    vector[n_gs] wts = freqs ./ det_p;
    real numerator = dot_product(wts, gs);
    real denominator = sum(wts);
    return numerator / denominator;
  }
}

data {
  int<lower=0> N;                      // number of observations (not used explicitly)
  int<lower=0> S;                      // number of species
  real delta;                          // bin width
  int<lower=1> n_site;                 // sites
  int<lower=1> n_distance_bins;        // distance bins
  int<lower=1> n_gs;                   // number of group sizes
  vector[n_gs] gs;                     // group sizes
  vector[n_distance_bins] midpts;      // midpoints
  real<lower=1> max_distance;          // truncation distance (m)
  int<lower=1> max_int_dist;           // max distance integer
  real<lower=0> theta_frac;            // fraction of camera view
  array[n_site] int effort;            // effort (integer)
  array[n_site, n_gs, S] int n_obs;    // counts by group size
  array[n_site, n_distance_bins, n_gs] int y; // DS histograms

  // presence–absence surveys mapped onto a flat trans index
  int<lower=1> trans;
  array[S, trans] int<lower=0, upper=1> y2;
  array[n_site] int<lower=0, upper=trans> start_idx;
  array[n_site] int<lower=0, upper=trans> end_idx;

  // number of (camera) surveys at each site (kept for compatibility)
  array[n_site] int<lower=0> n_survey;

  // availability prior
  real<lower=0> bshape;
  real<lower=0> bscale;

  // detection (distance) covariates
  int<lower=0> det_ncb;
  matrix[n_site, det_ncb] det_model_matrix;
  array[n_gs, n_distance_bins] real pa;

  // abundance/occupancy design
  int<lower=1> m_psi;
  matrix[n_site, m_psi] X_psi;       // first column is intercept

  int<lower=1> m_psi_abundance;
  matrix[n_site, m_psi_abundance] X_psi_abundance;       // first column is intercept

  // transect detection model for presence–absence
  int<lower=1> trans_det_ncb;
  matrix[trans, trans_det_ncb] trans_pred_matrix;

  // prediction data
  int<lower=1> npc;
  matrix[npc, m_psi] X_pred_psi;
  matrix[npc, m_psi_abundance] X_pred_psi_abundance;
  vector[npc] prop_pred;
  int<lower=1> np_reg;
  array[n_site] int<lower=1> site_reg;
  array[npc] int<lower=1> pred_reg;

  // bioregion RE
  int<lower=1> np_bioreg;
  array[n_site] int<lower=1> site_bioreg;
  array[npc] int<lower=1> pred_bioreg;

  // key function selector
  int keyfun;

  // ---- Regularized Horseshoe hyperparameters (not used now, kept for compatibility) ----
  real<lower=0> hs_df;            // slab df
  real<lower=0> hs_scale;         // slab scale
  real<lower=0> hs_global_scale;  // tau0
}

transformed data {
  vector[n_site] survey_area;
  array[S, n_site] real cam_seen;
  array[S] vector[n_gs] n_freqs;

  for (i in 1:n_site) {
    survey_area[i] = theta_frac * effort[i] * pi() * square(max_distance / 1000);
    for (s in 1:S) cam_seen[s, i] = sum(n_obs[i, , s]);
  }
  for (j in 1:n_gs) {
    for (s in 1:S) {
      real acc = 0;
      for (i in 1:n_site) acc += n_obs[i, j, s];
      n_freqs[s, j] = acc;
    }
  }
}

parameters {
  // ---- RHS for X_psi (first column = intercept, NOT regularised) ----
  array[S] vector[m_psi] beta_psi;
  array[S] vector[m_psi_abundance] beta_psi_abundance;
  real<lower=0> sigma_beta;  // global scale for Laplace on slopes

  // detection / abundance intercept
  vector[det_ncb] beta_det;                 // distance detection (sigma model)
  array[S] real beta_occ;                   // log-mean intercepts (for counts)
  real log_theta;                           // hazard parameter (if keyfun==1)

  // transect detection for presence–absence
  vector[trans_det_ncb] beta_trans_det;

  // temporal availability
  real<lower=0, upper=1> activ;

  // group-size mixture (logistic-normal over group sizes)
  array[S] vector[n_gs] zeta;
  array[S] matrix[n_site, n_gs] eps_raw;
  array[S] real<lower=0> grp_sd;

  // bioregion random effects
  array[S] real<lower=0> bioregion_sd;
  array[S] vector[np_bioreg] bioregion_raw;

  // NB2 dispersion (per species), on log-scale for regularisation
  vector[S] log_phi;
}

transformed parameters {
  // random effects
  array[S] vector[np_bioreg] eps_bioregion; // bioregion random effect

  // detection (distance)
  array[n_site] real log_sigma;
  array[n_site] real sigma;
  array[n_site] vector[n_distance_bins] p_raw;
  array[n_site, n_distance_bins, n_gs] real log_p_raw;
  array[n_site, n_gs] real log_p;
  array[n_site, n_gs] real<lower=0, upper=1> p;

  // NB2 mean (on log scale)
  array[S, n_site, n_gs] real log_lambda;

  // presence–absence detection
  vector[trans] logit_trans_p = trans_pred_matrix * beta_trans_det;
  vector<lower=0,upper=1>[trans] r = inv_logit(logit_trans_p);

  // site-level linear predictors
  array[S] vector[n_site] eta;             // log-mean (abundance) base
  array[S] vector[n_site] log_lambda_psi;  // NB2 log-mean (base)
  array[S] vector[n_site] psi;             // Pr(non-structural) for ZIP

  // group-size mixture per site
  array[S, n_site] simplex[n_gs] eps_ngs;

  real<lower=0> theta = exp(log_theta);
  real log_activ = log(activ);
  vector[n_gs] p_mean;
  vector<lower=0>[S] phi;                  // NB2 dispersion

  // bioregion RE & NB dispersion
  for (s in 1:S) {
    eps_bioregion[s] = bioregion_sd[s] * bioregion_raw[s];
    phi[s]           = exp(log_phi[s]);    // ensure phi > 0
  }

  // detection kernel
  for (n in 1:n_site) {
    log_sigma[n] = det_model_matrix[n,] * beta_det;
    sigma[n]     = exp(log_sigma[n]);
    p_raw[n]     = prob_dist(sigma[n], theta, keyfun, midpts);

    for (j in 1:n_gs) {
      for (i in 1:n_distance_bins) {
        // log(p_raw * pa)
        log_p_raw[n, i, j] = log(p_raw[n, i]) + log(pa[j, i]);
      }
      log_p[n, j] = log_sum_exp(log_p_raw[n, , j]);
      p[n, j]     = exp(log_p[n, j]);
    }
  }

  // site-level linear predictors, ZIP params, and group-size mixture
  for (s in 1:S) {
    for (n in 1:n_site) {
      // log-mean for counts (abundance component)
      eta[s,n]            = beta_occ[s] + eps_bioregion[s, site_bioreg[n]] + X_psi_abundance[n,] * beta_psi_abundance[s];
      log_lambda_psi[s,n] = eta[s,n];  // base NB2 log-mean

      // ZIP mixing prob
      psi[s,n] = inv_logit( X_psi[n,] * beta_psi[s] );

      {
        vector[n_gs] epsi_local;
        for (j in 1:n_gs) {
          // logistic-normal over group sizes
          epsi_local[j] = zeta[s, j] + grp_sd[s] * eps_raw[s, n, j];
        }
        epsi_local       = exp(epsi_local);
        eps_ngs[s, n]    = epsi_local / sum(epsi_local);
      }
    }
  }

  // per-group size log mean (NB2)
  for (s in 1:S)
    for (n in 1:n_site)
      for (j in 1:n_gs)
        log_lambda[s,n,j] =
          log_lambda_psi[s,n]
          + log_p[n,j]
          + log_activ
          + log(eps_ngs[s,n,j])
          + log(survey_area[n]);

  // p_mean for av_gs
  for (j in 1:n_gs) {
    real acc = 0;
    for (n in 1:n_site) acc += p[n, j];
    p_mean[j] = acc / n_site;
  }
}

model {
  // priors
  // global shrinkage for slopes
  sigma_beta ~ exponential(1);

  for (s in 1:S) {
    // group-size random effects
    to_vector(eps_raw[s]) ~ std_normal();
    grp_sd[s] ~ normal(0, 1);
    zeta[s]   ~ normal(0, 2);

    // intercept (first column of X_psi) not regularised
    beta_psi[s, 1] ~ normal(0, 3);
    beta_psi_abundance[s] ~ normal(0,2);

    // remaining coefficients with Laplace prior
    if (m_psi > 1) {
      beta_psi[s, 2:m_psi] ~ double_exponential(0, sigma_beta);
    }

    // bioregion random effects
    bioregion_raw[s] ~ std_normal();
    bioregion_sd[s]  ~ normal(0, 1);

    // NB2 dispersion (regularised to avoid extreme overdispersion)
    // log_phi ~ normal(log(5), 0.5) -> phi mostly ~ 2–12
    log_phi[s] ~ normal(log(5), 0.5);
  }

  beta_trans_det ~ normal(0, 2);
  beta_occ       ~ normal(0, 2);
  beta_det       ~ normal(0, 4);
  activ          ~ beta(bshape, bscale);
  log_theta      ~ normal(0, 2);

  // --- Likelihood ---

  // Distance histogram
  for (n in 1:n_site)
    for (j in 1:n_gs)
      y[n, , j] ~ multinomial_logit( to_vector( log_p_raw[n, , j] ) );

  // ZIP–NB2 for counts: mixture on each (s,n,j) using log_lambda and phi[s]
  for (s in 1:S) {
    for (n in 1:n_site) {
      real log_psi   = log(psi[s,n]);
      real log1m_psi = log1m(psi[s,n]);

      for (j in 1:n_gs) {
        int k = n_obs[n, j, s];
        if (k == 0) {
          real nb0 = neg_binomial_2_log_lpmf(0 | log_lambda[s,n,j], exp(log_phi[s]));
          target += log_sum_exp(log_psi + nb0, log1m_psi);
        } else {
          target += log_psi +
                    neg_binomial_2_log_lpmf(k | log_lambda[s,n,j], exp(log_phi[s]));
        }
      }
    }
  }

  // Presence–absence as binomial detection with imperfect detection
  // For site n, segment trans indices [start_idx[n] : end_idx[n]].
  for (s in 1:S) {
    for (n in 1:n_site) {
      int a = start_idx[n];
      int b = end_idx[n];
      if (b >= a && a >= 1) {
        // present: standard Bernoulli with detection prob r over that segment
        real present_ll = bernoulli_lpmf( y2[s, a:b] | r[a:b] );

        // absent: only possible if all zeros
        int sum_y = 0;
        for (t in a:b) sum_y += y2[s, t];
        real absent_ll = (sum_y == 0) ? 0 : negative_infinity();

        // mixture with same psi (Pr(non-structural))
        target += log_mix( psi[s,n], present_ll, absent_ll );
      }
    }
  }
}

generated quantities {
  array[S, n_site, n_gs] real n_obs_pred;
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

  array[S] real av_gs;
  array[S] simplex[n_gs] eps_gs_ave;

  array[S, npc] real pred;
  array[S, np_reg] real Nhat_reg;
  array[S] real Nhat;
  real Nhat_sum;
  array[S] int trunc_counter;

  // prep
  for (s in 1:S) {
    eps_gs_ave[s] = exp(zeta[s]) / sum(exp(zeta[s]));
    av_gs[s]      = weighted_prob_mean(n_gs, gs, n_freqs[s], p_mean);
    trunc_counter[s] = 0;
  }

  // log-lik pieces for LOO (counts + distance hist)
  for (n in 1:n_site) {
    for (j in 1:n_gs) {
      log_lik1[n, j] = multinomial_logit_lpmf(
        y[n, , j] | to_vector( log_p_raw[n, , j] )
      );

      for (s in 1:S) {
        int k = n_obs[n, j, s];
        if (k == 0) {
          real nb0 = neg_binomial_2_log_lpmf(0 | log_lambda[s,n,j], exp(log_phi[s]));
          log_lik2[s, n, j] = log_mix( psi[s,n], nb0, 0 );
        } else {
          log_lik2[s, n, j] =
            log(psi[s,n]) +
            neg_binomial_2_log_lpmf(k | log_lambda[s,n,j], exp(log_phi[s]));
        }
      }
      log_lik2_site[n, j] = log_sum_exp( log_lik2[, n, j] );
    }

    log_lik_det[n] = log_sum_exp( log_lik1[n,] );
    log_lik[n]     = log_sum_exp( log_lik_det[n],
                                  log_sum_exp( log_lik2_site[n,] ) );
  }

  // per-species, per-site log-liks
  for (s in 1:S) {
    for (n in 1:n_site) {
      log_lik2_species[s, n] = log_sum_exp( log_lik2[s, n,] );
    }
  }

  // Posterior predictive counts at sites (ZIP–NB2)
  for (n in 1:n_site) {
    for (s in 1:S) {
      real z_draw = bernoulli_rng( psi[s, n] );
      for (j in 1:n_gs) {
        real log_lambda_true = log_lambda_psi[s,n] + log(eps_ngs[s,n,j]);
        real log_lambda_pred = log_lambda[s,n,j];

        if (z_draw == 1) {
          n_obs_true[s, n, j] =
            gs[j] * neg_binomial_2_log_rng( log_lambda_true, exp(log_phi[s]) );
          n_obs_pred[s, n, j] =
            gs[j] * neg_binomial_2_log_rng( log_lambda_pred, exp(log_phi[s]) );
        } else {
          n_obs_true[s, n, j] = 0;
          n_obs_pred[s, n, j] = 0;
        }
      }
      N_site[s, n]      = sum( n_obs_true[s, n, ] );
      N_site_pred[s, n] = sum( n_obs_pred[s, n, ] );
    }
  }

  // Detection curve (for plotting)
  for (n in 1:n_site) {
    if (keyfun == 0) {
      for (j in 0:max_int_dist)
        DetCurve[n, j+1] =
          exp( - square(j + 0.5) / (2 * square(sigma[n])) );
    } else if (keyfun == 1) {
      for (j in 0:max_int_dist)
        DetCurve[n, j+1] =
          1 - exp( - pow((j + 0.5) / sigma[n], -exp(log_theta)) );
    }
  }

  // Predictions on grid (ZIP–NB2 with offsets), safe RNG guard
  for (s in 1:S) {
    for (i in 1:np_reg) Nhat_reg[s, i] = 0;

    for (i in 1:npc) {
      if (prop_pred[i] <= 0) {
        pred[s, i] = 0;
      } else {
        real eta_pred = beta_occ[s] + eps_bioregion[s, pred_bioreg[i]]+ X_pred_psi_abundance[i,] * beta_psi_abundance[s]
                        + log(prop_pred[i])
                        + log(av_gs[s]);
        real psi_pred =
          inv_logit(X_pred_psi[i,] * beta_psi[s]);
        real eta_safe = fmin(eta_pred, 20.7);
        int z_draw = bernoulli_rng(psi_pred);
        pred[s, i] = z_draw *
                     neg_binomial_2_log_rng(eta_safe, exp(log_phi[s]));
      }

      if (pred[s, i] > max(N_site[s,])) trunc_counter[s] += 1;

      // regional total
      Nhat_reg[s, pred_reg[i]] += pred[s, i];
    }
    Nhat[s] = sum(pred[s, ]);
  }
  Nhat_sum = sum(Nhat);
}

