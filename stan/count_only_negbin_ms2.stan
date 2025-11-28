functions {
  /* Half-normal function
   * Args:
   *   sigma: scale
   *   midpoints: distance bin midpoints
   * Returns:
   *   detection probability at each midpoint
   */
  vector halfnorm(real sigma, vector midpoints) {
    int bins = rows(midpoints);
    vector[bins] p_raw;
    p_raw = exp(- square(midpoints) / (2 * square(sigma)));
    return p_raw;
  }

  /* Hazard function
   * Args:
   *   sigma: scale
   *   theta: shape
   *   midpoints: distance bin midpoints
   * Returns:
   *   detection probability at each midpoint
   */
  vector hazard(real sigma, real theta, vector midpoints) {
    int bins = rows(midpoints);
    vector[bins] p_raw;
    p_raw = 1 - exp(- pow(midpoints / sigma, -theta));
    return p_raw;
  }

  // Wrapper to choose key function
  vector prob_dist(real sigma, real theta, int keyfun, vector midpoints) {
    vector[rows(midpoints)] out;
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
  vector[n_distance_bins] midpts;      // distance bin midpoints
  real<lower=1> max_distance;          // truncation distance (m)
  int<lower=1> max_int_dist;           // max distance as integer
  real<lower=0> theta_frac;            // fraction of camera view
  array[n_site] int effort;            // effort (integer)
  array[n_site, n_gs, S] int n_obs;    // counts by group size
  array[n_site, n_distance_bins, n_gs] int y; // DS histograms

  // summary of whether species is known to be present at each site
  array[S, n_site] int<lower=0, upper=1> any_seen;

  // number of surveys at each site
  array[n_site] int<lower=0> n_survey;

  // availability prior
  real<lower=0> bshape;
  real<lower=0> bscale;

  // detection (distance) covariates
  int<lower=0> det_ncb;
  matrix[n_site, det_ncb] det_model_matrix;
  array[n_gs, n_distance_bins] real pa;

  // abundance / occupancy design
  int<lower=1> m_psi;
  matrix[n_site, m_psi] X_psi;

  // prediction data
  int<lower=1> npc;
  matrix[npc, m_psi] X_pred_psi;
  vector[npc] prop_pred;

  // bioregion random effects
  int<lower=1> np_bioreg;
  array[n_site] int<lower=1> site_bioreg;
  array[npc] int<lower=1> pred_bioreg;

  // region ids
  int<lower=1> np_reg;
  array[n_site] int<lower=1> site_reg;
  array[npc] int<lower=1> pred_reg;

  // key function selector
  int keyfun;
}

transformed data {
  vector[n_site] survey_area;
  array[S, n_site] real cam_seen;
  array[S] vector[n_gs] n_freqs;

  for (i in 1:n_site) {
    survey_area[i] = theta_frac * effort[i] * pi() * square(max_distance / 1000);
    for (s in 1:S) {
      real acc = 0;
      for (j in 1:n_gs) acc += n_obs[i, j, s];
      cam_seen[s, i] = acc;
    }
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
  // abundance (per species)
  array[S] vector[m_psi] beta_psi;

  // detection (distance)
  vector[det_ncb] beta_det;
  real log_theta;

  // temporal availability
  real<lower=0, upper=1> activ;

  // bioregion RE
  array[S] real<lower=0> bioregion_sd;
  array[S] vector[np_bioreg] bioregion_raw;

  // group-size mixture (logistic-normal)
  array[S] vector[n_gs] zeta;
  array[S] matrix[n_site, n_gs] eps_raw;
  array[S] real<lower=0> grp_sd;

  // NB2 overdispersion
  array[S] real od_mu;   // species-specific deviation
}

transformed parameters {
  // bioregion RE
  array[S] vector[np_bioreg] eps_bioregion;

  // distance kernel
  array[n_site] real log_sigma;
  array[n_site] real sigma;
  array[n_site] vector[n_distance_bins] p_raw;
  array[n_site, n_distance_bins, n_gs] real<upper=0> log_p_raw;
  array[n_site, n_gs] real log_p;
  array[n_site, n_gs] real<lower=0, upper=1> p;

  // abundance
  real log_activ = log(activ);
  array[S] vector[n_site] log_lambda_psi; // base log mean for RN & counts
  array[S, n_site, n_gs] real<lower=0> lambda;

  // group size mixture
  array[S, n_site] simplex[n_gs] eps_ngs;
  array[S, n_site] vector[n_gs] epsi;
  array[S] matrix[n_site, n_gs] eps_site;

  // dispersion
  real<lower=0> theta = exp(log_theta);
  array[S] real od;      // NB2 overdispersion for counts

  // for av_gs
  vector[n_gs] p_mean;

  // bioregion RE & species-specific dispersion
  for (s in 1:S) {
    eps_bioregion[s] = bioregion_sd[s] * bioregion_raw[s];
    od[s] = exp(od_mu[s]);
  }

  // distance kernel and abundance
  for (n in 1:n_site) {
    // detection kernel
    log_sigma[n] = det_model_matrix[n, ] * beta_det;
    sigma[n]     = exp(log_sigma[n]);
    p_raw[n]     = prob_dist(sigma[n], theta, keyfun, midpts);

    for (j in 1:n_gs) {
      for (i in 1:n_distance_bins) {
        // pr(animal occurs and is detected in bin i) given group size j
        log_p_raw[n, i, j] = log(p_raw[n][i] * pa[j, i]);
      }
      log_p[n, j] = log_sum_exp(log_p_raw[n, , j]);
      p[n, j]     = exp(log_p[n, j]);
    }

    // species-specific parts
    for (s in 1:S) {
      // base log mean abundance at site n (per species, shared across group sizes)
      log_lambda_psi[s][n] = X_psi[n, ] * beta_psi[s] + eps_bioregion[s][site_bioreg[n]];

      // group-size mixture (logistic-normal)
      for (j in 1:n_gs) {
        eps_site[s, n, j] = grp_sd[s] * eps_raw[s, n, j];
        epsi[s, n, j]     = exp(zeta[s][j] + eps_site[s, n, j]);
      }
      eps_ngs[s, n] = epsi[s, n] / sum(epsi[s, n]);
    }

    // per-group size lambda
    for (j in 1:n_gs) {
      for (s in 1:S) {
        lambda[s, n, j] =
          exp(log_lambda_psi[s][n] + log_p[n, j] + log_activ + log(eps_ngs[s, n, j]))
          * survey_area[n];
      }
    }
  }

  // p_mean for av_gs
  for (j in 1:n_gs) {
    real acc = 0;
    for (n in 1:n_site) acc += p[n, j];
    p_mean[j] = acc / n_site;
  }
}

model {
  // Priors
  for (s in 1:S) {
    beta_psi[s]    ~ normal(0, 3);
    bioregion_sd[s] ~ normal(0, 2);
    bioregion_raw[s] ~ normal(0, 1);

    to_vector(eps_raw[s]) ~ std_normal();
    grp_sd[s] ~ normal(0, 1);
    zeta[s]   ~ normal(0, 2);

    od_mu[s] ~ normal(0, 1);
  }

  beta_det       ~ normal(0, 4);
  activ          ~ beta(bshape, bscale);
  log_theta      ~ normal(0, 2);

  // Likelihood
  for (n in 1:n_site) {
    for (j in 1:n_gs) {
      // distance histogram at site n, group size j
      {
        vector[n_distance_bins] lp;
        for (i in 1:n_distance_bins)
          lp[i] = log_p_raw[n, i, j];
        y[n, , j] ~ multinomial_logit(lp);
      }

      // counts per species
      for (s in 1:S) {
        target += neg_binomial_2_lpmf(n_obs[n, j, s] | lambda[s, n, j], od[s]);
      }
    }
  }
}

generated quantities {
  // posterior predictive
  array[S, n_site, n_gs] real n_obs_pred;
  array[S, n_site, n_gs] real n_obs_true;
  array[S, n_site] real N_site;
  array[S, n_site] real N_site_pred;

  // detection curves
  array[n_site, max_int_dist + 1] real DetCurve;

  // log-lik components
  array[n_site, n_gs] real log_lik1;
  array[S, n_site, n_gs] real log_lik2;
  array[n_site, n_gs] real log_lik2_site;
  array[S, n_site] real log_lik2_species;
  array[n_site] real log_lik;
  array[n_site] real log_lik_det;

  // group size summaries
  array[S] real av_gs;
  array[S] simplex[n_gs] eps_gs_ave;

  // predictions
  array[S, npc] real pred;
  array[S, np_reg] real Nhat_reg;
  array[S, np_bioreg] real Nhat_bioreg;
  array[S] real Nhat;
  real Nhat_sum;
  array[S] int trunc_counter;

  // ---- prep ----
  for (s in 1:S) {
    eps_gs_ave[s] = exp(zeta[s]) / sum(exp(zeta[s]));
    av_gs[s]      = weighted_prob_mean(n_gs, gs, n_freqs[s], p_mean);
    trunc_counter[s] = 0;
  }

  // log-lik & posterior predictive by site
  for (n in 1:n_site) {
    // distance histograms
    for (j in 1:n_gs) {
      vector[n_distance_bins] lp;
      for (i in 1:n_distance_bins)
        lp[i] = log_p_raw[n, i, j];

      log_lik1[n, j] = multinomial_logit_lpmf(y[n, , j] | lp);

      // counts
      for (s in 1:S) {
        log_lik2[s, n, j] =
          neg_binomial_2_lpmf(n_obs[n, j, s] | lambda[s, n, j], od[s]);

        // posterior predictive at site (true underlying & observed)
        n_obs_true[s, n, j] =
          gs[j] * neg_binomial_2_log_rng(
                      log_lambda_psi[s][n] + log(eps_ngs[s, n, j]),
                      od[s]);

        // ==============================
        // FIX 1: guard neg_binomial_2_rng
        // ==============================
        {
          real mu_counts = lambda[s, n, j];
          if (mu_counts <= 0) {
            n_obs_pred[s, n, j] = 0;
          } else {
            n_obs_pred[s, n, j] =
              gs[j] * neg_binomial_2_rng(mu_counts, od[s]);
          }
        }
        // ==============================
      }

      log_lik2_site[n, j] = 0;
      for (s in 1:S)
        log_lik2_site[n, j] += log_lik2[s, n, j];
    }

    // site-level log-lik pieces
    log_lik_det[n] = 0;
    for (j in 1:n_gs)
      log_lik_det[n] += log_lik1[n, j];

    // counts over species & group sizes
    real lp_counts = 0;
    for (j in 1:n_gs)
      lp_counts += log_lik2_site[n, j];

    // RN contribution

    log_lik[n] = log_lik_det[n] + lp_counts;

    // species-specific log-lik + N_site
    for (s in 1:S) {
      real lp_counts_sn = 0;
      for (j in 1:n_gs)
        lp_counts_sn += log_lik2[s, n, j];
      log_lik2_species[s, n] = lp_counts_sn;

      N_site[s, n]      = 0;
      N_site_pred[s, n] = 0;
      for (j in 1:n_gs) {
        N_site[s, n]      += n_obs_true[s, n, j];
        N_site_pred[s, n] += n_obs_pred[s, n, j];
      }
    }

    // detection curve for plotting
    if (keyfun == 0) {
      for (j in 0:max_int_dist) {
        DetCurve[n, j + 1] =
          exp(- square(j + 0.5) / (2 * square(sigma[n])));
      }
    } else if (keyfun == 1) {
      for (j in 0:max_int_dist) {
        DetCurve[n, j + 1] =
          1 - exp(- pow((j + 0.5) / sigma[n], -theta));
      }
    }
  }

  // predictions on grid (NB2 with offsets)
  for (s in 1:S) {
    for (i in 1:np_reg)    Nhat_reg[s, i]    = 0;
    for (i in 1:np_bioreg) Nhat_bioreg[s, i] = 0;

    for (i in 1:npc) {
      real eta_pred = X_pred_psi[i, ] * beta_psi[s]
                      + eps_bioregion[s][pred_bioreg[i]];
      real mu_pred  = exp(eta_pred) * prop_pred[i] * av_gs[s];

      // ==============================
      // FIX 2: guard neg_binomial_2_rng
      // ==============================
      if (mu_pred <= 0) {
        pred[s, i] = 0;
      } else {
        pred[s, i] = neg_binomial_2_rng(mu_pred, od[s]);
      }
      // ==============================

      if (pred[s, i] > max(N_site[s, ])) {
        trunc_counter[s] += 1;
      }

      // example species-specific restriction: Hog deer only in Gippsland (region 2)
      if (pred_reg[i] != 2 && s == 4) {
        pred[s, i] = 0;
      }

      Nhat_reg[s, pred_reg[i]]       += pred[s, i];
      Nhat_bioreg[s, pred_bioreg[i]] += pred[s, i];
    }

    Nhat[s] = 0;
    for (i in 1:npc)
      Nhat[s] += pred[s, i];
  }

  Nhat_sum = 0;
  for (s in 1:S) Nhat_sum += Nhat[s];
}


