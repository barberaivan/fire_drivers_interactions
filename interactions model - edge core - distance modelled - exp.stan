/*
  Conditional logistic model to classify core vs. edge pixels.
  This is like a categorical model with K = 2, where the parameters are
  fixed among categories but the predictors values vary.

  In this model the effect of distance is modelled as an atenuator of
  coefficients: b * x * exp(a * dist),
  where a < 0
*/

data {
  // Ns
  int<lower=0> N; // observations
  int<lower=0> N_pairs; // pairs of observations (N/2)
  int<lower=0> NF; // fires (groups)
  int<lower=0> N_fire[NF]; // observations by fire
  int<lower=0> K; // observation-level predictors dimension
  int<lower=0> K_dist; // number of predictor terms interacting with ditance (continuous ones)
  int<lower=0> J; // group-level predictors dimension

  // ids
  int starts[NF];  // positions in X where each fire id starts
  int ends[NF];    // or ends
  int core_ids[N_pairs]; // indices corresponding to core in X
  int edge_ids[N_pairs]; // indices corresponding to edge in X

  // predictors
  matrix[N, K] X;
  matrix[N, K_dist] X_dist; // X * distance
  vector[N] dist;  // distance
  matrix[NF, J] Z; // group-level design matrix (fwi or fire area)

  // priors
  real prior_corr_eta;
  real prior_sigma_sd; // sd for the sigma among fires
  real prior_b_sd;     // sd for all coefficients
}

transformed data {
  matrix[J, NF] Zt = Z';
}

parameters {
  // intercepts and slopes for the hiperparameters
  matrix[K, J] g_raw;
  // coefficients for distance-attenuation effect (no sign constraint!)
  vector[K_dist] d_raw;
  // random effects
  matrix[K, NF] error_raw;
  // sigma for random effects
  vector<lower = 0>[K] sigma_raw;
  // correlation matrix for hyperparameters
  cholesky_factor_corr[K] L_corr;
}

transformed parameters {
  // intercepts and slopes for the hiperparameters
  matrix[K, J] g;
  // coefficients for interaction terms with distance
  vector[K_dist] d;
  vector[K] d_extra; // extended vector, with ones
  // random effects
  matrix[K, NF] error;
  // sigma for random effects
  vector[K] sigma;
  // correlation matrix for random effects
  corr_matrix[K] Rho;
  // vcov matrix
  cov_matrix[K] V;
  matrix[K, K] V_chol;
  // coefficients matrix
  matrix[K, NF] b; // b_mean + raneff

  // linear predictor
  vector[N] linpred;
  // fitted prob
  vector[N_pairs] fitted_prob;

  g = g_raw * prior_b_sd;
  d = d_raw;// * prior_b_sd * 0.2; // normal(0, 1) prior
  d_extra = append_row(rep_vector(0, K - K_dist), d);
  sigma = sigma_raw * prior_sigma_sd;
  Rho = multiply_lower_tri_self_transpose(L_corr); // = Lcorr * Lcorr'
  V = quad_form_diag(Rho, sigma);
  V_chol = cholesky_decompose(V);
  error = V_chol * error_raw;
  b = g * Zt + error; // fixed effects + random effects

  // compute linear predictor
  for(f in 1:NF) {
    int S = N_fire[f];

    // matrix with products X * b, before summation
    matrix[S, K] Xb = X[starts[f] : ends[f], ] .* (rep_vector(1, S) * b[, f]');

    // attenuation factors to multiply with each term
    matrix[S, K] eta_log = dist[starts[f] : ends[f]] * d_extra';
    matrix[S, K] eta = exp(eta_log); // attenuation factors

    // linear predictor (the product with the ones vector is used for summation)
    linpred[starts[f] : ends[f]] = (Xb .* eta) * rep_vector(1, K);
    // linpred[starts[f] : ends[f]] = Xb * rep_vector(1, K);
  }

  // Compute probabilities
  {
    matrix[N_pairs, 2] linpred_mat;
    linpred_mat = append_col(linpred[core_ids], linpred[edge_ids]);
    for(i in 1:N_pairs) {
      fitted_prob[i] = softmax(linpred_mat[i, ]')[1]; // core probability
    }
  }

}

model {
  // priors
  to_vector(error_raw) ~ std_normal();
  to_vector(g_raw) ~ std_normal();
  d_raw ~ std_normal();
  sigma_raw ~ std_normal();
  L_corr ~ lkj_corr_cholesky(prior_corr_eta);

  // likelihood
  target += log(fitted_prob);
}