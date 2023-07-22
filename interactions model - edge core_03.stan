/*
  Conditional logistic model to classify core vs. edge pixels.
  This is like a categorical model with K = 2, where the parameters are
  fixed among categories but the predictors values vary.
*/

data {
  // Ns
  int<lower=0> N; // observations
  int<lower=0> N_pairs; // pairs of observations (N/2)
  int<lower=0> NF; // fires (groups)
  int<lower=0> K; // observation-level predictors dimension
  int<lower=0> J; // group-level predictors dimension

  // ids
  int starts[NF];  // positions in X where each fire id starts
  int ends[NF];    // or ends
  int core_ids[N_pairs]; // indices corresponding to core in X
  int edge_ids[N_pairs]; // indices corresponding to edge in X

  // predictors
  matrix[N, K] X;
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
}

transformed parameters {
  // intercepts and slopes for the hiperparameters
  matrix[K, J] g;
  // coefficients matrix
  matrix[K, NF] b; // b_mean

  // linear predictor
  vector[N] linpred;
  // fitted prob
  vector[N_pairs] fitted_prob;

  g = g_raw * prior_b_sd;
  b = g * Zt;

  // compute linear predictor
  for(f in 1:NF) {
    // make a matrix multiplication but taking only each fire's data and
    // coefficients
    linpred[starts[f] : ends[f]] = X[starts[f] : ends[f], ] * b[, f];
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
  to_vector(g_raw) ~ std_normal();

  // likelihood
  target += log(fitted_prob);
}