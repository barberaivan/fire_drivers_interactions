/*
  Conditional logistic model to classify core vs. edge pixels.
  This is like a categorical model with K = 2, where the parameters are
  fixed among categories but the predictors values vary.
*/

data {
  // Ns
  int<lower=0> N; // observations
  int<lower=0> N_pairs; // pairs of observations (N/2)
  int<lower=0> K; // observation-level predictors dimension

  // ids
  int core_ids[N_pairs]; // indices corresponding to core in X
  int edge_ids[N_pairs]; // indices corresponding to edge in X

  // predictors
  matrix[N, K] X;

  // priors
  real prior_b_sd;     // sd for all coefficients
}

parameters {
  vector[K] b_raw;
}

transformed parameters {
  vector[K] b = b_raw * prior_b_sd;
  vector[N] linpred = X * b;
  matrix[N_pairs, 2] linpred_mat;
  vector[N_pairs] fitted_prob;

  linpred_mat = append_col(linpred[core_ids], linpred[edge_ids]);
  for(i in 1:N_pairs) {
    fitted_prob[i] = softmax(linpred_mat[i, ]')[1]; // core probability
  }
}

model {
  // priors
  b_raw ~ std_normal();

  // likelihood
  target += log(fitted_prob);
}