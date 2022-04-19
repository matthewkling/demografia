functions {
  vector project_simplex(matrix x, int t, int k, int K) {
    row_vector[K] y = rep_row_vector(0, K);
    y[k] = 1;
    for(j in 1:t) {
      y = y * x;
    }
    return to_vector(y);
  }
}

data {
  int K; // number of stages
  int N; // number of sites
  int D; // number of predictor variables
  array[K, K, N] int y; // observed outcomes
  matrix[N, D] x; // predictor data
  int t[N]; // timesteps since prior survey
  array[D] matrix[K, K] beta_switch;
  matrix[K, K] alpha_switch;
}

parameters {
  matrix<upper=0>[K, K] alpha;
  array[D] matrix[K, K] beta;
}

transformed parameters {
  matrix[K, K] alpha_switched;
  array[D] matrix[K, K] beta_switched;
  for (k in 1:K) {
    alpha_switched[k] = alpha[k] - 20 * (1.0 - alpha_switch[k]); // this gets alpha very small... my attempts to make it -Inf caused sampling to fail
  }
  for (d in 1:D) {
    beta_switched[d] = beta[d] .* beta_switch[d]; // zero-out irrelevant betas
  }
}

model {
  
  // priors
  for (k in 1:K) {
    logit(exp(alpha[k])) ~ normal(0, 5);
  }
  for (d in 1:D) {
    to_vector(beta[d]) ~ normal(0, 1);
  }
  
  // likelihood
  for (n in 1:N) {
    matrix[K, K] p = alpha_switched;
    for (d in 1:D) {
      p = p + x[n, d] * beta_switched[d];
    }
    for (k in 1:K) { // has to be separate loop from the one below
      p[k] = to_row_vector(softmax(to_vector(p[k])));
    }
    for (k in 1:K) {
      y[k, , n] ~ multinomial(project_simplex(p, t[n], k, K));
    }
  }
  
}
