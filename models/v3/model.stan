data {
  int K; // number of stages
  int N; // number of sites
  int D; // number of predictor variables
  array[K, K, N] int y; // observed outcomes
  matrix[N, D] x; // predictor data
  array[K] matrix[D, K] beta_switch;
  array[K] vector[K] alpha_switch;
}

parameters {
  array[K] vector[K] alpha;
  array[K] matrix[D, K] beta;
}

transformed parameters {
  array[K] vector[K] alpha_switched;
  array[K] matrix[D, K] beta_switched;
  for(k in 1:K) {
    // alpha_switched[k] = alpha[k] + log(alpha_switch[k]); // log converts [0,1] to [-Inf,0]
    alpha_switched[k] = alpha[k] - 20 * (1.0 - alpha_switch[k]); // this gets alpha very small... my attempts to make it -Inf caused sampling to fail
    beta_switched[k] = beta[k] .* beta_switch[k]; // zero-out irrelevant betas
  }
}

model {
  
  for (k in 1:K) {
    alpha[k] ~ normal(0, 10);
    to_vector(beta[k]) ~ normal(0, 1);
    
    matrix[N, K] a_x_b = x * beta_switched[k];
    for(kk in 1:K) {
      a_x_b[, kk] = a_x_b[, kk] + alpha_switched[k][kk];
    }
    
    for (n in 1:N) {
      y[k, , n] ~ multinomial(softmax(a_x_b[n]'));
    }
  }
  
}
