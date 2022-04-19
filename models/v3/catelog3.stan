data {
  int K; // number of stages
  int N; // number of sites
  int D; // number of predictor variables
  array[K, N, K] int y; // observed outcomes
  matrix[N, D] x; // predictor data
}
parameters {
  array[K] matrix[D, K] beta;
}
model {
  
  for (k in 1:K) {
    matrix[N, K] x_beta = x * beta[k];
    to_vector(beta[k]) ~ normal(0, 5);
    
    for (n in 1:N) {
      y[k, n] ~ categorical_logit(x_beta[n]');
    }
  }
  
}