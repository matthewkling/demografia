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
  int yk0; // index of first stage for y data
  int yK0; // number of t0 stages for y data
  int yK1; // number of tt stages for y data
  array[yK0, yK1, N] int y; // observed individual-level outcomes
  array[2, K, N] real m; // real if Gauss; observed class counts at 2 time points
  // array[2, K, N] int m; // int if multiniomial; observed class counts at 2 time points
  int mk0; // indices of first stage for which m is observed
  int mk1; // indices of last stage for which m is observed
  matrix[N, D] x; // predictor data
  int t[N]; // timesteps since prior survey
  real F; // fecundity of each stage
  int Fk; // index of reproductive stage
  array[D] matrix[K, K] zero_beta;
  matrix[K, K] zero_alpha;
}

parameters {
  matrix<upper=0>[K, K] alpha;
  array[D] matrix[K, K] beta;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[K, K] zeroed_alpha;
  array[D] matrix[K, K] zeroed_beta;
  for (k in 1:K) {
    zeroed_alpha[k] = alpha[k] - 20 * (1.0 - zero_alpha[k]); // this gets alpha very small... my attempts to make it -Inf caused sampling to fail
  }
  for (d in 1:D) {
    zeroed_beta[d] = beta[d] .* zero_beta[d]; // zero-out irrelevant betas
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
    
    // apply predictor effects
    matrix[K, K] A = zeroed_alpha;
    for (d in 1:D) {
      A = A + x[n, d] * zeroed_beta[d];
    }
    
    // normalize to simplexes (note: can't be combined w loop below)
    for (k in 1:K) {
      A[k] = to_row_vector(softmax(to_vector(A[k])));
    }
    
    // individual-based transitions
    matrix[yK1, yK1] Ay = A[yk0:K, yk0:K];
    for (k in 1:yK0) {
      y[k, , n] ~ multinomial(softmax(log(project_simplex(Ay, t[n], k, K - yK1 + 1))));
    }
    
    // collective distribution changes
    row_vector[K] P = to_row_vector(m[1, , n]);
    matrix[K, K] AF = A;
    AF[Fk, 1] = F;
    for(T in 1:t[n]) {
        P = P * AF;
    }
    
    
    // multinomial approach -- i.e. change in stage distribution
    // P = P[mk0:mk1];
    // m[2, mk0:mk1, n] ~ multinomial(to_vector(P / sum(P)));
    
    // another approach -- Gaussian error on log counts
    log(m[2, mk0:mk1, n]) ~ normal(log(P[mk0:mk1]), sigma);
    
  }
  
}
