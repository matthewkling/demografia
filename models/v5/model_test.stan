functions {
  vector project_simplex(matrix x, int k, int t) {
    row_vector[rows(x)] y = rep_row_vector(0, rows(x));
    y[k] = 1;
    for(j in 1:t) {
      y = y * x;
    }
    return to_vector(y);
  }
}

data {
  int K; // number of stages
  row_vector[K] n0; // collective frequencies at time 0
  vector[K] nt; // collective frequencies at time t
  int nk; // index of earliest stage at which frequency is observed
  array[2, K+1] int io; // individual outcomes
  array[2] int iok; // indices of stages included in io data
  int t; // number of timesteps
  real F; // fecundity 
  int Fk; // index of reproductive stage
  array[K] vector[K+1] gamma;
  matrix[K, K+1] zero_alpha;
}

parameters {
  array[K] simplex[K+1] alpha;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[K, K+1] A; // projection matrix

  for (k in 1:K) {
    A[k, ] = to_row_vector(alpha[k]) .* zero_alpha[k, ];
    A[k, ] = A[k, ] / sum(A[k, ]);
  }
  A[Fk, 1] = F;
}


model {
  
  // priors
  for (k in 1:K) {
    alpha[k] ~ dirichlet(gamma[k]);
  }
  
  // collective distribution changes
  row_vector[K] P = n0;
  for(T in 1:t) {
    P = P * A[, 1:K];
  }
  log(nt[nk:K]) ~ normal(log(P[nk:K]), sigma);
  
  // individual outcomes
  matrix[K+1, K+1] S = append_row(A, rep_row_vector(0, K+1));
  S[K+1, K+1] = 1;
  S[Fk, 1] = 0;
  for(k in 1:size(iok)) {
    vector[K+1] p = project_simplex(S, iok[k], t);
    p = p + 1e-12;
    p = p / sum(p);
    io[k] ~ multinomial(p);
  }

}
