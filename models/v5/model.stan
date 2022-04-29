functions {
  vector project_simplex(matrix x, int k, int t) {
    row_vector[rows(x)] y = rep_row_vector(0, rows(x));
    y[k] = 1;
    for(j in 1:t) {
      y = y * x;
    }
    return to_vector(y);
  }
  
  array[] int index_ones(array[] int x) {// return a vector of the indices where boolean verctor x == 1
    array[sum(x)] int index;
    int j = 1;
    for (i in 1:size(x)) {
      if (x[i] == 1) {
        index[j] = i;
        j = j + 1;
      }
    }
    return(index);
  }
}

data {
  int K; // number of stages
  int B; // number of blocks
  array[B] row_vector[K] n0; // collective frequencies at time 0
  array[B] vector[K] nt; // collective frequencies at time t
  int nk; // index of earliest stage at which frequency is observed
  array[B, 2, K+1] int io; // individual outcomes
  array[2] int iok; // indices of stages included in io data
  int t; // number of timesteps
  real F; // fecundity 
  int Fk; // index of reproductive stage
  array[K, K+1] int alpha_i;
  int n_alpha;
  vector[n_alpha] alpha_mu; // prior means for alpha
  vector<lower=0>[n_alpha] alpha_sigma; // prior variance for alpha
}

parameters {
  vector[n_alpha] alpha; // on LOG scale
  real<lower=0> sigma;
}

transformed parameters {
  matrix[K, K+1] A; // projection matrix
  
  {
    int start = 0;
    for (k in 1:K) {
      A[k,] = rep_row_vector(0, K+1);
      int na = sum(alpha_i[k]); // number of alphas for this life stage
      array[na] int ai;
      for (i in 1:na) ai[i] = i + start;
      A[k, index_ones(alpha_i[k])] = to_row_vector(softmax(alpha[ai]));
      start = start + na;
    }
  }
  
  A[Fk, 1] = F;
}


model {
  
  // priors
  alpha ~ normal(alpha_mu, alpha_sigma);
  
  for(b in 1:B) {
    
    // local projection matrix
    
    
    // collective distribution changes
    row_vector[K] P = n0[b];
    for(T in 1:t) {
      P = P * A[, 1:K];
    }
    log(nt[b, nk:K]) ~ normal(log(P[nk:K]), sigma);
    
    // individual outcomes
    matrix[K+1, K+1] S = append_row(A, rep_row_vector(0, K+1));
    S[K+1, K+1] = 1;
    S[Fk, 1] = 0;
    for(k in 1:size(iok)) {
      vector[K+1] p = project_simplex(S, iok[k], t);
      p = p + 1e-12;
      p = p / sum(p);
      io[b, k] ~ multinomial(p);
    }
    
  }
  
  
  
}
