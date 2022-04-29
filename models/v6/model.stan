functions {
  vector project_simplex(matrix x, int k, int t) {
    row_vector[rows(x)] y = rep_row_vector(0, rows(x));
    y[k] = 1;
    for(j in 1:t) {
      y = y * x;
    }
    return to_vector(y);
  }
  
  // return a vector of the indices where boolean verctor x == 1
  array[] int index_ones(array[] int x) {
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
  int N; // number of blocks
  int D; // number of predictors
  array[N] vector[K] n0; // collective frequencies at time 0
  array[N] vector[K] nt; // collective frequencies at time t
  int nk; // index of earliest stage at which frequency is observed
  array[N, 2, K+1] int io; // individual outcomes
  array[2] int iok; // indices of stages included in io data
  int t; // number of timesteps
  real F; // fecundity 
  int Fk; // index of reproductive stage
  matrix[N, D] x; // predictor data
  array[K, K+1] int alpha_i; // boolean flags for nonzero alphas
  int n_alpha;
  vector[n_alpha] alpha_mu; // prior means for alpha
  vector<lower=0>[n_alpha] alpha_sigma; // prior variances for alphas
  int n_beta;
  array[D, K, K+1] int beta_i; // boolean flags nonzero betas
  real<lower=0> beta_sigma; // prior variance for betas
}

parameters {
  vector[n_alpha] alpha; // on LOG scale
  vector[n_beta] beta;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[K, K+1] A; // projection matrix intercepts
  array[D] matrix[K, K+1] B; // projection matrix slopes
  
  {
    // A
    int start = 0;
    for (k in 1:K) {
      A[k, ] = rep_row_vector(0, K+1);
      int na = sum(alpha_i[k]); // number of alphas for this life stage
      array[na] int ai;
      for (i in 1:na) ai[i] = i + start;
      A[k, index_ones(alpha_i[k])] = to_row_vector(softmax(alpha[ai]));
      start = start + na;
    }
    
    // int i = 1;
    // for(j in 1:K) {
    //   for(k in 1:(K+1)) {
    //     A[j,k] = log(0);
    //     if(alpha_i[j,k] == 1) {
    //       A[j,k] = alpha[i];
    //       i = i + 1;
    //     }
    //   }
    //   A[j,] = to_row_vector(softmax(to_vector(A[j,])));
    // }
    
    // B
    int i = 1;
    for(j in 1:D) {
      for(k in 1:K) {
        for(l in 1:(K+1)) {
          B[j,k,l] = 0;
          if(beta_i[j,k,l] == 1) {
            B[j,k,l] = beta[i];
            i = i + 1;
          }
        }
      }
    }
    
  }
  
}

model {
  // priors
  alpha ~ normal(alpha_mu, alpha_sigma);
  beta ~ normal(0, beta_sigma);
  
  
  for(n in 1:N) {
    
    // assemble local projection matrix, incorporating betas and fecundity
    matrix[K, K+1] L = log(A);
    for (d in 1:D) {
      L = L + x[n, d] * B[d];
    }
    for(k in 1:K) {
        L[k, ] = to_row_vector(softmax(to_vector(L[k, ])));
    }
    L[Fk, 1] = F;
    
    // symmetric version, without fecundity
    matrix[K+1, K+1] S = append_row(L, rep_row_vector(0, K+1));
    S[K+1, K+1] = 1;
    S[Fk, 1] = 0;
    
    // collective distribution changes
    row_vector[K] P = to_row_vector(n0[n]);
    for(T in 1:t) {
      P = P * L[, 1:K];
    }
    sqrt(nt[n, nk:K]) ~ normal(sqrt(P[nk:K]), sigma); // log won't ultimately work, as it breaks with zeros 
    
    // individual outcomes
    vector[K+1] p;
    for(k in 1:size(iok)) {
      p = project_simplex(S, iok[k], t);
      p = p + 1e-20; // because probs can't be 0 for multinomial
      p = p / sum(p);
      io[n, k] ~ multinomial(p);
    }
    
    
  }
  
  
  
}
