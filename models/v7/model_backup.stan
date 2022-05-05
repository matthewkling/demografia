functions {
  
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
  int T; // number of timesteps
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
  
  int n_gamma;
  array[K, K, K+1] int gamma_i; // boolean flags nonzero gammas
  real<lower=0> gamma_scale; // prior variance for gammas
}

parameters {
  vector[n_alpha] alpha; // on LOG scale
  vector[n_beta] beta;
  vector<lower=0>[n_gamma] gamma; // positive paramater, later made negative during transformation
  real<lower=0> sigma;
}

transformed parameters {
  matrix[K, K+1] A; // projection matrix intercepts
  array[D] matrix[K, K+1] B; // environmental effects (betas)
  array[K] matrix[K, K+1] G; // density effects (gammas)
  
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
    
    // G
    i = 1;
    for(j in 1:K) {
      for(k in 1:K) {
        for(l in 1:(K+1)) {
          G[j,k,l] = 0;
          if(gamma_i[j,k,l] == 1) {
            G[j,k,l] = 0 - gamma[i];
            i = i + 1;
          }
        }
      }
    }
    

  }
  
}

model {
  
  // local varibles
  matrix[K, K+1] AB;
  matrix[K, K+1] ABG; // projection matrix incorporating density
  matrix[K, K+1] ABG_c;
  matrix[K+1, K+1] ABG_i;
  row_vector[K] P; // total population in each stage
  matrix[size(iok), K+1] M; // individual outcomes
  
  // priors
  alpha ~ normal(alpha_mu, alpha_sigma);
  beta ~ normal(0, beta_sigma);
  gamma ~ exponential(gamma_scale); // exponential takes inverse scale param in stan (opposite of R)
  
  
  for(n in 1:N) {
    
    // environmental effects (beta)
    AB = log(A);
    for (d in 1:D) AB = AB + x[n, d] * B[d];
    
    // local variables
    P = to_row_vector(n0[n]); // total population in each stage
    M = rep_matrix(0, size(iok), K+1);
    for(i in 1:size(iok)) M[i, iok[i]] = 1;
    
    // multiyear population projection
    for(t in 1:T) {
      
      // density effects (gamma)
      ABG = AB;
      for(k in 1:K) ABG = ABG + P[k] * G[k];
      
      // convert to probabilities
      for(k in 1:K) ABG[k, ] = to_row_vector(softmax(to_vector(ABG[k, ])));
      
      // collective projection (incorporating fecundity)
      ABG_c = ABG;
      ABG_c[Fk, 1] = F;
      P = P * ABG_c[, 1:K];
      
      // individual outcomes (symmetric matrix without fecundity)
      ABG_i = append_row(ABG, rep_row_vector(0, K+1));
      ABG_i[K+1, K+1] = 1;
      for(i in 1:size(iok)) M[i, ] = M[i, ] * ABG_i;
    }
    
    // likelihood
    sqrt(nt[n, nk:K]) ~ normal(sqrt(P[nk:K]), sigma); 
    for(i in 1:size(iok)) io[n, i] ~ multinomial(to_vector((M[i, ] + 1e-20) / sum(M[i, ] + 1e-20)));
  }
  
}
