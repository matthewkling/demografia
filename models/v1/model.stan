data {
  int K; // number of stages
  array[K, K+1] int y; // observed transitions
  array[K] vector[K+1] alpha; // dirichlet priors
}

parameters {
  simplex[K+1] theta[K];
}

model {
  for(i in 1:K){
    y[i] ~ multinomial(theta[i]); // comment out to sample from prior
    theta[i] ~ dirichlet(alpha[i]);
  }
}
