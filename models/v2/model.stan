data {
  int K; // number of stages
  int t; // nubmer of timesteps
  array[K, K+1] int y; // observed transitions
  array[K] vector[K+1] alpha; // dirichlet priors
}

parameters {
  simplex[K+1] delta[K];
}

transformed parameters {
  simplex[K+1] theta[K]; // annualized version of delta. could implement in postprocessing instead of model.
  for(i in 1:K){
    theta[i] = delta[i];
    real p = pow(theta[i][i], 1.0 / t);
    theta[i][i] = 0;
    theta[i] = (1.0 - p) * theta[i] / sum(theta[i]);
    theta[i][i] = p;
  }
}

model {
  for(i in 1:K){
    y[i] ~ multinomial(delta[i]); // comment out to sample from prior
    delta[i] ~ dirichlet(alpha[i]);
  }
}
