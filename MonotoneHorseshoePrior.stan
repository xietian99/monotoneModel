// The input data is a vector 'y' of length 'N'.
data { // fixed 
  // Define N, L, D(k)  X  Y 
  int<lower=0> N;
  int<lower=1> L; // number of dose-level
  int<lower = 0, upper = 1> Y[N]; // outcome binary variable
  int<lower = 1, upper = L> D[N]; // dose level 
  vector[N] X;  // for one co-variate, for other matrix[N, M] X 
  
  // Define the local parameter needed 
  real<lower = 0> local_dof_stan; // ??degree of freedom of pi(lambda), = 1
  real<lower = 0> global_dof_stan; // ??degree of freedom  of pi(tau), = 1
  real<lower = 0> c_sq; // c^2  TIP: may vector[L] c for inequal C interval. 
  
  real<lower = 0> slab_precision_stan; // assumed to be 1 in the manuscript
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'alpha' and 'beta'.
// And hyperparameter of alpha/beta 

parameters {
  // Main Parameter Estimated
  vector<lower = 0.0>[L + 1] alpha_base; //?? why alpha_raw? 
  real beta_base; // if multi-dim covariae vector<M> beta
  
  // Hyper-paramter
  // tau: half-cauchy decomposed into chi-square (N^2) and inverse-gamma
  real<lower = 0.0> tau_base_sq;
  real<lower = 0.0> tau_scale_sq;
  // lambda: half-cauchy decomposed into chi-square (N^2) and inverse-gamma
  vector<lower = 0.0>[L+1] lambda_base_sq;
  vector<lower = 0.0>[L+1] lambda_scale_sq; 


}

transformed parameters {
  // Interested
  vector<lower = 0.0>[L+1] alpha; 
  real<lower = 0.0> beta;

  vector<lower = 0.0, upper = 1.0>[L+1] theta; // alpha = alpha_base * theta
  real<lower = 0.0> sigma_sq; // beta = beta * sqrt(sigma)
  real<lower = 0.0> tau_sq;
  vector<lower = 0.0>[L+1] lambda_sq;
  
  // Epsilon
  vector<lower = 0.0, upper = 1>[L] epsi;
  vector<lower = 0.0, upper = 1.0>[L+1] normalized_alpha;
  

  // transformation
  beta = beta_base * sqrt(sigma_sq);
  tau_sq = tau_base_sq * tau_scale_sq;
  lambda_sq = lambda_base_sq .* lambda_scale_sq;
  for(i in 1:(L+1)) {
    theta[i] = 1.0 / sqrt(slab_precision_stan + (1.0 / (c_sq * tau_sq * lambda_sq[i])));
  }
  alpha = (theta .* alpha_base);
  normalized_alpha = alpha / sum(alpha);
  epsi[1] = normalized_alpha[1];
  if(L > 1) {
    for(i in 2:L) {
      epsi[i] = epsi[i-1] + normalized_alpha[i];
    }
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  alpha_base ~ normal(0.0, 1.0);
  beta_base ~ normal(0.0, 1.0);
  tau_base_sq ~ chi_square(1.0);
  tau_scale_sq ~ inv_gamma(global_dof_stan/2.0, global_dof_stan/2.0);
  lambda_base_sq ~ chi_square(1.0);
  lambda_scale_sq ~ inv_gamma(local_dof_stan/2.0, local_dof_stan/2.0);
  real p;
  for (i in 1:N) {
    p = inv_logit(X[i] * beta + logit(epsi[D[i]]));
    Y[i] ~ bernoulli(p);
  }
}


