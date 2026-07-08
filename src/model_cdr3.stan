
data {
  int<lower=0> N;
  vector[N] y_pt;
  vector[N] y_sln;
}

transformed data {
  vector [N] y_delta = y_pt - y_sln;
}

parameters {
  real alpha;
  real<lower=0> sigma;
}

model {
  alpha ~ normal(0,1);
  sigma ~ normal(0,1);
  y_delta ~ normal(alpha, sigma);
}

generated quantities {
  real yhat_PT_vs_SLN;
  yhat_PT_vs_SLN = normal_rng(alpha, sigma);
}

