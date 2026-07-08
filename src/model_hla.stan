data {
  int<lower=1> K;                 // number of rows (cases)
  int<lower=0> y[K];              // hits
  int<lower=0> n[K];              // tries
  int<lower=0> Y[K];              // control hits
  int<lower=0> N[K];              // control tries
}

parameters {
  vector<lower=0,upper=1>[K] p;   // probabilities for cases
  vector<lower=0,upper=1>[K] P;   // probabilities for controls
}

model {

  // weak priors (uniform Beta(1,1))
  p ~ beta(1,1);
  P ~ beta(1,1);

  // likelihood
  y ~ binomial(n, p);
  Y ~ binomial(N, P);
}

generated quantities {

  vector[K] diff;
  vector[K] log_fold_change;

  for (k in 1:K) {
    diff[k] = p[k] - P[k];
    log_fold_change[k] = log(p[k]) - log(P[k]);
  }

}