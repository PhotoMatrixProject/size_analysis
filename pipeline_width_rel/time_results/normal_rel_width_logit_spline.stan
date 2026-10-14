
data {
  int<lower=1> N;
  vector[N] y;

  int<lower=1> K;
  matrix[N,K] B;

  int<lower=1> N_iss;
  int<lower=1> N_vol;
  int<lower=1> N_mag;

  array[N] int<lower=1,upper=N_iss> iss;
  array[N] int<lower=1,upper=N_vol> vol;
  array[N] int<lower=1,upper=N_mag> mag;
}

parameters {
  real a;
  vector[K] w;

  vector[N_iss] z_iss;
  vector[N_vol] z_vol;
  vector[N_mag] z_mag;

  real<lower=0> sd_iss;
  real<lower=0> sd_vol;
  real<lower=0> sd_mag;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu;
  mu = a + B*w + z_iss[iss]*sd_iss + z_vol[vol]*sd_vol + z_mag[mag]*sd_mag;
}

model {
  a ~ normal(0,2);
  w ~ normal(0,1);
  z_iss ~ normal(0,1);
  z_vol ~ normal(0,1);
  z_mag ~ normal(0,1);
  sd_iss ~ exponential(1);
  sd_vol ~ exponential(1);
  sd_mag ~ exponential(1);
  sigma ~ exponential(1);
  y ~ normal(mu,sigma);
}

generated quantities {
  vector[N] log_lik;
  for(n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu[n],sigma);
}

