
data {
  int<lower=1> N;
  vector[N] y;
  vector[N] time;
  real t_min;
  real t_max;
  real<lower=0> kappa;

  int<lower=1> N_iss;
  int<lower=1> N_vol;
  int<lower=1> N_mag;

  array[N] int<lower=1,upper=N_iss> iss;
  array[N] int<lower=1,upper=N_vol> vol;
  array[N] int<lower=1,upper=N_mag> mag;
}

parameters {
  vector[4] step_a;
  ordered[3] c_raw;

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
  vector[3] cut;
  vector[N] s1;
  vector[N] s2;
  vector[N] s3;

  for(j in 1:3) cut[j] = t_min + (t_max-t_min)*inv_logit(c_raw[j]);
  s1 = inv_logit(kappa*(time-cut[1]));
  s2 = inv_logit(kappa*(time-cut[2]));
  s3 = inv_logit(kappa*(time-cut[3]));
  mu = step_a[1]*(1-s1) + step_a[2]*(s1-s2) + step_a[3]*(s2-s3) + step_a[4]*s3 + z_iss[iss]*sd_iss + z_vol[vol]*sd_vol + z_mag[mag]*sd_mag;
}

model {
  step_a ~ normal(0,2);
  c_raw ~ normal(0,1.5);
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

