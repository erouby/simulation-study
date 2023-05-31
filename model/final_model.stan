data {
 int<lower = 1> n_ind_surv;                             // number of individuals for survival
 int<lower = 1> n_ind_repro;                            // number of individuals for reproduction
 int<lower = 1> n_ind_diff;                             // number of individuals added for reproduction
 int<lower = 1> n_ind_total;
 vector[n_ind_total] LOGTIME;                            // time to event (log scale) for survivorship
 vector<lower= 0.0>[n_ind_repro] TIME;                  // time event (normal scale) for reproduction
 int<lower = 1> n_year;                                 // number of years
 int<lower = 1> n_cov;                                  // number of covariates
 vector<lower = 1, upper = n_year>[n_year] x_trend;     // trend in the time series for both
 matrix[n_ind_total, n_cov] X;                           // covariate matrix for survivorship
 int<lower = 1, upper = n_year> YEAR[n_ind_total];       // Year index for survivorship data
 int<lower = 0, upper = 1> SURVIVAL[n_ind_total];        // Survivorship index (all are dead)
 int<lower = 0, upper = 2> CENSORING[n_ind_repro];      // Censoring index (depends on sexual maturity)
 vector<lower = 0>[2] prior_scale_for_intercept;        // depends on log TIME scale and TIME scale
 real<lower = 0> prior_scale_for_residual_sd;           // depends on log TIME scale and TIME scale
 vector<lower = 0>[2] prior_scale_for_slope;            // depends on log TIME scale and TIME scale
 real<lower = 0> prior_scale_for_sigma_frailty;         // depends on log TIME scale and TIME scale
 real<lower = 0> prior_scale_for_alpha;                 // depends on TIME scale 
 vector<lower = 0, upper = 1>[2] gamma_cov;                // parameter to add a covariate in the estimation or not
 vector<lower = 0, upper = 1>[2] gamma_random;             // parameter to add a random effect in the estimation or not
 vector<lower = 0, upper = 1>[2] gamma_trend;             // parameter to add a trend effect in the estimation or not
}

transformed data {
  vector[n_year] stdYEAR;
  for(i in 1:n_year) {
   stdYEAR[i] = (i - mean(x_trend)) / sd(x_trend);
  }
}

parameters {
 vector[2] unscaled_intercept;                // Intercept as real for both survivorship and reproduction
 vector[2] unscaled_slope[n_cov];             // Slope for covariate as real for both survivorship and reproduction
 real<lower = 0> unscaled_sigma2;             // Scale parameter for the frailty
 vector[2] u_year[n_year];                    // random effect parameter
 vector[2] trend;                             // yearly trend in hazard                           
 vector<lower = 0>[4] tau;                    // Parameter for beta2 distribution
 simplex[4] prop;                             // variance partitioning
 vector<lower = 0>[n_ind_surv] Z;             // Individual Frailty
 real log_alpha;                              // log of shape parameter for reproduction
 real<lower = -1, upper = 1> cor;
}

transformed parameters {
 vector[2] intercept;
 vector[2] sigma_year;
 real residual_sd;
 vector[2] R_sq;
 vector[n_ind_surv] linear_predictor_surv;
 vector[n_ind_surv] residuals_surv;
 vector[n_ind_repro] linear_predictor_repro;
 vector[n_ind_repro] residuals_repro;
 vector[2] slope[n_cov];
 real sigma_frailty;
 cov_matrix[2] Omega;
 real alpha;
 real beta;
// scale for frailty
 sigma_frailty = prior_scale_for_sigma_frailty * sqrt(prop[1] * unscaled_sigma2 / tau[1]);
 beta = 1 / sigma_frailty;
 // residual std. dev.
 residual_sd = prior_scale_for_residual_sd * sqrt(prop[2] * unscaled_sigma2 / tau[2]);
 // scale for year effect
 sigma_year[1] = prior_scale_for_residual_sd * sqrt(prop[3] * unscaled_sigma2 / tau[3]);
 sigma_year[2] = prior_scale_for_residual_sd * sqrt(prop[4] * unscaled_sigma2 / tau[4]);
 Omega[1, 1] = 1.0;
 Omega[2, 2] = 1.0;
 Omega[1, 2] = cor;
 Omega[2, 1] = cor;
 // slope
 for (n in 1:n_cov) {
 slope[n, 1] = unscaled_slope[n, 1] * prior_scale_for_slope[1];
 slope[n, 2] = unscaled_slope[n, 2] * prior_scale_for_slope[2];
}
 // intercept
 intercept[1] = unscaled_intercept[1] * prior_scale_for_intercept[1];
 intercept[2] = unscaled_intercept[2] * prior_scale_for_intercept[2];
 // shape param. (proportional to a precision)
 alpha = exp(log_alpha) * prior_scale_for_alpha;  // should be > 1 (hazard is increasing with time)
 // linear predictors
 for(j in 1:n_ind_repro) {
  // part without the covariates
  if (gamma_random[2] == 1) 
  linear_predictor_repro[j] = intercept[2] + sigma_year[2] * (gamma_trend[2] * trend[2] * stdYEAR[YEAR[j]] + u_year[YEAR[j], 2] * gamma_random[2]);
  else if (gamma_random[2] == 0)
  linear_predictor_repro[j] = intercept[2] + 1 * (gamma_trend[2] * trend[2] * stdYEAR[YEAR[j]] + u_year[YEAR[j], 2] * gamma_random[2]);
  for (n in 1:n_cov) {
   // add each covariate
   linear_predictor_repro[j] += gamma_cov[2] * X[j, n] * slope[n, 2];
   }
   // weibull regression
   linear_predictor_repro[j] = exp(-linear_predictor_repro[j]/alpha);
   // now linear_predictor is complete: compute residual
  residuals_repro[j] = TIME[j] - linear_predictor_repro[j];
 }
  for(i in 1:n_ind_surv){
  // part without the covariates
  if (gamma_random[1] == 1)
  linear_predictor_surv[i] = intercept[1] - (Z[i] * sigma_frailty) + sigma_year[1] * (gamma_trend[1] * trend[1] * stdYEAR[YEAR[i + n_ind_diff]] + u_year[YEAR[i + n_ind_diff], 1] * gamma_random[1]);
  else if (gamma_random[1] == 0)
  linear_predictor_surv[i] = intercept[1] - (Z[i] * sigma_frailty) + 1 * (gamma_trend[1] * trend[1] * stdYEAR[YEAR[i + n_ind_diff]] + u_year[YEAR[i + n_ind_diff], 1] * gamma_random[1]);
  for (n in 1:n_cov) {
   // add each covariate
   linear_predictor_surv[i] += gamma_cov[1] * X[i + n_ind_diff, n] * slope[n, 1];
  }
  // now linear_predictor is complete: compute residual
  residuals_surv[i] = LOGTIME[i + n_ind_diff] - linear_predictor_surv[i];
}
 R_sq[1] = 1 - variance(residuals_surv) / variance(LOGTIME);
 R_sq[2] = 1 - variance(residuals_repro) / variance(TIME);
}

model {
 // weakly informative priors
 unscaled_intercept[1] ~ normal(0.0, 1.0);
 unscaled_intercept[2] ~ normal(0.0, 1.0);
 unscaled_sigma2 ~ gamma(0.5, 1.0);
 for (n in 1:n_cov) {
  unscaled_slope[n] ~ multi_normal(rep_vector(0.0, 2), diag_matrix(rep_vector(1.0, 2)));
 }
 for(i in 1:n_year) {
  u_year[i] ~ multi_normal(rep_vector(0.0, 2), Omega);
 }
 trend[1] ~ normal(0.0, 1.0);
 trend[2] ~ normal(0.0, 1.0);
 tau ~ gamma(1.0, 1.0);
 Z ~ exponential(1.0);
 log_alpha ~ normal(0.0, 1.0);
 
 for (i in 1:n_ind_surv) {
  if(SURVIVAL[i + n_ind_diff] == 0) { 
    target += normal_lpdf(LOGTIME[i + n_ind_diff] | linear_predictor_surv[i], residual_sd);
  }
  else { // censored
    target += normal_lccdf(LOGTIME[i + n_ind_diff] | linear_predictor_surv[i], residual_sd);
  }
}
 for (j in 1:n_ind_repro) {
  if(CENSORING[j] == 0) { // observed
   target += weibull_lpdf(TIME[j]| alpha, linear_predictor_repro[j]);
  }
  else {
   if(CENSORING[j] == 1) { // left-censored
    target += weibull_lcdf(TIME[j]| alpha, linear_predictor_repro[j]);
   }
   else { // right-censored
    target += weibull_lccdf(TIME[j]| alpha, linear_predictor_repro[j]);
   }
  }
 }
}

generated quantities {
 vector[n_ind_surv] y_rep_surv;
 vector[n_ind_repro] y_rep_repro;
 vector[n_ind_surv] log_lik_surv;
 vector[n_ind_repro] log_lik_repro;
 real ppc;
 real TG;

for(j in 1:n_ind_repro) {
  y_rep_repro[j] = weibull_rng(alpha, linear_predictor_repro[j]);
  if(CENSORING[j] == 0) { // observed
    log_lik_repro[j] = weibull_lpdf(TIME[j]| alpha, linear_predictor_repro[j]);
  }
  else { // censored
   if(CENSORING[j] == 1) { // left-censored
     log_lik_repro[j] = weibull_lcdf(TIME[j]| alpha, linear_predictor_repro[j]);
   }
   else { // right-censored
     log_lik_repro[j] = weibull_lccdf(TIME[j]| alpha, linear_predictor_repro[j]);
   }
  }
 ppc = (max(y_rep_surv) > max(LOGTIME)) ? 1 : 0;
 TG = beta_rng(22.036479, 6.6720656);
 } 
 for(i in 1:n_ind_surv) {
  y_rep_surv[i] = normal_rng(linear_predictor_surv[i], residual_sd);
  if(SURVIVAL[i + n_ind_diff] == 0) { 
   log_lik_surv[i] = normal_lpdf(LOGTIME[i + n_ind_diff] | linear_predictor_surv[i], residual_sd);
  }
  else { // censored
   log_lik_surv[i] = normal_lccdf(LOGTIME[i + n_ind_diff] | linear_predictor_surv[i], residual_sd);
  }
 } 
}
