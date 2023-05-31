lapply(c("rstan", "loo"), library, character.only = TRUE)

options(mc.cores = 4)

load("data/Dolphin.RData")
writeLines("Data are loaded. Now fit the model and select the better one.")

### Stan parameters and model fit
load("data/stanmodel.RData")

## MCMC settings
n_chains <- 4
n_thin <- 1
n_warm <- 1000 * n_thin
n_iter <- n_warm + 2500 * n_thin

for (i in model_index) {
  stanfit = sampling(object = test,
                          data = list(n_ind_surv = n_ind_surv,
                                      n_ind_repro = n_ind_repro,
                                      n_ind_diff = n_ind_diff,
                                      n_ind_total = n_ind_total,
                                      LOGTIME = log(age_data),
                                      TIME =  SR_age_data,
                                      SURVIVAL = rep(0, n_ind_total),
                                      CENSORING = 2 - SR_data,
                                      YEAR = YEAR,
                                      n_year = n,
                                      x_trend = x_trend,
                                      n_cov = n_cov,
                                      X = X,
                                      prior_scale_for_intercept = rep(1.5, 2),
                                      prior_scale_for_sigma_frailty = log(2),
                                      prior_scale_for_residual_sd = 0.1,
                                      prior_scale_for_slope = rep(1, 2),
                                      prior_scale_for_alpha = 1.0,
                                      gamma_cov = possibility[[gammas$gamma_cov[i]]], # First is for survival, second is for reproduction
                                      gamma_random = possibility[[gammas$gamma_random[i]]], # First is for survival, second is for reproduction
                                      gamma_trend = possibility[[gammas$gamma_trend[i]]] # First is for survival, second is for reproduction
                          ),
                          pars = c("intercept", "residual_sd", "beta", "sigma_frailty", "slope", "alpha", "trend", "cor", "sigma_year", "u_year", "R_sq", "ppc", "TG", "log_lik_surv", "log_lik_repro"),
                          chains = n_chains, iter = n_iter, warmup = n_warm, thin = n_thin,
                          control = list(adapt_delta = 0.9, max_treedepth = 15)
                          )
  
  save(list = "stanfit", file = paste("output/model_", i, ".RData" , sep = ""))
  
  waic_m <- loo::waic(loo::extract_log_lik(stanfit = stanfit, 
                                              parameter_name = c("log_lik_surv", "log_lik_repro")
                                              )
                         )[3:8]
  
  save(list = "waic_m", file = paste("output/waic_", i, ".RData" , sep = ""))
}

