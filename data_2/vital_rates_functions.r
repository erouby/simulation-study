### Survival function
get_surv <- function(stanfit, logt = NULL, alpha = 0.20, y_sim = NULL, n_year = NULL, cov = NULL, what = "survival") {
  
  ### Useful function 
  surv_stan <- function(stanfit, logt = NULL, alpha = 0.20, y_sim = NULL, y = NULL, what = "survival") {
    ### useful functions
    get_ci <- function(x) {
      c(coda::HPDinterval(coda::as.mcmc(x), prob = 1 - alpha)[1], 
        mean(x, na.rm = TRUE), 
        coda::HPDinterval(coda::as.mcmc(x), prob = 1 - alpha)[2]
      )
    }
    
    # sanity checks
    if(is.null(logt)) {
      if(is.null(y_sim)) {
        stop("Must provide either a vector of logt or a dataset y_sim")
      }
      else{
        logt <- log(survival::survfit(survival::Surv(Value, Status) ~ 1, 
                                      data = data.frame(Value = y_sim,
                                                        Status = rep(1, length(y_sim))
                                      ))$time
        )
      }
    }
    
    if (covariable == "cov1_0"){
      
      ### extract posteriors
      if(any(stanfit@model_pars == "beta")) {
        beta <- drop(rstan::extract(stanfit, 'beta')$beta)
        mu <- drop(rstan::extract(stanfit, 'intercept')$intercept[,1]) + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,1])) * (y - mean(x_trend)) / sd(x_trend))
        sd <- drop(rstan::extract(stanfit, 'residual_sd')$residual_sd)
        alpha <- rep(0, length(mu))
        sigma <- rep(0, length(mu))
        
        # survivorship fct
        surv <- function(logt, mu, sd, beta, sigma, alpha) {
          (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * (1 - pnorm(beta * sd  + (logt-mu)/sd))
        }
        
        # annual survival rate fct
        ann_surv <- function(logt, mu, sd, beta, sigma, alpha) {
          (1 - pnorm(logt[i+1], mu[i+1], sd[i+1])) - exp(beta[i+1] * (logt[i+1] - mu[i+1]) + 0.5 * (beta[i+1] * sd[i+1])^2) * (1 - pnorm(beta[i+1] * sd  + (logt[i+1]-mu[i+1])/sd[i+1])) / (1 - pnorm(logt[i], mu[i], sd[i])) - exp(beta[i] * (logt[i] - mu[i]) + 0.5 * (beta[i] * sd[i])^2) * (1 - pnorm(beta[i] * sd[i]  + (logt[i]-mu[i])/sd[i]))
        }
        
        # hazard fct
        hazard <- function(logt, mu, sd, beta, sigma, alpha) {
          # pdf of a Reed distribution
          reed_pdf <- function(logt, mu, sd, beta) {
            beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
          }
          # survival fct of a Reed distribution
          reed_phi <- function(logt, mu, sd, beta) {
            (1 - pnorm(logt, mu, sd)) - 
              exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
              (1 - pnorm(beta * sd + (logt - mu) / sd))
          }
          # hazard of a reed distribution
          reed_hz <- function(logt, mu, sd, beta) {
            reed_pdf(logt, mu, sd, beta) / reed_phi(logt, mu, sd, beta)
          }
          return(reed_hz(logt, mu, sd, beta))
        }
      }
      
      ### survival
      if(what %in% c("survival", "hazard", "ann_surv")) {
        if(what == "survival") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 surv(logt = logt,
                                      mu = mu[row_i],
                                      sd = sd[row_i],
                                      beta = beta[row_i],
                                      alpha = alpha[row_i],
                                      sigma = sigma[row_i]
                                 )
                               }
          )
        }
        if(what == "ann_surv") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 ann_surv(logt = logt,
                                          mu = mu[row_i],
                                          sd = sd[row_i],
                                          beta = beta[row_i],
                                          alpha = alpha[row_i],
                                          sigma = sigma[row_i]
                                 )
                               }
          )
        }
        if(what == "hazard") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 hazard(logt = logt,
                                        mu = mu[row_i],
                                        sd = sd[row_i],
                                        beta = beta[row_i],
                                        alpha = alpha[row_i],
                                        sigma = sigma[row_i]
                                 )
                               }
          )
        }
      }
      else {
        stop("Must choose either 'survival' or 'hazard' for arg. 'what'")
      }
      
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Stan"
      output$year <- y
      return(output)
    }
    
    if (covariable == "cov2_0"){
      
      ### extract posteriors
      if(any(stanfit@model_pars == "beta")) {
        beta <- drop(rstan::extract(stanfit, 'beta')$beta)
        mu <- drop(rstan::extract(stanfit, 'intercept')$intercept[,1]) + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,1])) * (y - mean(x_trend)) / sd(x_trend))
        sd <- drop(rstan::extract(stanfit, 'residual_sd')$residual_sd)
        alpha <- rep(0, length(mu))
        sigma <- rep(0, length(mu))
        
        # survivorship fct
        surv <- function(logt, mu, sd, beta, sigma, alpha) {
          (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * (1 - pnorm(beta * sd  + (logt-mu)/sd))
        }
        # hazard fct
        hazard <- function(logt, mu, sd, beta, sigma, alpha) {
          # pdf of a Reed distribution
          reed_pdf <- function(logt, mu, sd, beta) {
            beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
          }
          # survival fct of a Reed distribution
          reed_phi <- function(logt, mu, sd, beta) {
            (1 - pnorm(logt, mu, sd)) - 
              exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
              (1 - pnorm(beta * sd + (logt - mu) / sd))
          }
          # hazard of a reed distribution
          reed_hz <- function(logt, mu, sd, beta) {
            reed_pdf(logt, mu, sd, beta) / reed_phi(logt, mu, sd, beta)
          }
          return(reed_hz(logt, mu, sd, beta))
        }
      }
      
      ### survival
      if(what %in% c("survival", "hazard")) {
        if(what == "survival") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 surv(logt = logt,
                                      mu = mu[row_i],
                                      sd = sd[row_i],
                                      beta = beta[row_i],
                                      alpha = alpha[row_i],
                                      sigma = sigma[row_i]
                                 )
                               }
          )
        }
        if(what == "hazard") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 hazard(logt = logt,
                                        mu = mu[row_i],
                                        sd = sd[row_i],
                                        beta = beta[row_i],
                                        alpha = alpha[row_i],
                                        sigma = sigma[row_i]
                                 )
                               }
          )
        }
      }
      else {
        stop("Must choose either 'survival' or 'hazard' for arg. 'what'")
      }
      
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Stan"
      output$year <- y
      return(output)
    }
    
    if (covariable == "cov1_1"){
      
      ### extract posteriors
      if(any(stanfit@model_pars == "beta")) {
        beta <- drop(rstan::extract(stanfit, 'beta')$beta)
        mu <- drop(rstan::extract(stanfit, 'intercept')$intercept[,1] + rstan::extract(stanfit, 'slope')$slope[,1,1]) + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,1])) * (y - mean(x_trend)) / sd(x_trend))
        sd <- drop(rstan::extract(stanfit, 'residual_sd')$residual_sd)
        alpha <- rep(0, length(mu))
        sigma <- rep(0, length(mu))
        
        # survivorship fct
        surv <- function(logt, mu, sd, beta, sigma, alpha) {
          (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * (1 - pnorm(beta * sd  + (logt-mu)/sd))
        }
        # hazard fct
        hazard <- function(logt, mu, sd, beta, sigma, alpha) {
          # pdf of a Reed distribution
          reed_pdf <- function(logt, mu, sd, beta) {
            beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
          }
          # survival fct of a Reed distribution
          reed_phi <- function(logt, mu, sd, beta) {
            (1 - pnorm(logt, mu, sd)) - 
              exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
              (1 - pnorm(beta * sd + (logt - mu) / sd))
          }
          # hazard of a reed distribution
          reed_hz <- function(logt, mu, sd, beta) {
            reed_pdf(logt, mu, sd, beta) / reed_phi(logt, mu, sd, beta)
          }
          return(reed_hz(logt, mu, sd, beta))
        }
      }
      
      ### survival
      if(what %in% c("survival", "hazard")) {
        if(what == "survival") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 surv(logt = logt,
                                      mu = mu[row_i],
                                      sd = sd[row_i],
                                      beta = beta[row_i],
                                      alpha = alpha[row_i],
                                      sigma = sigma[row_i]
                                 )
                               }
          )
        }
        if(what == "hazard") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 hazard(logt = logt,
                                        mu = mu[row_i],
                                        sd = sd[row_i],
                                        beta = beta[row_i],
                                        alpha = alpha[row_i],
                                        sigma = sigma[row_i]
                                 )
                               }
          )
        }
      }
      else {
        stop("Must choose either 'survival' or 'hazard' for arg. 'what'")
      }
      
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Stan"
      output$year <- y
      return(output)
    }
    
    if (covariable == "cov2_1"){
      
      ### extract posteriors
      if(any(stanfit@model_pars == "beta")) {
        beta <- drop(rstan::extract(stanfit, 'beta')$beta)
        mu <- drop(rstan::extract(stanfit, 'intercept')$intercept[,1] + rstan::extract(stanfit, 'slope')$slope[,2,1]) + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,1])) * (y - mean(x_trend)) / sd(x_trend))
        sd <- drop(rstan::extract(stanfit, 'residual_sd')$residual_sd)
        alpha <- rep(0, length(mu))
        sigma <- rep(0, length(mu))
        
        # survivorship fct
        surv <- function(logt, mu, sd, beta, sigma, alpha) {
          (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * (1 - pnorm(beta * sd  + (logt-mu)/sd))
        }
        # hazard fct
        hazard <- function(logt, mu, sd, beta, sigma, alpha) {
          # pdf of a Reed distribution
          reed_pdf <- function(logt, mu, sd, beta) {
            beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
          }
          # survival fct of a Reed distribution
          reed_phi <- function(logt, mu, sd, beta) {
            (1 - pnorm(logt, mu, sd)) - 
              exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
              (1 - pnorm(beta * sd + (logt - mu) / sd))
          }
          # hazard of a reed distribution
          reed_hz <- function(logt, mu, sd, beta) {
            reed_pdf(logt, mu, sd, beta) / reed_phi(logt, mu, sd, beta)
          }
          return(reed_hz(logt, mu, sd, beta))
        }
      }
      
      ### survival
      if(what %in% c("survival", "hazard")) {
        if(what == "survival") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 surv(logt = logt,
                                      mu = mu[row_i],
                                      sd = sd[row_i],
                                      beta = beta[row_i],
                                      alpha = alpha[row_i],
                                      sigma = sigma[row_i]
                                 )
                               }
          )
        }
        if(what == "hazard") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 hazard(logt = logt,
                                        mu = mu[row_i],
                                        sd = sd[row_i],
                                        beta = beta[row_i],
                                        alpha = alpha[row_i],
                                        sigma = sigma[row_i]
                                 )
                               }
          )
        }
      }
      else {
        stop("Must choose either 'survival' or 'hazard' for arg. 'what'")
      }
      
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Stan"
      output$year <- y
      return(output)
    }  
    
  }
  
  surv_table <- NULL
  hazard_table <- NULL
  
  if (what == "survival" & cov == "cov1_1") {
    
    for (y in 1:n_year) {
      
      ### Surv
      covariable <- "cov1_1"
      surv_df <- surv_stan(stanfit = stanfit, 
                                 logt = age_modelise_surv,
                                 what = "survival",
                                 y = y
                       )
      
      surv_df$cov <- covariable
      
      surv_table <- rbind(surv_table, surv_df)
      
    }
    
    return(surv_table)
    
  }
  
  if (what == "hazard" & cov == "cov1_1") { 
    for (y in 1:n_year){
      covariable <- "cov1_1"
      ### Hazard
      hazard_df <- surv_stan(stanfit = stanfit, 
                             logt = age_modelise_surv,
                             what = "hazard",
                             y = y
      )
      hazard_df$cov <- covariable
      
      hazard_table <- rbind(hazard_table, hazard_df)
      
    }
    
    return(hazard_table)
    
  }
  
  if (what == "survival" & cov == "cov1_0") {
    for (y in 1:n_year){
      ### Surv
      covariable <- "cov1_0"
      surv_df <- surv_stan(stanfit = stanfit, 
                                 logt = age_modelise_surv,
                                 what = "survival",
                                 y = y
      )
      
      surv_df$cov <- covariable
      
      surv_table <- rbind(surv_table, surv_df)
    }
    
    return(surv_table)
  }
  
  if (what == "hazard" & cov == "cov1_0") {
    for (y in 1:n_year){
      covariable <- "cov1_0"
      ### Hazard
      hazard_df <- surv_stan(stanfit = stanfit, 
                             logt = age_modelise_surv,
                             what = "hazard",
                             y = y
      )
      
      hazard_df$cov <- covariable
      
      hazard_table <- rbind(hazard_table, hazard_df)
      
    }
    return(hazard_table)
  }
  
  if (what == "survival" & cov == "cov2_1") {
    
    for (y in 1:n_year){
      
      ### Surv
      covariable <- "cov2_1"
      surv_df <- surv_stan(stanfit = stanfit, 
                                 logt = age_modelise_surv,
                                 what = "survival",
                                 y = y
      )
      
      surv_df$cov <- covariable
      
      surv_table <- rbind(surv_table, surv_df)
      
    }
    
    return(surv_table)
    
  }
  
  if (what == "hazard" & cov == "cov2_1") { 
    for (y in 1:n_year){
      covariable <- "cov2_1"
      ### Hazard
      hazard_df <- surv_stan(stanfit = stanfit, 
                             logt = age_modelise_surv,
                             what = "hazard",
                             y = y
      )
      hazard_df$cov <- covariable
      
      hazard_table <- rbind(hazard_table, hazard_df)
      
    }
    
    return(hazard_table)
    
  }
  
  if (what == "survival" & cov == "cov2_0") {
    for (y in 1:n_year){
      ### Surv
      covariable <- "cov2_0"
      surv_df <- surv_stan(stanfit = stanfit, 
                                 logt = age_modelise_surv,
                                 what = "survival",
                                 y = y
      )
      
      surv_df$cov <- covariable
      
      surv_table <- rbind(surv_table, surv_df)
    }
    
    return(surv_table)
  }
  
  if (what == "hazard" & cov == "cov2_0") {
    for (y in 1:n_year){
      covariable <- "cov2_0"
      ### Hazard
      hazard_df <- surv_stan(stanfit = stanfit, 
                             logt = age_modelise_surv,
                             what = "hazard",
                             y = y
      )
      
      hazard_df$cov <- covariable
      
      hazard_table <- rbind(hazard_table, hazard_df)
      
    }
    return(hazard_table)
  }
  
}

### Repro
get_repro <- function(stanfit, t = NULL, alpha = 0.20, n_year = NULL, type = "survival") {
  
  fec_table <- NULL
  
  ### Get useful function
  repro <- function(t, stanfit, type = "hazard", covariable = NULL, y = NULL) {
    # for getting confidence interval
    get_ci <- function(x) {
      c(coda::HPDinterval(coda::as.mcmc(x), prob = 0.80)[1], 
        mean(x), 
        coda::HPDinterval(coda::as.mcmc(x), prob = 0.80)[2])
    }
    
    ### hazard
    hz_fec <- function(t, mu, alpha) {
      (alpha/(exp(-mu/alpha))^alpha) * t^(alpha - 1)
    }
    
    ### Obtention taux fecondite
    taux_fec <- function(t, mu, TG, alpha) {
      (1-(exp(-(t/exp(-mu/alpha))^alpha))) * TG
    }
    
    ### survival
    surv_fec <- function(t, mu, alpha) {
      1 - (exp(-(t/exp(-mu/alpha))^alpha))
    }
    
    ### probability density function
    pdf_fec <- function(t, mu, alpha) {
      (alpha * t^(alpha - 1)/(exp(-mu/alpha))^alpha) * exp (-(t/exp(-mu/alpha))^alpha)
    }
    
      ### Asignements of values
      mu <- rstan::extract(stanfit, 'intercept')$intercept[,2] + drop(rstan::extract(stanfit, 'sigma_year')$sigma_year[,2]) * (drop(rstan::extract(stanfit, 'u_year')$u_year[, y,2]))
      alpha <- rstan::extract(stanfit, 'alpha')$alpha
      
      # compute hazard at t
      if(type == "hazard") {
        y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
                             function(row_i) {
                               hz_fec(t = t,
                                      mu = mu[row_i], 
                                      alpha = alpha[row_i]
                               )
                             }
        )
        y_pred <- do.call('rbind', y_pred)
        output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
        names(output) <- c("lower", "mean", "upper")
        output$age <- t
        output$year <- y
      }
      
      # compute cumulative survival at t
      if(type == "survival") {
        y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
                             function(row_i) {
                               surv_fec(t = t,
                                        mu = mu[row_i], 
                                        alpha = alpha[row_i]
                               )
                             }
        )
        y_pred <- do.call('rbind', y_pred)
        output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
        names(output) <- c("lower", "mean", "upper")
        output$age <- t
        output$year <- y
      }
      
      # compute fecondity rate
      if(type == "fecondite") {
        y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
                             function(row_i) {
                               taux_fec(t = t,
                                        mu = mu[row_i], 
                                        alpha = alpha[row_i],
                                        TG = rstan::extract(stanfit, 'TG')$TG[row_i]
                               )
                             }
        )
        y_pred <- do.call('rbind', y_pred)
        output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
        names(output) <- c("lower", "mean", "upper")
        output$age <- t
        output$year <- y
      }
      
      
      # compute pdf
      if(type == "pdf") {
        y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
                             function(row_i) {
                               pdf_fec(t = t,
                                       mu = rstan::extract(stanfit, 'intercept')$intercept[row_i,2] + drop(rstan::extract(stanfit, 'sigma_year')$sigma_year[,2]) * (drop(rstan::extract(stanfit, 'u_year')$u_year[,y,2])), 
                                       alpha = rstan::extract(stanfit, 'alpha')$alpha[row_i]
                               )
                             }
        )
        y_pred <- do.call('rbind', y_pred)
        output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
        names(output) <- c("lower", "mean", "upper")
        output$age <- t
        output$year <- y
      }
      return(output)
    
    
  }
  
    for (y in 1:n_year) {
      fec_df <- repro(stanfit = stanfit, 
                      t = age_modelise_fec,
                      type = "survival",
                      y = y
      )
      
      fec_df$age <- age_modelise_fec
      
      fec_table <- rbind(fec_table, fec_df)
      
    }
    return(fec_table)
  }

### Leslie matrixes
get_growth_rates <- function(data_surv = NULL, data_fec = NULL, n_year = NULL) {
  
  for (y in 1:n_year) {
    
    table <- data.frame("age" = filter(data_surv, method == "Stan" & year == y)$age, # Age-class
                        "lx" = filter(data_surv, method == "Stan" & year == y)$mean, # Cumulative survival rate
                        "fec_x" = filter(data_fec, year == y)$mean, # Proportion matures
                        "mx" = rep(NA),
                        "sx" = rep(NA)) # Age-Specific survival rate
    
    table_lower <- data.frame("age" = filter(data_surv, method == "Stan" & year == y)$age, # Age-class
                              "lx" = filter(data_surv, method == "Stan" & year == y)$lower, # Cumulative survival rate
                              "fec_x" = filter(data_fec, year == y)$lower, # Proportion matures
                              "mx" = rep(NA),
                              "sx" = rep(NA)) # Age-Specific survival rate
    
    table_upper <- data.frame("age" = filter(data_surv, method == "Stan" & year == y)$age, # Age-class
                              "lx" = filter(data_surv, method == "Stan" & year == y)$upper, # Cumulative survival rate
                              "fec_x" = filter(data_fec, year == y)$upper, # Proportion matures
                              "mx" = rep(NA),
                              "sx" = rep(NA)) # Age-Specific survival rate
    
    for (i in 1:nrow(table)) {
      
      table$sx[i] <- table$lx[i+1]/table$lx[i]
      
      if (table$age[i] >= 0 & table$age[i] < 4 ) {
        table$mx[i] <- table$fec_x[i] * (0) * table$sx[i]
      }
      if (table$age[i] >= 4 & table$age[i] < 9 ) {
        table$mx[i] <- table$fec_x[i] * (0.125) * table$sx[i]
      }
      if (table$age[i] >= 9 & table$age[i] < 13 ) {
        table$mx[i] <- table$fec_x[i] * (0.40) * table$sx[i]
      }
      if (table$age[i] >= 13 & table$age[i] < 18 ) {
        table$mx[i] <- table$fec_x[i] * (0.35) * table$sx[i]
      }
      if (table$age[i] >= 18) {
        table$mx[i] <- table$fec_x[i] * (0.30) * table$sx[i]
      }
    }
    
    for (i in 1:nrow(table_lower)) {
      
      table_lower$sx[i] <- table_lower$lx[i+1]/table_lower$lx[i]
      
      if (table_lower$age[i] >= 0 & table_lower$age[i] < 4 ) {
        table_lower$mx[i] <- table_lower$fec_x[i] * (0) * table_lower$sx[i]
      }
      if (table_lower$age[i] >= 4 & table_lower$age[i] < 9 ) {
        table_lower$mx[i] <- table_lower$fec_x[i] * (0.125) * table_lower$sx[i]
      }
      if (table_lower$age[i] >= 9 & table_lower$age[i] < 13 ) {
        table_lower$mx[i] <- table_lower$fec_x[i] * (0.40) * table_lower$sx[i]
      }
      if (table_lower$age[i] >= 13 & table_lower$age[i] < 18 ) {
        table_lower$mx[i] <- table_lower$fec_x[i] * (0.35) * table_lower$sx[i]
      }
      if (table_lower$age[i] >= 18) {
        table_lower$mx[i] <- table_lower$fec_x[i] * (0.30) * table_lower$sx[i]
      }
    }
    
    for (i in 1:nrow(table_upper)) {
      
      table_upper$sx[i] <- table_upper$lx[i+1]/table_upper$lx[i]
      
      if (table_upper$age[i] >= 0 & table_upper$age[i] < 4 ) {
        table_upper$mx[i] <- table_upper$fec_x[i] * (0) * table_upper$sx[i]
      }
      if (table_upper$age[i] >= 4 & table_upper$age[i] < 9 ) {
        table_upper$mx[i] <- table_upper$fec_x[i] * (0.125) * table_upper$sx[i]
      }
      if (table_upper$age[i] >= 9 & table_upper$age[i] < 13 ) {
        table_upper$mx[i] <- table_upper$fec_x[i] * (0.40) * table_upper$sx[i]
      }
      if (table_upper$age[i] >= 13 & table_upper$age[i] < 18 ) {
        table_upper$mx[i] <- table_upper$fec_x[i] * (0.35) * table_upper$sx[i]
      }
      if (table_upper$age[i] >= 18) {
        table_upper$mx[i] <- table_upper$fec_x[i] * (0.30) * table_upper$sx[i]
      }
    }
    
    matrix <- leslie.matrix(lx = table$lx, # vector of either age-specific cumulative survival or person-years lived in the interval
                            mx = table$mx, # age-specific fertility rates
                            L = FALSE, # if ’FALSE’, lx is taken to be cumulative survival to exact age x + n
                            peryear = 1, # Multiplier for fertility. Defaults to peryear = 5
                            one.sex = TRUE, # FALSE, # if ’TRUE’, fertility rates will be divided by (1+SRB)
                            infant.class = TRUE # ’TRUE’ if lx contains a value for the infant age-class
    )
    
    matrix_lower <- leslie.matrix(lx = table_lower$lx, # vector of either age-specific cumulative survival or person-years lived in the interval
                                  mx = table_lower$mx, # age-specific fertility rates
                                  L = FALSE, # if ’FALSE’, lx is taken to be cumulative survival to exact age x + n
                                  peryear = 1, # Multiplier for fertility. Defaults to peryear = 5
                                  one.sex = TRUE, # FALSE, # if ’TRUE’, fertility rates will be divided by (1+SRB)
                                  infant.class = TRUE # ’TRUE’ if lx contains a value for the infant age-class
    )
    
    matrix_upper <- leslie.matrix(lx = table_upper$lx, # vector of either age-specific cumulative survival or person-years lived in the interval
                                  mx = table_upper$mx, # age-specific fertility rates
                                  L = FALSE, # if ’FALSE’, lx is taken to be cumulative survival to exact age x + n
                                  peryear = 1, # Multiplier for fertility. Defaults to peryear = 5
                                  one.sex = TRUE, # FALSE, # if ’TRUE’, fertility rates will be divided by (1+SRB)
                                  infant.class = TRUE # ’TRUE’ if lx contains a value for the infant age-class
    )
    
    l <- eigen.analysis(matrix)$lambda
    l_low <- eigen.analysis(matrix_lower)$lambda
    l_upp <- eigen.analysis(matrix_upper)$lambda
    rho <- eigen.analysis(matrix)$rho
    
    proj_dd <- data.frame("Percentage" = l^(0:time_serie),
                          "Time" = seq(0,100,1),
                          "Cohort" = rep(y),
                          "Lambda" = l,
                          "Lambda_low" = l_low,
                          "Lambda_upp" = l_upp,
                          "Rho" = rho)
    
    proj_years <- rbind(proj_years, proj_dd)
    
  }
  
  table_lambda <- proj_years %>% 
    select(- Percentage, - Time) %>%
    unique()
  
  return(table_lambda)
  
}
