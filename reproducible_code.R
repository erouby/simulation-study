\section{Reproducible code}

\begin{tiny}
\begin{lstlisting}[language=R]
## ---------------------------
##
## Script name: Simulation 30 individuals (iterations 1 to 10 out of 100)
##
## Purpose of script: Perform a simulation study to evaluate the ability of a model to estimate
## decline of survival from age-at-death data in various scenarios of decline strength. 
##
## Author: Dr. Etienne Rouby, Dr. Matthieu Authier, Dr. Floriane Plard
##
## Date Created: 2023-06-20
##
## Email: etiennerouby.research@gmail.com
##
## ---------------------------
##
## Notes: This script allows to simulate 30 age-at-deaths each year for 30 years according to 5 scenarios
## of decline in survival over the time period. 
## Then it allows to fit 4 nimble on these data to estimate survival rates, including temporal effects. 
## If you want to modify the sampling size, you have to modify the parameter k at the first data simulation. 
## If you want to modify the number of years modify n_years at the first data simulation.
## The parameter l corresponds to the first simulation index. 
## We divided into 1 to 10 out of 100 because otherwise it is too long to run. 
## This script is made to run on a cluster (otherwise we would have include functions separately).
##
## ---------------------------

rm(list = ls())

### Prepare the environment
lapply(c("tidyverse", "dplyr", "nimble", "janitor"), library, character.only = TRUE)

WorkDir <- getwd() # Working directory
ModDir <- paste(WorkDir, "model", sep = "/") # Model directory
OutDir <- paste(WorkDir, "outputs", sep = "/") # Output directory
FunDir <- paste(WorkDir, "functions", sep = "/") # Output directory

# FUNCTIONS ------

### Simulate age at death
age_at_death <- function(n_year,
                         n_ind = 30,
                         effect_size = -0.2,
                         phi = c(seq(0.7, 0.9, length.out = 5), # juv
                                 seq(0.91, 0.95, length.out = 5), # ad
                                 seq(0.95, 0.75, length.out = 20) # senescence
                         ), 
                         frailty = 0
) {
  
  ### need to include a burn-in to reach a stationary distribution
  n_tot <-  length(phi) + n_year
  
  ### year effect
  upsilon <- exp(c(rep(0, length(phi)), # burn-in
                   seq(0, effect_size, length.out = n_year)
  )
  )
  
  ## Phi : age (row) by year (column) effect
  year_phi <- matrix(NA, nrow = n_tot, ncol = n_tot)
  for(t in 1:n_tot) {
    n <- min(length(phi), n_tot - t + 1)
    year_phi[t, t:n_tot] <- 0
    year_phi[t, t + 0:(n - 1)] <- phi[1:n]
    year_phi[t, ] <- year_phi[t, ] * upsilon
  }
  
  ### frailty
  ind_phi <- function(p) {
    p <- matrix(rep(log(p / (1 - p)), 
                    each = n_ind
    ), 
    ncol = length(p), 
    byrow = FALSE
    ) +
      matrix(rep(frailty * rnorm(n_ind), 
                 each = length(p)
      ),
      nrow = n_ind, byrow = TRUE
      )
    p <- plogis(p)
    return(p)
  }
  
  all_phi <- lapply(1:nrow(year_phi), 
                    function(i) {
                      ind_phi(year_phi[i, ])
                    }
  )
  
  ### survival
  longevity <- function(datamat) {
    out1 <- do.call('rbind',
                    apply(datamat, 1, 
                          function(row_i) {
                            row_ina = row_i[!is.na(row_i)]
                            S <- cumprod(rbinom(length(row_ina), 
                                                size = 1, 
                                                prob = row_ina
                            )
                            )
                            
                            out0 <- data.frame(cohort_birth = min(which(!is.na(row_i))) - length(phi),
                                               age_at_death = ifelse(prod(S) == 0, sum(S), NA)
                            )
                            return(out0)
                          }
                    )
    )
    return(out1)
  }
  out <- do.call('rbind',
                 lapply(all_phi, longevity)
  )
  return(out)
}

### Get MCMC results
get_chain <- function(model = NULL, chain = NULL, iter = NA) {
  
  if(model == "Mn"){
    fit <- M1_fit
  }
  if(model == "Mr"){
    fit <- M2_fit
  }
  if(model == "Mt"){
    fit <- M3_fit
  }
  if(model == "Mrt"){
    fit <- M4_fit
  }
  
  a <- pivot_longer(as.data.frame(fit$samples[[chain]][,str_detect(colnames(fit$samples[[chain]]), 
                                                                   "(Z|lambda|tau|unscaled_sigma2)", 
                                                                   negate = TRUE)]), 
                    cols = everything()) %>%
    mutate("chain" = rep(chain),
           "model" = rep(model),
           iter = iter)
  
  return(a)
  
}

### Survival functions

# Kaplan-Meier survival
surv_km <- function(y_sim) {
  with(survival::survfit(survival::Surv(Value, Status) ~ 1, 
                         data = data.frame(Value = y_sim, 
                                           Status = rep(1, length(y_sim))
                         )
  ),
  data.frame(lower = lower,
             mean = surv,
             upper = upper,
             age = time,
             method = "KM"
  )
  )
}

# Survival function without effects (Model Null: Mn)
surv_nimble <- function(rb = NULL, logt = NULL, alpha = 0.20, y_sim = NULL, what = "survival") {
  
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
  
  ### extract posteriors
  
  beta <- 1/rb[,,2]
  mu <- rb[,,3]
  sd <- rb[,,1]
  alpha <- rep(0, length(mu))
  sigma <- rep(0, length(mu))
  
  # survivorship fct
  surv <- function(logt, mu, sd, beta, sigma, alpha) {
    (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * 
      (1 - pnorm(beta * sd  + (logt-mu)/sd))
  }
  # hazard fct
  hazard <- function(logt, mu, sd, beta, sigma, alpha) {
    # pdf of a Reed distribution
    reed_pdf <- function(logt, mu, sd, beta) {
      beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
        (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
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
  
  ### survival
  if(what %in% c("survival", "hazard", "ann_surv")) {
    if(what == "survival" | what == "ann_surv") {
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
  
  if(what %in% c("survival", "hazard", "ann_surv")) {
    if(what == "survival") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Nimble"
      output$year <- "none"
      output$vital_rate <- "survivorship"
      return(output)
    }
    if(what == "hazard") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Nimble"
      output$year <- "none"
      output$vital_rate <- "hazard"
      return(output)
    }
    if(what == "ann_surv") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      
      a <- NULL
      b <- NULL
      new_y_pred <- NULL
      
      for (iter in 1:length(mu)) {
        for (s in 1:(length(logt)-1)) {
          a <- y_pred[iter, s+1]/y_pred[iter, s]
          
          b <- cbind(b, a)
        }
        new_y_pred <- rbind(new_y_pred, b)
        
        rm(a, b)
        
        a <- NULL
        b <- NULL
      }
      colnames(new_y_pred) <- seq(1:ncol(new_y_pred))
      output <- as.data.frame(t(apply(new_y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt[-length(logt)])
      output$method <- "Nimble"
      output$year <- "none"
      output$vital_rate <- "ann_surv"
      return(output)
    }
  }
}

# Survival function which includes a random year effect (Model Random: Mr)
surv_nimble_year <- function(rb = NULL, logt = NULL, alpha = 0.20, 
                             y_sim = NULL, y = NULL, what = "survival") {
  
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
  
  ### extract posteriors
  
  beta <- 1/rb[,,2]
  mu <- rb[,,3+y]
  sd <- rb[,,1]
  alpha <- rep(0, length(mu))
  sigma <- rep(0, length(mu))
  
  # survivorship fct
  surv <- function(logt, mu, sd, beta, sigma, alpha) {
    (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * 
      (1 - pnorm(beta * sd  + (logt-mu)/sd))
  }
  # hazard fct
  hazard <- function(logt, mu, sd, beta, sigma, alpha) {
    # pdf of a Reed distribution
    reed_pdf <- function(logt, mu, sd, beta) {
      beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
        (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
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
  
  
  ### survival
  if(what %in% c("survival", "hazard", "ann_surv")) {
    if(what == "survival" | what == "ann_surv") {
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
  
  if(what %in% c("survival", "hazard", "ann_surv")) {
    if(what == "survival") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Nimble"
      output$year <- y
      output$vital_rate <- "survivorship"
      return(output)
    }
    if(what == "hazard") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Nimble"
      output$year <- y
      output$vital_rate <- "hazard"
      return(output)
    }
    if(what == "ann_surv") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      
      a <- NULL
      b <- NULL
      new_y_pred <- NULL
      
      for (iter in 1:length(mu)) {
        for (s in 1:(length(logt)-1)) {
          a <- y_pred[iter, s+1]/y_pred[iter, s]
          
          b <- cbind(b, a)
        }
        new_y_pred <- rbind(new_y_pred, b)
        
        rm(a, b)
        
        a <- NULL
        b <- NULL
      }
      colnames(new_y_pred) <- seq(1:ncol(new_y_pred))
      output <- as.data.frame(t(apply(new_y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt[-length(logt)])
      output$method <- "Nimble"
      output$year <- y
      output$vital_rate <- "ann_surv"
      return(output)
    }
  }
}

# Survival function which includes a temporal trend effect (Model Trend: Mt)
surv_nimble_trend <- function(rb = NULL, logt = NULL, alpha = 0.20, 
                              y_sim = NULL, y = NULL, what = "survival") {
  
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
  
  ### extract posteriors
  
  beta <- 1/rb[,,2]
  mu <- rb[,,2+y] # because there is no sigma_year so the real mu is sigma_frailty as it is mu - Z/B 
  sd <- rb[,,1]
  alpha <- rep(0, length(mu))
  sigma <- rep(0, length(mu))
  
  # survivorship fct
  surv <- function(logt, mu, sd, beta, sigma, alpha) {
    (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * 
      (1 - pnorm(beta * sd  + (logt-mu)/sd))
  }
  # hazard fct
  hazard <- function(logt, mu, sd, beta, sigma, alpha) {
    # pdf of a Reed distribution
    reed_pdf <- function(logt, mu, sd, beta) {
      beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
        (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
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
  
  
  ### survival
  if(what %in% c("survival", "hazard", "ann_surv")) {
    if(what == "survival" | what == "ann_surv") {
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
  
  if(what %in% c("survival", "hazard", "ann_surv")) {
    if(what == "survival") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Nimble"
      output$year <- y
      output$vital_rate <- "survivorship"
      return(output)
    }
    if(what == "hazard") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Nimble"
      output$year <- y
      output$vital_rate <- "hazard"
      return(output)
    }
    if(what == "ann_surv") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      
      a <- NULL
      b <- NULL
      new_y_pred <- NULL
      
      for (iter in 1:length(mu)) {
        for (s in 1:(length(logt)-1)) {
          a <- y_pred[iter, s+1]/y_pred[iter, s]
          
          b <- cbind(b, a)
        }
        new_y_pred <- rbind(new_y_pred, b)
        
        rm(a, b)
        
        a <- NULL
        b <- NULL
      }
      colnames(new_y_pred) <- seq(1:ncol(new_y_pred))
      output <- as.data.frame(t(apply(new_y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt[-length(logt)])
      output$method <- "Nimble"
      output$year <- y
      output$vital_rate <- "ann_surv"
      return(output)
    }
  }
}

# Survival function which includes a random effect and a temporal trend effect (Model Random Trend: Mrt)
surv_nimble_year_trend <- function(rb = NULL, logt = NULL, alpha = 0.20, 
                                   y_sim = NULL, y = NULL, what = "survival") {
  
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
  
  ### extract posteriors
  
  beta <- 1/rb[,,2]
  mu <- rb[,,3+y]
  sd <- rb[,,1]
  alpha <- rep(0, length(mu))
  sigma <- rep(0, length(mu))
  
  # survivorship fct
  surv <- function(logt, mu, sd, beta, sigma, alpha) {
    (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * 
      (1 - pnorm(beta * sd  + (logt-mu)/sd))
  }
  # hazard fct
  hazard <- function(logt, mu, sd, beta, sigma, alpha) {
    # pdf of a Reed distribution
    reed_pdf <- function(logt, mu, sd, beta) {
      beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
        (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
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
  
  
  ### survival
  if(what %in% c("survival", "hazard", "ann_surv")) {
    if(what == "survival" | what == "ann_surv") {
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
  
  if(what %in% c("survival", "hazard", "ann_surv")) {
    if(what == "survival") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Nimble"
      output$year <- y
      output$vital_rate <- "survivorship"
      return(output)
    }
    if(what == "hazard") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Nimble"
      output$year <- y
      output$vital_rate <- "hazard"
      return(output)
    }
    if(what == "ann_surv") {
      ### formating
      y_pred <- do.call('rbind', y_pred)
      
      a <- NULL
      b <- NULL
      new_y_pred <- NULL
      
      for (iter in 1:length(mu)) {
        for (s in 1:(length(logt)-1)) {
          a <- y_pred[iter, s+1]/y_pred[iter, s]
          
          b <- cbind(b, a)
        }
        new_y_pred <- rbind(new_y_pred, b)
        
        rm(a, b)
        
        a <- NULL
        b <- NULL
      }
      colnames(new_y_pred) <- seq(1:ncol(new_y_pred))
      output <- as.data.frame(t(apply(new_y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt[-length(logt)])
      output$method <- "Nimble"
      output$year <- y
      output$vital_rate <- "ann_surv"
      return(output)
    }
  }
}

dIdataNorm <- nimbleFunction(
  run = function(x = double(1), u_year = double(1), Z = double(1),YEAR= double(1), 
                 sd = double(0), length = integer(), log = double()) {
    ll <- 0
    for(i in 1:length) {
      ll <- ll + dnorm(x[i], u_year[YEAR[i]]-Z[i], sd=sd, log=TRUE)
    }
    returnType(double())
    if(log) return(ll) else return(exp(ll))
  }
)

rIdataNorm <- nimbleFunction(
  run = function(n = integer(),  u_year = double(1), Z = double(1),YEAR= double(1), 
                 sd = double(0), length = integer()) {
    print('this should never run')
    ##x <- numeric(length)
    declare(x, double(1, length))
    returnType(double(1))
    return(x)
  }
)

dIdataNorm0 <- nimbleFunction(
  run = function(x = double(1), u_year = double(0), Z = double(1),YEAR= double(1), 
                 sd = double(0), length = integer(), log = double()) {
    ll <- 0
    for(i in 1:length) {
      ll <- ll + dnorm(x[i], u_year-Z[i], sd=sd, log=TRUE)
    }
    returnType(double())
    if(log) return(ll) else return(exp(ll))
  }
)

rIdataNorm0 <- nimbleFunction(
  run = function(n = integer(),  u_year = double(0), Z = double(1),YEAR= double(1),
                 sd = double(0), length = integer()) {
    print('this should never run')
    ##x <- numeric(length)
    declare(x, double(1, length))
    returnType(double(1))
    return(x)
  }
)

registerDistributions(list(
  dIdataNorm = list(
    BUGSdist = 'dIdataNorm(u_year, Z, YEAR, sd, length)',
    types    = c('value = double(1)', 'u_year = double(1)', 'Z = double(1)',
                 'YEAR = double(1)', 'sd = double(0)', 'length = integer()')
  ),
  dIdataNorm0 = list(
    BUGSdist = 'dIdataNorm0(u_year, Z, YEAR, sd, length)',
    types    = c('value = double(1)', 'u_year = double(0)', 'Z = double(1)',
                 'YEAR = double(1)', 'sd = double(0)', 'length = integer()')
  )
))

### Function to estimate Rhat (Gelman and Rubin convergence diagnostic)
Rhatfun <- function(rb, nch, m, long){
  vari <- apply(rb, c(2,3), var)
  mea <- apply(rb, c(2,3), mean)
  meag <- apply(rb, 3, mean)
  
  W <- apply(vari, 2, mean)
  W1 <- apply(vari, 2, var)
  B <- apply(mea, 2, var)
  cov1 <- cov2 <- rep(0, long)
  for (i in 1:long){
    cov1[i] <- cov(vari[,i], y = (mea^2)[,i])
    cov2[i] <- cov(vari[,i], y = mea[,i])
  }
  sig2 <- ((m-1)/m)*W+B
  V <- sqrt(sig2+B/3)^2
  varV <- ((m-1)/m)^2/3*W1+(4/3)^2*B^2+2*(m-1)*4/(9*m)*(cov1-2*meag*cov2)
  df <- 2*V^2/varV
  Rhat <- abs((V/W*df)/(df-2))
  return(Rhat)
}

sum_nim <- function(rb2){
  m <- dim(rb2)[1]
  nch <- dim(rb2)[2]
  long <- dim(rb2)[3]
  sumres <- matrix(NA,nrow=dim(rb2)[3],ncol=5)
  rownames(sumres) <- dimnames(rb2)[[3]]
  colnames(sumres) <- c('mean', 'sd', 'QI 2.5', 'QI 97.5', 'Rhat')
  sumres[,1] <- apply(rb2,3, mean)
  sumres[,2] <- apply(rb2,3, sd)
  sumres[,3] <- apply(rb2,3, quantile, 0.025, na.rm = TRUE)
  sumres[,4] <- apply(rb2,3, quantile, 0.975, na.rm = TRUE)
  sumres[,5] <- Rhatfun(rb2, nch, m, long)
  return(sumres)
}

# Record the start time for benchmarking
start_time <- Sys.time()

### Reproducibility
# Number of simulations
n_simul <- 10
# Set seed for reproducibility
set.seed(20220107)
seed_MC <- sample(1:1e6, size = n_simul, replace = FALSE) # One seed per iteration

# Set seed for main simulation
set.seed(20220107)

# Specify the effect size for simulation
effect_size <- -0
# Alternatively, you can use multiple effect sizes by uncommenting the line below
# effect_size <- c(-0.05, -0.10, -0.15, -0.20)

# Simulation range identifier
simulation_range <- "_1_10"

# Initialize data frames and tables
param_df <- NULL
vital_M1 <- NULL
vital_M2 <- NULL
vital_M3 <- NULL
vital_M4 <- NULL
WAIC_table <- NULL
chain_table <- NULL

# --------------------------------------------------------------------------------------------------------------
# FIRST DATA SIMULATION
# --------------------------------------------------------------------------------------------------------------

# Loop through each effect size
for (e in effect_size) {
  
  effect <- e
  
  # Invariant parameters
  n_year <- 30
  k <- 30
  l <- 1
  
  # Age at death parameters
  phi <- c(seq(0.7, 0.9, length.out = 5), # juvenile
           seq(0.91, 0.95, length.out = 5), # adult
           seq(0.95, 0.75, length.out = 20)) # senescence
  
  ntrue = 1100
  n_ind_surv = 1050
  
  ### Data simulation
  # Ensure that the number of simulated individuals meets the specified condition
  while (ntrue > n_ind_surv) {
    # Simulate age at death
    raw_simulated_dataset <- age_at_death(n_ind = k, 
                                          phi = phi,
                                          n_year = n_year, 
                                          effect_size = effect)
    
    # Process and filter the simulated dataset
    simulated_dataset <- raw_simulated_dataset %>%
      mutate(year_of_stranding = age_at_death + cohort_birth) %>%
      filter(!is.na(age_at_death),
             year_of_stranding > 0) %>%
      arrange(year_of_stranding, age_at_death) %>%
      select(year_of_stranding, age_at_death, cohort_birth)
    
    # Adjust age_at_death for individuals with age 0
    for (i in 1:nrow(simulated_dataset)) {
      if (simulated_dataset$age_at_death[i] == 0) {
        simulated_dataset$age_at_death[i] <- 0.001
      }
    }
    
    # Add additional columns to the dataset
    simulated_dataset$n_ind <- k
    simulated_dataset$iter <- l
    
    # Update the count of true individuals
    ntrue = nrow(simulated_dataset)
  }
  
  # Prepare vectors for log-transformed time and year
  LOGTIME = c(log(simulated_dataset$age_at_death), rep(NA, n_ind_surv - ntrue))
  YEAR = c(simulated_dataset$year_of_stranding - min(simulated_dataset$year_of_stranding) + 
             1, rep(NA, n_ind_surv - ntrue))
  
  # Specifications for Runs
  nthin <- 10
  niter <- 10000
  nburnin <- 5000
  m <- (niter - nburnin) / nthin
  nch <- 3
  
  # Define a trend sequence
  x_trend <- seq(1:n_year)
  
  # --------------------------------------------------------------------------------------------------------------
  # FIRST MODEL SPECIFICATIONS AND COMPILING
  # --------------------------------------------------------------------------------------------------------------
  
  # --- ### MODEL 1 NULL ### --- #
  
  # parameters to nimble
  constants <- list(n_ind_surv = n_ind_surv,
                    prior_scale = 0.1
  )
  
  
  data <- list(LOGTIME = LOGTIME,
               ntrue = ntrue,
               YEAR = rep(1,n_ind_surv)
  )
  
  inits <- list(intercept = 0,
                # prop = c(1/3,1/3,1/3),
                lambda = rgamma(3, 1, 1),
                tau = 1,
                Z = rep(1,n_ind_surv),
                unscaled_sigma2 = 1
                
  )
  
  #Nimble model
  NCode_M1 <- nimbleCode({    
    
    #priors
    intercept ~ dnorm(0, sd = 5)
    
    # prop[1:3] ~ ddirch(alpha[1:3])# variance partitioning
    for(i in 1:2) {
      lambda[i] ~ dgamma(1.0, 1.0)
      prop[i] <- lambda[i] / sum(lambda[1:2])
    }
    unscaled_sigma2 ~ dgamma(0.5, 1.0)
    tau ~ dgamma(1, 1)
    
    for (i in 1:n_ind_surv) {Z[i] ~ dexp(1/sigma_frailty)}
    
    #derived parameters
    sd_res <-  sqrt(prop[2] * unscaled_sigma2 / tau) * prior_scale
    sigma_frailty <- sqrt(prop[1] * unscaled_sigma2 / tau) * prior_scale
    
    # Likelihood
    LOGTIME[1:n_ind_surv]  ~ dIdataNorm0(intercept, Z = Z[1:n_ind_surv], 
                                         YEAR = YEAR[1:n_ind_surv], sd = sd_res, length = ntrue)
    
  })
  
  parameters_1 <- c('intercept', 'sd_res', 'sigma_frailty', 'lambda', 'unscaled_sigma2', 'tau', 'Z')
  
  Rmodel_1 <- nimbleModel(NCode_M1, constants, data, inits, check = T)
  ## configure MCMC
  conf_1 <- configureMCMC(Rmodel_1, monitors = parameters_1, thin = nthin,
                          enableWAIC = T, controlWAIC = list(online = F))
  ## build MCMC
  Rmcmc_1 <- buildMCMC(conf_1)
  ## compile model and MCMC
  Cmodel_1 <- compileNimble(Rmodel_1, showCompilerOutput = T)
  Cmcmc_1 <- compileNimble(Rmcmc_1, project = Rmodel_1)
  
  # --- ### MODEL 2 YEARS ### --- #
  
  # parameters to nimble
  constants <- list(n_year = n_year, 
                    n_ind_surv = n_ind_surv,
                    prior_scale = 0.1
                    
  )
  
  data <- list(LOGTIME = LOGTIME,
               YEAR = YEAR,
               ntrue = ntrue
  )
  
  inits <- list(intercept = 0,
                # prop = c(1/3,1/3,1/3),
                lambda = rgamma(3, 1, 1),
                tau = 1,
                Z = rep(1,n_ind_surv),
                u_year = rep(0,n_year),
                unscaled_sigma2 = 1
                
  )
  
  #Nimble model
  NCode_M2 <- nimbleCode({    
    
    #priors
    intercept ~ dnorm(0, sd = 5)
    
    # prop[1:3] ~ ddirch(alpha[1:3])# variance partitioning
    for(i in 1:3) {
      lambda[i] ~ dgamma(1.0, 1.0)
      prop[i] <- lambda[i] / sum(lambda[1:3])
    }
    unscaled_sigma2 ~ dgamma(0.5, 1.0)
    tau ~ dgamma(1, 1)
    
    for (i in 1:n_ind_surv) {Z[i] ~ dexp(1/sigma_frailty)}
    for (i in 1:n_year) {u_year[i] ~ dnorm(intercept, sd = sigma_year)}
    
    #derived parameters
    sd_res <-  sqrt(prop[2] * unscaled_sigma2 / tau) * prior_scale
    sigma_frailty <- sqrt(prop[1] * unscaled_sigma2 / tau) * prior_scale
    sigma_year <- sqrt(prop[3] * unscaled_sigma2 / tau) * prior_scale
    
    # Likelihood
    # for (i in 1:n_ind_surv) {
    #   LOGTIME[i]  ~ dnorm(u_year[YEAR[i]] - Z[i], sd = sd_res)
    # }
    LOGTIME[1:n_ind_surv]  ~ dIdataNorm(u_year[1:n_year], Z = Z[1:n_ind_surv], 
                                        YEAR = YEAR[1:n_ind_surv], sd = sd_res, length = ntrue)
    
  })
  
  parameters_2 <- c('intercept', 'sd_res', 'sigma_frailty', 'sigma_year', 
                    'u_year', 'lambda', 'unscaled_sigma2', 'tau', 'Z')
  
  Rmodel_2 <- nimbleModel(NCode_M2, constants, data, inits, check = T)
  ## configure MCMC
  conf_2 <- configureMCMC(Rmodel_2, monitors = parameters_2, thin = nthin,
                          enableWAIC = T, controlWAIC = list(online = F))
  ## build MCMC
  Rmcmc_2 <- buildMCMC(conf_2)
  ## compile model and MCMC
  Cmodel_2 <- compileNimble(Rmodel_2, showCompilerOutput = T)
  Cmcmc_2 <- compileNimble(Rmcmc_2, project = Rmodel_2)
  
  # --- ### MODEL 3 TREND ### --- #
  
  # parameters to nimble
  constants <- list(n_year = n_year, 
                    n_ind_surv = n_ind_surv,
                    prior_scale = 0.1
                    
  )
  
  data <- list(LOGTIME = LOGTIME,
               YEAR = YEAR,
               ntrue = ntrue,
               stdYEAR = (x_trend - mean(x_trend))/sd(x_trend)
  )
  
  inits <- list(intercept = 0,
                # prop = c(1/3,1/3,1/3),
                lambda = rgamma(3, 1, 1),
                tau = 1,
                b_trend = 0,
                Z = rep(1,n_ind_surv),
                u_year = rep(0,n_year),
                unscaled_sigma2 = 1
                
  )
  
  #Nimble model
  NCode_M3 <- nimbleCode({    
    
    #priors
    intercept ~ dnorm(0, sd = 5)
    
    # prop[1:3] ~ ddirch(alpha[1:3])# variance partitioning
    for(i in 1:2) {
      lambda[i] ~ dgamma(1.0, 1.0)
      prop[i] <- lambda[i] / sum(lambda[1:2])
    }
    unscaled_sigma2 ~ dgamma(0.5, 1.0)
    tau ~ dgamma(1, 1)
    
    for (i in 1:n_ind_surv) {Z[i] ~ dexp(1/sigma_frailty)}
    b_trend ~ dnorm(0, sd = 1)
    for (i in 1:n_year) {u_year[i] <- intercept + b_trend * stdYEAR[i]}
    
    #derived parameters
    sd_res <-  sqrt(prop[2] * unscaled_sigma2 / tau) * prior_scale
    sigma_frailty <- sqrt(prop[1] * unscaled_sigma2 / tau) * prior_scale
    
    # Likelihood
    # for (i in 1:n_ind_surv) {
    #   LOGTIME[i]  ~ dnorm(u_year[YEAR[i]] - Z[i], sd = sd_res)
    # }
    LOGTIME[1:n_ind_surv]  ~ dIdataNorm(u_year[1:n_year], Z = Z[1:n_ind_surv], 
                                        YEAR = YEAR[1:n_ind_surv], sd = sd_res, length = ntrue)
  })
  
  parameters_3 <- c('intercept', 'sd_res', 'sigma_frailty', 'u_year', 
                    'b_trend', 'lambda', 'unscaled_sigma2', 'tau', 'Z')
  
  Rmodel_3 <- nimbleModel(NCode_M3, constants, data, inits, check = T)
  ## configure MCMC
  conf_3 <- configureMCMC(Rmodel_3, monitors = parameters_3, thin = nthin,
                          enableWAIC = T, controlWAIC = list(online = F))
  ## build MCMC
  Rmcmc_3 <- buildMCMC(conf_3)
  ## compile model and MCMC
  Cmodel_3 <- compileNimble(Rmodel_3, showCompilerOutput = T)
  Cmcmc_3 <- compileNimble(Rmcmc_3, project = Rmodel_3)
  
  # --- ### MODEL 4 YEARS TREND ### --- #
  
  # Options parallelisation
  options(mc.cores = parallel::detectCores())
  
  # parameters to nimble
  constants <- list(n_year = n_year, 
                    n_ind_surv = n_ind_surv,
                    prior_scale = 0.1
                    
  )
  
  data <- list(LOGTIME = LOGTIME,
               YEAR = YEAR,
               ntrue = ntrue,
               stdYEAR = (x_trend - mean(x_trend))/sd(x_trend)
  )
  
  inits <- list(intercept = 0,
                # prop = c(1/3,1/3,1/3),
                lambda = rgamma(3, 1, 1),
                tau = 1,
                b_trend = 0,
                Z = rep(1,n_ind_surv),
                u_year = rep(0,n_year),
                unscaled_sigma2 = 1
                
  )
  
  #Nimble model
  NCode_M4 <- nimbleCode({    
    
    #priors
    intercept ~ dnorm(0, sd = 5)
    
    # prop[1:3] ~ ddirch(alpha[1:3])# variance partitioning
    for(i in 1:3) {
      lambda[i] ~ dgamma(1.0, 1.0)
      prop[i] <- lambda[i] / sum(lambda[1:3])
    }
    unscaled_sigma2 ~ dgamma(0.5, 1.0)
    tau ~ dgamma(1, 1)
    
    for (i in 1:n_ind_surv) {Z[i] ~ dexp(1/sigma_frailty)}
    b_trend ~ dnorm(0, sd = 1)
    for (i in 1:n_year) {u_year[i] ~ dnorm(intercept + b_trend * stdYEAR[i], sd = sigma_year)}
    
    #derived parameters
    sd_res <-  sqrt(prop[2] * unscaled_sigma2 / tau) * prior_scale
    sigma_frailty <- sqrt(prop[1] * unscaled_sigma2 / tau) * prior_scale
    sigma_year <- sqrt(prop[3] * unscaled_sigma2 / tau) * prior_scale
    
    # Likelihood
    LOGTIME[1:n_ind_surv]  ~ dIdataNorm(u_year[1:n_year], Z = Z[1:n_ind_surv], 
                                        YEAR = YEAR[1:n_ind_surv], sd = sd_res, length = ntrue)
  })
  
  parameters_4 <- c('intercept', 'sd_res', 'sigma_frailty', 'sigma_year', 'u_year', 'b_trend', 'lambda', 'unscaled_sigma2', 'tau', 'Z')
  
  Rmodel_4 <- nimbleModel(NCode_M4, constants, data, inits, check = T)
  ## configure MCMC
  conf_4 <- configureMCMC(Rmodel_4, monitors = parameters_4, thin = nthin,
                          enableWAIC = T, controlWAIC = list(online = F))
  ## build MCMC
  Rmcmc_4 <- buildMCMC(conf_4)
  ## compile model and MCMC
  Cmodel_4 <- compileNimble(Rmodel_4, showCompilerOutput = T)
  Cmcmc_4 <- compileNimble(Rmcmc_4, project = Rmodel_4)
  
  # --------------------------------------------------------------------------------------------------------------
  # FIRST MODEL FIT
  # --------------------------------------------------------------------------------------------------------------
  
  M1_fit <- runMCMC(Cmcmc_1, niter = niter, nburnin = nburnin, 
                    nchains = nch, progressBar = F, summary = F, WAIC = T)
  
  M2_fit <- runMCMC(Cmcmc_2, niter = niter, nburnin = nburnin, 
                    nchains = nch, progressBar = F, summary = F, WAIC = T)
  
  M3_fit <- runMCMC(Cmcmc_3, niter = niter, nburnin = nburnin, 
                    nchains = nch, progressBar = F, summary = F, WAIC = T)
  
  M4_fit <- runMCMC(Cmcmc_4, niter = niter, nburnin = nburnin, 
                    nchains = nch, progressBar = F, summary = F, WAIC = T)
  
  # M1
  
  rb_1 <- array(NA, dim = c(m, nch, ncol(M1_fit$samples$chain1)))
  rb_1[,1,] <- M1_fit$samples$chain1
  rb_1[,2,] <- M1_fit$samples$chain2
  rb_1[,3,] <- M1_fit$samples$chain3
  
  dimnames(rb_1)[[3]] <- colnames(M1_fit$samples$chain1)
  
  res_1 <- as.data.frame(sum_nim(rb_1)) %>%
    rownames_to_column(var = "param") %>%
    janitor::clean_names() 
  
  res_1$n_ind <- k
  res_1$sim <- l
  res_1$effect <- e
  res_1$model <- "M1"
  
  params <- matrix(c("sd_res", "sigma_frailty", "intercept"), ncol = 1)
  
  rb_1 <- rb_1[, , params[, 1]]
  
  # M2
  
  rb_2 <- array(NA, dim = c(m, nch, ncol(M2_fit$samples$chain1)))
  rb_2[,1,] <- M2_fit$samples$chain1
  rb_2[,2,] <- M2_fit$samples$chain2
  rb_2[,3,] <- M2_fit$samples$chain3
  
  dimnames(rb_2)[[3]] <- colnames(M2_fit$samples$chain1)
  
  res_2 <- as.data.frame(sum_nim(rb_2)) %>%
    rownames_to_column(var = "param") %>%
    janitor::clean_names() 
  
  res_2$n_ind <- k
  res_2$sim <- l
  res_2$effect <- e
  res_2$model <- "M2"
  
  u_years <- NULL
  
  for(i in 1:n_year) {
    u <- paste("u_year[", i, "]", sep = "")
    u_years <- c(u_years, u)
  }
  
  params <- matrix(c("sd_res", "sigma_frailty", "sigma_year", u_years), ncol = 1)
  
  rb_2 <- rb_2[, , params[, 1]]
  
  # M3
  
  rb_3 <- array(NA, dim = c(m, nch, ncol(M3_fit$samples$chain1)))
  rb_3[,1,] <- M3_fit$samples$chain1
  rb_3[,2,] <- M3_fit$samples$chain2
  rb_3[,3,] <- M3_fit$samples$chain3
  
  dimnames(rb_3)[[3]] <- colnames(M3_fit$samples$chain1)
  
  res_3 <- as.data.frame(sum_nim(rb_3)) %>%
    rownames_to_column(var = "param") %>%
    janitor::clean_names() 
  
  res_3$n_ind <- k
  res_3$sim <- l
  res_3$effect <- e
  res_3$model <- "M3"
  
  u_years <- NULL
  
  for(i in 1:n_year) {
    u <- paste("u_year[", i, "]", sep = "")
    u_years <- c(u_years, u)
  }
  
  params <- matrix(c("sd_res", "sigma_frailty", u_years), ncol = 1)
  
  rb_3 <- rb_3[, , params[, 1]]
  
  # M4
  
  rb_4 <- array(NA, dim = c(m, nch, ncol(M4_fit$samples$chain1)))
  rb_4[,1,] <- M4_fit$samples$chain1
  rb_4[,2,] <- M4_fit$samples$chain2
  rb_4[,3,] <- M4_fit$samples$chain3
  
  dimnames(rb_4)[[3]] <- colnames(M4_fit$samples$chain1)
  
  res_4 <- as.data.frame(sum_nim(rb_4)) %>%
    rownames_to_column(var = "param") %>%
    janitor::clean_names() 
  
  res_4$n_ind <- k
  res_4$sim <- l
  res_4$effect <- e
  res_4$model <- "M4"
  
  u_years <- NULL
  
  for(i in 1:n_year) {
    u <- paste("u_year[", i, "]", sep = "")
    u_years <- c(u_years, u)
  }
  
  params <- matrix(c("sd_res", "sigma_frailty", "sigma_year", u_years), ncol = 1)
  
  rb_4 <- rb_4[, , params[, 1]]
  
  param_df <- rbind(param_df, res_1, res_2, res_3, res_4)
  
  WAIC_df <- data.frame(model = c("M1", "M2", "M3", "M4"),
                        WAIC = c(M1_fit$WAIC, M2_fit$WAIC, M3_fit$WAIC, M4_fit$WAIC),
                        sim = rep(l),
                        n_ind = rep(k))
  
  WAIC_table <- rbind(WAIC_table, WAIC_df)
  
  rm(WAIC_df)
  
  chains <- rbind(rbind(get_chain(model = "Mn", chain = "chain1", iter = l), 
                        get_chain(model = "Mn", chain = "chain2", iter = l), 
                        get_chain(model = "Mn", chain = "chain3", iter = l)), 
                  rbind(get_chain(model = "Mr", chain = "chain1", iter = l), 
                        get_chain(model = "Mr", chain = "chain2", iter = l), 
                        get_chain(model = "Mr", chain = "chain3", iter = l)),
                  rbind(get_chain(model = "Mt", chain = "chain1", iter = l),
                        get_chain(model = "Mt", chain = "chain2", iter = l),
                        get_chain(model = "Mt", chain = "chain3", iter = l)),
                  rbind(get_chain(model = "Mrt", chain = "chain1", iter = l),
                        get_chain(model = "Mrt", chain = "chain2", iter = l),
                        get_chain(model = "Mrt", chain = "chain3", iter = l))
  )
  
  chain_table <- rbind(chain_table, chains)
  
  # --------------------------------------------------------------------------------------------------------------
  # VITAL RATES
  # --------------------------------------------------------------------------------------------------------------
  
  age_modelise_surv <- log(seq(0.001, 31))
  
  vital_M1_dd <- rbind(surv_nimble(rb = rb_1,
                                   logt = age_modelise_surv,
                                   what = "ann_surv"),
                       surv_nimble(rb = rb_1,
                                   logt = age_modelise_surv,
                                   what = "hazard"),
                       surv_nimble(rb = rb_1,
                                   logt = age_modelise_surv,
                                   what = "survival")
  )
  
  vital_M1_dd$sim <- l
  vital_M1_dd$n_ind <- k
  vital_M1_dd$effect <- e
  vital_M1_dd$model <- "M1"
  
  vital_M1 <- rbind(vital_M1, vital_M1_dd)
  
  for (y in 1:n_year) {
    
    ### M2
    vital_M2_dd <- rbind(surv_nimble_year(rb = rb_2,
                                          logt = age_modelise_surv,
                                          what = "ann_surv", 
                                          y = y),
                         surv_nimble_year(rb = rb_2,
                                          logt = age_modelise_surv,
                                          what = "hazard", 
                                          y = y),
                         surv_nimble_year(rb = rb_2,
                                          logt = age_modelise_surv,
                                          what = "survival", 
                                          y = y)
    )
    
    vital_M2_dd$sim <- l
    vital_M2_dd$n_ind <- k
    vital_M2_dd$effect <- e
    vital_M2_dd$model <- "M2"
    
    vital_M2 <- rbind(vital_M2, vital_M2_dd)
    
    ### M3
    vital_M3_dd <- rbind(surv_nimble_trend(rb = rb_3,
                                           logt = age_modelise_surv,
                                           what = "ann_surv", 
                                           y = y),
                         surv_nimble_trend(rb = rb_3,
                                           logt = age_modelise_surv,
                                           what = "hazard", 
                                           y = y),
                         surv_nimble_trend(rb = rb_3,
                                           logt = age_modelise_surv,
                                           what = "survival", 
                                           y = y)
    )
    
    vital_M3_dd$sim <- l
    vital_M3_dd$n_ind <- k
    vital_M3_dd$effect <- e
    vital_M3_dd$model <- "M3"
    
    vital_M3 <- rbind(vital_M3, vital_M3_dd)
    
    ### M4
    vital_M4_dd <- rbind(surv_nimble_year_trend(rb = rb_4,
                                                logt = age_modelise_surv,
                                                what = "ann_surv", 
                                                y = y),
                         surv_nimble_year_trend(rb = rb_4,
                                                logt = age_modelise_surv,
                                                what = "hazard", 
                                                y = y),
                         surv_nimble_year_trend(rb = rb_4,
                                                logt = age_modelise_surv,
                                                what = "survival", 
                                                y = y)
    )
    
    vital_M4_dd$sim <- l
    vital_M4_dd$n_ind <- k
    vital_M4_dd$effect <- e
    vital_M4_dd$model <- "M4"
    
    vital_M4 <- rbind(vital_M4, vital_M4_dd)
    
  }
  
  # --------------------------------------------------------------------------------------------------------------
  # LOOP
  # --------------------------------------------------------------------------------------------------------------
  
  start_sim <- l
  
  for (l in c(start_sim:n_simul)) {
    set.seed(seed_MC[l])
    if (l %% 10 == 0) {
      writeLines(paste("\t\tsimulation #", l, " out of ", n_simul, sep = ""))
    }
    
    # --------------------------------------------------------------------------------------------------------------
    # DATA SIMULATION
    # --------------------------------------------------------------------------------------------------------------
    
    ntrue = 1100
    n_ind_surv = 1050
    
    ### Data simulation
    while(ntrue >  n_ind_surv){
      raw_simulated_dataset <- age_at_death(n_ind = k, 
                                            phi = phi,
                                            n_year = n_year, 
                                            effect_size = effect)
      
      simulated_dataset <- raw_simulated_dataset %>%
        mutate(year_of_stranding = age_at_death + cohort_birth) %>%
        filter(!is.na(age_at_death),
               year_of_stranding > 0, # remove burn-in
        ) %>%
        arrange(year_of_stranding, age_at_death) %>%
        select(year_of_stranding, age_at_death, cohort_birth)
      
      for (i in 1:nrow(simulated_dataset)) {
        if (simulated_dataset$age_at_death[i] == 0){
          simulated_dataset$age_at_death[i] <- 0.001
        }
      }
      
      simulated_dataset$n_ind <- k
      simulated_dataset$iter <- l
      
      ntrue = nrow(simulated_dataset)
    }
    
    LOGTIME <- c(log(simulated_dataset$age_at_death), rep(NA, n_ind_surv-ntrue))
    YEAR <- c(simulated_dataset$year_of_stranding-
                min(simulated_dataset$year_of_stranding)+1, rep(NA, n_ind_surv-ntrue))
    
    # --------------------------------------------------------------------------------------------------------------
    # MODEL SPECIFICATIONS AND COMPILING
    # --------------------------------------------------------------------------------------------------------------
    
    # --- ### MODEL 1 NULL ### --- #
    
    data <- list(LOGTIME = LOGTIME,
                 ntrue = ntrue,
                 YEAR = rep(1,n_ind_surv)
    )
    
    Cmodel_1$LOGTIME <- LOGTIME
    Cmodel_1$ntrue <- ntrue
    Cmodel_1$YEAR <- rep(1,n_ind_surv)
    
    M1_fit <- runMCMC(Cmcmc_1, niter = niter, nburnin = nburnin, 
                      nchains = nch, progressBar = F, summary = F, WAIC = T)
    
    # --- ### MODEL 2 YEARS ### --- #
    
    data <- list(LOGTIME = LOGTIME,
                 YEAR = YEAR,
                 ntrue = ntrue
    )
    
    Cmodel_2$LOGTIME <- LOGTIME
    Cmodel_2$ntrue <- ntrue
    Cmodel_2$YEAR <- YEAR
    
    M2_fit <- runMCMC(Cmcmc_2, niter = niter, nburnin = nburnin, 
                      nchains = nch, progressBar = F, summary = F, WAIC = T)
    
    # --- ### MODEL 3 TREND ### --- #
    
    data <- list(LOGTIME = LOGTIME,
                 YEAR = YEAR,
                 ntrue = ntrue,
                 stdYEAR = (x_trend - mean(x_trend))/sd(x_trend)
    )
    
    Cmodel_3$LOGTIME <- LOGTIME
    Cmodel_3$ntrue <- ntrue
    Cmodel_3$YEAR <- YEAR
    
    M3_fit <- runMCMC(Cmcmc_3, niter = niter, nburnin = nburnin, 
                      nchains = nch, progressBar = F, summary = F, WAIC = T)
    
    # --- ### MODEL 4 YEARS TREND ### --- #
    
    data <- list(LOGTIME = LOGTIME,
                 YEAR = YEAR,
                 ntrue = ntrue,
                 stdYEAR = (x_trend - mean(x_trend))/sd(x_trend)
    )
    
    Cmodel_4$LOGTIME <- LOGTIME
    Cmodel_4$ntrue <- ntrue
    Cmodel_4$YEAR <- YEAR
    
    M4_fit <- runMCMC(Cmcmc_4, niter = niter, nburnin = nburnin, 
                      nchains = nch, progressBar = F, summary = F, WAIC = T)
    
    # --------------------------------------------------------------------------------------------------------------
    # CHECK FOR CONVERGENCE AND PARAMETERS
    # --------------------------------------------------------------------------------------------------------------
    
    # M1
    
    rb_1 <- array(NA, dim = c(m, nch, ncol(M1_fit$samples$chain1)))
    rb_1[,1,] <- M1_fit$samples$chain1
    rb_1[,2,] <- M1_fit$samples$chain2
    rb_1[,3,] <- M1_fit$samples$chain3
    
    dimnames(rb_1)[[3]] <- colnames(M1_fit$samples$chain1)
    
    res_1 <- as.data.frame(sum_nim(rb_1)) %>%
      rownames_to_column(var = "param") %>%
      janitor::clean_names() 
    
    res_1$n_ind <- k
    res_1$sim <- l
    res_1$effect <- e
    res_1$model <- "M1"
    
    params <- matrix(c("sd_res", "sigma_frailty", "intercept"), ncol = 1)
    
    rb_1 <- rb_1[, , params[, 1]]
    
    # M2
    
    rb_2 <- array(NA, dim = c(m, nch, ncol(M2_fit$samples$chain1)))
    rb_2[,1,] <- M2_fit$samples$chain1
    rb_2[,2,] <- M2_fit$samples$chain2
    rb_2[,3,] <- M2_fit$samples$chain3
    
    dimnames(rb_2)[[3]] <- colnames(M2_fit$samples$chain1)
    
    res_2 <- as.data.frame(sum_nim(rb_2)) %>%
      rownames_to_column(var = "param") %>%
      janitor::clean_names() 
    
    res_2$n_ind <- k
    res_2$sim <- l
    res_2$effect <- e
    res_2$model <- "M2"
    
    u_years <- NULL
    
    for(i in 1:n_year) {
      u <- paste("u_year[", i, "]", sep = "")
      u_years <- c(u_years, u)
    }
    
    params <- matrix(c("sd_res", "sigma_frailty", "sigma_year", u_years), ncol = 1)
    
    rb_2 <- rb_2[, , params[, 1]]
    
    # M3
    
    rb_3 <- array(NA, dim = c(m, nch, ncol(M3_fit$samples$chain1)))
    rb_3[,1,] <- M3_fit$samples$chain1
    rb_3[,2,] <- M3_fit$samples$chain2
    rb_3[,3,] <- M3_fit$samples$chain3
    
    dimnames(rb_3)[[3]] <- colnames(M3_fit$samples$chain1)
    
    res_3 <- as.data.frame(sum_nim(rb_3)) %>%
      rownames_to_column(var = "param") %>%
      janitor::clean_names() 
    
    res_3$n_ind <- k
    res_3$sim <- l
    res_3$effect <- e
    res_3$model <- "M3"
    
    u_years <- NULL
    
    for(i in 1:n_year) {
      u <- paste("u_year[", i, "]", sep = "")
      u_years <- c(u_years, u)
    }
    
    params <- matrix(c("sd_res", "sigma_frailty", u_years), ncol = 1)
    
    rb_3 <- rb_3[, , params[, 1]]
    
    # M4
    
    rb_4 <- array(NA, dim = c(m, nch, ncol(M4_fit$samples$chain1)))
    rb_4[,1,] <- M4_fit$samples$chain1
    rb_4[,2,] <- M4_fit$samples$chain2
    rb_4[,3,] <- M4_fit$samples$chain3
    
    dimnames(rb_4)[[3]] <- colnames(M4_fit$samples$chain1)
    
    res_4 <- as.data.frame(sum_nim(rb_4)) %>%
      rownames_to_column(var = "param") %>%
      janitor::clean_names() 
    
    res_4$n_ind <- k
    res_4$sim <- l
    res_4$effect <- e
    res_4$model <- "M4"
    
    u_years <- NULL
    
    for(i in 1:n_year) {
      u <- paste("u_year[", i, "]", sep = "")
      u_years <- c(u_years, u)
    }
    
    params <- matrix(c("sd_res", "sigma_frailty", "sigma_year", u_years), ncol = 1)
    
    rb_4 <- rb_4[, , params[, 1]]
    
    param_df <- rbind(param_df, res_1, res_2, res_3, res_4)
    
    WAIC_df <- data.frame(model = c("M1", "M2", "M3", "M4"),
                          WAIC = c(M1_fit$WAIC, M2_fit$WAIC, M3_fit$WAIC, M4_fit$WAIC),
                          sim = rep(l),
                          n_ind = rep(k))
    
    WAIC_table <- rbind(WAIC_table, WAIC_df)
    
    rm(WAIC_df)
    
    chains <- rbind(rbind(get_chain(model = "Mn", chain = "chain1", iter = l), 
                          get_chain(model = "Mn", chain = "chain2", iter = l), 
                          get_chain(model = "Mn", chain = "chain3", iter = l)), 
                    rbind(get_chain(model = "Mr", chain = "chain1", iter = l), 
                          get_chain(model = "Mr", chain = "chain2", iter = l), 
                          get_chain(model = "Mr", chain = "chain3", iter = l)),
                    rbind(get_chain(model = "Mt", chain = "chain1", iter = l),
                          get_chain(model = "Mt", chain = "chain2", iter = l),
                          get_chain(model = "Mt", chain = "chain3", iter = l)),
                    rbind(get_chain(model = "Mrt", chain = "chain1", iter = l),
                          get_chain(model = "Mrt", chain = "chain2", iter = l),
                          get_chain(model = "Mrt", chain = "chain3", iter = l))
    )
    
    chain_table <- rbind(chain_table, chains)
    
    # --------------------------------------------------------------------------------------------------------------
    # VITAL RATES
    # --------------------------------------------------------------------------------------------------------------
    
    age_modelise_surv <- log(seq(0.001, 31))
    
    vital_M1_dd <- rbind(surv_nimble(rb = rb_1,
                                     logt = age_modelise_surv,
                                     what = "ann_surv"),
                         surv_nimble(rb = rb_1,
                                     logt = age_modelise_surv,
                                     what = "hazard"),
                         surv_nimble(rb = rb_1,
                                     logt = age_modelise_surv,
                                     what = "survival")
    )
    
    vital_M1_dd$sim <- l
    vital_M1_dd$n_ind <- k
    vital_M1_dd$effect <- e
    vital_M1_dd$model <- "M1"
    
    vital_M1 <- rbind(vital_M1, vital_M1_dd)
    
    for (y in 1:n_year) {
      
      ### M2
      vital_M2_dd <- rbind(surv_nimble_year(rb = rb_2,
                                            logt = age_modelise_surv,
                                            what = "ann_surv", 
                                            y = y),
                           surv_nimble_year(rb = rb_2,
                                            logt = age_modelise_surv,
                                            what = "hazard", 
                                            y = y),
                           surv_nimble_year(rb = rb_2,
                                            logt = age_modelise_surv,
                                            what = "survival", 
                                            y = y)
      )
      
      vital_M2_dd$sim <- l
      vital_M2_dd$n_ind <- k
      vital_M2_dd$effect <- e
      vital_M2_dd$model <- "M2"
      
      vital_M2 <- rbind(vital_M2, vital_M2_dd)
      
      ### M3
      vital_M3_dd <- rbind(surv_nimble_trend(rb = rb_3,
                                             logt = age_modelise_surv,
                                             what = "ann_surv", 
                                             y = y),
                           surv_nimble_trend(rb = rb_3,
                                             logt = age_modelise_surv,
                                             what = "hazard", 
                                             y = y),
                           surv_nimble_trend(rb = rb_3,
                                             logt = age_modelise_surv,
                                             what = "survival", 
                                             y = y)
      )
      
      vital_M3_dd$sim <- l
      vital_M3_dd$n_ind <- k
      vital_M3_dd$effect <- e
      vital_M3_dd$model <- "M3"
      
      vital_M3 <- rbind(vital_M3, vital_M3_dd)
      
      ### M4
      vital_M4_dd <- rbind(surv_nimble_year_trend(rb = rb_4,
                                                  logt = age_modelise_surv,
                                                  what = "ann_surv", 
                                                  y = y),
                           surv_nimble_year_trend(rb = rb_4,
                                                  logt = age_modelise_surv,
                                                  what = "hazard", 
                                                  y = y),
                           surv_nimble_year_trend(rb = rb_4,
                                                  logt = age_modelise_surv,
                                                  what = "survival", 
                                                  y = y)
      )
      
      vital_M4_dd$sim <- l
      vital_M4_dd$n_ind <- k
      vital_M4_dd$effect <- e
      vital_M4_dd$model <- "M4"
      
      vital_M4 <- rbind(vital_M4, vital_M4_dd)
      
    }
    
    vital_df <- rbind(vital_M1, vital_M2, vital_M3, vital_M4)
    
    ### Save the tables
    
    write.table(vital_df, paste(WorkDir, paste("vital_", k, effect, 
                                               simulation_range, ".txt", sep = ""), sep = "/" ), quote = TRUE, sep = "\t",
                row.names = FALSE, col.names = TRUE
    )
    
    write.table(param_df, paste(WorkDir, paste("param_", k, effect, 
                                               simulation_range, ".txt", sep = ""), sep = "/" ), quote = TRUE, sep = "\t",
                row.names = FALSE, col.names = TRUE
    )
    
    write.table(WAIC_table, paste(WorkDir, paste("waic_", k, effect, 
                                                 simulation_range, ".txt", sep = ""), sep = "/" ), quote = TRUE, sep = "\t",
                row.names = FALSE, col.names = TRUE
    )
    
    
    write.table(chain_table, paste(WorkDir, paste("chain_table_", k, effect, 
                                                  simulation_range, ".txt", sep = ""), sep = "/" ), quote = TRUE, sep = "\t",
                row.names = FALSE, col.names = TRUE
    )
    
  }  
}

end_time <- Sys.time()

run_time <- end_time - start_time