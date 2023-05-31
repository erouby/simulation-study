




### WAIC ### --------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Make WAIC table

waic_df <- NULL

for (i in 1:64) {
  
load(paste("C:/Users/erouby/Desktop/PELAGIS/Demography/test/output/waic_", i, ".RData", sep = ""))

  waic_dd <- as.data.frame(waic_m)
  waic_dd$model <- i
  
  waic_df <- rbind(waic_df, waic_dd)
  
  rm(waic_m, waic_dd)
}

### Save WAIC table
write.csv(waic_df, "output/table_waic.csv", row.names=FALSE)

### VITAL RATES ### --------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(rstan)

### Load Best model
min_waic <- as.numeric(dplyr::filter(waic_df, waic == min(waic_df$waic))$model)
load(paste("C:/Users/erouby/Desktop/PELAGIS/Demography/test/output/model_", min_waic, ".RData", sep = ""))

source("data/vital_rates_functions.r")

### Survival
age_modelise_surv <- log(seq(0.001, 31, length = 100))
n_year <- 23
x_trend <- seq(1:n_year)

survival_rates <- NULL

for (what in c("survival", "hazard")) {
  for (cov in c("cov1_0", "cov1_1", "cov2_0", "cov2_1")) {
    table <- get_surv(stanfit = stanfit, 
                      logt = age_modelise_surv,
                      what = what,
                      n_year = n_year, 
                      cov = cov
    )
    
    table$vital_rate <- what
    table$model <- paste("Model", min_waic, sep = " ")
    
    survival_rates <- rbind(survival_rates, table)
  }
}

### Fecundity
age_modelise_fec <- seq(0.001, 31, length = 100)
prop_mature <- NULL

table <- get_repro(stanfit = stanfit, 
                   t = age_modelise_fec, 
                   alpha = 0.20, 
                   n_year = n_year, 
                   type = "survival"
  )
  
table$vital_rate <- "proportion_of_mature"
table$model <- paste("Model", min_waic, sep = " ")
  
prop_mature <- rbind(prop_mature, table)

write.table(survival_rates, paste("output/", "survival_rates.txt", sep = ""), quote = TRUE, sep = "\t",
            row.names = FALSE, col.names = TRUE
)

write.table(prop_mature, paste("output/", "fecundity_rates.txt", sep = ""), quote = TRUE, sep = "\t",
            row.names = FALSE, col.names = TRUE
)

### LESLIE MATRIXES ### ----------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(demogR)
time_serie <- 100

proj_years <- NULL
elas_years <- NULL
contrib_years <- NULL

female_surv <- survival_rates %>%
  filter(cov == "cov2_0" & vital_rate == "survival")

female_fec <- prop_mature

table_lambda <- get_growth_rates(data_surv = female_surv, data_fec = female_fec, n_year = 23)

write.table(table_lambda, paste("output/", "table_lambda.txt", sep = ""), quote = TRUE, sep = "\t",
            row.names = FALSE, col.names = TRUE
)
