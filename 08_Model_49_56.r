lapply(c("tidyverse", "rstan", "loo"), library, character.only = TRUE)

model_index <- 49:56

source("data/fit_models.r")

