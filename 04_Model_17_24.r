lapply(c("tidyverse", "rstan", "loo"), library, character.only = TRUE)

model_index <- 17:24

source("data/fit_models.r")

