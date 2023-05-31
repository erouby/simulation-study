lapply(c("tidyverse", "rstan", "loo"), library, character.only = TRUE)

model_index <- 33:40

source("data/fit_models.r")

