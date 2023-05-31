lapply(c("tidyverse", "rstan", "loo"), library, character.only = TRUE)

model_index <- 1:8

source("data/fit_models.r")

