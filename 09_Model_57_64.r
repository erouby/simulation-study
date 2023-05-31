lapply(c("tidyverse", "rstan", "loo"), library, character.only = TRUE)

model_index <- 57:64

source("data/fit_models.r")

