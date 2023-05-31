lapply(c("tidyverse", "rstan", "loo"), library, character.only = TRUE)

model_index <- 25:32

source("data/fit_models.r")

