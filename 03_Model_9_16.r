lapply(c("tidyverse", "rstan", "loo"), library, character.only = TRUE)

model_index <- 9:16

source("data/fit_models.r")

