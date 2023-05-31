lapply(c("tidyverse", "rstan", "loo"), library, character.only = TRUE)

model_index <- 41:48

source("data/fit_models.r")

