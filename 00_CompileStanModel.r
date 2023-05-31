
library(rstan)

### compile model

test <- stan_model(file = paste("model", "final_model.stan", sep = "/"),
                   model_name = "final_model"
                   )

save(list = "test", file = "data/stanmodel.RData")
