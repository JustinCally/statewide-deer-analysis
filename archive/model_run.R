# read in data
model_data_file <- "data/model_data.rds"
model_data <- readRDS(model_data_file)

# list deer species to model
deer_species_all <- c("Cervus unicolor", "Dama dama", "Cervus elaphus", "Axis porcinus")

# mcmc params
ni <- 750 # sampling iterations
nw <- 750 # warmup iterations
nc <- 8 # number of chains

# read in models
model_poisson <- cmdstan_model(here::here("stan", "count_det_nondet_poisson_hs.stan"))
model_poisson_co <- cmdstan_model(here::here("stan", "count_only_poisson_hs.stan"))

# specify which models to run (hog deer cant use transect info)
integrated_model <- c(TRUE, TRUE, TRUE, FALSE)
distributions <- c("poisson")
model_fits <- list()

for(j in 1:length(distributions)) {

  for(i in 1:length(deer_species_all)) {

    if(integrated_model[i]) {
      model_to_fit <- get(paste0("model_", distributions[j]))
    } else {
      model_to_fit <- get(paste0("model_", distributions[j], "_co"))
    }

    model_fits[[i]] <- model_to_fit$sample(data = model_data[[i]],
                                           chains = nc,
                                           parallel_chains = nc,
                                           init = 0.1,
                                           max_treedepth = 10,
                                           refresh = 20,
                                           step_size = 0.001,
                                           adapt_delta = 0.999,
                                           show_messages = TRUE,
                                           save_warmup = FALSE,
                                           iter_sampling = ni,
                                           iter_warmup = nw)

    model_fits[[i]]$save_object(paste0("outputs/models/fit_",
                                       distributions[j],
                                       "_",
                                       deer_species_all[i],".rds"))

  }

}
