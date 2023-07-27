distributions <- c("poisson")
model_fits <- list()

for(j in 1:length(distributions)) {

  for(i in 4:length(deer_species_all)) {

    if(evaluate_transects[i]) {
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
                                           step_size = 0.1,
                                           adapt_delta = 0.99,
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
