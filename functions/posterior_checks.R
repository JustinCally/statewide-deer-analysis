#### Posterior predictive checks ####
# Check number of zeros (first)
# Check number of obs pred (site by group)
# prop zero
prop_zero<- function(x) mean(x == 0, na.rm = T)
prop_low <- function(x, bound = 20) mean(x < bound, na.rm = T)

#' Run posterior checks for CTDS distance sampling models
#'
#' @param model cmdstanr model
#' @param model_data model data generated from prepare_model_data()
#' @param stat stat passed to ppc_stat
#' @param title title for plot
#' @param ... additional arguments passed to ppc_stat
#'
#' @return plot grid
#' @export
posterior_checks <- function(model, model_data, stat, title, integrated = F, only_det = F, ...) {
  # get n_obs from data
  n_obs <- model_data[["n_obs"]]
  n_obs_totals <- n_obs %*% model_data[["gs"]] %>% as.vector()

  if(only_det) {
    which_inc <- which(n_obs_totals > 0)
  } else {
    which_inc <- seq(1:length(n_obs_totals))
  }

  # if(integrated) {
  #   n_obs_totals <- pmax(n_obs_totals, model_data[["any_seen"]])
  # }

  model_draws <- model$draws("N_site_pred", format = "matrix")

  if(identical(stat, ppc_scatter_avg)) {
    ppc_plots <- bayesplot::ppc_scatter_avg(n_obs_totals[which_inc],
                                     model_draws[,which_inc],
                                     ...) +
      ggplot2::ggtitle(label = title) +
      ggplot2::scale_x_continuous(trans = "sqrt") +
      ggplot2::scale_y_continuous(trans = "sqrt") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")
  } else if(identical(stat, ppc_rootogram)) {
    ppc_plots <- bayesplot::ppc_rootogram(y = n_obs_totals[which_inc],
                                          yrep = round(model_draws[,which_inc]),
                                            ...) +
      ggplot2::ggtitle(label = title) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")
  } else {
    ppc_plots <- bayesplot::ppc_stat(n_obs_totals[which_inc],
                                     model_draws[,which_inc],
                                     stat=stat,
                                     ...) +
      ggplot2::ggtitle(label = title) +
      # ggplot2::scale_x_continuous(limits = c(qmin, qmax)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")
  }

    return(ppc_plots)
}

posterior_checks_multispecies <- function(model, model_data, species_index,
                                          stat, title,
                                          integrated = F, only_det = F, ...) {
  # get n_obs from data
  n_obs <- model_data[["n_obs"]][,,species_index]
  n_obs_totals <- n_obs %*% model_data[["gs"]] %>% as.vector()

  if(only_det) {
    which_inc <- which(n_obs_totals > 0)
  } else {
    which_inc <- seq(1:length(n_obs_totals))
  }

  if(integrated) {
    n_obs_totals <- pmax(n_obs_totals, model_data[["any_seen"]][species_index,])
  }

  model_draws <- model$draws("N_site_pred", format = "matrix")

  which_sp <- which(stringr::str_detect(string = colnames(model_draws),
                                         pattern = paste0("N_site_pred\\[", species_index)))

  model_draws <- model_draws[,which_sp]

  if(identical(stat, ppc_scatter_avg)) {
    ppc_plots <- bayesplot::ppc_scatter_avg(n_obs_totals[which_inc],
                                            model_draws[,which_inc],
                                            ...) +
      ggplot2::ggtitle(label = title) +
      ggplot2::scale_x_continuous(trans = "sqrt") +
      ggplot2::scale_y_continuous(trans = "sqrt") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")
  } else if(identical(stat, ppc_rootogram)) {
    ppc_plots <- bayesplot::ppc_rootogram(y = n_obs_totals[which_inc],
                                          yrep = round(model_draws[,which_inc]),
                                          ...) +
      ggplot2::ggtitle(label = title) +
      ggplot2::scale_x_continuous(trans = "sqrt") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")
  } else {
    ppc_plots <- bayesplot::ppc_stat(n_obs_totals[which_inc],
                                     model_draws[,which_inc],
                                     stat=stat,
                                     ...) +
      ggplot2::ggtitle(label = title) +
      # ggplot2::scale_x_continuous(limits = c(qmin, qmax)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")
  }

  return(ppc_plots)
}
