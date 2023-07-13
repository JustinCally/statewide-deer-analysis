#### Posterior predictive checks ####
# Check number of zeros (first)
# Check number of obs pred (site by group)
# prop zero
prop_zero<- function(x) mean(x == 0, na.rm = T)

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

  if(integrated) {
    n_obs_totals <- pmax(n_obs_totals, model_data[["any_seen"]])
  }

  model_draws <- model$draws("N_site_pred", format = "matrix")

  qmin <- min(quantile(apply(model_draws, 1, stat, simplify = T), 0.005, na.rm = T), do.call(stat, args = list(x = n_obs_totals[which_inc])))
  qmax <- max(quantile(apply(model_draws, 1, stat, simplify = T), 0.995, na.rm = T), do.call(stat, args = list(x = n_obs_totals[which_inc])))

    ppc_plots <- bayesplot::ppc_stat(n_obs_totals[which_inc],
                                     model_draws[,which_inc],
                                     stat=stat,
                                     ...) +
      ggplot2::ggtitle(label = title) +
      ggplot2::scale_x_continuous(limits = c(qmin, qmax)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")

    ppc_plots
}
