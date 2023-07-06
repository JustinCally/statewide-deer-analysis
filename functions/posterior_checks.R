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
posterior_checks <- function(model, model_data, stat, title, xlims = c(0,1), integrated = F, ...) {

  # get n_obs from data
  n_obs <- model_data[["n_obs"]]
  n_obs_totals <- rowSums(n_obs)

  if(integrated) {
    n_obs[,1] <- pmax(n_obs[,1], model_data[["any_seen"]])
  }

  model_draws <- model$draws("n_obs_pred", format = "matrix")

  split_groups <- list()
  ppc_zeros <- list()
  ppc_plots <- list()
  for(i in 1:ncol(n_obs)) {
    col_sel <- paste0("n_obs_pred[", 1:nrow(n_obs), ",", i, "]")
    split_groups[[i]] <- model_draws[,col_sel]

    ppc_plots[[i]] <- bayesplot::ppc_stat(n_obs[,i],
                                          split_groups[[i]],
                                          stat=stat,
                                          ...) +
      ggplot2::ggtitle(label = title, subtitle = paste0("Group size = ", i)) +
      ggplot2::xlim(xlims) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")

  }
  cowplot::plot_grid(plotlist = ppc_plots, ncol = 2)
}
