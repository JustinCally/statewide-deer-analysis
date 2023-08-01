marginal_effects_cmd <- function(draws, param, param_number, model_data, model_column, pwr = 1, log = FALSE, transition = FALSE, abundance = FALSE) {
  require(data.table)
  intercept = paste0(param, "[1]")
  if(class(model_data[[model_column]]) == "factor") {
    levels <- length(levels(model_data[[model_column]]))
    sr <- c(1)
    factor <- TRUE
  } else {
    sr <- seq(from = -1, to = 1, by = 0.1)
    factor <- FALSE
  }
  poi = paste0(param, "[", param_number, "]")
  model_draws <- draws[, c(intercept, poi)]
  # If predictors are scaled for occ just use intercept
  if(transition) {
  alpha <- draws[,"beta_phi[1]"]
  }

  if(param_number == 1) {
    mr <- 0
  } else {
    mr <- 1
  }




  effect <- matrix(ncol = length(sr),
                   nrow = nrow(model_draws))
  c <- 1
  for(i in sr) {

    if(transition) {
      effect[,c] <- plogis(alpha[[1]] + model_draws[[1]] * mr + model_draws[[2]] * i) - plogis(alpha[[1]]) # scaled so no need to add additional mean vals
    } else if(abundance) {
      effect[,c] <- model_draws[[2]] * i
    } else {
      effect[,c] <- plogis(model_draws[[1]] * mr + model_draws[[2]] * i) # scaled so no need to add additional mean vals
    }
    c = c+1
  }

  colnames(effect) <- sr
  effect <- data.table::as.data.table(effect)
  effect <- data.table::melt(effect, id.vars = NULL, measure.vars = colnames(effect))
  if(!class(model_data[[model_column]]) == "factor") {
    if(log) {
      effect[, variable := scales::rescale(log(as.numeric(variable))^pwr, to = c(min(model_data[[model_column]]),                                                                max(model_data[[model_column]])))]
    } else {
      effect[, variable := scales::rescale(as.numeric(variable)^pwr, to = c(min(model_data[[model_column]]),                                                                max(model_data[[model_column]])))]
    }

  }
  effect$param <- model_column
  return(effect)
}

marginal_effects_plot_cmd <- function(data, col, factor = FALSE, trans = "identity", ylab = "Occupancy Probability") {
  data_sum <- data %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(median = median(value),
              q5 = quantile(value, 0.05),
              q25 = quantile(value, 0.25),
              q75 = quantile(value, 0.75),
              q95 = quantile(value, 0.95))

  if(factor) {
    plot <- data %>%
      ggplot2::ggplot(aes(x = variable, y = value)) +
      ggplot2::geom_violin(draw_quantiles = c(0.05, 0.5, 0.95), fill = col, alpha = 0.5) +
      ggplot2::ylab(ylab) +
      # ylim(c(0,47)) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.title.y = ggtext::element_markdown())
  } else{

    plot <- data_sum %>%
      ggplot2::ggplot(aes(x = variable, y = median)) +
      ggplot2::geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.25, fill = col) +
      ggplot2::geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.25, fill = col) +
      ggplot2::geom_line(linewidth = 1.5, colour = col) +
      ggplot2::ylab(ylab) +
      ggplot2::scale_x_continuous(trans = trans) +
      # ylim(c(0,47)) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.title.y = ggtext::element_markdown())
  }
  return(plot)
}

marginal_effects_plot_cmd_all <- function(data, col, factor = FALSE, trans = "identity", ylab = "Occupancy Probability") {
  data_sum <- data %>%
    dplyr::group_by(variable, species, param) %>%
    dplyr::summarise(median = median(value),
                     q5 = quantile(value, 0.025),
                     q25 = quantile(value, 0.25),
                     q75 = quantile(value, 0.75),
                     q95 = quantile(value, 0.975))

  if(factor) {
    plot <- data %>%
      ggplot2::ggplot(aes(x = variable, y = value)) +
      ggplot2::geom_violin(draw_quantiles = c(0.05, 0.5, 0.95), fill = col, alpha = 0.5) +
      ggplot2::ylab(ylab) +
      # ylim(c(0,47)) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.title.y = ggtext::element_markdown())
  } else{

    plot <- data_sum %>%
      ggplot2::ggplot(aes(x = variable, y = median)) +
      ggplot2::geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.25, fill = "grey20") +
      ggplot2::geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.25, fill = "grey20") +
      ggplot2::geom_line(linewidth = 1.5, colour = "grey20") +
      ggplot2::ylab(ylab) +
      ggplot2::scale_x_continuous(trans = trans) +
      # ylim(c(0,47)) +
      ggplot2::theme(axis.title.y = ggtext::element_markdown())
  }
  return(plot)
}
