#' @export
plot_rf <- function(rf, hyperparameter = "pi"){
  #=========cross val plot==========
  cvp <- ggplot2::ggplot(rf[[hyperparameter]]$cross_validation,
                         ggplot2::aes(x = real, y = pred)) +
    ggplot2::geom_point(ggplot2::aes(color = as.logical(in_error))) +
    ggplot2::theme_bw() + ggplot2::labs(color = "In 95% C.I.")
  #========density=================
  pretty_breaks <- pretty(rf[[hyperparameter]]$parameter_density$pi, n = 3)
  pretty_breaks <- sort(unlist(c(pretty_breaks, rf[[hyperparameter]]$point_estimate[c(1, 4, 5)])))
  names(pretty_breaks) <- NULL
  hdp <- ggplot2::ggplot(rf[[hyperparameter]]$parameter_density,
                         ggplot2::aes_string(x = hyperparameter, y = "prob")) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(ggplot2::aes(xintercept = pred), color = "red", data = rf[[hyperparameter]]$point_estimate) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = lower_0.05), color = "red",
                        linetype = "dashed", data = rf[[hyperparameter]]$point_estimate) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = upper_0.05), color = "red",
                        linetype = "dashed", data = rf[[hyperparameter]]$point_estimate) +
    ggplot2::theme_bw() + ggplot2::ylab("Density")

  #========return=================
  return(list(cross_val = cvp, density = hdp))
}

#' @export
plot_reg <- function(reg){
  #=============cross validaiton path===========
  dep_conf_plot <- ggplot2::ggplot(reg$cross_val_data, ggplot2::aes(x = indep, y = dep)) + ggplot2::geom_point() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pred + clow, ymax = pred + chigh), alpha = 0.5) + ggplot2::theme_bw() +
    ggplot2::xlab(reg$independent[1]) + ggplot2::ylab(reg$dependent)

  #=============Joint quantile plot===============
  string <- "Plotting esitmates. Red polygon represents 95% prediction limits of independent variable along y axis, and 95% prediction limits for each of those independent varibale values along the y axis."
  tp <- ggplot2::ggplot(data = reg$joint_quantile_long, ggplot2::aes(x = independent_1, y = dependent, z = norm_joint_quantile)) + ggplot2::stat_summary_hex(bins = 100) +
    ggplot2::scale_fill_viridis_c() + ggplot2::scale_color_viridis_c() + ggplot2::theme_bw() +
    ggplot2::geom_polygon(data = reg$ci_path, ggplot2::aes(x = independent_1, y = dependent), color = "red", fill = NA, size = 1) +
    ggplot2::xlab(reg$independent[1]) + ggplot2::ylab(reg$dependent) + ggplot2::labs(fill = "Average joint probability")

  return(list(cross_val_plot = dep_conf_plot, joint_quantile_plot = tp, msg = string))
}
