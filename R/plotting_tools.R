#' @export
plot_sim_res <- function(dists, sims, acceptance_ratio = 0.01, viridis.option = "B",
                         real.sims = FALSE, real.alpha = 0.06, smooth = T, log_dist = T){
  # note accepted runs
  dists$hits <- 0
  if(acceptance_ratio > 1){
    acceptance_ratio <- acceptance_ratio/nrow(dists)
  }
  dists$hits <- ifelse(dists$dist <= quantile(dists$dist, acceptance_ratio), 1, 0)

  # prepare data
  ## real data if provided
  if(real.sims[1] != FALSE){
    real.sims$rID <- paste0(real.sims$iter, "_", real.sims$run)
    real.sims <- na.omit(real.sims)
  }

  ## sim data
  sims <- na.omit(sims)
  sims$dist <- dists$dist[match(sims$iter, dists$iter)]
  bh.sims <- sims[which(sims$iter %in% dists$iter[dists$hits == 1]),]
  bh.sims$rID <- paste0(bh.sims$iter, "_", bh.sims$run)

  # plot
  p <- ggplot2::ggplot()
  if(real.sims[1] != FALSE){
    p <- p + ggplot2::geom_line(data = real.sims, mapping = ggplot2::aes(x = gen, y = N, group = rID), alpha = real.alpha, color = "steelblue")
  }
  if(smooth){
    if(log_dist){
      p <- p + ggplot2::geom_smooth(data = bh.sims, mapping = ggplot2::aes(x = gen, y = N, group = iter, color = log10(dist)))
    }
    else{
      p <- p + ggplot2::geom_smooth(data = bh.sims, mapping = ggplot2::aes(x = gen, y = N, group = iter, color = dist))
    }
  }
  else{
    if(log_dist){
      p <- p + ggplot2::geom_line(data = bh.sims, mapping = ggplot2::aes(x = gen, y = N, group = rID, color = log10(dist)))
    }
    else{
      p <- p + ggplot2::geom_line(data = bh.sims, mapping = ggplot2::aes(x = gen, y = N, group = rID, color = dist))
    }
  }
  p <- p + ggplot2::scale_color_viridis_c(option = viridis.option) + ggplot2::theme_bw()

  return(p)
}

