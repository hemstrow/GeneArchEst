calc_distribution_stats <- function(x, meta, phenos, scheme = "gwas", chr = "chr",
                                    peak_delta = .5, peak_pcut = 0.0005, window_sigma = 50,
                                    burnin = NULL, thin = NULL, chain_length = NULL){
  #===========functions to run gwas or gp=============
  gwas <- function(x, phenos, meta, windows, peak_delta, peak_pcut, chr){
    x_pi <- pred(x, phenotypes = phenos,
                 prediction.program = "GMMAT",
                 maf.filt = F, runID = "gmmat_real")$e.eff$PVAL

    stats <- dist_desc(x_pi, meta, windows, peak_delta, peak_pcut, chr, pvals = T)
    return(stats = stats)
  }
  gp <- function(x,  phenos, meta, windows, peak_delta, peak_pcut, chr){
    cat("Beginning pseudo data", method, "run.\n")
    pseudo.pred <-pred(x, phenotypes = phenos,
                       burnin = burnin, thin = thin, chain_length = chain_length,
                       prediction.program = "BGLR", prediction.model = method,
                       runID = "gp_pseudo", verbose = F)

    stats <- dist_desc(pseudo.pred, meta, windows, peak_delta, peak_pcut, chr, pvals = F)

    return(stats = stats)
  }

  #=============get windows and run=========
  windows <- mark_windows(meta, window_sigma, chr)

  if(scheme == "gwas"){
    stats <- gwas(x, phenos, meta, windows, peak_delta, peak_pcut, chr)
  }
  else{
    stats <- gp(x, phenos, meta, windows, peak_delta, peak_pcut, chr)
  }

  return(stats)
}



#' Estimate pi from simulations.
#'
#' Estimate pi from distribution descriptions from simulated effect size distributions
#' using a random forest.
#'
#' @param x numeric. Descriptive statistics for the estimated effect sizes/associaion p-values
#'   from the observed data.
#' @param sims data.frame. Data.frame that matches that produced by \code{\link{sim_gen}}, containing descriptive statistics for the estimated effect sizes/associaion p-values
#'   from the simulted data as well as a column named 'pi' containing the simulated pi values. Each row should be a single
#'   simulation. May contain other columns titled 'df', 'scale', and 'h'.
#' @param p numeric < 1 and > 0, default .25. Proportion of sims to hold out from model estimation for use in cross-evalutation.
#' @param num.trees numeric, default 1000. Number of trees to grow during the random forest.
#' @param mtry numeric, default ncol(sims) - 1. Number of variables to possibly split at each node during random forest. See
#'   \code{\link[ranger]{ranger}} for details.
#' @param num.threads numeric, default NULL. Number of processing threads to use for tree growth and cross-evaluation.
#' @param pi_transform function or FALSE, default function(pi) log10(1 - pi). Transformation to use on pi. Usefull if pi values in simulations
#'   are heavily skewed, as is likely given that it is a non-linear parameter.
#' @param importance character, default "none". Determines how variable importance is computed, if it is at all. See \code{\link[ranger]{ranger}}
#'   for details. Note that "permutation", while accurate, can be very slow. By and large, this isn't needed for genetic architecture prediciton.
#' @param ... extra arguments passed to \code{\link[ranger]{ranger}}.
#'
#' @author William Hemstrom
#' @export
estimate_gen_arch_from_sims <- function(x, meta, phenos, sims, ABC_res, p = .25, num.trees = 1000, mtry = ncol(sims) - 1, num.threads = NULL,
                                        pi_transform = function(pi) log10(1 - pi), pi_back_transform = function(pi) 1 - 10^(pi),
                                        importance = "none",scheme = "gwas", chr = "chr",
                                        peak_delta = .5, peak_pcut = 0.0005, window_sigma = 50,
                                        burnin = NULL, thin = NULL, chain_length = NULL,
                                        ABC_acceptance_threshold = 0.005, ABC_dist_var = "ks",
                                        ABC_scale_transform = function(scale) log10(scale),
                                        ABC_scale_back_transform = function(scale) 10^scale, ...){

  #===========get stats for the real data==========
  cat("Getting descriptive statistics for the real data...\n")
  stats <- calc_distribution_stats(x, meta, phenos, scheme, chr,
                                   peak_delta, peak_pcut, window_sigma,
                                   burnin, thin, chain_length)
  cat("Done!\n")
  #===========prepare data=========
  # transform pi
  if(is.function(pi_transform)){
    cat("Pi transformed according to pi_transform. Please back-transform before estimating scale.\n")
    sims$pi <- pi_transform(sims$pi)
  }

  # remove extra columns from sims
  sims$df <- NULL
  sims$scale <- NULL
  sims$h <- NULL
  sims$iter <- NULL

  # check for NAs in the stats, if any then remove those collumns from the data
  na.stats <- which(is.na(stats))
  if(length(na.stats) > 0){
    stats <- x[-na.stats]
    sims <- sims[,-na.stats]
  }

  # remove nas
  sims <- na.omit(sims)

  # remove pi
  spi <- sims$pi
  sims$pi <- NULL

  #===========run the random forest=========
  cat("Generating a random forest model via ranger. This may take a while...\n")
  # grab training/test samples
  test_samps <- sample(nrow(sims), nrow(sims)*p, replace = F)
  dat_test <- sims[test_samps,]
  dat_train <- sims[-test_samps,]

  pi_test <- spi[test_samps]
  pi_train <- spi[-test_samps]
  dat_train$pi <- spi[-test_samps]

  # run ranger
  rf <- ranger::ranger(data = dat_train,
                       num.trees = num.trees, mtry = mtry,
                       dependent.variable.name = "pi",
                       num.threads = num.threads,
                       importance = importance,
                       keep.inbag = TRUE, ...)

  #===========do cross-evaluation=========
  cat("Done!\nPreforming cross-evaluation...\n")
  # get the predicted values and errors for each holdout
  pe <- forestError::quantForestError(forest = rf, X.train = dat_train[,-which(colnames(dat_train) == "pi")],
                                      X.test = dat_test,
                                      Y.train = dat_train$pi, n.cores = ifelse(is.null(num.threads), 1, num.threads))


  # add in real data and calculate errors
  pe$estimates$real <- pi_test
  pe$estimates$in_error <- ifelse(pe$estimates$real <= pe$estimates$upper_0.05 & pe$estimates$real >= pe$estimates$lower_0.05, 1, 0)
  pe$estimates$err <- pe$estimates$real - pe$estimates$pred
  pe$estimates$squared.err <- pe$estimates$err^2

  cat("Complete!\nr^2 of predictions vs true pi:", cor(pe$estimates$pred, pe$estimates$real)^2, ".\n")

  # make a diagnostic plot
  cv.plot <- ggplot2::ggplot(pe$estimates, ggplot2::aes(x = real, y = pred)) +
    ggplot2::theme_bw() + ggplot2::geom_point() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_0.05, ymax = upper_0.05, x = pred), fill = "grey", alpha = 0.5)
  cat("Creating diagnostic plot...\n")
  print(cv.plot)
  cat("Done!\n")

  #===========esitmate pi================
  cat("Estimating pi...\n")
  # fetch a pi value distribution for the real data
  pe2 <- forestError::quantForestError(forest = rf, X.train = dat_train[,-which(colnames(dat_train) == "pi")],
                                       X.test = rbind.data.frame(stats, dat_test),
                                       Y.train = dat_train$pi, n.cores = num.threads)

  # get the distribution of possible pi values around the point estimate.
  quantiles <- seq(0 + .0001, 1 - .0001, by = .0001)
  qs <- pe2$qerror(quantiles, 1) # quantiles for errors from 0.0001 to 0.9999 (. 1% to 99.99%)
  ed <- data.frame(q = as.numeric(qs), val = quantiles)
  ed$q <- ed$q + pe2$estimates$pred[1]
  ed$val[ed$val > .5] <- 1 - ed$val[ed$val > .5] # ed now contains: q: the pi values, val: the probabilities of those pi values.
  colnames(ed) <- c("pi", "prob")

  # back transform
  ed$pi <- pi_back_transform(ed$pi)

  cat("Done!\n")

  #===========estimate scale=============
  cat("Estimating scale...\n")
  scale <- estimate_scale_from_ABC(ABC_res, ed$pi, pi_transform = pi_transform,
                                   threshold = ABC_acceptance_threshold,
                                   dist_var = ABC_dist_var, quantiles = quantiles, scale_transform = ABC_scale_transform,
                                   scale_back_transform = ABC_scale_back_transform)

  browser()
  #==========return======================
  return(list(forest = rf, cross_validation = pe, pi_point = pe2$estimates[1,], pi_density = ed,
              descriptive_stats = stats))
}


#' Estimate scale ABC results
#'
#' Esitmates scale based on ABC results using gam smoothing.
#'
#' @param x data.frame. ABC results, such as those given by \code{\link{ABC_on_hyperparameters}}.
#' @param pi numeric. Vector of pi values for which to calculate scale.
#' @param pi_transform function or FALSE, default function(pi) log10(1 - pi). Transformation to use on pi. Usefull if pi values in simulations
#'   are heavily skewed, as is likely given that it is a non-linear parameter.
#' @param threshold numeric, default 0.005. Proportion of accepted runs.
#' @param dist_var character, default "ks". Name of the distance varaible to use in picking accepted runs.
#'
#' @export
estimate_scale_from_ABC <- function(x, pi, pi_transform = function(pi) log10(1 - pi), threshold = .005, dist_var = "ks",
                                    quantiles = seq(0 + .0001, 1 - .0001, by = .0001), scale_transform = function(scale) log10(scale),
                                    scale_back_transform = function(scale) 10^scale){
  #==================grab the accepted runs============
  hits <- which(x[,dist_var] <= quantile(x[,dist_var], threshold))

  #==================transform and esitmate===========
  # transform
  if(is.function(pi_transform)){
    x$pi <- pi_transform(x$pi)
    pi <- pi_transform(pi)
  }
  x$scale <- scale_transform(x$scale)
  y <- x[hits,]
  hold<- sample(nrow(y), nrow(y)*.25)
  y_train <- y[-hold,]
  y_test <- y[hold,]


  # smooth scale vs pi on the accepted runs and predict
  es <- data.frame(pi = pi) # for prediction
  gs_test <- mgcv::gam(scale ~ s(pi), data = y_train)
  gs <- mgcv::gam(scale ~ s(pi), data = y) # the model

  # get confidence limits
  get_conf_int <- function(mod, pis, quantiles){
    # get standard errors
    Designmat <- predict(mod, newdata=data.frame(pi = pis), type="lpmatrix")
    pr <- predict(mod, data.frame(pi = pis))
    diagMult <- compiler::cmpfun(function(m1, m2) sapply(seq_len(nrow(m1)), function(i) m1[i,] %*% m2[,i]))
    tmp <- Designmat %*% vcov(mod)
    predvar <-  diagMult(tmp, t(Designmat))
    SE <- sqrt(predvar)
    SE2 <- sqrt(predvar+mod$sig2)

    # get intervals for each quantile
    tfrach <- qt(quantiles, mod$df.residual)
    intervals <- sapply(tfrach, function(z) z * SE2)
    intervals <- matrix(intervals, nrow = length(SE2), ncol = length(tfrach)) # a row for each pi value, a column for each quantile
    colnames(intervals) <- quantiles
    intervals <- cbind(pi = pis, pr_scale = pr, intervals)
    return(as.data.frame(intervals))
  }

  # with test data for plot
  intervals <- get_conf_int(gs_test, y_test$pi, c(.025, .975))
  intervals <- cbind(scale = y_test$scale, intervals)
  colnames(intervals)[4:5] <- c("clow", "chigh")
  scale_conf_plot <- ggplot2::ggplot(intervals, ggplot2::aes(x = pi, y = scale)) + ggplot2::geom_point() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pr_scale + clow, ymax = pr_scale + chigh), alpha = 0.5) + ggplot2::theme_bw()


  cat("Plotting 95% confidence intervals for test values for fit evaluation of scale.\n")
  print(scale_conf_plot)

  # with real data
  intervals <- get_conf_int(gs, pi, quantiles)

  pr_scales <- intervals$pr_scale + intervals[,-c(1:2)]

  #=========back transform and return===========
  if(is.function(scale_transform)){
    pr_scales <- scale_back_transform(pr_scales)
    intervals$pr_scale <- scale_back_transform(intervals$pr_scale)
  }

  return(list(scale_quantiles = pr_scales, opt_scale = intervals$pr_scale, fit_plot = scale_conf_plot))
}
