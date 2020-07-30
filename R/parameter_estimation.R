calc_distribution_stats <- function(x, meta, phenos, center = T, scheme = "gwas", chr = "chr",
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
  if(center){
    phenos <- phenos - mean(phenos)
  }

  windows <- mark_windows(meta, window_sigma, chr)

  if(scheme == "gwas"){
    stats <- gwas(x, phenos, meta, windows, peak_delta, peak_pcut, chr)
  }
  else{
    stats <- gp(x, phenos, meta, windows, peak_delta, peak_pcut, chr)
  }

  return(stats)
}



#' Estimate hyperparameters from simulations using a random forest.
#'
#' Estimate one or more hyperparameters using descriptions of GWAS p-values distributions from
#' simulated effect size distributions, such as those produced by \code{\link{sim_gen}}
#' using a random forest.
#'
#' @param x numeric matrix. Input genotypes, SNPs as rows, columns as individuals. Genotypes formatted as 0,1,2 for the major homozygote, heterozygote, and minor homozygote, respectively.
#' @param meta data.frame. Metadata for SNPs. First two columns must hold chromosome ID and position. Futher columns ignored.
#' @param phenos numeric vector. Observed phenotypes, one per individual.
#' @param sims data.frame. Data.frame that matches that produced by \code{\link{sim_gen}}, containing descriptive statistics for the estimated effect sizes/associaion p-values
#'   from the simulted data as well as other columns containing the hyperparameters of interest for those simulations. Each row should be a single
#'   simulation. Columns not matching those expected and produced via \code{\link{sim_gen}} will be ignored.
#' @param hyperparameter_to_estimate character vector, default "pi". Names of the hyperparameters to estimate via random forest. Must match column names in sims.
#' @param center logical, default T. Determines if the phenotypes provided
#'   should be centered (have their means set to 0).
#'   This should match what was provided to \code{\link{sim_gen}}, as it does given the defaults for both functions.
#' @param hold_percent numeric < 1 and > 0, default .25. Proportion of sims to hold out from model estimation for use in cross-evalutation.
#' @param num_trees numeric, default 1000. Number of trees to grow during the random forest.
#' @param mtry function, default function(columns) columns. A function that, when given the number of columns containing distribution summary statistics,
#'   returns the number of variables to possibly split at each node during random forest. For example, function(columns) columns/2 would have an mtry
#'   equal to half the number of summary statistics. See \code{\link[ranger]{ranger}} for details about mtry.
#' @param num_threads numeric, default NULL. Number of processing threads to use for tree growth and cross-evaluation.
#' @param parameter_transforms. Named list of parameter transformations or NULL, default list(pi = function(pi) log10(1 - pi)).
#'   Transformations to use on any estimated parameters. Any estimated hyperparameters with names matching those
#'   in this list will be transformed as given. Usefull if pi or other hyperparameter values in simulations
#'   are heavily skewed, as is often likely.
#' @param parameter_back_transforms. Named list of parameter back transformations or NULL, default ist(pi = function(pi) 1 - 10^pi ).
#'   Back transformations to use on any estimated parameters. Any estimated hyperparameters with names matching those
#'   in this list will be back_transformed as given prior to being returned.
#' @param importance character, default "permutation". Determines how variable importance is computed, if it is at all. See \code{\link[ranger]{ranger}}
#'   for details. Note that "permutation", while accurate and thus the default, can be very slow.
#'   By and large, this isn't needed for genetic architecture prediciton, and can be set to "none" if not wanted.
#' @param peak_delta numeric, default 0.5. Value used to determine spacing between called peaks during peak identification for distribution description.
#' @param peak_pcut numeric, default 0.0005. Only p-values below this quantile will be used for peak detection during peak indentification for distribution description.
#' @param window_sigma numeric, default = 50. Size of the windows in megabases to be used during distribution description.
#' @param quantiles numeric, default seq(0 + 0.001, 1 - 0.001, .001). Density quantiles over which to estimate parameter values.
#' @param save_rf logical, default FALSE. If true, the raw ranger random forest object is returned. Can be extremely large, and
#'   not needed unless different quantiles/predictions/etc are needed.
#' @param ... Extra arguments passed to \code{\link[ranger]{ranger}}.
#'
#' @author William Hemstrom
#' @export
hyperparameter_random_forest <- function(x, meta, phenos, sims, hyperparameter_to_estimate = c("pi"),
                                         center = T, hold_percent = .25, num_trees = 1000,
                                         mtry = function(columns) columns, num_threads = NULL,
                                         parameter_transforms = list(pi = function(pi) log10(1 - pi)),
                                         parameter_back_transforms = list(pi = function(pi) 1 - 10^pi ),
                                         importance = "permutation", scheme = "gwas",
                                         peak_delta = .5, peak_pcut = 0.0005, window_sigma = 50,
                                         quantiles = seq(0 + .001, 1 - .001, by = .001), save_rf = FALSE,
                                         ...){
  #==========rf construction, evaluation, and prediction subfunction==========
  make_and_predict_rf <- function(dat_test, dat_train, param_train, param_test, hyperparameter,
                                  mtry, num_trees, num_threads, importance, stats, save_rf = FALSE, ...){


    #=================construct the random forest=================
    cat("Constructing random forest for:", hyperparameter, "\n")
    dat_train$pred_var <- param_train[[hyperparameter]]
    colnames(dat_train)[ncol(dat_train)] <- hyperparameter

    rf <- ranger::ranger(data = dat_train,
                         num.trees = num_trees,
                         mtry = mtry(ncol(dat_train) - 1),
                         dependent.variable.name = hyperparameter,
                         num.threads = num_threads,
                         importance = importance,
                         keep.inbag = TRUE, ...)



    #===========do cross-evaluation=========
    cat("Done!\nPreforming cross-evaluation...\n")
    # get the predicted values and errors for each holdout
    pe <- forestError::quantForestError(forest = rf, X.train = dat_train[,-which(colnames(dat_train) == hyperparameter)],
                                        X.test = dat_test,
                                        Y.train = dat_train[,hyperparameter],
                                        n.cores = ifelse(is.null(num_threads), 1, num_threads))

    # add in real data and calculate errors
    pe$estimates$real <- param_test[,hyperparameter]
    pe$estimates$in_error <- ifelse(pe$estimates$real <= pe$estimates$upper_0.05 & pe$estimates$real >= pe$estimates$lower_0.05, 1, 0)
    pe$estimates$err <- pe$estimates$real - pe$estimates$pred
    pe$estimates$squared.err <- pe$estimates$err^2

    cat("Complete!\nr^2 of predictions vs true", paste0(hyperparameter, ":"), cor(pe$estimates$pred, pe$estimates$real)^2, "\n")

    # make a diagnostic plot
    cv.plot <- ggplot2::ggplot(pe$estimates, ggplot2::aes(x = real, y = pred)) +
      ggplot2::theme_bw() + ggplot2::geom_point() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_0.05, ymax = upper_0.05, x = pred), fill = "grey", alpha = 0.5)
    cat("Creating diagnostic plot...\n")
    print(cv.plot)
    cat("Done!\n")

    #===========esitmate pi================
    cat("Estimating", hyperparameter, "...\n")
    # fetch a pi value distribution for the real data
    pe2 <- forestError::quantForestError(forest = rf, X.train = dat_train[,-which(colnames(dat_train) == hyperparameter)],
                                         X.test = rbind.data.frame(stats, dat_test),
                                         Y.train = dat_train[,hyperparameter], n.cores = ifelse(is.null(num_threads), 1, num_threads))

    # get the distribution of possible pi values around the point estimate.
    qs <- pe2$qerror(quantiles, 1) # quantiles for errors from 0.0001 to 0.9999 (. 1% to 99.99%)
    ed <- data.frame(q = as.numeric(qs), val = quantiles)
    ed$q <- ed$q + pe2$estimates$pred[1]
    ed$val[ed$val > .5] <- 1 - ed$val[ed$val > .5] # ed now contains: q: the pi values, val: the probabilities of those pi values.
    colnames(ed) <- c(hyperparameter, "prob")

    cat("Done!\n")


    return(list(rf = ifelse(save_rf, rf, FALSE), cross_validation = pe$estimates, point_estimate = pe2$estimates[1,], parameter_density = ed,
                descriptive_stats = stats,
                cross_val_plot = cv.plot))
  }





  #===========sanity checks========================
  msg <- character()
  # sim data
  ## colnames
  sim_cols <- colnames(sims)
  meta_sim_cols <- which(!colnames(sims) %in% names_descriptive_stats)

  # no predict variable data
  if(!all(hyperparameter_to_estimate %in% colnames(sims))){
    missing_hypers <- which(!hyperparameter_to_estimate %in% colnames(sims))
    msg <- c(msg, paste0("Some hyperparameters not found in simulated data:\n\t", paste0(missing_hypers, collapse = "\n\t")))
  }

  if(length(msg) > 0){
    stop(paste0(msg, collapse = "\n"))
  }

  sims <- as.data.frame(sims)
  meta <- as.data.frame(meta)

  #===========get stats for the real data==========
  cat("Getting descriptive statistics for the real data...\n")
  stats <- calc_distribution_stats(x, meta, phenos, center, scheme, colnames(meta)[1],
                                   peak_delta, peak_pcut, window_sigma)
  cat("Done!\n")
  #===========prepare data=========
  # pull parameters to estimate and transform if requested
  predict_params <- sims[, hyperparameter_to_estimate, drop = F]
  if(length(parameter_transforms) > 0){
    for(i in 1:length(parameter_transforms)){
      col.match <- which(colnames(predict_params) == hyperparameter_to_estimate[i])
      predict_params[,col.match] <- parameter_transforms[[i]](predict_params[,col.match])
    }
  }

  # remove extra columns from sims
  sims <- sims[,-meta_sim_cols]

  # check for NAs in the stats, if any then remove those collumns from the data
  na.stats <- which(is.na(stats))
  if(length(na.stats) > 0){
    stats <- x[-na.stats]
    sims <- sims[,-na.stats]
    predict_params[,-na.stats]
  }

  # remove nas
  sims <- na.omit(sims)
  predict_params <- predict_params[-attributes(sims)$na.action,, drop = F]


  #===========run the random forest prediction function=========
  cat("Generating  random forest models for each esitmated hyperparameter via ranger. This may take a while...\n")
  # grab training/test samples
  test_samps <- sample(nrow(sims), nrow(sims)*hold_percent, replace = F)
  dat_test <- sims[test_samps,]
  dat_train <- sims[-test_samps,]

  param_test <- predict_params[test_samps,,drop = F]
  param_train <- predict_params[-test_samps,, drop = F]

  rfl <- vector("list", length(hyperparameter_to_estimate))
  names(rfl) <- hyperparameter_to_estimate

  # run the rfs
  for(i in 1:length(rfl)){
    rfl[[i]] <- make_and_predict_rf(dat_test = dat_test, dat_train = dat_train,
                                    param_train = param_train, param_test = param_test,
                                    hyperparameter = hyperparameter_to_estimate[i],
                                    mtry = mtry,
                                    num_trees = num_trees,
                                    num_threads = num_threads,
                                    importance = importance,
                                    stats = stats,
                                    save_rf = save_rf,
                                    ...)
  }


  #==========return======================
  if(length(parameter_back_transforms) > 0){
    for(i in 1:length(parameter_transforms)){
      # point estimate:
      match <- rfl[[names(parameter_back_transforms)[i]]]
      point_est <- parameter_back_transforms[[i]](match$point_estimate)
      point_est_fix <- point_est[c(1,2,3,5,4)]
      point_est_fix$mspe <- abs(point_est_fix$mspe)
      names(point_est_fix) <- names(point_est)
      match$point_estimate <- point_est_fix

      # densities:
      match$parameter_density[1] <- parameter_back_transforms[[i]](match$parameter_density[1])

      # return
      rfl[[names(parameter_transforms)[i]]] <- match
    }
  }

  return(rfl)
}


#' Estimate a hyperparamter from ABC results via regression.
#'
#' Esitmates scale based on ABC results using gam smoothing or regression forest.
#' @export
hyperparameter_regression_on_ABC <- function(ABC, input_independent_parameters, formula = scale ~ pi,
                                             regression_method = "rf",
                                             num_trees = 10000,
                                             num_threads = NULL,
                                             parameter_transforms = list(pi = function(pi) log10(1 - pi),
                                                                         scale = function(scale) log10(scale)),
                                             parameter_back_transforms = list(pi = function(pi) 1 - 10^pi,
                                                                              scale = function(scale) 10^scale),
                                             acceptance_threshold = .005, dist_var = "ks",
                                             hold_percent = .25,
                                             quantiles = seq(0 + .001, 1 - .001, by = .001),
                                             independent_quantiles = seq(0 + .001, 1 - .001, by = .001)){

  #==================grab the accepted runs============
  hits <- which(ABC[,dist_var] <= quantile(ABC[,dist_var], acceptance_threshold, na.rm = T))
  ABC <- ABC[hits,]


  #==================fix input explanitories and transform any requested variables===========
  # grab the independent and dependent variables
  independent <- all.vars(formula)[-1]
  dependent <- all.vars(formula)[1]


  # fix input independent parameters if not passed as a data.frame
  if(!is.data.frame(input_independent_parameters)){
    input_independent_parameters <- data.frame(input_independent_parameters)
    colnames(input_independent_parameters) <- independent
  }

  # transform
  if(length(parameter_transforms) > 0){
    for(i in 1:length(parameter_transforms)){
      tparm <- names(parameter_transforms)[i]
      ABC[,tparm] <- parameter_transforms[[i]](ABC[,tparm])


      if(tparm %in% colnames(input_independent_parameters)){
        input_independent_parameters[,tparm] <-
          parameter_transforms[[i]](input_independent_parameters[,tparm])
      }
    }
  }

  #=================fit the regression===========================
  hold<- sample(nrow(ABC), nrow(ABC)*hold_percent)
  y_train <- ABC[-hold,]
  y_test <- ABC[hold,]

  # smooth scale vs pi on the accepted runs and predict
  if(regression_method == "gam"){
    gs_test <- mgcv::gam(formula, data = y_train)
    gs <- mgcv::gam(formula, data = ABC) # the model

    plot(gs)
    points(ABC[,independent[1]], ABC[,dependent])
  }
  else if (regression_method == "lm"){
    gs_test <- lm(formula, data = y_train)
    gs <- lm(formula, data = ABC)
  }
  else if (regression_method == "rf"){
    gs_test <- ranger::ranger(data = y_train[,c(independent, dependent)],
                              num.trees = num_trees,
                              dependent.variable.name = dependent,
                              num.threads = num_threads,
                              keep.inbag = TRUE)
    gs <- ranger::ranger(data = ABC[,c(independent, dependent)],
                         num.trees = num_trees,
                         num.threads = num_threads,
                         dependent.variable.name = dependent,
                         keep.inbag = TRUE)
  }



  # get confidence limits
  get_conf_int <- function(mod, independent_parms, quantiles, round = "test"){

    # get standard errors
    if(regression_method == "gam"){
      Designmat <- predict(mod, newdata = independent_parms, type="lpmatrix")
      pr <- predict(mod, independent_parms)
      diagMult <- compiler::cmpfun(function(m1, m2) sapply(seq_len(nrow(m1)), function(i) m1[i,] %*% m2[,i]))
      tmp <- Designmat %*% vcov(mod)
      predvar <-  diagMult(tmp, t(Designmat))
      SE <- sqrt(predvar)
      SE2 <- sqrt(predvar+mod$sig2)

      # get intervals for each quantile
      tfrach <- qt(quantiles, mod$df.residual)
      intervals <- sapply(tfrach, function(z) z * SE2)
      intervals <- matrix(intervals, nrow = length(SE2), ncol = length(tfrach)) # a row for each parameter value, a column for each quantile
      colnames(intervals) <- quantiles

      pr <- as.data.frame(pr)
      colnames(pr) <- paste0("pred_", dependent)
      intervals <- cbind(independent_parms, pr, intervals)
      return(as.data.frame(intervals))
    }
    else if(regression_method == "rf"){
      if(round == "test"){
        pe <- forestError::quantForestError(forest = mod, X.train = y_train[,independent, drop = F],
                                            X.test = independent_parms,
                                            Y.train = y_train[,dependent], n.cores = ifelse(is.null(num_threads), 1, num_threads))
      }
      else if(round == "prediction"){
        pe <- forestError::quantForestError(forest = mod, X.train = ABC[,independent, drop = F],
                                            X.test = independent_parms,
                                            Y.train = ABC[,dependent], n.cores = ifelse(is.null(num_threads), 1, num_threads))

      }
      intervals <- pe$qerror(quantiles) # quantiles for errors from 0.0001 to 0.9999 (. 1% to 99.99%)
      colnames(intervals) <- quantiles
      intervals <- cbind(pe$estimates$pred, intervals)
      colnames(intervals)[1] <- paste0("pred_", dependent)
      intervals <- cbind(independent_parms, intervals)
    }
  }

  # with test data for plot
  intervals <- get_conf_int(gs_test, y_test[,independent, drop = F], c(.025, .975), round = "test")
  intervals <- cbind(y_test[,dependent], intervals)
  colnames(intervals)[1] <- dependent
  colnames(intervals)[4:5] <- c("clow", "chigh")
  pdat <- intervals[,c(dependent, independent[1], "clow", "chigh", paste0("pred_", dependent))]
  colnames(pdat)[1:2] <- c("dep", "indep")
  colnames(pdat)[5] <- "pred"
  dep_conf_plot <- ggplot2::ggplot(pdat, ggplot2::aes(x = indep, y = dep)) + ggplot2::geom_point() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pred + clow, ymax = pred + chigh), alpha = 0.5) + ggplot2::theme_bw() +
    ggplot2::xlab(independent[1]) + ggplot2::ylab(dependent)


  cat("Plotting 95% confidence intervals for test values for fit evaluation of", dependent, ".\n")
  print(dep_conf_plot)

  # with real data
  intervals <- get_conf_int(gs, input_independent_parameters, quantiles, round = "prediction")

  predictions <- intervals[,paste0("pred_", dependent)] + intervals[,-c(1:(1 + length(independent)))]

  #=========back transform and return===========
  if(length(parameter_back_transforms) > 1){

    for(i in 1:length(parameter_back_transforms)){
      if(dependent == names(parameter_back_transforms)[i]){
        predictions <- parameter_back_transforms[[i]](predictions)
        match_cols <- which(!colnames(intervals) %in% independent) # need to back transform both the point and quantile predictions (everything except the independent vars)
        intervals[,match_cols] <- parameter_back_transforms[[i]](intervals[,match_cols])
      }
      else{
        intervals[,names(parameter_back_transforms)[i]] <-
          parameter_back_transforms[[i]](intervals[,names(parameter_back_transforms)[i]])
      }
    }
  }

  #===========clean, plot, and return====================
  # working point
  cat("Cleaning and preparing summary plot...\n")
  # get the probs for the quantiles
  probs <- independent_quantiles
  probs[probs > 0.5] <- 1 - probs[probs > 0.5]

  # add pi data and melt the scale data
  combined_predictions <- cbind(intervals[,which(colnames(intervals) %in% independent), drop = F],
                                independent_quantiles = probs,
                                as.data.table(predictions))

  qm <- reshape2::melt(combined_predictions, id.vars = c(independent, "independent_quantiles"))
  qm$variable <- as.numeric(as.character(qm$variable))
  colnames(qm)[which(colnames(qm) %in% c("variable", "value"))] <- c("dependent_quantile", "dependent")
  colnames(qm)[which(colnames(qm) %in% independent)] <- paste0("independent_", 1:length(independent))

  # fix scale quantiles, get joint quantile, and plot
  probs <- qm[,"dependent_quantile"]
  probs[probs > 0.5] <- 1 - probs[probs > 0.5]
  qm[,"dependent_quantile"] <- probs
  qm$joint_quantile <- qm$independent_quantiles * probs
  qm$norm_joint_quantile <- qm$joint_quantile/sum(qm$joint_quantile)
  ## prepareconfidence limits for geom_polygon
  upper_quant_name <- colnames(predictions)[which.min(abs(as.numeric(colnames(predictions)) - 0.975))] # find the quantiles closest to 97.5% and 2.5% (for 95% CIs. These will usually be those CIs, but might vary a bit depending on the quantiles set)
  lower_quant_name <- colnames(predictions)[which.min(abs(as.numeric(colnames(predictions)) - 0.025))]
  if(upper_quant_name != "0.975" | lower_quant_name != "0.025"){
    warning("Since the 97.5 and 2.5 quantiles are not in the requested quantiles, upper and lower CIs in plot instead reflect:\n\tupper:",
            as.numeric(upper_quant_name)*100, "\n\tlower:",
            as.numeric(lower_quant_name)*100)
  }
  dependent_upper_95 <- predictions[,upper_quant_name]
  dependent_lower_95 <- predictions[,lower_quant_name]
  dependent_cis <- data.frame(independent_1 = intervals[,independent[1]], upper = dependent_upper_95, lower = dependent_lower_95)
  path <- dependent_cis[which(dependent_cis$independent %in% intervals[,independent[1]][which(quantiles >= 0.025 & quantiles <= 0.975)]),]
  path <- reshape2::melt(path, id.vars = "independent_1")
  colnames(path)[3] <- "dependent"
  path$order <- 1:nrow(path)
  path <- path[c(which(path$variable == "lower"), rev(path$order[path$variable == "upper"])),]
  path$norm_joint_quantile <- NA
  ## plot
  tp <- ggplot2::ggplot(data = qm, ggplot2::aes(x = independent_1, y = dependent, z = norm_joint_quantile)) + ggplot2::stat_summary_hex(bins = 100) +
    ggplot2::scale_fill_viridis_c() + ggplot2::scale_color_viridis_c() + ggplot2::theme_bw() +
    ggplot2::geom_polygon(data = path, ggplot2::aes(x = independent_1, y = dependent), color = "red", fill = NA, size = 1) +
    ggplot2::xlab(independent[1]) + ggplot2::ylab(dependent) + ggplot2::labs(fill = "Average joint probability")

  cat("Plotting esitmates. Red polygon represents 95% prediction limits of independent variable along y axis, and 95% prediction limits for each of those independent varibale values along the y axis.\n")
  print(tp)

  ret_quants <- cbind(intervals[,independent], predictions)
  colnames(ret_quants)[1:length(independent)] <- independent
  return(list(optimal_fits = intervals[,c(independent, paste0("pred_", dependent))],
              quantiles = ret_quants,
              cross_val_fit_plot = dep_conf_plot,
              joint_quantile_plot = tp))
}
