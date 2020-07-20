#' @export
sim_gen <- function(x, meta, iters, center = T, scheme = "gwas",
                    effect_distribution = rbayesB,
                    parameter_distributions = list(pi = function(x) rbeta(x, 25, 1),
                                                   d.f = function(x) runif(x, 1, 100),
                                                   scale = function(x) rbeta(x, 1, 3)*100),
                    h_dist = function(x) rep(.5, x),
                    # burnin = burnin, thin = thin, chain_length = chain_length, method = "BayesB",
                    par = F, joint_res = NULL, joint_acceptance = NULL, joint_res_dist = "ks",
                    peak_delta = .5, peak_pcut = 0.0005, window_sigma = 50){

  #============schem functions for one simulation=============
  # gp <- function(x, pi, df, scale, method, t_iter, h, windows, center = center){
  #   pseudo <- generate_pseudo_effects(x, effect_distribution, parameters, h, center = center)
  #
  #   cat("Beginning pseudo data", method, "run.\n")
  #   pseudo.pred <-pred(x, phenotypes = pseudo$p,
  #                      burnin = burnin, thin = thin, chain_length = chain_length,
  #                      prediction.program = "BGLR", prediction.model = method,
  #                      runID = paste0(t_iter, "_pseudo"), verbose = F)
  #
  #   stats <- dist_desc(pseudo.pred, meta, windows, peak_delta, peak_pcut, chr, pvals = F)
  #
  #   return(list(stats = stats, e = pseudo$e))
  # }
  gwas <- function(x, effect_distribution, parameters, h, center = center,
                   t_iter, G, windows){

    pseudo <- generate_pseudo_effects(x, effect_distribution, parameters, h, center = center)

    pseudo_pi <- pred(x, phenotypes = pseudo$p,
                      prediction.program = "GMMAT",
                      maf.filt = F, runID = paste0(t_iter, "_gmmat"),
                      pass_G = G)$e.eff$PVAL

    stats <- dist_desc(pseudo_pi, meta, windows, peak_delta, peak_pcut, colnames(meta)[1], pvals = T)
    return(stats)
  }

  loop_func <- function(x, effect_distribution, parameters, scheme,
                        t_iter, G = NULL, h, windows, center = center){
    # if(scheme == "gp"){
    #   dist <- gp(x, pi, df, scale, method, t_iter, h, windows, center = center)
    # }
    if(scheme == "gwas"){
      dist <- gwas(x, effect_distribution, parameters = parameters, h = h, center = center,
                   t_iter = t_iter, windows = windows, G = G)
    }
    return(dist)
  }

  #============prep for simulations======================================
  # if any joint parameter priors, calculate and disambiguate
  joint_parms <- names(parameter_distributions)[which(parameter_distributions == "joint")]
  if(length(joint_parms) > 0){
    parms <- gen_parms(iters, joint_res, joint_acceptance, joint_parms, dist.var = joint_res_dist)

    if(ncol(parms) < length(parameter_distributions)){
      other_parms <- names(parameter_distributions)[which(!names(parameter_distributions) %in% colnames(parms))]
      run_parameters <- vector("list", length = length(other_parms))
      names(run_parameters) <- other_parms
      for(i in 1:length(run_parameters)){
        run_parameters[[i]] <- parameter_distributions[[other_parms[i]]](iters)
      }
      parms <- cbind(parms, as.data.frame(run_parameters))
    }
    run_parameters <- parms
    rm(parms)
  }
  h <- h_dist(iters)


  # can pass a g matrix forward once if doing gwas
  if(scheme == "gwas"){
    ind.genos <- convert_2_to_1_column(x)
    colnames(ind.genos) <- paste0("m", 1:ncol(ind.genos)) # marker names
    rownames(ind.genos) <- paste0("s", 1:nrow(ind.genos)) # ind IDS

    G <- AGHmatrix::Gmatrix(ind.genos, missingValue = -18, method = "Yang", maf = 0.05)
    colnames(G) <- rownames(ind.genos)
    rownames(G) <- rownames(ind.genos)
  }
  else{
    G <- NULL
  }

  # pre-run the window function
  windows <- mark_windows(meta, window_sigma, colnames(meta)[1])


  #============run the simulations===========================
  # initialize storage
  dist_output <-  matrix(0, iters, ncol = number_descriptive_stats)
  colnames(dist_output) <- names_descriptive_stats

  # run the simulations
  ## serial
  if(par == F){
    for(i in 1:iters){
      cat("Iter: ", i, ".\n")
      dist_output[i,] <- loop_func(x, effect_distribution,
                                   parameters = as.list(run_parameters[i,,drop = F]),
                                   scheme = scheme,
                                   t_iter = i, G = G, h = h[i], windows = windows, center = center)
    }
    ret <- cbind(as.data.table(run_parameters), h = h, as.data.table(dist_output))
    if(!is.list(ret)){
      ret <- as.data.frame(t(ret))
    }

    return(list(stats = ret, errors = list(parms = ret[-c(1:nrow(ret)),], msgs = character())))
  }


  # parallel
  else{
    cl <- snow::makeSOCKcluster(par)
    doSNOW::registerDoSNOW(cl)

    # divide up into ncore chunks
    chunks <- list(parms = split(run_parameters, (1:iters)%%par),
                   dist_output = split(as.data.frame(dist_output), (1:iters)%%par),
                   h = split(as.data.frame(h), (1:iters)%%par))


    # prepare reporting function
    progress <- function(n) cat(sprintf("Chunk %d out of", n), par, "is complete.\n")
    opts <- list(progress=progress)

    output <- foreach::foreach(q = 1:par, .inorder = FALSE,
                               .options.snow = opts, .packages = c("data.table", "GeneArchEst")
                               ) %dopar% {

                                 parm_chunk <- chunks$parm[[q]]
                                 h_chunk <- unlist(chunks$h[[q]])
                                 dist_chunk <- chunks$dist_output[[q]]
                                 is.err <- numeric(0)
                                 errs <- character(0)


                                 # run once per iter in this chunk
                                 for(i in 1:nrow(parm_chunk)){
                                   b <- try(loop_func(x, effect_distribution,
                                                      parameters = as.list(parm_chunk[i,,drop = F]),
                                                      scheme = scheme,
                                                      t_iter = paste0(q, "_", i), G = G, h = h_chunk[i],
                                                      windows = windows, center = center), silent = T)
                                   if(class(b) == "try-error"){
                                     is.err <- c(is.err, i)
                                     errs <- c(errs, b)
                                   }
                                   else{
                                     dist_chunk[i,] <- b
                                   }
                                 }
                                 out <- vector("list", 2)
                                 names(out) <- c("successes", "fails")
                                 if(length(is.err) > 0){
                                   out[[1]] <- cbind(parm_chunk[-is.err,], h = h_chunk[-is.err], dist_chunk[-is.err,])
                                   out[[2]] <- list(error_msg = errs,
                                                    error_parms = cbind(parm_chunk[is.err,,drop = F], h = h_chunk[is.err]))
                                 }
                                 else{
                                   out[[1]] <- cbind(parm_chunk, h = h_chunk, dist_chunk)
                                   out[[2]] <- list(error_msg = character(0),
                                                    error_parms = as.data.frame(matrix(NA, ncol = ncol(parm_chunk) + 1, nrow = 1)))
                                   colnames(out[[2]]$error_parms) <- c(colnames(parm_chunk), "h")
                                 }
                                 out
                               }


    parallel::stopCluster(cl)
    doSNOW::registerDoSNOW()
    gc();gc()
    dat <- dplyr::bind_rows(rvest::pluck(output, 1))
    errs <- rvest::pluck(output, 2)
    err_parms <- dplyr::bind_rows(rvest::pluck(errs, 2))
    err_msgs <- unlist(rvest::pluck(errs, 1))
    err_parms <- na.omit(err_parms)
    errs <- list(parms = err_parms, msgs = err_msgs)

    if(length(errs$msgs) > 0){
      warning("Errors occured on some simulations. See named element 'errors' in returned list for details.\n")
    }
    else{
      cat("Complete, no errors on any iterations.\n")
    }
    return(list(stats = dat, errors = errs))
  }
}
