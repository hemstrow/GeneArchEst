#' @export
sim_gen <- function(x, meta, iters, chr = "chr", pi_func = function(x) rbeta(x, 25, 1),
                    df_func = NULL,  scale_func = NULL, h_func = NULL, center = T, scheme = "gwas",
                    method = "BayesB",
                    burnin = burnin, thin = thin, chain_length = chain_length,
                    par = F, save_effects = T,
                    joint_res = NULL, joint_acceptance = NULL, joint_res_dist = "ks.D",
                    peak_delta = .5, peak_pcut = 0.0005, window_sigma = 50, run_number = NULL){

  #============general subfunctions=========================
  generate_pseudo_effects <- function(x, pi, df, scale, method, h = NULL, center = T){
    if(method == "BayesB"){

      pseudo_effects <- rbayesB(nrow(x), pi, df, scale)
      pseudo_phenos <- get.pheno.vals(x, pseudo_effects, h)$p
      if(center){
        pseudo_phenos <- pseudo_phenos - mean(pseudo_phenos)
      }
    }
    return(list(e = pseudo_effects, p = pseudo_phenos))
  }

  #============schem functions for one simulation=============
  gp <- function(x, pi, df, scale, method, t_iter, h, windows, center = center){
    pseudo <- generate_pseudo_effects(x, pi, df, scale, method, h, center = center)

    cat("Beginning pseudo data", method, "run.\n")
    pseudo.pred <-pred(x, phenotypes = pseudo$p,
                       burnin = burnin, thin = thin, chain_length = chain_length,
                       prediction.program = "BGLR", prediction.model = method,
                       runID = paste0(t_iter, "_pseudo"), verbose = F)

    stats <- dist_desc(pseudo.pred, meta, windows, peak_delta, peak_pcut, chr, pvals = F)

    return(list(stats = stats, e = pseudo$e))
  }
  gwas <- function(x, pi, df, scale, method, t_iter, G, h, windows, center = center){
    pseudo <- generate_pseudo_effects(x, pi, df, scale, method, h, center = center)
    pseudo_pi <- pred(x, phenotypes = pseudo$p,
                      prediction.program = "GMMAT",
                      maf.filt = F, runID = paste0(t_iter, "gmmat_pseudo"),
                      pass_G = G)$e.eff$PVAL

    stats <- dist_desc(pseudo_pi, meta, windows, peak_delta, peak_pcut, chr, pvals = T)
    return(list(stats = stats, e = pseudo$e))
  }

  loop_func <- function(x, pi, df, scale, method, scheme, t_iter, G = NULL, h, windows, center = center){
    if(scheme == "gp"){
      dist <- gp(x, pi, df, scale, method, t_iter, h, windows, center = center)
    }
    else if(scheme == "gwas"){
      dist <- gwas(x, pi, df, scale, method, t_iter, G = G, h = h, windows, center = center)
    }
    return(dist)
  }

  #============prep for simulations======================================
  # get the random values to run
  joint_params <- character()
  if(!is.null(df_func)){
    if(!is.function(df_func)){
      if(df_func == "joint"){
        joint_params <- "df"
      }
      else{
        stop("df func is neither a function or 'joint'.\n")
      }
    }
    else{
      run_dfs <- df_func(iters)
    }
  }
  else{
    run_dfs <- rep(NA, iters)
  }
  if(!is.null(scale_func)){
    if(!is.function(scale_func)){
      if(scale_func == "joint"){
        joint_params <- c(joint_params, "scale")
      }
      else{
        stop("scale func is neither a function or 'joint'.\n")
      }
    }
    else{
      run_scales <- scale_func(iters)
    }
  }
  else{
    run_scales <- rep(NA, iters)
  }
  if(!is.function(pi_func)){
    if(pi_func == "joint"){
      joint_params <- c(joint_params, "pi")
    }
    else{
      stop("pi func is neither a function or 'joint'.\n")
    }
  }
  else{
    run_pis <- pi_func(iters)
  }
  run_hs <- h_func(iters)

  # if any joint parameter priors, calculate and disambiguate
  if(length(joint_params) > 0){
    joint_params <- gen_parms(iters, joint_res, joint_acceptance, joint_params, dist.var = joint_res_dist)
    colnames(joint_params) <- paste0("run_", colnames(joint_params), "s")
    for(i in 1:ncol(joint_params)){
      assign(colnames(joint_params)[i], joint_params[,i])
    }
  }

  # can pass a g matrix forward once if doing gwas
  if(scheme == "gwas"){
    ind.genos <- convert_2_to_1_column(x)
    colnames(ind.genos) <- paste0("m", 1:ncol(ind.genos)) # marker names
    rownames(ind.genos) <- paste0("s", 1:nrow(ind.genos)) # ind IDS

    G <- AGHmatrix::Gmatrix(ind.genos, missingValue = NA, method = "Yang", maf = 0.05)
    colnames(G) <- rownames(ind.genos)
    rownames(G) <- rownames(ind.genos)
  }
  else{
    G <- NULL
  }

  # pre-run the window function
  windows <- mark_windows(meta, window_sigma, chr)


  #============run the simulations===========================
  # initialize storage
  out <- cbind(pi = run_pis, df = run_dfs, scale = run_scales, h = run_hs, matrix(NA, length(run_dfs), number_descriptive_stats))


  # run the ABC
  ## serial
  if(par == F){
    # initialize effects storage
    if(save_effects){
      out.effects <- data.table::as.data.table(matrix(NA, nrow = nrow(x), ncol = iters))
    }

    for(i in 1:iters){
      cat("Iter: ", i, ".\n")
      if(is.numeric(run_number)){rn <- run_number}
      else{rn <- i}

      tout <- loop_func(x = x, pi = out[i,"pi"], df = out[i,"df"], scale = out[i,"scale"],
                        method = method, scheme = scheme, t_iter = rn, G = G, h = out[i,"h"],
                        windows = windows, center = center)

      out[i, 5:ncol(out)] <- tout$stats

      if(save_effects){
        data.table::set(out.effects, j = i,  value = tout$e)
      }
    }
    colnames(out)[5:ncol(out)] <- names(tout$stats)
  }


  # parallel
  else{
    parms <- out[,-ncol(out)]
    cl <- snow::makeSOCKcluster(par)
    doSNOW::registerDoSNOW(cl)

    # divide up into ncore chunks
    chunks <- split(as.data.frame(out), (1:nrow(out))%%par)

    # prepare reporting function
    progress <- function(n) cat(sprintf("Chunk %d out of", n), par, "is complete.\n")
    opts <- list(progress=progress)

    output <- foreach::foreach(i = 1:par, .inorder = FALSE,
                               .options.snow = opts,
                               .export = c("rbayesB", "get.pheno.vals", "e.dist.func", "pred",
                                           "convert_2_to_1_column", "calc_dist_stats",
                                           "cucconi.stat", "lepage.stat"), .packages = c("data.table", "inline"),
                               .noexport = "weighted.colSums") %dopar% {

                                 # remake the weighted.colSums function, since the inline part doesn't work in packages
                                 # and writing the same code as a .cpp is actually much slower
                                 weighted.colSums <- function(data, weights){
                                   return(crossprod(t(data), weights))
                                 }


                                 out <- chunks[[i]]
                                 if(save_effects){
                                   out.effects <- data.table::as.data.table(matrix(NA, nrow = nrow(x), ncol = nrow(out)))
                                 }

                                 # run once per iter in this chunk
                                 for(j in 1:nrow(out)){
                                   if(is.numeric(run_number)){rn <- run_number}
                                   else{rn <- j}
                                   tout <- loop_func(x = x, pi = out[j,"pi"], df = out[j,"df"], scale = out[j,"scale"],
                                                     method = method, scheme = scheme, t_iter = rn, G = G, h = out[j,"h"],
                                                     windows = windows, center = center)
                                   out[j, 5:ncol(out)] <- tout$stats
                                   if(save_effects){
                                     data.table::set(out.effects, j = j,  value = tout$e)
                                   }
                                 }
                                 colnames(out)[5:ncol(out)] <- names(tout$stats)
                                 if(save_effects){
                                   out <- list(stats = out, effects = out.effects)
                                 }
                                 out
                               }


    # release cores and clean up
    parallel::stopCluster(cl)
    doSNOW::registerDoSNOW()
    gc();gc()

    # bind, relies on rvest::pluck to grab only the first or only the second part
    if(save_effects){
      out <- dplyr::bind_rows(rvest::pluck(output, 1))
      out.effects <- dplyr::bind_cols(rvest::pluck(output, 2))
    }
    else{
      out <- dplyr::bind_rows(output)
    }
  }

  if(save_effects){
    return(list(stats = out, effects = out.effects))
  }
  else{
    return(list(stats = out))
  }
}
