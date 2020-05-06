#' Conduct an ABC on a range of effect size distribution hyperparameters
#'
#' Runs Approximate Bayesian Computation across a range of marker effect size
#' distribution hyperparameters using one of three different schemes in order to determine
#' the hyperparamters that generate a distribution most like the real genomic architecture of the trait.
#'
#' ABC schemes: \itemize{
#'     \item{"A": }{Genomic Data -> prediction with prior hyperparameters -> compare predicted phenotypes to real phenotypes.}
#'     \item{"B": }{Genomic Data -> generate psuedo marker effects using distribution with prior hyperparameters -> prediction with defaults -> compare predicted phenotypes to real phenotypes.}
#'     \item{"C": }{Part 1: Genomic Data -> generate psuedo marker effects using distribution with prior hyperparameters -> "pseudo" predicted marker effects.
#'                  Part 2: Genomic Data -> prediction with defaults -> direct estimated marker effects.
#'                  Part 3: Compare direct to "pseudo" estimated marker effects.}
#' }
#'
#' @param x matrix. Input genotypes, SNPs as rows, columns as individuals. Genotypes formatted as 0,1,2 for the major homozygote, heterozygote, and minor homozygote, respectively.
#' @param phenotypes numeric vector. Observed phenotypes, one per individual.
#' @param iters numeric. Number of ABC permutations to run.
#' @param pi_func function, default function(x) rbeta(x, 25, 1). A distribution function for generating pi prior. Should take only one argument (n, the number of iters).
#' @param df_func function, default NULL. A distribution function for generating df prior. Should take only one argument (n, the number of iters).numeric vector of length 2, default NULL. Range (min, max) of degrees of freedom values to run.
#' @param scale_func function, default NULL. A distribution function for generating scale prior. Should take only one argument (n, the number of iters).numeric vector of length 2, default NULL. Range (min, max) of scale values to run.
#' @param h numeric, default NULL. Heritability to use. Will take a range in the future.
#' @param julia.path character, defualt "julia". File path to the julia executable, required for JWAS.
#' @param chain_length numeric, default 100000. Length of the MCMC chains used in each step of the ABC.
#' @param burnin numeric, default 5000. Number of MCMC chain steps discarded at the start of each MCMC.
#' @param thin numeric, default 100. Number of MCMC chain steps discarded between each sample used to form the posterior.
#' @param method character, default "bayesB". The marker effect size distribution/prediction method to use.
#' @param ABC_scheme character, default "A". The ABC_scheme to use. See details.
#' @param par numeric or FALSE, default FALSE. If numeric, the number of cores on which to run the ABC.
#' @param run_number numeric, default NULL. Controls how the itermediate output directories are named. If numeric, will be named for the number, otherwise, will be named for the iteration.
#'
#' @export
ABC_on_hyperparameters <- function(x, phenotypes, iters, center = T, pi_func = function(x) rbeta(x, 25, 1),
                                   df_func = NULL,  scale_func = NULL, h = NULL,
                                   julia.path = "julia", chain_length = 100000,
                                   burnin = 5000,
                                   thin = 100, method = "BayesB", ABC_scheme = "A",
                                   par = F, run_number = NULL, est_h = F, save_effects = T,
                                   joint_res = NULL, joint_acceptance = NULL, joint_res_dist = "ks.D",
                                   delta = .5, pcut = 0.0005){

  # ks <- which(matrixStats::rowSums2(x)/ncol(x) >= 0.05)
  #============general subfunctions=========================
  euclid.dist <- function(o, p){
    dist <- sqrt(sum((o - p)^2))
    return(dist)
  }
  euclid.distribution.dist <- function(o, p){
    if(sum(p) == 0){
      return(rep(NA, number_diff_stats))
    }
    dist <- compare_distributions(o, p)
    return(dist)
  }
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

  #============ABC_scheme functions for one rep=============
  scheme_A <- function(x, phenotypes, pi, method, t_iter){
    p <- pred(x, pi = pi, phenotypes = phenotypes, julia.path = julia.path,
              burnin = burnin, thin = thin, chain_length = chain_length,
              prediction.program = "JWAS", prediction.model = method, runID = t_iter, verbose = F)

    dist <- euclid.dist(phenotypes, p$est.phenos)
    return(dist)
  }
  scheme_B <- function(x, phenotypes, pi, df, scale, method, t_iter, center = center){
    pseudo <- generate_pseudo_effects(x, pi, df, scale, method, h, center = center)

    p <- pred(x, phenotypes = pseudo$p,
              burnin = burnin, thin = thin, chain_length = chain_length,
              prediction.program = "BGLR", prediction.model = method, runID = t_iter, verbose = F)

    p.phenos <- as.vector(convert_2_to_1_column(p$x)%*%p$output.model$mod$ETA[[1]]$b)

    dist <- euclid.distribution.dist(phenotypes, p.phenos)
    return(return(list(dist = dist, e = pseudo$e)))
  }
  scheme_C <- function(x, phenotypes, r.p.eff, pi, df, scale, method, t_iter, center = center){
    pseudo <- generate_pseudo_effects(x, pi, df, scale, method, h, center = center)

    cat("Beginning pseudo data", method, "run.\n")
    pseudo.pred <-pred(x, phenotypes = pseudo$p,
                       burnin = burnin, thin = thin, chain_length = chain_length,
                       prediction.program = "BGLR", prediction.model = method,
                       runID = paste0(t_iter, "_pseudo"), verbose = F)

    upper <- mean(pseudo.pred$output.model$mod$ETA[[1]]$b) + sd(pseudo.pred$output.model$mod$ETA[[1]]$b) * pcut
    lower <- mean(pseudo.pred$output.model$mod$ETA[[1]]$b) - sd(pseudo.pred$output.model$mod$ETA[[1]]$b) * pcut
    peaks.p <- findpeaks_multi(cbind(meta[pseudo.pred$kept.snps,], effect = pseudo.pred$output.model$mod$ETA[[1]]$b), delta,
                               pcut = c(lower, upper), "group", pvals = F)

    upper <- mean(r.p.eff) + sd(r.p.eff) * pcut
    lower <- mean(r.p.eff) - sd(r.p.eff) * pcut
    peaks.o <- findpeaks_multi(cbind(meta[pseudo.pred$kept.snps,], effect = r.p.eff), delta, pcut = c(upper, lower), "group", pvals = F)

    dist <- euclid.distribution.dist(r.p.eff, pseudo.pred$output.model$mod$ETA[[1]]$b)
    dist.peaks <- compare_peaks(peaks.o, peaks.p)
    names(dist.peaks)[-1] <- paste0("peak_", names(dist.peaks)[-1])


    return(list(dist = c(dist, dist.peaks), e = pseudo$e))

  }
  scheme_D <- function(x, phenotypes, pi, df, scale, method, center = center){
    pseudo <- generate_pseudo_effects(x, pi, df, scale, method, h, center = center)
    dist <- euclid.distribution.dist(phenotypes, pseudo$p)
    return(list(dist = dist, e = pseudo$e))
  }
  scheme_E <- function(x, phenotypes, real_pi_dist, pi, df, scale, method, t_iter, G, center = center){
    pseudo <- generate_pseudo_effects(x, pi, df, scale, method, h, center = center)

    pseudo_pi <- pred(x, phenotypes = pseudo$p,
                      prediction.program = "GMMAT",
                      maf.filt = F, runID = paste0(t_iter, "gmmat_pseudo"),
                      pass_G = G)$e.eff$PVAL

    peaks.p <- findpeaks_multi(cbind(meta, PVAL = pseudo_pi, logp = -log10(pseudo_pi)), delta, pcut = pcut, "group")
    peaks.o <- findpeaks_multi(cbind(meta, PVAL = real_pi_dist, logp = -log10(real_pi_dist)), delta, pcut = pcut, "group")

    dist <- euclid.distribution.dist(real_pi_dist, pseudo_pi)
    dist.peaks <- compare_peaks(peaks.o, peaks.p)
    names(dist.peaks)[-1] <- paste0("peak_", names(dist.peaks)[-1])

    return(list(dist = c(dist, dist.peaks), e = pseudo$e))
  }

  loop_func <- function(x, phenotypes, pi, df, scale, method, scheme, t_iter, r.p.phenos = NULL, r.p.eff = NULL, real_pi_dist = NULL, G = NULL, center = center){
    if(scheme == "A"){
      dist <- scheme_A(x, phenotypes, pi, method, t_iter)
    }
    else if(scheme == "B"){
      dist <- scheme_B(x, phenotypes, pi, df, scale, method, t_iter, center = center)
    }
    else if(scheme == "C"){
      dist <- scheme_C(x, phenotypes, r.p.eff, pi, df, scale, method, t_iter, center = center)
    }
    else if(scheme == "D"){
      dist <- scheme_D(x, phenotypes, pi, df, scale, method, center = center)
    }
    else if(scheme == "E"){
      dist <- scheme_E(x, phenotypes, real_pi_dist, pi, df, scale, method, t_iter, G = G, center = center)
    }
    return(dist)
  }

  #============initialization======================================
  # center the distribution if requested
  if(center){
    phenotypes <- phenotypes - mean(phenotypes)
  }


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
  # if any joint parameter priors, calculate and disambiguate
  if(length(joint_params) > 0){
    joint_params <- gen_parms(iters, joint_res, joint_acceptance, joint_params, dist.var = joint_res_dist)
    colnames(joint_params) <- paste0("run_", colnames(joint_params), "s")
    for(i in 1:ncol(joint_params)){
      assign(colnames(joint_params)[i], joint_params[,i])
    }
  }


  # initialize storage
  out <- cbind(pi = run_pis, df = run_dfs, scale = run_scales, matrix(0, length(run_pis), ncol = number_diff_stats))
  colnames(out)[-c(1:3)] <- names_diff_stats

  # if doing a method where prediction needs to be run on the real data ONCE, or if h should be estimated, do that now:
  if(ABC_scheme == "C" | est_h == T){
    cat("Beginning real data", method, "run.\n")
    real.pred <- pred(x, phenotypes = phenotypes,
                      burnin = burnin, thin = thin, chain_length = chain_length,
                      prediction.program = "BGLR", prediction.model = method, runID = "real_pred", verbose = F)
    if(ABC_scheme == "C"){
      #r.p.phenos <- as.vector(convert_2_to_1_column(real.pred$x)%*%real.pred$output.model$mod$ETA[[1]]$b)
      r.p.phenos <- NULL
      r.p.eff <- real.pred$output.model$mod$ETA[[1]]$b
    }
    if(est_h == T){
      h <- real.pred$h
    }
    if(ABC_scheme != "C"){
      rm(real.pred)
      r.p.phenos <- NULL
      r.p.eff <- NULL
    }
  }
  else{
    r.p.phenos <- NULL
  }

  # can pass a g matrix forward and get comparison p-values once if doing scheme E.
  if(ABC_scheme == "E"){
    ind.genos <- convert_2_to_1_column(x)
    colnames(ind.genos) <- paste0("m", 1:ncol(ind.genos)) # marker names
    rownames(ind.genos) <- paste0("s", 1:nrow(ind.genos)) # ind IDS

    G <- AGHmatrix::Gmatrix(ind.genos, missingValue = NA, method = "Yang", maf = 0.05)
    colnames(G) <- rownames(ind.genos)
    rownames(G) <- rownames(ind.genos)

    real_pi_dist <- pred(x, phenotypes = phenotypes, prediction.program = "GMMAT", maf.filt = F, runID = "gmmat_real",
                         pass_G = G)$e.eff$PVAL
  }
  else{
    real_pi_dist <- NULL
    G <- NULL
  }

  #============ABC loop===================
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
      tout <- loop_func(x, phenotypes, out[i,"pi"], out[i,"df"], out[i,"scale"], method,
                        ABC_scheme, t_iter = rn, r.p.phenos = r.p.phenos, r.p.eff = r.p.eff,
                        real_pi_dist = real_pi_dist, G = G, center = center)
      out[i, 4:ncol(out)] <- tout$dist
      if(save_effects){
        data.table::set(out.effects, j = i,  value = tout$e)
      }
    }
    colnames(out)[4:ncol(out)] <- names(tout$dist)
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



    output <- foreach::foreach(i = 1:par, .inorder = FALSE, .errorhandling = "pass",
                               .options.snow = opts, .packages = c("data.table", "GeneArchEst")
                               ) %dopar% {

                                 out <- chunks[[i]]
                                 if(save_effects){
                                   out.effects <- data.table::as.data.table(matrix(NA, nrow = nrow(x), ncol = nrow(out)))
                                 }

                                 # run once per iter in this chunk
                                 for(j in 1:nrow(out)){
                                   if(is.numeric(run_number)){rn <- run_number}
                                   else{rn <- j}
                                   tout <- loop_func(x, phenotypes, out[j,"pi"], out[j,"df"], out[j,"scale"], method, ABC_scheme, t_iter = rn,
                                                     r.p.phenos = r.p.phenos, r.p.eff = r.p.eff, real_pi_dist = real_pi_dist, G = G, center = center)
                                   out[j, 4:ncol(out)] <- tout$dist
                                   if(save_effects){
                                     data.table::set(out.effects, j = j,  value = tout$e)
                                   }
                                 }
                                 colnames(out)[4:ncol(out)] <- names(tout$dist)
                                 if(save_effects){
                                   out <- list(dists = out, effects = out.effects)
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
    return(list(dists = out, effects = out.effects))
  }
  else{
    return(list(dists = out))
  }
}


#' Draw parameter values from the posteriour distribution of an ABC using a density kernal
#' @param x numeric, the number of values to draw.
#' @param res data.frame, the results of an ABC run, must include a "dist" column.
#' @param num_accepted numeric, either the proportion of runs to accept or the number of runs to accept.
#' @param parameters character, the parameter names for which to draw values.
gen_parms <- function(x, res, num_accepted, parameters, grid = 2000, dist.var = "ks"){
  if(length(parameters) != 2){
    stop("Only two parameters accepted at the moment.\n")
  }

  # assign accepted runs
  if(num_accepted < 1){
    qcut <- num_accepted
  }
  else{
    qcut <- num_accepted/nrow(res)
  }

  # remove results with NA values for the dist var
  bads <- which(is.na(res[,dist.var]))
  if(length(bads) > 0){
    res <- res[-bads,]
  }


  res$hits <- ifelse(res[,dist.var] <= quantile(res[,dist.var], qcut), 1, 0)

  # generate kernal
  op <- GenKern::KernSur(res[res$hits == 1,which(colnames(res) == parameters[1])],
                         res[res$hits == 1,which(colnames(res) == parameters[2])],
                         range.x = c(min(res[,which(colnames(res) == parameters[1])]),
                                     max(res[,which(colnames(res) == parameters[1])])),
                         range.y = c(min(res[,which(colnames(res) == parameters[2])]),
                                     max(res[,which(colnames(res) == parameters[2])])),
                         ygridsize = grid, xgridsize = grid)

  # sample from kernal and readjust to rows/columns
  cells <- sample(1:length(op$zden), x, replace = T, prob = op$zden)
  rows <- cells %% ncol(op$zden)
  cols <- floor(cells/nrow(op$zden)) + 1
  if(any(rows == 0)){
    cols[rows == 0] <- cols[rows == 0] - 1
    rows[rows == 0] <- nrow(op$zden)
  }

  # grab parameter values and return
  vals <- data.frame(p1 = op$xords[rows], p2 = op$yords[cols])
  colnames(vals) <- parameters
  return(vals)
}
