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
ABC_on_hyperparameters <- function(x, phenotypes, iters,
                                   effect_distribution = rbayesB,
                                   parameter_distributions = list(pi = function(x) rbeta(x, 25, 1),
                                                                  d.f = function(x) runif(x, 1, 100),
                                                                  scale = function(x) rbeta(x, 1, 3)*100),
                                   h_dist = function(x) rep(.5, x),
                                   center = T,
                                   par = F,
                                   run_number = NULL){

  #============general subfunctions=========================
  euclid.dist <- function(o, p){
    dist <- sqrt(sum((o - p)^2))
    return(dist)
  }
  euclid.distribution.dist <- function(o, p){
    if(sum(p) == 0){
      return(rep(NA, length(names_diff_stats)))
    }
    dist <- compare_distributions(o, p)
    return(dist)
  }
  generate_pseudo_effects <- function(x, effect_distribution, parameters, h, center = T){
    pseudo_effects <- do.call(effect_distribution, c(list(n = nrow(x)), parameters))
    pseudo_phenos <- get.pheno.vals(x, pseudo_effects, h)$p
    if(center){
      pseudo_phenos <- pseudo_phenos - mean(pseudo_phenos)
    }
    return(list(e = pseudo_effects, p = pseudo_phenos))
  }

  #============ABC_scheme function for one rep=============
  scheme_D <- function(x, phenotypes, effect_distribution, parameters, h, center = center){
    pseudo <- generate_pseudo_effects(x, effect_distribution, parameters, h, center = center)
    dist <- euclid.distribution.dist(phenotypes, pseudo$p)
    return(dist)
  }

  #============initialization======================================
  # center the distribution if requested
  if(center){
    phenotypes <- phenotypes - mean(phenotypes)
  }


  # get the random values to run
  run_parameters <- vector("list", length = length(parameter_distributions))
  names(run_parameters) <- names(parameter_distributions)
  for(i in 1:length(run_parameters)){
    run_parameters[[i]] <- parameter_distributions[[i]](iters)
  }
  run_parameters <- as.data.frame(run_parameters)
  h <- h_dist(iters)

  # initialize storage
  dist_output <-  matrix(0, iters, ncol = length(names_diff_stats))
  colnames(dist_output) <- names_diff_stats

  #============ABC loop===================
  # run the ABC
  ## serial
  if(par == F){

    for(i in 1:iters){
      cat("Iter: ", i, ".\n")
      dist_output[i,] <- scheme_D(x = x,
                                  phenotypes = phenotypes,
                                  effect_distribution = effect_distribution,
                                  parameters = as.list(run_parameters[i,,drop = F]),
                                  h = h[i],
                                  center = center)
    }
    return(cbind(as.data.table(run_parameters), h = h, as.data.table(dist_output)))
  }


  # parallel
  else{
    # cl <- snow::makeSOCKcluster(par)
    # doSNOW::registerDoSNOW(cl)

    # divide up into ncore chunks
    chunks <- list(parms = split(run_parameters, (1:iters)%%par),
                   dist_output = split(as.data.frame(dist_output), (1:iters)%%par),
                   h = split(as.data.frame(h), (1:iters)%%par))

    # prepare reporting function
    progress <- function(n) cat(sprintf("Chunk %d out of", n), par, "is complete.\n")
    opts <- list(progress=progress)

    # output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
    #                            .options.snow = opts, .packages = c("data.table", "GeneArchEst")
    #                            ) %dopar% {
    #
    #                              parm_chunk <- chunks$parm[[q]]
    #                              h_chunk <- unlist(chunks$h[[q]])
    #                              dist_chunk <- chunks$dist_output[[q]]
    #
    #                              # run once per iter in this chunk
    #                              for(i in 1:nrow(parm_chunk)){
    #                                dist_chunk[i,] <- scheme_D(x = x,
    #                                                            phenotypes = phenotypes,
    #                                                            effect_distribution = effect_distribution,
    #                                                            parameters = as.list(parm_chunk[i,,drop = F]),
    #                                                            h = h_chunk[i],
    #                                                            center = center)
    #                              }
    #                              cbind(parm_chunk, h = h_chunk, dist_chunk)
    #                            }


    output <- vector("list", par)
    for(q in 1:par){

      parm_chunk <- chunks$parm[[q]]
      h_chunk <- unlist(chunks$h[[q]])
      dist_chunk <- chunks$dist_output[[q]]

      # run once per iter in this chunk
      for(i in 1:nrow(parm_chunk)){
        b <- try(scheme_D(x = x,
                          phenotypes = phenotypes,
                          effect_distribution = effect_distribution,
                          parameters = as.list(parm_chunk[i,,drop = F]),
                          h = h_chunk[i],
                          center = center), silent = T)
        if(class(b) == "try-error"){browser()}
        else{
          dist_chunk[i,] <- b
        }
      }
      output[[q]] <- cbind(parm_chunk, h = h_chunk, dist_chunk)
    }


    browser()
    # release cores and clean up
    parallel::stopCluster(cl)
    doSNOW::registerDoSNOW()
    gc();gc()
    output <- dplyr::bind_rows(output)
    return(output)
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
