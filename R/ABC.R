#' Conduct an ABC on a range of effect size distribution hyperparameters
#'
#' Runs Approximate Bayesian Computation across a range of marker effect size
#' distribution hyperparameters by comparing the distributions of phenotypes produced to those in the real data.
#'
#' Please note that no missing phenotypes or missing genotype data is permitted. Missing genotypes should be imputed
#' (such as with \code{\link{impute_and_phase_beagle}}) before this function is run.
#'
#' @param x matrix. Input genotypes, SNPs as rows, columns as individuals. Genotypes formatted as 0,1,2 for the major homozygote, heterozygote, and minor homozygote, respectively.
#' @param phenotypes numeric vector. Observed phenotypes, one per individual.
#' @param iters numeric. Number of ABC permutations to run.
#' @param effect_distribution function, default rbayesB. A function for a distribution of effect sizes. The first argument, n,
#'   should be the number of effects to draw. Other arguments must be named to match the names of the parameter distribution
#'   functions given in the parameter_distributions argument.
#' @param parameter_distributions List containing named functions, defualt
#'   list(pi = function(x) rbeta(x, 25, 1), d.f = function(x) runif(x, 1, 100), scale = function(x) rbeta(x, 1, 3)*100).
#'   A list containing functions for each parameter in the effect size distribution. The given functions should be named to
#'   match the parameter for which they produce distributions, and should take a single argument, x, which holds the number of
#'   parameter values to draw.
#' @param h_dist function, default function(x) rep(.5, x). A function for the distribution of heritability from which to draw
#'   h^2 values for each iteration. Must take a single argument, x, which holds the number of h values to draw. The default
#'   produces only heritabilities of .5.
#' @param center logical, default T. Determines if the phenotypes provided and generated during each iteration
#'   should be centered (have their means set to 0).
#' @param par numeric or FALSE, default FALSE. If numeric, the number of cores on which to run the ABC.
#'
#' @export
ABC_on_hyperparameters <- function(x, phenotypes, iters,
                                   effect_distribution = rbayesB,
                                   parameter_distributions = list(pi = function(x) rbeta(x, 25, 1),
                                                                  d.f = function(x) runif(x, 1, 100),
                                                                  scale = function(x) rbeta(x, 1, 3)*100),
                                   h_dist = function(x) rep(.5, x),
                                   center = TRUE,
                                   joint_res = NULL,
                                   joint_acceptance = NULL,
                                   joint_res_dist = NULL,
                                   par = FALSE, phased = FALSE,
                                   save_effects = FALSE, grid = 2000,
                                   save_effects_nkeep = ceiling(iters * .1),
                                   save_effects_keep_var = "ks"){


  if(!isFALSE(save_effects)){
    look_for_files <- paste0(save_effects, c(".bk", ".rds"))
    files_here <- file.exists(look_for_files)

    if(any(files_here)){
      stop(paste0("File(s) ", paste0(look_for_files[files_here], collapse = ", "), " already exist(s).\n"))
    }
  }

  #============ABC_scheme function for one rep=============
  scheme_D <- function(x, phenotypes, effect_distribution, parameters, h, center = center, phased = F, save_effects = FALSE){
    pseudo <- generate_pseudo_effects(x, effect_distribution, parameters, h, center = center, phased = phased)
    dist <- compare_distributions(phenotypes, pseudo$p)
    if(save_effects){
      return(list(dist = dist, effects = pseudo$e))
    }
    return(dist)
  }

  #============initialization======================================
  # center the distribution if requested
  if(center){
    phenotypes <- phenotypes - mean(phenotypes)
  }


  # get the random values to run
  run_parameters <- sample_parameters_from_distributions(parameter_distributions = parameter_distributions,
                                                         joint_res = joint_res, iters = iters,
                                                         joint_acceptance = joint_acceptance,
                                                         joint_res_dist = joint_res_dist,
                                                         reg_res = NULL, grid = grid)
  nparms <- ncol(run_parameters)

  # run_parameters <- vector("list", length = length(parameter_distributions))
  # names(run_parameters) <- names(parameter_distributions)
  # for(i in 1:length(run_parameters)){
  #   run_parameters[[i]] <- parameter_distributions[[i]](iters)
  # }
  # run_parameters <- as.data.frame(run_parameters)
  h <- h_dist(iters)

  # initialize storage
  dist_output <-  matrix(0, iters, ncol = length(names_diff_stats))
  colnames(dist_output) <- names_diff_stats

  #============ABC loop===================
  # run the ABC
  ## serial
  if(par == F | par == 1){

    for(i in 1:iters){
      cat("Iter:", i, "\n")

      d <- scheme_D(x = x,
                    phenotypes = phenotypes,
                    effect_distribution = effect_distribution,
                    parameters = as.list(run_parameters[i,,drop = F]),
                    h = h[i],
                    center = center,
                    phased = phased, save_effects = !isFALSE(save_effects))

      if(!isFALSE(save_effects)){
        data.table::fwrite(cbind(run_parameters[i,,drop = F], h = h[i], data.table::as.data.table(matrix(d$effects, nrow = 1))),
                           paste0(save_effects, ".tmp"), sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
        dist_output[i,] <- d$dist

      }
      else{
        dist_output[i,] <- d
      }
    }

    dist_output <- as.data.frame(dist_output)

    # filter saved effects
    # read effects into an FBM and save if requested (only the keepers!)
    if(!isFALSE(save_effects)){
      cat("Collating written effects into an FBM.\n")
      dist_output$keep <- FALSE
      dist_output$keep[order(dist_output[,save_effects_keep_var])[1:save_effects_nkeep]] <- TRUE # keep only the nkeep best runs


      effect_fbm <- bigstatsr::big_read(paste0(save_effects, ".tmp"),
                                        filter = which(dist_output$keep),
                                        select = 1:(nrow(x) + nparms + 1))

      file.remove(paste0(save_effects, ".tmp"))

      # save fbm
      effect_fbm <- effect_fbm$save()
    }

    return(cbind(as.data.table(run_parameters), h = h, as.data.table(dist_output)))
  }


  # parallel
  else{

    par <- min(par, iters)

    cl <- snow::makeSOCKcluster(par)
    doSNOW::registerDoSNOW(cl)

    # divide up into ncore chunks
    it_par <- (1:iters)%%par
    chunks <- list(parms = .smart_split(run_parameters, it_par),
                   dist_output = .smart_split(as.data.frame(dist_output), it_par),
                   h = .smart_split(as.data.frame(h), it_par))
    rm(run_parameters, dist_output)

    # prepare reporting function
    progress <- function(n) cat(sprintf("Chunk %d out of", n), par, "is complete.\n")
    opts <- list(progress=progress)

    need_removal <- paste0(save_effects, ".", 1:par, ".tmp.part")
    need_removal <- need_removal[which(file.exists(need_removal))]
    if(length(need_removal) > 0){
      file.remove(need_removal)

    }


    output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
                               .options.snow = opts, .packages = c("data.table", "GeneArchEst", "bigstatsr")
                               ) %dopar% {

                                 parm_chunk <- chunks$parm[[q]]
                                 h_chunk <- unlist(chunks$h[[q]])
                                 dist_chunk <- chunks$dist_output[[q]]
                                 is.err <- numeric(0)
                                 errs <- character(0)

                                 # run once per iter in this chunk
                                 for(i in 1:nrow(parm_chunk)){
                                   b <- try(scheme_D(x = x,
                                                     phenotypes = phenotypes,
                                                     effect_distribution = effect_distribution,
                                                     parameters = as.list(parm_chunk[i,,drop = F]),
                                                     h = h_chunk[i],
                                                     center = center,
                                                     phased = phased,
                                                     save_effects = !isFALSE(save_effects)), silent = T)


                                   if(class(b) == "try-error"){
                                     is.err <- c(is.err, i)
                                     errs <- c(errs, b)
                                   }
                                   else{
                                     if(!isFALSE(save_effects)){
                                       data.table::fwrite(cbind(parm_chunk[i,,drop = F], h = h_chunk[i], data.table::as.data.table(matrix(b$effects, nrow = 1))),
                                                          paste0(save_effects, ".", q, ".tmp.part"), sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
                                       dist_chunk[i,] <- b$dist
                                     }
                                     else{
                                       dist_chunk[i,] <- b
                                     }
                                   }
                                 }


                                 if(!isFALSE(save_effects)){
                                   dist_chunk$q_id <- q
                                   dist_chunk$iter_id <- 1:nrow(parm_chunk)
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

    # release cores and clean up
    parallel::stopCluster(cl)
    doSNOW::registerDoSNOW()
    gc();gc()

    cat("Run complete. Number of good runs per job:\n\t")
    dims <- unlist(lapply(purrr::map(output, 1), nrow))
    cat(paste0(1:length(dims), ": ", dims, "\n\t"))
    cat("Colnames:\n\t")
    cns <- lapply(purrr::map(output, 1), colnames)
    cns <- lapply(cns, paste0, collapse = "\t")
    cat(paste0(1:length(dims), ": ", cns, "\n\t"))



    dist_output <- dplyr::bind_rows(purrr::map(output, 1))
    errs <- purrr::map(output, 2)
    err_parms <- dplyr::bind_rows(purrr::map(errs, 2))
    err_msgs <- unlist(purrr::map(errs, 1))
    err_parms <- na.omit(err_parms)
    errs <- list(parms = err_parms, msgs = err_msgs)


    # read effects into an FBM and save if requested (only the keepers!)
    if(!isFALSE(save_effects)){
      cat("Collating written effects into an FBM.\n")
      dist_output$keep <- FALSE
      dist_output$keep[order(dist_output[,save_effects_keep_var])[1:save_effects_nkeep]] <- TRUE # keep only the nkeep best runs


      effect_fbm <- bigstatsr::FBM(save_effects_nkeep, nrow(x) + nparms + 1, backingfile = save_effects)


      progress <- 0
      for(i in 1:par){
        correct_iter <- which(dist_output$q_id == i)

        write_end <- sum(dist_output$keep[correct_iter])
        if(write_end == 0){
          next
        }
        tmpfile_bk <- stringi::stri_rand_strings(1, 20)
        while(file.exists(tmpfile_bk)){
          tmpfile_bk <- stringi::stri_rand_strings(1, 20)
        }
        tmp_fbm <- bigstatsr::big_read(paste0(save_effects, ".", i, ".tmp.part"),
                                       select = 1:(nrow(x) + nparms + 1),
                                       filter = which(dist_output$keep[correct_iter]),
                                       backingfile = tmpfile_bk)

        effect_fbm[(progress + 1):(progress + write_end),] <- tmp_fbm[]
        rm(tmp_fbm)
        gc();gc()



        file.remove(paste0(tmpfile_bk, ".bk"))
        file.remove(paste0(tmpfile_bk, ".rds"))

        file.remove(paste0(save_effects, ".", i, ".tmp.part"))

        progress <- progress + write_end
      }

      # clean dist_output
      dist_output <- dplyr::arrange(dist_output, q_id, iter_id) # sort by q and i to make sure our orders all match up.
      dist_output$q_id <- NULL
      dist_output$iter_id <- NULL

      # save fbm
      effect_fbm <- effect_fbm$save()
    }


    if(length(errs$msgs) > 0){
      warning("Errors occured on some ABC iterations. See named element 'errors' in returned list for details.\n")
      return(list(ABC_res = dist_output, errors = errs))
    }
    else{
      cat("Complete, no errors on any iterations.\n")
      return(dist_output)
    }
  }
}


#' Draw parameter values from the posteriour distribution of an ABC using a density kernal
#' @param x numeric, the number of values to draw.
#' @param res data.frame, the results of an ABC run, must include a "dist" column.
#' @param num_accepted numeric, either the proportion of runs to accept or the number of runs to accept.
#' @param parameters character, the parameter names for which to draw values.
#'
#' @export
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
