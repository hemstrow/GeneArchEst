#=======function to do growth and selection=======
# note: init means are we simply intiating a population under selection. We'll need to keep all markers if so.
#' Simulate selection and population growth
#'
#' @param genotypes unphased matrix of genotypes, one row per snp and one column
#'   per chromosome/gene copy. Copies from individuals must be sequential.
#' @param meta data.frame of snp metadata, with the first column containing
#'   chromosome info, the second containing position in base pairs. Effects can
#'   either be in one column, named "effects", if all effects are additive or in
#'   three columns named "effect_0", "effect_1", and "effect_2" for the effects
#'   for homozygote "0" individuals, heterozygotes, and homozygote "2"
#'   individuals, respectively.
#' @param chr.length numeric vector of chromosome lengths \emph{in the same
#'   order as would be given by unique(meta[,1])}.
#' @param h numeric. Narrow-sense heritability, between 0 and 1.
#' @param BVs numeric vector, breeding values for starting individuals. Either
#'   this or \code{phenotypes} is required.
#' @param phenotypes numeric vector, phenotypes for starting individuals. Either
#'   this or \code{BVs} is required.
#' @param gens numeric, maximum number of generations to simulate.
#' @param growth.function function, default \code{function(n) logistic_growth(n,
#'   500, 2)} the growth function to use during simulation, needs to take an
#'   argument \code{n}, the number of current individuals.
#' @param survival.function function, default \code{function(phenotypes,
#'   opt_pheno, ...) BL_survival(phenotypes, opt_pheno, omega = 1)}. Function
#'   for calculating survival probabilities of individuals. Needs two arguments:
#'   \code{phenotypes}: individual phenotypes and \code{opt_pheno}: optimum
#'   phenotype in a given generation.
#' @param selection.shift.function function, default \code{function(opt, iv)
#'   optimum_constant_slide(opt, iv, 0.3)}. Function for determining how much
#'   the fitness optimum changes each generation. Expects two arguments:
#'   \code{opt}: the current fitness optimum. \code{iv}: the initial genetic
#'   variance (the fitness optimum can slide as a function of the amount of
#'   initial genetic variance in the population).
#' @param thin logical, default TRUE. If TRUE, sites with no effect will be
#'   trimmed from analysis and reports unless all loci have no effects. This can
#'   massively reduce run time if few sites have effects. Set to FALSE if, for
#'   example, you wish to examine the effect of selection on linked, neutral
#'   loci.
#' @param fitnesses logical, default FALSE. If TRUE, effects are presumed to be
#'   fitness coefficients  rather than sheer effects. Net fitness is therefore
#'   the product off all effects. Fitnesses greater than 1 will be set to one,
#'   and less than zero will be set to zero. Must provide "effect_0",
#'   "effect_1", and "effect_2" columns to meta noting the fitness of each
#'   genotype.
#' @param stop_if_no_variance logical, default FALSE. If TRUE, will stop and
#'   return an error if no genetic variance remains in the population.
#' @export
gs <- function(genotypes,
               meta,
               phenotypes = NULL,
               BVs = NULL,
               h,
               gens,
               chr.length,
               migration = FALSE,
               starting.surv.opt = NULL,
               mutation = 0,
               mutation.effect.function = NULL,
               growth.function = function(n) logistic_growth(n, 500, 2),
               survival.function = function(phenotypes, opt_pheno, ...) BL_survival(phenotypes, opt_pheno, omega = 1),
               selection.shift.function = function(opt, iv) optimum_constant_slide(opt, iv, 0.3),
               rec.dist = function(n) rpois(n, lambda = 1),
               var.theta = 0,
               plot_during_progress = FALSE,
               facet = "group",
               do.sexes = TRUE,
               fitnesses = FALSE,
               effects = NULL,
               init = FALSE,
               thin = TRUE,
               thin_fixed = TRUE,
               verbose = FALSE,
               print.all.freqs = F,
               print.all.thinned.freqs = FALSE,
               sampling_point = "migrants",
               model = NULL,
               stop_if_no_variance = FALSE,
               K_thin_post_surv = NULL){

  #============subfuctions============
  # function completed each loop (done this way to allow for ease of multiple populations)
  one_gen <- function(genotypes, phenotypes,
                      BVs,
                      effects,
                      opt, fitnesses,
                      survival.function,
                      K_thin_post_surv,
                      meta,
                      rec.dist,
                      chr.length,
                      do.sexes,
                      h.av,
                      model,
                      selection.shift.function,
                      mutation,
                      pass_surv_genos = FALSE){

    #=========survival====
    if(length(unique(phenotypes)) == 1 & phenotypes[1] != 1 & stop_if_no_variance){stop("No genetic variance left.\n")}

    #survival:
    s <- rbinom(ncol(genotypes)/2, 1, #survive or not? Number of draws is the pop size in prev gen, survival probabilities are determined by the phenotypic variance and optimal phenotype in this gen.
                survival.function(phenotypes, opt))

    #if the population has died out, stop.
    if(sum(s) <= 1){
      #adjust selection optima
      return(NA)
    }

    # if doing carrying capacity on BREEDERS, not offspring, thin to K here.
    if(!is.null(K_thin_post_surv)){
      if(sum(s) > K_thin_post_surv){
        s[which(s == 1)][sample(sum(s), sum(s) - K_thin_post_surv, F)] <- 0
      }
    }

    #===============figure out next generation===============
    #what is the pop size after growth?
    offn <- round(growth.function(sum(s)))


    #make a new x with the survivors
    genotypes <- genotypes[, .SD, .SDcols = which(rep(s, each = 2) == 1)] #get the gene copies of survivors

    # save parents if needed
    if(pass_surv_genos){
      final_genotypes <- genotypes
      final_meta <- meta
      final_effects <- effects
    }

    #=============do random mating, adjust selection, get new phenotype scores, get ready for next gen====
    genotypes <- rand.mating(x = genotypes, N.next = offn, meta = meta, rec.dist = rec.dist, chr.length, do.sexes,
                             mutation = mutation)
    if(is.null(genotypes)){return(NA)} # return NA if empty (no individuals because everything was one sex)

    if(pass_surv_genos){
      return(list(genotypes = genotypes, final_genotypes = final_genotypes, final_meta = final_meta, final_effects = final_effects))
    }
    else{
      return(genotypes)
    }
  }

  do_mutation <- function(x.next, mutation, chr.length, meta, uf){
    #================mutation===================================
    # figure out number of mutations in each individual
    nmut <- lapply(x.next, function(z){
      if(!is.null(z)){return(rpois(ncol(z), mutation*sum(chr.length)))}
      else{return(NULL)}
    })
    pops <- unlist(lapply(nmut, length))
    nmut <- unlist(nmut)

    # exit if no mutations
    if(sum(nmut) == 0){
      meta$new <- FALSE
      return(list(genotypes = x.next, meta.next = meta, effects = effects))
    }

    # figure out positions for the mutations
    positions <- sample(sum(chr.length), sum(nmut), TRUE)
    chrs <- as.numeric(cut(positions, breaks = c(0, cumsum(chr.length))))
    positions <- positions - c(0, cumsum(chr.length)[-length(chr.length)])[chrs]
    inds <- rep(1:sum(pops), nmut) # which inds each mutation is in
    mut_info <- data.table::data.table(chrs = chrs, position = positions, ind = inds)
    # search for duplicates and re-sample them if they exist (unlikely)
    dups <- duplicated(mut_info[,1:3])
    loop_count <- 1
    while(any(dups)){
      if(loop_count > 20){
        stop("Tried to fix a double mutation at one locus in one individual 20 times without success--are your chromosomes too small for your mutation rate?\n")
      }
      npos <- sample(sum(chr.length), sum(dups), TRUE)
      nchrs <- as.numeric(cut(npos, breaks = c(0, cumsum(chr.length))))
      npos <- npos - c(0, cumsum(chr.length)[-length(chr.length)])[nchrs]
      chrs[dups] <- nchrs
      positions[dups] <- npos

      mut_info$chrs[dups] <- nchrs
      mut_info$position[dups] <- npos

      loop_count <- loop_count + 1
      dups <- duplicated(mut_info[,1:3])
    }

    # finish making the info df
    mut_info <- dplyr::arrange(mut_info, inds, chrs, positions)
    mut_info_index <- unique(mut_info[,1:2])
    mut_info_index$row <- 1:nrow(mut_info_index)
    mut_info <- merge(mut_info, mut_info_index)

    # get effects and thin if doing so
    mut.eff <- matrix(NA, nrow(mut_info_index), npops)
    if(length(mutation.effect.function) > 1){
      for(i in 1:npops){
        mut.eff[,i] <- mutation.effect.function[[i]](nrow(mut_info_index))
      }
    }
    else{
      mut.eff[,1:npops] <- mutation.effect.function(nrow(mut_info_index))
    }


    # thin if requested
    if(thin){
      zeros <- which(rowSums(mut.eff) == 0)

      # cut and return if empty
      if(length(zeros) == nrow(mut_info_index)){
        meta$new <- FALSE
        return(list(x.next = x.next, meta.next = meta, effects = effects))
      }

      # otherwise adjust row info
      mut_info_index <- mut_info_index[-zeros,]
      mut_info <- mut_info[-which(mut_info$row %in% zeros),]
      mut_info_index$row <- 1:nrow(mut_info_index)
      mut_info$row <- NULL
      mut_info <- merge(mut_info, mut_info_index)
      mut.eff <- mut.eff[-zeros,,drop=FALSE]
    }

    # make the data with mutations
    mut.x <- matrix(0, nrow(mut_info_index), sum(pops))
    mut.x[mut_info$row + ((mut_info$ind - 1) * nrow(mut.x))] <- 1
    mut.x <- as.data.table(mut.x)
    mut_info$chrs <- uf[mut_info$chrs]

    # determine if any new sites overlap the existing ones
    overlap_i_in_ref <- which(paste0(meta[,1], "_", meta[,2]) %in% paste0(mut_info$chrs, "_", mut_info$position))
    if(length(overlap_i_in_ref) > 0){
      # overlap_i_in_mut <- match(paste0(meta[overlap_i_in_ref,1], "_", meta[overlap_i_in_ref,2]),
      #                           paste0(mut_info$chrs, "_", mut_info$position))

      overlap_i_in_mut <- which(paste0(mut_info$chrs, "_", mut_info$position) %in%
                                  paste0(meta[overlap_i_in_ref,1], "_", meta[overlap_i_in_ref,2]))

      # flip any overlaps, then remove these from the data
      for(i in 1:npops){
        this_pops <- which(mut_info[overlap_i_in_mut,]$ind > c(0, cumsum(pops))[i] &
                             mut_info[overlap_i_in_mut,]$ind <= c(0, cumsum(pops))[i + 1])
        # move to next pop if none in this one
        if(length(this_pops) == 0){next}

        # find the overlaps in this pop
        t.overlaps <- mut_info[overlap_i_in_mut,][this_pops,]
        t.overlaps$ind <- t.overlaps$ind - c(0, cumsum(pops))[i]
        t.overlap.cols <- t.overlaps$ind

        # grab those rows, locate mutations, flip, and then update
        overlap_fill <- as.matrix(x.next[[i]][,..t.overlap.cols])
        overlap_indices <- t.overlaps$row + nrow(overlap_fill)*((1:nrow(t.overlaps)) - 1)
        overlap_fill[overlap_indices] <- ifelse(overlap_fill[overlap_indices] == 0, 1, 0)
        overlap_fill <- data.table::as.data.table(matrix(overlap_fill, nrow = nrow(x.next[[i]])))

        # fix the niche case where there are two mutation overlaps in one individual
        o_dups <- which(duplicated(t.overlap.cols) | duplicated(t.overlap.cols, fromLast = TRUE))
        if(length(o_dups) > 0){
          o_inds <- unique(t.overlap.cols[o_dups])
          rm_cols <- numeric(0)

          # for each dup, usually only one
          for(q in 1:length(o_inds)){
            toc <- o_inds[q]
            t_dup_fill <- as.matrix(x.next[[i]][,..toc])
            t_dup_fill[t.overlaps[ind == toc,]$row] <- ifelse(t_dup_fill[t.overlaps[ind == toc,]$row] == 0, 1, 0)
            data.table::set(overlap_fill, i = as.integer(1:nrow(overlap_fill)), j = which(t.overlap.cols == toc)[1], t_dup_fill)
            rm_cols <- c(rm_cols, which(t.overlap.cols == toc)[-1])
          }

          overlap_fill <- overlap_fill[,-..rm_cols]
          t.overlap.cols <- t.overlap.cols[-rm_cols]
        }

        # set
        data.table::set(x.next[[i]],
                        i = as.integer(1:nrow(x.next[[i]])),
                        j = as.integer(t.overlap.cols),
                        overlap_fill)
      }

      mut.eff <- mut.eff[-unique(mut_info$row[overlap_i_in_mut]),]
      mut_drop <- which(mut_info_index$row %in% mut_info$row[overlap_i_in_mut])
      mut.x <- mut.x[-mut_drop,, drop = FALSE]
      mut_info_index <- mut_info_index[-mut_drop,]
    }

    # split into pops
    empties <- which(pops == 0)
    mut.x <- split(data.table::transpose(mut.x), rep(1:npops, pops))
    mut.x <- lapply(mut.x, data.table::transpose)
    if(length(empties) > 0){
      rep.mut.x <- vector("list", npops)
      rep.mut.x[-empties] <- mut.x
      mut.x <- rep.mut.x; rm(rep.mut.x)
    }

    for(i in 1:npops){
      if(!is.null(x.next[[i]])){
        x.next[[i]] <- rbind(x.next[[i]], mut.x[[i]])
      }
    }
    effects <- rbind(effects, mut.eff)
    mut_info_index$chrs <- uf[mut_info_index$chrs]
    mut_info_index$new <- TRUE
    meta$new <- FALSE
    mut_info_index$row <- NULL
    colnames(mut_info_index) <- colnames(meta)
    meta.next <- rbind(meta, mut_info_index)



    return(list(x = x.next, meta.next = meta.next, effects = effects))
  }

  update_phenotypes <- function(genotypes, effects, h, h.av, fitnesses, model = NULL){
    for(j in 1:length(genotypes)){
      if(is.null(model)){
        if(additive){
          pa <- get.pheno.vals(genotypes[[j]], effect.sizes = effects[,j],
                               h = h,
                               hist.a.var = h.av[j],
                               phased = T,
                               fitnesses = fitnesses)
        }
        else{
          pa <- get.pheno.vals(genotypes, effect.sizes = meta[,c("effect_0", "effect_1", "effect_2")],
                               h = h,
                               hist.a.var = h.av[j],
                               phased = T,
                               fitnesses = fitnesses)
        }

      }
      else{
        if(class(model) == "ranger"){
          suppressWarnings(pa <- fetch_phenotypes_ranger(genotypes[[j]], model, h, h.av[j]))
        }
      }

      BVs[[j]] <- pa$a
      phenotypes[[j]] <- pa$p
    }

    return(list(BV = BVs, phenotypes = phenotypes))
  }

  #========================prep===========================
  if(verbose){cat("Initializing...\n")}
  if(!is.list(genotypes)){
    genotypes <- list(genotypes)
  }
  genotypes <- lapply(genotypes, data.table::as.data.table)
  if(length(unique(unlist(lapply(genotypes, nrow)))) != 1){
    stop("All genotype matrices must have the same number of loci.\n")
  }
  npops <- length(genotypes)

  if(is.null(model)){
    # thin genotypes if possible
    zeros <- numeric()
    if("effect" %in% colnames(meta)){
      additive <- TRUE
      if("effect" %in% colnames(meta) & any(meta$effect != 0) & any(meta$effect == 0)){
        if(fitnesses){
          stop("Cannot supply fitnesses with additive effects denoted by a single effect per allele. Supply three columns instead.\n")
        }
        zeros <- which(meta$effect == 0)

        effects <- as.matrix(meta[,"effect", drop = FALSE])
      }
    }
    else if(all(c("effect_0", "effect_1", "effect_2") %in% colnames(meta))){
      additive <- FALSE
      sum_effects <- rowSums(meta[,c("effect_0", "effect_1", "effect_2")])
      zeros <- which(sum_effects == ifelse(fitnesses, 3, 0))
    }
    else if(!is.null(effects)){
      if(!is.matrix(effects)){
        effects <- matrix(effects, ncol = 1)
      }
      if(nrow(effects) != nrow(genotypes[[1]])){
        stop("The number of locus effects is not equal to the number of loci.")
      }
      zeros <- which(rowSums(effects) == 0)
      additive <- TRUE
    }
    else{
      stop("Cannot locate effects in meta.\n")
    }

    if(thin & length(zeros) != 0 & length(zeros) != nrow(meta)){
      nmeta <- meta[-zeros,]
      # if chromosomes are missing, need to remove them from
      # the chr.length vector
      missing.chrs <- which(!unique(meta[,1]) %in% unique(nmeta[,1])) # which are now missing?
      if(length(missing.chrs) > 0){
        chr.length <- chr.length[-missing.chrs]
      }
      meta <- nmeta
      rm(nmeta)
      genotypes <- lapply(genotypes, function(x) x[-zeros,])
      effects <- effects[-zeros,,drop = FALSE]
    }
  }



  if(is.null(phenotypes) | is.null(BVs)){
    need_phenos <- is.null(phenotypes)
    need_BVs <- is.null(BVs)

    for(i in 1:npops){
      if(is.null(model)){ # fectch from provided effects
        if(additive){
          p <- get.pheno.vals(genotypes[[i]], effects[,i], h = h, phased = T, fitnesses = fitnesses)
        }
        else{
          p <- get.pheno.vals(genotypes[[i]], meta[,c("effect_0", "effect_1", "effect_2")], h = h, phased = T, fitnesses = fitnesses)
        }

      }
      else{ # fetch from model
        if(class(model) == "ranger"){
          p <- fetch_phenotypes_ranger(genotypes[[i]], model, h, phased = T)
        }
      }
      if(need_phenos){
        if(is.null(phenotypes)){
          phenotypes <- vector("list", npops)
        }
        phenotypes[[i]] <- p$p
      }
      if(need_BVs){
        if(is.null(BVs)){
          BVs <- vector("list", npops)
        }
        BVs[[i]] <- p$a
      }
      rm(p)
    }
  }

  colnames(effects) <- paste0("effects_", 1:npops)

  #================print out initial conditions, initialize final steps, and run===========
  #starting optimal phenotype, which is the starting mean additive genetic value.
  if(is.null(starting.surv.opt)){
    opt <- unlist(lapply(BVs, mean)) #optimum phenotype
    if(fitnesses){opt <- rep(1, length(opt))}
  }
  else if(length(starting.surv.opt) == 1){
    opt <- rep(starting.surv.opt, length(genotypes))
  }
  else{
    opt <- starting.surv.opt
  }


  if(verbose){
    cat("\n\n===============done===============\n\nStarting parms:\n\tstarting optimum phenotype:", opt,
        "\n\tmean phenotypic value:", paste0(unlist(lapply(phenotypes, mean)), collapse = ", "),
        "\n\taddative genetic variance:", paste0(unlist(lapply(BVs, var)), collapse = ", "),
        "\n\tphenotypic variance:", paste0(unlist(lapply(phenotypes, var)), collapse = ", "), "\n\th:", h, "\n")
  }

  #make output matrix and get initial conditions
  out <- array(NA, c(gens + 1, 8, npops))
  colnames(out) <- c("N", "mu_phenotypes", "mu_BVs", "opt", "diff", "var_BVs", "stochastic_opt", "gen")
  N <- unlist(lapply(genotypes, ncol))/2 #initial pop size
  h.av <- unlist(lapply(BVs, var)) #get the historic addative genetic variance.
  h.pv <- unlist(lapply(phenotypes, var)) #historic phenotypic variance.

  out[1,,] <- t(as.matrix(data.frame(N,
                                     unlist(lapply(phenotypes, mean)),
                                     unlist(lapply(BVs, mean)),
                                     opt, 0, h.av, opt, 0))) #add this and the mean initial additive genetic variance
  if(plot_during_progress){
    pdat <- apply(out, 3, function(x) reshape2::melt(as.data.frame(x), id.vars = "gen"))
    pdat <- data.table::rbindlist(pdat, idcol = "pop")
    colnames(pdat) <- c("pop", "Generation", "var", "val")
    pdat$pop <- as.factor(pdat$pop)
    pdat <- na.omit(pdat)
    print(ggplot2::ggplot(pdat, ggplot2::aes(Generation, val)) + ggplot2::geom_point(na.rm = T) +
            ggh4x::facet_grid2(pop~var, scales = "free_y", independent = "y") +
            ggplot2::theme_bw() +
            ggplot2::scale_x_continuous(limits = c(1, gens), breaks = scales::pretty_breaks()) +
            ggplot2::theme(strip.placement = "outside", axis.title.y = ggplot2::element_blank(),
                           strip.background = ggplot2::element_blank(),
                  strip.text = ggplot2::element_text(size = 11)))
  }

  #initialize matrix to return allele frequencies if requested.
  if(print.all.freqs){
    num <- c(NA, 0)
    a.fqs <- array(num[1], c(nrow(meta), gens + 1, npops))
    a.fqs[,1,] <- unlist(lapply(genotypes, function(x) rowSums(x)/ncol(x)))
  }



  if(!is.list(survival.function)){
    survival.function <- list(survival.function)[rep(1, npops)]
  }
  if(!is.list(rec.dist)){
    rec.dist <- list(rec.dist)[rep(1, npops)]
  }
  if(!is.list(selection.shift.function)){
    selection.shift.function <- list(selection.shift.function)[rep(1, npops)]
  }

  if(length(K_thin_post_surv) == 1){
    K_thin_post_surv <- rep(K_thin_post_surv, npops)
  }

  if(length(chr.length) == 1){
    chr.length <- rep(chr.length, length(unique(meta[,1])))
  }

  if(length(mutation) == 1){
    mutation <- rep(mutation, length(genotypes))
  }

  if(is.list(mutation.effect.function)){
    if(length(mutation.effect.function) != npops){
      stop("Either a single mutation effect function or one for each population must be provided.\n")
    }
  }

  if(length(var.theta) == 1){
    var.theta <- rep(var.theta, npops)
  }



  if(verbose){
    cat("\nBeginning run...\n\n================================\n\n")
  }

  # init thinned tracker, which won't be updated if not thinning fixed loci
  track_thinned_afs <- thin_fixed & print.all.freqs & print.all.thinned.freqs
  if(track_thinned_afs){
    thinned_a.fqs <- array(NA, c(0, ncol(a.fqs), npops))
    thinned_effects <- matrix(NA, 0, ncol(effects))
  }
  else{
    thinned_a.fqs <- NULL
    thinned_effects <- NULL
  }

  #================loop through each additional gen, doing selection, survival, and fisher sampling of survivors====

  for(i in 2:(gens+1)){

    if(thin_fixed){
      fixed <- lapply(genotypes, function(z) rowSums(z) == 0)
      fixed <- matrix(unlist(fixed), ncol = length(genotypes))
      fixed <- which(rowSums(fixed) == ncol(fixed))

      if(length(fixed) > 0){
        genotypes <- lapply(genotypes, function(z) z[-fixed,])
        meta <- meta[-fixed,,drop=FALSE]

        # pull thinned loci out of a.fq
        if(print.all.freqs){

          if(track_thinned_afs){
            new_thinned_a.fqs <- array(NA, c(nrow(thinned_a.fqs) + length(fixed), ncol(thinned_a.fqs), npops))
            if(nrow(thinned_a.fqs) > 0){new_thinned_a.fqs[1:nrow(thinned_a.fqs),,] <- thinned_a.fqs}
            new_thinned_a.fqs[(nrow(thinned_a.fqs)+1):nrow(new_thinned_a.fqs),,] <- a.fqs[fixed,,]
            thinned_a.fqs <- new_thinned_a.fqs
            rm(new_thinned_a.fqs)
            thinned_effects <- rbind(thinned_effects, effects[fixed,,drop=FALSE])
          }
          a.fqs <- a.fqs[-fixed,,]
        }

        effects <- effects[-fixed,,drop=FALSE]
      }
    }

    gen_res <- vector("list", npops)

    # get the optimum phenotype(s) this gen
    t.opt <- opt
    for(j in 1:npops){
      if(!fitnesses){t.opt[j] <- rnorm(1, opt[j], var.theta[j])}
      else{t.opt[j] <- 1}
    }

    # if(length(unique(unlist(lapply(genotypes, nrow)))) != 1){
    #   browser()
    # }

    for(j in 1:npops){
      if(!is.null(genotypes[[j]])){
        genotypes[[j]] <- one_gen(genotypes = genotypes[[j]],
                                phenotypes = phenotypes[[j]],
                                BVs = BVs[[j]],
                                effects = effects[,j],
                                opt = t.opt[j],
                                fitnesses = ifelse(isFALSE(fitnesses), FALSE, fitnesses[[j]]),
                                survival.function = survival.function[[j]],
                                K_thin_post_surv = K_thin_post_surv[j],
                                meta = meta,
                                rec.dist = rec.dist[[j]],
                                chr.length = chr.length,
                                do.sexes = do.sexes,
                                h.av = h.av[j],
                                model = model,
                                selection.shift.function = selection.shift.function[[j]],
                                mutation = mutation[j],
                                pass_surv_genos = ifelse(i == gens + 1 & sampling_point == "parents",
                                                         TRUE, FALSE)
                                )
      }
    }

    if(i == gens + 1 & sampling_point == "parents"){
      final_genotypes <- purrr::map(genotypes, "final_genotypes")
      final_meta <- purrr::map(genotypes, "final_meta")
      final_effects <- purrr::map(genotypes, "final_effects")
      genotypes <- purrr::map(genotypes, "genotypes")

      pa <- update_phenotypes(genotypes = genotypes, effects = effects, h = h, h.av = h.av,
                              fitnesses = fitnesses, model = model)

      final_BVs <- pa$BV; final_phenotypes <- pa$phenotypes

      final_genotypes <- lapply(final_genotypes, function(z){
        if(!is.data.table(z)){
          return(NULL)
        }
        else{
          return(z)
        }
      })
    }

    # replace NAs with NULLs, done this way to prevent list element removal
    genotypes <- lapply(genotypes, function(z){
      if(!is.data.table(z)){
        return(NULL)
      }
      else{
        return(z)
      }
    })

    if(all(unlist(lapply(genotypes, is.null)))){
      warning("All populations went extinct prior to designated number of generations.\n")
      res <- list(run_vars = out,
                  effects = effects, thinned_a.fqs = thinned_a.fqs, thinned_effects = thinned_effects,
                  meta = meta)
      if(print.all.freqs){
        res <- c(res, list(a.fqs = a.fqs))
      }

      return(res)
    }


    #====mutation======
    if(any(mutation > 0)){
      muts <- do_mutation(genotypes,
                          chr.length = chr.length, mutation = mutation,
                          meta = meta, uf = unique(meta[,1]))
      genotypes <- muts$x
      meta <- muts$meta.next
      meta$new <- NULL
      effects <- muts$effects

      # if tracking allele frequencies, add the new loci
      if(print.all.freqs & any(muts$meta.next$new)){
        new_a.fqs <- array(NA, c(nrow(meta), ncol(a.fqs), npops))
        new_a.fqs[1:nrow(a.fqs), 1:ncol(a.fqs),] <- a.fqs
        a.fqs <- new_a.fqs
        rm(new_a.fqs)
      }

      rm(muts); gc(FALSE)
    }

    #========migration================
    if(sampling_point == "offspring" & i == gens + 1){
      final_genotypes <- genotypes
      final_meta <- meta
      final_effects <- effects
      pa <- update_phenotypes(genotypes = genotypes, effects = effects, h = h, h.av = h.av,
                              fitnesses = fitnesses, model = model)

      final_BVs <- pa$BV; final_phenotypes <- pa$phenotypes
    }

    if(!isFALSE(migration)){
      dest <- vector("list", npops)
      new_genotypes <- dest

      # if(any(unlist(lapply(genotypes, is.null)))){browser()}


      # assign destinations according to migration probabilities
      for(k in 1:npops){
        if(!is.null(genotypes[[k]])){
          dest[[k]] <- rmultinom(1, ncol(genotypes[[k]])/2, migration[,k])
          dest[[k]] <- rep(1:nrow(dest[[k]]), dest[[k]])
          dest[[k]] <- sample(dest[[k]], length(dest[[k]]), FALSE)
        }
      }

      # move
      nsizes <- table(unlist(dest))
      for(k in 1:npops){ # into k
        new_genotypes[[k]] <- data.table::as.data.table(matrix(0, ncol = nsizes[k]*2, nrow = max(unlist(lapply(genotypes, nrow)))))
        prog <- 0

        if(ncol(new_genotypes[[k]]) > 0){
          for(j in 1:npops){ # from j
            if(!is.null(genotypes[[j]])){
              movers <- which(rep(dest[[j]], each = 2) == k)

              if(length(movers) > 0){
                data.table::set(new_genotypes[[k]], i = 1:nrow(new_genotypes[[k]]), j = (prog + 1):(prog + sum(dest[[j]] == k)*2),
                                genotypes[[j]][,..movers])

                prog <- (prog + (sum(dest[[j]] == k)*2))
              }
            }
          }
        }
        else{
          new_genotypes[[k]] <- 0
        }
      }

      genotypes <- new_genotypes
      rm(new_genotypes);gc(verbose = FALSE)
    }

    #========update phenotypes, etc==============
    # update genotypes, phenotypes, BVs, optima
    pa <- update_phenotypes(genotypes = genotypes, effects = effects, h = h, h.av = h.av,
                            fitnesses = fitnesses, model = model)

    BVs <- pa$BV; phenotypes <- pa$phenotypes

    if(sampling_point == "migrants" & i != gens + 1){
      final_BVs <- BVs; final_phenotypes <- phenotypes
      final_meta <- meta
      final_genotypes <- genotypes
      final_effects <- effects
    }

    #========progress report==========
    out[i,1,] <- unlist(lapply(genotypes, function(x) sum(ncol(x)/2)))
    out[i,2,] <- unlist(lapply(phenotypes, mean))
    out[i,3,] <- unlist(lapply(BVs, mean))
    out[i,4,] <- opt
    out[i,5,] <- opt - out[i,3,]
    out[i,6,] <- unlist(lapply(BVs, var))
    out[i,7,] <- t.opt
    out[i,8,] <- i - 1

    if(verbose){
      cat("gen:", i - 1,
          "\tfixed_opt:", paste0(round(out[i,4,],3), collapse = ","),
          "\tstoch_opt", paste0(round(out[i,7,],3), collapse = ","),
          "\tmean(BVs):", paste0(round(out[i,3,],3), collapse = ","),
          "\tvar(BVs):", paste0(round(out[i,6,],3), collapse = ","),
          "\tlag:", paste0(round(out[i,4,],3) - round(out[i,3,],3), collapse = ","),
          "\tN:", paste0(out[i,1,], collapse = ","),"\n")
    }

    if(plot_during_progress){
      pdat <- apply(out, 3, function(x) reshape2::melt(as.data.frame(x), id.vars = "gen"))
      pdat <- data.table::rbindlist(pdat, idcol = "pop")
      colnames(pdat) <- c("pop", "Generation", "var", "val")
      pdat$pop <- as.factor(pdat$pop)
      pdat <- na.omit(pdat)
      print(ggplot2::ggplot(pdat, ggplot2::aes(Generation, val)) + ggplot2::geom_point(na.rm = T) +
              ggh4x::facet_grid2(pop~var, scales = "free_y", independent = "y") +
              ggplot2::theme_bw() +
              ggplot2::scale_x_continuous(limits = c(1, gens), breaks = scales::pretty_breaks()) +
              ggplot2::theme(strip.placement = "outside", axis.title.y = ggplot2::element_blank(),
                             strip.background = ggplot2::element_blank(),
                             strip.text = ggplot2::element_text(size = 11)))
    }

    if(print.all.freqs){
      a.fqs[,i,] <- unlist(lapply(genotypes, function(z) rowSums(z)/ncol(z)))
    }

    #==========update selection optima===========
    if(!fitnesses){
      opt[j] <- selection.shift.function[[j]](opt[j], iv = sqrt(h.av[j]))
    }
  }

  if(!print.all.freqs){
    num <- c(NA, 0)
    a.fqs <- array(num[1], c(max(unlist(lapply(genotypes, nrow))), 1, npops))
    for(i in 1:npops){
      if(!is.null(genotypes[[i]])){
        a.fqs[,1,i] <- rowSums(genotypes[[i]])/ncol(genotypes[[i]])
      }
    }
  }

  return(list(run_vars = out, genotypes = final_genotypes, phenotypes = final_phenotypes, BVs = final_BVs, a.fqs = a.fqs,
              effects = final_effects, thinned_a.fqs = thinned_a.fqs, thinned_effects = thinned_effects,
              meta = final_meta))
}


#=======function to pull parameter values to run from a joint quantile result==========
#' Sample hyperparameter values from hyperparameter regression/random forest results.
#'
#' Samples hyperparameters from the joint quantile distribution resulting from a regression on ABC results
#' and a random forest.
#'
#' @param n numeric. The number of samples to draw.
#' @param reg regression results from \code{\link{hyperparameter_regression_on_ABC}}.
#' @return A data.frame containing the sampled hyperparameter values.
#'
#' @author William Hemstrom
#' @export
sample_joint_quantile <- function(n, reg){
  res <- reg$joint_quantile_long[sample(1:nrow(reg$joint_quantile_long), n, T, reg$joint_quantile_long$norm_joint_quantile), c(1, 4)]
  colnames(res) <- c(reg$independent, reg$dependent)
  return(res)
}

#=======distribution functions for survival simulations===============================
#' Calculate population size after logistic growth
#'
#' Given a starting population size, find the estiamted size after one generation of
#' logistic growth.
#'
#' @param n numeric. Starting population size.
#' @param K numeric. Carrying capacity.
#' @param r numeric. Logistic growth rate.
#'
#' @export
logistic_growth <- function(n, K, r, ...){
  return((K*n*exp(r))/(K + n*(exp(r*1) - 1)))
}

#' @export
BL_growth <- function(n, B, ...){
  return(B*n)
}

#' Calculate survivability probabilities using a scaled normal distribution.
#'
#' For a given optimum phenotype and surivival variance, calculate survival probabilities from a
#' normal distribution.
#'
#' Surival probabilities are given by \code{\link[stats]{pnorm}}. For values above the optimum phenotype,
#' the surivial is given by 1 - the resulting value. If a maxium survival other than 0.5 is given,
#' the surivivals are scaled such that individuals with the optimum phenotype will have the provided
#' survival probability.
#'
#' @param phenotypes numeric. Vector of phenotypes.
#' @param opt_pheno numeric. Optimum phenotype.
#' @param var numeric. Variance of surivial distribution.
#' @param max.survival numeric. Maximum possible survival probability.
#'
#' @export
normal_survival <- function(phenotypes, opt_pheno, var, max.survival = .5, ...){
  if(max.survival > 1){
    warning("Max surivial is above 1, will reset to 1.")
    max.survival <- 1
  }
  # get survival probabilities from normal curve
  x <- pnorm(phenotypes, opt_pheno, sqrt(var))
  x <- ifelse(x > .5, 1 - x, x)

  # scale for target max survival
  scale_factor <- max.survival/.5
  x <- x*scale_factor

  return(x)
}

#' Calculate survival using the Burger and Lynch method
#'
#' Surivial is caluclated using equation 1 from Burger and Lynch 1995.
#'
#' @param phenotypes numeric. Vector of phenotypes of individuals.
#' @param opt_pheno numeric. Optimum phenotype this generation.
#' @param omega numeric. Strength of stabilizing selection. Smaller numbers mean stronger selection.
#'
#' @export
BL_survival <- function(phenotypes, opt_pheno, omega, ...){
  exp(-((phenotypes - opt_pheno)^2)/(2*omega^2))
}

#' @export
prop_survival <- function(prop_high, opt_prop_high, omega){
  return(exp(-((prop_high - opt_prop_high)^2)/(2*omega^2)))
}

#' @export
prop_survival_scaled_omega <- function(prop_high, opt_prop_high, h.av, omega){
  omega <- omega*h.av
  return(exp(-((prop_high - opt_prop_high)^2)/(2*omega^2)))
}

#' Increase the selection optimum by a given percentage of starting variance each generation.
#'
#' @param opt numeric. Initial optimum phenotype.
#' @param iv numeric. Initial variance. The optimum will increase by some portion of this variance.
#' @param slide numeric. Proportion of the initial variance by which the optimum phenotype will slide.
#'
#' @export
optimum_constant_slide <- function(opt, iv, slide = 0.3, ...){
  return(opt + iv*slide)
}

#' @export
optimum_fixed_slide <- function(opt, slide = 0.3, ...){
  return(opt + slide)
}

#' @export
optimum_proportion_constant_slide <- function(opt, iv, slide, ...){
  return(min(1, opt + iv*slide))
}

#' @export
optimum_proportion_fixed_slide <- function(opt, slide, ...){
  return(min(1, opt + opt*slide))
}

#' Increases the selection optimum according to logistic growth.
#'
#' @param opt numeric. Current optimum phenotype.
#' @param init_opt numeric. Initial optimum phenotype
#' @param scale numeric. The proportion of the initial optimum phenotype at which the slide will max out.
#' @param r numeric. The growth rate of the optimum phenotype.
#'
#' @export
optimum_logistic_slide <- function(opt, init_opt, scale, r, ...){
  K <- init_opt*scale

  if(opt < 0){
    opt <- abs(opt)
    K <- -K
    return(-1*((K*opt*eoptp(r))/(K + opt*(eoptp(r*1) - 1))))
  }
  else{
    return((K*opt*eoptp(r))/(K + opt*(eoptp(r*1) - 1)))
  }
}



#==============quantitative genetics models of population trends==========
#' Estimate the average population sizes of a population under selection according to Burger and Lynch 1995.
#'
#' Uses the equations laid out in Burger and Lynch 1995 to estimate the average size of a population
#' overtime as the optimum phenotype changes. While some processes are \emph{modeled} as stochastic,
#' this returns the \emph{mean} population sizes over time, and so is not stochastic itself.
#'
#' If nloci, alpha, and mu are all provided, the stochastic house of cards
#' approximation of genetic variance is used to approximate changes in genetic variance
#' each generation based \emph{only} on changes in effective population size. Otherwise,
#' genetic variance is assumed to not change in each generation and is based solely on the
#' genetic variance calculated from the provided heritability and phenotypes.
#'
#' @param phenotypes numeric. Vector of initial phenotypes.
#' @param h numeric. Heritability of the trait.
#' @param B numeric. Number of offspring each surviving adult has. Not stochastic.
#' @param K numeric. Carrying capacity. Individuals surviving selection will be randomly culled to this number.
#' @param omega numeric. Width of the selection function. Smaller numbers equate to stronger selection around the optimum.
#' @param var.theta numeric. Variance of the stochastic selection optimum.
#' @param k numeric. Degree by which the optimum phenotype increases each generation (where mean(opt) = k*t).
#' @param gens numeric, default Inf. Number of generations to iterate.
#' @param n numeric, default NULL. Initial population size. If null, equal to the number of provided phenotypes.
#' @param nloci numeric, default NULL. Number of effective loci. If nloci, alpha, and mu are all provided, the stochastic house of cards approximation of genetic variance is used to approximate changes in genetic variance each generation based \emph{only} on changes in effective population size.
#' @param alpha numeric, default NULL. Standard deviation of the effect of new mutations. If nloci, alpha, and mu are all provided, the stochastic house of cards approximation of genetic variance is used to approximate changes in genetic variance each generation based \emph{only} on changes in effective population size.
#' @param mu numeric, default NULL. Mutation rate. If nloci, alpha, and mu are all provided, the stochastic house of cards approximation of genetic variance is used to approximate changes in genetic variance each generation based \emph{only} on changes in effective population size.
#'
#' @export
gs_BL <- function(phenotypes, h, K, omega, B, var.theta, k, gens = Inf, n = NULL, nloci = NULL, alpha = NULL, mu = NULL){
  # calculate variances
  p_var <- var(phenotypes)
  g_var <- h * p_var
  e_var <- p_var - g_var

  # standardize such that e_var = 1
  omega_shift_func <- function(gvar1, evar1, omega1, gvar2, evar2){
    omega1_var <- omega1^2
    omega2_var <- (omega1_var*(evar2 + gvar2))/(evar1 + gvar1)
    return(return(sqrt(omega2_var)))
  }

  ## standardize variances
  p_var.n <- p_var/e_var
  g_var.n <- g_var/e_var
  e_var.n <- e_var/e_var
  var.theta <- var.theta/e_var

  ## standardize omega
  omega.n <- omega_shift_func(g_var, e_var, omega, g_var.n, e_var.n)

  ## update
  p_var <- p_var.n
  g_var <- g_var.n
  e_var_original <- e_var
  e_var <- e_var.n
  omega <- omega.n

  ## standardize k
  k <- k/e_var

  # iterate across time to get pop sizes
  t <- 1
  Vs <- (omega^2) + e_var
  s <- g_var/(g_var + Vs)
  mean_phenos <- 0
  V_mean_phenos <- 0
  lambdas <- NA
  opt_pheno <- 0
  g_var_storage <- g_var*e_var_original
  if(is.null(n)){
    n <- length(phenotypes)
  }

  n <- min(n, k_thin)
  out_n <- n

  while(n[t] >= 2 & t <= gens){
    ne <- ((2*B)/((2*B) - 1))*n[t] # BL eq 13.

    # if nloci, alpha, and mu provided, use stochastic house of cards to estimate g_var this gen. Otherwise assume constant.
    if(!is.null(nloci) & !is.null(alpha) & !is.null(mu)){
      g_var <- (4*nloci*mu*ne*alpha^2)/(1 + ((ne*alpha^2)/Vs)) #  BL eq 14.
    }
    g_var_storage <- c(g_var_storage, g_var*e_var_original)

    #equation 7a
    Vgt <- (((((g_var + Vs)^2)/(ne*(g_var + 2*Vs))) + (g_var*var.theta)/(g_var + 2*Vs)) *
      (1 - (1 - s)^(2*t)))
    Vlt <- Vs + g_var + Vgt + var.theta
    Egt <- k*t - (k/s)*(1 - (1 - s)^t)

    # equation 9
    Bo <- B*omega/sqrt(Vlt)
    Lt <- Bo * exp(-((Egt - k*t)^2)/(2*Vlt))

    # growth
    nt <- Lt*n[t]
    out_n <- c(out_n, nt)
    if(nt > K){nt <- K}
    n <- c(n, nt)

    # update
    mean_phenos <- c(mean_phenos, Egt)
    V_mean_phenos <- c(V_mean_phenos, Vgt)
    lambdas <- c(lambdas, Lt)
    opt_pheno <- c(opt_pheno, k*t)
    if(t > 1 & lambdas[t] == lambdas[t + 1]){
      warning("Population sustainable: k is too small to cause extinction.\n")
      break
    }
    t <- t + 1
  }

  out <- data.frame(t = 1:length(n), n = n, mean_pheno = mean_phenos,
                    V_mean_pheno = V_mean_phenos, lambda = lambdas, opt_pheno = opt_pheno, g_var = g_var_storage)
  return(out)
}

#' Simulate population demographics under selection with the Breeder's Equation.
#'
#' Uses several equations from Burger and Lynch 1995 to calculate ne and additive genetic variance,
#' but otherwise simulates purely based on the Breeder's Equation.
#'
#' If nloci, alpha, and mu are all provided, the stochastic house of cards
#' approximation of genetic variance is used to approximate changes in genetic variance
#' each generation based \emph{only} on changes in effective population size. Otherwise,
#' genetic variance is assumed to not change in each generation and is based solely on the
#' genetic variance calculated from the provided heritability and phenotypes.
#'
#' @param phenotypes numeric. Vector of initial phenotypes.
#' @param h numeric. Heritability of the trait.
#' @param B numeric. Number of offspring each surviving adult has. Not stochastic.
#' @param K numeric. Carrying capacity. Individuals surviving selection will be randomly culled to this number.
#' @param omega numeric. Width of the selection function. Smaller numbers equate to stronger selection around the optimum.
#' @param var.theta numeric. Variance of the stochastic selection optimum.
#' @param k numeric. Amount by which the optimum phenotype increases each generation (where mean(opt) = k*t).
#' @param gens numeric, default Inf. Number of generations to iterate.
#' @param n numeric, default NULL. Initial population size. If null, equal to the number of provided phenotypes.
#' @param nloci numeric, default NULL. Number of effective loci. If nloci, alpha, and mu are all provided, the stochastic house of cards approximation of genetic variance is used to approximate changes in genetic variance each generation based \emph{only} on changes in effective population size.
#' @param alpha numeric, default NULL. Standard deviation of the effect of new mutations. If nloci, alpha, and mu are all provided, the stochastic house of cards approximation of genetic variance is used to approximate changes in genetic variance each generation based \emph{only} on changes in effective population size.
#' @param mu numeric, default NULL. Mutation rate. If nloci, alpha, and mu are all provided, the stochastic house of cards approximation of genetic variance is used to approximate changes in genetic variance each generation based \emph{only} on changes in effective population size.
#'
#' @export
gs_breeders <- function(phenotypes, h, B, K,
                        omega, var.theta, k, gens = Inf, n = NULL,
                        nloci = NULL, alpha = NULL, mu = NULL){
  # calculate variances
  p_var <- var(phenotypes)
  g_var <- h * p_var
  e_var <- p_var - g_var

  # iterate across time to get pop sizes
  t <- 1
  Vs <- (omega^2) + e_var
  s <- g_var/(g_var + Vs)
  opt_pheno <- 0
  phenotypes <- phenotypes - mean(phenotypes) # center
  mean_phenos <- 0
  Rs <- 0
  Ss <- 0
  opt <- 0
  g_var_storage <- g_var
  if(is.null(n)){
    n <- length(phenotypes)
  }

  n <- min(n, k_thin)



  while(n[t] >= 2 & t <= gens){

    # stochastic survival
    t_kt <- sum(k*t + rnorm(1, 0, sqrt(var.theta)))
    surv <- rbinom(length(phenotypes), 1, BL_survival(phenotypes, t_kt, omega))

    nsurv <- sum(surv)
    if(nsurv > K){surv[sample(which(surv == 1), nsurv - K)] <- 0; nsurv <- K} # if we are above K, kill some more at random.

    opt <- c(opt, t_kt)
    new_mean_phenos <- mean(phenotypes[which(surv == 1)])
    mean_phenos <- c(mean_phenos, new_mean_phenos)
    n <- c(n, nsurv)


    if(sum(surv) < 2){
      Rs <- c(Rs, NA)
      Ss <- c(Ss, NA)
      g_var_storage <- c(g_var_storage, NA)
      break
    }

    # apply the breeder's equation to get the response to selection (what portion of the change in phenotype is passed?)
    S <- new_mean_phenos - mean(phenotypes)
    R <- h*S

    # next gen's phenotypes



    # if nloci, alpha, and mu provided, use stochastic house of cards to estimate g_var this gen. Otherwise assume constant.
    if(!is.null(nloci) & !is.null(alpha) & !is.null(mu)){
      ne <- ((2*B)/((2*B) - 1))*nsurv # BL eq 13.
      g_var <- (4*nloci*mu*ne*alpha^2)/(1 + ((ne*alpha^2)/Vs)) # BL eq 14.

      phenotypes <- rnorm(nsurv*B, mean(phenotypes) + R, sqrt(g_var + e_var))
    }
    else{
      phenotypes <- rnorm(nsurv*B, mean(phenotypes) + R, sqrt(p_var))
    }
    g_var_storage <- c(g_var_storage, g_var)

    Rs <- c(Rs, R)
    Ss <- c(Ss, S)
    t <- t + 1
  }

  out <- data.frame(t = 1:length(n), n = n, mean_pheno = mean_phenos, opt_pheno = opt,
                    response_to_selection = Rs, selection_differential = Ss, g_var = g_var_storage)
  return(out)
}

#' @export
gs_outliers <- function(genotypes, pvals, scores, gens,
                        opt_prop_slide = function(opt) min(1, opt + opt*.5),
                        surivival.function = function(prop_high, opt_prop_high) prop_survival(prop_high, opt_prop_high, .1),
                        growth.function = function(n) logistic_growth(n, 500, .5),
                        n = NULL, pcrit = 1*10^-5,
                        K_thin_post_surv = NULL){
  # survival subfunciton
  surv_func <- surivival.function
  # maf function
  maf_func <- function(alleles){
    rowSums(alleles)/(2*ncol(alleles))
  }

  #===============initialize=============
  # find outliers and grab the genotypes
  outliers <- which(pvals <= pcrit)
  nloci <- length(outliers)
  out.meta <- data.frame(out_num = 1:length(outliers), pval = pvals[outliers])
  genotypes <- genotypes[outliers,]

  # flip everything to track the allele that has a positive effect on phenotype
  flip <- which(scores[outliers] < 0)
  genotypes[flip,] <- (genotypes[flip,]*-1) + 1

  # get ready to track things for sims
  alleles<- t(convert_2_to_1_column(genotypes))

  if(is.null(n)){
    n <- ncol(genotypes)/2
  }

  # fill out starting conditions
  out <- data.frame(n = numeric(gens + 1), mean_prop_high = numeric(gens + 1), opt_prop_high = numeric(gens + 1), var_prop_high = numeric(gens + 1), gen = 1:(gens + 1))
  out$n[1] <- n
  maf <- maf_func(alleles)
  prop_high <- colSums(alleles/(nloci*2))
  out$mean_prop_high[1] <- mean(prop_high)
  out$opt_prop_high[1] <- out$mean_prop_high[1]
  out$var_prop_high[1] <- var(prop_high)
  maf_output <- as.data.frame(matrix(0, gens, nloci))
  maf_output[1,] <- maf
  colnames(maf_output) <- paste0("locus_", 1:nloci)

  #===============run sim==================
  for(i in 2:(gens + 1)){
    # survival
    surv <- rbinom(out$n[i - 1], 1, surv_func(prop_high, out$opt_prop_high[i - 1], h.av = out$var_prop_high[1]))

    # if doing carrying capacity on BREEDERS, not offspring, thin to K here.
    if(!is.null(K_thin_post_surv)){
      if(sum(surv) > K_thin_post_surv){
        surv[which(surv == 1)][sample(sum(surv), K_thin_post_surv, F)] <- 0
      }
    }

    # growth
    out$n[i] <- round(growth.function(sum(surv)))

    if(out$n[i] < 2){
      out <- out[1:i,]
      maf_output <- maf_output[1:i,]
      break
    }

    # alleles for next generation
    maf <- maf_func(alleles[,which(surv == 1), drop = F])
    alleles <- rbinom(out$n[i]*nloci, 2, maf)
    alleles <- matrix(alleles, nloci)

    # update results
    prop_high <- colSums(alleles)/(nloci*2)
    out$mean_prop_high[i] <- mean(prop_high)
    out$var_prop_high[i] <- var(prop_high)
    out$opt_prop_high[i] <- opt_prop_slide(out$opt_prop_high[i - 1], iv = out$var_prop_high[1])
    maf_output[i,] <- maf
  }

  #==============return==========
  return(list(demographics = out, maf = maf_output, final_genotypes = alleles))
}
