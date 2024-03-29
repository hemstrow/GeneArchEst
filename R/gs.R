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
               growth.function = function(n) logistic_growth(n, 500, 2),
               survival.function = function(phenotypes, opt_pheno, ...) BL_survival(phenotypes, opt_pheno, omega = 1),
               selection.shift.function = function(opt, iv) optimum_constant_slide(opt, iv, 0.3),
               rec.dist = function(n) rpois(n, lambda = 1),
               var.theta = 0,
               plot_during_progress = FALSE,
               facet = "group",
               do.sexes = TRUE,
               fitnesses = FALSE,
               init = FALSE,
               thin = TRUE,
               verbose = FALSE,
               print.all.freqs = F,
               model = NULL,
               stop_if_no_variance = FALSE,
               K_thin_post_surv = NULL){

  #========================prep===========================
  if(verbose){cat("Initializing...\n")}
  genotypes <- data.table::as.data.table(genotypes)

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

      }
    }
    else if(all(c("effect_0", "effect_1", "effect_2") %in% colnames(meta))){
      additive <- FALSE
      sum_effects <- rowSums(meta[,c("effect_0", "effect_1", "effect_2")])
      zeros <- which(sum_effects == ifelse(fitnesses, 3, 0))
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
      genotypes <- genotypes[-zeros,]
    }
  }



  if(is.null(phenotypes) | is.null(BVs)){
    if(is.null(model)){ # fectch from provided effects
      if(additive){
        p <- get.pheno.vals(genotypes, meta$effect, h = h, phased = T, fitnesses = fitnesses)
      }
      else{
        p <- get.pheno.vals(genotypes, meta[,c("effect_0", "effect_1", "effect_2")], h = h, phased = T, fitnesses = fitnesses)
      }

    }
    else{ # fetch from model
      if(class(model) == "ranger"){
        p <- fetch_phenotypes_ranger(genotypes, model, h, phased = T)
      }
    }
    if(is.null(phenotypes)){
      phenotypes <- p$p
    }
    if(is.null(BVs)){
      BVs <- p$a
    }
    rm(p)
  }

  #================print out initial conditions, initialize final steps, and run===========
  #starting optimal phenotype, which is the starting mean additive genetic value.
  opt <- mean(BVs) #optimum phenotype
  if(fitnesses){opt <- 1}

  if(verbose){
    cat("\n\n===============done===============\n\nStarting parms:\n\tstarting optimum phenotype:", opt,
        "\n\tmean phenotypic value:", mean(phenotypes), "\n\taddative genetic variance:", var(BVs), "\n\tphenotypic variance:", var(phenotypes), "\n\th:", h, "\n")
  }

  #make output matrix and get initial conditions
  out <- matrix(NA, nrow = gens + 1, ncol = 8)
  colnames(out) <- c("N", "mu_phenotypes", "mu_BVs", "opt", "diff", "var_BVs", "stochastic_opt", "gen")
  N <- ncol(genotypes)/2 #initial pop size
  h.av <- var(BVs) #get the historic addative genetic variance.
  h.pv <- var(phenotypes) #historic phenotypic variance.

  out[1,] <- c(NA, mean(phenotypes), mean(BVs), opt, 0, h.av, opt, 1) #add this and the mean initial additive genetic variance
  if(plot_during_progress){
    library(ggplot2)
    pdat <- reshape2::melt(out)
    colnames(pdat) <- c("Generation", "var", "val")
    if(!fitnesses){
      ranges <- data.frame(var = c("N", "mu_phenotypes", "mu_BVs", "opt", "diff"),
                                        ymin = c(0, out[1,2]*2, out[1,3]*2, out[1,4]*2, -10),
                                        ymax = c(out[1,1]*1.05, 0, 0, 0, 10))
    }
    else{
      ranges <- data.frame(var = c("N", "mu_phenotypes", "mu_BVs", "opt", "diff"),
                           ymin = c(0, 0, 0, 0, 0),
                           ymax = c(out[1,1]*1.05, 1, 1, 1, 1))
    }
    pdat <- merge(pdat, ranges, by = "var")
    print(ggplot(pdat, aes(Generation, val)) + geom_point(na.rm = T) +
            facet_wrap(~var, ncol = 1, scales = "free_y", strip.position = "left") +
            geom_blank(aes(y = ymin)) +
            geom_blank(aes(y = ymax)) +
            theme_bw() + xlim(c(0, max(pdat$Generation))) +
            theme(strip.placement = "outside", axis.title.y = element_blank(), strip.background = element_blank(),
                  strip.text = element_text(size = 11)))
  }

  #initialize matrix to return allele frequencies if requested.
  if(print.all.freqs){
    a.fqs <- matrix(0, nrow(meta), gens + 1)
    a.fqs[,1] <- rowSums(x)/ncol(x)
  }

  #================loop through each additional gen, doing selection, survival, and fisher sampling of survivors====

  if(verbose){
    cat("\nBeginning run...\n\n================================\n\n")
  }

  offn <- N
  for(i in 2:(gens+1)){
    #=========survival====
    # get the optimum phenotype this gen
    if(!fitnesses){t.opt <- rnorm(1, opt, var.theta)}
    else{t.opt <- 1}

    if(length(unique(phenotypes)) == 1 & phenotypes[1] != 1 & stop_if_no_variance){stop("No genetic variance left.\n")}
    #survival:
    s <- rbinom(offn, 1, #survive or not? Number of draws is the pop size in prev gen, survival probabilities are determined by the phenotypic variance and optimal phenotype in this gen.
                survival.function(phenotypes, t.opt))

    #if the population has died out, stop.
    if(sum(s) <= 1){
      break
    }

    # if doing carrying capacity on BREEDERS, not offspring, thin to K here.
    if(!is.null(K_thin_post_surv)){
      if(sum(s) > K_thin_post_surv){
        s[which(s == 1)][sample(sum(s), sum(s) - K_thin_post_surv, F)] <- 0
      }
    }

    #===============save output============

    # update output and report
    out[i,1] <- sum(s)
    out[i,2] <- mean(phenotypes)
    out[i,3] <- mean(BVs)
    out[i,4] <- opt
    out[i,5] <- opt - mean(BVs)
    out[i,6] <- var(BVs)
    out[i,7] <- t.opt
    if(verbose){
      cat("gen:", i,
          "\tf_opt:", round(out[i,4],3),
          "\ts_opt", round(out[i,7],3),
          "\tmean(BVs):", round(out[i,3],3),
          "\tvar(BVs):", round(out[i,6],3),
          "\tlag:", round(out[i,4],3) - round(out[i,3],3),
          "\tN:", out[i,1],"\n")
    }

    # plot
    if(plot_during_progress){
      pdat <- reshape2::melt(out)
      colnames(pdat) <- c("Generation", "var", "val")
      if(!fitnesses){
        ranges <- data.frame(var = c("N", "mu_phenotypes", "mu_BVs", "opt", "diff"),
                             ymin = c(0, out[1,2]*2, out[1,3]*2, out[1,4]*2, -10),
                             ymax = c(out[1,1]*1.05, 0, 0, 0, 10))
      }
      else{
        ranges <- data.frame(var = c("N", "mu_phenotypes", "mu_BVs", "opt", "diff"),
                             ymin = c(0, 0, 0, 0, 0),
                             ymax = c(out[1,1]*1.05, 1, 1, 1, 1))
      }
      pdat <- merge(pdat, ranges, by = "var")
      print(ggplot(pdat, aes(Generation, val)) + geom_line(na.rm = T) + geom_point(na.rm = T) +
              facet_wrap(~var, ncol = 1, scales = "free_y", strip.position = "left") +
              geom_blank(aes(y = ymin)) +
              geom_blank(aes(y = ymax)) +
              theme_bw() + xlim(c(0, max(pdat$Generation))) +
              theme(strip.placement = "outside", axis.title.y = element_blank(), strip.background = element_blank(),
                    strip.text = element_text(size = 11)))
    }


    #add allele frequencies if requested
    if(print.all.freqs){
      a.fqs[,i] <- rowSums(genotypes)/ncol(genotypes)
    }


    #===============figure out next generation===============
    #what is the pop size after growth?
    offn <- round(growth.function(out[i,1]))


    #make a new x with the survivors
    genotypes <- genotypes[, .SD, .SDcols = which(rep(s, each = 2) == 1)] #get the gene copies of survivors

    #=============do random mating, adjust selection, get new phenotype scores, get ready for next gen====
    genotypes <- rand.mating(x = genotypes, N.next = offn, meta = meta, rec.dist = rec.dist, chr.length, do.sexes)
    # check that the pop didn't die due to every individual being the same sex (rand.mating returns NULL in this case.)
    if(is.null(genotypes)){
      break
    }

    #get phenotypic/genetic values
    if(is.null(model)){
      if(additive){
        pa <- get.pheno.vals(genotypes, effect.sizes = meta$effect,
                             h = h,
                             hist.a.var = h.av,
                             phased = T,
                             fitnesses = fitnesses)
      }
      else{
        pa <- get.pheno.vals(genotypes, effect.sizes = meta[,c("effect_0", "effect_1", "effect_2")],
                             h = h,
                             hist.a.var = h.av,
                             phased = T,
                             fitnesses = fitnesses)
      }

    }
    else{
      if(class(model) == "ranger"){
        suppressWarnings(pa <- fetch_phenotypes_ranger(genotypes, model, h, h.av))
      }
    }

    BVs <- pa$a
    phenotypes <- pa$p


    #adjust selection optima
    if(!fitnesses){opt <- selection.shift.function(opt, iv = sqrt(h.av))}

    gc()
  }

  #prepare stuff to return
  out[,"gen"] <- 1:nrow(out)
  out <- out[-nrow(out),, drop = F]

  if(print.all.freqs){
    a.fqs <- cbind(meta, a.fqs, stringsAsFactors = F)
    out <- list(summary = out, frequencies = a.fqs)
  }

  return(list(run_vars = out, genotypes = genotypes, phenotypes = phenotypes, BVs = BVs))
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
