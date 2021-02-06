#=======function to do growth and selection=======
# note: init means are we simply intiating a population under selection. We'll need to keep all markers if so.
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
               init = F,
               verbose = T,
               print.all.freqs = F,
               model = NULL, K_thin_post_surv = NULL){
  if(verbose){cat("Initializing...\n")}
  genotypes <- data.table::as.data.table(genotypes)

  if(is.null(phenotypes) | is.null(BVs)){
    if(is.null(model)){ # fectch from provided effects
      p <- get.pheno.vals(genotypes, meta$effect, h = h, phased = T)
    }
    else{ # fetch from model
      if(class(model) == "ranger"){
        p <- fetch_phenotypes_ranger(genotypes, model, h)
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

  #================print out initial conditions, intiallize final steps, and run===========
  #starting optimal phenotype, which is the starting mean addative genetic value.
  opt <- mean(BVs) #optimum phenotype

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

  out[1,] <- c(N, mean(phenotypes), mean(BVs), opt, 0, h.av, opt, 1) #add this and the mean initial additive genetic variance
  if(plot_during_progress){
    library(ggplot2)
    pdat <- reshape2::melt(out)
    colnames(pdat) <- c("Generation", "var", "val")
    ranges <- data.frame(var = c("N", "mu_phenotypes", "mu_BVs", "opt", "diff"),
                         ymin = c(0, out[1,2]*2, out[1,3]*2, out[1,4]*2, -10),
                         ymax = c(out[1,1]*1.05, 0, 0, 0, 10))
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

  for(i in 2:(gens+1)){
    #=========survival====
    # get the optimum phenotype this gen
    t.opt <- rnorm(1, opt, var.theta)

    #survival:
    s <- rbinom(out[(i-1),1], 1, #survive or not? Number of draws is the pop size in prev gen, surival probabilities are determined by the phenotypic variance and optimal phenotype in this gen.
                survival.function(phenotypes, t.opt))
    #if the population has died out, stop.
    if(sum(s) <= 1){
      break
    }

    # if doing carrying capacity on BREEDERS, not offspring, thin to K here.
    if(!is.null(K_thin_post_surv)){
      if(sum(s) > K_thin_post_surv){
        s[which(s == 1)][sample(sum(s), K_thin_post_surv, F)] <- 0
      }
    }

    #what is the pop size after growth?
    out[i,1] <- round(growth.function(sum(s)))

    #make a new x with the survivors
    genotypes <- genotypes[, .SD, .SDcols = which(rep(s, each = 2) == 1)] #get the gene copies of survivors

    #=============do random mating, adjust selection, get new phenotype scores, get ready for next gen====
    y <- rand.mating(x = genotypes, N.next = out[i,1], meta = meta, rec.dist = rec.dist, chr.length, do.sexes)
    # check that the pop didn't die due to every individual being the same sex (rand.mating returns NULL in this case.)
    if(is.null(y)){
      break
    }
    else{
      genotypes <- y
      rm(y)
    }

    #get phenotypic/genetic values
    if(is.null(model)){
      pa <- get.pheno.vals(genotypes, effect.sizes = meta$effect,
                           h = h,
                           hist.a.var = h.av,
                           phased = T)
    }
    else{
      if(class(model) == "ranger"){
        suppressWarnings(pa <- fetch_phenotypes_ranger(genotypes, model, h, h.av))
      }
    }

    BVs <- pa$a
    phenotypes <- pa$p


    #adjust selection optima
    opt <- selection.shift.function(opt, iv = sqrt(h.av))

    #save
    out[i,2] <- mean(phenotypes)
    out[i,3] <- mean(BVs)
    out[i,4] <- opt
    out[i,5] <- opt - mean(BVs)
    out[i,6] <- var(BVs)
    out[i,7] <- t.opt
    if(verbose){
      cat("gen:", i-1,
          "\tf_opt:", round(out[i-1,4],3),
          "\ts_opt", round(out[i-1,7],3),
          "\tmean(phenotypes):", round(out[i,2],3),
          "\tmean(BVs):", round(out[i,3],3),
          "\tvar(BVs):", round(var(BVs),3),
          "\tNs:", sum(s),
          "\tN(t+1):", out[i,1],"\n")
    }
    if(plot_during_progress){
      pdat <- reshape2::melt(out)
      colnames(pdat) <- c("Generation", "var", "val")
      ranges <- data.frame(var = c("N", "mu_phenotypes", "mu_BVs", "opt", "diff"),
                           ymin = c(0, out[1,2]*2, out[1,3]*2, out[1,4]*2, -10),
                           ymax = c(out[1,1]*1.05, 0, 0, 0, 10))
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

    gc()
  }

  #prepare stuff to return
  out[,"gen"] <- 1:nrow(out)
  out <- out[-nrow(out),]

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
#' @param k numeric. Proportion of the initial additive genetiv variance by which the optimum phenotype increases each generation (where mean(opt) = k*t).
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
  e_var <- e_var.n
  omega <- omega.n

  ## standardize k
  k <- k*g_var

  # iterate across time to get pop sizes
  t <- 1
  Vs <- (omega^2) + e_var
  s <- g_var/(g_var + Vs)
  mean_phenos <- 0
  V_mean_phenos <- 0
  lambdas <- NA
  opt_pheno <- 0
  if(is.null(n)){
    n <- length(phenotypes)
  }

  while(n[t] >= 2 & t <= gens){
    ne <- ((2*B)/((2*B) - 1))*n[t] # BL eq 13.

    # if nloci, alpha, and mu provided, use stochastic house of cards to estimate g_var this gen. Otherwise assume constant.
    if(!is.null(nloci) & !is.null(alpha) & !is.null(mu)){
      g_var <- (4*nloci*mu*ne*alpha^2)/(1 + ((ne*alpha^2)/Vs)) #  BL eq 14.
    }

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
    if(nt > K){nt <- K}

    # update
    n <- c(n, nt)
    mean_phenos <- c(mean_phenos, Egt)
    V_mean_phenos <- c(V_mean_phenos, Vgt)
    lambdas <- c(lambdas, Lt)
    opt_pheno <- c(opt_pheno, k*t)
    if(t > 1 & lambdas[t] == lambdas[t +1]){
      warning("Population sustainable: k is too small to cause extinction.\n")
      break
    }
    t <- t + 1
  }

  out <- data.frame(t = 1:length(n), n = n, mean_pheno = mean_phenos,
                    V_mean_pheno = V_mean_phenos, lambda = lambdas, opt_pheno = opt_pheno)
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
#' @param k numeric. Proportion of the initial additive genetiv variance by which the optimum phenotype increases each generation (where mean(opt) = k*t).
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

  ## standardize k
  k <- k*g_var

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
  if(is.null(n)){
    n <- length(phenotypes)
  }



  while(n[t] >= 2 & t <= gens){
    # stochastic survival
    t_kt <- sum(k*t + rnorm(1, 0, sqrt(var.theta)))
    surv <- rbinom(length(phenotypes), 1, BL_survival(phenotypes, t_kt, omega))

    nsurv <- sum(surv)
    if(nsurv > K){surv[sample(which(surv == 1), nsurv - K)] <- 0; nsurv <- K} # if we are above K, kill some more at random.

    n <- c(n, sum(surv))
    opt <- c(opt, t_kt)
    new_mean_phenos <- mean(phenotypes[which(surv == 1)])
    mean_phenos <- c(mean_phenos, new_mean_phenos)

    if(sum(surv) < 2){
      Rs <- c(Rs, NA)
      Ss <- c(Ss, NA)
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

    Rs <- c(Rs, R)
    Ss <- c(Ss, S)
    t <- t + 1
  }

  out <- data.frame(t = 1:length(n), n = n, mean_pheno = mean_phenos, opt_pheno = opt,
                    response_to_selection = Rs, selection_differential = Ss)
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
