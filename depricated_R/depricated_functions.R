



# Function to calculate the estimated time untill a population begins to crash (growth rate less than one) based on Burger and Lynch 1995.
#    g_var: addative genetic variance
#    e_var: environmental variance
#    omega: width of the fitness function, usually given as omega^2
#    k: rate of environmental change in phenotypic standard deviations
#    B: mean number of offspring per individual
#    Ne: effective population size
#    theta_var: environmental stochasticity
B_L_t1_func <- function(g_var, e_var, omega, k, B, Ne, theta_var){
  # calc Vs
  Vs <- (omega^2) + e_var

  # calc Vlam
  # simplified: Vlam = (Vs*(1+2*Ne))/2*Ne + (((1+2*Vs)*(g_var+theta_var))/2*Vs)
  V_gt <- (Vs/(2*Ne)) + (g_var*theta_var)/(2*Vs)
  Vlam <- Vs + g_var + V_gt + theta_var

  #calc kc
  Bo <- B*omega/sqrt(Vlam)
  if(Bo < 1){
    return(list(t1 = NA, kc = NA, Vs = Vs, Vlam = Vlam, Bo = Bo))
  }
  kc <- (g_var/(g_var + Vs))*sqrt(2*Vs*log(Bo))

  if(k<kc){
    t1 <- Inf
  }
  else{
    t1 <- -((g_var + Vs)/g_var)*log(1-(kc/k))
  }

  #calc t1
  return(list(t1 = t1, kc = kc, Vs = Vs, Vlam = Vlam, Bo = Bo))
}





# old gs function
gs <- function(x,
               gens,
               growth.function,
               survival.function,
               selection.shift.function,
               rec.dist,
               var.theta = 0,
               pred.method = "effects",
               plot_during_progress = FALSE,
               facet = "group", chr.length = 10000000,
               fgen.pheno = FALSE,
               intercept_adjust = FALSE,
               print.all.freqs = FALSE,
               adjust_phenotypes = FALSE,
               do.sexes = TRUE,
               init = F,
               verbose = T){
  if(verbose){cat("Initializing...\n")}
  #unpack x:
  if(pred.method == "effects"){ #unpack estimated effect sizes if provided.
    effect.sizes <- x$e.eff[,2]
  }
  if(fgen.pheno){ #unpack phenotypes if requested
    fgen.pheno <- x$phenotypes$p
  }
  h <- x$h
  meta <- x$meta
  if(pred.method != "real"){
    pred.mod <- x$output.model$mod
    pred.dat <- x$output.model$data
    model <- x$prediction.program
  }
  else{
    model <- "real"
    pred.method <- "effects" #since everything else works the same, just need to change inputs.
    effect.sizes <-  meta$effect
  }
  if(pred.method == "effects"){
    pred.mod <- NULL
    pred.dat <- NULL
  }
  if(pred.method == "model"){
    effect.sizes <- NULL
  }
  x <- x$x


  #=================checks========
  if(!pred.method %in% c("model", "effects")){
    stop("pred.method must be provided. Options:\n\tmodel: predict phenotypes directly from the model provided.\n\teffects: predict phenotypes from estimated effect sizes.\n")
  }

  if(!data.table::is.data.table(x)){
    x <- data.table::as.data.table(x)
  }

  if(pred.method == "effects"){
    if(nrow(x) != length(effect.sizes) | nrow(x) != nrow(meta)){
      stop("Provided x, effect sizes, and meta must all be of equal length!")
    }
  }
  else{
    if(nrow(x) != nrow(meta)){
      stop("Provided x and meta must be of equal length!")
    }
  }

  if(pred.method == "model"){
    if(!model %in% c("JWAS", "BGLR", "ranger")){
      stop("To predict from the model, a JWAS, BGLR, or ranger model must be provided.\n")
    }
  }
  else{
    if(model == "ranger"){
      stop("RF does not estimate effect sizes, so prediction must be done using the ranger model.\n")
    }
  }

  # before doing anything else, go ahead and remove any loci from those provided with no effect! Faster this way.
  # don't do this if initializing the population!
  if(pred.method == "effects" & !init){
    if(any(effect.sizes == 0)){
      n.eff <- which(effect.sizes == 0)
      x <- x[-n.eff,]
      meta <- meta[-n.eff,]
      effect.sizes <- effect.sizes[-n.eff]
    }
  }

  #=================get starting phenotypic values and BVs=========
  # get starting phenotypes and addative genetic values
  ## If first gen phenos aren't provided (should be uncommon)
  if(length(fgen.pheno) != ncol(x)/2){
    if(pred.method == "effects"){
      if(verbose){cat("Generating representative starting phenotypes from effect sizes.")}
      pheno <- get.pheno.vals(x, effect.sizes, h)
      a <- pheno$a # BVs
      pheno <- pheno$p # phenotypic values
    }
    else{
      if(verbose){cat("Generating representative starting phenotypes from model.")}
      a <- pred.BV.from.model(pred.mod, x, pred.method, model)
      pheno <- a +  e.dist.func(a, var(a), h) #add environmental effects
      #working here
    }
  }

  # otherwise use those, but still need to estimate BVs
  else{
    if(verbose){cat("Using provided phenotypic values.")}
    pheno <- fgen.pheno #provded phenotypic values.

    a <- pred.BV.from.model(pred.model = pred.mod, g = x, pred.method = pred.method,
                            model.source = model, h = h, h.av = "fgen", effect.sizes = effect.sizes)$a

  }

  #================set up BV variation adjustment to correct for drop in variance from GP methods============
  if(adjust_phenotypes){
    reorg_gcs <- rand.mating(x, ncol(x)/2, meta, rec.dist, chr.length, do.sexes, facet) # reorganize chrs once, since this causes one heck of a drop in var(a) in some GP results
    reorg_gcs <- rand.mating(reorg_gcs, ncol(x)/2, meta, rec.dist, chr.length, do.sexes, facet)
    re_p <- pred.BV.from.model(pred.model = pred.mod, g = reorg_gcs, pred.method = pred.method,
                               model.source = model, h = h, h.av = "fgen", effect.sizes = effect.sizes) # re-predict BVs
    re_a <- re_p$a
    re_p <- re_p$p
    adj.a.var <- var(re_a) #variance next gen


    # now need to adjust a and pheno to fit the variance a gen later
    # multiply future phenotypes by the square root of these values, then adjust the mean back to the correct mean.



    # old version which adjusts back to the starting phenotypic var every generation.
    # reorg_gcs <- rand.mating(x, ncol(x)/2, meta, rec.dist, chr.length, do.sexes, facet) # reorganize chrs once, since this causes one heck of a drop in var(a) in some GP results
    # re_p <- pred.BV.from.model(pred.model = pred.mod, g = reorg_gcs, pred.method = pred.method,
    #                                 model.source = model, h = h, h.av = "fgen", effect.sizes = effect.sizes) # re-predict BVs
    # re_a <- re_p$a
    # re_p <- re_p$p
    # ad.factor <- var(pheno)/(var(re_a)/h) # here's our adjustment factor
    # rm(re_a, re_p, reorg_gcs)
  }

  #if requested, get the amount to adjust phenotypes by in future gens.
  if(intercept_adjust){
    i.adj <- mean(pheno)
  }

  #================print out initial conditions, intiallize final steps, and run===========
  #starting optimal phenotype, which is the starting mean addative genetic value.
  opt <- mean(a) #optimum phenotype

  if(verbose){
    cat("\n\n===============done===============\n\nStarting parms:\n\tstarting optimum phenotype:", opt,
        "\n\tmean phenotypic value:", mean(pheno), "\n\taddative genetic variance:", var(a), "\n\tphenotypic variance:", var(pheno), "\n\th:", h, "\n")
  }

  #make output matrix and get initial conditions
  out <- matrix(NA, nrow = gens + 1, ncol = 8)
  colnames(out) <- c("N", "mu_pheno", "mu_a", "opt", "diff", "var_a", "stochastic_opt", "gen")
  N <- ncol(x)/2 #initial pop size
  h.av <- var(a) #get the historic addative genetic variance.
  h.pv <- var(pheno) #historic phenotypic variance.

  out[1,] <- c(N, mean(pheno), mean(a), opt, 0, h.av, opt, 1) #add this and the mean initial additive genetic variance
  if(plot_during_progress){
    library(ggplot2)
    pdat <- reshape2::melt(out)
    colnames(pdat) <- c("Generation", "var", "val")
    ranges <- data.frame(var = c("N", "mu_pheno", "mu_a", "opt", "diff"),
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
                survival.function(pheno, t.opt, hist.var = h.pv)) # calling the function in this way ensures that individuals with phenotypes at the optimum have a survival probability of whatever is set in the function.
    #if the population has died out, stop.
    if(sum(s) <= 1){
      break
    }

    #what is the pop size after growth?
    out[i,1] <- round(growth.function(sum(s)))

    #make a new x with the survivors
    x <- x[, .SD, .SDcols = which(rep(s, each = 2) == 1)] #get the gene copies of survivors

    # # check phenotypic variance...
    # temp <- get.pheno.vals(x, effect.sizes, h, hist.a.var = h.av)
    # ptemp <- data.frame(val = c(a, temp$a), class = c(rep("T0", length(a)), rep("T1", length(temp$a))))
    # temp <- tem$p
    # if(intercept_adjust){
    #   temp <- temp + i.adj
    # }
    # # adjust variance
    # if(adjust_phenotypes != FALSE){
    #   s.p.mean <- mean(temp)
    #   temp <- temp*sqrt(ad.factor)
    #   temp <- temp - (mean(temp) - s.p.mean)
    # }
    # print(var(temp))

    #=============do random mating, adjust selection, get new phenotype scores, get ready for next gen====
    y <- rand.mating(x, out[i,1], meta, rec.dist, chr.length, do.sexes, facet)
    # check that the pop didn't die due to every individual being the same sex (rand.mating returns NULL in this case.)
    if(is.null(y)){
      break
    }
    else{
      x <- y
      rm(y)
    }

    #get phenotypic/genetic values
    pa <- pred.BV.from.model(pred.model = pred.mod,
                             g = x,
                             pred.method = pred.method,
                             model.source = model,
                             h = h,
                             h.av = h.av,
                             effect.sizes = effect.sizes)
    a <- pa$a
    pheno <- pa$p

    #if requested, adjust the phenotypic values.
    # adjust intercept
    if(intercept_adjust){
      pheno <- pheno + i.adj
    }
    # adjust variance
    if(adjust_phenotypes != FALSE){
      s.p.mean <- mean(pheno)
      pheno <- pheno*sqrt(ad.factor)
      pheno <- pheno - (mean(pheno) - s.p.mean)
    }

    #adjust selection optima
    opt <- selection.shift.function(opt, iv = sqrt(h.av))

    #save
    out[i,2] <- mean(pheno)
    out[i,3] <- mean(a)
    out[i,4] <- opt
    out[i,5] <- opt - mean(a)
    out[i,6] <- var(a)
    out[i,7] <- t.opt
    if(verbose){
      cat("gen:", i-1,
          "\tf_opt:", round(out[i-1,4],3),
          "\ts_opt", round(out[i-1,7],3),
          "\tmean(pheno):", round(out[i,2],3),
          "\tmean(a):", round(out[i,3],3),
          "\tvar(a):", round(var(a),3),
          "\tNs:", sum(s),
          "\tN(t+1):", out[i,1],"\n")
    }
    if(plot_during_progress){
      pdat <- reshape2::melt(out)
      colnames(pdat) <- c("Generation", "var", "val")
      ranges <- data.frame(var = c("N", "mu_pheno", "mu_a", "opt", "diff"),
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
      a.fqs[,i] <- rowSums(x)/ncol(x)
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

  return(list(run_vars = out, x = x, phenos = pheno, BVs = a))
}

# from http://www2.univet.hu/users/jreiczig/locScaleTests/
lepage.stat=function(x1,x2){
  browser()
  enne1=as.numeric(length(x1))
  enne2=as.numeric(length(x2))
  enne=enne1+enne2
  e.w=enne1*(enne+1)/2
  v.w=enne1*enne2*(enne+1)/12
  e.a=enne1*(enne+2)/4
  v.a=enne1*enne2*(enne+2)*(enne-2)/48/(enne-1)
  w.o=as.numeric(wilcox.test(x1,x2,exact=FALSE)[1])+enne1*(enne1+1)/2
  a.o=as.numeric(ansari.test(x1,x2,exact=FALSE,alternative="two.sided")[1])
  wp.o=(w.o-e.w)^2/v.w
  ap.o=(a.o-e.a)^2/v.a
  return(wp.o+ap.o)
}
cucconi.stat=function(x1,x2){
  cuc=function(x1,x2){
    vett=c(x1,x2)
    enne1=as.numeric(length(x1))
    enne2=as.numeric(length(x2))
    enne=as.numeric(length(vett))
    ranghi=rank(vett)
    erre2=ranghi[(enne1+1):enne]
    media=enne2*(enne+1)*(2*enne+1)
    scarto=(enne1*enne2*(enne+1)*(2*enne+1)*(8*enne+11)/5)^0.5
    u=(6*sum(erre2^2)-media)/scarto
    v=(6*sum((enne+1-1*erre2)^2)-media)/scarto
    ro=2*(enne^2-4)/(2*enne+1)/(8*enne+11)-1
    cuc=(u^2+v^2-2*u*v*ro)/2/(1-ro^2)
  }
  return(.5*(cuc(x1,x2)+cuc(x2,x1)))
}



compare_peaks <- function(o, p){
  npdiff <- abs(nrow(o) - nrow(p))

  diffs <- rep(NA, 31)
  if(all(nrow(o) > 0 & nrow(p) > 0)){
    diffs <- calc_dist_stats(o$val, p$val)
  }

  return(c(npeak_diff = npdiff, diffs))
}
