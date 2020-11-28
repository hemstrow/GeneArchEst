#=======function to take in input data, use genomic prediction to estimate effect sizes based on phenotypes==============
#arugments:
#    x: Input data, matrix, df, or data.table. Converted internally to data.table. Columns are gene copies, rows are SNPs, formatted as 0 and 1.
#    effect.sizes: Vector of marker effect sizes. Will be optional eventually for prediction from phenotypes only.
#    ind.effects: Individual addative genetic values/BVs for individuals. Will eventually be optional for prediction from phenotypes only.
#    method: What method should we use to generate predictions? BGLR, JWAS, or ranger.
#    chain.length: How long should the MCMC chain in JWAS/BGLR be?
#    burnin: How many MCMC iterations should we discard at the start of the chain for JWAS/BGLR? Must be less than chain_length!
#    thin: How should the MCMC iterations be thinned? For BGLR only.
#    model: What model should we use? BGLR: FIXED, BL, BayesA-C, BRR. JWAS: G-BLUP, BayesA-C, RR-BLUP. ranger: RJ.
#    make.ig: should new input files for JWAS be created? If they don't exist already, set to TRUE.
#    sub.ig:
#    pass.resid: NULL or numeric >= 0. A numeric value tells the function to pass the
#                estimated residual variance in the model on to JWAS. A numeric value of 0 passes the exact variance,
#                a numeric value other than zero will fudge the variance number by up to the proportion given (1 fudges up to 100%).
#    pass.var: NULL or numeric >= 0. Like pass.resid, but for the true genetic variance.
#    standardize: Boolean. Should the addative genetic values be centered and scaled between -1 and 1 prior to entry into JWAS? Phenotypic values still won't be centered!
#    center: Boolean. Should the phenotypes be centered only?
#' @export
pred <- function(x, meta = NULL, effect.sizes = NULL, phenotypes = NULL,
                 prediction.program = "JWAS",
                 chain_length = 100000,
                 burnin = 5000,
                 thin = 100,
                 prediction.model = NULL,
                 make.ig = TRUE, sub.ig = FALSE, maf.filt = 0.05,
                 julia.path = "julia", runID = "r1", qtl_only = FALSE,
                 pass.resid = FALSE, pass.var = FALSE,
                 ntree = 50000,
                 mtry = 1,
                 h = NULL,
                 standardize = FALSE,
                 center = FALSE,
                 save.meta = TRUE, par = NULL, pi = NULL, pass_G = NULL, phased = F,
                 verbose = T){
  #============sanity checks================================
  # check that all of the required arguments are provided for the prediction.model we are running
  if(prediction.program %in% c("JWAS", "BGLR", "PLINK", "TASSEL", "ranger", "GMMAT")){

    # JWAS checks
    if(prediction.program == "JWAS"){
      # chain length and burnin
      if(!is.numeric(c(chain_length, burnin))){
        stop("Chain_length and burnin must be numeric values.")
      }
      if(chain_length <= burnin){
        stop("Chain_length must be larger than burnin in order to estimate effect sizes.")
      }

      # path to julia
      if(!is.character(julia.path)){
        stop("Invalid path to julia executable.")
      }
      if(!file.exists(julia.path)){
        stop("Invalid path to julia executable.")
      }

      # JWAS prediction.model-for now, only RR-BLUP.
      if(!prediction.model %in% c("RR-BLUP", "BayesB")){
        stop("Invalid JWAS prediction.model.")
      }

      # check that there is prior info to pass if pass resid is ture
      if((pass.resid != FALSE | pass.var != FALSE) & !is.null(effect.sizes)){
        stop("Marker effect sizes must be defined in order to pass prior residual and variance info to JWAS.")
      }

    }

    # BGLR checks:
    else if(prediction.program == "BGLR"){
      # chain length and burnin
      if(!is.numeric(c(chain_length, burnin))){
        stop("Chain_length and burnin must be numeric values.")
      }

      # check for accepted prediction.model.
      if(!prediction.model %in% c("BRR", "FIXED", "BayesA", "BayesB", "BayesC", "BL")){
        stop(paste0("prediction.model ", prediction.model, " not recognized by BGLR. Options: BRR, FIXED, BayesA-C, BL.\n"))
      }
    }

    # random forest checks:
    else if(prediction.program == "ranger"){
      # prediction.model
      if(!(prediction.model %in% c("RJ"))){
        stop("For random forest, the prediction.model must be RJ (random jungle) for now.")
      }

      # ntree
      if(!is.numeric(ntree)){
        stop("ntree must be an integer!")
      }
      if(ntree != floor(ntree)){
        stop("ntree must be an integer!")
      }

      if(mtry > 1 | mtry < 0){
        stop("mtry must be between 1 and 0.\n")
      }
    }

  }
  else{
    stop("Invalid prediction.program provided. Options: JWAS, BGLR, PLINK, TASSEL, or ranger.")
  }


  # check that we are either provided with input marker effects or with input phenotypes
  if(all(c(is.null(phenotypes), is.null(effect.sizes)))){
    stop("Individual phenotypes or marker effect sizes must be provided!")
  }
  else if(all(c(is.null(phenotypes), is.null(effect.sizes)))){
    warning("Both phenotypes and marker effect sizes provided. Input phenotypes will be ignored!")
  }

  # check that qtl_only filtering isn't requested if effect sizes aren't provided!
  if(qtl_only & is.null(effect.sizes)){
    stop("qtl_only filtering can only be performed if marker effect sizes are provieded.")
  }

  if(!is.null(effect.sizes)){
    if(!is.numeric(h)){
      stop("h must be provided if effect.sizes are used.\n")
    }
  }

  if(is.null(meta[1]) & save.meta){
    save.meta <- F
    warning("Since no SNP metadata provided, no SNP metadata will be saved.\n")
  }

  cat("Preparing model inputs...\n")

  #============subfunctions=================================
  # do filters if requested
  filter_snps <- function(x, qtl_only, sub.ig, maf.filt, effect.sizes){
    rejects <- numeric(nrow(x)) # track rejected snps
    # qtl_only
    if(qtl_only){
      if(sub.ig != FALSE | maf.filt != FALSE){
        warning("No subsampling (sub.ig) or maf filtering (maf.filt) will occur when qtl_only = TRUE.\n")
      }
      s.markers <- which(effect.sizes != 0)
      x <- x[s.markers,]
      rejects[-s.markers] <- 1
    }

    # otherwise, other filters to check
    else{
      #filter low minor allele frequencies if requested (like many prediction studies will do!)
      if(maf.filt != FALSE){
        af <- matrixStats::rowSums2(x)/ncol(x)
        m.keep <- af >= maf.filt & af <= (1-maf.filt)
        x <- x[m.keep,]
        rejects[-which(m.keep)] <- 1
      }

      #subset markers
      if(sub.ig != FALSE){
        if(nrow(x) > sub.ig){
          s.markers <- sort(sample(nrow(x), sub.ig))
          x <- x[s.markers,]
          rejects[rejects != 1][-s.markers] <- 1
        }
        else{
          warning("Fewer markers than sub.ig. Running all markers.\n")
        }
      }
    }
    return(list(x = x, snp_ids = which(rejects == 0)))
  }


  #============prepare directories and phenotypes=========
  if(!is.null(phenotypes)){
    # if standardization is requested, set the var(pheno) to 1
    if(standardize){
      phenotypes <- phenotypes/sd(phenotypes)
    }
    # centering, set mean to 0
    if(center | standardize){
      phenotypes <- phenotypes - mean(phenotypes)
    }
    r.ind.effects <- list(p = phenotypes) # backup the ind effects.
  }

  # get phenotypes if effect sizes are provided.
  if(!is.null(effect.sizes)){
    phenotypes <- get.pheno.vals(x, effect.sizes, h = h, standardize = standardize)
    r.ind.effects <- phenotypes
    phenotypes <- phenotypes$p
  }

  #============format data for prediction/GWAS==============
  #filter:
  x <- filter_snps(x, qtl_only, sub.ig, maf.filt, effect.sizes)
  kept.snps <- x$snp_ids
  if(!is.null(meta)){
    meta <- meta[kept.snps,]
  }
  x <- x$x


  if(prediction.program == "JWAS"){
    odw <- getwd()
    if(!dir.exists(runID)){
      dir.create(runID)
    }
    setwd(runID)

    # make an individual effect file.
    ind.effects <- cbind(samp = as.character(1:500), phenotypes = phenotypes)
    write.table(ind.effects, "ie.txt", quote = F, col.names = T, row.names = F)

    # make an individual genotype file if it isn't already constructed.
    if(make.ig){
      # convert format
      ind.genos <- convert_2_to_1_column(x)
      ind.genos <- cbind(samp = 1:nrow(ind.genos), ind.genos) # add sample info
      colnames(ind.genos) <- c("samp", paste0("m", 1:(ncol(ind.genos)-1)))
      data.table::fwrite(ind.genos, "ig.txt", sep = " ", col.names = T)
    }
  }

  else if(prediction.program == "BGLR"){
    # convert
    ind.genos <- convert_2_to_1_column(x) # rows are individuals, columns are SNPs
    colnames(ind.genos) <- paste0("m", 1:ncol(ind.genos)) # marker names
    rownames(ind.genos) <- paste0("s", 1:nrow(ind.genos)) # ind IDS

    # prepare ETA
    ETA <- list(list(X = ind.genos, model = prediction.model, saveEffects = T))
  }

  else if(prediction.program == "ranger"){
    t.x <- convert_2_to_1_column(x) # add sample info
    colnames(t.x) <- paste0("m", 1:ncol(t.x))

    t.eff <- data.frame(phenotype = phenotypes, stringsAsFactors = F)
    t.eff <- cbind(t.eff, t.x)
    colnames(t.eff)[-1] <- paste0("m", 1:ncol(t.x))
  }

  else if(prediction.program == "GMMAT"){
    ind.genos <- convert_2_to_1_column(x)
    colnames(ind.genos) <- paste0("m", 1:ncol(ind.genos)) # marker names
    rownames(ind.genos) <- paste0("s", 1:nrow(ind.genos)) # ind IDS

    if(is.null(pass_G)){
      G <- make_G(ind.genos, maf.filt, phased, par)
    }
    else{
      G <- pass_G
      rm(pass_G)
    }
  }


  #=========run genomic prediction or GWAS and return results==========
  # for JWAS:
  if(prediction.program == "JWAS"){
    cat("Calling JWAS.\n")
    options(scipen = 999)
    julia.call <- paste0(julia.path, " ", getwd(), "/analysis.jl ", chain_length, " ", burnin)
    # add the residual and genetic variance if requested.
    if(pass.resid != FALSE){
      rv <- var(r.ind.effects$p - r.ind.effects$a)
      rv <- rv + rv*runif(1, 0, pass.resid) #fudge according to factor provided
    }
    else{
      rv <- 1
    }
    if(pass.var != FALSE){
      gv <- var(r.ind.effects$a)
      gv <- gv + gv*runif(1, 0, pass.var) #fudge according to factor provided
    }
    else{
      gv <- 1
    }
    julia.call <- paste0(julia.call, " ", rv, " ", gv, " ", prediction.model)
    if(!is.null(pi)){
      julia.call <- paste0(julia.call, " ", pi)
    }
    else{
      julia.call <- paste0(julia.call, " ", "false")
    }

    # save the julia script
    browser()
    writeLines(analysis.jl, "analysis.jl")
    system(julia.call)

    #=========grab output and modify it to give the estimated effect size per locus=============
    e.eff <- read.table("est_effects.txt", header = F, sep = "\t")
    h <- read.table("h.txt")

    #save metadata for the selected markers if requested.
    if(save.meta){
      write.table(meta, "est_meta.txt", sep = "\t", quote = F, col.names = T, row.names = F)
    }

    setwd(owd)

    return(list(x = x, e.eff = e.eff, phenotypes = r.ind.effects, meta = meta, h = as.numeric(h), kept.snps = kept.snps))
  }

  # for BGLR:
  else if(prediction.program == "BGLR"){
    cat("Calling BGLR.\n")
    BGLR_mod <- BGLR::BGLR(y = phenotypes, ETA = ETA, nIter = chain_length, burnIn = burnin, thin = thin, verbose = verbose)

    # grab h2 estimate
    B <- BGLR::readBinMat('ETA_1_b.bin')
    h2 <- rep(NA,nrow(B))
    varU <- h2
    varE <- h2
    for(i in 1:length(h2)){
      u <- ind.genos%*%B[i,]
      varU[i] <- var(u)
      varE[i] <- var(phenotypes-u)
      h2[i] <- varU[i]/(varU[i] + varE[i])
    }
    h2 <- mean(h2)

    # return the values
    e.eff <- data.frame(V1 = BGLR_mod$ETA[[1]]$colNames, V2 = BGLR_mod$ETA[[1]]$b, stringsAsFactors = F)
    write.table(e.eff, "est_effects.txt", sep = "\t", col.names = F, row.names = F, quote = F)
    write.table(h2, "h.txt", quote = F)
    if(save.meta){
      write.table(meta, "est_meta.txt", sep = "\t", quote = F, col.names = T, row.names = F)
    }

    return(list(x = x, e.eff = e.eff, phenotypes = r.ind.effects, meta = meta, h = h2, prediction.program = "BGLR",
                prediction.model = prediction.model, output.model = list(mod = BGLR_mod, data = ETA), kept.snps = kept.snps))
  }

  # for randomForest
  else if(prediction.program == "ranger"){
    cat("Running randomforest (ranger rj implementaiton).\n")

    # figure out the mtry to use
    mtry <- nrow(x)*mtry

    # run the randomForest/jungle
    if(ncol(t.eff) - 1 >= 10000){
      rj <- ranger::ranger(dependent.variable.name = "phenotype", data = t.eff, mtry = mtry,
                           num.trees = ntree, verbose = T, save.memory = T, num.threads = par)
    }
    else{
      rj <- ranger::ranger(dependent.variable.name = "phenotype", data = t.eff,
                           mtry = mtry, num.trees = ntree, verbose = T, num.threads = par)
    }

    return(list(x = x, phenotypes = r.ind.effects, meta = meta, prediction.program = "ranger",
                prediction.model = "RJ", output.model = list(model = rj), kept.snps = kept.snps))
  }

  # for GMMAT
  else if(prediction.program == "GMMAT"){

    if(length(unique(phenotypes)) > 2){
      family <- gaussian(link = "identity")
    }
    else if(length(unique(phenotypes)) == 2){
      family <- binomial(link = "logit")
    }

    # run null model
    mod <- GMMAT::glmmkin(fixed = "phenotypes ~ 1",
                          data = data.frame(phenotypes = phenotypes, sampleID = rownames(ind.genos)),
                          kins = G,
                          id = "sampleID",
                          family = family, method.optim = "Brent")

    # run the test
    ## prepare infile
    infile <- paste0(runID, "_asso_in.txt")
    outfile <- paste0(runID, "_asso_out_score.txt")
    asso_in <- as.data.table(t(ind.genos))
    asso_in[, "snp_id" := colnames(ind.genos)]
    setcolorder(asso_in, c("snp_id", colnames(asso_in)[1:nrow(ind.genos)]))
    suppressMessages(data.table::fwrite(asso_in, infile, sep = "\t", col.names = T, row.names = F))
    ## run
    score.out <- GMMAT::glmm.score(obj = mod,
                                   infile = infile,
                                   outfile = outfile,
                                   infile.nrow.skip = 1,
                                   infile.ncol.skip = 1)
    score.out <- read.table(outfile, header = T, stringsAsFactors = F)
    ## clean
    file.remove(c(outfile, infile))

    # return
    return(list(x = x, e.eff = score.out, phenotypes = r.ind.effects, meta = meta, prediction.program = "GMMAT",
                prediction.model = prediction.model, kept.snps = kept.snps))
  }
}

#' Run a GMMAT GWAS on genotype/phenotype data.
#'
#' @param x genotype data, default NULL. Genotype
#'   with individuals in columns (either 1 or 2 columns per individual for unphased or phased data). Can either
#'   be coercable to a matrix or be a \code{\link[bigstatsr]{FBM}}.
#' @param phenotypes numeric. Phenotype data, one numeric value per individual, in the same order as the genotype data.
#' @param maf numeric, default 0.05. Minor allele frequency filter to be applied to both G-matrix and GWAS. Ignored if x is NULL.
#' @param pass_G numeric matrix, default NULL. An \emph{n x n} genetic relatedness matrix. If NULL, calculated from x via
#'   \code{\link[AGHmatrix]{Gmatrix}} using the Yang et al 2010 method.
#' @param GMMAT_infile character, default NULL. If provided, the file path to an infile for \code{\link[GMMAT]{glmm.score}}. If NULL,
#'   will be written from x.
#' @param phased logical, default F. Is the data in x phased (in two columns per individual)? Ignored if x is NULL.
#' @param par numeric, default 1. Number of parallel cores to use for G matrix creation. Ignored if x is NULL.
#' @param center logical, default T. If TRUE, the phenotypes will be centered (set to mean 0) prior to GWAS calculation
#'
#' @export
pred_gwas_FBM <- function(x = NULL, phenotypes, maf = 0.05, pass_G = NULL, GMMAT_infile = NULL,
                          phased = F, par = 1, center = T){
  #============transpose==================
  if(!is.null(x)){
    if(phased == T){
      xt <- convert_2_to_1_column(x)
    }
    else{
      xt <- bigstatsr::big_transpose(x)
    }
  }

  #============G prep or import===========
  if(is.null(pass_G)){
    G <- make_G(xt, maf, phased = F, par)
  }
  else{
    G <- pass_G
    rm(pass_G)
  }
  #============center phenotypes==========
  if(center){
    phenotypes <- phenotypes - mean(phenotypes)
  }

  #============prep gmmat infile or import=========
  if(is.null(GMMAT_infile)){
    if(!is.null(maf)){
      if(maf > 0){
        Frequency <- bigstatsr::big_colstats(xt)$sum/(2*nrow(xt))
        Frequency <- cbind(1 - Frequency, Frequency)
        Frequency <- cbind(Frequency, matrixStats::rowMins(Frequency))

        to.rm <- which(Frequency[,3] <= maf)
        if(length(to.rm) > 0){
          xt <- bigstatsr::FBM(nrow = nrow(xt), ncol = ncol(xt), type = "integer",
                                      xt[,-to.rm])
        }
        rm(Frequency)
      }
    }

    if("FBM" %in% class(xt)){
      bigstatsr::big_write(xt, file = "asso_in.txt", sep = "\t")
    }
    else{
      data.table::fwrite(xt, "asso_in.txt" ,sep = "\t")
    }
    GMMAT_infile <- "asso_in.txt"
  }

  #=========================run association==========
  if(length(unique(phenotypes)) > 2){
    family <- gaussian(link = "identity")
  }
  else if(length(unique(phenotypes)) == 2){
    family <- binomial(link = "logit")
  }

  # run null model
  colnames(G) <- 1:length(phenotypes)
  rownames(G) <- 1:length(phenotypes)
  mod <- GMMAT::glmmkin(fixed = "phenotypes ~ 1",
                        data = data.frame(phenotypes = phenotypes, sampleID = 1:length(phenotypes)),
                        kins = G,
                        id = "sampleID",
                        family = family, method.optim = "Brent")

  # run the GWAS
  outfile <- paste0("GMMAT_score.out")
  score.out <- GMMAT::glmm.score(obj = mod,
                                 infile = GMMAT_infile,
                                 outfile = outfile,
                                 infile.nrow.skip = 0,
                                 infile.ncol.skip = 0, infile.header.print = NULL,
                                 infile.ncol.print = NULL)
  score.out <- read.table(outfile, header = T, stringsAsFactors = F)

  # return
  return(list(e.eff = score.out))
}




