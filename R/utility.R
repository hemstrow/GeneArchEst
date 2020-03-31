#' @export
process_ms <- function(x, chr.length){
  infile <- x #infile
  lines <- readLines(x)
  lines <- lines[-which(lines == "")] #remove empty entries
  lines <- lines[-c(1,2)] #remove header info
  nss <- grep("segsites", lines) #get the number of segsites per chr
  chrls <- gsub("segsites: ", "", lines[nss]) #parse this to get the lengths
  chrls <- as.numeric(chrls)
  lines <- lines[-nss] #remove the segsites lines
  pos <- lines[grep("positions:", lines)] #find the positions
  lines <- lines[-grep("positions:", lines)] #remove the position
  div <- grep("//", lines) #find the seperators
  gc <- div[2] - div[1] - 1 #find the number of gene copies per chr
  if(is.na(gc)){gc <- length(lines) - 1} #if there's only one chr
  dat <- lines[-div] #get the data only
  dat <- strsplit(dat, "") #split the lines by individual snp calls
  x <- matrix(NA, nrow = sum(chrls), ncol = gc) #prepare output
  meta <- matrix(NA, nrow = sum(chrls), 2)

  #process this into workable data
  pchrls <- c(0, chrls)
  pchrls <- cumsum(pchrls)
  for(i in 1:length(chrls)){
    cat("\n\tChr ", i)
    tg <- dat[(gc*(i-1) + 1):(gc*i)] #get only this data
    tg <- unlist(tg) #unlist
    tg <- matrix(as.numeric(tg), ncol = chrls[i], nrow = gc, byrow = T) #put into a matrix
    tg <- t(tg) #transpose. rows are now snps, columns are gene copies
    tpos <- unlist(strsplit(pos[i], " ")) #grap and process the positions
    tpos <- tpos[-1]
    meta[(pchrls[i] + 1):pchrls[i + 1],] <- cbind(paste0(rep("chr", length = nrow(tg)), i), tpos)
    x[(pchrls[i] + 1):pchrls[i + 1],] <- tg #add data to output
  }

  meta <- as.data.frame(meta, stringsAsFactors = F)
  meta[,2] <- as.numeric(meta[,2])
  meta[,2] <- meta[,2] * chr.length

  colnames(meta) <- c("group", "position")
  colnames(x) <- paste0("gc_", 1:ncol(x))

  return(list(x = x, meta = meta))
}



# wrapper for gs() with no shift in surival and assuming real effect sizes (no model). x needs to contain phenotypes and meta with effects
#' @export
init_pop <- function(x,
                     init_gens,
                     growth.function,
                     survival.function,
                     rec.dist,
                     var.theta = 0,
                     plot_during_progress = FALSE,
                     facet = "group", chr.length = 10000000,
                     fgen.pheno = F,
                     print.all.freqs = FALSE,
                     do.sexes = TRUE){
  #=======set the parms for calling gs======
  s.shift.null <- function(x, ...){
    return(x)
  }

  if(!"effect" %in% colnames(x$meta)){
    stop("Effects must be provided in as a column in x$meta.\n")
  }
  if(!"h" %in% names(x)){
    stop("heritability must be provided in x$h.\n")
  }
  if(!is.numeric(x$h) | length(x$h) != 1){
    stop("heritability must be a single numeric value.\n")
  }
  #call gs correctly to initialize with no shift in the survival function and assuming real effects.
  out <- gs(x = x,
            gens = init_gens,
            growth.function = growth.function,
            survival.function = survival.function,
            selection.shift.function = s.shift.null,
            rec.dist = rec.dist,
            var.theta = var.theta,
            pred.method = "real",
            plot_during_progress = plot_during_progress,
            facet = facet,
            chr.length = chr.length,
            fgen.pheno = fgen.pheno,
            intercept_adjust = F,
            print.all.freqs = print.all.freqs,
            adjust_phenotypes = F,
            do.sexes = do.sexes,
            init = T)
  return(out)
}



#=======function to do a single generation of random mating===========
#' @export
rand.mating <- function(x, N.next, meta, rec.dist, chr.length, do.sexes = TRUE, facet = "group"){
  if(!data.table::is.data.table(x)){
    x <- data.table::as.data.table(x)
  }
  #=========get parents and assign gcs for the next gen====
  #make a new x with individuals in next gen
  ##find parents
  if(do.sexes){ # if there are two sexes
    sex <- rbinom(ncol(x)/2, 1, 0.5) #what are the parent sexes?
    if(sum(sex) == length(sex) | sum(sex) == 0){ # if every individual is the same sex, the population dies.
      return(NULL)
    }
    mates <- matrix(0, nrow = N.next, ncol = 2) # initialize, p1 and p2 are columns
    mates[,1] <- which(sex == 1)[sample(sum(sex), nrow(mates), T)] #get parents of sex a
    mates[,2] <- which(sex == 0)[sample(length(sex) - sum(sex), nrow(mates), T)] #get parents of sex b
  }
  else{ # if there is only one sex
    mates <- matrix(sample(ncol(x)/2, N.next*2, T), ncol = 2) #p1 and p2 are columns
    selfings <- which(mates[,1] == mates[,2]) #any selfing?
    while(length(selfings) > 0){ #correct selfing
      mates[selfings,] <- sample(ncol(x)/2, length(selfings)*2, T) #get new parents
      selfings <- which(mates[,1] == mates[,2]) #any selfing remaining?
    }
  }

  # table with the showing the distribution of the number of offspring for each adult:
  # table(c(table(mates), rep(0, (ncol(x)/2) - length(table(mates)))))
  x.next <- data.table::as.data.table(matrix(0, nrow(x), nrow(mates)*2)) #initialize x matrix for next gen

  #=========figure out which copy from each parent goes to offspring=====
  #randomly choose gene copies to push to individuals in the next gen.
  # for each individual, do they get copy 1 or copy 2 from the parent?
  uf <- unique(meta[,facet])
  chr.source.p1.i <- matrix(rbinom(nrow(mates)*length(uf), 1, .5), ncol = length(uf), byrow = T) + 1
  chr.source.p2.i <- matrix(rbinom(nrow(mates)*length(uf), 1, .5), ncol = length(uf), byrow = T) + 1

  # add the correct copies
  ##which column does the data come from?
  chr.source.p1 <- ifelse(chr.source.p1.i == 1, mates[,1] * 2 - 1,
                          mates[,1] * 2)

  chr.source.p2 <- ifelse(chr.source.p2.i == 1, mates[,2] * 2 - 1,
                          mates[,2] * 2)

  #the other copies?
  chr.nsource.p1 <- ifelse(chr.source.p1.i == 1, mates[,1] * 2,
                           mates[,1] * 2 - 1)
  chr.nsource.p2 <- ifelse(chr.source.p2.i == 1, mates[,2] * 2,
                           mates[,2] * 2 - 1)

  rm(chr.source.p1.i, chr.source.p2.i)

  #these now say which column in x to take the data from for each chr for each individual for bases that didn't recombine.

  #=========recombination and chromosome assignment======
  num.rec <- rec.dist(nrow(mates)*2*length(uf))
  rs <- sum(num.rec) #total number
  rec.pos <- sample(chr.length, rs, T) #get positions for each of these events. Assumes equal chr length, otherwise would need to put this into the per. chr loop. Could do later if necissary, but it'd be a bit slower.

  prog <- 0 #progress through num.rec tracker.
  pos.prog <- 0 #progress through recombination events tracker.

  #fill for each facet
  for(j in 1:length(unique(uf))){
    # cat(j,"\n")
    #overall approach: for each parent:
    # line up copy 1 and copy 2 chromosomes in a matrix, each column is a seperate chr. Copy one is the copy that is getting passed! Copy two is the one that isn't.
    # for each recombination event, on each chr, flip which chr we are taking from. For portions where we are taking from chr2, paste into chr 1, which will be the output.
    this.chr.pos <- meta$position[meta[,facet] == uf[j]]

    for(k in 1:2){

      trec <- c(prog + 1, prog + nrow(mates)) #which recombination events are we working with here?
      prog <- prog + nrow(mates)

      #get the number of recombination events per chr in this set
      tnr <- num.rec[trec[1]:trec[2]]

      #initialize matrix
      c1.mat <- data.table::as.data.table(matrix(0, length(this.chr.pos), ncol = nrow(mates)))
      c2.mat <- data.table::as.data.table(matrix(0, length(this.chr.pos), ncol = nrow(mates)))

      #paste in the values from x. c1.mat contains the "passed" chr, c2.mat contains the "unpassed" chr.
      if(k == 1){
        data.table::set(c1.mat, 1:nrow(c1.mat), as.integer(1:ncol(c1.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.source.p1[,j]])
        data.table::set(c2.mat, 1:nrow(c2.mat), as.integer(1:ncol(c2.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.nsource.p1[,j]])
      }
      else{
        data.table::set(c1.mat, 1:nrow(c1.mat), as.integer(1:ncol(c1.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.source.p2[,j]])
        data.table::set(c2.mat, 1:nrow(c2.mat), as.integer(1:ncol(c2.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.nsource.p2[,j]])
      }

      #only the recombining entries.
      wnz <- which(tnr != 0)
      if(length(wnz) == 0){
        #no recombination, mostly if pop size is VERY small...
        if(k == 1){
          data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(1, ncol(x.next), by = 2)), c1.mat)
        }
        else{
          data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(2, ncol(x.next), by = 2)), c1.mat)
        }
        next()
      }
      tnr_nz <- tnr[wnz]



      #get the positions of the recombination events
      trpos <- rec.pos[(pos.prog + 1):(pos.prog + sum(tnr_nz))]
      pos.prog <- pos.prog + sum(tnr_nz)


      # Now need to make and assign the actual output vector. I can't think of a good way to do this in an actually vectorized way, mostly because I have to look up positions to reference against for each chr.
      sort.prog.trpos <- 0 #how many positions have we searched through?
      for(m in 1:length(tnr_nz)){
        # cat("\t\t", m, "\n")
        sort.pos <-
          c(sort(trpos[(sort.prog.trpos + 1):(sort.prog.trpos + tnr_nz[m])])) #sort the correct recomb positions and add the ending positions.
        sort.prog.trpos <- sort.prog.trpos + tnr_nz[m] #update.

        #now figure out which chr each position will be drawn from.
        chr.ident <- numeric(length(this.chr.pos))
        for(q in 1:length(sort.pos)){
          chr.ident <- chr.ident + (this.chr.pos <= sort.pos[q])
        }
        chr.ident <- chr.ident %% 2 #if the number of crossing over events was even, 0s mean copy 1 and 1s mean copy 2. Otherwise reversed.

        #assign.
        if(length(sort.pos) %% 2 == 0){
          data.table::set(c1.mat, which(chr.ident == 1), j = wnz[m], value = c2.mat[which(chr.ident == 1), wnz[m], with = FALSE])
          #assign entries where the chr flips in c1 mat to the respecitve entries in c2 mat.
        }
        else{
          data.table::set(c1.mat, which(chr.ident == 0), j = wnz[m], value = c2.mat[which(chr.ident == 0), wnz[m], with = FALSE])
        }

      }

      #================assign chromosomes to x.next.
      if(k == 1){
        data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(1, ncol(x.next), by = 2)), c1.mat)
      }
      else{
        data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(2, ncol(x.next), by = 2)), c1.mat)
      }
    }
  }

  return(x.next)
}




# function to extract BV predictions from a model
# pred.model: The GWAS/GP/ECT provided model
# g: the genotype matrix, where columns are gene copies
# pred.method: Should the BVs be predicted directly off the model ("model") or off of loci effects ("effects")?
# model.source: Which program made the model? Currently supports BGLR, JWAS, or RJ (ranger).
# h: heritability estimate
# h.av: historic genetic varaince, for prediction from effect sizes.
# effect.sizes: marker effect sizes, for prediction from effect sizes.
#' @export
pred.BV.from.model <- function(pred.model, g, pred.method = NULL, model.source = NULL, h = NULL, h.av = "fgen", effect.sizes = NULL){
  if(pred.method == "effects"){
    pheno <- get.pheno.vals(g, effect.sizes, h, hist.a.var = h.av)
    a <- pheno$a #addative genetic values
    pheno <- pheno$p #phenotypic values
    return(list(a = a, p = pheno))
  }

  else{
    if(!is.matrix(g)){g <- as.matrix(g)}
    if(model.source == "BGLR"){
      g <- convert_2_to_1_column(g)
      a <- as.vector(g%*%pred.model$ETA[[1]]$b)
    }
    else if(model.source == "RJ"){
      g <- convert_2_to_1_column(g)
      colnames(g) <- pred.model$forest$independent.variable.names
      a <- predict(pred.model, as.data.frame(g))
    }
    else if(model.souce == "JWAS"){
      #in progress
    }
    else{
      #in progress
    }

    if(h.av == "fgen"){
      h.av <- var(a)
    }
    pheno <- a + e.dist.func(a, h.av, h)

    return(list(a = a, p = pheno))
  }
}





#function to do column sums faster
src <- '
  Rcpp::NumericMatrix dataR(data);
  Rcpp::NumericVector weightsR(weights);
  int ncol = dataR.ncol();
  Rcpp::NumericVector sumR(ncol);
  for (int col = 0; col<ncol; col++){
  sumR[col] = Rcpp::sum(dataR( _, col)*weightsR);
  }
  return Rcpp::wrap(sumR);'

#' @export
weighted.colSums <- inline::cxxfunction(
  signature(data="numeric", weights="numeric"), src, plugin="Rcpp")




#function to generate random environmental effects for a given set of BVs, addative genetic variance (either historic or current from BVs are typical), and heritability.
# note that standardization isn't perfect, but will result in data with a mean very close to 0 and var close to 1. The randomness of assigning random environmental effects everywhere will make it imperfect
# downstream standardization of the phenotypes will fix this if using estimated effect sizes!
e.dist.func <- function(A1, hist.a.var, h, standardize = F){
  esd <- sqrt((hist.a.var/h)-hist.a.var) # re-arrangement of var(pheno) = var(G) + var(E) and h2 = var(G)/var(pheno)
  env.vals <- rnorm(length(A1), 0, esd)

  #if it standardization is requested, do so
  if(standardize){
    env.vals <- env.vals/sqrt(var(env.vals)/(1-h)) # set variance to 1 - h
    env.vals <- env.vals - mean(env.vals) # set mean to 0. Should be close, but not perfect because of the random draws.
  }
  return(env.vals)
}

#get phenotypic values given genotypes, effect sizes, and heritabilities. If hist.a.var is true, uses the amount of genomic variability this gen and h to figure out how big of an env effect to add. Otherwise uses the provided value (probably that in the first generation).
#' @export
get.pheno.vals <- function(x, effect.sizes, h, hist.a.var = "fgen", standardize = FALSE){
  #get effect of each individual:
  a <- weighted.colSums(as.matrix(x), effect.sizes) # faster than t(x)%*%effect.sizes!
  a.ind <- a[seq(1, length(a), by = 2)] + a[seq(2, length(a), by = 2)] #add across both gene copies.

  #standardize the genetic variance if requested.
  if(standardize){
    a.ind <- a.ind/sqrt(var(a.ind)/h) # set the variance to h.
    a.ind <- a.ind - mean(a.ind) # set the mean to 0
  }

  #add environmental variance
  if(hist.a.var == "fgen"){
    pheno <- a.ind + e.dist.func(a.ind, var(a.ind), h, standardize)
    return(list(p = pheno, a = a.ind))
  }
  else{
    pheno <- a.ind + e.dist.func(a.ind, hist.a.var, h, standardize)
    return(list(p = pheno, a = a.ind))
  }
}

#converts 2 column to 1 column genotypes and transposes
convert_2_to_1_column <- function(x){
  if(!is.matrix(x)){x <- as.matrix(x)}
  ind.genos <- x[,seq(1,ncol(x), by = 2)] + x[,seq(2,ncol(x), by = 2)]
  ind.genos <- matrix(ind.genos, nrow = ncol(x)/2, byrow = T) # rematrix and transpose!
  return(ind.genos)
}
