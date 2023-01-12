#' @export
process_ms <- function(x, chr.length, fix_overlaps = T){
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

  if(fix_overlaps){
    meta <- offset_overlapping_positions(meta)
  }

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
rand.mating <- function(x, N.next, meta, rec.dist, chr.length, do.sexes = TRUE){

  #-=========prep========

  if(length(unique(meta[,1])) != length(chr.length)){
    stop("The number of unique chromosomes is not equal to the number of chromsome lengths provided.\n")
  }
  if(!data.table::is.data.table(x)){
    x <- data.table::as.data.table(x)
  }
  facet <- colnames(meta)[1]
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
  x.next <- data.table::as.data.table(matrix(0, 1, nrow(mates)*2))[rep(1, nrow(x))] #initialize x matrix for next gen

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
  n.rec.per.chr <- tapply(num.rec, rep(chr.length, each = nrow(mates)*2), sum) # figure out how many recombination events per chromosome
  rec.pos <- runif(rs, min = 0, max = rep(chr.length, n.rec.per.chr)) # get positions for each of these events.

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

      #================assign chromosomes to x.next.==========
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
weighted.colSums <- function(data, weights){
  if("FBM" %in% class(data)){
    return(bigstatsr::big_cprodVec(data, weights))
  }
  else{
    return(crossprod(as.matrix(data), weights))
  }
}
# weighted.colSums <- inline::cxxfunction(
#   signature(data="numeric", weights="numeric"), src, plugin="Rcpp")




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
get.pheno.vals <- function(x, effect.sizes, h, hist.a.var = "fgen", standardize = FALSE, phased = FALSE, fitnesses = FALSE){

  # additive
  if(is.null(ncol(effect.sizes))){
    if(fitnesses){
      stop("Cannot supply additive effects with one value per locus if providing fitnesses.\n")
    }

    # remove zeros
    zeros <- which(effect.sizes == 0)
    if(length(zeros) != 0 & length(zeros) != length(effect.sizes)){
      effect.sizes <- effect.sizes[-zeros]
      x <- x[-zeros,, drop = FALSE]
    }

    #get effect of each individual:
    a.ind <- weighted.colSums(x, effect.sizes) # faster than t(x)%*%effect.sizes!
    if(phased){
      a.ind <- a.ind[seq(1, length(a.ind), by = 2)] + a.ind[seq(2, length(a.ind), by = 2)] #add across both gene copies.
    }
    if(is.matrix(a.ind)){
      a.ind <- as.numeric(a.ind)
    }
  }

  # non-additive (three effects)
  else{
    effect.sizes <- as.matrix(effect.sizes)

    # remove zeros
    zeros <- which(rowSums(effect.sizes) == ifelse(fitnesses, 3, 0))
    if(length(zeros) != 0 & length(zeros) != nrow(effect.sizes)){
      effect.sizes <- effect.sizes[-zeros,, drop = FALSE]
      x <- x[-zeros,, drop = FALSE]
    }


    # fix if fitnesses
    if(fitnesses){
      effect.sizes[effect.sizes > 1] <- 1
      effect.sizes[effect.sizes < 0] <- 0
    }

    if(phased){
      x <- convert_2_to_1_column(x)
      x <- t(x)
    }

    if(!fitnesses){
      a.ind <- (x == 0)*effect.sizes[,1] + (x == 1)*effect.sizes[,2] + (x == 2)*effect.sizes[,3]
      a.ind <- matrixStats::colSums2(a.ind)
    }
    else{
      a.ind <- (x == 0)*effect.sizes[,1] + (x == 1)*effect.sizes[,2] + (x == 2)*effect.sizes[,3]
      a.ind <- matrixStats::colProds(a.ind)
    }
  }



  # make sure h isn't 0, less than 0, or greater than 1
  if(h <= 0){
    h <- 0.00000001
  }
  else if(h > 1){
    h <- 1
  }

  #standardize the genetic variance if requested.
  if(standardize){
    a.ind <- a.ind/sqrt(var(a.ind)/h) # set the variance to h.
    a.ind <- a.ind - mean(a.ind) # set the mean to 0
  }

  #add environmental variance
  if(hist.a.var == "fgen"){
    pheno <- a.ind + e.dist.func(a.ind, var(a.ind), h, standardize)
  }
  else{
    pheno <- a.ind + e.dist.func(a.ind, hist.a.var, h, standardize)
  }

  # set max to 1 if fitnesses
  if(fitnesses){
    pheno[pheno > 1] <- 1
  }

  return(list(p = pheno, a = a.ind))
}

#converts 2 column to 1 column genotypes and transposes
#' @export
convert_2_to_1_column <- function(x){
  if("FBM" %in% class(x)){
    x1 <- bigstatsr::FBM(nrow(x), ncol(x)/2, init =
                           x[,seq(1, ncol(x), by = 2)])
    x2 <- bigstatsr::FBM(nrow(x), ncol(x)/2, init =
                           x[,seq(2, ncol(x), by = 2)])
    x <- bigstatsr::big_apply(x1, a.FUN = function(x, inds) x1[,inds] + x2[,inds])
    return(bigstatsr::big_transpose(x))
  }
  else{
    if(!is.matrix(x)){x <- as.matrix(x)}
    ind.genos <- x[,seq(1,ncol(x), by = 2)] + x[,seq(2,ncol(x), by = 2)]
    ind.genos <- matrix(ind.genos, nrow = ncol(x)/2, byrow = T) # rematrix and transpose!
    return(ind.genos)
  }
}

#' Make bed files from genotype/phenotype data
#'
#' From input genotypes, metadata, and phenotypes, creates a bed, ped, map, and bim file
#' for use in PLINK, GCTA, or other programs. Files of each type will be written to the current
#' working directory with the "data" prefix. Note that file named "data.sh" will be written which
#' must be ran via bash command line to produce the .bed binary file.
#'
#' @param x matrix or object coercable to a matrix. Genotypes, in phased format (two columns per individual).
#'   Genotypes should be in single number format (0 or 1), such as produced by \code{\link{process_ms}}.
#' @param meta data.frame of object coercable to data.frame. Metadata for snp data, where column 1 is the chromsome
#'   and column 2 is the position on the chromosome in bp,such as produced by \code{\link{process_ms}}.
#' @param phenos character vector. A vector containing the phenotypes for each individual,
#'   sorted identically to individuals in x.
#' @param return_objects logical, default FALSE, If TRUE, returns a list with ped, bed, map, and bim data. Otherwise just saves
#'   files.
#' @param missing_genotypes numeric, default -9. Encoding for missing alleles/genotypes.
#'
#' @return Either NULL or a list containing ped, bed, map, and bim formatted data.
#'
#' @author William Hemstrom
#'
#' @references Source code pulled from the snpR package by the same author for use here.
#'
#' @export
make_bed <- function(x, meta, phenos, plink_path = "/usr/bin/plink.exe", return_objects = F, missing_genotypes = -9){
  # get allele names for map file down the line
  a.names <- matrix(c(1,2), nrow(x), 2, T)

  #===============make a fam file=================
  cat("Straightening sheets (making ped and fam)... ")
  ped <- data.frame(fam = 1,
                    ind = paste0("Ind_", 1:(ncol(x)/2)),
                    PatID = 0,
                    MatID = 0,
                    Sex = 0,
                    Phenotype = phenos)

  # save .fam
  fam <- ped

  # change missing data value and add a space between alleles.
  x <- x + 1
  x[x == missing_genotypes + 1] <- 0
  x.np <- paste0(x[,seq(1,ncol(x), by = 2)], " ", x[,seq(2,ncol(x), by = 2)])
  x.np <- matrix(x.np, nrow = ncol(x)/2, byrow = T)

  # rebind
  ped <- cbind(as.data.table(ped), as.data.table(x.np), stringsAsFactors = F)

  # clean
  rm(x.np)
  gc()
  cat("All settled!\n\n")

  #===============make an extended map file=================
  # without morgans
  bim <- data.frame(chr = meta[,1],
                    rs = 1:nrow(meta),
                    mr = 0,
                    bp = meta[,2],
                    a1 = a.names[,1],
                    a2 = a.names[,2])

  # recode chr
  bim$chr <- as.numeric(as.factor(bim$chr))

  # grab normal map file
  map <- bim[,1:4]

  #============save files============
  cat("Fluffing pillows (saving plain text map, fam, ped, and bim files)... ")
  # save easy files
  data.table::fwrite(map, "data.map", quote = F, col.names = F, sep = " ", row.names = F, na = NA)
  data.table::fwrite(fam, "data.fam", quote = F, col.names = F, sep = " ", row.names = F, na = NA)
  data.table::fwrite(ped, "data.ped", quote = F, col.names = F, sep = " ", row.names = F, na = NA)
  data.table::fwrite(bim, "data.bim", quote = F, col.names = F, sep = " ", row.names = F, na = NA)

  cat("Smells like fresh shampoo!\n\n")

  # use plink to make a .bed file. Too much of a pain in the butt to write binary with R.
  cat("Smoothing out wrinkles (using PLINK to make a .bed file)...\n\n")
  system(paste0(plink_path, " --file data --make-bed --out data --allow-no-sex --aec"))
  cat("All done! Looks nice and cozy.\n\n")

  # name output
  if(return_objects){
    if(file.exists("data.bed")){cat(".bed succesfully made.\n\n")}
    else{cat("No .bed file located. Time to go to a mattress store. Check log for plink errors.\n\n")}
    return(list(ped = ped, map = map, bim = bim, fam = fam))
  }
  else{
    if(file.exists("data.bed")){cat(".bed succesfully made.\n");return(TRUE)}
    else{cat("No .bed file located. Check log for plink errors.\n");return(FALSE)}
  }
}

#' @export
estim_h <- function(x, meta,
                    num_threads = 1, autosome_num = length(unique(meta[,1])), maf = 0.01,
                    gcta_path = "/usr/bin/gcta.exe", plink_path = "/usr/bin/plink.exe",
                    phenos = NULL,
                    gcta_extra_arg_make_grm = NULL,
                    gcta_extra_arg_est_h = NULL){
  #===========sanity checks=========
  if(!all(c(file.exists("data.bed"), file.exists("data.map"), file.exists("data.bim")))){
    if(is.null(phenos)){
      stop("Phenos must be provided if no .bed, .map, and .bim files exist.")
    }
    cat("Making .bed...\n\n")
    make_bed(x, meta, phenos, plink_path = plink_path)
  }

  #===========construct call=======
  call <- paste0(gcta_path, " --bfile data --autosome --autosome-num ", autosome_num, " --maf ", maf,
                 " --make-grm --out data --thread-num ", num_threads)
  if(!is.null(gcta_extra_arg_make_grm)){
    call <- paste0(call, " ", gcta_extra_arg_make_grm)
  }

  #===========call gcta to make a grm============
  cat("Calling gcta to make a grm...\n\n")
  system(call)

  #===========call gcta to estimate heritability===========
  # make a .phen file
  cat("Making a .phen file...\n\n")
  data.table::fwrite(
    data.table::fread("data.fam")[j = c(1, 2, 6)],
    file = "data.phen",
    sep = " ")


  call2 <- paste0(gcta_path,
                  " --grm data --pheno data.phen --reml --out data --thread-num ",
                  num_threads)
  if(!is.null(gcta_extra_arg_est_h)){
    call2 <- paste0(call2, " ", gcta_extra_arg_est_h)
  }

  system(call2)

  #==========parse========================================
  res <- readr::read_delim("data.hsq", delim = "\t")
  return(list(res = res, call_make_grm = call, call_est_h = call2))
}

offset_overlapping_positions <- function(meta){
  meta[,2] <- round(meta[,2])
  meta$ord <- 1:nrow(meta)

  # function that checks for duplicates, returns dups, sorted meta, and pos_id
  dup_check <- function(meta){
    pos_id <- paste0(meta[,1], "_", meta[,2])
    meta <- meta[order(pos_id),]
    pos_id <- sort(pos_id)
    dups <- duplicated(pos_id) | rev(duplicated(rev(pos_id)))
    return(list(dups = dups, meta = meta, pos_id = pos_id))
  }

  # function that takes a dup check result and adjusts any duplicate sites
  adjust_pos <- function(dc){
    # unpack dc
    meta <- dc$meta
    dups <- dc$dups
    pos_id <- dc$pos_id

    # get a dup count table
    tab <- table(pos_id[dups])

    # make a new vector of positions for each duplicate site and overwrite
    ## prepare vectors
    cent <- ceiling(tab/2)
    current_pos <- meta[match(names(tab), pos_id), 2]
    lower <- current_pos - (cent - 1)
    upper <- ifelse(tab %% 2 == 0, current_pos + cent, current_pos + cent - 1)
    vec <- apply(matrix(c(lower, upper), ncol = 2), 1, function(x){x[1]:x[2]})

    ## overwrite and return
    meta[dups,2] <- unlist(vec)

    return(meta)
  }

  # run adjustment and re-dup check until there are no duplicated sites
  dc <- dup_check(meta)
  while(any(dc$dups)){
    meta <- adjust_pos(dc) # adjust
    dc <- dup_check(meta) # dup check
  }

  meta <- meta[order(meta$ord),]
  meta$ord <- NULL

  return(meta)
}

#' @export
make_vcf <- function(x, meta, missing_genotypes = -9){
  #==================convert to 0/0, 0/1, 1/1, or ./.=============
  ind.genos <- x[,seq(1,ncol(x), by = 2)] + x[,seq(2,ncol(x), by = 2)]
  tab <- data.frame(gt = c(0,1,2,missing_genotypes + missing_genotypes), recode = c("0/0", "0/1", "1/1", "./."))
  ind.genos <- tab[match(ind.genos, tab$gt),2]
  ind.genos <- matrix(ind.genos, ncol = ncol(x)/2)

  #==================add metadata to genos============
  # NOTE:
  # (vcf expects bp, may edit later to allow a meta containing that info to be imported)
  # for now, 0s will be A, 1s will be T, -1s will be .
  # for now, going to do this via conversion to 0 1 2, then converting. In the future will have a skip for phased data.

  vcf <- data.table::data.table(CHROM = as.numeric(factor(meta[,1], levels = unique(meta[,1]))),
                                POS = meta[,2],
                                ID = paste0("snp", 1:nrow(meta)),
                                REF = "A",
                                ALT = "T",
                                QUAL = ".",
                                FILTER = "PASS",
                                INFO = ".",
                                FORMAT = "GT"
  )
  colnames(ind.genos) <- paste0("SAMP_", 1:ncol(ind.genos))
  vcf <- cbind(vcf, as.data.table(ind.genos))
  colnames(vcf)[1] <- '#CHROM'

  writeLines("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=Integer,Description='Genotype'>\n##FORMAT=<ID=GP,Number=G,Type=Float,Description='Genotype Probabilities'>\n##FORMAT=<ID=PL,Number=G,Type=Float,Description='Phred-scaled Genotype Likelihoods'>", "data.vcf")
  data.table::fwrite(vcf, "data.vcf", sep = "\t", append = T, col.names = T, row.names = F, scipen = 999)
  return(vcf)
}

#' @export
impute_and_phase_beagle <- function(x = NULL, meta = NULL,
                                    beagle_path = "/usr/bin/beagle.jar",
                                    num_threads = 1,
                                    ne = 1000000,
                                    additional_args = NULL){
  #===============sanity checks===========
  # make a vcf if it doesn't exist
  make_vcf(x, meta)

  #===============construct call==========
  old.scipen <- getOption("scipen")
  options(scipen = 999)
  call <- paste0("java -jar ", beagle_path, " gt=data.vcf out=data.gt nthreads=",
                 num_threads, " ne=", ne)
  if(!is.null(additional_args)){
    call <- paste0(call, " ", additional_args)
  }

  #===============call beagle============
  system(call)
  options(scipen = old.scipen)


  #==============parse results===========
  # read in vcf and convert back to input format
  #res <- process_vcf("data.gt.vcf.gz", filter_non_poly = F)

  #return(res$genotypes)
  return(TRUE)
}

#' @export
process_vcf <- function(file, phased = T, filter_non_poly = T){
  # read in and parse
  if(Sys.info()[1] != "Windows"){
    res <- data.table::fread(cmd = paste("grep -v ##", file, "| tr '|' '\t'"), sep = "\t") # doesn't work on windows because it refuses to escape the pipe
  }
  else{
    newname <- paste0(file, ".sep")
    shell(paste("grep -v '##'", file, "| tr '|' '\t' >", newname))
    file <- newname
    res <- data.table::fread(file, sep = "\t", header = T)
  }
  header <- colnames(res)[-c(1:9)]
  meta <- res[,1:2]
  res <- res[,-c(1:9)]
  res <- res[, lapply(.SD, function(x) cbind(substr(x, 1, 1), substr(x, 3, 3)))]
  gc()
  res <- res[, lapply(.SD, as.integer)]
  colnames(res) <- paste(rep(header, each = 2), 1:2, sep = "_")
  colnames(meta) <- c("chr", "position")

  # fix things where the 1 is a minor.
  # not really a problem, just makes every - allele effect essentially a + and vice versa, but cleaner this way.
  min_flip <- which(rowMeans(res, na.rm = T) > 0.5)
  data.table::set(res, min_flip, 1:ncol(res), value = as.data.table(ifelse(res[min_flip,] == 0, 1, 0)))

  if(filter_non_poly){
    np <- which(rowSums(res, na.rm = T) == 0)
    if(length(np) != 0){
      res <- res[-np,]
      meta <- meta[-np,]
    }
  }

  if(!phased){
    res[res == "."] <- -4.5
    res <- matrix(as.numeric(res), nrow(res), ncol(res))
    res <- res[,seq(1,ncol(res), by = 2)] + res[,seq(2,ncol(res), by = 2)]
    colnames(res) <- header
  }

  # return
  return(list(genotypes = as.data.table(res), meta = meta))
}


#' Adds missing data to a genotype dataset.
#' @param x object coercable to numeric matrix containing genotype calls. Rows are loci, columns are samples. Each sample
#'   should be the genotype for a gene copy
#' @export
sprinkle_missing_data <- function(x,
                                  missing_dist = function(x){
                                    missing_dist <- rnorm(x, mean = 0.05, sd = 0.02);
                                    missing_dist[missing_dist < 0] <- 0;
                                    return(missing_dist)
                                  },
                                  missing_code = -9){
  # get missingness
  missingness <- missing_dist(nrow(x))

  # do a binom draw for missing data in each individual at each site
  missing_logi <- rbinom((ncol(x)/2)*nrow(x), size = 1, rep(missingness, each = ncol(x)/2))

  # arrange missing_logi like the data
  missing_logi <- matrix(as.logical(missing_logi), nrow(x), ncol(x)/2, byrow = T)

  # double the missing_logi to represent two gene copies per ind
  col_seq <- rep(1:ncol(missing_logi), each = 2)
  missing_logi <- missing_logi[,col_seq]

  # add missing loci
  x[missing_logi] <- missing_code

  return(x)
}

#' @export
assess_imputation <- function(x, x_imp, x_missing){
  # prep input data
  ## transpose
  xmc <- convert_2_to_1_column(x_missing)
  xc <- convert_2_to_1_column(x)
  xic <- convert_2_to_1_column(x_imp)

  # select missing and get metadata and mafs
  missing_genos <- xmc == -18
  xc <- xc[missing_genos]
  xic <- xic[missing_genos]
  per_loci_missing <- colSums(missing_genos)
  chr_matrix <- matrix(rep(meta[,1], each = nrow(xmc)), nrow = nrow(xmc))
  pos_matrix <- matrix(rep(meta[,2], each = nrow(xmc)), nrow = nrow(xmc))
  mafs <- rowMeans(x)
  mafs[mafs > .5] <- 1 - mafs[mafs > .5]
  maf_matrix <- matrix(rep(mafs, each = nrow(xmc)), nrow = nrow(xmc))

  # bind
  eval <- data.table(true_geno = xc, imputed_geno = xic,
                     chr = chr_matrix[missing_genos], pos = pos_matrix[missing_genos],
                     maf = maf_matrix[missing_genos])

  # get per loci cor
  cors <- eval[, .(rsq = cor(true_geno, imputed_geno)^2), by=.(chr, pos)]
  cors$maf <- eval[match(paste0(cors$chr, "_", cors$pos), paste0(eval$chr, "_", eval$pos)),5]

  # plot
  p <- ggplot2::ggplot(cors, ggplot2::aes(x = maf, y = rsq)) + ggplot2::geom_point(alpha = 0.05) + ggplot2::theme_gray()
  print(p)

  return(list(overall_rsq = cor(xc, xic)^2, per_loci_rsq = cors, plot = p))
}

#' Remove genotypes without phenotypic data.
#'
#' @param genotypes numeric matrix, data.table, or data.frame. Input genotypes, two columns per individual.
#' @param phenos numeric vector. Phenotypes for each individual. Missing data should be marked with NA.
#' @export
clean_phenotypes <- function(genotypes, phenos){
  missing <- which(is.na(phenos))
  if(length(missing) > 0){
    phenos <- phenos[-missing]
    iids <- rep(1:(ncol(genotypes)/2), each = 2)
    matches <- which(iids %in% missing)
    if(is.data.table(genotypes)){
      genotypes <- genotypes[,-..matches]
    }
    else{
      genotypes <- genotypes[,-matches]
    }
  }

  return(list(genotypes = genotypes, phenos = phenos))
}

#' Make a G matrix
#'
#' Make a G matrix using the Yang et al (2010) method from genotypes.
#'
#' @param ind.genos Unphased genotypes, rows are \emph{individuals}.
#' @export
make_G <- function(ind.genos, maf = 0.05,par = 1){
  colnames(ind.genos) <- paste0("m", 1:ncol(ind.genos)) # marker names
  rownames(ind.genos) <- paste0("s", 1:nrow(ind.genos)) # ind IDS
  mig <- min(ind.genos)
  G <- AGHmatrix::Gmatrix(ind.genos, missingValue = ifelse(mig == 0, NA, mig), method = "Yang", maf = maf)
  colnames(G) <- rownames(ind.genos)
  rownames(G) <- rownames(ind.genos)
  return(G)
}


#' Yang method G matrix from a SNP big matrix
#'
#' Make a G matrix using the method of Yang et al (2010) from a bigstatsr
#' FBM (File-Backed Matrix). Uses code adapted from the AGH matrix package.
#' Slower than simple vectorized code or the AGHmatrix function, but runs
#' on huge genotype files in parallel without using huge amounts of memory.
#'
#' @param SNPmatrix File-backed matrix (FBM), where each row is an \emph{individual} and each
#'   column is a \emph{loci}. Values are the intergers 0, 1, or 2, where 0 and 2 are homozygotes and 1 is
#'   a heterozygote.
#' @param maf numeric, default 0.05. Loci with minor allele frequencies below this will be removed.
#' @param par numeric, default 1. Number of cores to use for calculations.
make_yang_G_FBM <- function(SNPmatrix, maf = 0.05, par = 1){
  nloci <- ncol(SNPmatrix)
  Frequency <- bigstatsr::big_colstats(SNPmatrix)$sum/(2*nrow(SNPmatrix))
  Frequency <- cbind(1 - Frequency, Frequency)
  Frequency <- cbind(Frequency, matrixStats::rowMins(Frequency))

  if(maf != 0){
    rm <- which(Frequency[,3] <= maf)
    if(length(rm) > 0){
      SNPmatrix <- bigstatsr::big_copy(SNPmatrix, ind.col = 1:ncol(SNPmatrix)[-rm], type = "integer")
      Frequency <- Frequency[-rm, -3]
    }
  }

  # FBM_rep_by_parts <- function(to.rep, each, par = 1){
  #   X <- bigstatsr::FBM(nrow(to.rep), ncol = each)
  # }

  FreqP <- bigstatsr::FBM(nrow(SNPmatrix), ncol(SNPmatrix), init =
                            rep(Frequency[, 2], each = nrow(SNPmatrix)))
  FreqPQ <- bigstatsr::FBM(nrow(SNPmatrix), ncol(SNPmatrix),
                           init = rep(2 * Frequency[, 1] * Frequency[,2], each = nrow(SNPmatrix)))
  G.all <- bigstatsr::FBM(nrow(SNPmatrix), ncol(SNPmatrix), init =
                            bigstatsr::big_apply(SNPmatrix, a.FUN = function(X, ind)
                              (X[,ind]^2 - (1 + 2 * FreqP[,ind]) * SNPmatrix[,ind] + 2 * (FreqP[,ind]^2))/FreqPQ[,ind],
                              a.combine = cbind, ncores = par))


  G.ii <- as.matrix(bigstatsr::big_colstats(bigstatsr::big_transpose(G.all))$sum)
  SNPmatrix <- bigstatsr::FBM(nrow(SNPmatrix), ncol(SNPmatrix), init =
                                bigstatsr::big_apply(SNPmatrix, a.FUN = function(X, ind){
                                  y <- (SNPmatrix[,ind] - 2* FreqP[,ind])/sqrt(FreqPQ[,ind])
                                  y[is.na(y)] <- 0
                                  return(y)
                                }, a.combine = cbind, par = par))
  G.ii.hat <- 1 + (G.ii)/nloci
  Gmatrix <- bigstatsr::big_tcrossprodSelf(SNPmatrix)[1:nrow(SNPmatrix), 1:nrow(SNPmatrix)]/nloci
  diag(Gmatrix) <- G.ii.hat
  return(Gmatrix)
}

smart_transpose_and_unphase <- function(genotypes, phased){
  if(phased){
    return(convert_2_to_1_column(genotypes))
  }

  if("FBM" %in% class(genotypes)){
    return(bigstatsr::big_transpose(genotypes))
  }
  else{
    return(t(genotypes))
  }
}

fetch_phenotypes_ranger <- function(genotypes, model, h, a.var = NULL, phased = T){
  tgt <- smart_transpose_and_unphase(genotypes, phased)
  colnames(tgt) <- model$forest$independent.variable.names
  warning("Loci in genotypes assumed to be in same order as provided to random forest model.\n")
  a <- predict(model, data = tgt)$predictions
  if(is.null(a.var)){
    p <- a + e.dist.func(a, var(a), h = h)
  }
  else{
    p <- a + e.dist.func(a, a.var, h = h)

  }
  p <- list(a = a, p = p)
  return(p)
}


#' Names of the statistics calculated during distribution comparisons
#'
#' @format A character vector of length 32.
"names_diff_stats"

#' Names of the descriptive statistics calculated from GWAS results.
#'
#' @format A character vector of length 729.
"names_descriptive_stats"



# given parameter distributions (some joint, some not), iters, and joint distribution results, returns sampled values
sample_parameters_from_distributions <- function(parameter_distributions, iters, joint_res, joint_acceptance, joint_res_dist, reg_res, grid){
  # if any joint parameter priors, calculate and disambiguate
  joint_parms <- names(parameter_distributions)[which(parameter_distributions == "joint")]
  parms <- as.data.frame(matrix(NA, iters, 0))
  if(length(joint_parms) > 0){
    parms <- cbind(parms, gen_parms(iters, joint_res, joint_acceptance, joint_parms, dist.var = joint_res_dist, grid = grid))
  }
  if(!is.null(reg_res)){
    parms <- cbind(parms, sample_joint_quantile(iters, reg_res))
  }
  if(ncol(parms) < length(parameter_distributions)){
    other_parms <- names(parameter_distributions)[which(!names(parameter_distributions) %in% colnames(parms))]
    run_parameters <- vector("list", length = length(other_parms))
    names(run_parameters) <- other_parms
    for(i in 1:length(run_parameters)){
      run_parameters[[i]] <- parameter_distributions[[other_parms[i]]](iters)
    }
    parms <- cbind(parms, as.data.frame(run_parameters))
  }
  run_parameters <- parms

  return(as.data.frame(run_parameters))
}

.smart_split<- function(x, f, ...){
  if(is.null(x)){
    return(NULL)
  }
  else{
    return(split(x, f, ...))
  }
}
