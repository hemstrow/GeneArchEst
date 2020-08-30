

# functions to find and compare peaks
findpeaks <- function(x, delta, pcut = .005, pvals = T){
  peakdet <- function(v, delta, x){
    # from: https://gist.github.com/dgromer/ea5929435b8b8c728193
    maxtab <- NULL
    mintab <- NULL

    if (length(v) != length(x))
    {
      stop("Input vectors v and x must have the same length")
    }

    if (!is.numeric(delta))
    {
      stop("Input argument delta must be numeric")
    }

    if (delta <= 0)
    {
      stop("Input argument delta must be positive")
    }

    mn <- Inf
    mx <- -Inf

    mnpos <- NA
    mxpos <- NA

    lookformax <- TRUE

    for(i in seq_along(v))
    {
      this <- v[i]

      if (this > mx)
      {
        mx <- this
        mxpos <- x[i]
      }

      if (this < mn)
      {
        mn <- this
        mnpos <- x[i]
      }

      if (lookformax)
      {
        if (this < mx - delta)
        {
          maxtab <- rbind(maxtab, data.frame(pos = mxpos, val = mx))

          mn <- this
          mnpos <- x[i]

          lookformax <- FALSE
        }
      }
      else
      {
        if (this > mn + delta)
        {
          mintab <- rbind(mintab, data.frame(pos = mnpos, val = mn))

          mx <- this
          mxpos <- x[i]

          lookformax <- TRUE
        }
      }
    }

    list(maxtab = maxtab, mintab = mintab)
  }

  if(pvals){
    lpcut <- -log10(pcut)
    x$effect <- -log10(x$effect)
    peaks <- x[x$effect >= lpcut,]
    peaks <- na.omit(peaks)
    xs <- smooth(peaks$effect)
    peaks$slogp <- as.numeric(xs)
  }
  else{
    peaks <- x[x$effect >= pcut[1] | x$effect <= pcut[2],]
    peaks <- na.omit(peaks)
    xs <- smooth(abs(peaks$effect))
    peaks$seffect <- as.numeric(xs)
  }
  tpk <- peakdet(as.numeric(xs), delta = delta, x = peaks$position)
  return(tpk$maxtab)
}

# wrapper for findpeaks to find across multiple chromosomes
#' @export
findpeaks_multi <- function(x, delta, pcut, chr = "chr", pvals = T){
  ug <- unique(x[,chr])
  lout <- vector("list", length(ug))
  for(i in 1:length(lout)){
    lout[[i]] <- as.data.table(findpeaks(x[x[,chr] == ug[i],], delta = delta, pcut = pcut, pvals = pvals))
    if(nrow(lout[[i]] > 0)){
      lout[[i]]$chr <- ug[i]
    }
  }
  return(data.table::rbindlist(lout))
}
