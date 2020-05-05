#' Marks unique windows for each snp.
#' @param x data.frame. Must contain a "position" column and a chromosome info column.
#' @param sigma numeric. Size of windows, in kb.
#' @param chr character, default "chr". Name of chromosome info column in x.
mark_windows <- function(x, sigma, chr = "chr"){
  sigma <- sigma * 1000

  # window start and end points
  unique.chr <- sort(unique(x[,chr]))
  chr.max <- tapply(x$position, x[,chr], max)
  chr.min <- tapply(x$position, x[,chr], min)
  chr.range <- matrix(c(chr.max, chr.min), ncol = 2)
  starts <- apply(chr.range, 1, function(x) seq(from = 0, to = x[1], by = sigma))
  if(is.matrix(starts)){
    starts <- as.list(as.data.frame(starts))
  }
  ends <- foreach::foreach(q = 1:length(chr.max), .inorder = T) %do%{
    c(starts[[q]][-1], chr.max[q])
  }

  # for each chr, assign windows to all snps
  ## function per chr
  assign_windows <- function(y, starts, ends){
    lmat <- outer(y$position, starts, function(pos, starts) pos > starts)
    lmat <- lmat * outer(y$position, ends, function(pos, ends) pos <= ends)
    lmat[lmat == 1] <- rep(1:ncol(lmat), each = nrow(lmat))[lmat == 1]
    if(any(colSums(lmat) == 0)){
      warning("Colsums are 0.")
    }
    lmat <- as.numeric(lmat[lmat != 0])
    return(lmat)
  }
  ## note: need to make sure that the window ids are unique!
  windows <- vector("list", length(unique.chr))
  tot.windows <- 0
  for(i in 1:length(chr.max)){
    windows[[i]] <- assign_windows(x[x[,chr] == unique.chr[i],], starts[[i]], ends[[i]]) + tot.windows
    tot.windows <- tot.windows + length(unique(windows[[i]]))
  }
  windows <- unlist(windows)

  return(windows)
}



#' Return discriptive statistics about effect sizes/p-values.
#'
#' Calculates a range of descriptive statistics both about the p-values as a whole,
#' effect size peaks, and sliding windows.
#'
#' Peaks are identified with a peak finding algorithm based on the delta and pcut values.
#' These values should be chosen based on the true p-value or effect size distribution. p-values are
#' -log10 transformed.
#'
#' Stats about windows are run through the basic stat calculation twice: once per window (means per
#' window, for example), then once *per stat* (so the mean of the mean values per window).
#'
#' The following stats are calculated: \itemize{\item{mean} \item{median} \item{sd} \item{variance}
#' \item{20 evenly spaced quantiles between 0 and 1} \item{kurtosis} \item{skewness}}.
#'
#' In addition, the number of unique peaks is also calculated.
#'
#' @param p numeric. p-values or effect sizes for each SNP.
#' @param meta. data.frame. The position and chromosome for each SNP. Needs a column named "position".
#' @param delta numeric. Value used to determine spacing between peaks.
#' @param pcut numeric, default 0.005. Only p-values/effect sizes above/below this quantile will be used for peak detection.
#' @param chr character,default "chr". Name of the meta column containing chromosome info.
#' @param pvals logical, default TRUE. Specifies if p-values are provided. p-values will be -log10 transformed.
dist_desc <- function(p, meta, windows, delta, pcut = .005, chr = "chr", pvals = T){
  #==============subfunctions================================
  basic <- function(p){
    # quantiles
    p <- na.omit(p)
    quantsp <- quantile(p, probs = seq(0, 1, length.out = 22)[-c(1, 22)])
    names(quantsp) <- paste0("Quantile_", names(quantsp))

    # shape
    out_desc_p <- c(kurtosis = e1071::kurtosis(p),
                    skewness = e1071::skewness(p),
                    mean = mean(p), median = median(p),
                    sd = sd(p), var = var(p))
    return(c(quantsp, out_desc_p))
  }

  #==============peaks=======================================
  # find peaks
  peaks <- findpeaks_multi(cbind(meta[,c(chr, "position")], effect = p), delta, pcut, chr, pvals)
  peak_stats <- basic(peaks$val)
  names(peak_stats) <- paste0("peak_", names(peak_stats))
  peak_stats <- c(peak_stats, peak_count = nrow(peaks))

  #===============stats for the basic distribution===========
  if(pvals){
    p <- -log10(p)
  }

  dist_stats <- basic(p)

  #==============windows======================================
  window_meta <- cbind(data.table::as.data.table(meta[,c(chr, "position")]), effect = p, window = windows)
  window_meta <- window_meta[,as.list(basic(effect)), by = "window"]
  colnames(window_meta)[-1] <- paste0("window_", colnames(window_meta)[-1])
  window_meta <- data.table::melt(window_meta, id.vars = "window")
  window_meta <- window_meta[,as.list(basic(value)), by = "variable"]
  colnames(window_meta)[1] <- "window_stat"
  window_meta <- data.table::melt(window_meta, id.vars = "window_stat")
  window_meta$stat <- paste0(window_meta$window_stat, "_summary_", window_meta$variable)
  window_stats <- window_meta[["value"]]
  names(window_stats) <- window_meta[["stat"]]

  return(c(dist_stats, peak_stats, window_stats))
}



compare_distributions <- function(o, p){
  # descriptive for p
  quantsp <- quantile(p, probs = seq(0.1, 0.9, length.out = 20))
  names(quantsp) <- paste0("Quantile_", names(quantsp))
  out_desc_p <- c(kurtosis = e1071::kurtosis(p),
                  skewness = e1071::skewness(p),
                  mean = mean(p), median = median(p), sd = sd(p), var = var(p),
                  quantsp)
  names(out_desc_p) <- paste0("sim_", names(out_desc_p))

  # descriptive for o
  quantso <- quantile(o, probs = seq(0.1, 0.9, length.out = 20))
  names(quantso) <- paste0("Quantile_", names(quantso))
  out_desc_o <- c(kurtosis = e1071::kurtosis(o),
                  skewness = e1071::skewness(o),
                  mean = mean(o), median = median(o), sd = sd(o), var = var(o),
                  quantso)
  names(out_desc_o) <- paste0("observed_", names(out_desc_o))


  # comparative
  capture.output(invisible((lepage <- nonpar::lepage.test(p, o))))
  capture.output(invisible(cucconi <- nonpar::cucconi.test(p, o, "permutation")))
  suppressMessages(out_dist <- c(ks = unlist(ks.test(p, o)$statistic),
                                 norm.ks = tryCatch(ks.test(p/sum(p), o/sum(o))$statistic, error=function(err) NA),
                                 lepage.p = lepage$p.value, lepage.s = lepage$obs.stat,
                                 cucconi.p = cucconi$p.value, cucconi.s = cucconi$C))
  out <- c(abs(out_desc_o - out_desc_p), out_dist)
  names(out) <- names_diff_stats
  return(out)
}
