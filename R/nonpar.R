# all of these functions come from tpepler/nonpar, which doesn't export them for some reason. Should re-implment.


cucconi.test <- function(x, y, method = c("permutation", "bootstrap")){

  # Implementation of the Cucconi test for the two-sample location-scale problem
  # A permutation/bootstrap distribution of the test statistic (C) under the
  # null hypothesis is used to calculate the p-value.
  # Reference: Marozzi (2013), p. 1302-1303

  m <- length(x)
  n <- length(y)
  C <- cucconi.teststat(x = x, y = y, m = m, n = n)

  if(method[1] == "permutation"){
    h0dist <- cucconi.dist.perm(x = x, y = y)
  }

  if(method[1] == "bootstrap"){
    h0dist <- cucconi.dist.boot(x = x, y = y)
  }

  p.value <- length(h0dist[h0dist >= C]) / length(h0dist)

  cat("\nCucconi two-sample location-scale test\n")
  cat("\nNull hypothesis: The locations and scales of the two population distributions are equal.\n")
  cat("Alternative hypothesis: The locations and/or scales of the two population distributions differ.\n")
  cat(paste("\nC = ", round(C, 3), ", p-value = ", round(p.value, 4), "\n\n", sep=""))

  return(list(C = C,
              method = method[1],
              p.value = p.value))
}



# from tpepler/nonpar. No longer built with package.
lepage.test <- function(x, y = NA, g = NA, method = NA, n.mc = 10000){

  ##Adapted from kruskal.test()##
  if (is.list(x)) {
    if (length(x) < 2L)
      stop("'x' must be a list with at least 2 elements")
    y <- x[[2]]
    x <- x[[1]]
  }
  else {
    if(min(is.na(y)) != 0){
      k <- length(unique(g))
      if (length(x) != length(g))
        stop("'x' and 'g' must have the same length")
      if (k < 2)
        stop("all observations are in the same group")
      y <- x[g == 2]
      x <- x[g == 1]
    }
  }
  #####################

  outp <- list()
  outp$m <- m <- length(x)
  outp$n <- n <- length(y)
  outp$n.mc <- n.mc
  N <- outp$m + outp$n
  outp$ties <- (length(c(x, y)) != length(unique(c(x, y))))
  even <- (outp$m + outp$n + 1)%%2
  outp$stat.name <- "Lepage D"


  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(choose(outp$m + outp$n, outp$n) <= 10000){
      method <- "Exact"
    }
    if(choose(outp$m + outp$n, outp$n) > 10000){
      method <- "Monte Carlo"
    }
  }
  #####################################################################
  outp$method <- method

  tmp.W <- rank(c(x, y))

  our.data <- rbind(c(x, y), c(rep(1, length(x)), rep(0, length(y))))
  sorted <- our.data[1, order(our.data[1, ]) ]
  x.labels <-our.data[2, order(our.data[1, ]) ]

  med <- ceiling(N / 2)
  if(even){no.ties <- c(1:med, med:1)}
  if(!even){no.ties <- c(1:med, (med - 1):1)}

  obs.group <- numeric(N)
  group.num <- 1
  for(i in 1:N){
    if(obs.group[i] == 0){
      obs.group[i] <- group.num
      for(j in i:N){
        if(sorted[i] == sorted[j]){
          obs.group[j] <- obs.group[i]
        }
      }
      group.num <- group.num + 1;
    }
  }

  group.ranks <- tapply(no.ties, obs.group, mean)

  tied.ranks <- numeric(N)
  for(i in 1:group.num){
    tied.ranks[which(obs.group == as.numeric(names(group.ranks)[i]))] <- group.ranks[i]
  }

  tmp.C <- c(tied.ranks[x.labels == 1], tied.ranks[x.labels == 0])

  ##Only needs to depend on y values
  D.calc <- function(C.vals, W.vals){

    if(even){
      exp.C <- n * (N + 2) / 4
      var.C <- m * n * (N + 2) * (N - 2) / (48 * (N - 1))
    }
    if(!even){
      exp.C <- n * (N + 1)^2 / (4 * N)
      var.C <- m * n * (N + 1) * (3 + N^2) / (48 * N^2)
    }
    W.obs <- sum(W.vals)
    W.star <- (W.obs - n * (N + 1) / 2) / sqrt(m * n * (N + 1) / 12)
    C.star <- (sum(C.vals) - exp.C) / sqrt(var.C)
    return(W.star^2 + C.star^2)
  }

  outp$obs.stat <- D.calc(tmp.C[(m + 1):N], tmp.W[(m + 1):N])

  if(outp$method == "Exact"){
    possible.orders <- gtools::combinations(outp$m + outp$n, outp$n)

    possible.C <- t(apply(possible.orders, 1, function(x) tmp.C[x]))
    possible.W <- t(apply(possible.orders, 1, function(x) tmp.W[x]))

    theor.dist <- numeric(nrow(possible.C))
    for(i in 1:nrow(possible.C)){
      theor.dist[i] <- D.calc(possible.C[i, ], possible.W[i, ])
    }

    outp$p.value <- mean(theor.dist >= outp$obs.stat)
  }

  if(outp$method == "Asymptotic"){
    outp$p.value <- (1 - pchisq(outp$obs.stat, 2))
  }

  if(outp$method == "Monte Carlo"){
    outp$p.value <- 0
    for(i in 1:n.mc){
      mc.sample <- sample(1:N, n)

      if(D.calc(tmp.C[mc.sample], tmp.W[mc.sample]) >= outp$obs.stat){
        outp$p.value = outp$p.value + 1 / n.mc
      }
    }
  }

  cat("\nLepage two-sample location-scale test\n")
  cat("\nNull hypothesis: The locations and scales of the two population distributions are equal.\n")
  cat("Alternative hypothesis: The locations and/or scales of the two population distributions differ.\n")
  cat(paste("\nD = ", round(outp$obs.stat, 3), ", p-value = ", round(outp$p.value, 4), "\n\n", sep=""))

  #class(outp)="NSM3Ch5p"
  return(outp)
}



cucconi.dist.perm <- function(x, y, reps = 1000){

  # Computes the distribution of the Cucconi test statistic using random permutations

  m <- length(x)
  n <- length(y)
  N <- m + n
  alldata <- c(x, y)

  bootFunc <- function(){
    permdata <- alldata[sample(1:N, size = N, replace = FALSE)]
    xperm <- permdata[1:m]
    yperm <- permdata[(m + 1):N]
    return(cucconi.teststat(x = xperm, y = yperm, m = m, n = n))
  }
  permvals <- replicate(reps, expr = bootFunc())

  #  permvals<-rep(NA,times=reps)
  #  for(r in 1:reps){
  #    permdata<-alldata[sample(1:N,size=N,replace=FALSE)]
  #    xperm<-permdata[1:m]
  #    yperm<-permdata[(m+1):N]
  #    permvals[r]<-cucconi.teststat(x=xperm,y=yperm, m=m, n=n)
  #  }
  return(permvals)
}


cucconi.dist.boot <- function(x, y, reps = 10000){

  # Computes the distribution of the Cucconi test statistic using bootstrap sampling

  m <- length(x)
  n <- length(y)
  x.s <- (x - mean(x)) / sd(x) # standardise the x-values
  y.s <- (y - mean(y)) / sd(y) # standardise the y-values

  bootFunc <- function(){
    xboot <- x.s[sample(1:m, size = m, replace = TRUE)]
    yboot <- y.s[sample(1:n, size = n, replace = TRUE)]
    return(cucconi.teststat(x = xboot, y = yboot, m = m, n = n))
  }
  bootvals <- replicate(reps, expr = bootFunc())

  #  bootvals <- rep(NA,times=reps)
  #  for(r in 1:reps){
  #    xboot<-x.s[sample(1:m,size=m,replace=TRUE)]
  #    yboot<-y.s[sample(1:n,size=n,replace=TRUE)]
  #    bootvals[r]<-cucconi.teststat(x=xboot,y=yboot, m=m, n=n)
  #  }
  return(bootvals)
}

cucconi.teststat <- function(x, y, m = length(x), n = length(y)){

  # Calculates the test statistic for the Cucconi two-sample location-scale test

  N <- m + n
  S <- rank(c(x, y))[(m + 1):N]
  denom <- sqrt(m * n * (N + 1) * (2 * N + 1) * (8 * N + 11) / 5)
  U <- (6 * sum(S^2) - n * (N + 1) * (2 * N + 1)) / denom
  V <- (6 * sum((N + 1 - S)^2) - n * (N + 1) * (2 * N + 1)) / denom
  rho <- (2 * (N^2 - 4)) / ((2 * N + 1) * (8 * N + 11)) - 1
  C <- (U^2 + V^2 - 2 * rho * U * V) / (2 * (1 - rho^2))
  return(C)
}
