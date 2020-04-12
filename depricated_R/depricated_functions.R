



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





calc_dist_stats <- function(o, p){
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

  return(c(abs(out_desc_o - out_desc_p), out_dist))
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
