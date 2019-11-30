ChowStatistic <- function(Sc, S1, S2, N1, N2, k){
  # get Chow test statistic
  # https://en.wikipedia.org/wiki/Chow_test
  numer <- ((Sc - (S1 + S2)) / k)
  denom <- (S1 + S2) / (N1 + N2 - 2 * k)
  return(numer  / denom)
}

DoChowTest <- function(fit.combined, fit.s1, fit.s2){
  # Do chow test: 
  # https://en.wikipedia.org/wiki/Chow_test
  Sc <- sum(fit.combined$residuals ^ 2)
  S1 <- sum(fit.s1$residuals ^ 2)
  S2 <- sum(fit.s2$residuals ^ 2)
  N1 <- nrow(fit.s1$model)
  N2 <- nrow(fit.s2$model)
  # check k is same for all fits
  ks <- sapply(list(fit.combined, fit.s1, fit.s2), function(fit.obj) fit.obj$rank)
  if (length(unique(ks)) != 1){
    stop("Number of parameters must be same for all fits")
  }
  k <- fit.combined$rank
  chow.t.stat <- ChowStatistic(Sc, S1, S2, N1, N2, k)
  # get p-value
  pval <- 1 - pf(q = chow.t.stat, df1 = k, df2 = N1 + N2 - 2 * k)
  return(list(t.stat = chow.t.stat, pval = pval))
}
