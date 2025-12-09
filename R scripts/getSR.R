getSR <- function(df.tmb, sas, method = 'mean') {
  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values <- rownames(sdrep)

  years <- df.tmb$years
  SSB <- getSSB(df.tmb, sas)
  R <- getR(df.tmb, sas)


  if(df.tmb$recmodel == 1){
    # log(alpha) and its SE
    SRstock <- sdrep[rep.values == "logalpha", 1]
    SE      <- sdrep[rep.values == "logalpha", 2]

    # method should be a character, e.g. method <- "mean" or "median"
    if (method == "mean") {
      # process variance on log scale for mean of lognormal
      sigmaR  <- exp(sdrep["logSDrec", "Estimate"])  # SD on log scale
      sigma2  <- sigmaR^2                             # variance on log scale
      adj        <- 0.5 * sigma2          # shift for mean curve
    } else {
      # median: no shift from lognormal
      adj <- 0
    }

    # CI for log(alpha)
    low  <- SRstock - 2 * SE
    high <- SRstock + 2 * SE

    # SSB range and beta
    SSBrange <- seq(1, max(SSB$SSB, na.rm = TRUE), length.out = nrow(R))

    # see if beta is estimated
    beta_est <- exp(sas$reps$par.fixed['logbeta'])
    if(is.na(beta_est)){
      betaSR   <- as.numeric(df.tmb$betaSR)
    }else{
      betaSR <- beta_est
    }
    # main hockey-stick curve (median if adj = 0, mean if adj = 0.5*sigma2_hat)
    SR <- exp(SRstock + log(SSBrange) + adj)
    SR[SSBrange > betaSR] <- exp(SRstock + log(betaSR) + adj)

    # curves using log(alpha) ± 2*SE (parameter uncertainty band)
    mins <- exp(low + log(SSBrange) + adj)
    mins[SSBrange > betaSR] <- exp(low + log(betaSR) + adj)

    maxs <- exp(high + log(SSBrange) + adj)
    maxs[SSBrange > betaSR] <- exp(high + log(betaSR) + adj)

  }

  if(df.tmb$recmodel == 2){

    SRstock <- sdrep[rep.values == "logalpha", 1]
    SE <- (sdrep[rep.values == "logalpha", 2])
    mins <- SRstock - 2 * SE
    maxs <- SRstock + 2 * SE


  }

  # Get the estimated recruitment

  SR.out <- data.frame(SR = SR, R = R$R,
                       low = R$low,
                       high = R$high,
                       minSR = mins, maxSR = maxs, SSB = SSBrange)

  return(SR.out)
}
