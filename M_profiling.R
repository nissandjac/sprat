library(smsR)

wd <- "C:/Users/nsja/Dropbox/DTU/BEBRIS/two_seasons/"


# Mortality profiling 


mort_prof <- function(sas, 
                      M_range = seq(0.5, 1.2, length.out = 10)){
  
  df.tmb <- sas$dat 
  
  nll <- rep(0, length.out = length(M_range))
  
  for(i in 1:length(M_range)){
    
  df.tmp <- df.tmb  
  df.tmp$M  <- df.tmp$M * M_range[i]  
    
  parms <- getParms(df.tmp)  
  sas_tmp <- runAssessment(df.tmp, parms)
    
  nll[i] <- sas_tmp$opt$objective
  
  
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}