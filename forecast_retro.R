# Simulate all recruitment options     
library(smsR)


dats.comb <- expand.grid(
  power = c(0, NA),
  SRR   = c(0.1, 1),
  SDR   = c(0, 2)
)

library(purrr)
library(dplyr)

library(doParallel)
library(foreach)

# Number of workers (leave 1 core free)
n_cores <- max(1, parallel::detectCores() - 1)

cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# Optional: preallocate name for each combo
combo_names <- apply(dats.comb, 1, function(x) {
  paste(paste(names(dats.comb), x, sep = "="), collapse = "_")
})

# Run getAssess() in parallel for each row
results_list <- foreach(i = 1:nrow(dats.comb),
                        .packages = c('smsR'),      # add any packages getAssess() needs
                        .export   = c("getAssess2", "dats.comb"),  # make sure these are visible on workers
                        .combine  = 'list',
                        .multicombine = TRUE) %dopar% {
                          # One row of settings
                          this_row <- dats.comb[i, ]
                          
                          # Run your assessment
                          res <- getAssess2(this_row)
                          
                          # (Optional) attach the settings to the result for easier tracking
                          res$settings <- this_row
                          
                          # (Optional) save each run to disk immediately:
                          # saveRDS(res, file = sprintf("assess_result_%03d.rds", i))
                          
                          res
                        }

# Give each result a meaningful name
names(results_list) <- combo_names

# Clean up cluster
parallel::stopCluster(cl)


saveRDS(results_list, file = 'recruitment_parameters/forecast_sims_two.RDS')


getAssess <- function(dats, npeels =15){
  
  wd <- "C:/Users/nsja/Dropbox/DTU/BEBRIS/four_quarters/"
  
  maxage <- 3
  years = 1974:2024
  #nyear <- length(years)
  seasons <- 1:4
  nseason <- length(seasons)
  
  dat <- getDataSMS(wd,
                    maxage = maxage,
                    survey.age = list(0:3, 1:3, 1:3), # Ages in the two surveys
                    survey.years = list(1982:2024, 1992:2024, 2006:2024),
                    survey.names = c('IBTS_Q1','IBTS_Q3', 'HERAS'),
                    survey.quarter = c(3,1, 1),
                    years = years,
                    seasons = seasons
                    
  )
  Qminage = c(0,1,1) # Qminage = c(0,1) minimum age in surveys
  Qmaxage = c(3,3,3) #Qmaxage = c(1,3)
  surveyStart = c(0,0,0) #c(0.75,0)
  surveyEnd =  c(1,1,0) # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = c(3,1,1) # c(2,1)Which seasons do the surveys occur in
  surveyCV =  list(c(0,1),
                   c(1,2),
                   c(1,2)) #c(1,2)),
  
  # Load packages and files #
  
  
  ages <- 0:maxage
  beta <- 90000
  
  # Fishing eff
  Surveyobs <- survey_to_matrix(dat[['survey']], years)
  # Rescale the recruitment stuff
  Catchobs <- df_to_matrix(dat[['canum']], season =  1:4)
  nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#')#, skip = 3)
  #Feffort[which(df.tmb$years == 2013),2] <- 0 # This is different two places in the input files
  
  # Save some data for package example
  dat$nocatch <- as.matrix(nocatch)#*0+1
  powerIN <- c(dats$power, NA, NA)  
  
  #dat$nocatch[1,4] <- 0
  SSBest <- rep(NA, npeels)
  
  for(peels in 0:npeels){
  
  
  df.tmb <- get_TMB_parameters(
    mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
    Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
    Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
    years = years, # Years to run
    endYear = years[length(years)-peels],
    nseason = nseason, # Number of seasons
    ages = ages, # Ages of the species
    recseason = 1, # Season where recruitment occurs
    CminageSeason = c(0,0,0,0),
    Fmaxage = 2, # Fully selected fishing mortality age
    Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
    Qmaxage = Qmaxage, #Qmaxage = c(1,3)
    # minSDcatch = sqrt(0.01),
    # maxSDcatch = sqrt(1.5),
    # minSDsurvey = sqrt(0.2),
    #penepsC = 1e-10,
    penepsCmax = 1e-8,
    # peneps = 1e-10,
    Fbarage = c(1,2),
    powers = powerIN,
    endFseason = 2, # which season does fishing stop in the final year of data
    nocatch = as.matrix(nocatch),
    surveyStart = surveyStart, #c(0.75,0)
    surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
    surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
    surveySD =  surveyCV, #c(1,2)),
    catchSD = list(c(0,1,2),
                   c(0,1,2),
                   c(0,1,2),
                   c(0,1,2)),
    estSD = c(0,0,dats$SDR), # Estimate
    beta = beta, # Hockey stick plateau
    nllfactor = c(1,1,dats$SRR) # Factor for relative strength of log-likelihood
    
  )
  
  # Get initial parameter structure
  parms <- getParms(df.tmb)
  # Get non-estimated parameters, based on info in df.tmb
  
  sas <- runAssessment(df.tmb, parms = parms, silent = TRUE)
  
  
  SSB <- getSSB(df.tmb, sas)
  SSBest[peels] <- SSB$SSB[nrow(SSB)]
  
  if(peels == 0){
    out <- getSummary(df.tmb, sas)
    }
  }
  
  info <- list(SSBest = SSBest,
               parameters = dats)
  
  
  return(list(out = out, info = info))
  
}


getAssess2 <- function(dats, npeels =15){
  
  wd <- "C:/Users/nsja/Dropbox/DTU/BEBRIS/two_seasons/"
  maxage <- 3
  years = 1974:2024
  #nyear <- length(years)
  seasons <- 1:2
  nseason <- length(seasons)
  
  dat <- getDataSMS(wd,
                    maxage = maxage,
                    survey.age = list(0:3, 1:3, 1:3), # Ages in the two surveys
                    survey.years = list(1982:2024, 1992:2024, 2006:2024),
                    survey.names = c('IBTS_Q1','IBTS_Q3', 'HERAS'),
                    survey.quarter = c(3,1, 1),
                    years = years,
                    seasons = seasons
                    
  )
  Qminage = c(0,1,1) # Qminage = c(0,1) minimum age in surveys
  Qmaxage = c(3,3,3) #Qmaxage = c(1,3)
  surveyStart = c(0,0.17,0) #c(0.75,0)
  surveyEnd =  c(.5,.33,0.17) # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = c(2,1,1) # c(2,1)Which seasons do the surveys occur in
  surveyCV =  list(c(0,1),
                   c(1,2),
                   c(1,2)) #c(1,2)),
  
  # Load packages and files #
  
  
  ages <- 0:maxage
  beta <- 90000
  
  # Fishing eff
  Surveyobs <- survey_to_matrix(dat[['survey']], years)
  # Rescale the recruitment stuff
  Catchobs <- df_to_matrix(dat[['canum']], season =  1:2)
  nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#')#, skip = 3)
  #Feffort[which(df.tmb$years == 2013),2] <- 0 # This is different two places in the input files
  
  # Save some data for package example
  dat$nocatch <- as.matrix(nocatch)#*0+1
  powerIN <- c(dats$power, NA, NA)  
  #dat$nocatch[1,4] <- 0
  SSBest <- rep(NA, npeels)
  
  for(peels in 0:npeels){
    
    
    df.tmb <- get_TMB_parameters(
      mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
      Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
      Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
      years = years, # Years to run
      endYear = years[length(years)-peels],
      nseason = nseason, # Number of seasons
      ages = ages, # Ages of the species
      recseason = 1, # Season where recruitment occurs
      CminageSeason = c(0,0),
      Fmaxage = 2, # Fully selected fishing mortality age
      Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
      Qmaxage = Qmaxage, #Qmaxage = c(1,3)
      # minSDcatch = sqrt(0.01),
      # maxSDcatch = sqrt(1.5),
      # minSDsurvey = sqrt(0.2),
      #penepsC = 1e-10,
      penepsCmax = 1e-8,
      # peneps = 1e-10,
      Fbarage = c(1,2),
#      isFseason = c(1,1,1,0), # Seasons to calculate fishing in
      powers = powerIN,
      endFseason = 2, # which season does fishing stop in the final year of data
      nocatch = as.matrix(nocatch),
      surveyStart = surveyStart, #c(0.75,0)
      surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
      surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
      surveySD =  surveyCV, #c(1,2)),
      catchSD = list(c(0,1,2),
                     c(0,1,2),
                     c(0,1,2),
                     c(0,1,2)),
      estSD = c(0,0,dats$SDR), # Estimate
      beta = beta, # Hockey stick plateau
      nllfactor = c(1,1,dats$SRR) # Factor for relative strength of log-likelihood
      
    )
    
    # Get initial parameter structure
    parms <- getParms(df.tmb)
    # Get non-estimated parameters, based on info in df.tmb
    
    sas <- runAssessment(df.tmb, parms = parms, silent = TRUE)
    
    
    SSB <- getSSB(df.tmb, sas)
    SSBest[peels] <- SSB$SSB[nrow(SSB)]
    
    if(peels == 0){
      out <- getSummary(df.tmb, sas)
    }
  }
  
  info <- list(SSBest = SSBest,
               parameters = dats)
  
  
  return(list(out = out, info = info))
  
}
