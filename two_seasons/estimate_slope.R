# Compare estimatng the slope vs inputting it 

library(smsR)

wd <- "C:/Users/nsja/Dropbox/DTU/BEBRIS/two_seasons/"

# M is not correct here 
maxage <- 3
years = 1974:2024
nyear <- length(years)
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

Surveyobs <- survey_to_matrix(dat[['survey']], years)
Catchobs <- df_to_matrix(dat[['canum']], season =  1:2)

nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#')#, skip = 3)
# Calculate the catch yield per season 
dat$nocatch <- as.matrix(nocatch)#*0+1

df.tmb <- get_TMB_parameters(
  mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  years = years, # Years to run
  #endYear = 2023,
  nseason = nseason, # Number of seasons
  ages = ages, # Ages of the species
  recseason = 1, # Season where recruitment occurs
  CminageSeason = c(0,0),
  Fmaxage = 2, # Fully selected fishing mortality age
  Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
  Qmaxage = Qmaxage, #Qmaxage = c(1,3)
  # minSDsurvey = .2,
  # randomF = 1,
  blocks = c(1974,2015),
  endFseason = 1, # which season does fishing stop in the final year of data
  nocatch = as.matrix(dat$nocatch),
  surveyStart = surveyStart, #c(0.75,0)
  surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
  surveySD =  surveyCV, #c(1,2)),
  catchSD = list(c(0,1,2,3),
                 c(0,1,2,3)),
  estSD = c(0,0,0), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,1) # Factor for relative strength of log-likelihood
  
)

parms <- getParms(df.tmb )
mps <-getMPS(df.tmb, parms)# Set boundaries
# This model works best if SDsurvey is mapped
sas <- runAssessment(df.tmb, parms = parms,mps = mps, silent = TRUE)

# 
mps$logbeta <- NULL
sas_beta <- runAssessment(df.tmb, parms = parms,mps = mps, silent = TRUE)





