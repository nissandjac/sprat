# Run the comp data # 

library(smsR)
library(tidyverse)
wd <- "C:/Users/nsja/Dropbox/DTU/BEBRIS/propcatch/"
source(file.path(wd,'getSurvey.R'))

maxage <- 3
years = 1974:2024
nyear <- length(years)
seasons <- 1:2
nseason <- length(seasons)
ages <- 0:maxage
beta <- 90000
nage <- length(ages)

Qminage = c(0,1,1) # Qminage = c(0,1) minimum age in surveys
Qmaxage = c(3,3,2) #Qmaxage = c(1,3)

surveyCV =  list(c(0,1,2),
                 c(1,2,3),
                 c(1)
) 

# Manually read the survey data 
survey <- read.table(file.path(wd, "fleet_catch.in"), 
                     fill = TRUE)
survey.in <- getSurvey(survey,
                       maxage = maxage,
                       survey.age = list(0:3, 1:3, 1:3), # Ages in the two surveys
                       survey.years = list(1982:2024, 1992:2024, 2006:2024),
                       survey.names = c('IBTS_Q1','IBTS_Q3', 'HERAS'),
                       survey.quarter = c(3,1, 1),
                       years = years,
                       seasons = seasons
)

surveyStart = c(0,0.17,0) #c(0.75,0)
surveyEnd =  c(.5,.33,0.17) # c(1,0) Does the survey last throughout the season it's conducted?
surveySeason = c(2,1,1) # c(2,1)Which seasons do the surveys occur in


# Catch 
ctot <- read.table(file.path(wd,'ton.in'))# %>% mutate(years = rep(years, each = nseason),
                                                  # season = rep(seasons, nyear))  
# Proportions 
props <- read.table(file.path(wd,'prop.in'))# %>% mutate(years = rep(years, each = nseason),
                                                    # season = rep(seasons, nyear))  
props[is.na(props)] <- -1

weca <- read.table(file.path(wd,'weca.in'))
mat <- read.table(file.path(wd,'propmat.in'))
M2 <- read.table(file.path(wd,'natmor.in'))
nsamples <- read.table(file.path(wd,'nsamples.in'))
# Order stuff correctly 
weca_m <- array(
  data = t(as.matrix(weca)),   # transpose so ages become the first dimension
  dim  = c(nage, nyear+1, nseason)         # 4 ages, 51 years, 2 seasons
)
mat_m <- array(
  data = t(as.matrix(mat)),   # transpose so ages become the first dimension
  dim  = c(nage, nyear, nseason)         # 4 ages, 51 years, 2 seasons
)

mat_M2 <- array(
  data = t(as.matrix(M2)),   # transpose so ages become the first dimension
  dim  = c(nage, nyear, nseason)         # 4 ages, 51 years, 2 seasons
)

props <- array(
  data = t(as.matrix(props)),   # transpose so ages become the first dimension
  dim  = c(nage, nyear, nseason)         # 4 ages, 51 years, 2 seasons
)

nsamples <- matrix(nsamples[,1], nrow = nyear, ncol = nseason)
mtrx<- list("M" = mat_M2,
            "mat" = mat_m,
            "west" = weca_m,
            "weca" = weca_m)

Surveyobs <- survey_to_matrix(survey.in, years)
Surveyobs[4,,3] <- -1#Surveyobs[3,years %in% c(2024),2]*2

Qminage = c(0,1,1) # Qminage = c(0,1) minimum age in surveys
Qmaxage = c(3,3,2) #Qmaxage = c(1,3)

surveyCV =  list(c(0,1,2),
                 c(1,2,3),
                 c(1)
) #c(1,2)),
# Now modify catch and age comps 
Catchobs <- array(ctot[,1], dim = c(1,nyear,nseason))

# now into get_tmb_parameters 
nocatch <- matrix(1, ncol = nseason, nrow = nyear)
nocatch[length(nocatch)] <- 0


df.tmb <- get_TMB_parameters(
  mtrx = mtrx, # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  CatchProportions = props,
  nsamples = nsamples,
  years = years, # Years to run
  #endYear = 2023,
  leavesurveyout = c(1,1,1),
  nseason = nseason, # Number of seasons
  ages = ages, # Ages of the species
  recseason = 1, # Season where recruitment occurs
  CminageSeason = c(0,0),
  Fmaxage = 2, # Fully selected fishing mortality age
  Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
  Qmaxage = Qmaxage, #Qmaxage = c(1,3)
  minSDcatch = .1,
  minSDsurvey = .2,
  #startYear = 1980,
  blocks = c(1974,2015),
  maxSDcatch = 10,
  Fbarage = c(1,2),
  endFseason = 1, # which season does fishing stop in the final year of data
  nocatch = as.matrix(nocatch),
  surveyStart = surveyStart, #c(0.75,0)
  surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
  surveySD =  surveyCV, #c(1,2)),
  catchSD = list(c(0),
                 c(0)),
  #csd_break = c(2006),
  estSD = c(0,0,0), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,1) # Factor for relative strength of log-likelihood
  
)

df.tmb$nsamples[df.tmb$nsamples == 0] <- 1


parms <- getParms(df.tmb)
# parms$logsdc <- log(.2)
# df.tmb$tuneStart <- which(years %in% 2006)-1
mps <-getMPS(df.tmb, parms, mapExtra ='logSDcatch')# Set boundaries

sas <- runAssessment(df.tmb, parms = parms,mps = mps, silent = TRUE, debug = TRUE)

sas$reps





