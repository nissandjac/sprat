## Forecast model ##

library(smsR)

wd <- "C:/Users/nsja/Dropbox/DTU/BEBRIS/two_seasons/"
#devtools::load_all()
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

surveyStart = c(0,0.17,0) #c(0.75,0)
surveyEnd =  c(.5,.33,0.17) # c(1,0) Does the survey last throughout the season it's conducted?
surveySeason = c(2,1,1) # c(2,1)Which seasons do the surveys occur in

# Load packages and files #


ages <- 0:maxage
beta <- 90000
Surveyobs <- survey_to_matrix(dat[['survey']], years)
Surveyobs[4,,3] <- -1

Qminage = c(0,1,1) # Qminage = c(0,1) minimum age in surveys
Qmaxage = c(3,3,2) #Qmaxage = c(1,3)

surveyCV =  list(c(0,1,2),
                 c(1,2),
                 c(1)
) #c(1,2)),
Catchobs <- df_to_matrix(dat[['canum']], season =  1:2)


nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#')#, skip = 3)


df.tmb <- get_TMB_parameters(
  mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  years = years, # Years to run
  #endYear = 2023,
  leavesurveyout = c(1,1,1),
  nseason = nseason, # Number of seasons
  ages = ages, # Ages of the species
  recseason = 1, # Season where recruitment occurs
  CminageSeason = c(0,0),
  Fmaxage = 3, # Fully selected fishing mortality age
  Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
  Qmaxage = Qmaxage, #Qmaxage = c(1,3)
  minSDcatch = .1,
  minSDsurvey = .3,
  startYear = 1980,
  blocks = c(1974,2015),
  maxSDcatch = 10,
  Fbarage = c(1,2),
  tuneCatch = 1,
  endFseason = 1, # which season does fishing stop in the final year of data
  nocatch = as.matrix(nocatch),
  surveyStart = surveyStart, #c(0.75,0)
  surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
  surveySD =  surveyCV, #c(1,2)),
  catchSD = list(c(0,1,2),
                 c(0,1,2)),
  csd_break = c(2006),
  estSD = c(0,0,0), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,1) # Factor for relative strength of log-likelihood

)

parms <- getParms(df.tmb)
mps <-getMPS(df.tmb, parms)# Set boundaries
mps$logbeta <- NULL
sas <- runAssessment(df.tmb, parms = parms,mps = mps, silent = TRUE, debug = TRUE)

sas$reps
plot(sas)

#mr <- mohns_rho(df.tmb, parms , peels = 5, mps =mps)
p2 <- plotDiagnostics(df.tmb,sas)
p2$SR

library(tidyverse)
yobs <- getYield(df.tmb)
yest <- getCatch(df.tmb, sas)
yest$catchobs <- yobs$Yield
ggplot(yest %>% filter(years > 1990), aes(x = years, y = (catchobs-Catch)/catchobs))+geom_line()+theme_classic()+
  geom_col()

ggplot(yest , aes(x = years, y = Catch))+geom_line()+theme_classic()+
  geom_point(aes(y = catchobs))+
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'red', alpha = .2)

plotBubbles(sas)

p2 <- plotDiagnostics(df.tmb, sas)
p2$survey

# Now one where the slope is estimated


parms <- getParms(df.tmb)
mps <-getMPS(df.tmb, parms)# Set boundaries
mps$logbeta <- NULL
sas_beta <- runAssessment(df.tmb, parms = parms,mps = mps, silent = TRUE, debug = TRUE)

sas$reps
plot(sas)



x <- getR(df.tmb, sas)
x2 <- getR(sas_beta$dat, sas_beta)

plot(x$R)
lines(x2$R)


library(patchwork)

plotBubbles(sas)/plotBubbles(sas_beta)
