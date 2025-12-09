library(smsR)

wd <- "C:/Users/nsja/Dropbox/DTU/BEBRIS/Final model/"
#devtools::load_all()
# M is not correct here
maxage <- 3
years = 1980:2024
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
# Load packages and files #


ages <- 0:maxage
beta <- 90000
Surveyobs <- survey_to_matrix(dat[['survey']], years)
#Surveyobs[4,,2] <- -1#Surveyobs[2,years %in% c(2023),2]*0.5
Surveyobs[4,,3] <- -1#Surveyobs[3,years %in% c(2024),2]*2

Qminage = c(0,1,1) # Qminage = c(0,1) minimum age in surveys
Qmaxage = c(3,3,2) #Qmaxage = c(3,3,2)

surveyCV =  list(c(0,1,2),
                 c(1,2,3),
                 c(1)
                 ) #c(1,2)),


library(tidyverse)


xx <- data.frame(age1 = Surveyobs[2,1:(nyear-1),2],
                 age2 = Surveyobs[3,2:nyear,2],
                 years = years[2:nyear])

ggplot(xx %>% filter(years > 1991), aes(x = age1, y = age2))+geom_point()+
  geom_text(aes(label = years), vjust = -0.5)
Catchobs <- df_to_matrix(dat[['canum']], season =  1:2)
#Catchobs[]

nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#')#, skip = 3)
powerIN <- c(NA, NA, NA)


df.tmb <- get_TMB_parameters(
  mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  years = years, # Years to run
 # endYear = 2022,
  leavesurveyout = c(1,1,1),
  nseason = nseason, # Number of seasons
  ages = ages, # Ages of the species
  recseason = 1, # Season where recruitment occurs
  CminageSeason = c(0,0),
  Fmaxage = 3, # Fully selected fishing mortality age
  Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
  Qmaxage = Qmaxage, #Qmaxage = c(1,3)
  minSDcatch = .2,
  maxSDcatch = 10,
  minSDsurvey = .3,
  #penepsC = 1e-3,
  blocks = c(1980,2015),
  Fbarage = c(1,2),
  tuneCatch = 1,
  tuneStart = 2012, # Put this into the last selectivity season
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
sas <- runAssessment(df.tmb, parms = parms,mps = mps, silent = TRUE, debug = TRUE)

sas$reps
plot(sas)

mr <- mohns_rho(df.tmb, parms , peels = 5, mps =mps)

plotBubbles(sas)
p2 <- plotDiagnostics(df.tmb, sas)
p2$catch
p2$survey

library(tidyverse)
yobs <- getYield(df.tmb)
yest <- getCatch(df.tmb, sas)

yest$catchobs <- yobs$Yield

ggplot(yest , aes(x = years, y = Catch))+geom_line()+theme_classic()+
  geom_point(aes(y = catchobs))+
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'red', alpha = .2)


ggplot(yest %>% filter(years > 1990), aes(x = years, y = (catchobs-Catch)/catchobs))+geom_line()+theme_classic()+
  geom_col()
  #geom_point(aes(y = catchobs))+


# Plot Catchobs weight per season 
# Calculate the catch yield per season
yield <- Catchobs * dat[['mtrx']]$weca[,1:nyear,]
dimnames(yield) <- list(ages = df.tmb$age,
                        years = df.tmb$years,
                        seasons = 1:df.tmb$nseason)
ydf <- as.data.frame.table(yield, responseName = "value") %>% mutate(years = as.numeric(years))

ggplot(ydf, aes(x = years, y = value, color = ages))+facet_wrap(~seasons, nrow = 2, scales = 'free_y')+
  geom_line()+geom_point()+theme_classic()

