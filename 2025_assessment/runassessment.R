library(smsR)
devtools::load_all()
wd <- "C:/Users/nsja/Dropbox/DTU/BEBRIS/2025_assessment/"

maxage <- 3
years = 1974:2024
nyear <- length(years)
seasons <- 1:4
nseason <- length(seasons)

  dat <- getDataSMS(wd,
                        maxage = maxage,
                        survey.age = list(0, 0:3, 1:3, 1:3), # Ages in the two surveys
                        survey.years = list(1983:2023 ,1982:2021, 1992:2021, 2006:2021),
                        survey.names = c('age0_ITBS','IBTS_Q1','IBTS_Q3', 'HERAS'),
                        survey.quarter = c(1,1,3, 2),
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
  
  powerIN <- list(0, NA, NA)
# Load packages and files #


ages <- 0:maxage
beta <- 90000

# Fishing eff

# Normalize effort to 1
# Format input data matrices into TMB style
# mtrx <- sumtable_to_matrix(sandeel.age)
Surveyobs <- survey_to_matrix(dat[['survey']], years)
Catchobs <- df_to_matrix(dat[['canum']], season =  1:4)

nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#', skip = 3)
#Feffort[which(df.tmb$years == 2013),2] <- 0 # This is different two places in the input files

# Save some data for package example
dat$nocatch <- nocatch#*0+1
#dat$nocatch[dat$effort == 0] <- 0



df.tmb <- get_TMB_parameters(
  mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  years = years, # Years to run
  nseason = nseason, # Number of seasons
  useEffort = 0,
  ages = ages, # Ages of the species
  recseason = 1, # Season where recruitment occurs
  CminageSeason = c(0,0,0,0),
  Fmaxage = 2, # Fully selected fishing mortality age
  Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
  Qmaxage = Qmaxage, #Qmaxage = c(1,3)
  minSDcatch = 0.1,
  Fbarage = c(1,2),
  isFseason = c(1,1,1,0), # Seasons to calculate fishing in
  powers = list(0,NA,NA),
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
  recmodel = 2, # Chose recruitment model (2 = estimated)
  estSD = c(0,0,2), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,0.1) # Factor for relative strength of log-likelihood

)



# Get initial parameter structure
parms <- getParms(df.tmb )

# Get non-estimated parameters, based on info in df.tmb
mps <- getMPS(df.tmb, parms)

# Set boundaries
# This model works best if SDsurvey is mapped
#parms[names(parms) == 'SDsurvey'][[1]] <- c(0.5, 0.3, 0.3,.3)
parms <- getSpratparms(wd, df.tmb, parms)


#parms$Fyear <- c(parms$Fyear,1)
# Try initial parameters in weird spots
sas <- runAssessment(df.tmb, parms = parms, 
                     mps = mps, silent = FALSE)
#
sas$reps

p1 <- smsPlots(df.tmb = df.tmb,sas)

p2 <- plotDiagnostics(df.tmb, sas)

#mr <- mohns_rho(df.tmb, peels = 5, parms, mps, plotfigure = TRUE)


# Compare all the results with ADMB

source(file.path(wd,"compare_TMB_admb.R"))

library(tidyverse)
compareTMB_ADMB(wd, sas, df.tmb)


# Compare everything with the correctly parameterized one 

## SSB 

sms <- read.table(file.path(wd, 'summary_table_raw.out'), header = TRUE)

plot(sms$Year, sms$SSB)
lines(c(df.tmb$years, max(df.tmb$years)+1), sas$x$SSB)

# Survey and catch residuals 

sresid <- read.table(file.path(wd, 'survey_residuals.out'), sep = ',', fill = TRUE)

Sresid.tmb <- sas$x$resid_survey


survey1 <- sresid %>% filter(V2 == 1) %>% select(-V1) %>% 
  rename('survey' = V2,
         'year' = V3,
         '0' = V4,
         '1' = V5,
         '2' = V6,
         '3' = V7) %>% 
  pivot_longer(3:6, names_to = 'age', values_to = 'residual') %>% mutate(model = 'admb')


survey1.tmb <- data.frame(survey = 1,year = df.tmb$years,
                          '0' = Sresid.tmb[1,,1],
                          '1' = Sresid.tmb[2,,1],
                          '2' = Sresid.tmb[3,,1],
                          '3' = Sresid.tmb[4,,1]) %>%  
  pivot_longer(3:6, names_to = 'age', values_to = 'residual')  %>% 
      mutate(age = str_replace(age,'X0', '0'),
         age = str_replace(age,'X1', '1'),
         age = str_replace(age,'X2', '2'),
         age = str_replace(age,'X3', '3')) %>% 
      mutate(model = 'tmb') %>% filter(residual > -98)

df.plot <- rbind(survey1, survey1.tmb)

ggplot(survey1,aes( x = year, y = residual, color = model))+
  geom_line()+
  geom_point(data = survey1.tmb, alpha =.5)+
  facet_wrap(~factor(age))


# Survey 2 

survey2 <- sresid %>% filter(V2 == 2) %>% select(-c(V1,V7)) %>% 
  rename('survey' = V2,
         'year' = V3,
         '1' = V4,
         '2' = V5,
         '3' = V6) %>% 
  pivot_longer(3:5, names_to = 'age', values_to = 'residual') %>% mutate(model = 'admb') 

survey2.tmb <- data.frame(survey = 2,year = df.tmb$years,
                          '1' = Sresid.tmb[2,,2],
                          '2' = Sresid.tmb[3,,2],
                          '3' = Sresid.tmb[4,,2]) %>%  
  pivot_longer(3:5, names_to = 'age', values_to = 'residual')  %>% 
  mutate(age = str_replace(age,'X1', '1'),
         age = str_replace(age,'X2', '2'),
         age = str_replace(age,'X3', '3')) %>% 
  mutate(model = 'tmb') %>% filter(residual > -98)

df.plot <- rbind(survey2, survey2.tmb)

ggplot(survey2 ,aes( x = year, y = residual, color = model))+geom_line()+
  geom_point(data = survey2.tmb, alpha = 0.5)+
  facet_wrap(~factor(age))


# Survey 2 

survey3 <- sresid %>% filter(V2 == 3) %>% select(-c(V1,V7)) %>% 
  rename('survey' = V2,
         'year' = V3,
         '1' = V4,
         '2' = V5,
         '3' = V6) %>% 
  pivot_longer(3:5, names_to = 'age', values_to = 'residual') %>% mutate(model = 'admb') 

survey3.tmb <- data.frame(survey = 3,year = df.tmb$years,
                          '1' = Sresid.tmb[2,,3],
                          '2' = Sresid.tmb[3,,3],
                          '3' = Sresid.tmb[4,,3]) %>%  
  pivot_longer(3:5, names_to = 'age', values_to = 'residual')  %>% 
  mutate(age = str_replace(age,'X1', '1'),
         age = str_replace(age,'X2', '2'),
         age = str_replace(age,'X3', '3')) %>% 
  mutate(model = 'tmb') %>% filter(residual > -98)

df.plot <- rbind(survey3, survey3.tmb)

ggplot(survey3 ,aes( x = year, y = residual, color = model))+geom_line()+
  geom_point(data = survey3.tmb, alpha = 0.5)+
  facet_wrap(~factor(age))




# Okay check catch residuals and solve survey after 


cresid <- read.table(file.path(wd, 'catch_residuals.out'), sep = ',', fill = TRUE)

Cresid.tmb <- sas$x$resid_catch


catch1 <- cresid  %>% select(-V1) %>% 
  rename('season' = V2,
         'year' = V3,
         '0' = V4,
         '1' = V5,
         '2' = V6,
         '3' = V7) %>% 
  pivot_longer(3:6, names_to = 'age', values_to = 'residual') %>% mutate(model = 'admb')
catch1$residual[catch1$residual < -50] <- NA


catch1.tmb <- data.frame(season = rep(c(1,2,3,4), each = df.tmb$nyears),
                          year = df.tmb$years,
                          '0' = as.numeric(Cresid.tmb[1,,]),
                          '1' = as.numeric(Cresid.tmb[2,,]),
                          '2' = as.numeric(Cresid.tmb[3,,]),
                          '3' = as.numeric(Cresid.tmb[4,,])) %>%  
  pivot_longer(3:6, names_to = 'age', values_to = 'residual')  %>% 
  mutate(age = str_replace(age,'X0', '0'),
         age = str_replace(age,'X1', '1'),
         age = str_replace(age,'X2', '2'),
         age = str_replace(age,'X3', '3')) %>% 
  mutate(model = 'tmb')
catch1.tmb$residual[catch1.tmb$residual  < -50] <- NA

df.plot <- rbind(catch1, catch1.tmb)

ggplot(catch1, aes( x = year, y = residual, color = model))+geom_line()+
  geom_point(data = catch1.tmb)+
  facet_grid(season~factor(age))








