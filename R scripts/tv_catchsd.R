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


# Read the survey data
# sd_q3 <- read.csv('data/spratindices-nohauldur/SDQ3.csv', sep = '')
# sd_q3$years <- as.numeric(rownames(sd_q3))
# Fishing eff
#dat$survey$eff[dat$survey$Survey == "IBTS_Q1_0"] <- 1e-5
# Normalize effort to 1
# Format input data matrices into TMB style
# mtrx <- sumtable_to_matrix(sandeel.age)
Surveyobs <- survey_to_matrix(dat[['survey']], years)
Surveyobs[1, years %in% c(2022,2023),1] <- median(Surveyobs[1, ,1])
#Surveyobs[1,years %in% 2024,1] <- -1
# Rescale the recruitment stuff
# Adjusted this to get it on the same scale as the other Qs
# Surveyobs[1,,1] <- Surveyobs[1,,1] * 0.001
# Surveyobs[Surveyobs < 0] <- -99

Catchobs <- df_to_matrix(dat[['canum']], season =  1:2)
nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#')#, skip = 3)
#Feffort[which(df.tmb$years == 2013),2] <- 0 # This is different two places in the input files
#nocatch[50,] <- 0
# Save some data for package example

# Calculate the catch yield per season
yield <- Catchobs * dat[['mtrx']]$weca[,1:nyear,]
yield_sum <- apply(yield, c(2, 3), sum)
rel_contrib <- prop.table(yield_sum, margin = 1)
names(rel_contrib) <- paste('season',1:nseason)

dat$nocatch <- as.matrix(nocatch)#*0+1
# dat$nocatch[rel_contrib < .2] <- 0
# #dat$nocatch[length(dat$nocatch)] <- 0
# dat$nocatch[20:nyear,2] <- 0
#dat$nocatch[dat$effort == 0] <- 0
powerIN <- list(NA, NA, NA)

#Surveyobs[,years %in% 2023,3] <- -1


df.tmb <- get_TMB_parameters(
  mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  years = years, # Years to run
  #startYear = 1997,
  #endYear = 2023,
  nseason = nseason, # Number of seasons
  ages = ages, # Ages of the species
  recseason = 1, # Season where recruitment occurs
  CminageSeason = c(0,0),
  Fmaxage = 2, # Fully selected fishing mortality age
  Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
  Qmaxage = Qmaxage, #Qmaxage = c(1,3)
  leavesurveyout = c(1,1,1),
  minSDcatch = .1,
  powers = powerIN,
  maxSDcatch =10,
  randomF = 1,
  blocks = c(1974,2015),
  endFseason = 1, # which season does fishing stop in the final year of data
  nocatch = as.matrix(dat$nocatch),
  surveyStart = surveyStart, #c(0.75,0)
  surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
  surveySD =  surveyCV, #c(1,2)),
  catchSD = list(c(0,1,2),
                 c(0,1,2)),
  #csd_break = c(2000),
  estSD = c(0,0,0), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,1) # Factor for relative strength of log-likelihood

)

# df.tmb$Cidx_CV[,2] <- df.tmb$Cidx_CV[,1]
# df.tmb$cscalar[1,] <- 0.2


parms <- getParms(df.tmb )
# parms$logSDcatch <- parms$logSDcatch*0+log(4)
# parms$logSDcatch[1:3] <- log(0.2)
mps <-getMPS(df.tmb, parms)#, mapExtra = 'logSDcatch')# Set boundaries

# mps$logSDcatch <- factor(c(0,1,2,NA, NA,NA))
# mps$logbeta <- NULL
# This model works best if SDsurvey is mapped
#lbound <- parms$logSDcatch*0+log(0.5)

sas <- runAssessment(df.tmb, parms = parms,mps = mps, silent = TRUE,
                     upr = list('logSDCatch'  = log(4)))
sas$reps
plot(sas)
plotBubbles(sas)
p2 <- plotDiagnostics(df.tmb, sas)
p2$catch
p2$survey

mr <- mohns_rho(df.tmb, parms, peels = 5)

# Compare yield 
library(tidyverse)

yobs <- getYield(df.tmb)
yest <- getCatch(df.tmb, sas)

yest$catchobs <- yobs$Yield

ggplot(yest, aes(x = years, y = Catch))+geom_line()+theme_classic()+
  geom_point(aes(y = catchobs))+
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'red', alpha = .2)

# And the catch per season 
library(tidyverse)
yseason <- getYieldSeason(df.tmb) %>% pivot_longer(paste('season',1:nseason, sep='-'))
yest_s <- getCatchSeason(df.tmb,sas)
yest_s$Catchobs <- yseason$value

ggplot(yest_s, aes(x = years, y = Catchseason))+geom_line()+
  geom_point(aes(y = Catchobs))+theme_classic()+
  facet_wrap(~season)
# See how wrong we are 
ggplot(yest, aes(x='1', y= (log(catchobs)-log(Catch))/log(catchobs)))+geom_violin()+
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "red")+
  theme_classic()



p2 <- plotDiagnostics(df.tmb, sas)
p2$catch

psurv <- removeSurvey(sas)
#saveRDS(psurv, file.path(wd, 'leavesurvey_2010.RDS'))
ggsave(file.path(wd,'leavesurvey_out_change.png'),psurv, width= 16, height = 10)

saveRDS(sas, file = file.path(wd,'two_seasons_csd.RDS'))
write.table(mr$df.save, file = file.path(wd,'mohns_table_csd.csv'), row.names = FALSE)
write.table(mr$mohns, file = file.path(wd,'mohns_table_result_csv.csv'), row.names = FALSE)


library(smsR)

wd <- "C:/Users/nsja/Dropbox/DTU/BEBRIS/four_quarters/"

maxage <- 3
years = 1974:2024
nyear <- length(years)
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

powerIN <- list(NA, NA, NA)
# Load packages and files #


ages <- 0:maxage
beta <- 90000

# Fishing eff
#dat$survey$eff[dat$survey$Survey == "IBTS_Q1_0"] <- 1e-5
# Normalize effort to 1
# Format input data matrices into TMB style
# mtrx <- sumtable_to_matrix(sandeel.age)
Surveyobs <- survey_to_matrix(dat[['survey']], years)
# Rescale the recruitment stuff
# Adjusted this to get it on the same scale as the other Qs
# Surveyobs[1,,1] <- Surveyobs[1,,1] * 0.001
# Surveyobs[Surveyobs < 0] <- -99

Catchobs <- df_to_matrix(dat[['canum']], season =  1:4)

nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#')#, skip = 3)
#Feffort[which(df.tmb$years == 2013),2] <- 0 # This is different two places in the input files
yield <- Catchobs * dat[['mtrx']]$weca[,1:nyear,]
yield_sum <- apply(yield, c(2, 3), sum)
rel_contrib <- as.data.frame(prop.table(yield_sum, margin = 1)) 
names(rel_contrib) <- paste('age',0:3)



nocatch <- as.matrix(nocatch)#*0+1
# nocatch[rel_contrib < .01] <- 0
# Save some data for package example
dat$nocatch <- as.matrix(nocatch)#*0+1
#dat$nocatch[1,4] <- 0
#dat$nocatch[dat$effort == 0] <- 0
powerIN <- list(NA, NA, NA)


df.tmb <- get_TMB_parameters(
  mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  years = years, # Years to run
  nseason = nseason, # Number of seasons
  ages = ages, # Ages of the species
  recseason = 1, # Season where recruitment occurs
  CminageSeason = c(0,0,0,0),
  Fmaxage = 2, # Fully selected fishing mortality age
  Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
  Qmaxage = Qmaxage, #Qmaxage = c(1,3)
  blocks = c(1974,2015),
  #maxSDcatch = sqrt(10),
  Fbarage = c(1,2),
  isFseason = c(1,1,1,0), # Seasons to calculate fishing in
  powers = powerIN,
  endFseason = 2,  # which season does fishing stop in the final year of data
  nocatch = as.matrix(dat$nocatch),
  surveyStart = surveyStart, #c(0.75,0)
  surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
  surveySD =  surveyCV, #c(1,2)),
  #randomF = 1,
  csd_break = 1995,
  catchSD = list(c(0,1,2),
                 c(0,1,2),
                 c(0,1,2),
                 c(0,1,2)),
  estSD = c(0,0,0), # Estimate,
  recmodel = 1,
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,1) # Factor for relative strength of log-likelihood
  
)
# Get initial parameter structure
parms <- getParms(df.tmb )
#parms$logSDcatch <- log(0*parms$logSDcatch + c(rep(0.2,8), rep(0.8,8)))
mps <-getMPS(df.tmb, parms)
sas <- runAssessment(df.tmb, parms = parms,mps = mps, 
                     silent = TRUE, 
                     debug = TRUE
)

plot(sas)
plotBubbles(sas)
mr <- mohns_rho(df.tmb, parms, peels = 5)

yobs <- getYield(df.tmb)
yest <- getCatch(df.tmb, sas)

yest$catchobs <- yobs$Yield

ggplot(yest, aes(x = years, y = Catch))+geom_line()+theme_classic()+
  geom_point(aes(y = catchobs))+
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'red', alpha = .2)

psurv <- removeSurvey(sas)
ggsave(file.path(wd,'leavesurvey_out_change.png'),psurv, width= 16, height = 10)

saveRDS(sas, file = file.path(wd,'four_seasons_csd.RDS'))
write.table(mr$df.save, file = file.path(wd,'mohns_table_csd.csv'), row.names = FALSE)
write.table(mr$mohns, file = file.path(wd,'mohns_table_result_csv.csv'), row.names = FALSE)

yseason <- getYieldSeason(df.tmb) %>% pivot_longer(paste('season',1:nseason, sep='-'))
yest_s <- getCatchSeason(df.tmb,sas)
yest_s$Catchobs <- yseason$value

ggplot(yest_s, aes(x = years, y = Catchseason))+geom_line()+
  geom_point(aes(y = Catchobs))+theme_classic()+
  facet_wrap(~season, scales = 'free_y')
# See how wrong we are 
ggplot(yest, aes(x='1', y= (log(catchobs)-log(Catch))/log(catchobs)))+geom_violin()+
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "red")+
  theme_classic()+scale_y_continuous('relative error\nlog')+
  xlab(label = '')

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

# Read the survey data
# sd_q3 <- read.csv('data/spratindices-nohauldur/SDQ3.csv', sep = '')
# sd_q3$years <- as.numeric(rownames(sd_q3))
# Fishing eff
#dat$survey$eff[dat$survey$Survey == "IBTS_Q1_0"] <- 1e-5
# Normalize effort to 1
# Format input data matrices into TMB style
# mtrx <- sumtable_to_matrix(sandeel.age)
Surveyobs <- survey_to_matrix(dat[['survey']], years)
# Rescale the recruitment stuff
# Adjusted this to get it on the same scale as the other Qs
# Surveyobs[1,,1] <- Surveyobs[1,,1] * 0.001
# Surveyobs[Surveyobs < 0] <- -99

Catchobs <- df_to_matrix(dat[['canum']], season =  1:2)
nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#')#, skip = 3)
#Feffort[which(df.tmb$years == 2013),2] <- 0 # This is different two places in the input files
#nocatch[50,] <- 0
# Save some data for package example

# Calculate the catch yield per season
yield <- Catchobs * dat[['mtrx']]$weca[,1:nyear,]
yield_sum <- apply(yield, c(2, 3), sum)
rel_contrib <- prop.table(yield_sum, margin = 1)
names(rel_contrib) <- paste('season',1:nseason)

dat$nocatch <- as.matrix(nocatch)#*0+1
# dat$nocatch[rel_contrib < .2] <- 0
# #dat$nocatch[length(dat$nocatch)] <- 0
# dat$nocatch[20:nyear,2] <- 0
#dat$nocatch[dat$effort == 0] <- 0
powerIN <- list(NA, NA, NA)

#Surveyobs[,years %in% 2023,3] <- -1


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
  minSDcatch = .1,
  # maxSDcatch = .6,
  # randomF = 1,
  blocks = c(1974,2015),
  endFseason = 1, # which season does fishing stop in the final year of data
  nocatch = as.matrix(dat$nocatch),
  surveyStart = surveyStart, #c(0.75,0)
  surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
  surveySD =  surveyCV, #c(1,2)),
  catchSD = list(c(0,1,2),
                 c(0,1,2)),
  csd_break = c(1995),
  estSD = c(0,0,0), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,1) # Factor for relative strength of log-likelihood
  
)

parms <- getParms(df.tmb )
parms$logSDcatch <- parms$logSDcatch*0 + log(0.5)
parms$logSDcatch[,2] <- parms$logSDcatch[2:3,1]*0 + log(0.2)

mps <-getMPS(df.tmb, parms, mapExtra = 'logSDcatch')# Set boundaries
#mps$logbeta <- NULL


sas <- runAssessment(df.tmb, parms = parms,mps = mps, silent = TRUE)
sas$reps
plot(sas)
plotBubbles(sas)

mr <- mohns_rho(df.tmb, parms, peels = 5)

# Compare yield 
library(tidyverse)

yobs <- getYield(df.tmb)
yest <- getCatch(df.tmb, sas)

yest$catchobs <- yobs$Yield

ggplot(yest, aes(x = years, y = Catch))+geom_line()+theme_classic()+
  geom_point(aes(y = catchobs))+
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'red', alpha = .2)

# And the catch per season 
library(tidyverse)
yseason <- getYieldSeason(df.tmb) %>% pivot_longer(paste('season',1:nseason, sep='-'))
yest_s <- getCatchSeason(df.tmb,sas)
yest_s$Catchobs <- yseason$value

ggplot(yest_s, aes(x = years, y = Catchseason))+geom_line()+
  geom_point(aes(y = Catchobs))+theme_classic()+
  facet_wrap(~season)
# See how wrong we are 
ggplot(yest, aes(x='1', y= (log(catchobs)-log(Catch))/log(catchobs)))+geom_violin()+
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "red")+
  theme_classic()



p2 <- plotDiagnostics(df.tmb, sas)
p2$catch

psurv <- removeSurvey(sas)
#saveRDS(psurv, file.path(wd, 'leavesurvey_2010.RDS'))
ggsave(file.path(wd,'leavesurvey_out_change.png'),psurv, width= 16, height = 10)

saveRDS(sas, file = file.path(wd,'two_seasons_fixed.RDS'))
write.table(mr$df.save, file = file.path(wd,'mohns_table_fixed.csv'), row.names = FALSE)
write.table(mr$mohns, file = file.path(wd,'mohns_table_result_fixed.csv'), row.names = FALSE)






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

# Read the survey data
# sd_q3 <- read.csv('data/spratindices-nohauldur/SDQ3.csv', sep = '')
# sd_q3$years <- as.numeric(rownames(sd_q3))
# Fishing eff
#dat$survey$eff[dat$survey$Survey == "IBTS_Q1_0"] <- 1e-5
# Normalize effort to 1
# Format input data matrices into TMB style
# mtrx <- sumtable_to_matrix(sandeel.age)
Surveyobs <- survey_to_matrix(dat[['survey']], years)
# Rescale the recruitment stuff
# Adjusted this to get it on the same scale as the other Qs
# Surveyobs[1,,1] <- Surveyobs[1,,1] * 0.001
# Surveyobs[Surveyobs < 0] <- -99

Catchobs <- df_to_matrix(dat[['canum']], season =  1:2)
Catchobs[,years > 2015,2] <- 0
nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#')#, skip = 3)
#Feffort[which(df.tmb$years == 2013),2] <- 0 # This is different two places in the input files
#nocatch[years>2021,] <- 0
# Save some data for package example

# Calculate the catch yield per season
yield <- Catchobs * dat[['mtrx']]$weca[,1:nyear,]
yield_sum <- apply(yield, c(2, 3), sum)
rel_contrib <- prop.table(yield_sum, margin = 1)
names(rel_contrib) <- paste('season',1:nseason)

dat$nocatch <- as.matrix(nocatch)#*0+1
# dat$nocatch[rel_contrib < .2] <- 0
# #dat$nocatch[length(dat$nocatch)] <- 0
# dat$nocatch[20:nyear,2] <- 0
#dat$nocatch[dat$effort == 0] <- 0
powerIN <- list(NA, NA, NA)

#Surveyobs[,years %in% 2023,3] <- -1


df.tmb <- get_TMB_parameters(
  mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  years = years, # Years to run
  endYear = 2023,
  nseason = nseason, # Number of seasons
  ages = ages, # Ages of the species
  recseason = 1, # Season where recruitment occurs
  CminageSeason = c(0,0),
  Fmaxage = 2, # Fully selected fishing mortality age
  Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
  Qmaxage = Qmaxage, #Qmaxage = c(1,3)
  minSDcatch = .1,
  leavesurveyout = c(1,1,1),
  # maxSDcatch = .6,
  # randomF = 1,
  blocks = c(1974,2015),
  endFseason = 1, # which season does fishing stop in the final year of data
  nocatch = as.matrix(dat$nocatch),
  surveyStart = surveyStart, #c(0.75,0)
  surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
  surveySD =  surveyCV, #c(1,2)),
  catchSD = list(c(0,1,2),
                 c(0,1,2)),
  csd_break = c(1995),
  estSD = c(0,0,0), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,1) # Factor for relative strength of log-likelihood
  
)

parms <- getParms(df.tmb )
#parms$logSDcatch <- parms$logSDcatch*0 + log(0.01)
#parms$logSDcatch[,2] <- parms$logSDcatch[2:3,1]*0 + log(0.01)

mps <-getMPS(df.tmb, parms)#, mapExtra = 'logSDcatch')# Set boundaries
#mps$logbeta <- NULL


sas <- runAssessment(df.tmb, parms = parms,mps = mps, silent = TRUE)
sas$reps
plot(sas)
plotBubbles(sas, CVscale = FALSE)

mr <- mohns_rho(df.tmb, parms, peels = 5)


p2 <- plotDiagnostics(df.tmb, sas)
p2$catch


removeSurvey(sas)
# Compare yield 
library(tidyverse)

yobs <- getYield(df.tmb)
yest <- getCatch(df.tmb, sas)

yest$catchobs <- yobs$Yield

ggplot(yest, aes(x = years, y = Catch))+geom_line()+theme_classic()+
  geom_point(aes(y = catchobs))+
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'red', alpha = .2)

# And the catch per season 
library(tidyverse)
yseason <- getYieldSeason(df.tmb) %>% pivot_longer(paste('season',1:nseason, sep='-'))
yest_s <- getCatchSeason(df.tmb,sas)
yest_s$Catchobs <- yseason$value

ggplot(yest_s, aes(x = years, y = Catchseason))+geom_line()+
  geom_point(aes(y = Catchobs))+theme_classic()+
  facet_wrap(~season)
# See how wrong we are 
ggplot(yest, aes(x='1', y= (log(catchobs)-log(Catch))/log(catchobs)))+geom_violin()+
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "red")+
  theme_classic()



p2 <- plotDiagnostics(df.tmb, sas)
p2$catch

# psurv <- removeSurvey(sas)
# #saveRDS(psurv, file.path(wd, 'leavesurvey_2010.RDS'))
# ggsave(file.path(wd,'leavesurvey_out_change.png'),psurv, width= 16, height = 10)

saveRDS(sas, file = file.path(wd,'two_seasons_sdc2015.RDS'))
write.table(mr$df.save, file = file.path(wd,'mohns_table_sdc2015.csv'), row.names = FALSE)
write.table(mr$mohns, file = file.path(wd,'mohns_table_result_sdc2015.csv'), row.names = FALSE)
