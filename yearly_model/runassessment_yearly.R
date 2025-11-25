library(smsR)

wd <- "C:/Users/nsja/Dropbox/DTU/BEBRIS/yearly_model/"

# M is not correct here 


maxage <- 3
years = 1974:2024
nyear <- length(years)
seasons <- 1
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
  surveyStart = c(0,0.5,0.6) # Survey starts in jan, juli, august 
  surveyEnd =  c(.2,.58,0.7) # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = c(1,1,1) # c(2,1)Which seasons do the surveys occur in
  surveyCV =  list(c(0,1),
                   c(1,2),
                   c(1,2)) #c(1,2)),

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

Catchobs <- df_to_matrix(dat[['canum']], season =  1)

nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#')#, skip = 3)
#Feffort[which(df.tmb$years == 2013),2] <- 0 # This is different two places in the input files

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
  CminageSeason = c(0),
  Fmaxage = 2, # Fully selected fishing mortality age
  Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
  Qmaxage = Qmaxage, #Qmaxage = c(1,3)
  # minSDcatch = sqrt(0.01),
  # minSDsurvey = sqrt(0.01),
  # penepsC = 1e-10,
  # penepsCmax = 1e-8,
  # peneps = 1e-10,
  # maxSDcatch = sqrt(2),
  Fbarage = c(1,2),
  #isFseason = c(1,1,1,0), # Seasons to calculate fishing in
  powers = powerIN,
  endFseason = 2, # which season does fishing stop in the final year of data
  nocatch = as.matrix(nocatch),
  surveyStart = surveyStart, #c(0.75,0)
  surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
  surveySD =  surveyCV, #c(1,2)),
  catchSD = list(c(0,1,2),
                 c(0,1,2)),
  estSD = c(0,0,0), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,1) # Factor for relative strength of log-likelihood

)
# Get initial parameter structure
parms <- getParms(df.tmb )

# Get non-estimated parameters, based on info in df.tmb
#parms$logalpha <- log(1287.509)
#parms$logalpha <- log(1269.427)
mps <-getMPS(df.tmb, parms)
#mps$logalpha <- factor(parms$logalpha * NA)
# Set boundaries
# This model works best if SDsurvey is mapped
sas <- runAssessment(df.tmb, parms = parms,mps = mps, silent = TRUE)
sas$reps
sas$opt$time_to_eval

getForecastTable(df.tmb, sas, TACold = 74000, Btarget = 125000, Flimit =  .69)

plot(sas)

mr <- mohns_rho(df.tmb, parms = parms, peels = 5, plotfigure = TRUE) # 1 season has high mohns rho 

saveRDS(sas, file = file.path(wd,'yearly_model.RDS'))

# Save the table 
write.table(mr$df.save, file = file.path(wd,'mohns_table.csv'), row.names = FALSE)


source(file.path('C:/Users/nsja/Dropbox/DTU/BEBRIS',"compare_TMB_admb.R"))

library(tidyverse)
wd_dat <- 'C:/Users/nsja/Dropbox/DTU/SPRAT/SMS_2025/Sprat-div-4_plus_IIIa_S3_19'
x <- compareTMB_ADMB(wd_dat, sas, df.tmb)

ggsave(file.path(wd,'smsR_admb.png'), x)

xx <- list(df.tmb = df.tmb,
           sas = sas)

saveRDS(xx, file = 'sprat_2023.RDS')
# Compare everything with the correctly parameterized one

## SSB
library(tidyverse)

sms <- read.table(file.path('C:/Users/nsja/Dropbox/DTU/BEBRIS', 'compare'), header = TRUE)

# plot(sms$Year, sms$SSB)
# lines(c(df.tmb$years, max(df.tmb$years)+1), sas$x$SSB)
#
# Survey and catch residuals

sresid <- read.table(file.path(wd, 'survey_residuals.out'), sep = ',', fill = TRUE)
cresid <- read.table(file.path(wd, 'catch_residuals.out'), sep = ',', fill = TRUE)
names(cresid) <- c('Species','season','years','0','1','2','3')
cresid <- cresid %>% select(-Species) %>% pivot_longer(3:6, names_to = 'ages', values_to = 'ResidCatch') %>% mutate(model = 'admb', ages = as.numeric(ages))

TB <- getTSB(df.tmb, sas) %>% mutate(Area = 'North Sea', Species = 'sprat')
SSB <- getSSB(df.tmb, sas)

TB$SSB <- SSB$SSB[1:df.tmb$nyears]

# write.csv(TB, file = 'C:/Users/nsja/Dropbox/DTU/Myndighedsbetjening/Økosystem status 2023/summed data/sprat_ns.csv', row.names = FALSE)


# Compare the survey residuals

S1 <- sresid[1:length(1982:2022),4]

S.tmb <- getResidSurvey(df.tmb, sas)
C.tmb <- getResidCatch(df.tmb, sas) %>% mutate(model = 'tmb')


residC <- bind_rows(C.tmb, cresid)

plot(S1/S.tmb$ResidSurvey)

ggplot(residC %>% filter(season == 4, ResidCatch> -10),
       aes(x = years, y = ResidCatch, color = factor(model)))+geom_line()+facet_wrap(~ages)



# All resids


cresid <- read.table(file.path(wd, 'catch_survey_residuals.out'), sep = '', fill = TRUE, header = TRUE)

x0 <- cresid %>% filter(Age == 0, data == 'catch')


ggplot(x0 %>% filter(residual > -99), aes(x = Year, y=  residual))+geom_line()+facet_wrap(~Quarter)+theme_classic()#+geom_point(data = x0 %>% filter(model > 0))

# First survey
surveyresid <- cresid %>% filter(data == 'survey', fleet == 1)

plot(surveyresid$Year, surveyresid$residual/S.tmb$ResidSurvey[S.tmb$survey == 'IBTS_Q1_0'])
points(S.tmb$years[S.tmb$survey == 'IBTS_Q1_0'], S.tmb$ResidSurvey[S.tmb$survey == 'IBTS_Q1_0'], col ='red')

# Try the second survey
snames <- unique(S.tmb$survey)

# Second survey

# First survey
surveyresid <- cresid %>% filter(data == 'survey', fleet == 4) %>% rename(ResidSurvey = residual, ages = Age, years = Year, season = Quarter)

x<- bind_rows(surveyresid %>% mutate(model = 'admb'), S.tmb %>% filter(survey == snames[4]) %>% mutate(model = 'tmb'))

ggplot(x, aes(x= years, y = ResidSurvey, color = model))+geom_line()+facet_grid(~ages)


parms <- scan(file.path(wd,'sms.par'), comment.char = '#')
rec.scale <- scan(file.path(wd,'rec_scale.out'), comment.char = "#", skip = 1)[2]
R <- exp(parms[1:df.tmb$nyears]+rec.scale)
R.tmb <- getR(df.tmb, sas)

plot(R.tmb$R/R)
lines(rep(1, length(R)), lty = 2)

Ninit <- exp(parms[(df.tmb$nyears+1):(df.tmb$nyears+3)]+rec.scale)
Ninit.tmb <- exp(sas$reps$par.fixed[names(sas$reps$par.fixed) == 'logNinit'])

plot(Ninit/Ninit.tmb*100)
#lines(Ninit.tmb)

# For Ole

ggsave(p2$survey,filename = file.path(wd,'survey_fit_tmb.png'), width = 16, height = 12, units = 'cm')
ggsave(p2$catch,filename = file.path(wd,'catch_fit_tmb.png'), width = 20, height = 20, units = 'cm')
