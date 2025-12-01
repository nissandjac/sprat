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

  powerIN <- list(c(0,1), NA, NA)
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

# Save some data for package example
dat$nocatch <- as.matrix(nocatch)#*0+1
#dat$nocatch[1,4] <- 0
#dat$nocatch[dat$effort == 0] <- 0
powerIN <- list(0, NA, NA)


df.power <- get_TMB_parameters(
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
  minSDcatch = sqrt(0.01),
  maxSDcatch = sqrt(1.5),
  # minSDsurvey = sqrt(0.2),
  #penepsC = 1e-10,
  penepsCmax = 1e-8,
  # peneps = 1e-10,
  Fbarage = c(1,2),
  isFseason = c(1,1,1,0), # Seasons to calculate fishing in
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
  estSD = c(0,0,0), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,.1) # Factor for relative strength of log-likelihood

)
df.power$powers[2,1] <- 1
# Get initial parameter structure
parms.power <- getParms(df.power)
# Get non-estimated parameters, based on info in df.tmb

sas.power <- runAssessment(df.power, parms = parms.power, silent = TRUE)
plot(sas.power)



df.nopower <- get_TMB_parameters(
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
  minSDcatch = sqrt(0.01),
  maxSDcatch = sqrt(1.5),
  # minSDsurvey = sqrt(0.2),
  #penepsC = 1e-10,
  penepsCmax = 1e-8,
  # peneps = 1e-10,
  Fbarage = c(1,2),
  isFseason = c(1,1,1,0), # Seasons to calculate fishing in
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
  estSD = c(0,0,2), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,.1) # Factor for relative strength of log-likelihood
  
)
parms <- getParms(df.nopower)
sas.nopower <- runAssessment(df.nopower, parms = parms, silent = TRUE)


df.estRSD <- get_TMB_parameters(
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
  minSDcatch = sqrt(0.01),
  maxSDcatch = sqrt(1.5),
  # minSDsurvey = sqrt(0.2),
  #penepsC = 1e-10,
  penepsCmax = 1e-8,
  # peneps = 1e-10,
  Fbarage = c(1,2),
  isFseason = c(1,1,1,0), # Seasons to calculate fishing in
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
  estSD = c(0,0,0), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,.1) # Factor for relative strength of log-likelihood
  
)

parms <- getParms(df.estRSD)
sas_nopower_estRSD <- runAssessment(df.estRSD, parms = parms, silent = TRUE)

df.power.estRSD <- get_TMB_parameters(
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
  minSDcatch = sqrt(0.01),
  maxSDcatch = sqrt(1.5),
  # minSDsurvey = sqrt(0.2),
  #penepsC = 1e-10,
  penepsCmax = 1e-8,
  # peneps = 1e-10,
  Fbarage = c(1,2),
  isFseason = c(1,1,1,0), # Seasons to calculate fishing in
  powers = c(0, NA, NA),
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
  estSD = c(0,0,0), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,.1) # Factor for relative strength of log-likelihood
  
)
parms <- getParms(df.power.estRSD)
sas.power.estRSD <- runAssessment(df.power.estRSD, parms = parms, silent = TRUE)


df.estRSD_SR <- get_TMB_parameters(
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
  minSDcatch = sqrt(0.01),
  maxSDcatch = sqrt(1.5),
  # minSDsurvey = sqrt(0.2),
  #penepsC = 1e-10,
  penepsCmax = 1e-8,
  # peneps = 1e-10,
  Fbarage = c(1,2),
  isFseason = c(1,1,1,0), # Seasons to calculate fishing in
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
  estSD = c(0,0,0), # Estimate
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,1) # Factor for relative strength of log-likelihood
  
)

parms <- getParms(df.estRSD_SR)
sas_nopower_estRSD_SR <- runAssessment(df.estRSD_SR, parms = parms, silent = TRUE)

library(dplyr)
library(ggplot2)
library(tidyr)
sas_2024 <- readRDS('sprat_2025.RDS')[[1]]

models <- list(
  "no power - estRSD - SR" = sas_nopower_estRSD_SR,
  "nopower estRSD"    = sas_nopower_estRSD,
  "power - estRSD"      = sas.power.estRSD,
  "power"             = sas.power,
  "2024 model"              = sas_2024,
  "no power"           = sas.nopower
)

## 2. Extract SSB for each model and bind into one data frame ----
ssb_all <- lapply(names(models), function(mname) {
  mod <- models[[mname]]
  
  # getSSB(df, sas) where df = mod$dat, sas = model object
  tmp <- getSSB(mod$dat, mod)
  
  # Add model name
  tmp$model <- mname
  
  tmp
}) |> 
  bind_rows()

# Optional: make sure years is numeric (if it's not already)
ssb_all <- ssb_all |>
  mutate(years = as.numeric(years))

## 3. Plot SSB trajectories for all models ----
p <- ggplot(ssb_all, aes(x = years, y = SSB, colour = model)) +
  geom_line(linewidth = .8) +
  # Optional: add uncertainty ribbons if you like
  # geom_ribbon(aes(ymin = low, ymax = high, fill = model),
  #             alpha = 0.15, colour = NA) +
  labs(x = "Year",
       y = "Spawning stock biomass (SSB)",
       colour = "Model",
       title = "") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank()
  )

print(p)

# Now compare R 


R_all <- lapply(names(models), function(mname) {
  mod <- models[[mname]]
  
  # getSSB(df, sas) where df = mod$dat, sas = model object
  tmp <- getR(mod$dat, mod)
  
  # Add model name
  tmp$model <- mname
  
  tmp
}) |> 
  bind_rows()

# Optional: make sure years is numeric (if it's not already)
R_all <- R_all |>
  mutate(years = as.numeric(years))

## 3. Plot SSB trajectories for all models ----
p.R <- ggplot(R_all, aes(x = years, y = R, colour = model)) +
  geom_line(linewidth = .8) +
  # Optional: add uncertainty ribbons if you like
  # geom_ribbon(aes(ymin = low, ymax = high, fill = model),
  #             alpha = 0.15, colour = NA) +
  labs(x = "Year",
       y = "Recruitment",
       colour = "Model",
       title = "") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank()
  )

print(p.R)


ggsave('compare_SSB_inputs.png',p, width = 16, height = 12, units = 'cm')
ggsave('compare_R_inputs.png',p.R, width = 16, height = 12, units = 'cm')


# Make a table of AIC values 
## 2. Extract AIC values into a table
AIC_table <- lapply(names(models), function(mname) {
  data.frame(
    model = mname,
    AIC   = AIC(models[[mname]])
  )
}) |> 
  dplyr::bind_rows() |> dplyr::arrange(AIC)

AIC_table

write.table(AIC_table,'aic_0year.csv', row.names = FALSE)
