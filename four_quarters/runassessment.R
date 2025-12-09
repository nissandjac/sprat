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
rel_contrib <- prop.table(yield_sum, margin = 1)

nocatch <- as.matrix(nocatch)#*0+1
# nocatch[rel_contrib < .01] <- 0
# Save some data for package example
dat$nocatch <- as.matrix(nocatch)#*0+1
#dat$nocatch[1,4] <- 0
#dat$nocatch[dat$effort == 0] <- 0
powerIN <- list(NA, NA, NA)
#Surveyobs[,years < 2000,2] <- -1

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
  blocks = c(1974,2018),
  Fbarage = c(1,2),
  endFseason = 2, # which season does fishing stop in the final year of data
  nocatch = as.matrix(dat$nocatch),
  surveyStart = surveyStart, #c(0.75,0)
  surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
  surveySD =  surveyCV, #c(1,2)),
  #randomF = 1,
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
#parms$Fseason <- matrix(1, nrow = 1, ncol = max(df.tmb$bidx)+1)
# Get non-estimated parameters, based on info in df.tmb
#parms$logalpha <- log(1287.509)
#parms$logalpha <- log(1269.427)
mps <-getMPS(df.tmb, parms)
#mps$logalpha <- factor(parms$logalpha * NA)
# Set boundaries
# This model works best if SDsurvey is mapped
# df.tmb$randomR = 1
# df.tmb$randomF <- 1
sas <- runAssessment(df.tmb, parms = parms,mps = mps, 
                     silent = TRUE)
sas$reps
sas$opt$time_to_eval

# getForecastTable(df.tmb, sas, TACold = 74000, Btarget = 125000, Flimit =  .69)

# plot(sas)

plotBubbles(sas)
# Save
# p2 <- plotDiagnostics(df.tmb, sas)
# p2$cresids_scaled

mr <- mohns_rho(df.tmb, peels = 5, parms, mps, plotfigure = TRUE)


saveRDS(sas, file.path(wd,'four_seasons.RDS'))
write.table(mr$df.save, file = file.path(wd,'mohns_table.csv'), row.names = FALSE)
write.table(mr$mohns, file = file.path(wd,'mohns_table_tot.csv'), row.names = FALSE)


g   <- sas$obj$gr(sas$opt$par)
ord <- order(abs(g), decreasing = TRUE)
head(cbind(name = names(sas$opt$par)[ord], grad = g[ord]), 10)
problem_par  <- names(g)[ord[1]]
print(problem_par)


