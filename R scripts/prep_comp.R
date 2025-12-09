# Load Catch and write to files #
library(tidyverse)

path <- 'propcatch'
fls <- "no_prop_jul_jun_halfyear_IV_with_Q2_no_SMS_no_SEcoast_v29.csv"
fl_names <- c('propcatch')

wd <- 'C:/Users/nsja/Dropbox/DTU/BEBRIS/'

# Age is independent on input data
maxage <- 3
ages <- 0:maxage
nage <- length(ages)
start_year <- 1973
cyear <- 2024

i <-1
  dat <- read.csv(file.path(wd,path,fls[i]))
  names(dat)[3] <- 'season'
  dat <- dat %>% arrange(year, season) %>% filter(year > 1973)
  ifelse(!dir.exists(file.path(wd, fl_names[i])), dir.create(file.path(file.path(wd, fl_names[i]))), FALSE)

  
  styr <- min(dat$year)
  if(start_year > styr)styr <- start_year
  endyr <- cyear
  years <- styr:endyr
  nyear <- length(styr:endyr)
  nseason <- length(unique(dat$season))

  # Proportions in catch
  tmp.prop <- dat %>% select(paste('n',0:maxage,'_prop', sep=''), year, season) %>% 
    mutate(hashtag = '#') %>%
    filter(year <= cyear ) %>%
    relocate(paste('n',0:maxage,'_prop', sep=''), hashtag) %>% arrange(year, season)

  write.table(tmp.prop, file = file.path(wd,fl_names[i],'prop.in'), 
              row.names = FALSE, quote = FALSE, col.names = FALSE)
 
  # canum

  tmp.ton <- dat %>% select(ton, year, season) %>% mutate(hashtag = '#') %>%
    relocate(ton,hashtag) %>% arrange(year, season) %>% filter(year <= cyear)


  write.table(tmp.ton, 
              file = file.path(wd,fl_names[i],'ton.in'), row.names = FALSE, quote = FALSE, col.names = FALSE)

  # Number of observations 
  tmp.nobs <- dat %>% select(n_samples, year, season) %>% 
    mutate(hashtag = '#') %>%
    filter(year <= cyear ) %>%
    relocate(n_samples, hashtag) %>% arrange(year, season)
  tmp.nobs$n_samples[is.na(tmp.nobs$n_samples)] <- 0
  write.table(tmp.nobs, 
              file = file.path(wd,fl_names[i],'nsamples.in'), row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  # Maturity
  # Fix maturity
  mat <- c(0, 0.4, 0.9, 1) # age 4 is = 1, but is not used in the assessment


  mat.mtrx <- matrix(mat, nrow = (nyear+1)*nseason, ncol  = length(0:maxage), byrow = TRUE)
  mat.mtrx <- as.data.frame(mat.mtrx) %>% mutate(hashtag = '#',year = rep(c(years,endyr+1), each = nseason),
                                                 season = rep(1:nseason, nyear+1))
  write.table(mat.mtrx,
              file = file.path(wd,fl_names[i],'propmat.in'), row.names = FALSE, quote = FALSE, col.names = FALSE)


  # Fix Natural mortality for two seasons
  M2 <- read.table('C:/Users/nsja/Dropbox/DTU/SPRAT/SMS_2025/Sprat-div-4_plus_IIIa_age3_S3/natmor.in')
  M2 <- M2 %>% mutate(
    years = rep((start_year+1):cyear, each = 4),
    season = rep(1:4, length((start_year+1):cyear)))
  M2$season[M2$season %in% c(1,2)] <- 1
  M2$season[M2$season %in% c(3,4)] <- 2
  M2_sum <- M2 %>%
      pivot_longer(
        cols = starts_with("V"),      # V1–V4
        names_to = "ageclass",
        values_to = "value"
      ) %>%
      group_by(years, season, ageclass) %>%
      summarise(M2 = sum(value), .groups = "drop")

  M2_wide <- M2_sum %>%
      pivot_wider(names_from = ageclass, values_from = M2) %>%
      mutate(hashtag = '#') %>%
      relocate(paste('V',1:4, sep =''), hashtag)
    names(M2_wide)[1:nage] <- paste('Age',0:maxage, sep ='')

  write.table(M2_wide, file = file.path('C:/Users/nsja/Dropbox/DTU/BEBRIS/',fl_names[i],'natmor.in'),
                row.names = FALSE, quote = FALSE, col.names = FALSE)


  
  Q1 <- 'data/spratindices-nohauldur/NSspratw3a-IBTS-Q1.txt'
  dat <-read.table(Q1, skip = 1, fill = TRUE)
  idx_years <- dat[1,1]:dat[1,2] - 1 # Subtract one because of the timing mismatch
  idx_season <- c(dat[2,3],dat[2,4])
  idx_age <- dat[3,1]:dat[3,2] - 1 # Subtract one because of the timing mismatch

  Q1.export <- data.frame(dat[4:nrow(dat),2:ncol(dat)]) %>%
    mutate(hashtag = '#', years = idx_years, eff = 0.1) %>%
    relocate(eff)
  names(Q1.export)[2:(length(idx_age)+1)] <- paste('age',idx_age, sep='')


  Q3 <- 'data/spratindices-nohauldur/NSspratw3a-IBTS-Q3.txt'
  dat <-read.table(Q3, skip = 1, fill = TRUE)
  idx_years <- dat[1,1]:dat[1,2] # Dont subtract for Q3
  idx_season <- c(dat[2,3],dat[2,4])
  idx_age <- dat[3,1]:dat[3,2]# Don't Subtract one because of the timing mismatch

Q3.export <- data.frame(dat[4:nrow(dat),2:ncol(dat)]) %>%
    mutate(hashtag = '#', years = idx_years, eff = 0.1)  %>% # Revisit this value
    relocate(eff)
  names(Q3.export)[2:(length(idx_age)+1)] <- paste('age',idx_age, sep='')

  # remove age 0 (why) and group age 3 and 4
  Q3.export$age3plus <- Q3.export$age3+Q3.export$age4
  Q3.export <- Q3.export %>% select(-c(age0,age3,age4)) %>%
    relocate(eff, paste('age',1:2, sep=''),'age3plus')
  #  age_str <- paste('Ages =', )

  # Load the heras data #
  heras <- read.table(file.path(wd, 'data/heras.txt'))
  names(heras) <- c('eff',paste('age',1:3, sep=''))
  heras <- heras %>% mutate(hashtag = '#', years = 2006:2024)


  # Output all the survey files
  outfile <- file.path(wd,fl_names[i],'fleet_catch.in')

  # Open a file connection
  con <- file(outfile, open = "wt")

  # ---- Write the first data frame ----
  writeLines("# ---- Q1 IBTS ----", con)
  writeLines(paste("#", paste(names(Q1.export), collapse = " ")), con)
  write.table(Q1.export, con, row.names = FALSE, col.names = FALSE ,quote = FALSE)
  writeLines("", con)  # blank line

  # ---- Write the second data frame ----
  writeLines("# ---- Q3 IBTS ----", con)
  # Write column names as commented line
  writeLines(paste("#", paste(names(Q3.export), collapse = " ")), con)
  # Write the data (no header this time)
  write.table(Q3.export, con, row.names = FALSE, col.names = FALSE, quote = FALSE)
  writeLines("", con)


  # ---- Write the third data frame ----
  writeLines("# ---- HERAS ----", con)
  writeLines(paste("#", paste(names(heras), collapse = " ")), con)
  write.table(heras, con, row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Close the file
  close(con)




