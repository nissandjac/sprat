# Load Catch and write to files # 
library(tidyverse)

path <- 'data/catches_vs2/'
fls <- dir(path)
fl_names <- c('old_method','four_quarters','two_seasons','yearly_model')

wd <- 'C:/Users/nsja/Dropbox/DTU/BEBRIS/'

# Age is independent on input data 
maxage <- 3
ages <- 0:maxage
nage <- length(ages)
start_year <- 1973
cyear <- 2024



for(i in 1:length(fl_names)){

  
  dat <- read.csv(file.path(wd,path,fls[i]))
  names(dat)[3] <- 'season'
  dat <- dat %>% arrange(year, season) %>% filter(year > 1973) 
  ifelse(!dir.exists(file.path(wd, fl_names[i])), dir.create(file.path(file.path(wd, fl_names[i]))), FALSE)
  
  if(i == 4){
    names(dat)[3] <- 'n0'
    dat$season <- 1
  }
  
  styr <- min(dat$year)
  if(start_year > styr)styr <- start_year
  endyr <- cyear
  years <- styr:endyr
  nyear <- length(styr:endyr)
  nseason <- length(unique(dat$season))  
  
  tmp.weca <- dat %>% select(paste('mw',0:maxage, sep=''), year, season) %>% mutate(hashtag = '#') %>% 
    filter(year <= cyear ) %>% 
    relocate(paste('mw',0:maxage, sep=''), hashtag) %>% arrange(year, season)

  # Add 2025 (5 year average) (fix to whatever is used for sprat)
  add_rows <- lapply(1:nseason, function(j) {
    tmp.weca %>%
      dplyr::filter(year %in% (endyr-4):endyr, season == j) %>%
      dplyr::summarise(dplyr::across(1:nage, ~ mean(.x, na.rm = TRUE))) %>%
      dplyr::mutate(hashtag = "#", year = endyr + 1, season = j)
  })
  
  tmp.weca.out <- dplyr::bind_rows(tmp.weca, dplyr::bind_rows(add_rows)) %>%
    dplyr::select(dplyr::all_of(names(tmp.weca)))  # keep original column orde
  
  
  
  write.table(tmp.weca.out, file = file.path(wd,fl_names[i],'weca.in'), row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(tmp.weca.out, file = file.path(wd,fl_names[i],'west.in'), row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  
  # canum
  
  tmp.canum <- dat %>% select(paste('n', 0:maxage, sep=''), year, season) %>% mutate(hashtag = '#') %>%
    relocate(paste('n',0:maxage,sep=''),hashtag) %>% arrange(year, season) %>% filter(year <= cyear)
  
  
  write.table(tmp.canum, file = file.path(wd,fl_names[i],'canum.in'), row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  # Create a new nocatch matrix 
  
  mincatch <- 537518.7
  mincatch <- 100000
  
  # Create matrix of 0s and 1s
  
  nocatch <- as.matrix(rowSums(tmp.canum[, c("n0", "n1", "n2", "n3")])) > mincatch
  nocatch <- 1 * nocatch  # Convert TRUE/FALSE to 1/0  
  nocatch <- as.data.frame(nocatch) %>% mutate(hashtag = '#', seasons = rep(1:nseason, nyear),
                                               years = rep(years, each = nseason)) %>%
    pivot_wider(
      names_from = seasons,     # make one column per season
      values_from = V1,        # values come from V1
      names_prefix = "season"  # optional: prefix column names
    ) %>% relocate(paste('season',1:nseason, sep=''))
  
  write.table(nocatch, file = file.path(wd,fl_names[i],'zero_catch_year_season.in'), 
              row.names = FALSE, quote = FALSE, col.names = FALSE)
  
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
  
  if( i == 3){
    
  
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
              row.names = FALSE, quote = FALSE, col.names = FALSE)}

  
  if(i == 4){
    M2_sum <- M2 %>%
      pivot_longer(
        cols = starts_with("V"),      # V1–V4
        names_to = "ageclass",
        values_to = "value"
      ) %>%
      group_by(years, ageclass) %>%
      summarise(M2 = sum(value), .groups = "drop")
    M2_wide <- M2_sum %>%
      pivot_wider(names_from = ageclass, values_from = M2) %>% 
      mutate(hashtag = '#') %>% 
      relocate(paste('V',1:4, sep =''), hashtag) 
    
    names(M2_wide)[1:nage] <- paste('Age',0:maxage, sep ='')
    
    write.table(M2_wide, file = file.path('C:/Users/nsja/Dropbox/DTU/BEBRIS/',fl_names[i],'natmor.in'),
                row.names = FALSE, quote = FALSE, col.names = FALSE)
    
    
}

  # Read the natural mortality 
  # source('parse_sms_mor.R')
  # natmor <- parse_sms_mor('data/natmorm1m2.out') %>% filter(age > -1)
  # # Get sprat 
  # spr_m2 <- natmor %>% filter(species_id == 25, year > (styr-1)) %>% 
  #   select(-c(species_id, species_name))
  # 
  # # Change the quarters and years
  # 
  # 
  # 
  # # # Now correct the seasons to sprat seasons 
  # # natmor_sprat <- spr_m2 %>%
  # #   # 1) go back one year for ORIGINAL seasons 1 and 2
  # #   mutate(year = if_else(season %in% c(1L, 2L), year - 1L, year)) %>%
  # #   # 2) remap seasons: 1→3, 2→4, 3→1, 4→2
  # #   mutate(season = recode(as.integer(season),
  # #                          `1` = 3, `2` = 4, `3` = 1, `4` = 2)) %>%  
  # #   arrange(year,season,age) 
  # # 
  # natmor_sprat <- spr_m2 %>%
  #   mutate(season_orig = as.integer(season)) %>%
  #   # move YEAR back for original seasons 1 & 2
  #   mutate(year = if_else(season_orig %in% c(1L, 2L), year - 1L, year)) %>%
  #   # AGE shift: if we moved to previous year (orig seasons 1 & 2), decrement age by 1
  #   mutate(age = if_else(season_orig %in% c(1L, 2L), pmax(age - 1L, 0L), age)) %>%
  #   # (optional) keep plus group at max_age if you have one
  #   #mutate(age = if (!is.null(maxage)) pmin(age, maxage) else age) %>%
  #   # remap seasons: 1→3, 2→4, 3→1, 4→2
  #   mutate(season = recode(season_orig, `1` = 3L, `2` = 4L, `3` = 1L, `4` = 2L)) %>%
  #   select(-season_orig) %>%
  #   arrange(year, season, age) %>% filter(year >= styr)
  # 
  # 
  # 
  # natmor_roll <- natmor_sprat %>%
  #   arrange(age, season, year) %>%
  #   group_by(age, season) %>%
  #   mutate(value_roll3 = zoo::rollmean(value, k = 3, fill = 'extend')) %>%
  #   ungroup() %>% arrange(year,season )
  # 
  # # Now add season 1 and 2 in the first year 
  # #  first_year <- min(natmor_roll$year)
  # # rows to add: take seasons 1 & 2 from (first_year + 1), make them belong to first_year
  # rows_to_add <- natmor_roll %>%
  #   filter(year == first_year + 1, season %in% c(1, 2)) %>%
  #   mutate(year = first_year)
  # 
  # # avoid duplicating if those rows somehow already exist
  # keys <- c("year", "season", "age")
  # rows_to_add <- rows_to_add %>%
  #   anti_join(natmor_roll, by = keys)
  # 
  # # Left side fixed 
  # natmor_roll_filled <- bind_rows(natmor_roll, rows_to_add) %>%
  #   arrange(year, season, age)
  # 
  # # Now fix the years after the Norht Sea sms has been run 
  # df <- natmor_roll_filled
  # last_year     <- max(df$year)
  # years_to_add  <- (last_year + 1):(endyr+1)
  # 
  # tmpl <- df %>% filter(year == last_year)
  # 
  # if(max(tmpl$season) != nseason){
  #   
  #   # Get the missing season from two years ago 
  #   idx <- which(!1:nseason %in% unique(tmpl$season))
  #   tmpl_old <- df %>% filter(year == (last_year-1),
  #                             season %in% idx)
  #   
  #   # Add  to data.frame 
  #   df <- rbind(df, tmpl_old %>% mutate(year = last_year))
  #   
  #   tmpl <- rbind(tmpl_old, tmpl) %>% mutate(year = max(tmpl$year)) %>% arrange(season)
  #   
  # }
  # 
  # 
  # # Create copies of the last year for each new year
  # new_rows <- do.call(dplyr::bind_rows,
  #                     lapply(years_to_add, function(y) dplyr::mutate(tmpl, year = y)))
  # 
  # # If some of these years already exist, avoid duplicates
  # new_rows <- new_rows %>% anti_join(df, by = c("year","season","age"))
  # 
  # # Final output 
  # natmor_roll_ext <- bind_rows(df, new_rows) %>%
  #   arrange(year, season, age)
  # 
  # # Now write as data files 
  # 
  # 
  # if(nseason == 4){
  # m2.mtrx <- natmor_roll_ext %>% select(-value) %>% pivot_wider(names_from = age, values_from = value_roll3) %>% 
  #   mutate(hashtag = '#') %>% relocate(paste(0:maxage), hashtag, season, year)
  # }else{ # Only works for two seasons 
  #   m2.mtrx <- natmor_roll_ext %>%
  #     mutate(season2 = if_else(season %in% c(1, 2), 1, 2)) %>%
  #     group_by(year, season2, age) %>%
  #     summarise(
  #       value_roll3 = sum(value_roll3, na.rm = TRUE),
  #       .groups = "drop"
  #     ) %>%
  #     arrange(year, season2, age)  %>% pivot_wider(names_from = age, values_from = value_roll3) %>% 
  #     mutate(hashtag = '#') %>% relocate(paste(0:maxage), hashtag, season2, year)
  #   
  #   
  # }
  # 
  # write.table(m2.mtrx, 
  #             file = file.path(wd,fl_names[i],'natmor.in'), row.names = FALSE, quote = FALSE, col.names = FALSE)
  # 
  
# Now do the surveys 
  
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
  
  
  
}

