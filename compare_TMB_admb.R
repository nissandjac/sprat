compareTMB_ADMB <- function(wd, sas, df.tmb ){


  # Read general stuff in

  admb <- read.table(file.path(wd, 'summary_table.out'), skip = 4, fill = TRUE)
  names(admb) <- c('year', 'R','SSB','TSB','Chat', 'Fbar' )

  # Fix some stuff due to formatting
  admb$SSB[nrow(admb)] <- admb$year[nrow(admb)]
  admb$year[nrow(admb)] <- admb$year[nrow(admb)-1]+1
  admb$model <- 'admb'

  # Create the same data frame for smsR

  R <- getR(df.tmb, sas)
  SSB <- getSSB(df.tmb, sas)
  Fbar <- getFbar(df.tmb, sas)
  TSB <- getBiomass(df.tmb, sas)
  TSB <- TSB %>% filter(season == 1) %>% group_by(years) %>% summarise(Biomass = sum(Biomass))
  Yield <- getCatch(df.tmb, sas)
  #
  admb.tot <- read.table(file.path(wd, 'summary_table_raw.out'), fill = TRUE, header = TRUE)

  # Do estimated catch instead
  admb$Yield <- admb.tot$Yield.hat


  tmb.fit <- data.frame(year = c(df.tmb$years, max(df.tmb$years)+1),
                        R = c(R$R[1:df.tmb$nyears], NA),
                        SSB = SSB$SSB,
                        TSB = c(TSB$Biomass,NA),
                        Yield = c(Yield$Catch, NA),
                        Fbar = c(Fbar$Fbar,NA),
                        model = 'TMB')


  admb.plot <- admb %>% pivot_longer(c(R,SSB, TSB, Yield, Fbar))
  tmb.plot <- tmb.fit %>% pivot_longer(c(R,SSB, TSB, Yield, Fbar))


  p.all <- ggplot(admb.plot, aes(x = year, y = value, color = model))+geom_point()+
    geom_line(data = tmb.plot)+
    facet_wrap(~name, scales = 'free_y')+
    theme_classic()

  print(p.all)


  # Uncertainty between the two models
  admb.tot <- read.table(file.path(wd, 'summary.out'), fill = TRUE, header = TRUE)

return(p.all)

}
