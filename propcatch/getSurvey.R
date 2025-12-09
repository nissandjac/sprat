getSurvey <- function(survey_df,
                      maxage,
                      survey.age, # Ages in the two surveys
                      survey.years,
                      survey.names,
                      survey.quarter,
                      years = years,
                      seasons = seasons){

survey <- survey[1:length(unlist(survey.years)), ]


ages <- 0:maxage


survey.out <- matrix(NA, nrow = nrow(survey), ncol = length(ages) + 4) # Plus efficiency, name, quarter and years
survey.out <- as.data.frame(survey.out)
colnames(survey.out) <- c("eff", ages, "Survey", "Quarter", "year")
survey.out$eff <- survey[, 1]


for (k in 1:length(survey.names)) {
  sname <- rep(survey.names[k], length(survey.years[[k]]))
  syears <- survey.years[[k]]
  
  
  if (k == 1) {
    idx <- 1:length(syears)
  } else {
    idx <- (idx[length(idx)] + 1):((idx[length(idx)]) + length(syears))
  }
  
  
  survey.out$Survey[idx] <- sname
  survey.out$year[idx] <- syears
  survey.out$Quarter[idx] <- rep(survey.quarter[k], length(idx))
  
  
  for (i in 1:length(survey.age[[k]])) {
    survey.out[idx, names(survey.out) == survey.age[[k]][i]] <- survey[idx, i + 1]
  }
}

survey.out <- survey.out %>% tidyr::pivot_longer(2:(length(ages) + 1), names_to = "Age", values_to = "cpue")

survey <- survey.out
survey$cpue[is.na(survey$cpue)] <- -1

return(survey)
}
