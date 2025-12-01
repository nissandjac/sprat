### Plot all the results ### 
library(tidyverse)

dat <- readRDS('recruitment_parameters/all_sims.RDS')
scenarios <- 1:length(dat)


out.tot <- dat[[1]]$out %>% mutate(scenarios = scenarios[1])
parms <- as.data.frame(dat[[1]]$info$parameters) %>% mutate(scenarios = scenarios[1])
mr <- as.data.frame(dat[[1]]$info$mohns) %>% 
  mutate(scenarios = scenarios[1]) %>% mutate(AIC = dat[[1]]$info$AIC)


for(i in 2:length(dat)){
out.tot <- rbind(out.tot, dat[[i]]$out %>% mutate(scenarios = scenarios[i]))  
parms <- rbind(parms, as.data.frame(dat[[i]]$info$parameters) %>% mutate(scenarios = scenarios[i]))  
mr <- rbind(mr, as.data.frame(dat[[i]]$info$mohns%>% 
                                mutate(AIC = dat[[i]]$info$AIC)) %>% mutate(scenarios = scenarios[i])) 
  }


df.plot <- left_join(out.tot, parms, by = 'scenarios')

df.plot <- df.plot %>%
  mutate(
    SDR = case_when(
      SDR == 0 ~ "estimated",
      SDR == 2 ~ "tuned",
      TRUE     ~ as.character(SDR)   # fallback if other values appear
    )
  )

df.plot$label <- paste('PL=',df.plot$power, '- SDR'=df.plot$SDR, '- nll weight=',df.plot$SRR)
# Plot things 
ggplot(df.plot, aes(x = years, y= SSB, color = label))+
  geom_line()+theme_classic()+facet_wrap(~SDR)+
  theme(legend.position = 'top', legend.title = element_blank())

mr <- mr %>%
  left_join(parms, by = 'scenarios') %>%
  mutate(
    power = ifelse(is.na(power), "PLoff", "PLon"),
    SDR   = ifelse(SDR == 0, "SDRest", "SDRtun"),
    SRR   = paste0("SRR", SRR),
    label = paste(power, SDR, SRR, sep = " | ")
  ) 
mr.plot <- mr %>%
  pivot_longer(1:3, values_to = 'Mohns')

ggplot(mr.plot, aes(x = scenarios, y = Mohns, fill = label))+
  geom_col()+facet_wrap(~name)+theme_minimal()+labs(y = expression("Mohn's " * rho))+
  theme(legend.position = 'top',legend.title = element_blank())

ggplot(mr %>% filter(SRR =='SRR1'), aes(x = label, y = AIC, color = power))+
  geom_point()+theme_minimal()+labs(y = expression("AIC"))+
  theme(legend.position = 'none',legend.title = element_blank())+
  labs(x =NULL)


# Check the forecasting of the 8 different runs 

dat <- readRDS('recruitment_parameters/forecast_sims.RDS')





# Check one model
ssbest <- dat[[1]]$info$SSBest
yrs <- (2025-length(ssbest)+1):2025
ggplot(dat[[1]]$out, aes(x= years, y= SSB))+geom_line()+
  geom_line(data = data.frame(years = yrs, SSB = ssbest), linetype = 2)
  


df.compare <- data.frame(SSB_real = dat[[1]]$out$SSB[dat[[1]]$out$years %in% yrs], years = yrs,
                         SSB_est = ssbest,scenarios = 1)

for(i in 2:length(dat)){
tmp <- dat[[i]]
ssbest <- tmp$info$SSBest
df.compare <- rbind(
  df.compare, 
  data.frame(SSB_real = tmp$out$SSB[tmp$out$years %in% yrs], years = yrs,
             SSB_est = ssbest,scenarios = i)
  )
  
}  
  

df.compare <- left_join(df.compare, parms, by = 'scenarios')

df.compare <- df.compare %>%
  mutate(
    SDR = case_when(
      SDR == 0 ~ "estimated",
      SDR == 2 ~ "tuned",
      TRUE     ~ as.character(SDR)   # fallback if other values appear
    )
  )
df.compare$label <- paste('PL=',df.compare$power, '- SDR'=df.compare$SDR, '- nll weight=',df.compare$SRR)

df.compare2 <- df.compare %>%
  mutate(
    power_ord = ifelse(power == 0, 0, 1)  # 0 first, others second
  ) %>%
  arrange(power_ord, label) %>%
  distinct(label, power_ord) %>%     # unique rows for levels
  pull(label) -> label_levels
df.compare2 <- df.compare %>%
    mutate(
      label = factor(label, levels = label_levels)
    )
  
  
p1 <- ggplot(df.compare2, aes(x = scenarios, y = (SSB_real-SSB_est)/SSB_real, fill = label))+
  geom_violin()+theme_classic()+geom_hline(aes(yintercept = 0), linetype = 2)+
  stat_summary(fun = median, geom = "crossbar", width = 0.5, 
               fatten = 0, color = "black") +
  theme(legend.position = 'top',
        legend.title = element_blank())+
  facet_wrap(~power)
p1

p2 <- ggplot(df.compare, aes(x = years, y = SSB_real, color = factor(scenarios)))+
  geom_line(show.legend = FALSE)+
  geom_line(aes(y  = SSB_est), linetype = 2,show.legend = FALSE)+
  theme_classic()+
  #geom_hline(aes(yintercept = 0), linetype = 2)+
  facet_wrap(~label, nrow =2)+ylab('SSB')
p2
