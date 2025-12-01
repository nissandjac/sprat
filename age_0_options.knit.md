---
title: "0 group options (four quarters)"
author: "Nis Sand Jacobsen"
format:
  html:
    embed-resources: true
---



## Introduction

This document serves to analyse the 0 group parameters. There are three essential parameters here

- Survey power law (on or off for age 0)
- SRR weighting parameter (historically set to 0.1)
- Estimation of the recruitment standard deviation (removes need to weigh SRR manually).

Here I combine all different options in the four season model, leading to 8 different models. All the models converge, but I have not looked in detail into stability or 
The power law is a somewhat troublesome parameter because it interferes with the 0 group catchability. Estimating both those parameters lead to longer estimation times, and makes the Q0 parameter hardly identifiable. 



::: {.cell}

:::



# Differences in estimated quantities 
Neither of the parameters make a big difference in the historical estimates of SSB (@fig-ssb) or recruitment (@fig-R). However, it makes a big difference for the in year advice in the final year. The model without power law and low weight to the stock recruitment relationship predicts an incredibly high recruitment event in the final year. Both the power law and letting the model estimate the weight of the SRR relationship tunes this down. 



::: {.cell}
::: {.cell-output-display}
![Spawning stock biomass (SSB) trajectories by SDR scenario.](age_0_options_files/figure-html/fig-ssb-1.png){#fig-ssb width=100%}
:::
:::

::: {.cell}
::: {.cell-output-display}
![Recruitment (R) trajectories by SDR scenario.](age_0_options_files/figure-html/fig-R-1.png){#fig-R width=100%}
:::
:::


# Retrospective patterns
The power law is often implemented to reduce the retrospective patterns in recruitment. However in this case the retros are only slightly different. The power law model performs slightly better than without the power law for SBB, but performs worse for recruitment (@fig-mohns). 




::: {.cell}
::: {.cell-output-display}
![Absolute values of Mohns rho for the different recruitment configurations](age_0_options_files/figure-html/fig-mohns-1.png){#fig-mohns width=100%}
:::
:::



# AIC 
The AIC values can not really be compared for the cases where the SRR likelihood is weighted by 0.1, as those likelihoods will be articially deflated. If we compare the other options they are very close, however the power law with internally tuned has a slightly lower AIC than the other options. The tuned model without powerlaw performs the worst (@fig-aic). 



::: {.cell}
::: {.cell-output-display}
![AIC of power law and estimated SDR models](age_0_options_files/figure-html/fig-aic-1.png){#fig-aic width=100%}
:::
:::




# Most accurate forecast
Since this choice has large consequences for the in year advice as the recruits will be a part of next years SSB, we tested the forecast abilities of the models by peeling of 5 years and forecasting into the next year. We can then compare the estimated forecast SSB against the SSB estimated in the models. This is almost the same as a Mohns rho analysis, however I am looking only at the in year projection of SSB here, as that value is affected by both the terminal year recruitment and the age 2-4 available SSB. 
For each run I removed a year back in time (a *peel*) and predicted the spawning biomass in the following year (i.e., the advice year). 



::: {.cell}

:::


All the models tended to overestimate the SSB in the advice year, with the exception of the model with tuned Stock recruitment variance, no power law, and no weighting on the likelihood function, which had a median close to 0 of the relative observed minus forecasted values (@fig-ssbest). The models that had a weighting factor on the stock recruitment relationship all tended to overestimate the forecasted SSB (@fig-ssbline). Te trend could be a sign that the natural mortality on 0 year olds is higher than expected. 



::: {.cell}
::: {.cell-output-display}
![In year (SSB-SSBforecast)/SSB for the different recruitment models. The figure shows a violin plot of the 15 peels](age_0_options_files/figure-html/fig-ssbest-1.png){#fig-ssbest width=100%}
:::
:::


123



::: {.cell}
::: {.cell-output-display}
![In year forecast versus the estimated SSB the following year. Dashed line shows the forecasted SSB](age_0_options_files/figure-html/fig-ssbline-1.png){#fig-ssbline width=100%}
:::
:::
