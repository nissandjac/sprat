# Potential assessment models for North Sea Sprat
Nis Sand Jacobsen and Ole Henriksen

## Introduction

This repository serves to compare different potential model
configurations for the North Sea sprat. All the models we will show are
based on the seasonal $smsR$ model. Until 2025 the model has run with
ADMB, but here we use the smsR package available at

> [!NOTE]
>
> <https://github.com/nissandjac/smsr>

$smsR$ is an R package which can be used to set up seasonal assessment
through few R functions. Each of the functions offer simple to complex
model setups, which can both respect historical assessment methods and
more modern implementation of e.g., random effects. The package is
currently used as assessments for the four North Sea Sandeel stocks, and
a thorough comparison of the TMB and ADMB performance can be seen in the
WKSANDEEL report.

We have chosen model setups that are close to the previous ones, but all
the models we will provide here has a transparent data flow, and try to
follow best statistical practices. Over the history of the stock many
small modifications have been applied to the model setup and input data
to a point where some of the changes are not traceable anymore. By
resetting all these changes we provide options for a model with annual,
half-year, and quarterly time steps.

## Current model and issues

### Ole add general description about the seasons and data here

Put the text here Ole

### Issues with the current assessment

The current sprat model originated in the 2018(?) sprat benchmark. In
the past years the model has had convergence issues, particularly with
the maximum gradient being above an acceptable level. The issue has been
solved by taking the 0-group from the Q1 survey and scaling it
differently than the other surveys. $smsR$ does not internally scale the
surveys, so this method does not change that issue. The core of the
issue lies elsewhere; there is large confounding between the parameter
estimating survey density dependence (the so-called *power law*) and the
survey catchability parameter *Q*. The survey *power law* is implemented
for two reasons that are linked

1)  Avoid retrospective patterns in recruitment, as the in-year
    assessment can tend to overestimate the recruitment when the only
    data point is one observation in Q1

2)  Reduce the risk that one observation point leads to a (too) high in
    year advice, since the 0 year olds in the spring become part of the
    SSB

$$
N_{i,\mathrm{survey}} = q_i \, N_i^{p_i}
$$

Another issue has been massive residual patterns, which has been solved
by moving the catches into seasons where they did not originate. The
catches from season 4 has since the benchmark been moved into season 3,
except for the inital model year (as ADMB sms was then not able to
converge). There is some confusion to where the catches in season 4 in
the initial model year come from in the current data files, as the data
has been added one year at a time.

Additionally, several of the models parameters are estimated on the
boundaries, leading to slow model convergence in TMB. ADMB has a
particular feature that ‘nudges’ parameters away from boundary
conditions leading to the perception that a parameter has been correctly
estimated, when in reality it is stuck in an infinitely small difference
between the estimated parameter and the set boundary.

The

## Comparing the number of seasons
