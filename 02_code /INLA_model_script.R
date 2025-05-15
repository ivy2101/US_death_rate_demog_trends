library(dplyr)
library(here)
library(haven)
library(stringr)
library(readr)
library(purrr)
library(readxl)
library(tidyr)
library(posterior)
library(ggplot2)
library(coda)
library(spdep)
library(INLA)
library(ggplot2)
library(ggforce)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(tidyverse)

all_causes_summary <- read_rds("/Users/irisyang/Downloads/summary_1970_2022.rds")
unique(all_causes_summary$age)
age_levels <- c("<1", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
                "35-39", "40-44", "45-49", "50-54", "55-59", "60-64",
                "65-69", "70-74", "75-79", "80-84", "85+")
race_levels <- c("White", "Black or African American", "Asian or Pacific Islander",
                 "American Indian or Alaska Native", "All other races")

subset_test_f <- all_causes_summary %>%
  filter(year >= 1973, year <= 2022) %>%
  filter(sex == 'Female') %>%
  filter(population != 0) %>%
  filter(deaths <= population) %>%
  filter(!is.na(population)) %>%
  filter(!is.na(age)) %>% 
  filter(!is.na(state)) %>%
  filter(!is.na(year)) %>% 
  filter(!is.na(race)) %>%
  #slice_sample(n = 50000) %>% 
  mutate(
    x = year - 1972,
    age = factor(age, levels = age_levels, ordered = TRUE),
    race = factor(race, levels = race_levels, ordered = TRUE)
  )  %>%
  mutate(
    
    id_3way = as.integer(interaction(state, age, race, drop = TRUE)),
    id_3way2 = id_3way,
    
    state2 = state,
    state3 = as.integer(factor(state)),
    state4 = as.integer(factor(state)),
    state5 = as.integer(factor(state)),
    state6 = as.integer(factor(state)),
    age2 = age,
    age3 = age,
    age4 = age,
    age5 = age,
    age6 = age,
    race2 = race,
    race3 = as.integer(factor(race)),
    race4 = as.integer(factor(race)),
    race5 = race,
    race6 = race,
    year2 = x,
    e = 1:n()) %>%
  arrange(state, age, race, year) 

subset_test_m <- all_causes_summary %>%
  filter(year >= 1973, year <= 2022) %>%
  filter(sex == 'Male') %>%
  filter(population != 0) %>%
  filter(deaths <= population) %>%
  filter(!is.na(population)) %>%
  filter(!is.na(age)) %>% 
  filter(!is.na(state)) %>%
  filter(!is.na(year)) %>% 
  filter(!is.na(race)) %>%
  #slice_sample(n = 50000) %>% 
  mutate(
    x = year - 1972,
    age = factor(age, levels = age_levels, ordered = TRUE),
    race = factor(race, levels = race_levels, ordered = TRUE)
  )  %>%
  mutate(
    
    id_3way = as.integer(interaction(state, age, race, drop = TRUE)),
    id_3way2 = id_3way,
    
    state2 = state,
    state3 = as.integer(factor(state)),
    state4 = as.integer(factor(state)),
    state5 = as.integer(factor(state)),
    state6 = as.integer(factor(state)),
    age2 = age,
    age3 = age,
    age4 = age,
    age5 = age,
    age6 = age,
    race2 = race,
    race3 = as.integer(factor(race)),
    race4 = as.integer(factor(race)),
    race5 = race,
    race6 = race,
    year2 = x,
    e = 1:n()) %>%
  arrange(state, age, race, year)  

## INLA Code for national seasonal model
# INLA requires unique names for the first of each in f(). So for example you set:
# dat$age2 = dat$age
# etc. because f(AGE,...) and f(AGE2,...) need to be different
# The second term in each f() (e.g., f(age2, YEAR2,...)) can be repeated as the second term
# To test, you can try all 'rw1' as 'iid' first
# Please try most basic model to check that INLA runs first,
# i.e., fml <- deaths ~ 1 + year.month
## dat - data frame in long form
## year.month, year.month2, ... - month and year over time, values 1-N_yearmonth
## month1, month2, ... - indices for month, values 1-12
## e - column of 1...N_rows_dataframe
# hyperparameter value and hyperprior
loggamma_prior = function(gamma.shape = 0.001, gamma.rate = 0.001) {
  as.character(glue::glue(
    'list(prec = list(prior = "loggamma", param = ',
    'c({gamma.shape}, {gamma.rate})))'
  ))
}
hyperprior = loggamma_prior()

#INLA formula
fml <- deaths ~ 1 + # global intercept
  x + # global annual slope
  # state specific terms (iid)
  f(state, model='iid', hyper = hyperprior, constr = TRUE) + # state specific intercept
  f(state2, x, model='iid', hyper = hyperprior, constr = TRUE) + # state specific slope
  # age group specific terms (random walk)
  f(age, model='rw1', hyper = hyperprior, constr = TRUE) + # age group specific intercept
  f(age2, x, model='rw1', hyper = hyperprior, constr = TRUE) + # age group specific slope
  # race specific terms (iid)
  f(race, model='iid', hyper = hyperprior, constr = TRUE) + # race specific intercept
  f(race2, x, model='iid', hyper = hyperprior, constr = TRUE) + # race specific slope
  # age-race specific intercept and slope
  f(age3, model = 'rw1', group = race3, control.group = list(model = "iid", hyper = hyperprior),
    hyper = hyperprior, fixed = FALSE, constr = TRUE) +
  f(age4, x, model = 'rw1', group = race4, control.group = list(model = "iid", hyper = hyperprior),
    hyper = hyperprior, fixed = FALSE, constr = TRUE) +
  # age-state specific intercept and slope
  f(age5, model = 'rw1', group = state3,
    control.group = list(model = "iid", hyper = hyperprior),
    hyper = hyperprior, fixed = FALSE, constr = TRUE) +
  f(age6, x, model = 'rw1', group = state4,
    control.group = list(model = "iid", hyper = hyperprior),
    hyper = hyperprior, fixed = FALSE, constr = TRUE) +
  # race-state specific intercept and slope
  #f(id_racest, model='iid', hyper = hyperprior, constr = TRUE) + 
  #f(id_racest2, year, model='iid', hyper = hyperprior, constr = TRUE) +
  f(race5, model = 'iid', group = state5,
    control.group = list(model = "iid", hyper = hyperprior),
    hyper = hyperprior, fixed = FALSE, constr = TRUE) +
  f(race6, x, model = 'iid', group = state6,
    control.group = list(model = "iid", hyper = hyperprior),
    hyper = hyperprior, fixed = FALSE, constr = TRUE) +
  # TO DO: age-race-state interaction
  f(id_3way, model='iid', hyper = hyperprior, constr = TRUE) + # state-age-race specific intercept
  f(id_3way2, x, model='iid', hyper = hyperprior, constr = TRUE) +
  f(year2, model="rw1", hyper = hyperprior,  constr = TRUE) + # rw1 over time; scale = TRUE,
  f(e, model = "iid", hyper = hyperprior, constr = TRUE) # overdispersion term

# INLA model rough
mod.rough =
  inla(formula = fml,
       family = "poisson",
       data = subset_test_m,
       E = population,
       control.compute = list(dic=TRUE, openmp.strategy="huge"), # openmp.strategy="pardiso.parallel" ?
       control.predictor = list(link = 1),
       control.inla = list(diagonal=10000, int.strategy='eb',strategy='gaussian'),
       # control.fixed=list(prec=1),
       num.threads = 8,
       verbose=FALSE
  )
# INLA model proper
mod =
  inla(formula = fml,
       family = "poisson",
       data = subset_test_m,
       E = population,
       control.compute = list(config=TRUE, dic=TRUE, openmp.strategy="huge"), # openmp.strategy="pardiso.parallel" ?
       control.predictor = list(link = 1),
       control.inla=list(diagonal=0, cmin = 0, int.strategy = 'eb',strategy='gaussian'),
       control.mode = list(result = mod.rough, restart = FALSE),
       #control.fixed=list(prec=1),
       num.threads = 8,
       verbose=FALSE
  )
plot(mod) 
summary(mod)