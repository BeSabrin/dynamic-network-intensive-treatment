########################################################################
###            Preparatory Analyses for the manuscript               ###
###  Testing implications of the network theory of mental disorders  ###
###               in intensive short-term treatment                  ###
########################################################################

# Author: Sabrina Benz (LMU Munich, Division of Clinical Psychology & 
# Psychological Treatment) & Lea Schumacher (University Medical Center 
# Hamburg-Eppendorf). Large parts of the code are adapted from the preliminary 
# study by Schumacher et al. (2023): https://osf.io/fhqmk/files/kxd7f

# This script fits and evaluates six competing multilevel logistic model
# variants predicting each of six symptoms (mood, loss of interest, decrease
# in activity, nervousness, tension, arousal) from lagged symptom presence:
#   Model 1  — Random intercept only
#   Model 1b — Model 1 + treatment as a fixed interaction effect
#   Model 2  — Model 1 + random slopes for lagged symptom presence
#   Model 2b — Model 2 + treatment as a fixed interaction effect
#   Model 3  — Model 2 + random slopes for the linear change rate
#   Model 3b — Model 3 + treatment as a fixed interaction effect
# Models are fit with brms (Bayesian multilevel models) and compared via
# WAIC/LOO/model weights in the "Model fit evaluation" section below.

# last updated: September 2026

####################################################
####            1. Model Evaluation             ####
####################################################

library("brms")
library("ggplot2")
library("dplyr")
library("tidyr")
library("ggpubr")
library("knitr")
library("shinystan")

# read network data
dat_N <- read.csv2("data.csv", sep= ";")

#### Defining priors ####

# priors of the population intercept is a normal distribution
# which grants relatively high weight to values between 0 and 3 
ggplot(data = data.frame(x = c(-4,4)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 3)) + ylim(0,0.3) 

# priors of the population slopes is a normal distribution 
# which grants relatively high weight to values between -3 and 3
ggplot(data = data.frame(x = c(-4,4)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 3)) + ylim(0,0.3) 

# prior of the sd of population estimates is a HalfCauchy distribution 
# granting higher probabilities to smaller values 
# while also allowing higher (positive) numbers 
ggplot(data = data.frame(x = c(0, 20)), aes(x)) +
  stat_function(fun = dcauchy, args = list( scale = 3)) +   ylim(0,0.15) 


# setting priors 
prior1<- c(
  prior(normal(0,3), class= Intercept), # prior on population intercept 
  prior(normal(0,3), class= b), # prior on population slopes 
  prior(cauchy (0,3), class = sd)) # prior on sd of group-level effects 


#################### Model 1: Random Intercept Model ###########################
# includes random effects for the intercept

## mood 
model1_mood<-brm(
  formula = mood ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1_mood, type = "trace")
mcmc_plot(model1_mood, type = "rhat_hist")
save(model1_mood, file = "Model1_mood.rda")
launch_shinystan(model1_mood) 
summary(model1_mood)

# 3000 iterations
model1_mood_3000 <- update(model1_mood, chains = 4, warmup = 3000,  iter = 9000, 
                           thin=8, cores = 8) # should result in post. sample of 3000 data points
mcmc_plot(model1_mood_3000, type = "trace")
mcmc_plot(model1_mood_3000, type = "rhat_hist")
save(model1_mood_3000, file = "Model1_mood_3000.rda")

## loss_interest 
model1_loss_interest<-brm(
  formula = loss_interest ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1_loss_interest, type = "trace")
mcmc_plot(model1_loss_interest, type = "rhat_hist")
save(model1_loss_interest, file = "Model1_loss_interest.rda")
launch_shinystan(model1_loss_interest) 
summary(model1_loss_interest)


## decrease_activity 
model1_decrease_activity<-brm(
  formula = decrease_activity ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1_decrease_activity, type = "trace")
mcmc_plot(model1_decrease_activity, type = "rhat_hist")
save(model1_decrease_activity, file = "Model1_decrease_activity.rda")
launch_shinystan(model1_decrease_activity) 
summary(model1_decrease_activity)


## nervous 
model1_nervous<-brm(
  formula = nervous ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1_nervous, type = "trace")
mcmc_plot(model1_nervous, type = "rhat_hist")
save(model1_nervous, file = "Model1_nervous.rda")
launch_shinystan(model1_nervous) 
summary(model1_nervous)


## tension 
model1_tension<-brm(
  formula = tension ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1_tension, type = "trace")
mcmc_plot(model1_tension, type = "rhat_hist")
save(model1_tension, file = "Model1_tension.rda")
launch_shinystan(model1_tension) 
summary(model1_tension)


## arousal 
model1_arousal<-brm(
  formula = arousal ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1_arousal, type = "trace")
mcmc_plot(model1_arousal, type = "rhat_hist")
save(model1_arousal, file = "Model1_arousal.rda")
launch_shinystan(model1_arousal) 
summary(model1_arousal)


########### Model 1b: Random Intercept Model with Treatment Effect ############
# make sure treatment effect is a factor
dat_N$treatment <- as.factor(dat_N$treatment)

# treatment effect added as fixed interaction term
## mood 
model1b_mood<-brm(
  formula = mood ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1b_mood, type = "trace")
mcmc_plot(model1b_mood, type = "rhat_hist")
save(model1b_mood, file = "model1b_mood.rda")
# launch_shinystan(model1b_mood) 
summary(model1b_mood)


## loss_interest 
model1b_loss_interest<-brm(
  formula = loss_interest ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1b_loss_interest, type = "trace")
mcmc_plot(model1b_loss_interest, type = "rhat_hist")
save(model1b_loss_interest, file = "model1b_loss_interest.rda")
# launch_shinystan(model1b_loss_interest) 
summary(model1b_loss_interest)


## decrease_activity 
model1b_decrease_activity<-brm(
  formula = decrease_activity ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1b_decrease_activity, type = "trace")
mcmc_plot(model1b_decrease_activity, type = "rhat_hist")
save(model1b_decrease_activity, file = "model1b_decrease_activity.rda")
# launch_shinystan(model1b_decrease_activity) 
summary(model1b_decrease_activity)


## nervous 
model1b_nervous<-brm(
  formula = nervous ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1b_nervous, type = "trace")
mcmc_plot(model1b_nervous, type = "rhat_hist")
save(model1b_nervous, file = "model1b_nervous.rda")
# launch_shinystan(model1b_nervous) 
summary(model1b_nervous)


## tension 
model1b_tension<-brm(
  formula = tension ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1b_tension, type = "trace")
mcmc_plot(model1b_tension, type = "rhat_hist")
save(model1b_tension, file = "model1b_tension.rda")
# launch_shinystan(model1b_tension) 
summary(model1b_tension)


## arousal 
model1b_arousal<-brm(
  formula = arousal ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1b_arousal, type = "trace")
mcmc_plot(model1b_arousal, type = "rhat_hist")
save(model1b_arousal, file = "model1b_arousal.rda")
# launch_shinystan(model1b_arousal) 
summary(model1b_arousal)


########### Model 2: Random Baseline Symptom Association Model #################
# includes random effects for the intercept and lagged-symptom presence

## mood 
model2_mood<-brm(
  formula = mood ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2_mood, type = "trace")
mcmc_plot(model2_mood, type = "rhat_hist")
save(model2_mood, file = "model2_mood.rda")
launch_shinystan(model2_mood) 
summary(model2_mood)


## loss_interest 
model2_loss_interest<-brm(
  formula = loss_interest ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2_loss_interest, type = "trace")
mcmc_plot(model2_loss_interest, type = "rhat_hist")
save(model2_loss_interest, file = "model2_loss_interest.rda")
launch_shinystan(model2_loss_interest) 
summary(model2_loss_interest)


## decrease_activity 
model2_decrease_activity<-brm(
  formula = decrease_activity ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID),  
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2_decrease_activity, type = "trace")
mcmc_plot(model2_decrease_activity, type = "rhat_hist")
save(model2_decrease_activity, file = "model2_decrease_activity.rda")
launch_shinystan(model2_decrease_activity) 
summary(model2_decrease_activity)


## nervous 
model2_nervous<-brm(
  formula = nervous ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2_nervous, type = "trace")
mcmc_plot(model2_nervous, type = "rhat_hist")
save(model2_nervous, file = "model2_nervous.rda")
launch_shinystan(model2_nervous) 
summary(model2_nervous)


## tension 
model2_tension<-brm(
  formula = tension ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2_tension, type = "trace")
mcmc_plot(model2_tension, type = "rhat_hist")
save(model2_tension, file = "model2_tension.rda")
launch_shinystan(model2_tension) 
summary(model2_tension)


## arousal 
model2_arousal<-brm(
  formula = arousal ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2_arousal, type = "trace")
mcmc_plot(model2_arousal, type = "rhat_hist")
save(model2_arousal, file = "model2_arousal.rda")
launch_shinystan(model2_arousal) 
summary(model2_arousal)



####### Model 2b: Random Baseline Association Model with Treatment Effect ######
# includes random effects for the intercept and lagged-symptom presence and 
# treatment effect as fixed interaction term

## mood 
model2b_mood<-brm(
  formula = mood ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2b_mood, type = "trace")
mcmc_plot(model2b_mood, type = "rhat_hist")
save(model2b_mood, file = "model2b_mood.rda")
launch_shinystan(model2b_mood) 
summary(model2b_mood)


## loss_interest 
model2b_loss_interest<-brm(
  formula = loss_interest ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2b_loss_interest, type = "trace")
mcmc_plot(model2b_loss_interest, type = "rhat_hist")
save(model2b_loss_interest, file = "model2b_loss_interest.rda")
launch_shinystan(model2b_loss_interest) 
summary(model2b_loss_interest)


## decrease_activity 
model2b_decrease_activity<-brm(
  formula = decrease_activity ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID),  
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2b_decrease_activity, type = "trace")
mcmc_plot(model2b_decrease_activity, type = "rhat_hist")
save(model2b_decrease_activity, file = "model2b_decrease_activity.rda")
launch_shinystan(model2b_decrease_activity) 
summary(model2b_decrease_activity)


## nervous 
model2b_nervous<-brm(
  formula = nervous ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2b_nervous, type = "trace")
mcmc_plot(model2b_nervous, type = "rhat_hist")
save(model2b_nervous, file = "model2b_nervous.rda")
launch_shinystan(model2b_nervous) 
summary(model2b_nervous)


## tension 
model2b_tension<-brm(
  formula = tension ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2b_tension, type = "trace")
mcmc_plot(model2b_tension, type = "rhat_hist")
save(model2b_tension, file = "model2b_tension.rda")
launch_shinystan(model2b_tension) 
summary(model2b_tension)


## arousal 
model2b_arousal<-brm(
  formula = arousal ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model2b_arousal, type = "trace")
mcmc_plot(model2b_arousal, type = "rhat_hist")
save(model2b_arousal, file = "model2b_arousal.rda")
launch_shinystan(model2b_arousal) 
summary(model2b_arousal)



######### Model 3: Random Baseline & Change Symptom Association Model ##########
# includes random effects for all parameters: intercept, lagged-symptom presence, 
# and interaction with timepoint 

## mood 
model3_mood<-brm(
  formula = mood ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3_mood, type = "trace")
mcmc_plot(model3_mood, type = "rhat_hist")
save(model3_mood, file = "model3_mood.rda")
launch_shinystan(model3_mood) 
summary(model3_mood)


## loss_interest 
model3_loss_interest<-brm(
  formula = loss_interest ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3_loss_interest, type = "trace")
mcmc_plot(model3_loss_interest, type = "rhat_hist")
save(model3_loss_interest, file = "model3_loss_interest.rda")
launch_shinystan(model3_loss_interest) 
summary(model3_loss_interest)


## decrease_activity 
model3_decrease_activity<-brm(
  formula = decrease_activity ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3_decrease_activity, type = "trace")
mcmc_plot(model3_decrease_activity, type = "rhat_hist")
save(model3_decrease_activity, file = "model3_decrease_activity.rda")
launch_shinystan(model3_decrease_activity) 
summary(model3_decrease_activity)


## nervous 
model3_nervous<-brm(
  formula = nervous ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3_nervous, type = "trace")
mcmc_plot(model3_nervous, type = "rhat_hist")
save(model3_nervous, file = "model3_nervous.rda")
launch_shinystan(model3_nervous) 
summary(model3_nervous)


## tension 
model3_tension<-brm(
  formula = tension ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3_tension, type = "trace")
mcmc_plot(model3_tension, type = "rhat_hist")
save(model3_tension, file = "model3_tension.rda")
launch_shinystan(model3_tension) 
summary(model3_tension)


## arousal 
model3_arousal<-brm(
  formula = arousal ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3_arousal, type = "trace")
mcmc_plot(model3_arousal, type = "rhat_hist")
save(model3_arousal, file = "model3_arousal.rda")
launch_shinystan(model3_arousal) 
summary(model3_arousal)


######### Model 3b: Random Baseline & Change Symptom Association Model with Treatment Effect ##########
# includes random effects for all parameters: intercept, lagged-symptom presence, 
# and interaction with timepoint

## mood 
model3b_mood<-brm(
  formula = mood ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3b_mood, type = "trace")
mcmc_plot(model3b_mood, type = "rhat_hist")
save(model3b_mood, file = "model3b_mood.rda")
launch_shinystan(model3b_mood) 
summary(model3b_mood)


## loss_interest 
model3b_loss_interest<-brm(
  formula = loss_interest ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3b_loss_interest, type = "trace")
mcmc_plot(model3b_loss_interest, type = "rhat_hist")
save(model3b_loss_interest, file = "model3b_loss_interest.rda")
launch_shinystan(model3b_loss_interest) 
summary(model3b_loss_interest)


## decrease_activity 
model3b_decrease_activity<-brm(
  formula = decrease_activity ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3b_decrease_activity, type = "trace")
mcmc_plot(model3b_decrease_activity, type = "rhat_hist")
save(model3b_decrease_activity, file = "model3b_decrease_activity.rda")
launch_shinystan(model3b_decrease_activity) 
summary(model3b_decrease_activity)


## nervous 
model3b_nervous<-brm(
  formula = nervous ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3b_nervous, type = "trace")
mcmc_plot(model3b_nervous, type = "rhat_hist")
save(model3b_nervous, file = "model3b_nervous.rda")
launch_shinystan(model3b_nervous) 
summary(model3b_nervous)


## tension 
model3b_tension<-brm(
  formula = tension ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3b_tension, type = "trace")
mcmc_plot(model3b_tension, type = "rhat_hist")
save(model3b_tension, file = "model3b_tension.rda")
launch_shinystan(model3b_tension) 
summary(model3b_tension)


## arousal 
model3b_arousal<-brm(
  formula = arousal ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    mood_lag:week:treatment + loss_interest_lag:week:treatment + decrease_activity_lag:week:treatment + 
    nervous_lag:week:treatment + tension_lag:week:treatment + arousal_lag:week:treatment +
    (1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
       nervous_lag + tension_lag + arousal_lag +
       mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
       nervous_lag:week + tension_lag:week + arousal_lag:week|ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 8000, thin=10, cores = 6,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model3b_arousal, type = "trace")
mcmc_plot(model3b_arousal, type = "rhat_hist")
save(model3b_arousal, file = "model3b_arousal.rda")
launch_shinystan(model3b_arousal) 
summary(model3b_arousal)


######### Save Models ##########

save(model1_mood, model1_loss_interest, model1_decrease_activity, model1_nervous, model1_tension, model1_arousal, file = "models1.rda")
save(model1b_mood, model1b_loss_interest, model1b_decrease_activity, model1b_nervous, model1b_tension, model1b_arousal, file = "models1b.rda")
save(model2_mood, model2_loss_interest, model2_decrease_activity, model2_nervous, model2_tension, model2_arousal, file = "models2.rda")
save(model2b_mood, model2b_loss_interest, model2b_decrease_activity, model2b_nervous, model2b_tension, model2b_arousal, file = "models2b.rda")
save(model3_mood, model3_loss_interest, model3_decrease_activity, model3_nervous, model3_tension, model3_arousal, file = "models3.rda")
save(model3b_mood, model3b_loss_interest, model3b_decrease_activity, model3b_nervous, model3b_tension, model3b_arousal, file = "models3b.rda")

####################################################
####            2. Model Fit Evaluation         ####
####################################################

library("Rcpp")
library("brms")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggmcmc")

# load models
load("models1.rda")
load("models1b.rda")
load("models2.rda")
load("models2b.rda")
load("models3.rda")
load("models3b.rda")

## model 1
waic(model1_mood)
loo(model1_mood)
bayes_R2(model1_mood)

waic(model1_loss_interest)
loo(model1_loss_interest, cores = 2)

waic(model1_decrease_activity)
loo(model1_decrease_activity, cores = 2, moment_match = TRUE)
bayes_R2(model1_decrease_activity)

waic(model1_nervous)
loo(model1_nervous, cores = 2)

waic(model1_tension)
loo(model1_tension, cores = 2)

waic(model1_arousal)
loo(model1_arousal, cores = 2, moment_match = TRUE)

## model 1b
waic(model1b_mood)
loo(model1b_mood,cores = 2)

waic(model1b_loss_interest)
loo(model1b_loss_interest, cores = 2, moment_match = TRUE)

waic(model1b_decrease_activity)
loo(model1b_decrease_activity, cores = 2, moment_match = TRUE)

waic(model1b_nervous)
loo(model1b_nervous, cores = 2)

waic(model1b_tension)
loo(model1b_tension, cores = 2)

waic(model1b_arousal)
loo(model1b_arousal, cores = 2, moment_match = TRUE)
bayes_R2(model1b_arousal)

## model 2
waic(model2_mood)
loo(model2_mood,cores = 2, reloo = TRUE)

waic(model2_loss_interest)
loo(model2_loss_interest, cores = 2, moment_match = TRUE) 
kfold(model2_loss_interest, k = 10)

waic(model2_decrease_activity)
loo(model2_decrease_activity, cores = 2, reloo = TRUE)

waic(model2_nervous)
loo(model2_nervous, cores = 2) 

waic(model2_tension)
loo(model2_tension, cores = 6, reloo = TRUE) 

waic(model2_arousal)
loo(model2_arousal, cores = 2, moment_match = TRUE) 
kfold(model2_arousal, k = 10)
  
## model 2b
waic(model2b_mood)
loo(model2b_mood,cores = 2, reloo = TRUE) 
#bayes_R2(model2b_mood)

waic(model2b_loss_interest)
loo(model2b_loss_interest, cores = 8)
#bayes_R2(model2b_loss_interest)

waic(model2b_decrease_activity)
loo(model2b_decrease_activity, cores = 8)
#bayes_R2(model2b_decrease_activity)

waic(model2b_nervous)
loo(model2b_nervous, cores = 8, moment_match = TRUE)
kfold(model2b_nervous, k = 10)
#bayes_R2(model2b_nervous)

waic(model2b_tension)
loo(model2b_tension, cores = 2, reloo = TRUE) 
#bayes_R2(model2b_tension)

waic(model2b_arousal)
loo(model2b_arousal, cores = 2, moment_match = TRUE)
#bayes_R2(model2b_arousal)

## model 3
waic(model3_mood)
loo(model3_mood,cores = 2)
#bayes_R2(model3_mood)

waic(model3_loss_interest)
loo(model3_loss_interest, cores = 2)
#bayes_R2(model3_loss_interest)

waic(model3_decrease_activity)
loo(model3_decrease_activity, cores = 2)
#bayes_R2(model3_decrease_activity)

waic(model3_nervous)
loo(model3_nervous, cores = 2, moment_match = TRUE)
#bayes_R2(model3_nervous)

waic(model3_tension)
loo(model3_tension, cores = 2, moment_match = TRUE)
bayes_R2(model3_tension)

waic(model3_arousal)
loo(model3_arousal, cores = 2)
#bayes_R2(model3_arousal)

## model 3b
waic(model3b_mood)
loo(model3b_mood, cores = 2)
#bayes_R2(model3b_mood)

summary(model3b_loss_interest)
waic(model3b_loss_interest)
loo(model3b_loss_interest, cores = 2) 
#bayes_R2(model3b_loss_interest)

waic(model3b_decrease_activity)
loo(model3b_decrease_activity, cores = 2) 
#bayes_R2(model3b_decrease_activity)

waic(model3b_nervous)
loo(model3b_nervous, cores = 2) 
#bayes_R2(model3b_nervous)

waic(model3b_tension)
loo(model3b_tension, cores = 2)
#bayes_R2(model3b_tension)

waic(model3b_arousal)
loo(model3b_arousal, cores = 2)
#bayes_R2(model3b_arousal)

### model comparison 
# mood 
model_weights(model1_mood, model1b_mood, model2_mood, model2b_mood, model3_mood, model3b_mood, weights = "waic")
model_weights(model1_mood, model1b_mood, model2_mood, model2b_mood, model3_mood, model3b_mood, weights = "loo")

# loss_interest
model_weights(model1_loss_interest, model1b_loss_interest, model2_loss_interest, 
              model2b_loss_interest, model3_loss_interest, model3b_loss_interest, weights = "waic")
model_weights(model1_loss_interest, model1b_loss_interest, model2_loss_interest, 
              model2b_loss_interest, model3_loss_interest, model3b_loss_interest, weights = "loo")

# decrease_activity
model_weights(model1_decrease_activity, model1b_decrease_activity, model2_decrease_activity, 
              model2b_decrease_activity, model3_decrease_activity, model3b_decrease_activity, weights = "waic")
model_weights(model1_decrease_activity, model1b_decrease_activity, model2_decrease_activity, 
              model2b_decrease_activity, model3_decrease_activity, model3b_decrease_activity, weights = "loo")

# nervous
model_weights(model1_nervous, model1b_nervous, model2_nervous, model2b_nervous, 
              model3_nervous, model3b_nervous, weights = "waic")
model_weights(model1_nervous, model1b_nervous, model2_nervous, model2b_nervous, 
              model3_nervous, model3b_nervous, weights = "loo")

# tension
model_weights(model1_tension, model1b_tension, model2_tension, model2b_tension, 
              model3_tension, model3b_tension, weights = "waic")
model_weights(model1_tension, model1b_tension, model2_tension, model2b_tension, 
              model3_tension, model3b_tension, weights = "loo")

# arousal
model_weights(model1_arousal, model1b_arousal, model2_arousal, model2b_arousal, 
              model3_arousal, model3b_arousal, weights = "waic")
model_weights(model1_arousal, model1b_arousal, model2_arousal, model2b_arousal, 
              model3_arousal, model3b_arousal, weights = "loo")
