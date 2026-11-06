########################################################################
###                Main Analysis for the manuscript                  ###
###  Testing implications of the network theory of mental disorders  ###
###               in intensive short-term treatment                  ###
########################################################################

# Author: Sabrina Benz (LMU Munich, Division of Clinical Psychology & 
# Psychological Treatment) & Lea Schumacher (University Medical Center 
# Hamburg-Eppendorf). Large parts of the code are adapted from the preliminary 
# study by Schumacher et al. (2023): https://osf.io/fhqmk/files/rdcv7

# This script fits the final model (Model 1: random intercept only, selected 
# via WAIC/LOO comparison in the preparatory analyses script) for each of six 
# symptoms, then uses the posterior distributions to construct and evaluate 
# time-varying symptom networks across 8 weekly timepoints of inpatient 
# treatment. It contains five parts:

#   1. Model estimation      — refit Model 1 with 3000 post-warmup draws 
#                               (4 chains, thin = 8) and check convergence

#   2. Results inspection    — fixed effects as odds ratios; random-intercept 
#                               variance components

#   3. Model fit evaluation  — Bayesian R^2; model-implied vs. observed 
#                               symptom prevalence over time

#   4. Network analysis      — cross-lagged edges reconstructed from the full 
#                               posterior per week; node- and network-level 
#                               parameters (in-/outstrength, expected 
#                               influence, global connectivity)

#   5. Plotting              — network graphs (all 8 weeks and key 
#                               timepoints), node-level and global 
#                               network-parameter trajectories over time

#   6. Sensitivity Analysis  - adding time as main effect: estimations, inspection,
#                               network construction & plotting

# Assumes the preparatory analyses script has already been run to select 
# Model 1.

# last updated: September 2026


################################
###### 1. Model estimation #####
################################

# load best fitting model (model 1) with 3000 iterations
load("models1_3000.rda")

library("brms")
library("ggplot2")
library("ggmcmc")
library("dplyr")

### Defining priors ###

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
  stat_function(fun = dcauchy, args = list( scale = 3)) +
  ylim(0,0.15) 

# setting priors 
prior1<- c(
  prior(normal(0,3), class= Intercept), # prior on population intercept 
  prior(normal(0,3), class= b), # prior on population slopes 
  prior(cauchy (0,3), class = sd)) # prior on sd of group-level effects


#### Model estimation ###
# Model 1: Random Intercept Only Model

## mood
model1_mood_3000 <-brm(
  formula = mood ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

save(model1_mood_3000, file = "model1_mood_3000.rda")
summary(model1_mood_3000)

# check convergence 
# load("model1_mood_3000.rda")
launch_shinystan(model1_mood_3000) 
mcmc_plot(model1_mood_3000, type = "rhat_hist")

# make traceplots including warm up phase 
posterior <- ggs(model1_mood_3000)
par <- variables(model1_mood_3000)[1:14] # Intercept, 6x lag, 6x change rate, sd(Intercept)
trace_list<- list()

for(i in par){
  posterior1<-posterior[which(posterior$Parameter == i),]
  posterior1$Chain<- as.factor(posterior1$Chain)
  
  trace_list[[i]]<-ggplot(data= posterior1, aes(x = Iteration, y = value, color = Chain)) +
    geom_rect(aes(xmin=0,xmax=375,ymin=-Inf,ymax=Inf),
              fill="white")+
    geom_line() +
    theme(panel.grid = element_blank())+
    scale_x_continuous(limits = c(0, 1125)) +
    ggtitle(i)
}

pdf("trace_plot_mood.pdf", width = 5.8 , height = 4.5)
for (i in par){
  plot(trace_list[[i]])}
dev.off()

mcmc_plot(model1_mood_3000, type = "trace")

## loss_interest 
model1_loss_interest_3000 <-brm(
  formula = loss_interest ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

save(model1_loss_interest_3000, file = "model1_loss_interest_3000.rda")
summary(model1_loss_interest_3000)

# check convergence
# load("model1_loss_interest_3000.rda")
launch_shinystan(model1_loss_interest_3000) 
mcmc_plot(model1_loss_interest_3000, type = "rhat_hist")

# make traceplots including warm up phase 
posterior<- ggs(model1_loss_interest_3000)
par<- variables(model1_loss_interest_3000)[1:14]
trace_list<- list()

for(i in par){
  posterior1<-posterior[which(posterior$Parameter == i),]
  posterior1$Chain<- as.factor(posterior1$Chain)
  
  trace_list[[i]]<-ggplot(data= posterior1, aes(x = Iteration, y = value, color = Chain)) +
    geom_rect(aes(xmin=0,xmax=375,ymin=-Inf,ymax=Inf),
              fill="white")+
    geom_line() +
    theme(panel.grid = element_blank())+
    scale_x_continuous(limits = c(0, 1125)) +
    ggtitle(i)
}

pdf("trace_plot_loss_interest.pdf", width = 5.8 , height = 4.5)
for (i in par){
  plot(trace_list[[i]])}
dev.off()

mcmc_plot(model1_loss_interest_3000, type = "trace")

## decrease_activity 
model1_decrease_activity_3000 <-brm(
  formula = decrease_activity ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

save(model1_decrease_activity_3000, file = "model1_decrease_activity_3000.rda")
summary(model1_decrease_activity_3000)

# check convergence
# load("model1_decrease_activity_3000.rda")
launch_shinystan(model1_decrease_activity_3000) 
mcmc_plot(model1_decrease_activity_3000, type = "rhat_hist")

# make traceplots including warm up phase 
posterior<- ggs(model1_decrease_activity_3000)
par<- variables(model1_decrease_activity_3000)[1:14]
trace_list<- list()

for(i in par){
  posterior1<-posterior[which(posterior$Parameter == i),]
  posterior1$Chain<- as.factor(posterior1$Chain)
  
  trace_list[[i]]<-ggplot(data= posterior1, aes(x = Iteration, y = value, color = Chain)) +
    geom_rect(aes(xmin=0,xmax=375,ymin=-Inf,ymax=Inf),
              fill="white")+
    geom_line() +
    theme(panel.grid = element_blank())+
    scale_x_continuous(limits = c(0, 1125)) +
    ggtitle(i)
}

pdf("trace_plot_decrease_activity.pdf", width = 5.8 , height = 4.5)
for (i in par){
  plot(trace_list[[i]])}
dev.off()

mcmc_plot(model1_decrease_activity_3000, type = "trace")

## nervous 
model1_nervous_3000 <-brm(
  formula = nervous ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

save(model1_nervous_3000, file = "model1_nervous_3000.rda")
summary(model1_nervous_3000)

# check convergence
# load("model1_nervous_3000.rda")
launch_shinystan(model1_nervous_3000) 
mcmc_plot(model1_nervous_3000, type = "rhat_hist")

# make traceplots including warm up phase 
posterior<- ggs(model1_nervous_3000)
par<- variables(model1_nervous_3000)[1:14]
trace_list<- list()

for(i in par){
  posterior1<-posterior[which(posterior$Parameter == i),]
  posterior1$Chain<- as.factor(posterior1$Chain)
  
  trace_list[[i]]<-ggplot(data= posterior1, aes(x = Iteration, y = value, color = Chain)) +
    geom_rect(aes(xmin=0,xmax=375,ymin=-Inf,ymax=Inf),
              fill="white")+
    geom_line() +
    theme(panel.grid = element_blank())+
    scale_x_continuous(limits = c(0, 1125)) +
    ggtitle(i)
}

pdf("trace_plot_nervous.pdf", width = 5.8 , height = 4.5)
for (i in par){
  plot(trace_list[[i]])}
dev.off()

mcmc_plot(model1_nervous_3000, type = "trace")

## tension 
model1_tension_3000 <-brm(
  formula = tension ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

save(model1_tension_3000, file = "model1_tension_3000.rda")
summary(model1_tension_3000)

# check convergence
# load("model1_tension_3000.rda")
launch_shinystan(model1_tension_3000) 
mcmc_plot(model1_tension_3000, type = "rhat_hist")

# make traceplots including warm up phase 
posterior<- ggs(model1_tension_3000)
par<- variables(model1_tension_3000)[1:14]
trace_list<- list()

for(i in par){
  posterior1<-posterior[which(posterior$Parameter == i),]
  posterior1$Chain<- as.factor(posterior1$Chain)
  
  trace_list[[i]]<-ggplot(data= posterior1, aes(x = Iteration, y = value, color = Chain)) +
    geom_rect(aes(xmin=0,xmax=375,ymin=-Inf,ymax=Inf),
              fill="white")+
    geom_line() +
    theme(panel.grid = element_blank())+
    scale_x_continuous(limits = c(0, 1125)) +
    ggtitle(i)
}

pdf("trace_plot_tension.pdf", width = 5.8 , height = 4.5)
for (i in par){
  plot(trace_list[[i]])}
dev.off()

mcmc_plot(model1_tension_3000, type = "trace")

## arousal 
model1_arousal_3000 <-brm(
  formula = arousal ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

save(model1_arousal_3000, file = "model1_arousal_3000.rda")
summary(model1_arousal_3000)

# check convergence
# load("model1_arousal_3000.rda")
launch_shinystan(model1_arousal_3000) 
mcmc_plot(model1_arousal_3000, type = "rhat_hist")

# make traceplots including warm up phase 
posterior<- ggs(model1_arousal_3000)
par<- variables(model1_arousal_3000)[1:14]
trace_list<- list()

for(i in par){
  posterior1<-posterior[which(posterior$Parameter == i),]
  posterior1$Chain<- as.factor(posterior1$Chain)
  
  trace_list[[i]]<-ggplot(data= posterior1, aes(x = Iteration, y = value, color = Chain)) +
    geom_rect(aes(xmin=0,xmax=375,ymin=-Inf,ymax=Inf),
              fill="white")+
    geom_line() +
    theme(panel.grid = element_blank())+
    scale_x_continuous(limits = c(0, 1125)) +
    ggtitle(i)
}

pdf("trace_plot_arousal.pdf", width = 5.8 , height = 4.5)
for (i in par){
  plot(trace_list[[i]])}
dev.off()

mcmc_plot(model1_arousal_3000, type = "trace")

# save all models in one file
save(model1_mood_3000, model1_loss_interest_3000, model1_decrease_activity_3000, 
     model1_nervous_3000, model1_tension_3000, model1_arousal_3000, 
     file = "models1_3000.rda")



##################################
###### 2. Results inspection #####
##################################

#### create large table with model parameters of all models ####

library("brms")
library("ggplot2")
library("openxlsx")

load("models1_3000.rda")

## extract fixed effects and convert to odds ratios ##
# first 13 variables of each model include fixed effects (1 Intercept, 6 slopes, 6 change slopes)

extract_OR_table <- function(model, suffix, n_fixed = 13) {
  fit <- as.data.frame(posterior_summary(model)[1:n_fixed, ])
  fit$OR   <- exp(fit$Estimate)
  fit$CI_l <- exp(fit$Q2.5)
  fit$CI_h <- exp(fit$Q97.5)
  fit <- round(fit, 2)
  fit$CI    <- paste(fit$CI_l, fit$CI_h, sep = "-")
  fit$dif_1 <- ifelse(fit$CI_l > 1 | fit$CI_h < 1, "*", "")
  out <- fit[, c("OR", "CI", "dif_1")]
  colnames(out) <- paste(c("OR", "CI", "dif_1"), suffix, sep = "_")
  out
}

# one model per symptom, with the suffix used in the original column names
models_list <- list(
  m = model1_mood_3000,
  i = model1_loss_interest_3000,
  d = model1_decrease_activity_3000,
  n = model1_nervous_3000,
  t = model1_tension_3000,
  a = model1_arousal_3000
)

final_table <- do.call(cbind, Map(extract_OR_table, models_list, names(models_list)))

# save table to csv file
write.csv2(final_table, "final_table.csv")


## extract random effect variance
# please note that the only random parameter in the final model is the intercept
options(digits = 4)

# mood
re_mood <- summary(model1_mood_3000)[["random"]][["ID"]][, 1:4]
re_mood$CI <-paste(round(re_mood$`l-95% CI`, 3),
                   round(re_mood$`u-95% CI`, 3), sep = "-")
re_mood$Estimate <-round(re_mood$Estimate, 3)
write.csv2(re_mood, "random_effects_mood.csv")

# loss interest
re_loss_interest <-summary(model1_loss_interest_3000)[["random"]][["ID"]][, 1:4]
re_loss_interest$CI <-paste(round(re_loss_interest$`l-95% CI`, 3),
                   round(re_loss_interest$`u-95% CI`, 3), sep = "-")
re_loss_interest$Estimate <-round(re_loss_interest$Estimate, 3)
write.csv2(re_loss_interest, "random_effects_loss_interest.csv")

# decrease activity
re_decrease_activity <-summary(model1_decrease_activity_3000)[["random"]][["ID"]][,1:4]
re_decrease_activity$CI <-paste(round(re_decrease_activity$`l-95% CI`, 3),
                      round(re_decrease_activity$`u-95% CI`, 3), sep = "-")
re_decrease_activity$Estimate <-round(re_decrease_activity$Estimate, 3)
write.csv2(re_decrease_activity, "random_effects_decrease_activity.csv")

# nervous
re_nervous <-summary(model1_nervous_3000)[["random"]][["ID"]][, 1:4]
re_nervous$CI <-paste(round(re_nervous$`l-95% CI`, 3),
                      round(re_nervous$`u-95% CI`, 3), sep = "-")
re_nervous$Estimate <-round(re_nervous$Estimate, 3)
write.csv2(re_nervous, "random_effects_nervous.csv")

# tension
re_tension <-summary(model1_tension_3000)[["random"]][["ID"]][, 1:4]
re_tension$CI <-paste(round(re_tension$`l-95% CI`, 3),
                  round(re_tension$`u-95% CI`, 3), sep = "-")
re_tension$Estimate <-round(re_tension$Estimate, 3)
write.csv2(re_tension, "random_effects_tension.csv")

# arousal
re_arousal <-summary(model1_arousal_3000)[["random"]][["ID"]][, 1:4]
re_arousal$CI <-paste(round(re_arousal$`l-95% CI`, 3),
                         round(re_arousal$`u-95% CI`, 3), sep = "-")
re_arousal$Estimate <-round(re_arousal$Estimate, 3)
write.csv2(re_arousal, "random_effects_arousal.csv")


###################################
###### 3. Model fit evaluation ####
###################################

library("dplyr")
library("brms")
library("dvmisc")
library("ggplot2")
library("ggpubr")
library("jtools")

# load data and models
load("models1_3000.rda")

# load data

#### calculate Bayesian pseudo R2 
bayes_R2(model1_mood_3000)
bayes_R2(model1_loss_interest_3000)
bayes_R2(model1_decrease_activity_3000)
bayes_R2(model1_nervous_3000)
bayes_R2(model1_tension_3000)
bayes_R2(model1_arousal_3000)

#### compare the estimated and observed symptom prevalence 
symptoms<- c("Depressive mood", "Loss of interest", "Decrease in activity", 
             "Nervousness", "Physical tension", "Autonomic arousal")

# parameter list with fixed effects
par_list <-list(
 as.data.frame(posterior_summary(model1_mood_3000)[1:13,]),
 as.data.frame(posterior_summary(model1_loss_interest_3000)[1:13,]),
 as.data.frame(posterior_summary(model1_decrease_activity_3000)[1:13,]),
 as.data.frame(posterior_summary(model1_nervous_3000)[1:13,]),
 as.data.frame(posterior_summary(model1_tension_3000)[1:13,]),
 as.data.frame(posterior_summary(model1_arousal_3000)[1:13,]))

names(par_list)<- symptoms

# create function for transformation from log-ratios to probabilities
# (mathematically identical to plogis())
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# make a dataframe for all estimated symptom probabilities
symp <- as.data.frame(matrix(nrow=8, ncol = 7)) # nrow = no. of measurement points, ncol = no of symptoms + 1
colnames(symp) <- c("week", symptoms)
symp$week <- 0:7

# prevalence for the 0. session is estimated by taking the log odds at session 0
# (log odds = log(P(symptom present)/P(symptom not present))) -> log(P(symptom)/1-P(symptom))
# P(symptom) = mean(symptom) -> binary coded
dat_N <- as.data.frame(dat_N)

symp[1,]<- c(0,c(
              (log(mean(dat_N[dat_N$week == 0, "mood"], na.rm = TRUE)/
              (1-mean(dat_N[dat_N$week == 0, "mood"], na.rm = TRUE)))),
              
              (log(mean(dat_N[dat_N$week == 0, "loss_interest"], na.rm = TRUE)/
              (1-mean(dat_N[dat_N$week == 0, "loss_interest"], na.rm = TRUE)))),
              
              (log(mean(dat_N[dat_N$week == 0, "decrease_activity"], na.rm = TRUE)/
              (1-mean(dat_N[dat_N$week == 0, "decrease_activity"], na.rm = TRUE)))),
              
              (log(mean(dat_N[dat_N$week == 0, "nervous"], na.rm = TRUE)/
              (1-mean(dat_N[dat_N$week == 0, "nervous"], na.rm = TRUE)))),
              
              (log(mean(dat_N[dat_N$week == 0, "tension"], na.rm = TRUE)/
              (1-mean(dat_N[dat_N$week == 0, "tension"], na.rm = TRUE)))),
              
              (log(mean(dat_N[dat_N$week == 0, "arousal"], na.rm = TRUE)/
              (1-mean(dat_N[dat_N$week == 0, "arousal"], na.rm = TRUE))))
              ))

symp[1,-1]<- logit2prob(symp[1,-1])

# calculate the probability of each symptom during each week of therapy
for(i in 1:7){
  for(k in symptoms){
    symp[which(symp$week == i), k] <- 
    logit2prob(
      par_list[[k]]["b_Intercept", "Estimate"]  +
      par_list[[k]]["b_mood_lag", "Estimate"]*
        symp[which(symp$week == i-1),"Depressive mood"]+
      par_list[[k]]["b_loss_interest_lag", "Estimate"]*
        symp[which(symp$week == i-1),"Loss of interest"]+
      par_list[[k]]["b_decrease_activity_lag", "Estimate"]*
        symp[which(symp$week == i-1),"Decrease in activity"] + 
      par_list[[k]]["b_nervous_lag", "Estimate"]*
        symp[which(symp$week == i-1),"Nervousness"] +
      par_list[[k]]["b_tension_lag", "Estimate"]*
        symp[which(symp$week == i-1),"Physical tension"] +
      par_list[[k]]["b_arousal_lag", "Estimate"]*
        symp[which(symp$week == i-1),"Autonomic arousal"]+
        
      par_list[[k]]["b_mood_lag:week", "Estimate"]* i *
        symp[which(symp$week == i-1),"Depressive mood"]+
      par_list[[k]]["b_loss_interest_lag:week", "Estimate"]* i *
        symp[which(symp$week == i-1),"Loss of interest"]+
      par_list[[k]]["b_decrease_activity_lag:week", "Estimate"]* i *
        symp[which(symp$week == i-1),"Decrease in activity"] + 
      par_list[[k]]["b_nervous_lag:week", "Estimate"]* i *
        symp[which(symp$week == i-1),"Nervousness"] +
      par_list[[k]]["b_tension_lag:week", "Estimate"]* i *
        symp[which(symp$week == i-1),"Physical tension"] +
      par_list[[k]]["b_arousal_lag:week", "Estimate"]* i *
        symp[which(symp$week == i-1),"Autonomic arousal"]
      )
  }
}

# check
symp

## calculate the observed relative frequencies 
obs_symp <- as.data.frame(matrix(nrow=8, ncol = 7)) # nrow = no of time points; ncol= no of symps + 1
colnames(obs_symp) <- c("week", symptoms)
obs_symp$week <- 0:7

for(i in 0:7){
  obs_symp[i+1,2:7]<-c(mean(dat_N[dat_N$week == i, "mood"], na.rm = T),
                          mean(dat_N[dat_N$week == i, "loss_interest"], na.rm = T),
                          mean(dat_N[dat_N$week == i, "decrease_activity"], na.rm = T),
                          mean(dat_N[dat_N$week == i, "nervous"], na.rm = T),
                          mean(dat_N[dat_N$week == i, "tension"], na.rm = T),
                          mean(dat_N[dat_N$week == i, "arousal"], na.rm = T))
  }


colnames(symp)[2:7]<- gsub(" ", "_",x= paste(symptoms, "est"))
colnames(obs_symp)[2:7]<-  gsub(" ", "_",x= paste(symptoms, "obs"))

all_symp <- cbind(symp, obs_symp[,-1])

save(all_symp, file = "observedpredictedfrequencies.rda")


# plot the observed and predicted frequencies over time
colnames(all_symp)

obs_pred_Symp <- ggarrange(
# mood
ggplot(all_symp) + 
  geom_point(aes(x = week , y= Depressive_mood_obs))+
  geom_smooth(aes(x = week , y= Depressive_mood_obs), se=F)+
  
  geom_point(aes(x = week , y= Depressive_mood_est), shape = 1)+
  geom_line(aes(x = week , y= Depressive_mood_est), linetype ="dashed")+

  labs(title = "Depressive mood", x = "week", y= "Observed frequency/Predicted probability")+

  scale_x_continuous(breaks=c(seq(0,7,1)))+
  scale_y_continuous(limits = c(0,1))+
  theme_apa(),

# loss interest
ggplot(all_symp)+ 
  geom_point(aes(x = week , y= Loss_of_interest_obs))+
  geom_smooth(aes(x = week , y= Loss_of_interest_obs), se=F)+
  
  geom_point(aes(x = week , y= Loss_of_interest_est), shape = 1)+
  geom_line(aes(x = week , y= Loss_of_interest_est), linetype ="dashed")+
  
  labs(title = "Loss of interest", x = "week", y= "Observed frequency/Predicted probability")+
  
  scale_x_continuous(breaks=c(seq(0,7,1)))+
  scale_y_continuous(limits = c(0,1))+
  theme_apa(),

# decrease interest
ggplot(all_symp)+ 
  geom_point(aes(x = week , y= Decrease_in_activity_obs))+
  geom_smooth(aes(x = week , y= Decrease_in_activity_obs), se=F)+
  
  geom_point(aes(x = week , y= Decrease_in_activity_est), shape = 1)+
  geom_line(aes(x = week , y= Decrease_in_activity_est),linetype ="dashed")+

  labs(title = "Decrease in activity", x = "week", y= "Observed frequency/Predicted probability")+
  
  scale_x_continuous(breaks=c(seq(0,7,1)))+
  scale_y_continuous(limits = c(0,1))+
  theme_apa(),

# Nervousness
ggplot(all_symp)+ 
  geom_point(aes(x = week , y= Nervousness_obs))+
  geom_smooth(aes(x = week , y=Nervousness_obs), se=F)+
  
  geom_point(aes(x = week , y=Nervousness_est), shape = 1)+
  geom_line(aes(x = week , y=Nervousness_est), linetype ="dashed")+
  
  labs(title = "Nervousness", x = "week", y= "Observed frequency/Predicted probability")+
  
  scale_x_continuous(breaks=c(seq(0,7,1)))+
  scale_y_continuous(limits = c(0,1))+
  theme_apa(),

# Physical tension
ggplot(all_symp)+ 
  geom_point(aes(x = week , y= Physical_tension_obs))+
  geom_smooth(aes(x = week , y= Physical_tension_obs), se=F)+
  
  geom_point(aes(x = week , y= Physical_tension_est), shape = 1)+
  geom_line(aes(x = week , y= Physical_tension_est), linetype ="dashed")+

  labs(title = "Physical tension", x = "week", y= "Observed frequency/Predicted probability")+
  
  scale_x_continuous(breaks=c(seq(0,7,1)))+
  scale_y_continuous(limits = c(0,1))+
  theme_apa(),

# Autonomic arousal
ggplot(all_symp)+ 
  geom_point(aes(x = week , y= Autonomic_arousal_obs))+
  geom_smooth(aes(x = week , y= Autonomic_arousal_obs), se=F)+
  
  geom_point(aes(x = week , y=Autonomic_arousal_est), shape = 1)+
  geom_line(aes(x = week , y=Autonomic_arousal_est), linetype ="dashed")+
  
  labs(title = "Autonomic arousal", x = "week", y= "Observed frequency/Predicted probability")+
  
  scale_x_continuous(breaks=c(seq(0,7,1)))+
  scale_y_continuous(limits = c(0,1))+
  theme_apa(), 

  nrow= 2, ncol = 3)

pdf("observed_and_predicted_frequencies.pdf", width = 12 , height = 10) 
obs_pred_Symp
dev.off()


################################
###### 4. Network Analysis #####
################################

library("qgraph")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("brms")

# load models 
load("models1_3000.rda")

symptoms<- c("Depressive mood", "Loss of interest", "Decrease in activity", 
             "Nervousness", "Physical tension", "Autonomic arousal")

week<- as.character(0:7)
lines<- 1:3000

################# Edge calculation #############################################
# to be able to form a CI for each edge, edges are calculated from the whole 
# posterior distribution 

## collect all posterior samples in a list  
post_samp <- list(as_draws_df(model1_mood_3000)[,1:13],
                as_draws_df(model1_loss_interest_3000)[,1:13],
                as_draws_df(model1_decrease_activity_3000)[,1:13],
                as_draws_df(model1_nervous_3000)[,1:13],
                as_draws_df(model1_tension_3000)[,1:13],
                as_draws_df(model1_arousal_3000)[,1:13])

names(post_samp) <- symptoms

save(post_samp, file = "posterior_coef_sample.RData")

## calculate edges for every session 

# create 
# - matrices for the incoming edges of each symptom for the whole posterior 
#   distribution (3000 x 6)
# - a list with 6 matrices for each symptom 
# - a list of 8 lists for each week which includes 6 matrices for each symptom
#   including the whole posterior distribution 

edges <- matrix(nrow = 3000, ncol = 6)
colnames(edges) <- symptoms
symp_edges <- replicate(6, edges, simplify=FALSE)
names(symp_edges) <- symptoms
time_edges <- replicate(8,symp_edges, simplify=FALSE) 
names(time_edges) <- week

for(k in 0:7){
  for (i in symptoms){
    time_edges[[as.character(k)]][[i]] <- as.data.frame(cbind(
      
      post_samp[[i]][,"b_mood_lag"] + post_samp[[i]][,"b_mood_lag:week"]*k, 
        
      post_samp[[i]][,"b_loss_interest_lag"] + post_samp[[i]][,"b_loss_interest_lag:week"]*k,
      
      post_samp[[i]][,"b_decrease_activity_lag"] + post_samp[[i]][,"b_decrease_activity_lag:week"]*k,
      
      post_samp[[i]][,"b_nervous_lag"] + post_samp[[i]][,"b_nervous_lag:week"]*k,
        
      post_samp[[i]][,"b_tension_lag"] + post_samp[[i]][,"b_tension_lag:week"]*k,
      
      post_samp[[i]][,"b_arousal_lag"] + post_samp[[i]][,"b_arousal_lag:week"]*k))
    
      colnames(time_edges[[as.character(k)]][[i]]) <- paste(symptoms,i, sep = "-")
  }
}


## construction of edge matrices
## create an edge matrix for the complete posterior distribution of edges 
edge_m <- as.data.frame(matrix(nrow = 6, ncol = 6))
colnames(edge_m) <- symptoms
post_edge_m <- replicate(3000,edge_m, simplify=FALSE)
time_post_edge_m <- replicate(8,post_edge_m, simplify=FALSE)
names(time_post_edge_m) <- week

for(i in week){
  for(l in lines){
     for(k in symptoms){
      time_post_edge_m[[i]][[l]][,k]<- as.numeric(
                                    time_edges[[i]][[k]][l,])
    }
  }
}

save(time_post_edge_m, file = "edge_matrices.rda")

## create an edge matrix for the posterior mean edges and its upper and 
## lower limit of the 95% credibility interval
mean_edge_m <- as.data.frame(matrix(nrow = 6, ncol = 6))
colnames(mean_edge_m)<- symptoms
time_mean_edge_m<- replicate(8,mean_edge_m, simplify=FALSE)
names(time_mean_edge_m)<- week
time_cih_edge_m <- time_cil_edge_m <- time_mean_edge_m # CI upper and lower limit

for(i in week){
    for(k in symptoms){
      for(l in 1:6){ 
        time_mean_edge_m[[i]][l,k] <- mean(
        time_edges[[i]][[k]][,l])
        
        time_cil_edge_m[[i]][l,k] <- sort(time_edges[[i]][[k]][,l])[76] # for lower 2.5%
        time_cih_edge_m[[i]][l,k] <- sort(time_edges[[i]][[k]][,l])[2926] # for upper 2.5%
    }
  }
}

save(time_mean_edge_m, file = "mean_edge_matrix.rda")
save(time_cil_edge_m, file = "cil_edge_matrix.rda")
save(time_cih_edge_m, file = "cih_edge_matrix.rda")

write.xlsx(time_mean_edge_m, file = "mean_edge_matrix.xlsx")
write.xlsx(time_cil_edge_m, file = "cil_edge_matrix.xlsx")
write.xlsx(time_cih_edge_m, file = "cih_edge_matrix.xlsx")


## Edge parameters across timepoints
# create one table with edge weight and 95% cred interval per timepoint

library(dplyr)
library(openxlsx) 

symptoms <- c("mood", "loss of interest", "decrease in activity", 
              "nervousness", "physical tension", "autonomic arousal")

final_tables <- list()

for (t in 1:8) {
  est <- time_mean_edge_m[[t]]
  ci_l <- time_cil_edge_m[[t]]
  ci_h <- time_cih_edge_m[[t]]
  
  colnames(est) <- colnames(ci_l) <- colnames(ci_h) <- symptoms
  
  est <- est %>% mutate(across(all_of(symptoms), ~round(.x,2)))
  colnames(est) <- paste0(colnames(est), "_est")
  
  ci_l <- ci_l %>% mutate(across(all_of(symptoms), ~round(.x,2)))
  colnames(ci_l) <- paste0(colnames(ci_l), "_l")
  
  ci_h <- ci_h %>% mutate(across(all_of(symptoms), ~round(.x,2)))
  colnames(ci_h) <- paste0(colnames(ci_h), "_h")
  
  final_table <- cbind(est, ci_l, ci_h)
  
  for (sym in symptoms) {
    final_table[[paste0(sym, "_CI")]] <- paste0(final_table[[paste0(sym, "_l")]], 
                                                " - ", 
                                                final_table[[paste0(sym, "_h")]])
    
    
    final_table[[paste0(sym, "_rel")]] <- ifelse(final_table[[paste0(sym, "_l")]] > 0 | 
                                                   final_table[[paste0(sym, "_h")]] < 0, "*", "")
  }
  all_columns <- unlist(lapply(symptoms, function(sym) {
    c(paste0(sym, "_est"), paste0(sym, "_CI"), paste0(sym, "_rel"))
  }))
  
  final_table <- final_table %>% select(all_of(all_columns))
  
  final_tables[[paste0("t", t)]] <- final_table
}

write.xlsx(final_tables, file = "Edge_parameters_across_timepoints.xlsx")



############### Calculation of network parameters ##############################
load("edge_matrices.rda")

symptoms<- c("Depressive mood", "Loss of interest", "Decrease in activity", 
             "Nervousness", "Physical tension", "Autonomic arousal")

week <- as.character(0:7)
lines<- 1:3000

#### calculate instrength and outstrength ####
# calculate instrength and outstrength for every symptom and every timepoint
# for the whole posterior distribution of edges

outstr <- as.data.frame(matrix(ncol = 8, nrow = 3000)) # ncol = no of symps + "lines" + "week
colnames(outstr) <- c("lines", "week", symptoms)
outstr$lines <- lines
time_outstr <- replicate(8, outstr, simplify=FALSE) # replicate no. of timepoints
names(time_outstr) <- week
time_instr <- time_outstr

for(i in week){ 
  for(k in lines){
    for(l in 1:6){
      time_outstr[[i]][,"week"]<- i
      
      time_outstr[[i]][k,l+2]<-
        sum(abs(time_post_edge_m[[i]][[k]][l,-l]))
      
      time_instr[[i]][,"week"]<- i
      
      time_instr[[i]][k,l+2]<-
        sum(abs(time_post_edge_m[[i]][[k]][-l,l])) 
    }
  }
}

save(time_outstr, file="outstrength.rda")
save(time_instr, file="instrength.rda")


#### calculate expected influence ####
out_ei <- as.data.frame(matrix(ncol = 8, nrow = 3000)) # ncol = no of symps + "lines" + "week
colnames(out_ei) <- c("lines", "week", symptoms)
out_ei$lines <- lines
time_out_ei <- replicate(8,out_ei, simplify=FALSE) # replicate no. of timepoints
names(time_out_ei) <- week
time_in_ei <- time_out_ei

for(i in week){
  for(k in lines){
    for(l in 1:6){
      time_out_ei[[i]][,"week"]<- i
      
      time_out_ei[[i]][k,l+2]<-
        sum(time_post_edge_m[[i]][[k]][l,-l])
      
      time_in_ei[[i]][,"week"]<- i
      
      time_in_ei[[i]][k,l+2]<-
        sum(time_post_edge_m[[i]][[k]][-l,l])
    }
  }
}

save(time_out_ei, file="out_EI.rda")
save(time_in_ei, file="in_EI.rda")


#### calculate global network parameters ####
# calculate the sum of positive edge weights, the sum of negative edge weights
# and the raw and absolute average edge weight for the complete posterior  
# (without autoregressive effects)

global_net_par <- as.data.frame(matrix(nrow =3000, ncol = 6)) # ncol = see names next line
colnames(global_net_par) <- c("line", "week", "sum_pos_edges", "sum_neg_edges", "abs_average_edge_weight", "raw_average_edge_weight")
global_net_par$line <- lines
time_global_net_par <- replicate(8, global_net_par, simplify=FALSE) # replicate for time points (8)
names(time_global_net_par) <- week
egdes_wo_auto <- time_post_edge_m


for(i in week){
  for(k in lines){
    time_global_net_par[[i]][k, "week"] <- i
    diag(egdes_wo_auto[[i]][[k]]) <- 0
    
    time_global_net_par[[i]][k, "sum_pos_edges"] <- sum(
      egdes_wo_auto[[i]][[k]][egdes_wo_auto[[i]][[k]] > 0])
    
    time_global_net_par[[i]][k, "sum_neg_edges"] <- sum(
      egdes_wo_auto[[i]][[k]][egdes_wo_auto[[i]][[k]] < 0])
    
    time_global_net_par[[i]][k, "abs_average_edge_weight"] <- mean(
     (abs(as.matrix(egdes_wo_auto[[i]][[k]]))))
    
    time_global_net_par[[i]][k, "raw_average_edge_weight"] <- mean(
      as.matrix(egdes_wo_auto[[i]][[k]]))
  }
}

save(time_global_net_par, file = "global_net_par.rda")


## get the mean absolute and raw average edge weight before and after treatment
load("global_net_par.rda")

## absolute average edge weight
# treatment start
mean(time_global_net_par[[1]]$abs_average_edge_weight) # 0.5478
sort(time_global_net_par[[1]]$abs_average_edge_weight)[76]# 0.4651
sort(time_global_net_par[[1]]$abs_average_edge_weight)[2924] # 0.6401

# treatment end
mean(time_global_net_par[[8]]$abs_average_edge_weight) # 0.4601
sort(time_global_net_par[[8]]$abs_average_edge_weight)[76]# 0.3822
sort(time_global_net_par[[8]]$abs_average_edge_weight)[2924] # 0.5443


## raw average edge weight 
# treatment start
mean(time_global_net_par[[1]]$raw_average_edge_weight) # 0.4453
sort(time_global_net_par[[1]]$raw_average_edge_weight)[76] # 0.3788
sort(time_global_net_par[[1]]$raw_average_edge_weight)[2924] #  0.5127

# treatment end
mean(time_global_net_par[[8]]$raw_average_edge_weight) # 0.3664
sort(time_global_net_par[[8]]$raw_average_edge_weight)[76] #  0.2993
sort(time_global_net_par[[8]]$raw_average_edge_weight)[2924] # 0.4336


################################
######## 5. Plotting  ##########
################################

library("ggplot2")
library("qgraph")
library("ggpubr")
library("jtools")

load("edge_matrices.rda")
load("mean_edge_matrix.rda") 

symptoms <- c("Depressive mood", "Loss of interest", "Decrease in activity", 
             "Nervousness", "Physical tension", "Autonomic arousal")

week <- as.character(0:7)
lines <- 1:3000

#### Plot networks for all sessions ####

# find the maximum edge weight without autoregressive effects 
time_mean_edge_m_wa <- time_mean_edge_m # matrix without autoregressive effects
for(i in week){
  diag(time_mean_edge_m_wa[[i]]) <- 0
}

lapply(lapply(time_mean_edge_m_wa,FUN=abs), FUN = max) # 1.37
# set the max to 1.4

# create a list with all network plots 
graph_List<- list()

for (i in 1:8) {
  graph_List[[i]]<- qgraph(time_mean_edge_m[[i]],layout = "circle",
                           vsize=9, edge.width = 1.4, labels = c(1:6), 
                           asize = 6, minimum =0.1, maximum= 1.4,  legend = F,
                           nodeNames = symptoms, title = paste(i), theme = "colorblind")
}

# plot all network plots in a pdf (quick plotting, no further editing)
pdf("all_symptom_networks.pdf", width = 4.5 , height = 5.8)
par(mfrow= c(3,2))

for (i in 1:8){
  plot(graph_List[[i]])
}
dev.off()

# include markers for substantial edges and those that changed substantially during treatment with 95% certainty
p <- 6

subst_t0 <- c(1, 2, 7, 8, 9, 14, 15, 22, 23, 28, 29, 36)
subst_t1 <- c(1, 2, 7, 8, 9, 14, 15, 22, 23, 25, 28, 29, 36)
subst_t2 <- c(1, 2, 7, 8, 9, 13, 14, 15, 22, 23, 24, 25, 28, 29, 30, 36)
subst_t3 <- c(1, 2, 3, 6, 7, 8, 9, 12, 13, 14, 15, 22, 23, 24, 25, 28, 29, 30, 31, 36)
subst_t4 <- c(1, 2, 3, 7, 8, 9, 13, 14, 15, 22, 23, 24, 25, 27, 28, 29, 30, 31, 36)
subst_t5 <- c(1, 2, 3, 7, 8, 9, 13, 14, 15, 22, 23, 24, 25, 27, 29, 30, 31, 36)
subst_t6 <- c(1, 2, 7, 8, 9, 13, 14, 15, 22, 23, 24, 25, 27, 29, 30, 31, 36)
subst_t7 <- c(1, 2, 9, 22, 23, 29, 30)

# Color function: blue/red for substantial, grey for non-substantial
make_edgecolor <- function(graph_matrix, subst_idx, p = 6){
  base_colors <- ifelse(graph_matrix > 0, "blue3",
                        ifelse(graph_matrix < 0, "red", NA))
  col_mat <- matrix("grey60", p, p)
  col_mat[subst_idx] <- base_colors[subst_idx]
  col_mat
}

# lty function: 
#   - substantial edges (positive or negative) -> solid
#   - non-substantial edges -> solid if positive, dashed if negative
make_lty <- function(graph_matrix, subst_idx, p = 6){
  L <- matrix(1, p, p)                 # solid by default
  is_subst <- matrix(FALSE, p, p)
  is_subst[subst_idx] <- TRUE
  L[!is_subst & graph_matrix < 0] <- 2 # dashed ONLY for non-substantial negative edges
  L
}

# Generate color and lty matrices for all 8 timepoints
col_t0 <- make_edgecolor(time_mean_edge_m[[1]], subst_t0, p)
col_t1 <- make_edgecolor(time_mean_edge_m[[2]], subst_t1, p)
col_t2 <- make_edgecolor(time_mean_edge_m[[3]], subst_t2, p)
col_t3 <- make_edgecolor(time_mean_edge_m[[4]], subst_t3, p)
col_t4 <- make_edgecolor(time_mean_edge_m[[5]], subst_t4, p)
col_t5 <- make_edgecolor(time_mean_edge_m[[6]], subst_t5, p)
col_t6 <- make_edgecolor(time_mean_edge_m[[7]], subst_t6, p)
col_t7 <- make_edgecolor(time_mean_edge_m[[8]], subst_t7, p)

lty_t0 <- make_lty(time_mean_edge_m[[1]], subst_t0, p)
lty_t1 <- make_lty(time_mean_edge_m[[2]], subst_t1, p)
lty_t2 <- make_lty(time_mean_edge_m[[3]], subst_t2, p)
lty_t3 <- make_lty(time_mean_edge_m[[4]], subst_t3, p)
lty_t4 <- make_lty(time_mean_edge_m[[5]], subst_t4, p)
lty_t5 <- make_lty(time_mean_edge_m[[6]], subst_t5, p)
lty_t6 <- make_lty(time_mean_edge_m[[7]], subst_t6, p)
lty_t7 <- make_lty(time_mean_edge_m[[8]], subst_t7, p)

# Edge labels — mark edges that changed substantially during treatment
lable_subst <- c(rep(NA, 36))
lable_subst[c(8, 30, 36)] <- "O"

# Layout: 4 rows x 2 columns
layout.n <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4, byrow = TRUE)

pdf("all_symptom_networks_marked.pdf", width = 12, height = 21)
# png("all_symptom_networks_marked.png", width = 12, height = 21, units = "cm", res = 600)
layout(layout.n)
par(oma = c(1, 2, 2, 2), mar = c(1, 1, 1, 1))
cex_panel <- par("pin")[2] / 3

# Week 0
qgraph(time_mean_edge_m[[1]], layout = "circle",
       vsize = 13, edge.width = 1.4,
       labels = c("DM", "IL", "DA", "NER", "TEN", "ARO"), label.cex = 1,
       asize = 6, minimum = 0.1, maximum = 1.4,
       legend = FALSE, nodeNames = symptoms,
       edge.labels = lable_subst, edge.label.cex = 1.8,
       edge.label.color = "black", edge.label.position = 0.55,
       edge.color = col_t0, lty = lty_t0)
mtext("t0", side = 3, line = 0.2, adj = 0.05, cex = cex_panel)

# Week 1
qgraph(time_mean_edge_m[[2]], layout = "circle",
       vsize = 13, edge.width = 1.4,
       labels = c("DM", "IL", "DA", "NER", "TEN", "ARO"), label.cex = 1,
       asize = 6, minimum = 0.1, maximum = 1.4,
       legend = FALSE, nodeNames = symptoms,
       edge.labels = lable_subst, edge.label.cex = 1.8,
       edge.label.color = "black", edge.label.position = 0.55,
       edge.color = col_t1, lty = lty_t1)
mtext("t1", side = 3, line = 0.2, adj = 0.05, cex = cex_panel)

# Week 2
qgraph(time_mean_edge_m[[3]], layout = "circle",
       vsize = 13, edge.width = 1.4,
       labels = c("DM", "IL", "DA", "NER", "TEN", "ARO"), label.cex = 1,
       asize = 6, minimum = 0.1, maximum = 1.4,
       legend = FALSE, nodeNames = symptoms,
       edge.labels = lable_subst, edge.label.cex = 1.8,
       edge.label.color = "black", edge.label.position = 0.55,
       edge.color = col_t2, lty = lty_t2)
mtext("t2", side = 3, line = 0.2, adj = 0.05, cex = cex_panel)

# Week 3
qgraph(time_mean_edge_m[[4]], layout = "circle",
       vsize = 13, edge.width = 1.4,
       labels = c("DM", "IL", "DA", "NER", "TEN", "ARO"), label.cex = 1,
       asize = 6, minimum = 0.1, maximum = 1.4,
       legend = FALSE, nodeNames = symptoms,
       edge.labels = lable_subst, edge.label.cex = 1.8,
       edge.label.color = "black", edge.label.position = 0.55,
       edge.color = col_t3, lty = lty_t3)
mtext("t3", side = 3, line = 0.2, adj = 0.05, cex = cex_panel)

# Week 4
qgraph(time_mean_edge_m[[5]], layout = "circle",
       vsize = 13, edge.width = 1.4,
       labels = c("DM", "IL", "DA", "NER", "TEN", "ARO"), label.cex = 1,
       asize = 6, minimum = 0.1, maximum = 1.4,
       legend = FALSE, nodeNames = symptoms,
       edge.labels = lable_subst, edge.label.cex = 1.8,
       edge.label.color = "black", edge.label.position = 0.55,
       edge.color = col_t4, lty = lty_t4)
mtext("t4", side = 3, line = 0.2, adj = 0.05, cex = cex_panel)

# Week 5
qgraph(time_mean_edge_m[[6]], layout = "circle",
       vsize = 13, edge.width = 1.4,
       labels = c("DM", "IL", "DA", "NER", "TEN", "ARO"), label.cex = 1,
       asize = 6, minimum = 0.1, maximum = 1.4,
       legend = FALSE, nodeNames = symptoms,
       edge.labels = lable_subst, edge.label.cex = 1.8,
       edge.label.color = "black", edge.label.position = 0.55,
       edge.color = col_t5, lty = lty_t5)
mtext("t5", side = 3, line = 0.2, adj = 0.05, cex = cex_panel)

# Week 6
qgraph(time_mean_edge_m[[7]], layout = "circle",
       vsize = 13, edge.width = 1.4,
       labels = c("DM", "IL", "DA", "NER", "TEN", "ARO"), label.cex = 1,
       asize = 6, minimum = 0.1, maximum = 1.4,
       legend = FALSE, nodeNames = symptoms,
       edge.labels = lable_subst, edge.label.cex = 1.8,
       edge.label.color = "black", edge.label.position = 0.55,
       edge.color = col_t6, lty = lty_t6)
mtext("t6", side = 3, line = 0.2, adj = 0.05, cex = cex_panel)

# Week 7
qgraph(time_mean_edge_m[[8]], layout = "circle",
       vsize = 13, edge.width = 1.4,
       labels = c("DM", "IL", "DA", "NER", "TEN", "ARO"), label.cex = 1,
       asize = 6, minimum = 0.1, maximum = 1.4,
       legend = FALSE, nodeNames = symptoms,
       edge.labels = lable_subst, edge.label.cex = 1.8,
       edge.label.color = "black", edge.label.position = 0.55,
       edge.color = col_t7, lty = lty_t7)
mtext("t7", side = 3, line = 0.2, adj = 0.05, cex = cex_panel)

dev.off()



#### Plot for main timepoints (i.e., at baseline, 3, and 7 weeks) ####
layout.m <- matrix(c(0, 1, 1, 2, 2, 0,
                     0, 0, 3, 3, 0, 0), nrow = 2, byrow = TRUE)

pdf("symptom_networks_main_timepoints.pdf", width = 14, height = 10)

#pdf Layout
layout(layout.m, widths = c(1, 2, 2, 2, 2, 1), heights = c(2, 2))

par(oma = c(2, 2, 2, 2), mar = c(3, 3, 4, 2))

qgraph(time_mean_edge_m[[1]], layout = "circle",
       vsize = 13, edge.width = 1.4, theme = "colorblind",
       labels = c("DM", "IL", "DA", "NER", "TEN", "ARO"), label.cex = 1.2,
       asize = 6, minimum = 0.1, maximum = 1.4,
       legend = FALSE, nodeNames = symptoms,
       edge.labels = lable_subst, edge.label.cex = 1.8,
       edge.label.color = "black", edge.label.position = 0.55, 
       edge.color = col_t0, lty = lty_t0)
mtext(text = "T0", outer = FALSE, line = 2, adj = 0.05, cex = 1.5)

# Plot B
qgraph(time_mean_edge_m[[4]], layout = "circle",
       vsize = 13, edge.width = 1.4,
       labels = c("DM", "IL", "DA", "NER", "TEN", "ARO"), label.cex = 1.2,
       asize = 6, minimum = 0.1, maximum = 1.4,
       legend = FALSE, nodeNames = symptoms,
       edge.labels = lable_subst, edge.label.cex = 1.8,
       edge.label.color = "black", edge.label.position = 0.55, 
       edge.color = col_t3, lty = lty_t3)
mtext(text = "T3", outer = FALSE, line = 2, adj = 0.05, cex = 1.5)

# Plot C
qgraph(time_mean_edge_m[[8]], layout = "circle",
       vsize = 13, edge.width = 1.4,
       labels = c("DM", "IL", "DA", "NER", "TEN", "ARO"), label.cex = 1.2,
       asize = 6, minimum = 0.1, maximum = 1.4,
       legend = FALSE, nodeNames = symptoms,
       edge.labels = lable_subst, edge.label.cex = 1.8,
       edge.label.color = "black", edge.label.position = 0.55, 
       edge.color = col_t7, lty = lty_t7)
mtext(text = "T7", outer = FALSE, line = 2, adj = 0.05, cex = 1.5)

dev.off()


#### Plot symptom specific network parameters ####
load("outstrength.rda")
load("instrength.rda")

## strength
# combine the whole posterior of all timepoints in one large matrix
all_outstr <- do.call("rbind",time_outstr)
all_instr <- do.call("rbind",time_instr)

all_outstr$week <- as.numeric(all_outstr$week)
all_instr$week<- as.numeric(all_instr$week)

# calculate mean in- and outstrength of posterior
mean_outstr <- as.data.frame(matrix(ncol = 7, nrow = 8)) # ncol = no of symptoms + 1, nrow = no of time points
colnames(mean_outstr)<- c("week", symptoms)
mean_outstr$week <- 0:7
mean_instr <- mean_outstr 

for(i in week){
  for(k in symptoms){ 
    mean_outstr[which(mean_outstr$week == i),k]<-mean(
    all_outstr[which(all_outstr$week == i), k])
    
    mean_instr[which(mean_instr$week == i),k]<-mean(
      all_instr[which(all_instr$week == i), k])
    }
}


## expected influence 
load("in_EI.rda")
load("out_EI.rda")

# combine the whole posterior of all time points in one large matrix
all_out_ei <- do.call("rbind",time_out_ei)
all_in_ei <- do.call("rbind",time_in_ei)

all_out_ei$week <- as.numeric(all_out_ei$week)
all_in_ei$week <- as.numeric(all_in_ei$week)

# calculate mean in- and out_strength of posterior
mean_out_ei <- as.data.frame(matrix(ncol = 7, nrow = 8)) # ncol = no of symptoms + 1, nrow = no of time points
colnames(mean_out_ei) <- c("week", symptoms)
mean_out_ei$week<- 0:7
mean_in_ei <- mean_out_ei

for(i in week){
  for(k in symptoms){ 
    mean_out_ei[which(mean_out_ei$week == i),k]<-mean(
      all_out_ei[which(all_out_ei$week == i), k])
    
    mean_in_ei[which(mean_in_ei$week == i),k]<-mean(
      all_in_ei[which(all_in_ei$week == i), k])
  }
}

# plot all symptom specific network parameters
net_par_plot1 <- ggarrange(
# Depressive mood
ggplot(all_instr)+
  geom_line(aes(x=week,y=all_instr[,"Depressive mood"],group= lines),color = "orange",alpha =0.4)+

  geom_line(data=mean_instr,aes(x=week, y= mean_instr[,"Depressive mood"], color="tomato2"), linewidth=1.5)+

  labs( x = "Week", y= "",
        title= paste("Depressive Mood", "- instrength "))+
  scale_x_continuous(breaks=c(seq(0,7,1)))+
  scale_y_continuous(limits = c(0, 10), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_outstr)+
  geom_line(aes(x=week,y=all_outstr[,"Depressive mood"],group= lines),color = "orange",alpha =0.4)+

  geom_line(data=mean_outstr,aes(x=week, y= mean_outstr[,"Depressive mood"], color="tomato2"), linewidth=1.5)+

  labs( x = "Week", y= "",
        title= paste("Depressive mood", "- outstrength "))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(0, 10), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_in_ei)+
  geom_line(aes(x=week,y=all_in_ei[,"Depressive mood"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_in_ei,aes(x=week, y= mean_in_ei[,"Depressive mood"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Depressive mood", "- ingoing EI "))+
  scale_x_continuous(breaks = seq(0,7,1))+         
  scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_out_ei)+
  geom_line(aes(x=week,y=all_out_ei[,"Depressive mood"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_out_ei,aes(x=week, y= mean_out_ei[,"Depressive mood"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Depressive mood", "- outgoing EI "))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
  theme_apa(legend.pos = "none"),


# Loss of interest
ggplot(all_instr)+
  geom_line(aes(x=week,y=all_instr[,"Loss of interest"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_instr,aes(x=week, y= mean_instr[,"Loss of interest"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Loss of interest", "- instrength "))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(0, 11), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_outstr)+
  geom_line(aes(x=week,y=all_outstr[,"Loss of interest"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_outstr,aes(x=week, y= mean_outstr[,"Loss of interest"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Loss of interest", "- outstrength "))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(0, 11), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_in_ei)+
  geom_line(aes(x=week,y=all_in_ei[,"Loss of interest"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_in_ei,aes(x=week, y= mean_in_ei[,"Loss of interest"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Loss of interest", "- ingoing EI "))+
  scale_x_continuous(breaks = seq(0,7,1))+         
  scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_out_ei)+
  geom_line(aes(x=week,y=all_out_ei[,"Loss of interest"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_out_ei,aes(x=week, y= mean_out_ei[,"Loss of interest"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Loss of interest", "- outgoing EI "))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
  theme_apa(legend.pos = "none"),
labels = c(1:8),
ncol = 2, nrow = 4)
net_par_plot1

net_par_plot2<-ggarrange(
# Decrease in activity
  ggplot(all_instr)+
    geom_line(aes(x=week,y=all_instr[,"Decrease in activity"],group= lines),color = "orange",alpha =0.4)+
    
    geom_line(data=mean_instr,aes(x=week, y= mean_instr[,"Decrease in activity"], color="tomato2"), linewidth=1.5)+
    
    labs( x = "Week", y= "",
          title= paste("Decrease in activity", "- instrength "))+
    scale_x_continuous(breaks = seq(0,7,1))+
    scale_y_continuous(limits = c(0, 11), n.breaks = 10) +
    theme_apa(legend.pos = "none"),
  
  ggplot(all_outstr)+
    geom_line(aes(x=week,y=all_outstr[,"Decrease in activity"],group= lines),color = "orange",alpha =0.4)+
    
    geom_line(data=mean_outstr,aes(x=week, y= mean_outstr[,"Decrease in activity"], color="tomato2"), linewidth=1.5)+
    
    labs( x = "Week", y= "",
          title= paste("Decrease in activity", "- outstrength "))+
    scale_x_continuous(breaks = seq(0,7,1))+
    scale_y_continuous(limits = c(0, 11), n.breaks = 10) +
    theme_apa(legend.pos = "none"),
  
  ggplot(all_in_ei)+
    geom_line(aes(x=week,y=all_in_ei[,"Decrease in activity"],group= lines),color = "orange",alpha =0.4)+
    
    geom_line(data=mean_in_ei,aes(x=week, y= mean_in_ei[,"Decrease in activity"], color="tomato2"), linewidth=1.5)+
    
    labs( x = "Week", y= "",
          title= paste("Decrease in activity", "- ingoing EI "))+
    scale_x_continuous(breaks = seq(0,7,1))+         
    scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
    theme_apa(legend.pos = "none"),
  
  ggplot(all_out_ei)+
    geom_line(aes(x=week,y=all_out_ei[,"Decrease in activity"],group= lines),color = "orange",alpha =0.4)+
    
    geom_line(data=mean_out_ei,aes(x=week, y= mean_out_ei[,"Decrease in activity"], color="tomato2"), linewidth=1.5)+
    
    labs( x = "Week", y= "",
          title= paste("Decrease in activity", "- outgoing EI "))+
    scale_x_continuous(breaks = seq(0,7,1))+
    scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
    theme_apa(legend.pos = "none"),

# Nervousness
ggplot(all_instr)+
  geom_line(aes(x=week,y=all_instr[,"Nervousness"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_instr,aes(x=week, y= mean_instr[,"Nervousness"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Nervousness", "- instrength "))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(0, 11), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_outstr)+
  geom_line(aes(x=week,y=all_outstr[,"Nervousness"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_outstr,aes(x=week, y= mean_outstr[,"Nervousness"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Nervousness", "- outstrength "))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(0, 11), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_in_ei)+
  geom_line(aes(x=week,y=all_in_ei[,"Nervousness"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_in_ei,aes(x=week, y= mean_in_ei[,"Nervousness"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Nervousness", "- ingoing EI "))+
  scale_x_continuous(breaks = seq(0,7,1))+         
  scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_out_ei)+
  geom_line(aes(x=week,y=all_out_ei[,"Nervousness"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_out_ei,aes(x=week, y= mean_out_ei[,"Nervousness"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Nervousness", "- outgoing EI "))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
  theme_apa(legend.pos = "none"),
labels = c(9:16),
ncol = 2, nrow = 4)


net_par_plot3<- ggarrange(
# Physical tension
  ggplot(all_instr)+
    geom_line(aes(x=week,y=all_instr[,"Physical tension"],group= lines),color = "orange",alpha =0.4)+
    
    geom_line(data=mean_instr,aes(x=week, y= mean_instr[,"Physical tension"], color="tomato2"), linewidth=1.5)+
    
    labs( x = "Week", y= "",
          title= paste("Physical tension", "- instrength "))+
    scale_x_continuous(breaks = seq(0,7,1))+
    scale_y_continuous(limits = c(0, 11), n.breaks = 10) +
    theme_apa(legend.pos = "none"),
  
  ggplot(all_outstr)+
    geom_line(aes(x=week,y=all_outstr[,"Physical tension"],group= lines),color = "orange",alpha =0.4)+
    
    geom_line(data=mean_outstr,aes(x=week, y= mean_outstr[,"Physical tension"], color="tomato2"), linewidth=1.5)+
    
    labs( x = "Week", y= "",
          title= paste("Physical tension", "- outstrength "))+
    scale_x_continuous(breaks = seq(0,7,1))+
    scale_y_continuous(limits = c(0, 11), n.breaks = 10) +
    theme_apa(legend.pos = "none"),
  
  ggplot(all_in_ei)+
    geom_line(aes(x=week,y=all_in_ei[,"Physical tension"],group= lines),color = "orange",alpha =0.4)+
    
    geom_line(data=mean_in_ei,aes(x=week, y= mean_in_ei[,"Physical tension"], color="tomato2"), linewidth=1.5)+
    
    labs( x = "Week", y= "",
          title= paste("Physical tension", "- ingoing EI "))+
    scale_x_continuous(breaks = seq(0,7,1))+         
    scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
    theme_apa(legend.pos = "none"),
  
  ggplot(all_out_ei)+
    geom_line(aes(x=week,y=all_out_ei[,"Physical tension"],group= lines),color = "orange",alpha =0.4)+
    
    geom_line(data=mean_out_ei,aes(x=week, y= mean_out_ei[,"Physical tension"], color="tomato2"), linewidth=1.5)+
    
    labs( x = "Week", y= "",
          title= paste("Physical tension", "- outgoing EI "))+
    scale_x_continuous(breaks = seq(0,7,1))+
    scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
    theme_apa(legend.pos = "none"),
  
  
# Autonomic arousal
ggplot(all_instr)+
  geom_line(aes(x=week,y=all_instr[,"Autonomic arousal"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_instr,aes(x=week, y= mean_instr[,"Autonomic arousal"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Autonomic arousal", "- instrength "))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(0, 11), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_outstr)+
  geom_line(aes(x=week,y=all_outstr[,"Autonomic arousal"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_outstr,aes(x=week, y= mean_outstr[,"Autonomic arousal"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Autonomic arousal", "- outstrength "))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(0, 11), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_in_ei)+
  geom_line(aes(x=week,y=all_in_ei[,"Autonomic arousal"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_in_ei,aes(x=week, y= mean_in_ei[,"Autonomic arousal"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Autonomic arousal", "- ingoing EI "))+
  scale_x_continuous(breaks = seq(0,7,1))+         
  scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
  theme_apa(legend.pos = "none"),

ggplot(all_out_ei)+
  geom_line(aes(x=week,y=all_out_ei[,"Autonomic arousal"],group= lines),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_out_ei,aes(x=week, y= mean_out_ei[,"Autonomic arousal"], color="tomato2"), linewidth=1.5)+
  
  labs( x = "Week", y= "",
        title= paste("Autonomic arousal", "- outgoing EI "))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(-8, 10), n.breaks = 10) +
  theme_apa(legend.pos = "none"),
labels = c(17:24),
ncol = 2, nrow = 4)


# save the plots in a pdf
pdf("all_node_parameters.pdf", width = 8 , height = 11)
net_par_plot1
net_par_plot2
net_par_plot3

dev.off()


#### Plot global network parameters ####
## NOTE: Scaling tried to be held constant between different variations: 
## sum pos edges: lamda = 17
## mean strength: lamda = 0.8

load("global_net_par.rda")

# combine global network parameters in one dataframe
all_global_net_par <- do.call("rbind",time_global_net_par)
all_global_net_par$week <- as.numeric(all_global_net_par$week)

# calculate the posterior mean for each global network parameter
mean_global_net_par <- as.data.frame(matrix(nrow = 8, ncol = 5)) # nrow = no of time points; ncol = see next line
colnames(mean_global_net_par)<- c("week", "sum_pos_edges", "sum_neg_edges", "abs_average_edge_weight", "raw_average_edge_weight")
mean_global_net_par$week<- 0:7

for(i in 0:7){
  mean_global_net_par[which(mean_global_net_par$week ==i), "sum_pos_edges"]<- mean(
    all_global_net_par[which(all_global_net_par$week == i),"sum_pos_edges"])
  
  mean_global_net_par[which(mean_global_net_par$week ==i), "sum_neg_edges"]<- mean(
    all_global_net_par[which(all_global_net_par$week == i),"sum_neg_edges"])
  
  mean_global_net_par[which(mean_global_net_par$week ==i), "abs_average_edge_weight"]<- mean(
    all_global_net_par[which(all_global_net_par$week == i),"abs_average_edge_weight"])  
  
  mean_global_net_par[which(mean_global_net_par$week ==i), "raw_average_edge_weight"]<- mean(
    all_global_net_par[which(all_global_net_par$week == i),"raw_average_edge_weight"])
}

pdf("global_network_parameters_final.pdf", width = 15, height = 10 )
ggarrange(
  
# plot sum of positive edges 
ggplot(all_global_net_par)+
  geom_line(aes(x=week,y=sum_pos_edges,group= line),color = "orange",alpha =0.4)+
  
  geom_line(data=mean_global_net_par,aes(x=week, y= sum_pos_edges, color="tomato2"), linewidth= 1.5)+

  labs( x = "Week", y= "",
        title= "Sum of cross-lagged positive edges")+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(5.5, 22.5))+
  theme_apa(legend.pos = "none")+
  theme(text = element_text(size = 14)),

# plot sum of negative edges 
ggplot(all_global_net_par)+
  geom_line(aes(x=week,y=sum_neg_edges,group= line),color = "orange",alpha =0.4)+

  geom_line(data=mean_global_net_par,aes(x=week, y= sum_neg_edges, color="tomato2"), linewidth= 1.5)+

  labs( x = "Week", y= "",
        title= "Sum of cross-lagged negative edges")+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(-10,0))+
  theme_apa(legend.pos = "none"),

# plot absolute average_edge_weight
ggplot(all_global_net_par)+
  geom_line(aes(x=week,y=abs_average_edge_weight,group= line),color = "orange",alpha =0.4)+

  geom_line(data=mean_global_net_par,aes(x=week, y= abs_average_edge_weight, color="tomato2"), linewidth= 1.5)+

  labs( x = "Week", y= "", title= "Absolute overall strength of symptom interactions (valence-ignorant)")+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(0.15,0.95))+
  theme_apa(legend.pos = "none"),

# plot raw average_edge_weight
ggplot(all_global_net_par)+
  geom_line(aes(x=week,y=raw_average_edge_weight,group= line),color = "orange",alpha =0.4)+

  geom_line(data=mean_global_net_par,aes(x=week, y= raw_average_edge_weight, color="tomato2"), linewidth= 1.5)+

  labs( x = "Week", y= "",
        title= "Raw overall strength of symptom interactions (valence-sensitive) ")+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(limits = c(0.15,0.95))+
  theme_apa(legend.pos = "none"),
labels = c(LETTERS[1:4]),
ncol = 2, nrow = 2)
dev.off()

# only mean strength 
pdf("mean_strength.pdf", width = 6, height = 8)
ggarrange(
  # plot absolute average_edge_weight
  ggplot(all_global_net_par)+
    geom_line(aes(x=week,y=abs_average_edge_weight,group= line),color = "orange",alpha =0.4)+
    
    geom_line(data=mean_global_net_par,aes(x=week, y= abs_average_edge_weight, color="tomato2"), linewidth= 1.5)+
    
    labs( x = "Week", y= "logit scale", title= "Average absolute strength (valence-ignorant)")+
    scale_x_continuous(breaks = seq(0,7,1))+
    scale_y_continuous(limits = c(0.2,0.8), breaks = seq(0.2, 0.8, by = 0.1))+
    theme_apa(legend.pos = "none"),
  
  # plot raw average_edge_weight
  ggplot(all_global_net_par)+
    geom_line(aes(x=week,y=raw_average_edge_weight,group= line),color = "orange",alpha =0.4)+
    
    geom_line(data=mean_global_net_par,aes(x=week, y= raw_average_edge_weight, color="tomato2"), linewidth= 1.5)+
    
    labs( x = "Week", y= "logit scale",
          title= "Average raw strength (valence-sensitive) ")+
    scale_x_continuous(breaks = seq(0,7,1))+
    scale_y_continuous(limits = c(0.2,0.8), breaks = seq(0.2, 0.8, by = 0.1))+
    theme_apa(legend.pos = "none"),
  labels = c(LETTERS[1:2]),
  ncol = 1, nrow = 2)
dev.off()  


############################################################################
##### 6. Sensitivity Analysis: Adding Time as Main Effect to the Model #####
############################################################################
library("brms")
library("ggplot2")
library("ggmcmc")
library("dplyr")

# load data

### Defining priors ###

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
  stat_function(fun = dcauchy, args = list( scale = 3)) +
  ylim(0,0.15) 

# setting priors 
prior1<- c(
  prior(normal(0,3), class= Intercept), # prior on population intercept 
  prior(normal(0,3), class= b), # prior on population slopes 
  prior(cauchy (0,3), class = sd)) # prior on sd of group-level effects

# model1s = model 1, but sensitivity check (for time trend as fixed effect)

## mood 
model1s_mood<-brm(
  formula = mood ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + week + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1s_mood, type = "trace")
mcmc_plot(model1s_mood, type = "rhat_hist")
save(model1s_mood, file = "model1s_mood.rda")
summary(model1s_mood)


## loss_interest 
model1s_loss_interest<-brm(
  formula = loss_interest ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + week + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1s_loss_interest, type = "trace")
mcmc_plot(model1s_loss_interest, type = "rhat_hist")
save(model1s_loss_interest, file = "model1s_loss_interest.rda")
summary(model1s_loss_interest)


## decrease_activity 
model1s_decrease_activity<-brm(
  formula = decrease_activity ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + week +
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1s_decrease_activity, type = "trace")
mcmc_plot(model1s_decrease_activity, type = "rhat_hist")
save(model1s_decrease_activity, file = "model1s_decrease_activity.rda")
summary(model1s_decrease_activity)


## nervous 
model1s_nervous<-brm(
  formula = nervous ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + week + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1s_nervous, type = "trace")
mcmc_plot(model1s_nervous, type = "rhat_hist")
save(model1s_nervous, file = "model1s_nervous.rda")
summary(model1s_nervous)


## tension 
model1s_tension<-brm(
  formula = tension ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + week + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1s_tension, type = "trace")
mcmc_plot(model1s_tension, type = "rhat_hist")
save(model1s_tension, file = "model1s_tension.rda")
summary(model1s_tension)


## arousal 
model1s_arousal<-brm(
  formula = arousal ~ 1 + mood_lag + loss_interest_lag + decrease_activity_lag + 
    nervous_lag + tension_lag + arousal_lag + week + 
    mood_lag:week + loss_interest_lag:week + decrease_activity_lag:week + 
    nervous_lag:week + tension_lag:week + arousal_lag:week +
    (1 |ID), 
  data = dat_N,
  family = bernoulli(link="logit"), 
  prior = prior1,
  chains = 4, warmup = 3000,  iter = 9000, thin=8, cores = 8,
  control = list(adapt_delta = 0.99), 
  save_pars = save_pars(all = TRUE))

mcmc_plot(model1s_arousal, type = "trace")
mcmc_plot(model1s_arousal, type = "rhat_hist")
save(model1s_arousal, file = "model1s_arousal.rda")
summary(model1s_arousal)

save(model1s_mood, model1s_loss_interest, model1s_decrease_activity, 
     model1s_nervous, model1s_tension, model1s_arousal, file = "models1s_3000.rda")


### extract model parameters from sensitivity analysis (time included as fixed effect)
#### create large table with model parameters of all models

library("brms")
library("ggplot2")
library("openxlsx")

load("models1s_3000.rda")

## extract fixed effects and change to odd-ratios
variables(model1s_mood) # first 14 variables include fixed effects (1 Intercept, 7 slopes (week included), 6 change slopes)

extract_OR_table <- function(model, suffix, n_fixed = 13) {
  fit <- as.data.frame(posterior_summary(model)[1:n_fixed, ])
  fit$OR   <- exp(fit$Estimate)
  fit$CI_l <- exp(fit$Q2.5)
  fit$CI_h <- exp(fit$Q97.5)
  fit <- round(fit, 2)
  fit$CI    <- paste(fit$CI_l, fit$CI_h, sep = "-")
  fit$dif_1 <- ifelse(fit$CI_l > 1 | fit$CI_h < 1, "*", "")
  out <- fit[, c("OR", "CI", "dif_1")]
  colnames(out) <- paste(c("OR", "CI", "dif_1"), suffix, sep = "_")
  out
}

models_list_s <- list(
  m = model1s_mood, i = model1s_loss_interest, d = model1s_decrease_activity,
  n = model1s_nervous, t = model1s_tension, a = model1s_arousal
)

# 14 fixed effects: 1 Intercept, 6 lag slopes, 1 week main effect, 6 lag:week interactions
final_table_s <- do.call(cbind, Map(extract_OR_table, models_list_s, names(models_list_s), n_fixed = 14))
write.csv2(final_table_s, "final_table_sensitivity.csv")


#### Model Fit for Sensitivity Analysis for Model 1a ####
library("dplyr")
library("brms")
library("dvmisc")
library("ggplot2")
library("ggpubr")
library("jtools")

# load data and models
load("models1s_3000.rda")
load("preprocessing_workspace.RData")

#### calculate Bayesian pseudo R2 
bayes_R2(model1s_mood)
bayes_R2(model1s_loss_interest)
bayes_R2(model1s_decrease_activity)
bayes_R2(model1s_nervous)
bayes_R2(model1s_tension)
bayes_R2(model1s_arousal)

summary(model1s_mood)
#### compare the estimated and observed symptom prevalence 
symptoms<- c("Depressive mood", "Loss of interest", "Decrease in activity", 
             "Nervousness", "Physical tension", "Autonomic arousal")

# parameter list with fixed effects
par_list <-list(
  as.data.frame(posterior_summary(model1s_mood)[1:14,]), # 14 parameters estimated (1 Intercept, 7 slopes (week included), 6 change slopes)
  as.data.frame(posterior_summary(model1s_loss_interest)[1:14,]),
  as.data.frame(posterior_summary(model1s_decrease_activity)[1:14,]),
  as.data.frame(posterior_summary(model1s_nervous)[1:14,]),
  as.data.frame(posterior_summary(model1s_tension)[1:14,]),
  as.data.frame(posterior_summary(model1s_arousal)[1:14,]))

names(par_list)<- symptoms

# create function for transformation from log-ratios to probabilities
# mathematically identical to plogis()
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# make a dataframe for all estimated symptom probabilities
symp <- as.data.frame(matrix(nrow=8, ncol = 7)) # nrow = no. of measurement points, ncol = no of symptoms + 1
colnames(symp) <- c("week", symptoms)
symp$week <- 0:7

# prevalence for the 0. session is estimated by taking the logodds at session 0
# (logodds = log(P(symptom present)/P(symptom not present))) -> log(P(symptom)/1-P(symptom))
# P(symptom) = mean(symptom) -> binary coded
dat_N <- as.data.frame(dat_N)

symp[1,]<- c(0,c(
  (log(mean(dat_N[dat_N$week == 0, "mood"], na.rm = TRUE)/
         (1-mean(dat_N[dat_N$week == 0, "mood"], na.rm = TRUE)))),
  
  (log(mean(dat_N[dat_N$week == 0, "loss_interest"], na.rm = TRUE)/
         (1-mean(dat_N[dat_N$week == 0, "loss_interest"], na.rm = TRUE)))),
  
  (log(mean(dat_N[dat_N$week == 0, "decrease_activity"], na.rm = TRUE)/
         (1-mean(dat_N[dat_N$week == 0, "decrease_activity"], na.rm = TRUE)))),
  
  (log(mean(dat_N[dat_N$week == 0, "nervous"], na.rm = TRUE)/
         (1-mean(dat_N[dat_N$week == 0, "nervous"], na.rm = TRUE)))),
  
  (log(mean(dat_N[dat_N$week == 0, "tension"], na.rm = TRUE)/
         (1-mean(dat_N[dat_N$week == 0, "tension"], na.rm = TRUE)))),
  
  (log(mean(dat_N[dat_N$week == 0, "arousal"], na.rm = TRUE)/
         (1-mean(dat_N[dat_N$week == 0, "arousal"], na.rm = TRUE))))
))

symp[1,-1]<- logit2prob(symp[1,-1])

# calculate the probability of each symptom during each week of therapy
for(i in 1:7){
  for(k in symptoms){
    symp[which(symp$week == i), k] <- 
      logit2prob(
        par_list[[k]]["b_Intercept", "Estimate"]  +
          par_list[[k]]["b_mood_lag", "Estimate"]*
          symp[which(symp$week == i-1),"Depressive mood"]+
          par_list[[k]]["b_loss_interest_lag", "Estimate"]*
          symp[which(symp$week == i-1),"Loss of interest"]+
          par_list[[k]]["b_decrease_activity_lag", "Estimate"]*
          symp[which(symp$week == i-1),"Decrease in activity"] + 
          par_list[[k]]["b_nervous_lag", "Estimate"]*
          symp[which(symp$week == i-1),"Nervousness"] +
          par_list[[k]]["b_tension_lag", "Estimate"]*
          symp[which(symp$week == i-1),"Physical tension"] +
          par_list[[k]]["b_arousal_lag", "Estimate"]*
          symp[which(symp$week == i-1),"Autonomic arousal"]+
          par_list[[k]]["b_week", "Estimate"]*i+ # added for sensitivity check

          par_list[[k]]["b_mood_lag:week", "Estimate"]* i *
          symp[which(symp$week == i-1),"Depressive mood"]+
          par_list[[k]]["b_loss_interest_lag:week", "Estimate"]* i *
          symp[which(symp$week == i-1),"Loss of interest"]+
          par_list[[k]]["b_decrease_activity_lag:week", "Estimate"]* i *
          symp[which(symp$week == i-1),"Decrease in activity"] + 
          par_list[[k]]["b_nervous_lag:week", "Estimate"]* i *
          symp[which(symp$week == i-1),"Nervousness"] +
          par_list[[k]]["b_tension_lag:week", "Estimate"]* i *
          symp[which(symp$week == i-1),"Physical tension"] +
          par_list[[k]]["b_arousal_lag:week", "Estimate"]* i *
          symp[which(symp$week == i-1),"Autonomic arousal"]
      )
  }
}

# check
symp

## calculate the observed relative frequencies 
obs_symp <- as.data.frame(matrix(nrow=8, ncol = 7)) # nrow = no of time points; ncol= no of symps + 1
colnames(obs_symp) <- c("week", symptoms)
obs_symp$week <- 0:7

for(i in 0:7){
  obs_symp[i+1,2:7]<-c(mean(dat_N[dat_N$week == i, "mood"], na.rm = T),
                       mean(dat_N[dat_N$week == i, "loss_interest"], na.rm = T),
                       mean(dat_N[dat_N$week == i, "decrease_activity"], na.rm = T),
                       mean(dat_N[dat_N$week == i, "nervous"], na.rm = T),
                       mean(dat_N[dat_N$week == i, "tension"], na.rm = T),
                       mean(dat_N[dat_N$week == i, "arousal"], na.rm = T))
}


colnames(symp)[2:7]<- gsub(" ", "_",x= paste(symptoms, "est"))
colnames(obs_symp)[2:7]<-  gsub(" ", "_",x= paste(symptoms, "obs"))

all_symp <- cbind(symp, obs_symp[,-1])

save(all_symp, file = "observedpredictedfrequencies_sensitivity.rda")


# plot the observed and predicted frequencies over time
colnames(all_symp)

obs_pred_Symp <- ggarrange(
  # mood
  ggplot(all_symp) + 
    geom_point(aes(x = week , y= Depressive_mood_obs))+
    geom_smooth(aes(x = week , y= Depressive_mood_obs), se=F)+
    
    geom_point(aes(x = week , y= Depressive_mood_est), shape = 1)+
    geom_line(aes(x = week , y= Depressive_mood_est), linetype ="dashed")+
    
    labs(title = "Depressive mood", x = "week", y= "Observed frequency/Predicted probability")+
    
    scale_x_continuous(breaks=c(seq(0,7,1)))+
    scale_y_continuous(limits = c(0,1))+
    theme_apa(),
  
  # loss interest
  ggplot(all_symp)+ 
    geom_point(aes(x = week , y= Loss_of_interest_obs))+
    geom_smooth(aes(x = week , y= Loss_of_interest_obs), se=F)+
    
    geom_point(aes(x = week , y= Loss_of_interest_est), shape = 1)+
    geom_line(aes(x = week , y= Loss_of_interest_est), linetype ="dashed")+
    
    labs(title = "Loss of interest", x = "week", y= "Observed frequency/Predicted probability")+
    
    scale_x_continuous(breaks=c(seq(0,7,1)))+
    scale_y_continuous(limits = c(0,1))+
    theme_apa(),
  
  # decrease interest
  ggplot(all_symp)+ 
    geom_point(aes(x = week , y= Decrease_in_activity_obs))+
    geom_smooth(aes(x = week , y= Decrease_in_activity_obs), se=F)+
    
    geom_point(aes(x = week , y= Decrease_in_activity_est), shape = 1)+
    geom_line(aes(x = week , y= Decrease_in_activity_est),linetype ="dashed")+
    
    labs(title = "Decrease in activity", x = "week", y= "Observed frequency/Predicted probability")+
    
    scale_x_continuous(breaks=c(seq(0,7,1)))+
    scale_y_continuous(limits = c(0,1))+
    theme_apa(),
  
  # Nervousness
  ggplot(all_symp)+ 
    geom_point(aes(x = week , y= Nervousness_obs))+
    geom_smooth(aes(x = week , y=Nervousness_obs), se=F)+
    
    geom_point(aes(x = week , y=Nervousness_est), shape = 1)+
    geom_line(aes(x = week , y=Nervousness_est), linetype ="dashed")+
    
    labs(title = "Nervousness", x = "week", y= "Observed frequency/Predicted probability")+
    
    scale_x_continuous(breaks=c(seq(0,7,1)))+
    scale_y_continuous(limits = c(0,1))+
    theme_apa(),
  
  # Physical tension
  ggplot(all_symp)+ 
    geom_point(aes(x = week , y= Physical_tension_obs))+
    geom_smooth(aes(x = week , y= Physical_tension_obs), se=F)+
    
    geom_point(aes(x = week , y= Physical_tension_est), shape = 1)+
    geom_line(aes(x = week , y= Physical_tension_est), linetype ="dashed")+
    
    labs(title = "Physical tension", x = "week", y= "Observed frequency/Predicted probability")+
    
    scale_x_continuous(breaks=c(seq(0,7,1)))+
    scale_y_continuous(limits = c(0,1))+
    theme_apa(),
  
  # Autonomic arousal
  ggplot(all_symp)+ 
    geom_point(aes(x = week , y= Autonomic_arousal_obs))+
    geom_smooth(aes(x = week , y= Autonomic_arousal_obs), se=F)+
    
    geom_point(aes(x = week , y=Autonomic_arousal_est), shape = 1)+
    geom_line(aes(x = week , y=Autonomic_arousal_est), linetype ="dashed")+
    
    labs(title = "Autonomic arousal", x = "week", y= "Observed frequency/Predicted probability")+
    
    scale_x_continuous(breaks=c(seq(0,7,1)))+
    scale_y_continuous(limits = c(0,1))+
    theme_apa(), 
  
  nrow= 2, ncol = 3)

pdf("observed_and_predicted_frequencies_sensitivity.pdf", width = 12 , height = 10) 
obs_pred_Symp
dev.off()

