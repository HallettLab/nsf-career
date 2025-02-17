## Model 
setwd("~/Desktop/career_r")
model.dat <- read.csv("taca_mod_dat.csv")
#model.dat.filt <- read.csv("data/model_dat.csv")

library(tictoc)
library(beepr)
library(tidyverse)
library(bayesplot)
library(rstan)
library(here)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

date <- Sys.Date()
species <- c("TACA")

model.output <- list()
warnings <- list()


for(i in species){
  tic("model time")
  dat <- subset(model.dat, phyto == i)
  
  ## create vectors of the various data inputs
  Fecundity <- as.integer(round(dat$num_seeds)) ## seeds out, is it the number of seeds or the average?
  N_blocks <- as.integer(length(unique(dat$block))) ## number of blocks
  Blocks <- as.integer(dat$block) ## vector of block vals
  
  N <- as.integer(length(Fecundity)) ## number of observations
  N_i <- as.integer(dat$ph_n_indiv) ## stem # of focal species
  trt <- as.integer(dat$trt) ## treatment (many)
  g_i<-dat$mean_germ
  s_i<-dat$mean_surv
  ## stems data
  #taca<-as.integer(dat$TACA)
  acam<-as.integer(dat$ACAM)
  avba<-as.integer(dat$AVBA)
  gitr<-as.integer(dat$GITR)
  lomu<-as.integer(dat$LOMU)
  pler<-as.integer(dat$PLER)
  weeds <- as.integer(dat$total_weeds)
  
  ## make a vector of data inputs to model
  
  data_vec <- c("N", "Fecundity", "N_i","g_i","s_i", "N_blocks", "Blocks", "trt", "weeds","acam","avba","gitr","lomu","pler")
  
  print(i)
  
  ## create initials for epsilon and sigma
  initials <- list(epsilon=rep(1,N_blocks), sigma = 1,lambda_base=2,lambda_dev=1,alpha_acam_base=2,alpha_avba_base=2,alpha_gitr_base=2,alpha_lomu_base=2,
                   alpha_pler_base=2,alpha_weeds_base=2,alpha_acam_dev=1,alpha_avba_dev=1,alpha_gitr_dev=1,alpha_lomu_dev=1,alpha_pler_dev=1,
                   alpha_weed_dev=1)
  initials1<- list(initials, initials, initials)
  
  # Model ####
  
PrelimFit <- stan(file = 'beverton_holt_negbi.stan', 
                   data = data_vec, init = initials1, iter = 20000, chains = 3, cores=3,, 
                   control = list(adapt_delta = 0.9, max_treedepth = 15)) 
  
  ## save model output
  save(PrelimFit, file = paste0(i, "_posteriors_BH_negbin", date ,".rdata"))
  beep(8)
  toc()
}
######## diag
launch_shinystan(PrelimFit)