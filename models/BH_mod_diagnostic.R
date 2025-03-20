##script to run models for initial diagnostics, different conditions from "full" model
library(tictoc)
library(beepr)
library(tidyverse)
library(rstan)
setwd("~/Desktop/career_r/data")
model.dat <- read.csv("career_model_data.csv")%>%
  filter(!(is.na(num_seeds)))##### removes samples that are in datasheet but weren't found during processing
setwd("~/Desktop/career_repo")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

date <- Sys.Date()
species <- c("ACAM","ERBO","GITR","LOMU","PLER","TACA")
#treatment<-c("D","A","AC","AG","ANG","ING","IN","DNG","IG","AN","DN","ACN","DG","I")
treatment<-c("D","A","I")### only using water trts for quick models

model.output <- list()
warnings <- list()

for(i in species){
  for(j in treatment){
    tic("model time")
    dat <- subset(model.dat, phyto == i & treatment==j)
    
    ## create vectors of the various data inputs
    Fecundity <- as.integer(round(dat$num_seeds)) ## seeds out, is it the number of seeds or the average?
    N <- as.integer(length(Fecundity)) ## number of observations
    if (N!=0){print(paste(i,j,N,sep = "_"))}
    else next ####### I think this logic loop runs through itself entirely before moving on
    
    N_blocks <- as.integer(length(unique(dat$block))) ## number of blocks
    block <- as.integer(unique(dat$block))## vector of block vals, might need to be changed if subsetting by treatment
    # for each unique value in Blocks, reorder starting at 1 for the lowest value in vector, then 2 for second lowest etc.
    if (length(block)==4){
      block_new<-c(1,2,3,4)
    } else if (length(block)==3){
      block_new<-c(1,2,3)
    } else if (length(block)==2){
      block_new<-c(1,2)
    }else if (length(block)==1){
      block_new<-1
    }
    bl_rep<-as.data.frame(cbind(block,block_new))
    dat<-left_join(dat,bl_rep,by="block")
    Blocks<-as.integer(dat$block_new)
    
    N_i<- as.integer(dat$ph_n_indiv) ## stem # of focal species
    g_i<-dat$mean_germ
    s_i<-dat$mean_surv
    
    ## stems data
    taca<-as.integer(dat$TACA)
    acam<-as.integer(dat$ACAM)
    avba<-as.integer(dat$AVBA)
    erbo<-as.integer(dat$ERBO)
    gitr<-as.integer(dat$GITR)
    lomu<-as.integer(dat$LOMU)
    pler<-as.integer(dat$PLER)
    
    if (i== "TACA"){
      weeds <- as.integer(dat$total_weeds_n) 
    } else if(i=="LOMU"){
      weeds <- as.integer(dat$total_weeds_n)
    } else if(i == "GITR"){
      weeds<-as.integer(dat$total_weed_cover)
    } else if(i== "PLER"){
      weeds<-as.integer(dat$total_weed_cover)
    }else if (i=="ERBO"){
      weeds<-as.integer(dat$total_weed_cover)
    }
    else if (i=="ACAM"){
      weeds<-as.integer(dat$neigh_weeds)
    }
    
    
    ## make a vector of data inputs to model
    
    data_vec <- c("N","Fecundity","N_i","g_i","s_i","N_blocks","Blocks","weeds","acam","avba","erbo","gitr","lomu","pler","taca")
    
    
    ## create initials for epsilon and sigma
    initials <- list(epsilon=rep(1,N_blocks),sigma = 1,alpha_acam=1,alpha_avba=1,alpha_erbo=1,
                     alpha_gitr=1,alpha_lomu=1,alpha_pler=1,alpha_taca=1,alpha_weeds=1)
    
    initials1<- list(initials, initials, initials,initials)
    
    # Model ####
    print(paste("running",i,j,sep = "_"))
    PrelimFit <- stan(file="~/Desktop/career_repo/models/beverton_holt_negbi.stan", model_name="beverton_holt_negbi",
                      data = data_vec, init = initials1, iter = 2000, chains = 4, cores=4, thin=1,
                      control = list(adapt_delta = 0.9, max_treedepth = 12)) 
    
    
    ## save model output
    save(PrelimFit, file = paste0("~/Desktop/test_mods/",i,"_",j,"_posteriors_BH_erbo_update", date ,".rdata"))
    beep(8)
    toc()
  }
}
