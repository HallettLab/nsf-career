## Model script

library(tidyverse)
library(bayesplot)
library(rstan)
#library(here)
setwd("~/Data")
model.dat <- read.csv("career_model_data_Aug20_25.csv")%>%
  filter(!(is.na(num_seeds)))%>%##### removes samples that are in datasheet but weren't found during processing
  mutate_at((vars("weeds","sp_mean_lambda_mc")),~replace_na(.,0))%>%
  mutate_at("sp_sd_lambda_mc",~replace_na(.,1)) #### some Mono samples (and TRWI bkgs) have weeds as NA

setwd("~/R_scripts")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

date <- Sys.Date()
species <- c("ACAM","AVBA","ERBO","GITR","LOMU","PLER","TACA","TRWI")
treatment<-c("D","A","AC","AG","ANG","ING","IN","DNG","IG","AN","DN","ACN","DG","I")
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
    # for each unique value in blocks, reorder starting at 1 for the lowest value in vector, then 2 for second lowest etc.
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
    mc_mean_lambda <- unique(dat$sp_mean_lambda_mc)
    mc_sd_lambda <- unique(dat$sp_sd_lambda_mc)
    
    ## stems data
    taca<-as.integer(dat$TACA)
    acam<-as.integer(dat$ACAM)
    avba<-as.integer(dat$AVBA)
    erbo<-as.integer(dat$ERBO)
    gitr<-as.integer(dat$GITR)
    lomu<-as.integer(dat$LOMU)
    pler<-as.integer(dat$PLER)
    avba<-as.integer(dat$AVBA)
    trwi<-as.integer(dat$TRIW)
    
    
    if (i== "TACA"){
      weeds <- as.integer(dat$total_weeds_n) 
    } else if(i=="LOMU"){
      weeds <- as.integer(dat$total_weeds_n)
    } else if(i=="AVBA"){
      weeds <- as.integer(dat$total_weeds_n) 
    } else if(i == "GITR"){
      weeds<-as.integer(dat$total_weed_cover)
    } else if(i== "PLER"){
      weeds<-as.integer(dat$total_weed_cover)
    }else if (i=="ERBO"){
      weeds<-as.integer(dat$total_weed_cover)
    }else if (i=="TRWI"){
      weeds<-as.integer(dat$total_weed_cover)
    }
    else if (i=="ACAM"){
      weeds<-as.integer(dat$neigh_weeds)
    }
    
    
    ## make a vector of data inputs to model
    
    data_vec <- c("N","Fecundity","N_i","g_i","s_i","mc_mean_lambda","mc_sd_lambda","N_blocks","Blocks","weeds","acam","avba","erbo","gitr","lomu","pler","taca","trwi")
    
    
    ## create initials for epsilon and sigma
    initials <- list(epsilon=rep(1,N_blocks), sigma = 1,alpha_acam=0.1,alpha_avba=0.1,alpha_erbo=0.1,
                     alpha_gitr=0.1,alpha_lomu=0.1,alpha_pler=0.1,alpha_taca=0.1,alpha_trwi=0.1,alpha_weeds=0.1)
    
    initials1<- list(initials, initials, initials, initials)
    
    # Model ####
    print(paste("running",i,j,sep = "_"))
    PrelimFit <- stan(file="~/R_scripts/beverton_holt_negbi.stan", model_name="beverton_holt_negbi",
                      data = data_vec, init = initials1, iter = 20000, chains = 4, cores=4, thin=1,
                      control = list(adapt_delta = 0.9, max_treedepth = 12)) 
    
    
    ## save model output
    save(PrelimFit, file = paste0("~/model_outputs/",i,"_",j,"_posteriors_BH_negbin", date ,".rdata"))

  }
}


###################
########## now extract posteriors
species2<-c("ACAM","AVAB","ERBO","GITR","LOMU","PLER","TACA","TRWI")
treatment2<-c("D","A","AC","AG","ANG","ING","IN","DNG","IG","AN","DN","ACN","DG","I")

#c("D","A","AC","IN","AN","DN","I")
#"D","A","AC","AG","ANG","ING","IN","DNG","IG","AN",
for (i in species2){
  for (j in treatment2){
    if(file.exists(paste0("~/Desktop/talapas_model_outputs/Feb_09_26/",i,"_",j,"_posteriors_BH_negbin2026-02-09.rdata"))==TRUE){
      load(file=paste0("~/Desktop/talapas_model_outputs/Feb_09_26/",i,"_",j,"_posteriors_BH_negbin2026-02-09.rdata"))
    } else next
    
    PrelimPosteriors <- extract(PrelimFit) ###
    
    ##### Diagnostic plots
    # First check the distribution of Rhats and effective sample sizes
    #hist(summary(PrelimFit)$summary[,"Rhat"])
    rhat_plot<-stan_rhat(PrelimFit)+labs(title = paste0(i,"_",j))
    ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_rhat.pdf"), plot=rhat_plot,device="pdf",width = 10, height = 8)
    
    ### effective sample size
    ess_plot<-stan_ess(PrelimFit)+labs(title = paste0(i,"_",j))### neff to total sample size
    ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_ess.pdf"), plot=ess_plot,device="pdf",width = 10, height = 8)
    
    #hist(summary(PrelimFit)$summary[,"n_eff"])
    # Next check the correlation among key model parameters and identify any
    #       divergent transitions
    parms<-nuts_params(PrelimFit)
    mcmc_plot<-mcmc_pairs(PrelimFit,pars = c("lambda","alpha_acam","alpha_avba","alpha_erbo","alpha_gitr","alpha_lomu","alpha_pler","alpha_taca","alpha_trwi"),max_treedepth = 15,np=parms)
    #pair_plot<-pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"),main=paste0(i,"_",j))
    ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_pairs.pdf"), plot=mcmc_plot,device="pdf",width = 10, height = 8)
    
    # Finally, check for autocorrelation in the posteriors of key model parameters
    ac_plot<-stan_ac(PrelimFit,pars = c("lambda","alpha_acam","alpha_avba","alpha_erbo","alpha_gitr","alpha_lomu","alpha_pler","alpha_taca","alpha_trwi"))+labs(title = paste0(i,"_",j))
    ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_acf.pdf"), plot=ac_plot,device="pdf",width = 10, height = 8)
  }
}
