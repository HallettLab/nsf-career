###### sparse modeling script, adapted from Weiss-Lehman et al.2022
# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript

rm(list = ls())
library(beepr)
library(tidyverse)
library(rstan)
library(here)
library(HDInterval)
library(tictoc)
library(bayesplot)
library(ggsci)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("~/Desktop/career_r/data")

model.dat <- read.csv("career_model_data_Aug20_25.csv")%>%
  filter(!(is.na(num_seeds)))%>%##### removes samples that are in datasheet but weren't found during processing
  mutate_at((vars("weeds","sp_mean_lambda_mc")),~replace_na(.,0))%>%
  mutate_at("sp_sd_lambda_mc",~replace_na(.,1)) #### some Mono samples (and TRWI bkgs) have weeds as NA
####### set priors for AVBA and ERBO as 0,1 normal
  setwd("~/Desktop/career_repo")

date <- Sys.Date()
#species <- c("ACAM","AVBA","ERBO","GITR","LOMU","PLER","TACA","TRWI")
species <- c("LOMU")

#treatment<-c("D","A","AC","AG","ANG","ING","IN","DNG","IG","AN","DN","ACN","DG","I")
treatment<-c("IN","DNG","IG","AN","DN","ACN","DG","I")### testing sparse model with TACA in three trts

#model.output <- list()
#warnings <- list()
ptm<-proc.time()
for(i in species){
  for(j in treatment){
    
    tic("model time")
    dat <- subset(model.dat, phyto == i & treatment==j)
    FocalSpecies<-i
    ## create vectors of the various data inputs
    Fecundity <- as.integer(round(dat$num_seeds)) ## seeds out, is it the number of seeds or the average?
    N <- as.integer(length(Fecundity)) ## number of observations
    if (N!=0){print(paste(i,j,N,sep = "_"))}
    else next ####### I think this logic loop runs through itself entirely before moving on
    
    AllSpAbunds <- select(dat,TRWI:ACAM,weeds)
    AllSpNames <-  names(select(dat,TRWI:ACAM,weeds))
    SpTotals <- colSums(AllSpAbunds)
    SpToKeep <- SpTotals > 0
    S <- sum(SpToKeep)
    SpMatrix <- matrix(NA, nrow = N, ncol = S)
    f <- 1
    
    for(s in 1:ncol(AllSpAbunds)){
      if(SpToKeep[s] == 1){
        SpMatrix[,f] <- AllSpAbunds[,s]
        f <- f + 1
      }
    }
    SpNames <- AllSpNames[SpToKeep]
    Intra <- ifelse(SpNames == FocalSpecies, 1, 0)
    
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
    
# Set the parameters defining the regularized horseshoe prior, as described in
#       the "Incorporating sparsity-inducing priors" section of the manuscript.
tau0 <- 1
slab_scale <- sqrt(2)
slab_df <- 4

data_vec <- c("N","N_i","S","SpMatrix","Intra","tau0","slab_scale","slab_df", "Fecundity","g_i","s_i","N_blocks","Blocks","mc_mean_lambda","mc_sd_lambda")

initials <- list(epsilon=rep(1,N_blocks), sigma = 1)
                 #lambdas=c(5,5),alpha_generic_tilde=2,
                # alpha_intra_tilde=2,alpha_hat_ij_tilde=rep(0,S),c2_tilde=1.25,
                 #local_shrinkage_ij= rep(5, S), tau_tilde = 15,
                # alpha_intra_tilde = 1)

initials1<- list(initials,initials,initials,initials)

# Model ####
print(paste("running",i,j,sep = "_"))
PrelimFit <- stan(file="~/Desktop/career_repo/models/BH_sparse_preliminary.stan", model_name="sparse_BH_negbi",
                  data = data_vec, iter = 5000,chains = 4, cores=4, init = initials1, 
                  control = list(adapt_delta = 0.95, max_treedepth = 15)) 


## save model output
save(PrelimFit, file = paste0("~/Desktop/sparse_model_outputs/Aug21_negbi/",i,"_",j,"_posteriors_BH_sparse", date ,".rdata"))
#beep(2)
toc() #### should output how long each run takes
  }
}
proc.time()-ptm
#### for checking what line error messages correspond to because stan can't count
stanmod <- readLines("~/Desktop/career_repo/models/BH_sparse_preliminary.stan")
stanmod <- stanmod[stanmod != ""]
stanmod[87]

########## now extract posteriors
species2<-c("LOMU")
treatment2<-c("D","A","AC","AG","ANG","ING","IN","DNG","IG","AN","DN","ACN","DG","I")

#c("D","A","AC","IN","AN","DN","I")
#"D","A","AC","AG","ANG","ING","IN","DNG","IG","AN",
for (i in species2){
  for (j in treatment2){
    if(file.exists(paste0("~/Desktop/sparse_model_outputs/Aug21_negbi/",i,"_",j,"_posteriors_BH_sparse2025-08-21.rdata"))==TRUE){
      load(file=paste0("~/Desktop/sparse_model_outputs/Aug21_negbi/",i,"_",j,"_posteriors_BH_sparse2025-08-21.rdata"))
    } else next
    
PrelimPosteriors <- extract(PrelimFit) ###

##### Diagnostic plots
# First check the distribution of Rhats and effective sample sizes
#hist(summary(PrelimFit)$summary[,"Rhat"])
rhat_plot<-stan_rhat(PrelimFit)+labs(title = paste0(i,"_",j))
ggsave(paste0("~/Desktop/sparse_model_outputs/diagnostic_figures/",i,"_",j,"_rhat.pdf"), plot=rhat_plot,device="pdf",width = 10, height = 8)

### effective sample size
ess_plot<-stan_ess(PrelimFit)+labs(title = paste0(i,"_",j))### neff to total sample size
ggsave(paste0("~/Desktop/sparse_model_outputs/diagnostic_figures/",i,"_",j,"_ess.pdf"), plot=ess_plot,device="pdf",width = 10, height = 8)

#hist(summary(PrelimFit)$summary[,"n_eff"])
# Next check the correlation among key model parameters and identify any
#       divergent transitions
parms<-nuts_params(PrelimFit)
mcmc_pairs(PrelimFit,pars = c("lambdas[1]","lambdas[2]","alpha_generic","alpha_intra"),max_treedepth = 15,np=parms)
#pair_plot<-pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"),main=paste0(i,"_",j))
ggsave(paste0("~/Desktop/sparse_model_outputs/diagnostic_figures/",i,"_",j,"_pairs.pdf"), plot=pair_plot,device="pdf",width = 10, height = 8)

# Finally, check for autocorrelation in the posteriors of key model parameters
ac_plot<-stan_ac(PrelimFit,pars = c("lambdas","alpha_generic","alpha_intra"))+labs(title = paste0(i,"_",j))
ggsave(paste0("~/Desktop/sparse_model_outputs/diagnostic_figures/",i,"_",j,"_acf.pdf"), plot=ac_plot,device="pdf",width = 10, height = 8)
  }
}


#### If the diagnostic plots don't reveal any problems with the model fit, now
#       move on to determining which parameters warrant inclusion in the final
#       model (i.e. the data pulled their posteriors away from 0). The final model
#       will then be run with only these species-specific parameters, but without
#       the regularized horseshoe priors.
species3<-c("LOMU","TACA","GITR","PLER")
#treatment3<-c("D","A","I")
#treatment<-c("D","A","AC","IN","AN","DN","ACN","DG","I")
treatment3<-c("D","A","AC","AG","ANG","ING","IN","DNG","IG","AN","ACN","DG","I")

ptm<-proc.time()

for (i in species3){
  for (j in treatment3){
    tic("model time")
    if(file.exists(paste0("~/Desktop/sparse_model_outputs/Aug21_negbi/",i,"_",j,"_posteriors_BH_sparse2025-08-21.rdata"))==TRUE){
      load(file=paste0("~/Desktop/sparse_model_outputs/Aug21_negbi/",i,"_",j,"_posteriors_BH_sparse2025-08-21.rdata"))
    } else next
    
PrelimPosteriors <- extract(PrelimFit) ###

dat <- subset(model.dat, phyto == i & treatment==j)
    FocalSpecies<-i
    ## create vectors of the various data inputs
    Fecundity <- as.integer(round(dat$num_seeds)) ## seeds out, is it the number of seeds or the average?
    N <- as.integer(length(Fecundity)) ## number of observations
    if (N!=0){print(paste(i,j,N,sep = "_"))}
    else next ####### I think this logic loop runs through itself entirely before moving on
    
    AllSpAbunds <- select(dat,TRWI:ACAM,weeds)
    AllSpNames <-  names(select(dat,TRWI:ACAM,weeds))
    SpTotals <- colSums(AllSpAbunds)
    SpToKeep <- SpTotals > 0
    S <- sum(SpToKeep)
    SpMatrix <- matrix(NA, nrow = N, ncol = S)
    f <- 1
    
    for(s in 1:ncol(AllSpAbunds)){
      if(SpToKeep[s] == 1){
        SpMatrix[,f] <- AllSpAbunds[,s]
        f <- f + 1
      }
    }
    SpNames <- AllSpNames[SpToKeep]
    Intra <- ifelse(SpNames == FocalSpecies, 1, 0)
    
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

    ############ parms to include
    Inclusion_ij <- rep(0,S)
    #Inclusion_eij <- matrix(data = 0, nrow = 2, ncol = S)
    IntLevel <- 0.5 #0.5 usually, 0.75 for Waitzia, shade
    for(s in 1:S){
      Ints_ij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
      # Ints_eij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_eij[,i,s], credMass = IntLevel)
      if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
        Inclusion_ij[s] <- 1
      }
    }
    
    # Set the parameters defining the regularized horseshoe prior, as described in
   #       the "Incorporating sparsity-inducing priors" section of the manuscript.
tau0 <- 1
slab_scale <- sqrt(2)
slab_df <- 4
#sum(Inclusion_ij)

############# Final fit

data_vec <- c("N","N_i","S","SpMatrix","Intra","Fecundity","g_i","s_i","N_blocks","Blocks",
              "mc_mean_lambda","mc_sd_lambda","Inclusion_ij")

initials <- list(epsilon=rep(1,N_blocks), sigma = 1)

initials1<- list(initials,initials,initials)

print(paste("running",i,j,sep = "_"))
FinalFit <- stan(file="~/Desktop/career_repo/models/BH_sparse_final.stan", model_name="sparse_BH_negbi",
                  data = data_vec, iter = 10000,chains = 3, cores=8, init = initials1, 
                  control = list(adapt_delta = 0.95, max_treedepth = 15)) 


## save model output
FileName <- paste0("~/Desktop/sparse_model_outputs/Aug22_negbi_final/",i,"_",j,"_posteriors_BH_sparse", date ,".rdata")
save(FinalFit, SpNames, N, S, Fecundity, SpMatrix, Inclusion_ij,
     tau0, slab_scale, slab_df, Intra, file = FileName)
#beep(2)
toc() #### should output how long each run takes
  }
}
proc.time()-ptm
#beep(8)

################
species4<-c("GITR","TACA","LOMU")
treatment4<-c("D","A","AC","AG","ANG","ING","IN","DNG","IG","AN","ACN","DG","I","DN")

# Final Diagnostic figures
for (i in species4){
  for (j in treatment4){
    if(file.exists(paste0("~/Desktop/sparse_model_outputs/Aug22_negbi_final/",i,"_",j,"_posteriors_BH_sparse2025-08-22.rdata"))==TRUE){
      load(file=paste0("~/Desktop/sparse_model_outputs/Aug22_negbi_final/",i,"_",j,"_posteriors_BH_sparse2025-08-22.rdata"))
    } else next
    
    FinalPosteriors <- extract(FinalFit) ###
    
    ##### Diagnostic plots
    # First check the distribution of Rhats and effective sample sizes
    #hist(summary(PrelimFit)$summary[,"Rhat"])
    rhat_plot<-stan_rhat(FinalFit)+labs(title = paste0(i,"_",j))
    ggsave(paste0("~/Desktop/sparse_model_outputs/diagnostic_figures_final/",i,"_",j,"_rhat.pdf"), plot=rhat_plot,device="pdf",width = 10, height = 8)
    
    ### effective sample size
    ess_plot<-stan_ess(FinalFit)+labs(title = paste0(i,"_",j))### neff to total sample size
    ggsave(paste0("~/Desktop/sparse_model_outputs/diagnostic_figures_final/",i,"_",j,"_ess.pdf"), plot=ess_plot,device="pdf",width = 10, height = 8)
    
    #hist(summary(PrelimFit)$summary[,"n_eff"])
    # Next check the correlation among key model parameters and identify any
    #       divergent transitions
    #parms<-nuts_params(FinalFit)
    #mcmc_pairs(FinalFit,pars = c("lambdas[1]","lambdas[2]","alpha_generic","alpha_intra"),max_treedepth = 15,np=parms)
    #pair_plot<-pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"),main=paste0(i,"_",j))
    #ggsave(paste0("~/Desktop/sparse_model_outputs/diagnostic_figures_final/",i,"_",j,"_pairs.pdf"), plot=pair_plot,device="pdf",width = 10, height = 8)
    
    # Finally, check for autocorrelation in the posteriors of key model parameters
    ac_plot<-stan_ac(FinalFit,pars = c("lambdas","alpha_generic","alpha_intra"))+labs(title = paste0(i,"_",j))
    ggsave(paste0("~/Desktop/sparse_model_outputs/diagnostic_figures_final/",i,"_",j,"_acf.pdf"), plot=ac_plot,device="pdf",width = 10, height = 8)

mod_diagnostics <- rstan::get_sampler_params(FinalFit) %>% 
  set_names(1:3) %>% 
  map_df(as_data_frame,.id = 'chain') %>% 
  group_by(chain) %>% 
  mutate(iteration = 1:length(chain)) %>% 
  mutate(warmup = iteration <= 1500)

### divergent chains
div<-mod_diagnostics %>% 
  group_by(warmup, chain) %>% 
  summarise(percent_divergent = mean(divergent__ >0)) %>% 
  ggplot() +
  labs(title = paste0(i,"_",j))+
  theme_classic()+
  geom_col(aes(chain, percent_divergent, fill = warmup), position = 'dodge', color = 'black') + 
  scale_y_continuous(labels = scales::percent, name = "% Divergent Runs")  + 
  scale_fill_npg()

ggsave(paste0("~/Desktop/sparse_model_outputs/diagnostic_figures_final/",i,"_",j,"_div.pdf"), plot=div,device="pdf",width = 10, height = 8)

#### max treedepth
tree<-mod_diagnostics %>% 
  ggplot(aes(iteration, treedepth__, color = chain)) + 
  geom_line() + 
  labs(title = paste0(i,"_",j))+
  geom_hline(aes(yintercept = 15), color = 'red') + 
  scale_color_locuszoom()+
  theme_classic()

ggsave(paste0("~/Desktop/sparse_model_outputs/diagnostic_figures_final/",i,"_",j,"_tree.pdf"), plot=tree,device="pdf",width = 10, height = 8)
  }
}
beep(3)

