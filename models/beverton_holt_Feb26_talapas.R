## Model script
library(ggsci)
library(tidyverse)
library(bayesplot)
library(rstan)
#library(here)
setwd("~/Desktop/career_r/data")
model.dat <- read.csv("career_model_data_Aug20_25.csv")%>%
  filter(!(is.na(num_seeds)))%>%##### removes samples that are in datasheet but weren't found during processing
  mutate_at((vars("weeds","sp_mean_lambda_mc")),~replace_na(.,0))%>%
  mutate_at("sp_sd_lambda_mc",~replace_na(.,1)) #### some Mono samples (and TRWI bkgs) have weeds as NA, the priors have been updated for each species

setwd("~/Desktop/career_repo/models")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

date <- Sys.Date()
species <- c("ACAM","AVBA","ERBO","GITR","LOMU","PLER","TACA") ## dropping TRWI
treatment<-c("D","A","AC","AG","ANG","ING","IN","DNG","IG","AN","DN","ACN","DG","I")
model.output <- list()
warnings <- list()
ptm<-proc.time()
for(i in species){
  for(j in treatment){
    #tic("model time")
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
    mc_mean_lambda <- unique(dat$sp_mean_lambda_mc) #### meagcomp priors for GITR, PLER, TACA, LOMU, ACAM, 
                                                    #priors from Hallett et al. 2019 for AVBA, ERBO, mean across all environments, no priors for TRWI
    mc_sd_lambda <- unique(dat$sp_sd_lambda_mc)
    
    ## stems data
    taca<-as.integer(dat$TACA)
    acam<-as.integer(dat$ACAM)
    avba<-as.integer(dat$AVBA)
    erbo<-as.integer(dat$ERBO)
    gitr<-as.integer(dat$GITR)
    lomu<-as.integer(dat$LOMU)
    pler<-as.integer(dat$PLER)
    weeds<-as.integer(dat$weeds)
    #trwi<-as.integer(dat$TRIW)
    
    # if (i== "TACA"){
    #   weeds <- as.integer(dat$total_weeds_n) 
    # } else if(i=="LOMU"){
    #   weeds <- as.integer(dat$total_weeds_n)
    # } else if(i=="AVBA"){
    #   weeds <- as.integer(dat$total_weeds_n) 
    # } else if(i == "GITR"){
    #   weeds<-as.integer(dat$total_weed_cover)
    # } else if(i== "PLER"){
    #   weeds<-as.integer(dat$total_weed_cover)
    # }else if (i=="ERBO"){
    #   weeds<-as.integer(dat$total_weed_cover)
    # #}else if (i=="TRWI"){
    #   #weeds<-as.integer(dat$total_weed_cover)
    # }
    # else if (i=="ACAM"){
    #   weeds<-as.integer(dat$neigh_weeds)
    # }
    
    ## make a vector of data inputs to model
    
    data_vec <- c("N","Fecundity","N_i","g_i","s_i","mc_mean_lambda","mc_sd_lambda","N_blocks","Blocks","weeds","acam","avba","erbo","gitr","lomu","pler","taca")
    
    
    ## create initials for epsilon and sigma
    initials <- list(epsilon=rep(1,N_blocks), sigma = 1,alpha_acam=0.1,alpha_avba=0.1,alpha_erbo=0.1,
                     alpha_gitr=0.1,alpha_lomu=0.1,alpha_pler=0.1,alpha_taca=0.1,alpha_weeds=0.1)
    
    initials1<- list(initials, initials, initials, initials)
    
    # Model ####
    print(paste("running",i,j,sep = "_"))
    PrelimFit <- stan(file="~/Desktop/career_repo/models/beverton_holt_negbi.stan", model_name="beverton_holt_negbi",
                      data = data_vec, init = initials1, iter = 25000, chains = 4, cores=4,
                      control = list(adapt_delta = 0.95, max_treedepth = 15)) 
    
    
    ## save model output
    save(PrelimFit, file = paste0("~/model_outputs/Feb-20/",i,"_",j,"_posteriors_BH_negbin", date ,".rdata"))

  }
}
#### if you get an error message and to check specific lines
stanmod <- readLines("~/Desktop/career_repo/models/BH_sparse_preliminary.stan")
stanmod <- stanmod[stanmod != ""]
stanmod[87] ### line number of the error message, the number given usually doesn't match the line of code in the model

###################
########## now extract posteriors
species2<-c("ACAM","AVBA","ERBO","GITR","LOMU","PLER","TACA")
species2<-c("LOMU")
#treatment2<-c("ACN")
treatment2<-c("D","A","AC","AG","ANG","ING","IN","DNG","IG","AN","DN","ACN","DG","I")

#c("D","A","AC","IN","AN","DN","I")
#"D","A","AC","AG","ANG","ING","IN","DNG","IG","AN",
for (i in species2){
  for (j in treatment2){
    if(file.exists(paste0("~/Desktop/talapas_model_outputs/Feb_09_26/",i,"_",j,"_posteriors_BH_negbin2026-02-09.rdata"))==TRUE){
      load(file=paste0("~/Desktop/talapas_model_outputs/Feb_09_26/",i,"_",j,"_posteriors_BH_negbin2026-02-09.rdata"))
    } else next
    print(paste("running",i,j,sep = "_"))
    
    # PrelimPosteriors <- extract(PrelimFit) ###
# 
#     ##### Diagnostic plots
#     # First check the distribution of Rhats and effective sample sizes
#     #hist(summary(PrelimFit)$summary[,"Rhat"])
#     rhat_plot<-stan_rhat(PrelimFit)+labs(title = paste0(i,"_",j))
#     ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_rhat.pdf"), plot=rhat_plot,device="pdf",width = 10, height = 8)
# 
#     ### effective sample size
#     ess_plot<-stan_ess(PrelimFit)+labs(title = paste0(i,"_",j))### neff to total sample size
#     ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_ess.pdf"), plot=ess_plot,device="pdf",width = 10, height = 8)
# 
#     #hist(summary(PrelimFit)$summary[,"n_eff"])
#     # Next check the correlation among key model parameters and identify any
#     #       divergent transitions
#     parms<-nuts_params(PrelimFit)
#     mcmc_plot<-mcmc_pairs(PrelimFit,pars = c("lambda","alpha_acam","alpha_avba","alpha_erbo","alpha_gitr","alpha_lomu","alpha_pler","alpha_taca","alpha_trwi"),max_treedepth = 15,np=parms)
#     #pair_plot<-pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"),main=paste0(i,"_",j))
#     ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_pairs.pdf"), plot=mcmc_plot,device="pdf",width = 10, height = 8)
# 
#     # Finally, check for autocorrelation in the posteriors of key model parameters
#     ac_plot<-stan_ac(PrelimFit,pars = c("lambda","alpha_acam","alpha_avba","alpha_erbo","alpha_gitr","alpha_lomu","alpha_pler","alpha_taca","alpha_trwi"))+labs(title = paste0(i,"_",j))
#     # ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_acf.pdf"), plot=ac_plot,device="pdf",width = 10, height = 8)
    
    output<-summary(PrelimFit)$summary%>%
      as.data.frame() %>% 
      mutate(variable = rownames(.)) %>% 
      select(variable, everything()) %>% 
      as.data.frame()
    
    ### neff is 40000, number of non warmup iters x chains
    #### divergence plots
    
    mod_diagnostics <- rstan::get_sampler_params(PrelimFit) %>% 
      set_names(1:4) %>% 
      map_df(as_data_frame,.id = 'chain') %>% 
      group_by(chain) %>% 
      mutate(iteration = 1:length(chain)) %>% 
      mutate(warmup = iteration <= 10000)
    
    # ## divergent chains
    # div<-mod_diagnostics %>%
    #   group_by(warmup, chain) %>%
    #   summarise(percent_divergent = mean(divergent__ >0)) %>%
    #   ggplot() +
    #   labs(title = paste0(i,"_",j))+
    #   theme_classic()+
    #   geom_col(aes(chain, percent_divergent, fill = warmup), position = 'dodge', color = 'black') +
    #   scale_y_continuous(labels = scales::percent, name = "% Divergent Runs")  +
    #   scale_fill_npg()
    # 
    # ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_div.pdf"), plot=div,device="pdf",width = 10, height = 8)

    #### trace
    trace1<-stan_trace(PrelimFit,pars = c("alpha_weeds","alpha_acam","alpha_avba","alpha_erbo","alpha_gitr","alpha_lomu","	
                            alpha_pler","alpha_taca"))+labs(title = paste0(i,"_",j))### trace plot, specify pars
    ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_trace1.pdf"), plot=trace1,device="pdf",width = 10, height = 8)
    
    if(sum(str_detect(output$variable, 'epsilon'))==4){
      trace2<-stan_trace(PrelimFit,pars = c("lambda","disp_dev","sigma","epsilon[1]","epsilon[2]","epsilon[3]","epsilon[4]"))+labs(title = paste0(i,"_",j))
    } else if (sum(str_detect(output$variable, 'epsilon'))==3){
      trace2<-stan_trace(PrelimFit,pars = c("lambda","disp_dev","sigma","epsilon[1]","epsilon[2]","epsilon[3]"))+labs(title = paste0(i,"_",j))
    }
    ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_trace2.pdf"), plot=trace2,device="pdf",width = 10, height = 8)
    
    ###### parameter estimates, bio and non bio
    dens_plots<-stan_dens(PrelimFit,pars = c("lambda","alpha_weeds","alpha_acam","alpha_avba","alpha_erbo","alpha_gitr","alpha_lomu","	
                            alpha_pler","alpha_taca"))+labs(title = paste0(i,"_",j))## density, specify pars
    ggsave(paste0("~/Desktop/talapas_model_outputs/diagnostic_figures/",i,"_",j,"_dens.pdf"), plot=dens_plots,device="pdf",width = 10, height = 8)
    print(paste("finished",i,j,sep = "_"))
   # rm(PrelimFit)
  }
}
