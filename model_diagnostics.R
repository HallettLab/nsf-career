####### diagnostic plots
library(readr)
library(tidyverse)
library(rstan)
library(bayesplot)
library(ggsci)
setwd("~/Desktop/career_repo/model_outputs/March_3_25")### deleted ACAM DN in file because it has no samples

#### ACAM_DN mod file does not contain samples, in mod data it is one ACAM_TRWI and two ACAM_ACAM
###### splitting up species due to memory issues, trying to save as plots are made instead of lists
species <- c("ACAM","ERBO","GITR","LOMU","PLER","TACA")
treatment<-c("D","A","AC","AG","ANG","ING","IN","DNG","IG","AN","DN","ACN","DG","I")
mod_list<-

for (i in species){
  for (j in treatment){
    if(file.exists(paste0(i,"_",j,"_posteriors_BH_negbin2025-03-03.rdata"))==TRUE){
      load(file=paste0(i,"_",j,"_posteriors_BH_negbin2025-03-03.rdata"))
    } else next
    
    output<-summary(PrelimFit)$summary%>%
      as.data.frame() %>% 
      mutate(variable = rownames(.)) %>% 
      select(variable, everything()) %>% 
      as.data.frame()
    
    mod_list[[i]][[j]]<-output
    
    ##### filter out rhat values > 1.1?
    
    ### neff is 40000, number of non warmup iters x chains
    #### divergence plots
    
    mod_diagnostics <- rstan::get_sampler_params(PrelimFit) %>% 
      set_names(1:4) %>% 
      map_df(as_data_frame,.id = 'chain') %>% 
      group_by(chain) %>% 
      mutate(iteration = 1:length(chain)) %>% 
      mutate(warmup = iteration <= 10000)
    
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
    
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_div.pdf"), plot=div,device="pdf",width = 10, height = 8)
    
    #### max treedepth
    tree<-mod_diagnostics %>% 
      ggplot(aes(iteration, treedepth__, color = chain)) + 
      geom_line() + 
      labs(title = paste0(i,"_",j))+
      geom_hline(aes(yintercept = 15), color = 'red') + 
      scale_color_locuszoom()+
      theme_classic()
    
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_tree.pdf"), plot=tree,device="pdf",width = 10, height = 8)
    
    ### stepsize
    step<-mod_diagnostics %>% 
      ggplot(aes(iteration, stepsize__, color = chain)) + 
      geom_line()  + 
      labs(title = paste0(i,"_",j))+
      theme_classic()+
      scale_color_locuszoom()
    
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_step.pdf"), plot=step,device="pdf",width = 10, height = 8)
    
    ###effective sample size, sum of the effectively independent sampling iterations across all chains, num chains x non warmup iters
    #output %>% 
    #ggplot(aes(n_eff)) + 
    #geom_histogram() + 
    #geom_vline(aes(xintercept = 40000), color = 'red')+
    #theme_classic()
    
    ess_plot<-stan_ess(PrelimFit)+labs(title = paste0(i,"_",j))### neff to total sample size
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_ess.pdf"), plot=ess_plot,device="pdf",width = 10, height = 8)
    
    ###rhat, tells whether chains have reached a stable posterior distribution
    #output %>% 
    #ggplot(aes(Rhat)) + 
    #geom_histogram() + 
    #geom_vline(aes(xintercept = 1.1), color = 'red')+
    #theme_classic()+
    #labs(title=paste0(i,"_",j))
    
    rhat_plot<-stan_rhat(PrelimFit)+labs(title = paste0(i,"_",j))### MCMC convergence
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_rhat.pdf"), plot=rhat_plot,device="pdf",width = 10, height = 8)
    
    ### autocorrelation
    ac_plot<-stan_ac(PrelimFit)+labs(title = paste0(i,"_",j))### autocorrelation
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_auto.pdf"), plot=ac_plot,device="pdf",width = 10, height = 8)
    
    ### mcmc se
    mcse_plot<-stan_mcse(PrelimFit)+labs(title = paste0(i,"_",j))### standard error
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_mcse.pdf"), plot=mcse_plot,device="pdf",width = 10, height = 8)
    
    #### trace
    trace1<-stan_trace(PrelimFit,pars = c("alpha_weeds_base","alpha_acam_base","alpha_avba_base","alpha_erbo_base","alpha_gitr_base","alpha_lomu_base","	
                            alpha_pler_base","alpha_taca_base"))+labs(title = paste0(i,"_",j))### trace plot, specify pars
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_trace1.pdf"), plot=trace1,device="pdf",width = 10, height = 8)
    
    if(sum(str_detect(output$variable, 'epsilon'))==4){
      trace2<-stan_trace(PrelimFit,pars = c("lambda_base","disp_dev","sigma","epsilon[1]","epsilon[2]","epsilon[3]","epsilon[4]"))+labs(title = paste0(i,"_",j))
    } else if (sum(str_detect(output$variable, 'epsilon'))==3){
      trace2<-stan_trace(PrelimFit,pars = c("lambda_base","disp_dev","sigma","epsilon[1]","epsilon[2]","epsilon[3]"))+labs(title = paste0(i,"_",j))
    }
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_trace2.pdf"), plot=trace2,device="pdf",width = 10, height = 8)
    
    ###pairs
    print("pairs")
    pairs<-mcmc_pairs(PrelimFit,pars = c("lambda_base","alpha_weeds_base","alpha_acam_base","alpha_avba_base","alpha_erbo_base","alpha_gitr_base","alpha_lomu_base","alpha_pler_base","alpha_taca_base"))
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_pairs.pdf"), plot=pairs,device="pdf",width = 10, height = 8)
    print("pairs 1")
    if(sum(str_detect(output$variable, 'epsilon'))==4){
      pairs2<-mcmc_pairs(PrelimFit,pars = c("lambda_base","disp_dev","sigma","epsilon[1]","epsilon[2]","epsilon[3]","epsilon[4]"))
    } else if (sum(str_detect(output$variable, 'epsilon'))==3){
      pairs2<-stan_trace(PrelimFit,pars = c("lambda_base","disp_dev","sigma","epsilon[1]","epsilon[2]","epsilon[3]"))
    }
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_pairs2.pdf"), plot=pairs2,device="pdf",width = 10, height = 8)
    print("pairs 2")
    ## parallel coordinates
    parcoord<-mcmc_parcoord(PrelimFit)+theme_classic()+labs(title = paste0(i,"_",j))
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_parcoord.pdf"), plot=parcoord,device="pdf",width = 10, height = 8)
    
    ###### parameter estimates, bio and non bio
    dens_plots<-stan_dens(PrelimFit,pars = c("lambda_base","alpha_weeds_base","alpha_acam_base","alpha_avba_base","alpha_erbo_base","alpha_gitr_base","alpha_lomu_base","	
                            alpha_pler_base","alpha_taca_base"))+labs(title = paste0(i,"_",j))## density, specify pars
    ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_dens.pdf"), plot=dens_plots,device="pdf",width = 10, height = 8)
    
    rm(PrelimFit)
    print(paste0(i,"_",j))
  }
}

############
########## diagnostic plots, trimmed down, for testing with erbo updated models, 03_20_2025
setwd("~/Desktop/test_mods")

############
species2 <- c("ACAM","ERBO","GITR","LOMU","PLER","TACA")
#treatment<-c("D","A","AC","AG","ANG","ING","IN","DNG","IG","AN","DN","ACN","DG","I")
treatment2<-c("D","A","I")

mod_list2<-list()  

for (i in species2){
  for (j in treatment2){
    if(file.exists(paste0(i,"_",j,"_posteriors_BH_erbo_update2025-03-20.rdata"))==TRUE){
      load(file=paste0(i,"_",j,"_posteriors_BH_erbo_update2025-03-20.rdata"))
    } else next
    
    output<-summary(PrelimFit)$summary%>%
      as.data.frame() %>% 
      mutate(variable = rownames(.)) %>% 
      select(variable, everything()) %>% 
      as.data.frame()
    
    mod_list2[[i]][[j]]<-output
    
    ##### filter out rhat values > 1.1?
    
    ### neff is 40000, number of non warmup iters x chains
    #### divergence plots
    
    mod_diagnostics <- rstan::get_sampler_params(PrelimFit) %>% 
      set_names(1:4) %>% 
      map_df(as_data_frame,.id = 'chain') %>% 
      group_by(chain) %>% 
      mutate(iteration = 1:length(chain)) %>% 
      mutate(warmup = iteration <= 1000)
    
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
    
    ggsave(paste0("~/Desktop/test_mods/",i,"_",j,"_div.pdf"), plot=div,device="pdf",width = 10, height = 8)
    
    #### max treedepth
    tree<-mod_diagnostics %>% 
      ggplot(aes(iteration, treedepth__, color = chain)) + 
      geom_line() + 
      labs(title = paste0(i,"_",j))+
      geom_hline(aes(yintercept = 15), color = 'red') + 
      scale_color_locuszoom()+
      theme_classic()
    
    ggsave(paste0("~/Desktop/test_mods/",i,"_",j,"_tree.pdf"), plot=tree,device="pdf",width = 10, height = 8)
    
    ### stepsize
    #step<-mod_diagnostics %>% 
    #ggplot(aes(iteration, stepsize__, color = chain)) + 
    #geom_line()  + 
    # labs(title = paste0(i,"_",j))+
    # theme_classic()+
    # scale_color_locuszoom()
    
    #ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_step.pdf"), plot=step,device="pdf",width = 10, height = 8)
    
    ###effective sample size, sum of the effectively independent sampling iterations across all chains, num chains x non warmup iters
    #output %>% 
    #ggplot(aes(n_eff)) + 
    #geom_histogram() + 
    #geom_vline(aes(xintercept = 40000), color = 'red')+
    #theme_classic()
    
    ess_plot<-stan_ess(PrelimFit)+labs(title = paste0(i,"_",j))### neff to total sample size
    ggsave(paste0("~/Desktop/test_mods/",i,"_",j,"_ess.pdf"), plot=ess_plot,device="pdf",width = 10, height = 8)
    
    ###rhat, tells whether chains have reached a stable posterior distribution
    #output %>% 
    #ggplot(aes(Rhat)) + 
    #geom_histogram() + 
    #geom_vline(aes(xintercept = 1.1), color = 'red')+
    #theme_classic()+
    #labs(title=paste0(i,"_",j))
    
    rhat_plot<-stan_rhat(PrelimFit)+labs(title = paste0(i,"_",j))### MCMC convergence
    ggsave(paste0("~/Desktop/test_mods/",i,"_",j,"_rhat.pdf"), plot=rhat_plot,device="pdf",width = 10, height = 8)
    
    ### autocorrelation
    #ac_plot<-stan_ac(PrelimFit)+labs(title = paste0(i,"_",j))### autocorrelation
    #ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_auto.pdf"), plot=ac_plot,device="pdf",width = 10, height = 8)
    
    ### mcmc se
    #mcse_plot<-stan_mcse(PrelimFit)+labs(title = paste0(i,"_",j))### standard error
    #ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_mcse.pdf"), plot=mcse_plot,device="pdf",width = 10, height = 8)
    
    #### trace
    trace1<-stan_trace(PrelimFit,pars = c("alpha_weeds","alpha_acam","alpha_avba","alpha_erbo","alpha_gitr","alpha_lomu","	
                            alpha_pler","alpha_taca"))+labs(title = paste0(i,"_",j))### trace plot, specify pars
    ggsave(paste0("~/Desktop/test_mods/",i,"_",j,"_trace1.pdf"), plot=trace1,device="pdf",width = 10, height = 8)
    
    if(sum(str_detect(output$variable, 'epsilon'))==4){
      trace2<-stan_trace(PrelimFit,pars = c("lambda_base","disp_dev","sigma","epsilon[1]","epsilon[2]","epsilon[3]","epsilon[4]"))+labs(title = paste0(i,"_",j))
    } else if (sum(str_detect(output$variable, 'epsilon'))==3){
      trace2<-stan_trace(PrelimFit,pars = c("lambda_base","disp_dev","sigma","epsilon[1]","epsilon[2]","epsilon[3]"))+labs(title = paste0(i,"_",j))
    }
    ggsave(paste0("~/Desktop/test_mods/",i,"_",j,"_trace2.pdf"), plot=trace2,device="pdf",width = 10, height = 8)
    
    ###pairs
    print("pairs")
    pairs<-mcmc_pairs(PrelimFit,pars = c("lambda_base","alpha_weeds","alpha_acam","alpha_avba","alpha_erbo","alpha_gitr","alpha_lomu","alpha_pler","alpha_taca"))
    ggsave(paste0("~/Desktop/test_mods/",i,"_",j,"_pairs.pdf"), plot=pairs,device="pdf",width = 10, height = 8)
    print("pairs 1")
    if(sum(str_detect(output$variable, 'epsilon'))==4){
      pairs2<-mcmc_pairs(PrelimFit,pars = c("lambda_base","disp_dev","sigma","epsilon[1]","epsilon[2]","epsilon[3]","epsilon[4]"))
    } else if (sum(str_detect(output$variable, 'epsilon'))==3){
      pairs2<-stan_trace(PrelimFit,pars = c("lambda_base","disp_dev","sigma","epsilon[1]","epsilon[2]","epsilon[3]"))
    }
    ggsave(paste0("~/Desktop/test_mods/",i,"_",j,"_pairs2.pdf"), plot=pairs2,device="pdf",width = 10, height = 8)
    print("pairs 2")
    ## parallel coordinates
    #parcoord<-mcmc_parcoord(PrelimFit)+theme_classic()+labs(title = paste0(i,"_",j))
    #ggsave(paste0("~/Desktop/career_repo/model_outputs/March_3_25/plots/",i,"_",j,"_parcoord.pdf"), plot=parcoord,device="pdf",width = 10, height = 8)
    
    ###### parameter estimates, bio and non bio
    dens_plots<-stan_dens(PrelimFit,pars = c("lambda_base","alpha_weeds","alpha_acam","alpha_avba","alpha_erbo","alpha_gitr","alpha_lomu","	
                            alpha_pler","alpha_taca"))+labs(title = paste0(i,"_",j))## density, specify pars
    ggsave(paste0("~/Desktop/test_mods/",i,"_",j,"_dens.pdf"), plot=dens_plots,device="pdf",width = 10, height = 8)
    
    rm(PrelimFit)
    print(paste0(i,"_",j))
  }
}

########using stan plot functions from r stan
#stan_diag(PrelimFit,information = "sample")+labs(title=paste0(i,"_",j))
#stan_par(PrelimFit,par=("sigma"))### single par lp, acceptance plots
#stan_scat(PrelimFit,pars = c("alpha_acam_base","alpha_taca_base"))##scatterplot, teo only

#### plots to show people, dens, plot, scat, 
### plots purely diagnostic: rhat, ac, sample, treedepth, divergence, ess, trace


##### looking at model input data
setwd("~/Desktop/career_r/data")

model.dat <- read.csv("career_model_data.csv")%>%
  filter(!(is.na(num_seeds)))##### removes samples that are in datasheet but weren't found during processing

ggplot(data=model.dat,aes(x=num_seeds))+
  geom_density()+
  facet_wrap(~phyto,scales = "free")+
  theme_classic()

ggplot(data=model.dat,aes(x=bg_n_indiv))+
  geom_histogram()+
  facet_wrap(~background,scales = "free")+
  theme_classic()

ggplot(data=model.dat,aes(x=ERBO))+
  geom_histogram()+
  facet_wrap(~phyto,scales = "free")+
  theme_classic()
