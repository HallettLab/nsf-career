######## figures for lab retreat
library(tidyverse)
library(data.table)
library(here)
rm(list = ls())

setwd("~/Desktop/career_r/data")

model.dat <- read.csv("career_model_data_Aug20_25.csv")%>%
  filter(!(is.na(num_seeds)))%>%##### removes samples that are in datasheet but weren't found during processing
  mutate_at((vars("weeds","sp_mean_lambda_mc")),~replace_na(.,0))%>%
  mutate_at("sp_sd_lambda_mc",~replace_na(.,1)) #### some Mono samples (and TRWI bkgs) have weeds as NA

setwd("~/Desktop/sparse_model_outputs/Aug22_negbi_final")

files = list.files(pattern="AVBA*")
dataset = do.call(rbind, lapply(files, fread))
rm(files)

file_names<-as.list(dir(pattern="AVBA_*"))
lapply(file_names,load,.GlobalEnv)
filenames <- list.files(file.path(getwd()), pattern='AVBA*', full.names=T)
results <- sapply(filenames, function(x) mget(load(x)), simplify = TRUE)

output<-summary(FinalFit)$summary%>%
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as.data.frame()
output2<-(FinalFit)%>%
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as.data.frame()
treatment<-c("D") #"A","AC","AG","ANG","ING","IN","DNG","IG","AN","DN","ACN","DG","I")

################ extract alphas
i<-"PLER"
j<-"IN"
mod<-paste0(i,"_",j)
rm(FinalFit, SpNames, N, S, Fecundity, SpMatrix, Inclusion_ij,
   tau0, slab_scale, slab_df, Intra, Post)
#for (i in species){
  #for (j in treatment){

load(file=paste0("~/Desktop/sparse_model_outputs/Aug22_negbi_final/",i,"_",j,"_posteriors_BH_sparse2025-08-22.rdata"))

Post <- rstan::extract(FinalFit)

GenericPlot <- array(NA, dim = c(3, N))
IntraPlot <- array(NA, dim = c(3, N))

# Create an object to hold the alpha values for any species identified by the sparse
#  modeling approach as deviating from the generic value
Inclusion_ij # if all 0 none included
SpecificPlot <- array(NA, dim = c(2, 3, N))

# Now calculate the alpha_e,i,j values for intraspecific alphas (IntraPlot), any heterospecifics that
#  affect the focal species as an average heterospecific neighbor (GenericPlot), and for any species
#  that differ from that generic heterospecific effect (SpecificPlot)

  # Calculate the posterior distribution for the relevant alpha_i,j
  GenericPost <- exp(Post$alpha_generic)
  #Sp1Post <- exp(Post$alpha_generic + Post$alpha_hat_ij[,9])
  #Sp2Post <- exp(Post$alpha_generic + Post$alpha_hat_ij[,6])
  IntraPost <- exp(Post$alpha_intra)
  # Now store the mean and 95% CI for these posteriors
  for(i in 1:N){
  GenericPlot[1,i] <- mean(GenericPost)
  GenericPlot[2:3,i] <- HDInterval::hdi(GenericPost)
  IntraPlot[1,i] <- mean(IntraPost)
  IntraPlot[2:3,i] <- HDInterval::hdi(IntraPost)
  #SpecificPlot[1,1,i] <- mean(Sp1Post)
  #SpecificPlot[1,2:3,i] <- HDInterval::hdi(Sp1Post)
  # SpecificPlot[2,1,i] <- mean(Sp2Post)
  # SpecificPlot[2,2:3,i] <- HDInterval::hdi(Sp2Post)
  }
  PLER_IN_dat <- data.frame(t(GenericPlot), t(IntraPlot))%>% #t(SpecificPlot[1,,]), t(SpecificPlot[2,,]
    rename(gen.mean = X1, gen.low = X2, gen.high = X3,
           intra.mean = X1.1, intra.low = X2.1, intra.high = X3.1)
           
            #sp1.mean = X1.2, sp1.low = X2.2, sp1.high = X3.2,
  PLER_IN_dat$mod<-mod                                                                  #sp2.mean = X1.3, sp2.low = X2.3, sp2.high = X3.3)
  PLER_IN_dat<-PLER_IN_dat%>%pivot_longer(!mod,names_to = "sp", values_to = "alpha")%>%
    unique()%>%
    separate(sp, into = c("species", "value")) %>%
    pivot_wider(names_from = "value", values_from = "alpha")
  PLER_IN_dat$trt<-j 

TACA_AN_dat$species[TACA_AN_dat$species=="sp1"]="weeds"

PLER_dat<-rbind(PLER_A_dat,PLER_AN_dat,PLER_AC_dat,PLER_ACN_dat,PLER_D_dat,PLER_DN_dat,PLER_I_dat,PLER_IN_dat)
#PLER_dat$phyto<-"PLER"
rm(FinalFit, SpNames, N, S, Fecundity, SpMatrix, Inclusion_ij,
   tau0, slab_scale, slab_df, Intra, Post)
#output2$Generic_Post<-exp(Post$alpha_generic)
#output2$IntraPost <- exp(Post$alpha_intra)

#dat_name<-paste0(i,"_",j)
#dat_name<- list(Post = Post, SpNames = SpNames, N = N, S = S, Fecundity = Fecundity,
            #SpMatrix = SpMatrix, Inclusion_ij = Inclusion_ij,
            #Intra = Intra)
alpha_dat<-list(AVBA_dat,GITR_dat,PLER_dat,TACA_dat)%>%reduce(full_join)
alpha_dat$species[alpha_dat$species=="sp1"]="weeds"

alpha_dat[1,6]="A"
alpha_dat[2,6]="A"
alpha_dat[33,6]="A"
alpha_dat[34,6]="A"
alpha_dat[49,6]="AN"
alpha_dat[50,6]="AN"
alpha_dat[51,6]="AN"
################ getting lambda values
lambda_list<-list()
rm(i,j)
lambda_list[[i]][[j]]<-list()
species<-"AVBA"
treatment<-c("A","AC","AN","ACN","D","DN","I","IN")

for (i in species){
  for (j in treatment){
mod<-paste0(i,"_",j)
rm(FinalFit, SpNames, N, S, Fecundity, SpMatrix, Inclusion_ij,
   tau0, slab_scale, slab_df, Intra, Post)

load(file=paste0("~/Desktop/sparse_model_outputs/Aug22_negbi_final/",i,"_",j,"_posteriors_BH_sparse2025-08-22.rdata"))

Post <- rstan::extract(FinalFit)
N<-(length(Post$lambdas)/2)
LambdaPlotVals <- array(NA, dim = c(N,1)) #metric (mean, lwr, upr)
for(k in 1:N){
    LPost <- exp(Post$lambdas[k,1] + Post$lambdas[k,2])
    LambdaPlotVals[k,1] <- (LPost)
}
cfi<-hdi(LambdaPlotVals) 
lambda_dat<-data.frame(mean_lambda=mean(LambdaPlotVals))
lambda_dat$upper<-cfi[2,1]
lambda_dat$lower<-cfi[1,1]
lambda_dat$mod<-mod
lambda_dat$trt<-j
lambda_dat$phyto<-i
lambda_list[[i]][[j]]<-lambda_dat
  }
}
PLER_lambda_dat<-rbind(PLER_A_lambda,PLER_AN_lambda,PLER_AC_lambda,PLER_ACN_lambda,PLER_D_lambda,PLER_DN_lambda,PLER_I_lambda,PLER_IN_lambda)
reduce(full_join,lambda_list)
AVBA_lambda_dat<-lambda_list[[1]] %>% reduce(full_join)
df_list<-list(AVBA_lambda_dat,GITR_lambda_dat,PLER_lambda_dat,TACA_lambda_dat)
lambda_data<-df_list %>% reduce(full_join)

########### FIGURES
jake_theme <- function () { 
  theme_linedraw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.border = element_blank(),
      panel.background  = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA), 
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0, vjust = 0.5, 
                                margin = margin(b = 0.2, unit = "cm")),
      plot.subtitle = element_text(size = 10, hjust = 0, vjust = 0.5, 
                                   margin = margin(b = 0.4, unit = "cm")),
      plot.caption = element_text(size = 7, hjust = 1, face = "italic", 
                                  margin = margin(t = 0.1, unit = "cm")),
      axis.text = element_text(size = 14),
      
      
    )
}
###### alphas
species<-c("AVBA","GITR","PLER")
for (s in species){
  fig_dat<-filter(alpha_dat,phyto==s)
ggplot(fig_dat,aes(x=mod,y=mean))+
  geom_point(aes(color=factor(species)),size=5)+
  geom_errorbar(aes(ymin=mean-low,ymax=mean+high),width=0.25,alpha=0.66,color="gray7")+
  guides(color=guide_legend(title="Interaction Coefficient"))+
  labs(x="",y="Mean Alpha")+
  scale_color_manual(values = c("dodgerblue", "seagreen"),
                     labels = c("Generic Interspecific", "Intraspecific"))+
  scale_x_discrete(guide = guide_axis(angle=45))+
  jake_theme()

ggsave(paste0("~/Desktop/sparse_model_outputs/figures/",s,"_alphas_fig.pdf"),plot=last_plot(),device=cairo_pdf,width = 10,height = 6)
}

fig_dat<-filter(alpha_dat,phyto=="TACA")
ggplot(fig_dat,aes(x=mod,y=mean))+
  geom_point(aes(color=factor(species)),size=5)+
  geom_errorbar(aes(ymin=mean-low,ymax=mean+high),width=0.25,alpha=0.66,color="gray7")+
  guides(color=guide_legend(title="Interaction Coefficient"))+
  labs(x="",y="Mean Alpha")+
  scale_color_manual(values = c("dodgerblue", "seagreen","salmon"),
                     labels = c("Generic Interspecific", "Intraspecific","Weeds"))+
  scale_x_discrete(guide = guide_axis(angle=45))+
  jake_theme()
ggsave(paste0("~/Desktop/sparse_model_outputs/figures/TACA_alphas_fig.pdf"),plot=last_plot(),device=cairo_pdf,width = 10,height = 6)

#########  by treatment not species

ggplot(alpha_dat,aes(x=trt,y=mean))+
  geom_point(aes(color=phyto,shape=species),size=5)+
  geom_errorbar(aes(ymin=mean-low,ymax=mean+high),width=0.25,alpha=0.66,color="gray7")+
  guides(color=guide_legend(title="Species"))+
  labs(x="",y="Mean Alpha") 

  scale_color_manual(values = c("dodgerblue", "seagreen","red"),
                     labels = c("Generic Interspecific", "Intraspecific","Weeds"))+
  jake_theme()


######### lambdas
ggplot(lambda_data,aes(x=trt,y=mean_lambda))+
  geom_point(aes(color=phyto),size=5)+
  geom_errorbar(aes(ymin=mean_lambda-lower,ymax=mean_lambda+upper),width=0.25,alpha=0.66,color="gray7")+
  guides(color=guide_legend(title="Species"))+
  labs(x="",y="Mean Lambda")+ 
scale_color_manual(values = c("goldenrod", "darkorchid","indianred4","palegreen4"),
                   labels = c("Avena barbata", "Gilia tricolor","Plantago erecta","Taeniatherum caput-medusae"))+
  jake_theme()+
  theme(legend.text = element_text(face="italic"))
ggsave(paste0("~/Desktop/sparse_model_outputs/figures/all_lambdas.pdf"),plot=last_plot(),device=cairo_pdf,width = 10,height = 6)

### by species
species<-c("AVBA","GITR","PLER","TACA")
#for (s in species){
  fig_dat<-filter(lambda_data,phyto=="AVBA")
  ggplot(lambda_data,aes(x=trt,y=mean_lambda))+
    facet_wrap(~phyto,scales="free")+
    geom_point(aes(color=phyto),size=5)+
    geom_errorbar(aes(ymin=mean_lambda-lower,ymax=mean_lambda+upper),width=0.25,alpha=0.66,color="gray7")+
    guides(color=guide_legend(title="Species"))+
    labs(x="",y="Mean Lambda (Seed Output)")+ 
    scale_color_manual(values = c("goldenrod", "darkorchid","indianred4","palegreen4"),
                       labels = c("Avena barbata", "Gilia tricolor","Plantago erecta","Taeniatherum caput-medusae"))+
    jake_theme()+
    scale_x_discrete(guide = guide_axis(angle=45))+
    theme(legend.text = element_text(face="italic"))+
    #theme(strip.text = element_text(color="black",size=16))+
    theme(strip.text = element_blank())+
    theme(strip.background = element_rect(fill="white",color="white"))
    
ggsave(paste0("~/Desktop/sparse_model_outputs/figures/all_lambdas_fig2.pdf"),plot=last_plot(),device=cairo_pdf,width = 10,height = 6)

####### Lambdas and alphas
la_dat<-left_join(alpha_dat,lambda_data)
ggplot(la_dat,aes(x=mean,y=mean_lambda))+
  facet_wrap(~phyto,scales="free",labeller=labeller(phyto=c("AVBA"="Avena barbata","GITR"="Gilia tricolor","PLER"="Plantago erecta","TACA"="Taeniatherum caput-medusae")))+
  geom_point(aes(color=species),size=5)+
  geom_smooth(se=FALSE)+
  labs(x="Mean Alpha",y="Mean Lambda (Seed Output)")+ 
  guides(color=guide_legend(title="Species"))+
  scale_color_manual(values = c("dodgerblue", "seagreen","red"),
                     labels = c("Generic Interspecific", "Intraspecific","Weeds"))+
  jake_theme()+
  theme(strip.text = element_text(color="black",face="italic",size=16))+
  theme(strip.background = element_rect(fill="white",color="white"))

ggsave(paste0("~/Desktop/sparse_model_outputs/figures/lambdas_alphas.pdf"),plot=last_plot(),device=cairo_pdf,width = 10,height = 6)
