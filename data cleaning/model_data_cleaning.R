##### Data cleaning for stan models
library(readr)
library(data.table)
library(tidyverse)
#######
setwd("~/Desktop/career_r/data") #haven't decided the best way to upload data so it's easily accessible so just local to me for now
########## germ and seed data
germ_dat<-read.csv("20220218_Germination-Data_full.csv")

mean_germ<-germ_dat%>%group_by(Species)%>%
  summarise(mean_germ=mean(p.germ))%>%
  ungroup()%>%
  filter(Species %in% c("ACMAME","AVEBAR","ELYCAP","EROBOT","FESPER","GILTRI","PLAERE","TRIWIL"))%>%
  mutate(Species=ifelse(Species=="ACMAME","ACAM",Species),
         Species=ifelse(Species=="AVEBAR","AVBA",Species),
         Species=ifelse(Species=="ELYCAP","TACA",Species),
         Species=ifelse(Species=="FESPER","LOMU",Species),
         Species=ifelse(Species=="EROBOT","ERBO",Species),
         Species=ifelse(Species=="GILTRI","GITR",Species),
         Species=ifelse(Species=="PLAERE","PLER",Species),
         Species=ifelse(Species=="TRIWIL","TRWI",Species),
  )%>%
  rename(phyto=Species)
### write_csv(mean_germ,"CAREER_mean_germ.csv)

### seed surv data
seed_dat<-read_csv("CAREER_mean_seed_surv.csv") ### erbo surv taken from Rice 1985 paper

########### species data
taca<-read_csv("taca_phyto_biomass.csv")
acam<-read_csv("acam_phyto_biomass.csv")
pler<-read_csv("pler_phyto_biomass.csv")
gitr<-read_csv("gitr_phyto_biomass.csv")
lomu<-read_csv("lomu_phyto_biomass.csv")
erbo<-read_csv("erbo_phyto_counts.csv")
avba<-read_csv("avba_phyto_counts.csv")
trwi<-read_csv("trwi_phyto_biomass.csv")
bl_key<-read_csv("block_treatments.csv")

############# clean up taca
taca<-taca %>%left_join(bl_key,by="block")%>%
  #filter(biomass!="NA")%>%
  #mutate_at(vars("biomass"),as.numeric)%>%
  mutate(avg_biomass=biomass/ph_n_indiv)%>%
  mutate(num_seeds=biomass*59.82549283)%>%
  mutate(avg_seeds=avg_biomass*59.82549283)%>%
  select(phyto,block,plot,sub,background,density,ph_n_indiv,biomass,avg_biomass,num_seeds,avg_seeds,treatment,water,grazing,nitrogen)

################
taca_bkg_n<-read_csv("taca_phyto_final.csv")

taca_zeros<-taca%>%filter(!(block==4 & plot ==21 & sub==5 & biomass ==0.596))%>%
  left_join(taca_bkg_n)

taca_mono<-read_csv("taca_mono_final.csv")

all_taca<-left_join(taca_zeros,taca_mono,by=join_by(block,plot,sub,background,density))%>%
  select(-c("ph_n_indiv.y","ph_date.y","ph_initials.y","ph_notes.y"))%>%
  rename(ph_n_indiv=ph_n_indiv.x,ph_date=ph_date.x,ph_notes=ph_notes.x,ph_initials=ph_initials.x)%>%
  mutate(background=ifelse(background=="TACA_mono","Mono",background))

############
taca_bkg_col<-read_csv("taca_bkg_collections.csv")%>%
  filter(bg_n_indiv_collected !=0)%>%
  mutate_at(vars("avg_biomass"),as.numeric)%>%
  mutate(num_seeds=biomass * 59.82549283,
         avg_seeds=avg_biomass * 59.82549283)%>%
  rename(phyto=species,ph_n_indiv=bg_n_indiv_collected)%>%
  select(-c("bg_collection_notes","processing_notes","scale_id","date","bg_date_collected","bg_initials"))

########### getting estimates of bkg plots for bkg collections
bkg_plots<-read_csv("bkg_estimates_all.csv")
na_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(is.na(bg_n_indiv))%>%
  filter(is.na(ph_date) & is.na(ph_initials))

taca_bkg_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(background!="LOMU_mono" & background!="BRHO" & background !="Hirtum")%>%
  filter(background=="TACA")%>%
  filter(!(is.na(bg_n_indiv)))%>%
  group_by(block, plot, background)%>%
  summarise(avg_bkg_n=mean(bg_n_indiv))%>%
  rename(bg_n_indiv=avg_bkg_n)

taca_bkg_plots$bg_n_indiv<-round(taca_bkg_plots$bg_n_indiv) ### avg num of taca in bkg plots from phyto collections
taca_bkg_plots<-left_join(taca_bkg_col,taca_bkg_plots)
taca_bkg_plots<-left_join(taca_bkg_plots,bl_key)


all_taca<-all_taca%>%mutate_at(vars("bg_n_indiv"),as.numeric)  
taca_all<-full_join(all_taca,taca_bkg_plots)
taca_fil<-taca_all%>%filter(ph_n_indiv!=0)%>%
  filter(background!="Mono")
taca_mono_dat<-taca_all%>%filter(ph_n_indiv!=0)%>%
  filter(background=="Mono")
  
################# phyto data without monos
mod_dat_taca<-taca_fil%>%mutate(phyto2=phyto,
                            background2=background,
                            background3=background,
                            bg_n_indiv2=bg_n_indiv,
                            treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("TRWI":"TACA")),~replace_na(.,0))%>%
  mutate_at(vars(c("TRWI":"TACA")),as.numeric)%>%
  mutate_at(vars("wd_forb_indiv"),~replace_na(.,0))%>%
  mutate(total_weeds_n=wd_grass_indiv+wd_forb_indiv,
         trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("total_weeds_n"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
 # filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("biomass","avg_biomass","water","grazing","nitrogen","ph_unique","ph_date","ph_initials","ph_notes","bg_date","bg_initials","bg_notes","seed_notes","unique.ID","wd_grass_indiv","wd_forb_indiv","wd_grass_cover","wd_forb_cover","bare_cover","cover_notes"))
mod_dat_taca<-left_join(mod_dat_taca,mean_germ)
mod_dat_taca<-left_join(mod_dat_taca,seed_dat) #### gets background and phyto collections in, no monos, are those only used for setting lambda priors?
#############
########## taca bkgs don't have weed counts so is set to 0, can estimate like bkg plots?
#### taca mono data, run separate models to get posteriors in each environment and then use those as priors for big model run?
taca_mono_dat2<-taca_mono_dat%>%#mutate(phyto2=phyto,
                                #background2=background,
                                #background3=background,
                                #bg_n_indiv2=bg_n_indiv,
                                #treatment2=treatment)%>%
 # pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
 # unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  #mutate_at(vars("Mono"),~replace_na(.,0))%>%
  #mutate_at(vars("Mono"),as.numeric)%>%
  #mutate_at(vars("wd_forb_indiv"),~replace_na(.,0))%>%
  mutate(total_weed_cover=wd_grass_cover+wd_forb_cover,
         trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
  # filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("bg_n_indiv","biomass","avg_biomass","water","grazing","nitrogen","ph_unique","ph_date","ph_initials","ph_notes","bg_date","bg_initials","bg_notes","seed_notes","unique.ID","wd_grass_indiv","wd_forb_indiv","wd_grass_cover","wd_forb_cover","bare_cover","cover_notes"))
taca_mono_dat2<-left_join(taca_mono_dat2,mean_germ)
taca_mono_dat2<-left_join(taca_mono_dat2,seed_dat)


################## Pler data clean
#pler<-read_csv("pler_phyto_biomass.csv") #has bg indivs in it and weed cover data

pler<-pler%>%left_join(bl_key,by="block")%>%
  #filter(biomass!="NA")%>%
  mutate(avg_biomass=biomass/ph_n_indiv)%>%
  mutate(num_seeds=biomass*346.8854513)%>%
  mutate(avg_seeds=avg_biomass*346.8854513)

########### getting estimates of bkg plots for bkg collections
bkg_plots<-read_csv("bkg_estimates_all.csv")
na_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(is.na(bg_n_indiv))%>%
  filter(is.na(ph_date) & is.na(ph_initials))

pler_bkg_col<-read_csv("pler_bkg_mass.csv")%>%
  filter(bg_n_indiv_collected !=0)%>%
  mutate_at(vars("avg_biomass"),as.numeric)%>%
  mutate(num_seeds=biomass *346.8854513 ,
         avg_seeds=avg_biomass *346.8854513 )%>%
  rename(phyto=species,ph_n_indiv=bg_n_indiv_collected)%>%
  select(-c("bg_collection_notes","processing_notes","scale_id","date","bg_date_collected","bg_initials","other_notes","seed_notes"))


pler_bkg_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(background!="LOMU_mono" & background!="BRHO" & background !="Hirtum")%>%
  filter(background=="PLER")%>%
  filter(!(is.na(bg_n_indiv)))%>%
  group_by(block, plot, background)%>%
  summarise(avg_bkg_n=mean(bg_n_indiv))%>%
  rename(bg_n_indiv=avg_bkg_n)

pler_bkg_plots$bg_n_indiv<-round(pler_bkg_plots$bg_n_indiv) ### avg num of taca in bkg plots from phyto collections
pler_bkg_plots<-left_join(pler_bkg_col,pler_bkg_plots)
pler_bkg_plots<-left_join(pler_bkg_plots,bl_key)

pler_all<-full_join(pler,pler_bkg_plots)
  
pler_fil<-pler_all%>%filter(ph_n_indiv!=0)%>%
  mutate(background=ifelse(background=="PLER_mono","Mono",background))%>%
  filter(background!="Mono" & background!="BRHO" & background !="Hirtum")

pler_mono_dat<-pler_all%>%filter(ph_n_indiv!=0)%>%
  filter(background=="Mono")
##### need to add in 0 for mono bkg indiv (weeds?) and estimated bkg densities for L and H density bkg collections
######## model data including 0 indivs and all mono indivs
mod_dat_pler<-pler_fil%>%mutate(phyto2=phyto,
                            background2=background,
                            background3=background,
                            bg_n_indiv2=bg_n_indiv,
                            treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("TRWI":"PLER")),~replace_na(.,0))%>%
  mutate_at(vars(c("TRWI":"PLER")),as.numeric)%>%
  mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>%
  mutate(trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("biomass","avg_biomass","water","grazing","nitrogen","phyto_unique","phyto_date","initials_ph","phyto_notes","background_date","initials_bg","background_notes","seed_notes","scale_id","disturbance_cover","bare_cover","weed_notes","processing_notes","other_notes","date"))
mod_dat_pler<-left_join(mod_dat_pler,mean_germ)
mod_dat_pler<-left_join(mod_dat_pler,seed_dat) 
########### weed cover for pler bkgs doesn't exist, can estimate from collections like bkgs?

######### pler mono data
pler_mono_dat2<-pler_mono_dat%>%#mutate(phyto2=phyto,
  #background2=background,
  #background3=background,
  #bg_n_indiv2=bg_n_indiv,
  #treatment2=treatment)%>%
  # pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  # unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  #mutate_at(vars("Mono"),~replace_na(.,0))%>%
  #mutate_at(vars("Mono"),as.numeric)%>%
  #mutate_at(vars("wd_forb_indiv"),~replace_na(.,0))%>%
  mutate(trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
  # filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("biomass","avg_biomass","water","grazing","nitrogen","phyto_unique","phyto_date","initials_ph","phyto_notes","background_date","initials_bg","background_notes","seed_notes","scale_id","disturbance_cover","bare_cover","weed_notes","processing_notes","other_notes","date"))
pler_mono_dat2<-left_join(pler_mono_dat2,mean_germ)
pler_mono_dat2<-left_join(pler_mono_dat2,seed_dat)
########## ACAM
acam$ph_n_indiv<-as.numeric(acam$ph_n_indiv)

acam<-acam%>%left_join(bl_key,by="block")%>%
  #filter(biomass!="NA")%>%
  mutate(avg_biomass=biomass/ph_n_indiv)%>%
  mutate(num_seeds=biomass*34.93414486)%>%
  mutate(avg_seeds=avg_biomass*34.93414486)

acam<-select(acam,phyto,block,plot,sub,background,density,ph_n_indiv,biomass,avg_biomass,num_seeds,avg_seeds,treatment,water,grazing,nitrogen)

acam_bkg_n<-read_csv("acam_phyto_final.csv")
acam_bkg_n$ph_n_indiv<-as.numeric(acam_bkg_n$ph_n_indiv)

acam_zeros<-left_join(acam,acam_bkg_n,relationship = "many-to-many")

acam_mono<-read_csv("acam_mono_final.csv")

all_acam<-left_join(acam_zeros,acam_mono,by=join_by(block,plot,sub,background,density))%>%
  select(-c("ph_n_indiv.y","ph_date.y","ph_initials.y","ph_notes.y"))%>%
  rename(ph_n_indiv=ph_n_indiv.x,ph_date=ph_date.x,ph_notes=ph_notes.x,ph_initials=ph_initials.x)%>%
  mutate(background=ifelse(background=="ACAM_mono","Mono",background))

############
acam_bkg_col<-read_csv("acam_bkg_collections.csv")%>%
  filter(bg_n_indiv_collected !=0)%>%
  mutate_at(vars("avg_biomass"),as.numeric)%>%
  mutate(num_seeds=biomass * 34.93414486,
         avg_seeds=avg_biomass * 34.93414486)%>%
  rename(phyto=species,ph_n_indiv=bg_n_indiv_collected)%>%
  select(-c("bg_collection_notes","processing_notes","scale_id","date","bg_date_collected","bg_initials"))

########### getting estimates of bkg plots for bkg collections
bkg_plots<-read_csv("bkg_estimates_all.csv")
na_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(is.na(bg_n_indiv))%>%
  filter(is.na(ph_date) & is.na(ph_initials))

acam_bkg_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(background!="LOMU_mono" & background!="BRHO" & background !="Hirtum")%>%
  filter(background=="ACAM")%>%
  filter(!(is.na(bg_n_indiv)))%>%
  group_by(block, plot, background)%>%
  summarise(avg_bkg_n=mean(bg_n_indiv))%>%
  rename(bg_n_indiv=avg_bkg_n)

acam_bkg_plots$bg_n_indiv<-round(acam_bkg_plots$bg_n_indiv) ### avg num of taca in bkg plots from phyto collections
acam_bkg_plots<-left_join(acam_bkg_col,acam_bkg_plots)
acam_bkg_plots<-left_join(acam_bkg_plots,bl_key)


acam_all<-full_join(all_acam,acam_bkg_plots)

acam_fil<-acam_all%>%filter(ph_n_indiv!=0)%>%
  filter(background!="Mono")

acam_mono_dat<-acam_all%>%filter(ph_n_indiv!=0)%>%
  filter(background=="Mono")

#################
mod_dat_acam<-acam_fil%>%mutate(phyto2=phyto,
                                background2=background,
                                background3=background,
                                bg_n_indiv2=bg_n_indiv,
                                treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("TRWI":"ACAM")),~replace_na(.,0))%>%
  mutate_at(vars(c("TRWI":"ACAM")),as.numeric)%>%
  mutate_at(vars(c("neighbor_grass_cover","neighbor_forb_cover")),~replace_na(.,0))%>%
  mutate(neigh_weeds=neighbor_grass_cover+neighbor_forb_cover,
         trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("neigh_weeds"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
  filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("biomass","avg_biomass","water","grazing","nitrogen","ph_date","ph_initials","ph_notes","bare_cover.x","bare_cover.y","litter_cover","bg_date","bg_initials","bg_notes","seed_notes","unique.ID","wd_grass_cover","wd_forb_cover","cover_notes"))
mod_dat_acam<-left_join(mod_dat_acam,mean_germ)
mod_dat_acam<-left_join(mod_dat_acam,seed_dat) #### gets background and phyto collections in, no monos, are those only used for setting lambda priors?

########## ACAM mono data
acam_mono_dat2<-acam_mono_dat%>%#mutate(phyto2=phyto,
  #background2=background,
  #background3=background,
  #bg_n_indiv2=bg_n_indiv,
  #treatment2=treatment)%>%
  # pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  # unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  #mutate_at(vars("Mono"),~replace_na(.,0))%>%
  #mutate_at(vars("Mono"),as.numeric)%>%
  #mutate_at(vars("wd_forb_indiv"),~replace_na(.,0))%>%
  mutate(total_weed_cover=wd_grass_cover+wd_forb_cover,
         trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
  # filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("biomass","avg_biomass","water","grazing","nitrogen","ph_date","ph_initials","ph_notes","bare_cover.x","bare_cover.y","litter_cover","bg_date","bg_initials","bg_notes","seed_notes","unique.ID","wd_grass_cover","wd_forb_cover","cover_notes","neighbor_forb_cover","neighbor_grass_cover","forb_notes","grass_notes","...10"))
acam_mono_dat2<-left_join(acam_mono_dat2,mean_germ)
acam_mono_dat2<-left_join(acam_mono_dat2,seed_dat)

########### LOMU
lomu$ph_n_indiv<-as.numeric(lomu$ph_n_indiv)

lomu<-lomu%>%left_join(bl_key,by="block")%>%
  #filter(biomass!="NA")%>%
  mutate(avg_biomass=biomass/ph_n_indiv)%>%
  mutate(num_seeds=biomass*192.008604)%>%
  mutate(avg_seeds=avg_biomass*192.008604)%>%
  select(phyto,block,plot,sub,background,density,ph_n_indiv,biomass,avg_biomass,num_seeds,avg_seeds,treatment,water,grazing,nitrogen)

################
lomu_bkg_n<-read_csv("lomu_phyto_final.csv")
lomu_bkg_n$ph_n_indiv<-as.numeric(lomu_bkg_n$ph_n_indiv)

lomu_zeros<-left_join(lomu,lomu_bkg_n)

lomu_mono<-read_csv("lomu_mono_final.csv")

all_lomu<-left_join(lomu_zeros,lomu_mono,by=join_by(block,plot,sub,background,density))%>%
  select(-c("ph_n_indiv.y","ph_date.y","ph_initials.y","ph_notes.y"))%>%
  rename(ph_n_indiv=ph_n_indiv.x,ph_date=ph_date.x,ph_notes=ph_notes.x,ph_initials=ph_initials.x)%>%
  mutate(background=ifelse(background=="LOMU_mono","Mono",background))

############
lomu_bkg_col<-read_csv("lomu_bkg_collections.csv")%>%
  filter(bg_n_indiv_collected !=0)%>%
  mutate_at(vars("avg_biomass"),as.numeric)%>%
  mutate(num_seeds=biomass * 192.008604,
         avg_seeds=avg_biomass * 192.008604)%>%
  rename(phyto=species,ph_n_indiv=bg_n_indiv_collected)%>%
  select(-c("bg_collection_notes","processing_notes","scale_id","date","bg_date_collected","bg_initials"))
########### getting estimates of bkg plots for bkg collections
bkg_plots<-read_csv("bkg_estimates_all.csv")
na_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(is.na(bg_n_indiv))%>%
  filter(is.na(ph_date) & is.na(ph_initials))

lomu_bkg_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(background!="LOMU_mono" & background!="BRHO" & background !="Hirtum")%>%
  filter(background=="LOMU")%>%
  filter(!(is.na(bg_n_indiv)))%>%
  group_by(block, plot, background)%>%
  summarise(avg_bkg_n=mean(bg_n_indiv))%>%
  rename(bg_n_indiv=avg_bkg_n)

lomu_bkg_plots$bg_n_indiv<-round(lomu_bkg_plots$bg_n_indiv) ### avg num of taca in bkg plots from phyto collections
lomu_bkg_plots<-left_join(lomu_bkg_col,lomu_bkg_plots)
lomu_bkg_plots<-left_join(lomu_bkg_plots,bl_key)


all_lomu<-all_lomu%>%mutate_at(vars("bg_n_indiv"),as.numeric) 

lomu_all<-full_join(all_lomu,lomu_bkg_plots)
  
lomu_fil<-lomu_all%>%filter(ph_n_indiv!=0)%>%
  filter(background!="Mono")

lomu_mono_dat<-lomu_all%>%filter(ph_n_indiv!=0)%>%
  filter(background=="Mono")
#################
mod_dat_lomu<-lomu_fil%>%mutate(phyto2=phyto,
                                background2=background,
                                background3=background,
                                bg_n_indiv2=bg_n_indiv,
                                treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("TRWI":"LOMU")),~replace_na(.,0))%>%
  mutate_at(vars(c("TRWI":"LOMU")),as.numeric)%>%
  mutate_at(vars(c("wd_forb_indiv","wd_grass_indiv")),as.numeric)%>%
  mutate_at(vars("wd_forb_indiv"),~replace_na(.,0))%>%
  mutate(total_weeds_n=wd_grass_indiv+wd_forb_indiv,
         trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("total_weeds_n"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
  filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("biomass","avg_biomass","water","grazing","nitrogen","ph_date","ph_initials","ph_notes","bg_date","bg_initials","bg_notes","seed_notes","unique.ID","wd_grass_indiv","wd_forb_indiv","wd_grass_cover","wd_forb_cover","bare_cover","cover_notes"))
mod_dat_lomu<-left_join(mod_dat_lomu,mean_germ)
mod_dat_lomu<-left_join(mod_dat_lomu,seed_dat) #### gets background and phyto collections in, no monos, are those only used for setting lambda priors?


########### lomu mono data
lomu_mono_dat2<-lomu_mono_dat%>%#mutate(phyto2=phyto,
  #background2=background,
  #background3=background,
  #bg_n_indiv2=bg_n_indiv,
  #treatment2=treatment)%>%
  # pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  # unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  #mutate_at(vars("Mono"),~replace_na(.,0))%>%
  #mutate_at(vars("Mono"),as.numeric)%>%
  #mutate_at(vars("wd_forb_indiv"),~replace_na(.,0))%>%
  mutate(total_weed_cover=wd_grass_cover+wd_forb_cover,
         trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
  # filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("biomass","avg_biomass","water","grazing","nitrogen","ph_date","ph_initials","ph_notes","bg_date","bg_initials","bg_notes","seed_notes","unique.ID","wd_grass_indiv","wd_forb_indiv","wd_grass_cover","wd_forb_cover","bare_cover","cover_notes"))
lomu_mono_dat2<-left_join(lomu_mono_dat2,mean_germ)
lomu_mono_dat2<-left_join(lomu_mono_dat2,seed_dat)

########## GITR
##################
#gitr<-read_csv("gitr_phyto_biomass.csv") #has bg indivs in it and weed cover data

gitr<-gitr%>%left_join(bl_key,by="block")%>%
  #filter(biomass!="NA")%>%
  mutate(avg_biomass=biomass/ph_n_indiv)%>%
  mutate(num_seeds=biomass*55.5053883)%>%
  mutate(avg_seeds=avg_biomass*55.5053883)

########### getting estimates of bkg plots for bkg collections
bkg_plots<-read_csv("bkg_estimates_all.csv")
na_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(is.na(bg_n_indiv))%>%
  filter(is.na(ph_date) & is.na(ph_initials))

gitr_bkg_col<-read_csv("gitr_bkg_collections.csv")%>%
  filter(bg_n_indiv!=0)%>%
  mutate_at(vars("avg_biomass"),as.numeric)%>%
  mutate(num_seeds=biomass *55.5053883 ,
         avg_seeds=avg_biomass *55.5053883 )%>%
  rename(phyto=species,ph_n_indiv=bg_n_indiv)%>%
  select(-c("bg_collection_notes","processing_notes","scale_id","date","bg_date_collected","bg_initials"))


gitr_bkg_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(background!="LOMU_mono" & background!="BRHO" & background !="Hirtum")%>%
  filter(background=="GITR")%>%
  filter(!(is.na(bg_n_indiv)))%>%
  group_by(block, plot, background)%>%
  summarise(avg_bkg_n=mean(bg_n_indiv))%>%
  rename(bg_n_indiv=avg_bkg_n)

gitr_bkg_plots$bg_n_indiv<-round(gitr_bkg_plots$bg_n_indiv) ### avg num of taca in bkg plots from phyto collections
gitr_bkg_plots<-left_join(gitr_bkg_col,gitr_bkg_plots)
gitr_bkg_plots<-left_join(gitr_bkg_plots,bl_key)

gitr_all<-full_join(gitr,gitr_bkg_plots)

gitr_fil<-gitr_all%>%filter(ph_n_indiv!=0)%>%
  mutate(background=ifelse(background=="GITR_mono","Mono",background))%>%
  filter(background!="Mono" & background!="BRHO" & background !="Hirtum")

gitr_mono_dat<-gitr_all%>%filter(ph_n_indiv!=0)%>%
  filter(background=="Mono")
##### need to add in 0 for mono bkg indiv (weeds?) and estimated bkg densities for L and H density bkg collections
######## model data including 0 indivs and all mono indivs
#which(duplicated(gitr_fil))
mod_dat_gitr<-gitr_fil%>%mutate(phyto2=phyto,
                                background2=background,
                                background3=background,
                                bg_n_indiv2=bg_n_indiv,
                                treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("PLER":"GITR")),~replace_na(.,0))%>%
  mutate_at(vars(c("PLER":"GITR")),as.numeric)%>%
  mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>%
  mutate(trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("sub_total_n","neighbor_total_n","biomass","avg_biomass","water","grazing","nitrogen","phyto_unique","phyto_date","initials_ph","phyto_notes","background_date","initials_bg","background_notes","scale_id","disturbance_cover","bare_cover","weed_notes","processing_notes","other_notes","date"))
mod_dat_gitr<-left_join(mod_dat_gitr,mean_germ)
mod_dat_gitr<-left_join(mod_dat_gitr,seed_dat) 

###### GITR Mono
gitr_mono_dat2<-gitr_mono_dat%>%#mutate(phyto2=phyto,
  #background2=background,
  #background3=background,
  #bg_n_indiv2=bg_n_indiv,
  #treatment2=treatment)%>%
  # pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  # unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  #mutate_at(vars("Mono"),~replace_na(.,0))%>%
  #mutate_at(vars("Mono"),as.numeric)%>%
  #mutate_at(vars("wd_forb_indiv"),~replace_na(.,0))%>%
  mutate(trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
  # filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("sub_total_n","neighbor_total_n","biomass","avg_biomass","water","grazing","nitrogen","phyto_unique","phyto_date","initials_ph","phyto_notes","background_date","initials_bg","background_notes","scale_id","disturbance_cover","bare_cover","weed_notes","processing_notes","other_notes","date"))
gitr_mono_dat2<-left_join(gitr_mono_dat2,mean_germ)
gitr_mono_dat2<-left_join(gitr_mono_dat2,seed_dat)


#####################
########### ERBO
erbo<-erbo%>%left_join(bl_key,by="block")%>%
#filter(biomass!="NA")%>%
mutate(total_flowers=n_flowers_mature+n_flowers_immature,
        avg_flowers=total_flowers/ph_n_indiv,
        num_seeds=total_flowers*5,
        avg_seeds=avg_flowers*5,
       avg_biomass=biomass/ph_n_indiv)


########### getting estimates of bkg plots for bkg collections
bkg_plots<-read_csv("bkg_estimates_all.csv")
na_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(is.na(bg_n_indiv))%>%
  filter(is.na(ph_date) & is.na(ph_initials))

erbo_bkg_col<-read_csv("erbo_bkg_counts.csv")%>%
  filter(bg_n_indiv_collected !=0)%>%
  mutate(total_flowers=n_flowers_mature+n_flowers_immature,
         num_seeds=total_flowers *5 ,
         avg_flowers=total_flowers/bg_n_indiv_collected,
         avg_seeds=avg_flowers *5,
         avg_biomass=biomass/bg_n_indiv_collected)%>%
  rename(phyto=species,ph_n_indiv=bg_n_indiv_collected)%>%
  select(-c("bg_notes","processing_notes","scale_id","bg_date","bg_initials","initials","processing_notes"))


erbo_bkg_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(background!="LOMU_mono" & background!="BRHO" & background !="Hirtum")%>%
  filter(background=="ERBO")%>%
  filter(!(is.na(bg_n_indiv)))%>%
  group_by(block, plot, background)%>%
  summarise(avg_bkg_n=mean(bg_n_indiv))%>%
  rename(bg_n_indiv=avg_bkg_n)

erbo_bkg_plots$bg_n_indiv<-round(erbo_bkg_plots$bg_n_indiv) ### avg num of taca in bkg plots from phyto collections
erbo_bkg_plots<-left_join(erbo_bkg_col,erbo_bkg_plots)
erbo_bkg_plots<-left_join(erbo_bkg_plots,bl_key)

erbo_all<-full_join(erbo,erbo_bkg_plots)

erbo_fil<-erbo_all%>%filter(ph_n_indiv!=0)%>%
  mutate(background=ifelse(background=="ERBO_mono","Mono",background))%>%
  filter(background!="Mono" & background!="BRHO" & background !="Hirtum")

erbo_mono_dat<-erbo_all%>%filter(ph_n_indiv!=0)%>%
  filter(background=="Mono")

##### need to add in 0 for mono bkg indiv (weeds?) and estimated bkg densities for L and H density bkg collections
######## model data including 0 indivs and all mono indivs
mod_dat_erbo<-erbo_fil%>%mutate(phyto2=phyto,
                                background2=background,
                                background3=background,
                                bg_n_indiv2=bg_n_indiv,
                                treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("TRWI":"ERBO")),~replace_na(.,0))%>%
  mutate_at(vars(c("TRWI":"ERBO")),as.numeric)%>%
  mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>%
  mutate(trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("biomass","water","grazing","nitrogen","phyto_unique","phyto_date","initials_ph","phyto_notes","background_date","Initials_bg","background_notes","scale_id","disturbance_cover","bare_cover","weed_notes","processing_notes","initials"))
mod_dat_erbo<-left_join(mod_dat_erbo,mean_germ)
mod_dat_erbo<-left_join(mod_dat_erbo,seed_dat) 


###### ERBO Mono dat
erbo_mono_dat2<-erbo_mono_dat%>%#mutate(phyto2=phyto,
  #background2=background,
  #background3=background,
  #bg_n_indiv2=bg_n_indiv,
  #treatment2=treatment)%>%
  # pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  # unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  #mutate_at(vars("Mono"),~replace_na(.,0))%>%
  #mutate_at(vars("Mono"),as.numeric)%>%
  #mutate_at(vars("wd_forb_indiv"),~replace_na(.,0))%>%
  mutate(trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
  # filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("biomass","avg_biomass","water","grazing","nitrogen","phyto_unique","phyto_date","initials_ph","phyto_notes","background_date","Initials_bg","background_notes","scale_id","disturbance_cover","bare_cover","weed_notes","processing_notes","initials"))
erbo_mono_dat2<-left_join(erbo_mono_dat2,mean_germ)
erbo_mono_dat2<-left_join(erbo_mono_dat2,seed_dat)


#AVBA
avba<-avba %>%left_join(bl_key,by="block")%>%
  #filter(biomass!="NA")%>%
  mutate_at(vars(c("ph_n_indiv","num_glumes")),as.numeric)%>%
  mutate(avg_n_glumes=num_glumes/ph_n_indiv,
                      num_seeds=num_glumes*2,
                      avg_seeds=num_seeds/ph_n_indiv)%>%
  select(phyto,block,plot,sub,background,density,ph_n_indiv,num_seeds,avg_seeds,treatment,water,grazing,nitrogen)

################
avba_bkg_n<-read_csv("avba_phyto_final.csv")%>%
  mutate_at(vars("ph_n_indiv"),as.numeric)

avba_zeros<-left_join(avba,avba_bkg_n)

avba_mono<-read_csv("avba_mono_final.csv")

all_avba<-left_join(avba_zeros,avba_mono,by=join_by(block,plot,sub,background,density,ph_notes))%>%
  select(-c("ph_n_indiv.y","ph_date.y","ph_initials.y"))%>%
  rename(ph_n_indiv=ph_n_indiv.x,ph_date=ph_date.x,ph_initials=ph_initials.x)%>%
  mutate(background=ifelse(background=="AVBA_mono","Mono",background))

############
avba_bkg_col<-read_csv("avba_bkg_counts.csv")%>%
  filter(bg_n_indiv_collected !=0)%>%
  #mutate_at(vars("avg_biomass"),as.numeric)%>%
  mutate(avg_n_glumes=num_glumes/bg_n_indiv_collected)%>%
  mutate(num_seeds=num_glumes*2)%>%
  mutate(avg_seeds=num_seeds/bg_n_indiv_collected)%>%
  rename(phyto=species,ph_n_indiv=bg_n_indiv_collected)%>%
  select(-c("total_seeds","bg_collection_notes","processing_notes","scale_id","date","bg_date_collected","bg_initials","num_glumes","avg_n_glumes"))

########### getting estimates of bkg plots for bkg collections
bkg_plots<-read_csv("bkg_estimates_all.csv")
na_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(is.na(bg_n_indiv))%>%
  filter(is.na(ph_date) & is.na(ph_initials))

avba_bkg_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(background!="LOMU_mono" & background!="BRHO" & background !="Hirtum")%>%
  filter(background=="AVBA")%>%
  filter(!(is.na(bg_n_indiv)))%>%
  group_by(block, plot, background)%>%
  summarise(avg_bkg_n=mean(bg_n_indiv))%>%
  rename(bg_n_indiv=avg_bkg_n)

avba_bkg_plots$bg_n_indiv<-round(avba_bkg_plots$bg_n_indiv) ### avg num of taca in bkg plots from phyto collections
avba_bkg_plots<-left_join(avba_bkg_col,taca_bkg_plots)
avba_bkg_plots<-left_join(avba_bkg_plots,bl_key)


all_avba<-all_avba%>%mutate_at(vars(c("bg_n_indiv","ph_n_indiv")),as.numeric)  
avba_all<-full_join(all_avba,avba_bkg_plots)

avba_fil<-avba_all%>%filter(ph_n_indiv!=0)%>%
  filter(background!="Mono")
avba_mono_dat<-avba_all%>%filter(ph_n_indiv!=0)%>%
  filter(background=="Mono")

################# phyto data without monos
mod_dat_avba<-avba_fil%>%mutate(phyto2=phyto,
                                background2=background,
                                background3=background,
                                bg_n_indiv2=bg_n_indiv,
                                treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("PLER":"AVBA")),~replace_na(.,0))%>%
  mutate_at(vars(c("PLER":"AVBA")),as.numeric)%>%
  mutate_at(vars("wd_forb_indiv"),~replace_na(.,0))%>%
  mutate(total_weeds_n=wd_grass_indiv+wd_forb_indiv,
         trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("total_weeds_n"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
  # filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("initials","biomass","avg_biomass","water","grazing","nitrogen","ph_unique","ph_date","ph_initials","ph_notes","bg_date","bg_initials","bg_notes","seed_notes","unique.ID","wd_grass_indiv","wd_forb_indiv","wd_grass_cover","wd_forb_cover","bare_cover","cover_notes"))
mod_dat_avba<-left_join(mod_dat_avba,mean_germ)
mod_dat_avba<-left_join(mod_dat_avba,seed_dat) #### gets background and phyto collections in, no monos, are those only used for setting lambda priors?

############# avba mono
avba_mono_dat2<-avba_mono_dat%>%#mutate(phyto2=phyto,
  #background2=background,
  #background3=background,
  #bg_n_indiv2=bg_n_indiv,
  #treatment2=treatment)%>%
  # pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  # unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  #mutate_at(vars("Mono"),~replace_na(.,0))%>%
  #mutate_at(vars("Mono"),as.numeric)%>%
  #mutate_at(vars("wd_forb_indiv"),~replace_na(.,0))%>%
  mutate(total_weed_cover=wd_grass_cover+wd_forb_cover,
         trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
  # filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("bg_n_indiv","biomass","avg_biomass","water","grazing","nitrogen","ph_unique","ph_date","ph_initials","ph_notes","bg_date","bg_initials","bg_notes","seed_notes","unique.ID","wd_grass_indiv","wd_forb_indiv","wd_grass_cover","wd_forb_cover","bare_cover","cover_notes"))
avba_mono_dat2<-left_join(avba_mono_dat2,mean_germ)
avba_mono_dat2<-left_join(avba_mono_dat2,seed_dat)

########## TRWI
trwi<-trwi%>%left_join(bl_key,by="block")%>%
  #filter(biomass!="NA")%>%
  mutate(avg_biomass=biomass/ph_n_indiv)%>%
  mutate(num_seeds=biomass*235.9913615)%>%
  mutate(avg_seeds=avg_biomass*235.9913615)

trwi_bkg_n<-read_csv("trwi_phyto_final.csv")
trwi_bkg_n$ph_n_indiv<-as.numeric(trwi_bkg_n$ph_n_indiv)

trwi_zeros<-left_join(trwi,trwi_bkg_n)

trwi_mono<-read_csv("trwi_mono_final.csv")

all_trwi<-left_join(trwi_zeros,trwi_mono)
  #by=join_by(species,block,plot,background,density,phyto_unique,bg_n_indiv,ph_n_indiv,ph_notes,phyto_date))%>%
  #select(-c("initials_ph.y","initials._ph.x","ph_notes.y"))%>%
  #rename(ph_n_indiv=ph_n_indiv.x,ph_date=ph_date.x,ph_notes=ph_notes.x,ph_initials=ph_initials.x)%>%
  #mutate(background=ifelse(background=="TRWI_mono","Mono",background))


########### getting estimates of bkg plots for bkg collections
bkg_plots<-read_csv("bkg_estimates_all.csv")
na_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(is.na(bg_n_indiv))%>%
  filter(is.na(ph_date) & is.na(ph_initials))

trwi_bkg_col<-read_csv("trwi_bkg_biomass.csv")%>%
  filter(bg_n_indiv !=0)%>%
  #mutate_at(vars("avg_biomass"),as.numeric)%>%
  mutate(num_seeds=biomass *235.9913615 ,
         avg_biomass=biomass/bg_n_indiv,
         avg_seeds=avg_biomass *235.9913615)%>%
  rename(phyto=species,ph_n_indiv=bg_n_indiv)%>%
  select(-c("background_notes","processing_notes","scale_id","date","background_date","bg_initials"))


trwi_bkg_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(background!="LOMU_mono" & background!="BRHO" & background !="Hirtum")%>%
  filter(background=="TRWI")%>%
  filter(!(is.na(bg_n_indiv)))%>%
  group_by(block, plot, background)%>%
  summarise(avg_bkg_n=mean(bg_n_indiv))%>%
  rename(bg_n_indiv=avg_bkg_n,)

trwi_bkg_plots$bg_n_indiv<-round(trwi_bkg_plots$bg_n_indiv) ### avg num of taca in bkg plots from phyto collections
trwi_bkg_plots<-left_join(trwi_bkg_col,trwi_bkg_plots)
trwi_bkg_plots<-left_join(trwi_bkg_plots,bl_key)

trwi_all<-full_join(all_trwi,trwi_bkg_plots)

trwi_fil<-trwi_all%>%filter(ph_n_indiv!=0)%>%
  mutate(background=ifelse(background=="TRWI_mono","Mono",background))%>%
  filter(background!="Mono" & background!="BRHO" & background !="Hirtum")

trwi_mono_dat<-trwi_all%>%filter(ph_n_indiv!=0)%>%
  filter(background=="Mono")
##### need to add in 0 for mono bkg indiv (weeds?) and estimated bkg densities for L and H density bkg collections
######## model data including 0 indivs and all mono indivs
mod_dat_trwi<-trwi_fil%>%mutate(phyto2=phyto,
                                background2=background,
                                background3=background,
                                bg_n_indiv2=bg_n_indiv,
                                treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("ERBO":"TRWI")),~replace_na(.,0))%>%
  mutate_at(vars(c("ERBO":"TRWI")),as.numeric)%>%
  #mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>%
  mutate(trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("biomass","avg_biomass","water","grazing","nitrogen","phyto_unique","phyto_date","initials_ph","ph_notes","bg_date","initials_bg","bg_notes","scale_id","disturbance_cover","bare_cover","weed_notes","processing_notes","date"))
mod_dat_trwi<-left_join(mod_dat_trwi,mean_germ)
mod_dat_trwi<-left_join(mod_dat_trwi,seed_dat) 
########### weed cover for pler bkgs doesn't exist, can estimate from collections like bkgs?

#########  mono data
trwi_mono_dat2<-trwi_mono_dat%>%#mutate(phyto2=phyto,
  #background2=background,
  #background3=background,
  #bg_n_indiv2=bg_n_indiv,
  #treatment2=treatment)%>%
  # pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  # unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  #mutate_at(vars("Mono"),~replace_na(.,0))%>%
  #mutate_at(vars("Mono"),as.numeric)%>%
  #mutate_at(vars("wd_forb_indiv"),~replace_na(.,0))%>%
  mutate(trt=ifelse(treatment=="A",1,treatment),
         trt=ifelse(treatment=="AC",2,trt),
         trt=ifelse(treatment=="ACN",3,trt),
         trt=ifelse(treatment=="AG",4,trt),
         trt=ifelse(treatment=="AN",5,trt),
         trt=ifelse(treatment=="ANG",6,trt),
         trt=ifelse(treatment=="D",7,trt),
         trt=ifelse(treatment=="DG",8,trt),
         trt=ifelse(treatment=="DN",9,trt),
         trt=ifelse(treatment=="DNG",10,trt),
         trt=ifelse(treatment=="I",11,trt),
         trt=ifelse(treatment=="IG",12,trt),
         trt=ifelse(treatment=="IN",13,trt),
         trt=ifelse(treatment=="ING",14,trt))%>%
  mutate_at(vars("total_weed_cover"),~replace_na(.,0))%>% ### setting mono weeds as 0 for now until I think of solution to percent cover problem
  # filter(ph_n_indiv!=0)%>%
  ##### add in bkg indivs for monos, should be set to 0?
  select(-c("biomass","avg_biomass","water","grazing","nitrogen","phyto_unique","phyto_date","initials_ph","ph_notes","bg_date","initials_bg","bg_notes","scale_id","disturbance_cover","bare_cover","weed_notes","processing_notes","date"))
trwi_mono_dat2<-left_join(trwi_mono_dat2,mean_germ)
trwi_mono_dat2<-left_join(trwi_mono_dat2,seed_dat)


#####################

########merging species data
df_list<-list(mod_dat_acam,mod_dat_avba,mod_dat_erbo,mod_dat_gitr,mod_dat_lomu,mod_dat_pler,mod_dat_taca,mod_dat_trwi)
model_dat<-df_list %>% reduce(full_join)%>%
  select(-c("initials","avg_biomass","...10","n_flowers_immature","n_flowers_mature","total_flowers","avg_flowers"))
seed_na<-filter(model_dat,is.na(num_seeds)) ### check, 16 missing samples that are accounted for
model_dat<-model_dat%>%filter(!(is.na(num_seeds))) ## filter them out
model.dat<-model.dat%>%mutate_at(vars(c("neigh_weeds","total_weeds_n","total_weed_cover")),as.numeric)%>%
  mutate(weeds=((coalesce(neigh_weeds,total_weeds_n,total_weed_cover)))) ## make a single weed column, since each species is a spearate model should be OK

#write_csv(seed_na,"missings_seeds_Aug0125.csv")
#write_csv(model_dat,"career_model_data_Aug0125.csv")
#############################
# mono data
df_list2<-list(acam_mono_dat2,avba_mono_dat2,erbo_mono_dat2,gitr_mono_dat2,lomu_mono_dat2,pler_mono_dat2,taca_mono_dat2,trwi_mono_dat2)
model_mono_dat<-df_list2 %>% reduce(full_join)%>%
  select(-c("grass_cover","forb_cover","initials","n_flowers_immature","n_flowers_mature","total_flowers","avg_flowers","bg_n_indiv"))
#mono_seed_na<-filter(model_mono_dat,is.na(num_seeds)) ### check, 16 missing samples that are accounted for


########### using posteriors from mnegacomp as priors for lambda estimates
mc_posts<-read_csv("megacomp_posteriors_20240714_models.csv")%>% rowwise()%>%
  mutate(mean_lambda=mean(c(lambda_d,lambda_c)))

mc_posts2<-mc_posts%>%group_by(phyto)%>%
           mutate(sp_mean_lambda_mc=mean(mean_lambda),
                  sp_sd_lambda_mc=sd(mean_lambda))%>%
        select(c(phyto,sp_mean_lambda_mc,sp_sd_lambda_mc))%>%
  filter(phyto %in% c("ACAM","GITR","LOMU","PLER","TACA","TWIL"))%>%
  unique()
  
#write_csv(model_mono_dat,"career_model_mono_data_Aug0425.csv")
#### version with mono data together with all other data
all_dat<-full_join(model.dat,mono.dat)%>%  
  mutate_at(vars(c("TRWI":"ACAM")),~replace_na(.,0))%>%
  mutate_at(vars(c("weeds","total_weed_cover")),as.numeric)%>%
  mutate(weeds=((coalesce(weeds,total_weed_cover))))

all_dat2<-left_join(all_dat,mc_posts2)
all_dat3<-select(all_dat2,c(phyto,block,plot,sub,background,density,ph_n_indiv,num_seeds,avg_seeds,treatment,bg_n_indiv,TRWI:ACAM,weeds,mean_germ,mean_surv,sp_mean_lambda_mc,sp_sd_lambda_mc))
#write_csv(all_dat3,"career_model_data_Aug20_25.csv")


############## Creating datasheet of missing samples
taca_miss<-read_csv("taca_phyto_biomass.csv")
taca_miss$ph_n_indiv<-as.numeric(taca_miss$ph_n_indiv)
taca_bkg_miss<-read_csv("taca_bkg_collections.csv")


acam_miss<-read_csv("acam_phyto_biomass.csv")
acam_miss$ph_n_indiv<-as.numeric(acam_miss$ph_n_indiv)
acam_bkg_miss<-read_csv("acam_bkg_collections.csv")

pler_miss<-read_csv("pler_phyto_biomass.csv")%>%
  select(-"other_notes")%>%
  mutate_at(vars(c("ph_n_indiv","bg_n_indiv")),as.numeric)

pler_bkg_miss<-read_csv("pler_bkg_mass.csv")%>%
  select(-"other_notes")

gitr_miss<-read_csv("gitr_phyto_biomass.csv")%>%
  select(-"other_notes")%>%
  mutate_at(vars(c("ph_n_indiv","bg_n_indiv")),as.numeric)

gitr_bkg_miss<-read_csv("gitr_bkg_collections.csv")%>%
  select(-"other_notes")

lomu_miss<-read_csv("lomu_phyto_biomass.csv")%>%
  mutate_at(vars(c("ph_n_indiv","bg_n_indiv")),as.numeric)
lomu_bkg_miss<-read_csv("lomu_bkg_collections.csv")

######## Cases where phyto is NA, most likely not entered correctly, check before data cleaning
tph<-taca_miss%>%filter(is.na(ph_n_indiv)) # 8
pph<-pler_miss%>%filter(is.na(ph_n_indiv)) # 31
aph<-acam_miss%>%filter(is.na(ph_n_indiv)) # 0
lph<-lomu_miss%>%filter(is.na(ph_n_indiv)) # 7
gph<-gitr_miss%>%filter(is.na(ph_n_indiv)) # 79


######## Cases where phyto biomass is missing, most likely sample is missing
tpm<-taca_miss%>%filter(is.na(biomass)) # 95
ppm<-pler_miss%>%filter(is.na(biomass)) # 119
apm<-acam_miss%>%filter(is.na(biomass)) # 584
lpm<-lomu_miss%>%filter(is.na(biomass)) # 137
gpm<-gitr_miss%>%filter(is.na(biomass)) # 450
#### includes where phyto is 0

####### Cases where bg_indiv is missing, most likely not entered correctly
tbg<-taca_bkg_miss%>%filter(is.na(bg_n_indiv_collected)) # 0
pbg<-pler_bkg_miss%>%filter(is.na(bg_n_indiv_collected)) # 0
abg<-acam_bkg_miss%>%filter(is.na(bg_n_indiv_collected)) # 0
lbg<-lomu_bkg_miss%>%filter(is.na(bg_n_indiv_collected)) # 0
gbg<-gitr_bkg_miss%>%filter(is.na(bg_n_indiv_collected)) # 1

######## Cases where bkg biomass is missing, most likely sample is missing
tbm<-taca_bkg_miss%>%filter(is.na(biomass)) # 1
pbm<-pler_bkg_miss%>%filter(is.na(biomass)) # 1
abm<-acam_bkg_miss%>%filter(is.na(biomass)) # 32
lbm<-lomu_bkg_miss%>%filter(is.na(biomass)) # 1
gbm<-gitr_bkg_miss%>%filter(is.na(biomass)) # 4

miss_list<-list(tph,tpm,tbg,tbm,aph,apm,abg,abm,pph,ppm,pbg,pbm,lph,lpm,lbg,gbm,gph,gpm,gbg,gbm)
miss_list<-lapply(miss_list, function(x){
  mutate_at(x,vars("processing_notes"),as.character)
  })

miss_list<-lapply(miss_list, function(x){
  mutate_at(x,vars("avg_biomass"),as.numeric)
})


miss_dat<-miss_list %>% reduce(full_join)%>%
filter(background!="BRHO" & background!="Hirtum")

miss_dat2<-filter(miss_dat,ph_n_indiv!=0 | is.na(ph_n_indiv))%>%
  filter()
write_csv(miss_dat2,"missing_data_feb17_25.csv")

######### Checking missing data
acam_all%>%filter(is.na(biomass)) # 3
taca_all%>%filter(is.na(biomass)) # 1
pler_all%>%filter(is.na(biomass)) # 9
lomu_all%>%filter(is.na(biomass)) # 9
gitr_all%>%filter(is.na(biomass)) # 5

acam_all%>%filter(is.na(bg_n_indiv)) # 26
taca_all%>%filter(is.na(bg_n_indiv)) # 2
pler_all%>%filter(is.na(bg_n_indiv)) # 0
lomu_all%>%filter(is.na(bg_n_indiv)) # 6
gitr_all%>%filter(is.na(bg_n_indiv)) # 4

acam_all%>%filter(is.na(ph_n_indiv)) # 0
taca_all%>%filter(is.na(ph_n_indiv)) # 0
pler_all%>%filter(is.na(ph_n_indiv)) # 0
lomu_all%>%filter(is.na(ph_n_indiv)) # 0
gitr_all%>%filter(is.na(ph_n_indiv)) # 0
