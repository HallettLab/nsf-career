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
bl_key<-read_csv("block_treatments.csv")

############# clean up taca
taca<-taca %>%left_join(bl_key,by="block")%>%
  #filter(biomass!="NA")%>%
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
taca_all<-full_join(all_taca,taca_bkg_plots)%>%
  filter(ph_n_indiv!=0)%>%
  filter(background!="Mono")
  
#################
mod_dat_taca<-taca_all%>%mutate(phyto2=phyto,
                            background2=background,
                            background3=background,
                            bg_n_indiv2=bg_n_indiv,
                            treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("Wildenovii":"TACA")),~replace_na(.,0))%>%
  mutate_at(vars(c("Wildenovii":"TACA")),as.numeric)%>%
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
  select(-c("biomass","avg_biomass","water","grazing","nitrogen","ph_unique","ph_date","ph_initials","ph_notes","bg_date","bg_initials","bg_notes","seed_notes","unique.ID","wd_grass_indiv","wd_forb_indiv","wd_grass_cover","wd_forb_cover","bare_cover","cover_notes"))
mod_dat_taca<-left_join(mod_dat_taca,mean_germ)
mod_dat_taca<-left_join(mod_dat_taca,seed_dat) #### gets background and phyto collections in, no monos, are those only used for setting lambda priors?

########## taca bkgs don't have weed counts so is set to 0, can estimate like bkg plots?

################## Pler data clean
pler<-read_csv("pler_phyto_biomass.csv") #has bg indivs in it and weed cover data

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

pler_all<-full_join(pler,pler_bkg_plots)%>%
  filter(ph_n_indiv!=0)%>%
  mutate(background=ifelse(background=="PLER_mono","Mono",background))%>%
  filter(background!="Mono" & background!="BRHO" & background !="Hirtum")


##### need to add in 0 for mono bkg indiv (weeds?) and estimated bkg densities for L and H density bkg collections
######## model data including 0 indivs and all mono indivs
mod_dat_pler<-pler_all%>%mutate(phyto2=phyto,
                            background2=background,
                            background3=background,
                            bg_n_indiv2=bg_n_indiv,
                            treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("Wildenovii":"PLER")),~replace_na(.,0))%>%
  mutate_at(vars(c("Wildenovii":"PLER")),as.numeric)%>%
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


acam_all<-full_join(all_acam,acam_bkg_plots)%>%
  filter(ph_n_indiv!=0)%>%
  filter(background!="Mono")


#################
mod_dat_acam<-acam_all%>%mutate(phyto2=phyto,
                                background2=background,
                                background3=background,
                                bg_n_indiv2=bg_n_indiv,
                                treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("Wildenovii":"ACAM")),~replace_na(.,0))%>%
  mutate_at(vars(c("Wildenovii":"ACAM")),as.numeric)%>%
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

lomu_all<-full_join(all_lomu,lomu_bkg_plots)%>%
  filter(ph_n_indiv!=0)%>%
  filter(background!="Mono")


#################
mod_dat_lomu<-lomu_all%>%mutate(phyto2=phyto,
                                background2=background,
                                background3=background,
                                bg_n_indiv2=bg_n_indiv,
                                treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("Wildenovii":"LOMU")),~replace_na(.,0))%>%
  mutate_at(vars(c("Wildenovii":"LOMU")),as.numeric)%>%
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

########## GITR
##################
gitr<-read_csv("gitr_phyto_biomass.csv") #has bg indivs in it and weed cover data

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
  filter(bg_n_indiv_collected !=0)%>%
  mutate_at(vars("avg_biomass"),as.numeric)%>%
  mutate(num_seeds=biomass *55.5053883 ,
         avg_seeds=avg_biomass *55.5053883 )%>%
  rename(phyto=species,ph_n_indiv=bg_n_indiv_collected)%>%
  select(-c("bg_collection_notes","processing_notes","scale_id","date","bg_date_collected","bg_initials","other_notes"))


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

gitr_all<-full_join(gitr,gitr_bkg_plots)%>%
  filter(ph_n_indiv!=0)%>%
  mutate(background=ifelse(background=="GITR_mono","Mono",background))%>%
  filter(background!="Mono" & background!="BRHO" & background !="Hirtum")


##### need to add in 0 for mono bkg indiv (weeds?) and estimated bkg densities for L and H density bkg collections
######## model data including 0 indivs and all mono indivs
mod_dat_gitr<-gitr_all%>%mutate(phyto2=phyto,
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

########### ERBO
erbo<-erbo%>%left_join(bl_key,by="block")%>%
#filter(biomass!="NA")%>%
mutate(avg_flowers=total_flowers/ph_n_indiv)%>%
  mutate(num_seeds=total_flowers*5)%>%
  mutate(avg_seeds=avg_flowers*5)

########### getting estimates of bkg plots for bkg collections
bkg_plots<-read_csv("bkg_estimates_all.csv")
na_plots<-bkg_plots%>%filter(background!="Mono")%>%
  filter(is.na(bg_n_indiv))%>%
  filter(is.na(ph_date) & is.na(ph_initials))

erbo_bkg_col<-read_csv("erbo_bkg_counts.csv")%>%
  filter(bg_n_indiv_collected !=0)%>%
  mutate(num_seeds=total_flowers *5 ,
         avg_flowers=total_flowers/bg_n_indiv_collected,
         avg_seeds=avg_flowers *5 )%>%
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

erbo_all<-full_join(erbo,erbo_bkg_plots)%>%
  filter(ph_n_indiv!=0)%>%
  mutate(background=ifelse(background=="ERBO_mono","Mono",background))%>%
  filter(background!="Mono" & background!="BRHO" & background !="Hirtum")


##### need to add in 0 for mono bkg indiv (weeds?) and estimated bkg densities for L and H density bkg collections
######## model data including 0 indivs and all mono indivs
mod_dat_erbo<-erbo_all%>%mutate(phyto2=phyto,
                                background2=background,
                                background3=background,
                                bg_n_indiv2=bg_n_indiv,
                                treatment2=treatment)%>%
  pivot_wider(names_from = background2,values_from = bg_n_indiv2)%>%
  unite("combo",c("phyto2","background3","treatment2"),sep = "_")%>%
  mutate_at(vars(c("GITR":"ERBO")),~replace_na(.,0))%>%
  mutate_at(vars(c("GITR":"ERBO")),as.numeric)%>%
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
  select(-c("biomass","water","grazing","nitrogen","phyto_unique","phyto_date","Initials_ph","Phyto_notes","background_date","Initials_bg","Background_notes","scale_id","disturbance_cover","bare_cover","weed_notes","processing_notes","initials"))
mod_dat_erbo<-left_join(mod_dat_erbo,mean_germ)
mod_dat_erbo<-left_join(mod_dat_erbo,seed_dat) 

########merging species data
df_list<-list(mod_dat_acam,mod_dat_erbo,mod_dat_gitr,mod_dat_lomu,mod_dat_pler,mod_dat_taca)
model_dat<-df_list %>% reduce(full_join)
#write_csv(model_dat,"career_model_data.csv")
#############################


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
