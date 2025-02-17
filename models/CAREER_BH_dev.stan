// Try to write a BH model
// should be able to handle data with two separate environmental treatments
//input data
data{
  int <lower=1> N;//number of observations
  int Fecundity [N];// at time t+1
  int <lower=1> N_blocks;// groups
  int Blocks[N];// block column
  
  vector[N] N_i;// pop size at time t
  vector[N] trt;// treatment
  vector[N] g_i;// mean germ rate
  vector[N] s_i;//mean seed surv
  
// population sizes of interacting species at time t
vector [N] acam;
vector [N] avba;
vector [N] gitr;
vector [N] lomu;
vector [N] pler;
vector [N] taca;
vector [N] weeds;

}

//model parameters 
parameters{
  //random effects
  real <lower=0>sigma;
  real epsilon [N_blocks];
  real<lower=0> disp_dev; //might only be required with the negative binomial?
  
  //lambda
  real <lower=0,upper=15000> lambda_base;//change to megacomp priors if ever using this model form again
  real lambda_dev;
  
  //alphas
  real  alpha_weeds_base;
  real  alpha_acam_base;
  real  alpha_avba_base;
  real  alpha_gitr_base;
  real  alpha_lomu_base;
  real  alpha_pler_base;
  real  alpha_taca_base;
  
  //alpha dev from treatments
  real alpha_acam_dev;
  real alpha_avba_dev;
  real alpha_gitr_dev;
  real alpha_lomu_dev;
  real alpha_pler_dev;
  real alpha_weeds_dev;
  real alpha_taca_dev;
  
}

// model block
model{
// prediction vectors
  vector [N] F_hat;
  vector [N] F_hat2;
  
//priors
  sigma~gamma(1,1);
  epsilon~gamma (sigma,sigma);
  lambda_base~exponential(0.0009);//check on this for my data, uninformative prior?
  lambda_dev~normal(0,100);//check on this
  
  alpha_acam_base~normal(0,5);
  alpha_acam_dev~normal(0,5);
  alpha_avba_base~normal(0,5);
  alpha_avba_dev~normal(0,5);
  alpha_gitr_base~normal(0,5);
  alpha_gitr_dev~normal(0,5);
  alpha_lomu_base~normal(0,5);
  alpha_lomu_dev~normal(0,5);
  alpha_pler_base~normal(0,5);
  alpha_pler_dev~normal(0,5);
  alpha_taca_base~normal(0,5);
  alpha_taca_dev~normal(0,5);
  alpha_weeds_base~normal(0,5);
  alpha_weeds_dev~normal(0,5);
  
//BH model

for(i in 1:N){

F_hat[i]=s_i[i]*(1-g_i[i])+(N_i[i]*g_i[i]*(lambda_base + lambda_dev*trt[i])/(1+
acam[i]*(alpha_acam_base*g_i[i] + alpha_acam_dev*trt[i])+
avba[i]*(alpha_avba_base*g_i[i] + alpha_avba_dev*trt[i])+
gitr[i]*(alpha_gitr_base*g_i[i] + alpha_gitr_dev*trt[i])+
lomu[i]*(alpha_lomu_base*g_i[i] + alpha_lomu_dev*trt[i])+
pler[i]*(alpha_pler_base*g_i[i] + alpha_pler_dev*trt[i])+
taca[i]*(alpha_taca_base*g_i[i] + alpha_taca_dev*trt[i])+
weeds[i]*(alpha_weeds_base*g_i[i] + alpha_weeds_dev*trt[i])));

F_hat2[i]=F_hat[i]*epsilon[Blocks[i]];//block effect
  
}

//likelihood
Fecundity~neg_binomial_2(F_hat2,disp_dev);//check on what distribution to use here

}
