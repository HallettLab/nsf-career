// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html
data{
  int <lower=1> N;//number of observations
  int Fecundity [N];// at time t+1
  int<lower = 1> S; // Number of species
  matrix[N,S] SpMatrix; // Matrix of abundances for each species 
  int<lower = 0> Intra[S]; // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations
  int <lower=1> N_blocks;// groups
  int Blocks[N];// block column
  vector[N] g_i;// mean germ rate
  vector[N] s_i;//mean seed surv
  vector[N] N_i;//indivs
  real mc_mean_lambda;
  real mc_sd_lambda;
  
  // The below values define the regularized horseshoe priors used for species-specific parameters
real tau0; 		// determines the scale of the global shrinkage parameter (tau)
real slab_scale;	// scale for significant alpha_sp values
real slab_df;		// effective degrees of freedom for significant alpha_sp values
}

transformed data{
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5*slab_df;
}

parameters{
  vector [2] lambdas;
  real alpha_generic_tilde;
  real alpha_intra_tilde;
  vector [S] alpha_hat_ij_tilde;
  vector<lower = 0>[S] local_shrinkage_ij;
  real<lower = 0> c2_tilde;
  real<lower = 0> tau_tilde;
  real<lower=0> disp_dev; //might only be required with the negative binomial?

  //random effects
  real <lower=0> sigma;
  real epsilon [N_blocks];
  //real<lower=0> disp_dev; //might only be required with the negative binomial?
}
transformed parameters{
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the normalized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  real c2;
  real tau;
  vector[S] alpha_hat_ij;
  vector[S] local_shrinkage_ij_tilde;
  real alpha_generic; 
  real alpha_intra;
  
  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
    for(s in 1:S){
      local_shrinkage_ij_tilde[s] = sqrt( c2 * square(local_shrinkage_ij[s]) / (c2 + square(tau) * square(local_shrinkage_ij[s])) );
      alpha_hat_ij[s] = tau * local_shrinkage_ij_tilde[s] * alpha_hat_ij_tilde[s];
    }
  

  // scale the lambdas and alphas values
  alpha_generic = 3 * alpha_generic_tilde - 6;
  alpha_intra = 3 * alpha_intra_tilde - 6;
}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a vector of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] F_hat2;
  vector[N] interaction_effects;
  row_vector[S] alpha_ij;
  vector[N] lambda_i;

  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  //print("target = ", target(), " when ag_tilde = ", alpha_generic_tilde);

  alpha_intra_tilde ~ normal(0,1);
  //print("target = ", target(), " when ai_tilde = ", alpha_intra_tilde);

  lambdas ~ normal(mc_mean_lambda, mc_sd_lambda);//using megacomp posteriors, except for AVBA and ERBO
//print("target = ", target(), " when lambdas = ", lambdas);

  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

    alpha_hat_ij_tilde ~ normal(0,1);
 // print("target = ", target(), " when ahij_tilde = ", alpha_hat_ij_tilde);

    local_shrinkage_ij ~ cauchy(0,1);
   // print("target = ", target(), " when ls_ij = ", local_shrinkage_ij);

  tau_tilde ~ cauchy(0,1);
   // print("target = ", target(), " when tau_tilde = ", tau_tilde);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
     // print("target = ", target(), " when c2_tilde = ", c2_tilde);

 sigma~gamma(1,1);
     // print("target = ", target(), " when sigma = ", sigma);

 epsilon~gamma(sigma,sigma);
   // print("target = ", target(), " when epsilon_ij = ", epsilon);

  // implement the biological model

  for(i in 1:N){
    lambda_i[i] = exp(lambdas[1] + lambdas[2]);
    for(s in 1:S){
        alpha_ij[s] = exp((1-Intra[s]) * alpha_generic + Intra[s] * alpha_intra + (1-Intra[s]) * alpha_hat_ij[s]);
    }
    interaction_effects[i] = sum(alpha_ij.* SpMatrix[i,]);
    
    F_hat[i] = (s_i[i]*(1-g_i[i]))+(N_i[i]*g_i[i]*lambda_i[i]) / (1 + interaction_effects[i]); //(s_i[i]*(1-g_i[i]))* removing this 
    F_hat2[i]=F_hat[i]*epsilon[Blocks[i]];//block effect
}
//likelihood
//Fecundity~poisson(F_hat2);//check on what distribution to use here
Fecundity~neg_binomial_2(F_hat2,disp_dev);//check on what distribution to use here

}
