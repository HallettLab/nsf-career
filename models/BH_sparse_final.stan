// final fits for BH model after determining parameters to be included post preliminary fits
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
  int Inclusion_ij[S];  // Boolean indicator variables to identify the species x reserve intercept parameters identified for inclusion in the final model
  real mc_mean_lambda;
  real mc_sd_lambda;
}

parameters{
  vector [2] lambdas;
  real alpha_generic_tilde;
  real alpha_intra_tilde;
  vector[S] alpha_hat_ij;
  
  //random effects
  real <lower=0> sigma;
  real epsilon [N_blocks];
  real<lower=0> disp_dev; //might only be required with the negative binomial?
}
transformed parameters{
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the normalized (and thus easier to sample)
  // 	counterparts declared in the parameters block

  real alpha_generic; 
  real alpha_intra;
  
  // scale the alpha values
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
  alpha_intra_tilde ~ normal(0,1);
  lambdas ~ normal(mc_mean_lambda, mc_sd_lambda);
  alpha_hat_ij ~ normal(0,1);
  sigma~gamma(1,1);
  epsilon~gamma(sigma,sigma);

  // implement the biological model
  
  for(i in 1:N){
    lambda_i[i] = exp(lambdas[1] + lambdas[2]);
    for(s in 1:S){
        alpha_ij[s] = exp((1-Intra[s]) * alpha_generic + Intra[s] * alpha_intra + Inclusion_ij[s] * alpha_hat_ij[s]);
    }   

    interaction_effects[i] = sum(alpha_ij.* SpMatrix[i,]);
    
  F_hat[i] = (s_i[i]*(1-g_i[i]))+(N_i[i]*g_i[i]*lambda_i[i]) / (1 + interaction_effects[i]);
  F_hat2[i]=F_hat[i]*epsilon[Blocks[i]];//block effect
}
//likelihood
Fecundity~neg_binomial_2(F_hat2,disp_dev);
}
