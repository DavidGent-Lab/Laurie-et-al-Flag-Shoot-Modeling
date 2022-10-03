// Mildew data logistic regression
// December 2017


data {
   // Define variables in data
   // Number of observations (an integer)
   int<lower=0> N;
   // Number of parameters
   int<lower=0> p;
   // Variables
   int<lower=0,upper=1> FS_Binary[N];
   int<lower=0,upper=1> Second_Year[N];
 }
 
parameters {
 // Define parameters to estimate
 real beta[p];
}
 
transformed parameters  {
 // Probability trasformation from linear predictor
 real<lower=0> odds[N];
 real<lower=0, upper=1> prob[N];

 for (i in 1:N) {
   odds[i] = exp(beta[1] + beta[2]*Second_Year[i]);
   prob[i] = odds[i] / (odds[i] + 1);
 }
}

model {
  // Prior part of Bayesian inference (flat if unspecified)

  // Likelihood part of Bayesian inference
 FS_Binary ~ bernoulli(prob);
}

