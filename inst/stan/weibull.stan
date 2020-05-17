#include /chunks/logliks.stan

data{
  // setting observed data:
  int<lower=1> n;
  int<lower=1> p;
  matrix[n, p] x;
  vector<lower=0>[n] time;
  real<lower=0> tau;
  vector<lower=0, upper=1>[n] status;
  int<lower=0, upper=1> approach;
  int<lower=0, upper=1> type;

  // setting hyperparameters:
  real<lower=0> a_lambda;
  real<lower=0> b_lambda;
  real<lower=0> a_gamma;
  real<lower=0> b_gamma;
  real mu_beta;
  real<lower=0> sigma_beta;

}


parameters {
  vector[p] beta;
  real<lower=0> lambda;
  real<lower=0> gamma;
}

transformed parameters{
  real<lower=0> lambda_transf;

  lambda_transf = lambda*(tau^gamma);
}



model{

  vector[n] H0t;
  vector[n] ht;
  vector[n] h0t;
  vector[n] lambda_PO;

  for(i in 1:n){
      H0t[i] = lambda_transf*exp(gamma*log(time[i]))*exp(x[i,]*beta);
      ht[i] = lambda_transf*gamma*exp((gamma-1)*log(time[i]))*exp(x[i,]*beta);
      h0t[i] = lambda_transf*gamma*exp((gamma-1)*log(time[i]));
      lambda_PO[i] = exp(x[i,]*beta);
  }

  if(type ==0){
  vector[num_elements(status)] loglik;
  loglik = loglik_weibull(H0t, ht, status);
  target += sum(loglik);}

  if(type ==1){
  vector[num_elements(status)] loglik;
  loglik = loglik_weibull_PO(H0t, h0t, status, lambda_PO);
  target += sum(loglik);}


  if(approach==1){
  // setting the prior distributions:
  lambda ~ gamma(a_lambda, b_lambda);
  gamma ~ gamma(a_gamma, b_gamma);
  beta ~ normal(mu_beta, sigma_beta);
  }

}





