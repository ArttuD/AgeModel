
functions {
   vector polyFit(vector x, real alpha, real beta) {
        return alpha+beta*x;
   }
   vector AgeFit(vector x, real a, real b, real alpha, real beta){
      return a/((b-1)*beta)*(pow(alpha,(1-b)) - pow((beta*x + alpha),(1-b))); //+ e 
      //return a*(c*x+d) + b
   }
   vector coFit(vector Fe_p, vector Mn_p, vector Co_p){
      return 50*Co_p ./ (Fe_p + Mn_p); //
   }
}
data {
   int n_i; // Number of Co points
   int n_j; // Number of Co points

   vector[n_i] Co_;
   vector[n_i] Co_e;
   vector[n_i] Fe_;
   vector[n_i] Fe_e;
   vector[n_i] Mn_;
   vector[n_i] Mn_e;

   //vector[n_i] Co;
   //vector[n_i] Co_error;

   vector[n_i] Co_d; 
   vector[n_i] Co_d_error;
   
   vector[n_j] Ages; 
   vector[n_j] Ages_error;
   vector[n_j] Ages_d; 
   vector[n_j] Ages_d_error; 

   array[n_i] int n_Co_n_indices;
   array[n_i] int n_Co_indices;
   array[n_j] int n_j_indices;
}
parameters {
    real<lower = 0> alpha;
    real<upper = 0> beta;

    real a;
    real b;

   vector<lower = 0>[n_i] Co_mu;
   vector<lower = 0>[n_i] Fe_mu;
   vector<lower = 0>[n_i] Mn_mu;

   vector<lower=0>[n_i] Co_p_sigma;
   vector[n_i] Co_n_mu;

    real<lower = 0> Co_d_mu;
    real<lower = 0> Co_d_sigma;
    vector<lower = 0>[n_i] d_Co_mu;

    real<lower=0>Co_sigma_mu;
    real<lower=0>Co_sigma_sigma;
    vector<lower=0>[n_i]Co_error_;

    real<lower = 0> Ages_d_mu;
    real<lower = 0> Ages_d_sigma;
    vector<lower = 0>[n_j] d_Ages_mu;

    real<lower=0>Ages_sigma_mu;
    real<lower=0>Ages_sigma_sigma;
    vector<lower=0>[n_j]Ages_error_;

    real<lower=0> nu;
    real<lower=0> nu2;
    real<lower=0> nu3;

   real Co_mu_mu;
   real<lower = 0> Co_mu_sigma;
   real Fe_mu_mu;
   real<lower = 0> Fe_mu_sigma;
   real Mn_mu_mu;
   real<lower = 0> Mn_mu_sigma;

}
transformed parameters {
   vector[n_i] Co_p = 10*coFit(Fe_mu, Mn_mu, Co_mu);
   vector[n_i] Co_fit = polyFit(d_Co_mu, alpha, beta);
   vector[n_j] Ages_fit = AgeFit(d_Ages_mu, a, b, alpha, beta); //e
}
model {

   alpha ~ normal(0.35,0.1);
   beta ~ normal(-0.3,0.1);

   a ~ normal(0.1,0.1);
   b ~ normal(1,0.1);

   Co_mu_mu ~ normal(0, 0.5);
   Co_mu_sigma ~ inv_gamma(4,1);
   Co_mu ~ normal(Co_mu_mu,Co_mu_sigma+1e-4);

   Fe_mu_mu ~ normal(5, 1);
   Fe_mu_sigma ~ inv_gamma(4,1);
   Fe_mu ~ normal(Fe_mu_mu,Fe_mu_sigma+1e-4);

   Mn_mu_mu ~ normal(5, 1);
   Mn_mu_sigma ~ inv_gamma(4,1);
   Mn_mu ~ normal(Mn_mu_mu,Mn_mu_sigma+1e-4);

   Co_  ~ normal(Co_mu, Co_e);
   Mn_  ~ normal(Mn_mu, Mn_e);
   Fe_  ~ normal(Fe_mu, Fe_e);

   Co_p_sigma ~ inv_gamma(4,1);
   nu3 ~ gamma(2,0.1);

   Co_n_mu ~ student_t(nu3, Co_p, Co_p_sigma+1e-4);

   Co_d_mu ~ normal(0.5,0.3);
   Co_d_sigma ~ inv_gamma(4,1);
   d_Co_mu ~ normal(Co_d_mu, Co_d_sigma);

   Co_sigma_mu ~ normal(0.01,0.05); 
   Co_sigma_sigma ~ inv_gamma(4,1);
   Co_error_ ~ normal(Co_sigma_mu, Co_sigma_sigma+1e-4);

   nu2 ~ gamma(2,0.1);
   Co_d ~ normal(d_Co_mu, Co_d_error+1e-4);

   Co_n_mu ~ student_t(nu2, Co_fit, Co_error_+1e-4);

   Ages_d_mu ~ normal(0.5, 0.3);
   Ages_d_sigma ~ inv_gamma(4,1);
   d_Ages_mu ~ normal(Ages_d_mu, Ages_d_sigma+1e-4);
   Ages_d ~ normal(d_Ages_mu, Ages_d_error+1e-4);

   Ages_sigma_mu ~ normal(0.1,0.1); 
   Ages_sigma_sigma ~ inv_gamma(4,1);
   Ages_error_ ~ normal(Ages_sigma_mu, Ages_sigma_sigma+1e-4);
   Ages_error ~ normal(Ages_error_, 0.1);
   nu ~ gamma(2,0.1);
   Ages ~ student_t(nu, Ages_fit, Ages_error_+1e-4);
}
generated quantities {
vector[n_j+n_i+n_i+n_i] log_lik;
array[n_j+n_i+n_i] real y_hat;

for (i in 1:n_i){
     log_lik[i] = normal_lpdf(Co_[i] | Co_mu[i], Co_e[i]);
     }
for (i in 1:n_i){
    log_lik[i+n_i] = normal_lpdf(Fe_[i] | Fe_mu[i], Fe_e[i]);
}
for (i in 1:n_i){
     log_lik[i+n_i+n_i] = normal_lpdf(Mn_[i] | Mn_mu[i], Mn_e[i]);
 }
for(i in 1:n_j) {
    log_lik[i+n_i+n_i+n_i] = student_t_lpdf(Ages[i]| nu ,Ages_fit[i], Ages_error_[i]);
}

 //y_hat = student_t_rng(5,intAge(d_Co, a, b, c, alpha, beta);
 y_hat[n_Co_n_indices] = normal_rng(Co_p, Co_p_sigma);
 y_hat[n_Co_indices] = student_t_rng(nu2,Co_fit,Co_error_); 
 y_hat[n_j_indices] = student_t_rng(nu,Ages_fit,Ages_error_);
}
