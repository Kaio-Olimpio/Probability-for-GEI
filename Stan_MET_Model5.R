########################################################################################################################
##                                                                                                                    ##                                        
##                                                                                                                    ##
##  Bayesian Analysis for MET data                                                                                    ##
##  Date: Feb, 2022                                                                                                   ##
##                                                                                                                    ##
##  Leveraging probability concepts for cultivar recommendation in multi‑environment trials                           ##
##  Dias et al. 2022. Model 5 (M5)                                                                                    ##
##                                                                                                                    ##
##  Authors:    KOG Dias        <kaio.o.dias@ufv.br>                                                                  ##
##              JPR dos Santos  <jhowpd@gmail.com>                                                                    ##
##              MD Krause       <krause.d.matheus@gmail.com>                                                          ##
##                                                                                                                    ##
########################################################################################################################

# Loading Library
library(dplyr)
library(rstan)

# Loading data
load("wheat_data.RData")

#----------------------------------------------Data for stan-------------------------------------------------#

# Number of observations
n <- dim(data)[1]
p_env <- nlevels(data$Env)

# SNPs matrix
Z <- wheat.X

# SNP matrix dimensions
p_z <- ncol(Z)
n_z <- nrow(Z)

# Subset response variable
y <- wheat.Y

# Create the known global hyperparameter:
phi <- max(y) * 10

# Create enviromental index
index_env <- as.integer(data$Env)

#  Create genotypic index
index_g <- as.integer(data$Line)

# Create a list to store data for stan
RRBLUP_stan <- list(n = n,
                    p_z = p_z,
                    n_z = n_z,
                    p_env = p_env,
                    Z = Z,
                    y = y,
                    index_env = index_env,
                    index_g = index_g,
                    phi = phi) 

#---------------------------------------------Stan model code------------------------------------------------#

stan_RRBLUP <- "

        data {
        //Number of row entries of the matrices or vectors:
        int<lower=1> n;
        int<lower=1> p_z;
        int<lower=1> n_z;
        int<lower=1> p_env;
     
        // Feature matrices:
        matrix[n_z, p_z] Z;

        // Phenotypic vector:
        matrix[n_z, p_env] y;

        // Global known hyperparameter:
        real phi;

        // Indexation vector for alpha
        int index_env[n];
        int index_g[n];
        }

        parameters {
        // Population mean parameter/hyperparameters:
        real<lower=0> s_mu;
        vector[p_env] mu;

        // Features parameter/hyperparameters:
        real<lower=0> s_alpha;
        matrix[p_z, p_env] alpha;  

        // Residual parameter/hyperparameters:
        real<lower=0> s_sigma;
        real<lower=0> sigma[p_env];

        // Defining variable to generate data from the model:
        real y_gen[n];
        }

        transformed parameters {
        // Declaring variables to receive input:
        matrix[n_z, p_env] expectation;

        // Computing the expectation of the likelihood function:
        for (i in 1:p_env){
        expectation[,i] = mu[i] + (Z * alpha[,i]);}
        }

        model {
        // Conditional probabilities distributions for mean:
        s_mu ~ cauchy(0, phi);
        mu ~ normal(0, s_mu);
        
        // Conditional probabilities distributions for features:
        s_alpha ~ cauchy(0, phi);
        to_vector(alpha) ~ normal(0, s_alpha);

        // Conditional probabilities distributions for residuals:
        s_sigma ~ cauchy(0, phi);
        sigma ~ cauchy(0, s_sigma);

        // Specifying the likelihood:
        for(i in 1:p_env){
                y[,i] ~ normal(expectation[,i], sigma[i]);
                }

        // Generating data from the model:
        to_vector(y_gen) ~ normal(to_vector(expectation), sigma[index_env]);

} " 


#----------------------------------------------Run stan code-------------------------------------------------#

# Compile into C++ the stan models
stan_RRBLUP_comp <- stan_model(model_code = stan_RRBLUP)

# Fit the stan models
ptm <- proc.time()
RRBLUP <- sampling(stan_RRBLUP_comp,
                   data = RRBLUP_stan,
                   iter = 2000,
                   cores = 8,
                   chain = 4)
timeRRblup_MET <- proc.time() - ptm

# Codes to extract stan results: please look the R codes from models M1, M2, M3, and M4.

# Save 
save.image("Model5.RData")

