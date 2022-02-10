########################################################################################################################
##                                                                                                                    ##                                        
##                                                                                                                    ##
##  Bayesian Analysis for MET data                                                                                    ##
##  Date: Feb, 2022                                                                                                   ##
##                                                                                                                    ##
##  Leveraging probability concepts for cultivar recommendation in multi‑environment trials                           ##
##  Dias et al. 2022. Model 3 (M3)                                                                                    ##
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
df <- read.csv("maize_dataset.csv", h = TRUE)

# Defining Factors 
df$M <- df$Region %>% as.factor        # Mega-Region
df$B <- df$Block  %>% as.factor        # Block
df$R <- df$Rep %>% as.factor           # Replicates
df$H <- df$Hybrid %>% as.factor        # Genotypes
df$L <- df$Location %>% as.factor      # Locations

#----------------------------------------------Data for stan-------------------------------------------------#

# Number of observations
n <- df %>% nrow

# Designs Matrices
Z_1 <- model.matrix( ~ -1 + R:L, data = df)  # Replication effects
Z_2 <- model.matrix( ~ -1 + B:L, data = df)  # Block effects
Z_3 <- model.matrix( ~ -1 + H, data = df)    # Genotype effects
Z_4 <- model.matrix( ~ -1 + L, data = df)    # Location effects
Z_5 <- model.matrix( ~ -1 + H:L, data = df)  # GxL effects
Z_6 <- model.matrix( ~ -1 + M, data = df)    # Region effects
Z_7 <- model.matrix( ~ -1 + H:M, data = df)  # GxM effects

# Number of columns for each design matrix
p_1 <- ncol(Z_1)
p_2 <- ncol(Z_2)
p_3 <- ncol(Z_3)
p_4 <- ncol(Z_4)
p_5 <- ncol(Z_5)
p_6 <- ncol(Z_6)
p_7 <- ncol(Z_7)

# Subset response variable
y <- df$GY

# Create the known global hyperparameter:
phi <- max(y) * 10

# Create a list to store data for stan
df_stan <- list(n = n,
               p_1 = p_1,
               p_2 = p_2,
               p_3 = p_3,
               p_4 = p_4,
               p_5 = p_5,
               p_6 = p_6,
               p_7 = p_7,
               Z_1 = Z_1,
               Z_2 = Z_2,
               Z_3 = Z_3,
               Z_4 = Z_4,
               Z_5 = Z_5,
               Z_6 = Z_6,
               Z_7 = Z_7,
               y = y,
               phi = phi) 

#---------------------------------------------Stan model code------------------------------------------------#

stan_df <- "

  data{
    // Number of observations
    int<lower=1> n;
  
    // Number of parameters
    int<lower=1> p_1;
    int<lower=1> p_2;
    int<lower=1> p_3;
    int<lower=1> p_4;
    int<lower=1> p_5;
    int<lower=1> p_6;
    int<lower=1> p_7;
    
    // Designs matrices
    matrix[n, p_1] Z_1;
    matrix[n, p_2] Z_2;
    matrix[n, p_3] Z_3;
    matrix[n, p_4] Z_4;
    matrix[n, p_5] Z_5;
    matrix[n, p_6] Z_6;
    matrix[n, p_7] Z_7;
  
    // Phenotype vector
    real y[n];
    
    // Global hyperparameter
    real phi;
  }

    parameters{
    // Residual standard deviation parameter/hyperparameters
    real<lower=0> s_sigma;
    real<lower=0> sigma;
    
    // Mean parameter/hyperparameters
    real<lower=0> s_mu;
    real mu;

    // Replication parameter/hyperparameters
    real<lower=0> s_r;
    vector[p_1] r;
    
    // Block parameter/hyperparameters
    real<lower=0> s_b;
    vector[p_2] b;
    
    // Genotype parameter/hyperparameters
    real<lower=0> s_g;
    vector[p_3] g; 
    
    // Location parameter/hyperparameters
    real<lower=0> s_l;
    vector[p_4] l; 
    
    // Hybrid by Location parameter/hyperparameters
    real<lower=0> s_gl;
    vector[p_5] gl; 
    
    // Region parameter/hyperparameters
    real<lower=0> s_m;
    vector[p_6] m; 
    
    // Hybrid by Region parameter/hyperparameters
    real<lower=0> s_gm;
    vector[p_7] gm; 
    
    // Defining variable to generate data from the model
    real y_gen[n];
  }
  
  transformed parameters{
  
    // Declaring variables to receive input
    vector[n] expectation;
    
    // Computing the expectation of the likelihood function
    expectation = mu + Z_1*r + Z_2*b + Z_3*g + Z_4*l + Z_5*gl + Z_6*m + Z_7*gm;
  } 

  model{
  
    // Conditional prior probabilities distributions for residual standard deviation
    s_sigma ~ cauchy(0, phi);
    sigma ~ cauchy(0, s_sigma); 
    
    // Conditional prior probabilities distributions for the mean
    s_mu ~ cauchy(0, phi);
    mu ~ normal(0, s_mu);
    
    // Conditional prior probabilities distributions for replications
    s_r ~ cauchy(0, phi);
    r ~ normal(0, s_r);
    
    // Conditional prior probabilities distributions for blocks
    s_b ~ cauchy(0, phi);
    b ~ normal(0, s_b);   
    
    // Conditional prior probabilities distributions  for genotypes
    s_g ~ cauchy(0, phi);
    g ~ normal(0, s_g);
    
    // Conditional prior probabilities distributions  for locations
    s_l ~ cauchy(0, phi);
    l ~ normal(0, s_l);

    // Conditional prior probabilities distributions  for genotype by location
    s_gl ~ cauchy(0, phi);
    gl ~ normal(0, s_gl);
    
    // Conditional prior probabilities distributions  for Regions
    s_m ~ cauchy(0, phi);
    m ~ normal(0, s_m);
    
    // Conditional prior probabilities distributions  for genotype by Region
    s_gm ~ cauchy(0, phi);
    gm ~ normal(0, s_gm);
    
    // Specifying the likelihood
    y ~ normal(expectation, sigma);

    // Generating data from the model
    y_gen ~ normal(expectation, sigma);
  }
    
    generated quantities {
    real y_log_like[n];
      for (j in 1:n) {
      // Computing log-likelihood of the observed data
        y_log_like[j] = cauchy_lpdf(y[j] | expectation[j], sigma);
      }
      
} "

#----------------------------------------------Run stan code-------------------------------------------------#

# Compile into C++ the stan models
stan_df_comp <- stan_model(model_code = stan_df)

# Fit the stan models
Model3 <- sampling(stan_df_comp,
                     data = df_stan, 
                     iter = 4000,
                     cores = 4,
                     chain = 4)

# Extract stan results
out <- rstan::extract(Model3, permuted = TRUE)

# Subset replication posterior effects
r_post <- out$r

# Subset blocks posterior effects
b_post <- out$b

# Subset hybrids posterior effects
g_post <- out$g
dim(g_post)

# Subset genotype by location posterior effects
gl_post <- out$gl 
dim(gl_post)

# Subset data generated by the model
y_gen_post <- out$y_gen

# Subset the replicatin variance
s2_r_post <- (out$s_r)^2
mean(s2_r_post) # replication variance

# Subset the block variance
s2_b_post <- (out$s_b)^2
mean(s2_b_post) # Block variance

# Subset the genotype (genetic) variance
s2_g_post <- (out$s_g)^2
mean(s2_g_post) # Genetic variance

# Subset the genotype by location variance
s2_gl_post <- (out$s_gl)^2
mean(s2_gl_post) # Genetic by location variance

# Subset the location variance
s2_l_post <- (out$s_l)^2
mean(s2_l_post) # location variance

# Subset the Region variance
s2_m_post <- (out$s_m)^2
mean(s2_m_post) # Region variance

# Subset the genotype by Region variance
s2_gm_post <- (out$s_gm)^2
mean(s2_gm_post) # Genetic by Region variance

# Subset the error variance
sigma <- (out$sigma)^2
mean(sigma)

# Getting maximum a posteriori values (MAP)
source('get_map.R') # Function to obtain the maximum a posteriori (MAP) value

# For the replication effects
r_map <- get_map(r_post)

# For the block effects
b_map <- get_map(b_post)

# For the genotype effects
g_map <- get_map(g_post)

# For the genotype effects
gl_map <- get_map(gl_post)


# Save Image
save.image("Model3.RData")







