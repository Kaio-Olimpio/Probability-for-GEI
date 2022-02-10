
##  Leveraging probability concepts for cultivar recommendation in multi‑environment trials                                             
##  Dias et al. 2022.               
## 
##  Function to perform the posterior predictive check                                                     
##                                                                                                                     
##  Authors:    KOG Dias        <kaio.o.dias@ufv.br>                                                                   
##              JPR dos Santos  <jhowpd@gmail.com>                                                                     
##              MD Krause       <krause.d.matheus@gmail.com>   

pp_check <- function(fit, data_list) {
  out = rstan::extract(fit, permuted = TRUE) # Extracting posterior values from the stan fit object
  ns = length(out$mu) # Counting the number of simulations
  y = data_list$y     # Indexing response variable
  N = length(y)       # number of observations of the response variable
  temp = apply(out$y_gen, 1, function(x) {
    c(
      max = max(x) > max(y),
      min = min(x) > min(y),
      median = quantile(x, 0.5) > quantile(y, 0.5),
      mean = mean(x) > mean(y),
      sd = sd(x) > sd(y)
    )
  })
  p_value_max = sum(temp["max",]) / ns
  p_value_min = sum(temp["min",]) / ns
  p_value_median = sum(temp["median.50%",]) / ns
  p_value_mean = sum(temp["mean",]) / ns
  p_value_sd = sum(temp["sd",]) / ns
  # WAIC information criterion
  temp_v = apply(out$y_log_like, 2, function(x) {
    c(val = log((1 / ns) * sum(exp(x))),
      var = var(x))
  })
  lppd = sum(temp_v["val", ])
  p_WAIC2 = sum(temp_v["var", ]) # Effective number of parameters
  elppd_WAIC2 = lppd - p_WAIC2
  WAIC2 = -2 * elppd_WAIC2
  output_p_check = round(t(
    cbind(
      p_value_max = p_value_max,
      p_value_min = p_value_min,
      p_value_median = p_value_median,
      p_value_mean = p_value_mean,
      p_value_sd = p_value_sd,
      Eff_N_of_parameters = p_WAIC2,
      WAIC2 = WAIC2
    )
  ), 4)
  colnames(output_p_check) <-
    "Posterior predictive check statistics"
  return(output_p_check)
}