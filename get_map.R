

##  Unraveling genotype-by-environment interaction with probability theory                                             
##  Dias et al. 2022.              
## 
##  Function to obtain the maximum a posteriori (MAP) value                                                      
##                                                                                                                     
##  Authors:    KOG Dias        <kaioolimpio@usp.br>                                                                   
##              JPR dos Santos  <jhowpd@gmail.com>                                                                     
##              MD Krause       <krause.d.matheus@gmail.com>                                                           
                                                                                                                   

get_map <- function(posterior) {
  posterior <- as.matrix(posterior) 
  if (ncol(posterior) > 1) {
    den = apply(posterior, 2, density)
    map = unlist(lapply(den, function(den)
      den$x[which.max(den$y)]))
  }
  else {
    den = density(posterior)
    map = den$x[which.max(den$y)]
  }
  return(map)
}
