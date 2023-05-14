# add noise to proportions
add_noise_to_pi <- function(pi_m, b, kappa, seed){
  set.seed(seed)
  n <- nrow(pi_m)
  pi_errors <- rnorm(n, b*mean(pi_m[,1]), kappa^2*0.2)
  # set sigma=0.2 such that
  # 0.25*pi_errors: sd=0.05
  # 0.5*pi_errors: sd=0.1
  # 1*pi_errors: sd=0.2
  pi_vec <- pi_m[,1]
  pi_vec <- pi_vec + pi_errors
  # bound the estimates
  pi_vec <- sapply(pi_vec, function(x){
    if(x > 1){
      return(1)
    }else if(x < 0){
      return(0)
    }else{
      return(x)
    }
  })
  # get cell type two estimates
  new_pi_m <- pi_m
  new_pi_m[,1] <- pi_vec
  new_pi_m[,2] <- 1-new_pi_m[,1]
  return(new_pi_m)
}
