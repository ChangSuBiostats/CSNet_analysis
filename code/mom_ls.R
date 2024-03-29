#
# mom_ls.R
# 
# least squares (either ordinary least squares or iteratively re-weighted least squares)
# for moment-based regression
# that estimates cell-type-specific covariances
#

# 04/26: consider two modes

# least squares for moment-based regressions
mom_ls <- function(P, X, 
	methods = list(var='nnls', covar = 'wls'),
	weight_th = NULL){
	require(nnls)
  p <- ncol(X)
  K <- ncol(P)
  P_2 <- P^2
  # mean regression
  mu <- apply(X, 2, function(x) nnls::nnls(P, x)$x)
  M <- P %*% mu
  # variance regression
  Sigma_array <- array(NA, c(p, p, K))
  X_centered <- X - M
  Y <- X_centered^2
  if(methods$var == 'nnls'){
  	sigma_var <- apply(Y, 2, function(y) nnls(P_2, y)$x)
  }else if(methods$var == 'ols'){
  	sigma_var <- solve(t(P_2) %*% P_2) %*% (t(P_2) %*% Y)
  }
  for(k in 1:K){
  	diag(Sigma_array[, , k]) <- sigma_var[k, ]
  }
  # covariance regression
  # weights
  if(methods$covar == 'ols'){
  	n <- nrow(P)
  	w <- rep(1, n)
  }else if(methods$covar == 'wls'){
  	obs_w <- P_2 %*% sigma_var
  }
  for(i in 1:(p-1)){
  	for(j in (i+1):p){
  		if(methods$covar == 'wls'){
  			w <- sqrt(obs_w[,i] * obs_w[,j])
			# set minimal weights
			# to avoid numerical issues with evaluating the inverse of design matrix
			min_w <- ifelse(is.null(weight_th), 0, weight_th)
			w[w < min_w] <- min_w
  		}
  		P_2_w <- P_2 / w
  		y_w <- X_centered[, i] * X_centered[, j] / w

		design_m <- t(P_2_w) %*% P_2_w		
		inv_dm <- solve(design_m)
  		Sigma_array[j, i, ] <- Sigma_array[i, j, ] <- inv_dm %*% t(P_2_w) %*% y_w
  	}
  }
  return(lapply(1:K, function(k) Sigma_array[ , , k]))
}
