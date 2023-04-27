# 
# evaluation_helper.R
#
# helper functions for evaluating the errors from different co-expression estimators
# in simulation

eval_errors <- function(R_hat, R, metrics, dense_estimate = F){
	error_vec <- numeric(length(metrics))
	names(error_vec) <- metrics
	for(m in metrics){
		if(m == 'F_norm'){
			error_vec[m] <- F_norm(R_hat, R)
		}else if(m == 'op_norm'){
			error_vec[m] <- op_norm(R_hat, R)
		}else if(m == 'TPR'){
			error_vec[m] <- eval_TPR(R_hat, R, off_diagonal = T)
		}else if(m == 'FPR'){
			error_vec[m] <- eval_FPR(R_hat, R, off_diagonal = T)
		}
	}
	if(dense_estimate) error_vec[c('TPR', 'FPR')] <- rep(NA, 2)
	return(error_vec)
}
