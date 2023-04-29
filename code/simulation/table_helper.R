#
# table_helper.R
#
# helper functions for printing tables
print_a_setting <- function(error_list, coexp_methods, n, p){
  n_methods <- length(coexp_methods)
  if(p == 100){
    cat('\\\\\\hline') 
    cat('\n')
    cat(sprintf('\\multirow{%i}{*}{%i} & \\multirow{%i}{*}{\\begin{tabular}[c]{@{}l@{}}$%i$\\end{tabular}}',
        4 * n_methods, n, 2 * n_methods, p))
  }else{
    cat('\\\\\\cline{2-11}')
    cat('\n')
    cat(sprintf('& \\multirow{%i}{*}{\\begin{tabular}[c]{@{}l@{}}$%i$\\end{tabular}}',
        2 * n_methods, p))
  }
  cat('\n')
  for(coexp_method in coexp_methods){
    print_a_method(error_list[[coexp_method]], coexp_method)
  }
  cat('\\\\\\hline')
}

print_a_method <- function(error_mat, coexp_method){
  K <- length(error_mat)
  if(coexp_method != 'Bulk'){
    cat('\\\\\\cline{3-11}')
    cat('\n')
    cat(sprintf('&& \\texttt{%s}', coexp_method))
  }else{
    cat('& \\texttt{Bulk}')
  }
  cat('\n')
  #for(j in 1:ncol(error_mat[[1]])){
  for(k in 1:K){
    for(j in 1:ncol(error_mat[[1]])){
    err_mean <- mean(error_mat[[k]][,j])
    err_sd <- sd(error_mat[[k]][,j])
    if(coexp_method == 'CSNet'){
      cat('& \\begin{tabular}[c]{@{}c@{}}',
          ifelse(is.na(err_mean), '{-}', sprintf('\\bf{%.2f}', err_mean)),
          '\\\\\\small{',
          ifelse(is.na(err_sd), '{-}', sprintf('\\bf{(%.2f)}', err_sd)),
          '}\\end{tabular}', sep='')
    }else{
      cat('& \\begin{tabular}[c]{@{}c@{}}',
          ifelse(is.na(err_mean), '{-}', sprintf('%.2f', err_mean)),
          '\\\\\\small{',
          ifelse(is.na(err_sd), '{-}', sprintf('(%.2f)', err_sd)),
          '}\\end{tabular}', sep='')
    }
    }
  }
  cat('\n')
}
