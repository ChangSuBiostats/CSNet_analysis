library(dplyr)
library(optparse)
library(ROCR)
library(scales)
source('../helper_function.R')

option_list = list(
  make_option(c("--result_prefix"), type="character",
              default="results/MA_rho_0.5_0.5/K_2/log_var_8.0_equal_strength_TRUE/n_150_p_100",
              help="which experiment to print", metavar="character"),
  make_option(c("--setting"), type="character",
              default="standard",
              help="which set of methods to print", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
result_prefix <- opt$result_prefix
setting <- opt$setting

n_p_setting <- strsplit(result_prefix, '/')[[1]][5]
n <- strsplit(n_p_setting, '_')[[1]][2] %>% as.integer
p <- strsplit(n_p_setting, '_')[[1]][4] %>% as.integer
K <- strsplit(strsplit(result_prefix, '/')[[1]][3], '_')[[1]][2] %>% as.numeric
n_rep <- 200

fig_prefix <- paste0('figures', strsplit(result_prefix, 'results')[[1]][2])

# -
# load ground truths
# -
sim_setting <- readRDS(sprintf('%s/sim_setting.rds', result_prefix))

# -
# load co-expression estimates
# -
if(setting == 'full'){
  coexp_methods <- c('d-Bulk', 'd-CSNet', 'bMIND', 'ENIGMA')
  names(coexp_methods) <- c('Bulk', 'CSNet', 'bMIND', 'ENIGMA')
}else if(setting == 'standard'){
  coexp_methods <- c('d-Bulk', 'd-CSNet', 'bMIND')
  names(coexp_methods) <- c('Bulk', 'CSNet', 'bMIND')
}
n_methods <- length(coexp_methods)

pred_list <- lapply(1:K, function(k) lapply(1:length(coexp_methods), function(i) list()))
for(k in 1:K) names(pred_list[[k]]) <- coexp_methods
CSNet_est <- list()

for(i_rep in 1:n_rep){
  est_tmp <- readRDS(sprintf('%s/n_rep_%i_i_rep_%i_est.rds', result_prefix, 200, i_rep))
  for(coexp_m in coexp_methods){
    for(k in 1:K){
      coexp_est <- est_tmp[[coexp_m]][[k]]
      pred_list[[k]][[coexp_m]][[i_rep]] <- coexp_est[upper.tri(coexp_est)] %>% abs
    }
  }
  CSNet_est[[i_rep]] <- est_tmp[['CSNet']]
}

perf_list <- list()
selected_points <- matrix(NA, nrow = n_rep, ncol = 2)
for(k in 1:K){
  # extract entries with true co-expression !=0
  R_star <- sim_setting$R[[k]]
  truths <- list()
  for(i_rep in 1:n_rep) truths[[i_rep]] <- as.numeric(R_star[upper.tri(R_star)] > 0)
  # evaluate the ROC curves for different methods
  for(coexp_m in coexp_methods){
    #print(sapply(pred_list[[k]][[coexp_m]], length))
    #print(table(truths[[1]]))
    #print(coexp_m)
    pred <- prediction(pred_list[[k]][[coexp_m]], truths)
    perf_list[[coexp_m]] <- performance(pred,
                            'tpr', 'fpr')
  }
  # for CSNet, evaluate the TPR and FPR at the selected threshold
  seleced_points <- matrix(NA, nrow = n_rep, ncol = 2)
  for(i_rep in 1:n_rep){
    selected_points[i_rep, 1] <- eval_FPR(CSNet_est[[i_rep]][[k]], R_star, T)
    selected_points[i_rep, 2] <- eval_TPR(CSNet_est[[i_rep]][[k]], R_star, T)
  } 
  # make AUC plot
  # https://mycolor.space/?hex=%23CA6277&sub=1
  color_palette <- c('#0090DB', '#CA6277', '#52922C', '#835B8A')
  names(color_palette) <- c('Bulk', 'CSNet', 'bMIND', 'ENIGMA')
  lty_list <- c('dotted', 'solid', 'dotdash', 'longdash')

  if(n_methods < 4){
    method_match <- match(names(coexp_methods), names(color_palette))
    color_palette <- color_palette[method_match]
    lty_list <- lty_list[method_match]
  }

  # Code reference: https://mlr.mlr-org.com/articles/tutorial/roc_analysis.html
  # https://preludeinr.com/graphs-and-plots/adding-elements-to-an-existing-graph/
  pdf(sprintf('%s/ROC_curves_ct_%i_%s.pdf', fig_prefix, k, setting), width = 5, height = 5)
  par(cex.axis=1, cex.lab=1.2, cex.main=1.8, mar=c(4.5,4.5,5,1))
  # Bulk
  plot(perf_list[[1]],
        avg='vertical',
        spread.estimate='stderror',
        lwd=3,
        # main=sprintf('(n,p)=(%i,%i), m=%s', n, p, ifelse(k==1, '2/3', '1/3')),
        main = sprintf('p=%i', p),
        lty = lty_list[1],
        ylab='True Positive Rate',
        xlab='False Positive Rate',
        #ylab = 'Average true positive rate (1 standard error)',
        col=color_palette[1])  
  
  plot(perf_list[[2]],
        avg='vertical',
        spread.estimate='stddev', #'stderror'
        lty = lty_list[2],
        lwd=3,
        col=color_palette[2], ylab='',
        add = TRUE)

  plot(perf_list[[3]],
          avg='vertical',
          spread.estimate='stderror',
          lty = lty_list[3],
          lwd=3,
          col=color_palette[3], ylab='',
          add = TRUE)
  
  if(setting == 'full'){
    plot(perf_list[[4]],
         avg='vertical',
         spread.estimate='stderror',
         lty = lty_list[4],
         lwd=3,
         col=color_palette[4], ylab='',
         add = TRUE)
  }
  legend("bottomright",
           names(color_palette),
           lty=lty_list,
           col=color_palette,
           bty="o",               # The box type for the whole legend
           lwd=4,
           pt.cex=1.4,            # Expansion factor for symbols
           cex=1.55)

  points(selected_points[,1],
         selected_points[,2],
         cex=1.25,
         cex.axis=1.5,
         pch=4, col=alpha(color_palette[2],0.5))
  
  n_setting <- sprintf('n=%i', n)
  text(0.8, ifelse(setting == 'standard', 0.45, 0.55), n_setting, cex = 1.8, font = 2) 
  dev.off()

}

print(sprintf('%s/ROC_curves_%s.pdf', fig_prefix, setting))
