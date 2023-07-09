library("optparse")
option_list = list(
  make_option(c("--i_geneset"), type="integer", default=1,
              help="index of GO gene set", metavar="integer")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

i_geneset <- opt$i_geneset

library(DirichletReg)
library(nnls)
library(ggplot2)
library(dplyr)
source('sensitivity_helper_functions.R')
source('/home/cs2629/project/CSNet/JASA_codes/CSNet_analysis/code/mom_ls.R')
source('/home/cs2629/project/CSNet/JASA_codes/CSNet_analysis/code/helper_function.R')
source('/home/cs2629/project/CSNet/JASA_codes/CSNet_analysis/code/cross_validation.R')
source('/home/cs2629/project/CSNet/JASA_codes/CSNet_analysis/code/simulation/visualization_helper.R')
source('real_data_helper_functions.R')

# -
# parameters
# -

# gene set of interests
ct_genesets <- c('GOCC_EXCITATORY_SYNAPSE',
                 'GOCC_MYELIN_SHEATH',
                 'GOBP_ASTROCYTE_DIFFERENTIATION',
		 'AD_risk_genes')
geneset_names <- c('excitatory synapse',
		   'myelin sheath',
		   'astrocyte differentiation',
		   'AD risk genes')

figure_dir <- 'figures/coexp_estimates'

# -
# load precomputed results from GO_gene_sets.R and AD_risk_genes.R
# -

geneset <- ct_genesets[i_geneset]
res <- readRDS(sprintf('output/%s_data.rds', geneset))
rownames(res$full_data_list[[1]]$P) <- rownames(res$full_data_list[[1]]$X)

CSNet_est <- readRDS(sprintf('../../data/ROSMAP/%s_CSNet_est.rds', geneset))
res$CSNet_est <- CSNet_est

# -
# generate Figure S12
# -

fn <- 'figures/sensitivity/noisy_pi_minus_pi.pdf'
if(!file.exists(fn)){
P <- res$full_data_list[[1]]$P
# add noises to observed proportions
kappas <- c(10,100,1000)
noisy_P_list <- lapply(kappas, function(kappa) sample_level_perturb(kappa, P))
# visualize boxplots of noises
P_diff_list <- lapply(noisy_P_list, function(nP) nP - P)
P_df <- do.call(rbind, P_diff_list) %>% as.data.frame
P_df$kappa <- rep(sprintf('kappa=%i', kappas), each = nrow(P))
P_df <- P_df[, c(1:4, 8)]
long_P_df <- melt(P_df, id = 'kappa')

ggplot(long_P_df) +
geom_boxplot(aes(x = variable, y = value)) +
facet_wrap(~kappa) +
labs(x = 'Cell type', y = 'Noisy pi - pi') +
theme_classic(base_size = 16.5)

ggsave('figures/sensitivity/noisy_pi_minus_pi.pdf', width = 6, height = 4)
}

# -
# generate Figure S13
# -

kappa_range <- c(10,100,1000)
sens_res <- real_data_sens(res, kappa_range, 100)
make_fig(sens_res, geneset, geneset_names[i_geneset])


# -
# generate Figure S14
# -

# randomly permute the proportions
data_list <- res$full_data_list[[1]]
set.seed(1)
n <- nrow(data_list$P)
data_list$P <- data_list$P[sample(n, n, replace = F), ]

# run CSNet with permuted proportions

# NNLS for variances and WLS for covariances
methods_settings <- list(var='nnls', covar='wls')
# threshold for small weights in WLS
th <- 0.01

dCSNet_cov_est <- mom_ls(data_list$P, data_list$X, methods_settings, th)
dCSNet_est <- lapply(dCSNet_cov_est, function(x) get_cor_from_cov(x))

cv_res <- cross_validation(data_list,
                           coexp_method = 'CSNet',
                           list(methods_settings = methods_settings, th = th),
                           ct_names = colnames(data_list$P),
                           to_plot = F,
                           seed = 1,
                          n_splits = 1)

CSNet_est <- lapply(1:length(dCSNet_est), function(i){
                            generalized_th(dCSNet_est[[i]], cv_res$th[i], 'SCAD', F) })
plot_four_heatmaps(CSNet_est, res$bulk_cl_hc, geneset, 'CSNet_permuted', cw = NULL)

