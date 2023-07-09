# organized from cross_cell_type_term_on_estimation.ipynb
# -
# load codes and data
# -

library(Matrix)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(nnls)

# load helper functions for this analysis
source('cross_ct_helper_functions.R')

data_dir <- '/gpfs/gibbs/pi/zhao/cs2629/ROSMAP/GeneExpression/snRNAseq_DLPFC_experiment2/pseudo_bulk'
source('../mom_ls.R')
source('../helper_function.R')

# load pseudo-bulk data
major_cts <- c('Excitatory Neurons', 'Oligodendrocytes',
              'Astrocyte', 'Microglia')

data_list <- lapply(major_cts, function(ct) load_data(ct))
names(data_list) <- major_cts


# -
# load cell type proportions
# -

# sample size per cell type
annot <- read.csv(sprintf('%s/../cell-annotation.csv', data_dir))
ct_by_ind_tab <- table(annot$individualID, annot$cell.type)
main_ct_tab <- ct_by_ind_tab[, colnames(ct_by_ind_tab) %in% major_cts]
main_ct_prop_tab <- main_ct_tab / rowSums(main_ct_tab)

main_ct_prop_tab <- main_ct_prop_tab[!is.na(main_ct_prop_tab[,1]), ]

features <- read.table(sprintf('%s/%s_features.txt', data_dir, major_cts[1]))[[1]]

# -
# compare two versions of CSNet
# -
gene_set_name_titles <- c('AD risk genes',
                         'Excitatory synapse genes',
                          'Myelin sheath genes',
                          'Astrocyte differentiation')
names(gene_set_name_titles) <- c('AD_risk_genes', 
                                'GOCC_EXCITATORY_SYNAPSE',
                                'GOCC_MYELIN_SHEATH',
                                'GOBP_ASTROCYTE_DIFFERENTIATION')
compare_two_CSNet('AD_risk_genes', 'ols', T)
compare_two_CSNet('AD_risk_genes', 'nnls', T)

compare_two_CSNet('GOCC_EXCITATORY_SYNAPSE', 'ols', T)
compare_two_CSNet('GOCC_EXCITATORY_SYNAPSE', 'nnls', T)

compare_two_CSNet('GOCC_MYELIN_SHEATH', 'ols', T)
compare_two_CSNet('GOCC_MYELIN_SHEATH', 'nnls', T)

compare_two_CSNet('GOBP_ASTROCYTE_DIFFERENTIATION', 'ols', T)
compare_two_CSNet('GOBP_ASTROCYTE_DIFFERENTIATION', 'nnls', T)
