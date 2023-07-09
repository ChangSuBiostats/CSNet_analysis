# -
# load packages and codes
# -

library(dplyr)
library(nnls)

# for visualization
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)

source('/home/cs2629/project/CSNet/JASA_codes/CSNet_analysis/code/mom_ls.R')
source('/home/cs2629/project/CSNet/JASA_codes/CSNet_analysis/code/helper_function.R')
source('/home/cs2629/project/CSNet/JASA_codes/CSNet_analysis/code/simulation/visualization_helper.R')

genesets <- c('GOCC_EXCITATORY_SYNAPSE',
              'GOCC_MYELIN_SHEATH',
              'GOBP_ASTROCYTE_DIFFERENTIATION')

GO_data <- list()
for(geneset in genesets){
    GO_data[[geneset]] <- readRDS(sprintf('output/%s_data.rds', geneset))$full_data_list[[1]]
}

# -
# combine three gene sets
# -
print('1-2')
intersect(colnames(GO_data[[1]]$X), colnames(GO_data[[2]]$X))
print('2-3')
intersect(colnames(GO_data[[3]]$X), colnames(GO_data[[2]]$X))
print('1-3')
intersect(colnames(GO_data[[1]]$X), colnames(GO_data[[3]]$X))

data_combined <- cbind(GO_data[[1]]$X, GO_data[[2]]$X, GO_data[[3]]$X)
dim(data_combined)
gn_table <- table(colnames(data_combined))
redun_genes <- names(gn_table)[which(gn_table > 1)]
redun_genes

data_combined <- cbind(GO_data[[1]]$X,
                       GO_data[[2]]$X[, !colnames(GO_data[[2]]$X) %in% redun_genes],
                       GO_data[[3]]$X[, !colnames(GO_data[[3]]$X) %in% redun_genes])
dim(data_combined)

geneset_memb <- rep(c('Ex', 'Oli', 'Ast'), c(ncol(GO_data[[1]]$X),
                                             ncol(GO_data[[2]]$X)-length(redun_genes),
                                             ncol(GO_data[[3]]$X)-length(redun_genes)))

geneset_memb_num <- rep(1:3,
                        c(ncol(GO_data[[1]]$X),
                          ncol(GO_data[[2]]$X)-length(redun_genes),
                          ncol(GO_data[[3]]$X)-length(redun_genes)))
geneset_memb_list <- sapply(1:3, function(i) which(geneset_memb_num == i))


# -
# evaluate d-CSNet and d-Bulk
# -

ROSMAP_data_list <- list(X = data_combined,
                        P = GO_data[[1]]$P)
methods_settings <- list(var='nnls', covar='wls')

dCSNet_cov_est <- mom_ls(ROSMAP_data_list$P, ROSMAP_data_list$X, methods_settings)
dCSNet_est <- lapply(dCSNet_cov_est, function(x) get_cor_from_cov(x))

# obtain bulk estimates
dbulk_est <- cor(ROSMAP_data_list$X)


# -
# visualize the estimates
# -

# reorder the genes based on clustering for better visualization
geneset_memb_list_new <- list()
for(i in 1:3){
    hc_res <- hclust(dist(dCSNet_est[[i]][geneset_memb_list[[i]], geneset_memb_list[[i]]]),
                    method = 'average')
    geneset_memb_list_new[[i]] <- geneset_memb_list[[i]][hc_res$order]
}

ct_names <- c('Excitatory neuron', 'Oligodendrocyte', 'Astrocyte')

g_list <- list()
for(i in 1:3){
    new_ordering <- c(geneset_memb_list_new[[1]],
                     geneset_memb_list_new[[2]],
                     geneset_memb_list_new[[3]])
    g_list[[i]] <- plot_heatmap(dCSNet_est[[i]][new_ordering, new_ordering],
                                geneset_memb_list,
                                ct_names[i],
                         legend = T,#ifelse(i == 3, T, F),
                         annotation_legend = F,
                         annotation_names = F,
                                width = 5,
                         silent = T)[[4]]
}


g_list[[4]] <- plot_heatmap(dbulk_est[new_ordering, new_ordering],
                            geneset_memb_list,
                            'Bulk brain',
                         legend = ifelse(i == 3, T, F),
                         annotation_legend = F,
                         annotation_names = F,
                            width = 5,
                         silent = T)[[4]]

# make the legend

legend_df <- data.frame(GO = c('Excitatory synapse', 'Myelin sheath', 'Astrocyte differentiation'), values = 1:3)
legend_df$GO <- factor(legend_df$GO,
                       levels = c('Excitatory synapse', 'Myelin sheath', 'Astrocyte differentiation'),
                       labels = c('Excitatory synapse', 'Myelin sheath', 'Astrocyte differentiation'))

legend_plot <- ggplot(legend_df) +
geom_point(aes(x = GO, y = values, color = GO), shape = 15, size = 10) +
theme_classic(base_size = 25) +
labs(color = 'GO gene sets   ') +
theme(legend.position="bottom") +
scale_color_manual(values = c('#E69F00', '#56B4E9', '#009E73'), breaks = legend_df$GO) +
guides(color=guide_legend(nrow=1,byrow=TRUE)) +
#theme(legend.box.margin = margin(0, 0, 0, 12))
theme(legend.spacing.x = unit(0.55, 'cm'))

g_list[[5]] <- as_ggplot(get_legend(legend_plot))

pw <- 4.5
options(repr.plot.width = (pw*1.1)*4, repr.plot.height =pw+0.6)
g <- grid.arrange(grobs = g_list, nrow = 2,
                  widths=rep(pw, 4),
                  heights = c(pw+0.5,0.8),
                  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 5, 5, 5)))

ggsave('figures/E_P3_three_gene_sets_combined_reordered_with_bulk.pdf', g,
      width = (pw*1.1)*4, height = pw+0.6)
