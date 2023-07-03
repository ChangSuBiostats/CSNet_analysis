#
# visualization_helper.R
#
# plot heatmaps of co-expression matrices

plot_heatmap <- function(cor_m, sub_cl, title,
                         legend = T,
                         annotation_legend = T,
			 annotation_names = T,
			 silent = F,
			 cw = 2){
  require(pheatmap)
  require(RColorBrewer)
  # color for heatmap
  palPos <- colorRampPalette(c("white", "red"), space = "rgb")
  palNeg <- colorRampPalette(c("blue", "white"), space = "rgb")
  coul <- c(palNeg(20), palPos(20))
  
  # label for clusters
  g_K <- length(sub_cl)
  p <- nrow(cor_m)
  K <- length(sub_cl)
  labels <- sapply(1:g_K, function(k) rep(k, length(sub_cl[[k]]))) %>% unlist
  labels <- c(labels, rep(g_K + 1, p - length(labels))) 
  
  # name corrlation matrix
  if(is.null(rownames(cor_m))){
      rownames(cor_m) <- colnames(cor_m) <- paste0('Gene_', 1:p)
  }
  
  # http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#:~:text=The%20default%20colors%20in%20ggplot2,not%20friendly%20for%20colorblind%20viewers.
  my_colors = list(
    cluster = c('1'='#E69F00',
                '2'='#56B4E9',
                '3'='#009E73',
                '4'='#F0E442',
                '5'='#0072B2',
                '6'='#D55E00',
                '7'='#CC79A7')[1:length(unique(labels))]
  )
  gene_cl <- data.frame(cluster = labels %>% as.character)
  rownames(gene_cl) <- rownames(cor_m)

  h <- pheatmap(cor_m[rev(1:p), ], color = coul,cluster_rows = F, cluster_cols = F,
                breaks = seq(-1, 1, by = 0.05),
  main = (title),
  fontsize = 16,
  na_col = 'grey',
  show_rownames = F, show_colnames = F,
  border_color = FALSE,
  annotation_row = gene_cl,
                    annotation_col = gene_cl,
                    annotation_colors = my_colors,
                    annotation_legend = annotation_legend,
                    legend = legend,
		    annotation_names_row = annotation_names, annotation_names_col = annotation_names,
                    cellwidth = cw, cellheight = cw,
		    silent = silent)
  h$gtable$grobs[[1]]$gp <- gpar(fontsize = 25, fontfamily = 'sans')
  return(h)
}
