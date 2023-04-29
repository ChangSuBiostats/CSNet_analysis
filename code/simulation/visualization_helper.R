#
# visualization_helper.R
#
# plot heatmaps of co-expression matrices

plot_heatmap <- function(cor_m, sub_cl, title,
                         legend = T,
                         annotation_legend = T){
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

  my_colors = list(
    cluster = c('1'='#65C1E8',
                '2'='#D85B63',
                '3'='#D680AD',
                '4'='#5C5C5C',
                '5'='#C0BA80',
                '6'='#FDC47D',
                '7'='#EA3B46')[1:length(unique(labels))]
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
                    cellwidth = 2, cellheight = 2)
  return(h)
}
