# load pseudo-bulk samples constructed from single cell data
load_bulk_sc_data <- function(gene_set, ct = 'Ast'){
    # -
    # load pre-computed pseudo-bulk data
    # -
    sc_pseudo_bulk <- readRDS(sprintf('%s/output/pseudo_bulk/bulk_scRNAseq_%s_nor_pseudo_bulk.rds', rosmap_data_dir, ct))
    subset_gene_exp <- sc_pseudo_bulk[match(gene_set, rownames(sc_pseudo_bulk)),]
    # set genes not available to exp=0
    print(sprintf('%i genes requested by the gene set are not available', sum(is.na(subset_gene_exp[,1]))))
    subset_gene_exp[is.na(subset_gene_exp)] <- 0
    rownames(subset_gene_exp) <- gene_set
    return(subset_gene_exp)
}

# load data
load_data <- function(genes){
    data_list <- list()
    data_list$X <- ROSMAP_bulk[['train']][match(genes, bulk_genes$external_gene_name),] %>% t
    colnames(data_list$X) <- genes
    data_list$X[is.na(data_list$X)] <- 0

    data_list$P <- train_props[, -which(colSums(train_props) == 0)]
    gene_mapping <- sapply(genes,function(g) sum(g == bulk_genes$external_gene_name, na.rm=T))
    print(table(gene_mapping))
    return(data_list)
}

load_GO_geneset <- function(geneset, low_filter = F, parse_style = 1, subset = NULL){
    # load genes in the gene set
    genes <- read.table(sprintf('%s/data/GO_pathway_genesets/%s', rosmap_data_dir, geneset),
                        skip=1, sep = '\n')[-1,] %>% as.vector()

    # subset expression data
    geneset_dat <- load_data(genes)
    print(dim(geneset_dat$X))

    # filter genes with low expression
    if(low_filter){
        gene_detected_count <- apply(geneset_dat$X, 2, function(x) sum(x > 0))

        print(sprintf('Cutoff for filtering: %f', nrow(geneset_dat$X) * 0.25))
        print(summary(gene_detected_count))

        filter_ind <- which(gene_detected_count < nrow(geneset_dat$X) * 0.25)

        if(length(filter_ind) > 0){
            print(sprintf('%i gene(s) were filtered', length(filter_ind)))
            geneset_dat$X <- geneset_dat$X[,-filter_ind]
            genes <- genes[-filter_ind]
        }
    }

    return(list(geneset_dat, genes))
}


# run WGCNA to obtain gene clusters and re-order the genes

run_wgcna <- function(data_list, geneset, power, to_plot = F){
    # pick power parameter
    log_nor_train_dat <- log(data_list$X + 1)

    # Choose a set of soft-thresholding powers
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    # Call the network topology analysis function
    sft = pickSoftThreshold(log_nor_train_dat, powerVector = powers, verbose = 5)

    if(to_plot){
    ## png(sprintf('%s/wgcna_bulk_softpower_%s.png', fig_prefix, ct_genesets[1]), width = 480, height = 480)
    # Plot the results:
    # sizeGrWindow(9, 5)
    par(mfrow = c(1,2));
    cex1 = 0.9;
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.90,col="red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    ## dev.off()
    }

    # perform clustering
    # one-step network construction
    net = blockwiseModules(log_nor_train_dat, power = power,
                           deepSplit=2,
    TOMType = "signed", minModuleSize = 3,
    reassignThreshold = 0, mergeCutHeight = 0.1,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = F,
    verbose = 3)

mergedColors = labels2colors(net$colors,
                                 colorSeq = c('#65C1E8',
                                              '#D85B63',
                                              '#D680AD',
                                              '#5C5C5C',
                                              '#C0BA80',
                                              '#FDC47D',
                                              '#EA3B46'))
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05)
    dev.off()

    moduleLabels = net$colors
    moduleColors = labels2colors(net$colors)
    MEs = net$MEs;
    geneTree = net$dendrograms[[1]];

    # output clustering results
    dendro_ordered_colors <- net$colors[net$dendrograms[[1]]$order]
    print(dendro_ordered_colors)
    cl_levels <- unique(dendro_ordered_colors) + 1
    bulk_cl_hc <- list(reordering = (net$dendrograms[[1]]$order)[order(dendro_ordered_colors)],
                       cluster = lapply(1:max(cl_levels), function(i) which(dendro_ordered_colors == i - 1)))
                     #cluster = dendro_ordered_colors[order(dendro_ordered_colors)] %>%
                     #  factor(levels = 0:max(dendro_ordered_colors), labels = 1:(max(dendro_ordered_colors)+1)))
    #saveRDS(bulk_cl_hc, sprintf('%s/%s_WGCNA_cl.rds', output_prefix, geneset))

    return(bulk_cl_hc)
}


# visualize the co-expression estimates in four most abundant brain cell types
plot_four_heatmaps <- function(est_list,
                               bulk_cl_hc,
                               geneset,
                               method_name,
                               ct_full_names=c('Excitatory neuron', 'Oligodendrocyte', 'Astrocyte', 'Microglia'),
                               plot_width=4.5,
                              cw = 5){
    g_list <- list()
    for(k in 1:4){
        g_list[[k]] <- plot_heatmap(est_list[[k]][bulk_cl_hc$reordering, bulk_cl_hc$reordering],
                     bulk_cl_hc$cluster,
                     ct_full_names[k],
                     legend = T, #ifelse(k %in% c(2,4), T, F),
                     annotation_legend = F,
                     annotation_names = F,
                     cw = cw,
                    width = plot_width,
                     silent = T)[[4]]
    }

    options(repr.plot.width = plot_width*2, repr.plot.height = plot_width*1.85)
    g <- grid.arrange(grobs = g_list, nrow = 2)
    ggsave(sprintf('%s/%s_%s.pdf', figure_dir, geneset, method_name), g,
          width = plot_width * 1.9, height = plot_width * 1.7)
}

real_data_run_bMIND <- function(data_list, ind_bulk, ind_profile){
    K <- ncol(data_list$P)
    genes <- colnames(data_list$X)

    bulk <- data_list$X %>% t
    rownames(data_list$P) <- colnames(bulk)

    sprintf('Cell type names match between cell type proportions and priors: %s',
            all(colnames(data_list$P) == colnames(prior$profile))) %>% print

    # for genes with valid prior inferred from single cell data
    # run bMIND with prior for these genes
    deconv_inf = bMIND(bulk = log2(1+bulk[ind_bulk, ]),
                       frac = data_list$P,
                       profile = prior$profile[ind_profile,],
                       covariance = prior$covariance[ind_profile,,],
                       ncore = 1)
    # for genes without valid prior inferred from single cell data
    # run bMIND with default mode without prior information
    deconv_np = bMIND(bulk = log2(1+bulk[-ind_bulk, ]),
                      frac = data_list$P,
                      ncore = NULL,
                      profile = NULL, # do not supply prior
                      covariance = NULL # do not supply prior
                      )
    # combine the estimates
    CTS_est <- array(0, dim = c(length(genes), K, ncol(bulk)))
    CTS_est[ind_bulk,,] <- deconv_inf$A
    CTS_est[-ind_bulk,,] <- deconv_np$A

    bMIND_est <- list()
    for(k in 1:K){
        bMIND_est[[k]] <- cor(CTS_est[,k,] %>% t)
        constant_exp_genes <- which(diag(bMIND_est[[k]]) == 0)
        bMIND_est[[k]][constant_exp_genes, ] <- NA
        bMIND_est[[k]][, constant_exp_genes] <- NA
    }
    return(bMIND_est)
}


