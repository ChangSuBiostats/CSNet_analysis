# load gene annotations for bulk RNA-seq data from Mostafavi et al.
bulk_genes <- readRDS(sprintf('%s/output/ROSMAP_annotLookup_new.rds', data_rosmap_dir))

# AD gene list from a 2020 Nature neuroscienc publication
# https://www.nature.com/articles/s41593-020-0599-5/tables/1

MS4A_cluster <- bulk_genes$external_gene_name[grepl('MS4A', bulk_genes$external_gene_name)]
# MS4A_cluster[!MS4A_cluster %in% rownames(ROSMAP_sc)]

AD_genes <- c('APOE', 'EPHA1', 'CLU', 'INPP5D',
              'HLA-DRB5', 'HLA-DRB1', #or
             'CR1', 'TREM2', 'CD33',
              MS4A_cluster,
              'ABI3', 'PLCG2',
              'ZCWPW1','PILRA',
              'MEF2C', 'CD2AP', 'BIN1', 'PICALM', 'CASS4',
              'CELF1', 'SPI1',
             'FERMT2', 'NME8', 'SORL1', 'ABCA7', 'SLC24A4–RIN3', 'PTK2B',
             'ADAM10', 'IGHV1-67', 'PPARGC1A', 'TP53INP1', 'ECHDC3', 'ACE',
             'ADAMTS1', 'IQCK', 'TRIP4', 'RORA', 'ZNF423', 'APP', 'IGHG3',
             'AC099552.4', 'ZNF655',
              'HBEGF', 'AFDN1',
              'BZRAP1-AS1', 'TPBG','DSG2',
              'CLNK', 'HS3ST1',
              'SCIMP',
              'PRKD3', 'NDUFAF7',
              'TREML2','SHARPIN',
              'MAPT', 'KANSL1',
              'CHRNE', 'C17orf107',
              'IL34',
             'CNTNAP2', 'ALPK2', 'ADAMTS4', 'APH1B', 'KAT8', 'SPPL2A', 'HESX1')
# AD genes that have keyword 'immune' or 'inflam' in their functions
immune_AD_genes <- c('APOE', 'EPHA1', 'CLU', 'INPP5D', 'HLA-DRB5', 'HLA-DRB1',
             'CR1', 'TREM2', 'CD33', MS4A_cluster, 'ABI3', 'PLCG2', 'MEF2C',
                    'PTK2B', 'RORA', 'ADAMTS1', 'IGHG3', 'CLNK', 'HS3ST1',
                  'TREML2', 'SHARPIN')

length(AD_genes)
length(immune_AD_genes)
length(MS4A_cluster)

81-15-8

# ----------------------------------------------------------------
# check if genes at the same locus are covered by ROSMAP bulk data
# ----------------------------------------------------------------
c('HLA-DRB5', 'HLA-DRB1') %in% bulk_genes$external_gene_name
c('ZCWPW1', 'PILRA') %in% bulk_genes$external_gene_name
c('CELF1', 'SPI1') %in% bulk_genes$external_gene_name
c('HBEGF', 'AFDN1') %in% bulk_genes$external_gene_name
c('CLNK', 'HS3ST1') %in% bulk_genes$external_gene_name
c('PRKD3', 'NDUFAF7') %in% bulk_genes$external_gene_name
c('MAPT','KANSL1') %in% bulk_genes$external_gene_name
c('CHRNE', 'C17orf107') %in% bulk_genes$external_gene_name

# ----------------------------------------------------------------
# check if all AD genes are covered by the bulk ROSMAP data
# if not, try to substitute with equivalent names that are available
# ----------------------------------------------------------------
print(AD_genes %>% length)
#print(AD_genes[!AD_genes %in% rownames(ROSMAP_sc)])
print(AD_genes[!AD_genes %in% bulk_genes$external_gene_name])

# SLC24A4–RIN3
'SLC24A4' %in% bulk_genes$external_gene_name
# "AC099552.4"
'ENSG00000278998' %in% bulk_genes$external_gene_name
# not found
# AFDN1
# not found
# BZRAP1-AS1
'TSPOAP1-AS1' %in% bulk_genes$external_gene_name

AD_genes[AD_genes == "SLC24A4–RIN3"] <- "SLC24A4"
AD_genes[AD_genes == "BZRAP1-AS1"] <- "TSPOAP1-AS1"

# ----------------------------------------------------------------
# subset AD genes that are covered by bulk ROSMAP data
# ----------------------------------------------------------------
print(AD_genes[!AD_genes %in% bulk_genes$external_gene_name])
avail_AD_genes <- AD_genes[AD_genes %in% bulk_genes$external_gene_name]
# avail_AD_genes <- AD_genes[AD_genes %in% rownames(ROSMAP_sc)]

print(sprintf('#AD genes: %i', length(AD_genes)))
print(sprintf('#AD genes available in ROSMAP bulk data: %i', length(avail_AD_genes)))
