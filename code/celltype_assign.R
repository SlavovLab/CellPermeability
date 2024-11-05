## For getting file path
# install.packages("rstudioapi")
library(rstudioapi)
source(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),'/func.R'))


################################################################################################################################
# Read in Fresh+Frozen protein and mRNA data and preprocess for integration
################################################################################################################################


# Read in protein data
convert_mouse <- Proc_fasta('https://drive.google.com/uc?export=download&id=1B3q5YBQ2mxvoWivoj_brj58GQLSgSW3N')


# Read in Fresh data
protein_norm_noimp_fresh <- read.prot_gene('https://drive.google.com/uc?export=download&id=1-f-gtu8cdC4OuWphLIxpJzLQ8AWDILU2',convert_mouse)
protein_norm_imp_fresh <- read.prot_gene('https://drive.google.com/uc?export=download&id=12QcmfUXbHsFAQMzMjfkwLqwyWgFRBU1t',convert_mouse)
meta_fresh <- read.csv('https://drive.google.com/uc?export=download&id=1u6pqYc8ADxmPnmOxeZP4rappFjGi4uEo')


# Read in Frozen data 
protein_norm_noimp_frozen <- read.prot_gene('https://drive.google.com/uc?export=download&id=18SIiEbI3yINHh0oWKRSF4fHda1E_aH2Z',convert_mouse)
protein_norm_imp_frozen <- read.prot_gene('https://drive.google.com/uc?export=download&id=1uuSxzHfgPcNNYjpaAZIVPQ48i5R-6Gmq',convert_mouse)
meta_frozen <- read.csv('https://drive.google.com/uc?export=download&id=1rAJ_nU7XMY2LvXxYYhWNrWwQuDQmzNFM')



# **** Download the following data from this link: https://drive.google.com/file/d/1R9-utw-_rsZsNSxi1LccHCpmKF4ons45/view?usp=drive_link
trach_rnaseq <- readRDS('/Users/andrewleduc/Desktop/Senescense/seurat_integrated_filered_700_named.rds')

counts <- (trach_rnaseq@assays$RNA@counts)

# Transform counts and filter out genes with low # of observations
dim(counts)
counts <- as.matrix(counts)
counts <- counts[rowSums(counts != 0) > 100,]
dim(counts)



# Get intersection of genes across fresh and frozen protein data and mRNA data
intersect_genes <- intersect(rownames(protein_norm_imp_frozen),intersect(rownames(protein_norm_imp_fresh),rownames(counts)))
length(intersect_genes) # 1132 genes shared

# filter for intersected
protein_norm_imp_frozen <- protein_norm_imp_frozen[intersect_genes,]
protein_norm_imp_fresh <- protein_norm_imp_fresh[intersect_genes,]
counts_filt <- counts_filt[intersect_genes,]

intersect_genes


# make new scaledata
trach_rnaseq <- ScaleData(trach_rnaseq, do.scale = TRUE ,features = intersect_genes)
scale_data <- trach_rnaseq@assays[["RNA"]]@scale.data
scale_data <- scale_data[intersect_genes,]

# unload heavy file from RAM
rm(trach_rnaseq)


# This is correlation vector analysis, where we look for genes that have similar geneXgene correlation patterns
# which gives us a better subset of genes to feed in to LIGER for integration

imputed_cor_vect_frozen <- GetCorVect(protein_norm_imp_frozen,scale_data)
imputed_cor_vect_fresh <- GetCorVect(protein_norm_imp_fresh,scale_data)


# Check protein order scrambled null dist is 0 centered
#imputed_cor_vect_null <- GetCorVect(protein_norm_imp[sample(rownames(protein_norm_imp)),],scale_data)

# Check results 
ggplot(imputed_cor_vect_frozen, aes(x = cors)) + geom_histogram()  +
  ylab('# of Genes') + xlab('Correlations') + ggtitle('Correlation vector comparison')
ggplot(imputed_cor_vect_fresh, aes(x = cors)) + geom_histogram()  +
  ylab('# of Genes') + xlab('Correlations') + ggtitle('Correlation vector comparison')
#ggplot(imputed_cor_vect_null, aes(x = cors)) + geom_histogram()


# Compare fresh and frozen -- Supplemental figure 1a
imputed_cor_vect_frozen$sample <- 'TMT Frozen'
imputed_cor_vect_fresh$sample <- 'TMT Fresh'
imputed_cor_vect_compare <- rbind(imputed_cor_vect_frozen,imputed_cor_vect_fresh)
ggplot(imputed_cor_vect_compare,aes(x = sample, y = cors)) + geom_boxplot()+ xlab('') +
  ylab('Correlations') + ggtitle('Correlation vector analysis (intersected Genes)')+
  dot_plot


# Filter for genes that covary 
imputed_cor_vect_good <- imputed_cor_vect_frozen %>% filter(cors > .2)
keep_corV <- imputed_cor_vect_good$genes #intersect(unimputed_cor_vect_good$genes,imputed_cor_vect_good$genes)


protein_norm_imp_fresh_filt <- protein_norm_imp_fresh[keep_corV,]
protein_norm_imp_frozen_filt <- protein_norm_imp_frozen[keep_corV,]

scale_dat_filt <- scale_data[keep_corV,]
counts <- counts[keep_corV,]




#### Integration plots showing mixing, cell type amounts per clust, and age per clust ####

# Generating confidence interval for integration score by boostrapping data
# score_fresh <- c()
# score_frozen <- c()
# for(i in 1:2){
#   
#   integrated_fresh <- Integrate_liger(counts[,sample(ncol(counts),2000)],protein_norm_imp_fresh_filt[,sample(ncol(protein_norm_imp_fresh_filt),500)],10)
#   integrated_frozen <- Integrate_liger(counts[,sample(ncol(counts),2000)],protein_norm_imp_frozen_filt[,sample(ncol(protein_norm_imp_frozen_filt),500)],10)
#   
#   score_fresh <- c(score_fresh,integrated_fresh[[1]])
#   score_frozen <- c(score_frozen,integrated_frozen[[1]])
#   
# }
# 
# median(score_fresh)
# median(score_frozen)



################################################################################################################################
# Integration and cell type assignment fresh samples
################################################################################################################################


integrated_fresh <- Integrate_liger(counts_filt,protein_norm_imp_fresh_filt,10)
integrated_fresh[[1]]
umap_liger_fresh <-integrated_fresh[[2]]
umap_liger_prot <- umap_liger_fresh %>% filter(Dataset == 'prot')
umap_liger_rna <- umap_liger_fresh %>% filter(Dataset == 'mRNA')


ggplot(umap_liger, aes(x =Dim1,y = Dim2,color = Dataset )) + geom_point()+dot_plot
ggplot(umap_liger_prot, aes(x =Dim1,y = Dim2,color = Cluster )) + geom_point()+dot_plot


## Transfer cell type assignment for co-clusters from preannotated rna seq data


  # Some useful marker proteins for additional confidence
  # Epithelial: Sult1d1, Gsta4
  # Basal: Ckmt1, Cyp2f2
  # Secretory: Aldh1a7, Cyp2f2, Cbr2
  # Ciliated: Clic1
  # Fib/immune: Vim, S100a4
  # Immune: Lcp1
  # Chondrocyte: Eno1

FeaturePlot(trach_rnaseq, features = 'Krt8' )

umap_liger_prot <- umap_liger_prot[colnames(protein_norm_imp_fresh_filt),]
umap_liger_prot$ID <- rownames(umap_liger_prot)
umap_liger_prot <- umap_liger_prot %>% left_join(meta,by = c('ID'))

umap_liger_prot$prot <- protein_norm_imp_fresh_filt['Krt8',]
ggplot(umap_liger_prot, aes(x =Dim1,y = Dim2,color = prot )) + geom_point() +dot_plot +
  scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white',high = 'red')


umap_liger_prot$dead <- 'no'
umap_liger_prot$dead[umap_liger_prot$Stain_Diameter > 10] <- 'yes'
ggplot(umap_liger_prot, aes(x =Dim1,y = Dim2,color = dead )) + geom_point() +dot_plot


umap_liger_prot$Cell_type <- NA
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 0] <- 'Basal'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 1] <- 'Fibroblast'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 2] <- 'Fibroblast'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 3] <- 'Club'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 4] <- 'Chondrocyte'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 5] <- 'Immune'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 6] <- 'Ciliated'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 7] <- 'Fibroblast'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 8] <- 'Club'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 9] <- 'not sure'


umap_liger_prot <- umap_liger_prot %>% filter(Cell_type != 'not sure')

write.csv(umap_liger_prot,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/celltype_annotations_fresh.csv')








################################################################################################################################
# Integration and cell type assignment Frozen samples
################################################################################################################################


integrated_fresh <- Integrate_liger(counts_filt,protein_norm_imp_frozen_filt,10)
integrated_fresh[[1]]
umap_liger_fresh <-integrated_fresh[[2]]
umap_liger_prot <- umap_liger_fresh %>% filter(Dataset == 'prot')
umap_liger_rna <- umap_liger_fresh %>% filter(Dataset == 'mRNA')


ggplot(umap_liger, aes(x =Dim1,y = Dim2,color = Dataset )) + geom_point()+dot_plot
ggplot(umap_liger_prot, aes(x =Dim1,y = Dim2,color = Cluster )) + geom_point()+dot_plot


## Transfer cell type assignment for co-clusters from preannotated rna seq data


# Some useful marker proteins for additional confidence
# Epithelial: Sult1d1, Gsta4
# Basal: Ckmt1, Cyp2f2
# Secretory: Aldh1a7, Cyp2f2, Cbr2
# Ciliated: Clic1
# Fib/immune: Vim, S100a4
# Immune: Lcp1
# Chondrocyte: Eno1

FeaturePlot(trach_rnaseq, features = 'Krt8' )

umap_liger_prot <- umap_liger_prot[colnames(protein_norm_imp_frozen_filt),]
umap_liger_prot$ID <- rownames(umap_liger_prot)
umap_liger_prot <- umap_liger_prot %>% left_join(meta,by = c('ID'))

umap_liger_prot$prot <- protein_norm_imp_frozen_filt['Krt8',]
ggplot(umap_liger_prot, aes(x =Dim1,y = Dim2,color = prot )) + geom_point() +dot_plot +
  scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white',high = 'red')


umap_liger_prot$dead <- 'no'
umap_liger_prot$dead[umap_liger_prot$Stain_Diameter > 10] <- 'yes'
ggplot(umap_liger_prot, aes(x =Dim1,y = Dim2,color = dead )) + geom_point() +dot_plot


umap_liger_prot$Cell_type <- NA
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 0] <- 'Basal'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 1] <- 'Fibroblast'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 2] <- 'Fibroblast'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 3] <- 'Club'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 4] <- 'Chondrocyte'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 5] <- 'Immune'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 6] <- 'Ciliated'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 7] <- 'Fibroblast'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 8] <- 'Club'
umap_liger_prot$Cell_type[umap_liger_prot$Cluster == 9] <- 'not sure'




write.csv(umap_liger_prot,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/celltype_annotations_frozen.csv')




