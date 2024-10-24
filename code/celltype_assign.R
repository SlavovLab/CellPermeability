source('/Users/andrewleduc/Desktop/Senescense/Trachea/TMT/Joint_mRNA_prot_functions.R')
library(stringr)
library(sva)
library(dplyr)
library(ggplot2)
# Read in protein data

convert_mouse <- Proc_fasta('/Users/andrewleduc/Desktop/Senescense/Mouse.fasta')


# Fresh
protein_norm_noimp <- read.prot_gene('/Users/andrewleduc/Desktop/Projects/AliveDead/Fresh/protein.csv',convert_mouse)
protein_norm_imp <- read.prot_gene('/Users/andrewleduc/Desktop/Projects/AliveDead/Fresh/protein.imputed.csv',convert_mouse)
meta <- read.csv()


# Frozen
protein_norm_noimp <- read.prot_gene('/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/protein.csv',convert_mouse)
protein_norm_imp <- read.prot_gene('/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/protein.imputed.csv',convert_mouse)
meta <- read.csv('/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/meta.csv')




protein_norm_imp <- as.matrix(protein_norm_imp)
#protein_norm_noimp <- as.matrix(protein_norm_noimp)
trach_gergana <- readRDS('/Users/andrewleduc/Desktop/Senescense/seurat_integrated_filered_700_named.rds')

counts <- (trach_gergana@assays$RNA@counts)#[intersect(rownames(protein_norm_imp),rownames(counts)),drop = FALSE]
counts <- counts[intersect(rownames(protein_norm_imp),rownames(counts)),]

dim(counts)
counts_filt <- as.matrix(counts)




counts_filt <- counts_filt[rowSums(counts_filt != 0) > 40,]
#protein_norm_noimp <- filt.mat.cr(protein_norm_noimp,.95,.95)#Filter_missingness(protein_norm_noimp,clust)

intersect_genes <- intersect(rownames(protein_norm_imp),rownames(counts_filt))#intersect(rownames(protein_norm_noimp),rownames(counts_filt))

# filter for intersected
#protein_norm_noimp <- protein_norm_noimp[intersect_genes,]
protein_norm_imp <- protein_norm_imp[intersect_genes,]
counts_filt <- counts_filt[intersect_genes,]

intersect_genes


# make new scaledata
trach_gergana <- ScaleData(trach_gergana, do.scale = TRUE ,features = intersect_genes)
scale_data <- trach_gergana@assays[["RNA"]]@scale.data
scale_data <- scale_data[intersect_genes,]


# I just ran the script twice and loaded each data set seperately, now saving each for comparison
imputed_cor_vect1 <- GetCorVect(protein_norm_imp,scale_data)
imputed_cor_vect2 <- GetCorVect(protein_norm_imp,scale_data)

#unimputed_cor_vect <- GetCorVect(protein_norm_noimp,scale_data)
#unimputed_cor_vect <- GetCorVect(protein_norm_noimp,scale_data)

# Check protein order scrambled null dist is 0 centered
imputed_cor_vect_null <- GetCorVect(protein_norm_imp[sample(rownames(protein_norm_imp)),],scale_data)

# Check results
ggplot(imputed_cor_vect2, aes(x = cors)) + geom_histogram()  +
  ylab('# of Genes') + xlab('Correlations') + ggtitle('Correlation vector comparison')
ggplot(imputed_cor_vect_null, aes(x = cors)) + geom_histogram()


# Compare how cell death affects overall agreement to mRNA data
imputed_cor_vect1$sample <- 'TMT Frozen'
imputed_cor_vect2$sample <- 'TMT Fresh'
imputed_cor_vect_compare <- rbind(imputed_cor_vect1,imputed_cor_vect2)
ggplot(imputed_cor_vect_compare,aes(x = sample, y = cors)) + geom_boxplot()+ xlab('') +
  ylab('Correlations') + ggtitle('Correlation vector analysis (intersected Genes)')+
  dot_plot




#unimputed_cor_vect_good <- unimputed_cor_vect %>% filter(cors > .3)
imputed_cor_vect_good <- imputed_cor_vect1 %>% filter(cors > .2)


keep_corV <- imputed_cor_vect_good$genes #intersect(unimputed_cor_vect_good$genes,imputed_cor_vect_good$genes)

protein_norm_imp_filt <- protein_norm_imp[keep_corV,]
scale_dat_filt <- scale_data[keep_corV,]
counts_filt <- counts_filt[keep_corV,]


#### Integration plots showing mixing, cell type amounts per clust, and age per clust ####

integrated <- Integrate_liger(counts_filt,protein_norm_imp_filt,10)
integrated[[1]]
umap_liger <-integrated[[2]]
umap_liger_prot <- umap_liger %>% filter(Dataset == 'prot')
umap_liger_rna <- umap_liger %>% filter(Dataset == 'mRNA')


ggplot(umap_liger, aes(x =Dim1,y = Dim2,color = Dataset )) + geom_point()+dot_plot
ggplot(umap_liger_prot, aes(x =Dim1,y = Dim2,color = Cluster )) + geom_point()+dot_plot


## Manual cell type assignment for co-clusters from preannotated rna seq data

#Find marker proteins for cell type anno

# Epithelial: Sult1d1, Gsta4
    # Basal: Ckmt1, Cyp2f2
    # Secretory: Aldh1a7, Cyp2f2, Cbr2
    # Ciliated: Clic1
# Fib/immune: Vim, S100a4
# Immune: Lcp1
# Chondrocyte: Eno1


FeaturePlot(trach_gergana, features = 'Krt8' )

umap_liger_prot <- umap_liger_prot[colnames(protein_norm_imp_filt),]
umap_liger_prot$ID <- rownames(umap_liger_prot)
umap_liger_prot <- umap_liger_prot %>% left_join(meta,by = c('ID'))

umap_liger_prot$prot <- protein_norm_imp_filt['Krt8',]
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

write.csv(umap_liger_prot,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/celltype_annotations.csv')




