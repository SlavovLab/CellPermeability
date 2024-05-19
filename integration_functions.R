library(rliger)
library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(openxlsx)
library(HGNChelper)
library(seqinr)
library(ComplexHeatmap)
library(Hmisc)

TLS <- function(vect1,vect2){

  vect1[vect1 == -Inf] <- NA
  vect1[vect1 == Inf] <- NA



  vect2[vect2 == -Inf] <- NA
  vect2[vect2 == Inf] <- NA

  int_x <- mean(vect1,na.rm=T)
  int_y <- mean(vect2,na.rm=T)

  vect1 <- vect1-int_x
  vect2 <- vect2-int_y

  mat <- cbind(vect1,vect2)

  mat <- mat[complete.cases(mat),]

  TLS_mat <- svd(mat)$v

  slope <- TLS_mat[1,1]/TLS_mat[2,1]

  int <- c(int_x,int_y)

  return(list(slope,int))

}


Filter_missingness <- function(data,cluster){

}

GetCorVect <- function(prot,mRNA){
  prot_cor <- rcorr(t(prot))$r
  rna_cor <- rcorr(t(mRNA))$r

  diag(prot_cor) <- NA
  diag(rna_cor) <- NA

  genes <- c()
  cors <- c()
  sd_diff <- c()

  for(i in 1:ncol(rna_cor)){
    genes <- c(genes,rownames(rna_cor)[i])
    cors <- c(cors,cor(prot_cor[,i],rna_cor[,i], use = 'pairwise.complete.obs'))
    sd_diff <- c(sd_diff,(sd(rna_cor[,i],na.rm = T)-sd(prot_cor,na.rm = T)))
  }

  df <- as.data.frame(cors)
  df$genes <- genes

  return(df)


}


Integrate_liger <- function(mRNA_counts,prot_dat,k_in){

  sc.batch_cor_p <- 2^prot_dat

  liger_epi <- rliger::createLiger(list(mRNA = mRNA_counts, prot = sc.batch_cor_p))

  liger_epi <- rliger::normalize(liger_epi)

  liger_epi <- selectGenes(liger_epi,var.thresh = .01)
  liger_epi@var.genes  <- rownames(mRNA_counts) #intersect(vgs$split_gene,rownames(RNA_Seq_dat))

  liger_epi <- scaleNotCenter(liger_epi)

  liger_epi@scale.data$prot <- t(sc.batch_cor_p[liger_epi@var.genes,])

  liger_epi <- optimizeALS(liger_epi, k=k_in) #, use.unshared = TRUE)
  liger_epi <- quantile_norm(liger_epi, ref_dataset= "mRNA")
  liger_epi <- louvainCluster(liger_epi,resolution = .3)
  liger_epi <- rliger::runUMAP(liger_epi)

  umap_plots <-plotByDatasetAndCluster(liger_epi, axis.labels = c("UMAP1","UMAP2"), return.plots = TRUE)


  calcAlignment(liger_epi)
  umap_liger <- as.data.frame(umap_plots[[1]][["data"]])
  return(list(calcAlignment(liger_epi),umap_liger))


}


filt.mat.cr<-function(mat, pct.r,pct.c){
  mat <- as.matrix(mat)
  kc<-c()
  for(k in 1:ncol(mat)){

    pct.na<-length(which(is.na(mat[,k]))) / length(mat[,k])
    if(pct.na <= pct.c){ kc<-c(kc,k)}
    #print(pct.na)


  }

  mat<-mat[,kc]

  kr<-c()
  for(k in 1:nrow(mat)){

    pct.na<-length(which(is.na(mat[k,]))) / length(mat[k,])
    if(pct.na <= pct.r){ kr<-c(kr,k)}
    #print(pct.na)


  }

  mat<-mat[kr,]



  return(mat)

}


Auto_Annotate <- function(mRNA_Seurat){
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  tissue = c("Lung") # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus


  # prepare gene sets
  gs_list = gene_sets_prepare(db_, tissue)
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = mRNA_Seurat[["RNA"]]@scale.data, scaled = TRUE,
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

  # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix.
  # In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
  # or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(mRNA_Seurat@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(mRNA_Seurat@meta.data[mRNA_Seurat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(mRNA_Seurat@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])

  mRNA_Seurat@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,];
    mRNA_Seurat@meta.data$customclassif[mRNA_Seurat@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }


  return(mRNA_Seurat)

}


Proc_fasta <- function(path){
  convert_mouse <- read.fasta(path,set.attributes = T,whole.header = T)
  convert_mouse <- names(convert_mouse)
  parse_row<-grep("GN=",convert_mouse, fixed=T)
  split_prot<-str_split(convert_mouse[parse_row], pattern = fixed("GN="))
  gene<-unlist(split_prot)[seq(2,2*length(split_prot),2)]
  prot <- unlist(split_prot)[seq(1,2*length(split_prot),2)]
  prot_parse <- grep("|",prot, fixed=T)
  gene_parse <- grep(" ",gene, fixed=T)
  split_gene<-str_split(gene[parse_row], pattern = fixed(" "))
  split_gene<-unlist(split_gene)[seq(1,3*length(split_gene),3)]
  split_prot<-str_split(prot[parse_row], pattern = fixed("|"))
  split_prot<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  convert_mouse  <- as.data.frame(cbind(split_prot,split_gene))

  return(convert_mouse)
}


read.prot_gene <- function(data.path,fasta_conv){
  fasta_conv <- convert_mouse

  data <- read.csv(data.path,row.names = 1)

  fasta_conv <- fasta_conv %>% distinct(split_gene,.keep_all = T)

  data <- data %>% dplyr::filter(rownames(data) %in% fasta_conv$split_prot)
  fasta_conv <- fasta_conv %>% filter(split_prot %in% rownames(data))



  data <- data[fasta_conv$split_prot,]

  rownames(data) <- fasta_conv$split_gene
  data <- as.matrix(data)
  return(data)

}
