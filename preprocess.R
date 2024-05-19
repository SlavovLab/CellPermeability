library(QuantQC)
library(ggpointdensity)
library(viridis)
library(ggpubr)
library(dplyr)
library(matrixStats)
library(reshape2)
library(stringr)



# link gene names and uniprot
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



## -----------------------------------------------
#pSCoPE Analysis


#searched SC data, find on massive repo
data_path_fresh <- "/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users/aleduc/AliveDead_project/plate2/evidence_p2.txt"
data_path_frozen <- "/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users/aleduc/AliveDead_project/plate1/evidence_p1.txt"

# link raw file name to well plate, find in the github repo /QuantQC/AnalysisFromPaper/pSCoPE/
linker_fresh <- "/Users/andrewleduc/Desktop/Projects/AliveDead/Fresh/linker_p2.csv"
linker_frozen <- "/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/linker_p1.csv"


# Read in cell isolation files from CellenONE and assign cell type, find in folder /QuantQC/AnalysisFromPaper/pSCoPE/
fresh <-"/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users/aleduc/AliveDead_project/plate2/Alive_isolated.xls"
Fresh_cells <- list(Fresh = fresh)

dmso <- "/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users/aleduc/AliveDead_project/plate1/dmso_deadttt.xls"
commer <- "/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users/aleduc/AliveDead_project/plate1/commerc_solated.xls"
Frozen_cells <- list(Commercial = commer,
                      DMSO = dmso)



# Generate the HTML Report - Fresh
output_path_fresh <- '/Users/andrewleduc/Desktop/Projects/AliveDead/report_fresh.html'

Gen_QQC_report_DDA(data_path = data_path_fresh,
                   linker_path = linker_fresh,
                   isolation = Fresh_cells,
                   output_path = output_path_fresh,
                   CV_thresh = .42,
                   plex = 29)


# Generate the HTML Report - Frozen
output_path_frozen <- '/Users/andrewleduc/Desktop/Projects/AliveDead/report_frozen.html'
Gen_QQC_report_DDA(data_path = data_path_frozen,
                   linker_path = linker_frozen,
                   isolation = Frozen_cells,
                   output_path = output_path_frozen,
                   CV_thresh = .45,
                   plex = 29)



######------------------------------------------------------------------------
# Processing of single cell data, Ran on both preps to generate protein X cell matricies for downstream analysis


#Generate nPOP object from raw data
AppNote <- MQ_to_QQC(data_path_fresh,linker_fresh, plex = 29,PIF_in = .7, PEP_in = 1)


# Normalize single cell runs to reference channel,
# filter out data points over twice reference
# Generate turn long table format to Peptide X Cell matrix

AppNote <- TMT_Reference_channel_norm(AppNote)


## Mapping cellenONE meta data to raw data
AppNote <- link_cellenONE_Raw(AppNote,Fresh_cells)

#plot exp design on glass slide
PlotSlideLayout_celltype(AppNote)

PlotSlideLayout_label(AppNote)




## Calculate statistics to look for decrease in LC/MS performance
# over time
AppNote <- Calculate_run_order_statistics(AppNote)

# Plot MS Intensity Statistics
PlotIntensityDrift(AppNote)

#Plot Retention time Statistics
PlotRTDrift(AppNote)


# Test negative controls, i.e. samples with no cell
AppNote <- EvaluateNegativeControls(AppNote)



# Make the classic neg ctrl plots
PlotNegCtrl(AppNote,CV_thresh = .45)


# filter bad cells based off above, put in log10 intensity
AppNote <- FilterBadCells(AppNote, CV_thresh = .45)



# Compute the size relative to the carrier of each cell
PlotSCtoCarrierRatio(AppNote)


# Plot cell size vs intensity in MS, options to color code by" "Run order" or "sample"
PlotCellSizeVsIntensity(AppNote, type = "sample")


# Collapse and give statistics options:
#1 = median relative peptide, 2 = maxLFQ... soon to be directLFQwhich is much faster

AppNote <- CollapseToProtein(AppNote, 1)



# plot protein numbers
PlotProtAndPep(AppNote)



PlotDataComplete(AppNote)



# Compute correlations between peptides mapping to prot
AppNote <- SharedPeptideCor(AppNote)

PlotPepCor(AppNote)


# Impute, still old imputation very slow, soon will update with C++ implementation
AppNote <- KNN_impute(AppNote)


## currently does label and LC/MS run but options comming soon
AppNote <- BatchCorrect(AppNote,run = T,labels = T)


AppNote <- ComputePCA(AppNote,imputed = F)

## plot PCA options are "Run order" "Total protein" "Condition" "Label"
PlotPCA(AppNote, by = "Run order")

## also assigns louvain clusters
AppNote <- ComputeUMAP(AppNote)

## plots by cluster
PlotUMAP(AppNote)
PlotUMAP(AppNote, by = 'Condition')
PlotUMAP(AppNote, by = 'Total protein')
PlotUMAP(AppNote, by = 'Label')
PlotUMAP(AppNote, by = 'Run order')


FeatureUMAP(AppNote, prot = 'P24472', imputed = F)
FeatureUMAP(AppNote, prot = 'P00329', imputed = F)
FeatureUMAP(AppNote, prot = 'P56395', imputed = F)

FeatureUMAP(AppNote, prot = 'P20152', imputed = F) #Fib
FeatureUMAP(AppNote, prot = 'Q61233', imputed = F) #immune

FeatureUMAP(AppNote, prot = 'P33267', imputed = F)

View(AppNote@pep.cor[[1]])

FeaturePlot(trach_gergana, features = 'Cyp2f2')


test2 <- AppNote@reductions$UMAP
test2$ID <- rownames(test2)

testt <- AppNote@meta.data %>% dplyr::select(ID,Intensity,Stain_Diameter)

test <- test %>% left_join(testt, by = c('ID'))
test$dead <- 'no'
test$dead[test$Stain_Diameter > 10] <- 'yes'

sum(test$dead == 'no')/nrow(test)

ggplot(test, aes(x = umap_1, y = umap_2, color = dead))+ geom_point() + dot_plot + theme_classic()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())

ggplot(test, aes(x = umap_1, y = umap_2, color = cluster))+ geom_point() + dot_plot + theme_classic()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())


write.csv(AppNote@matricies@protein,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/protein.csv')
write.csv(AppNote@matricies@protein.imputed,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/protein.imputed.csv')
write.csv(AppNote@meta.data,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/meta.csv')




