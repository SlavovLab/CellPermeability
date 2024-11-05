## For getting file path
# install.packages("rstudioapi")
library(rstudioapi)

# Source all the functions and packages needed
source(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),'/func.R'))

# you need to install QuantQC as well 
devtools::install_github("https://github.com/Andrew-Leduc/QuantQC")
library(QuantQC)




## -----------------------------------------------
#pSCoPE Analysis


#searched SC data, find on massive repo
data_path_fresh <- 'https://drive.google.com/uc?export=download&id=181SaMsT8ZjaqsAUUnKOkIFMsqHa8bbVu'
data_path_frozen <- 'https://drive.google.com/uc?export=download&id=1n3QphpaqomX0v6W0q7CJBSbvb6WtG32s'

# link raw file name to well plate, find in the github repo /QuantQC/AnalysisFromPaper/pSCoPE/
linker_fresh <- "https://drive.google.com/uc?export=download&id=1cx6iu0zJMGr8leS5k3Mhcew3BoVhCFnd"
linker_frozen <- "https://drive.google.com/uc?export=download&id=1n_LKjkyfT2_GqQZRbTafx1p7MLW-qxHH"

# Read in cell isolation files from CellenONE and assign cell type, find in folder /QuantQC/AnalysisFromPaper/pSCoPE/
fresh <-"https://drive.google.com/uc?export=download&id=1e1uDfu3d6YeSBlg3wXOmvOWdHBdZyEiH"
Fresh_cells <- list(Fresh = fresh)

frozen_sort1 <- "https://drive.google.com/uc?export=download&id=1DI4Il0_PvkV3CvGDhq9vSnwkR_3AIR_n"
frozen_sort2 <- "https://drive.google.com/uc?export=download&id=1mLJBRouNYqUxvGwmSXuDFChS5DXciCLt"
Frozen_cells <- list(Sort1 = frozen_sort1,
                     Sort2 = frozen_sort2)



# Generate the HTML Report - Fresh 
# Input desired out path
output_path_fresh <- 'report_fresh.html'

Gen_QQC_report_DDA(data_path = data_path_fresh,
                   linker_path = linker_fresh,
                   isolation = Fresh_cells,
                   output_path = output_path_fresh,
                   CV_thresh = .42,
                   plex = 29)


# Generate the HTML Report - Frozen
# Input desired out path
output_path_frozen <- 'report_frozen.html'
Gen_QQC_report_DDA(data_path = data_path_frozen,
                   linker_path = linker_frozen,
                   isolation = Frozen_cells,
                   output_path = output_path_frozen,
                   CV_thresh = .45,
                   plex = 29)



######------------------------------------------------------------------------
# Processing of single cell data, Ran on both preps to generate protein X cell matricies for downstream analysis


#------------------------------FRESH------------------------------------------

#Generate nPOP object from raw data
ProtLeakage <- MQ_to_QQC(data_path_fresh,linker_fresh, plex = 29,PIF_in = .7, PEP_in = 1)


# Normalize single cell runs to reference channel,
# filter out data points over twice reference
# Generate turn long table format to Peptide X Cell matrix

ProtLeakage <- TMT_Reference_channel_norm(ProtLeakage)


## Mapping cellenONE meta data to raw data
ProtLeakage <- link_cellenONE_Raw(ProtLeakage,Fresh_cells)

#plot exp design on glass slide
PlotSlideLayout_celltype(ProtLeakage)

PlotSlideLayout_label(ProtLeakage)




## Calculate statistics to look for decrease in LC/MS performance
# over time
ProtLeakage <- Calculate_run_order_statistics(ProtLeakage)

# Plot MS Intensity Statistics
PlotIntensityDrift(ProtLeakage)

#Plot Retention time Statistics
PlotRTDrift(ProtLeakage)


# Test negative controls, i.e. samples with no cell
ProtLeakage <- EvaluateNegativeControls(ProtLeakage)



# Make the classic neg ctrl plots
PlotNegCtrl(ProtLeakage,CV_thresh = .45)


# filter bad cells based off above, put in log10 intensity
ProtLeakage <- FilterBadCells(ProtLeakage, CV_thresh = .45)



# Compute the size relative to the carrier of each cell
PlotSCtoCarrierRatio(ProtLeakage)


# Plot cell size vs intensity in MS, options to color code by" "Run order" or "sample"
PlotCellSizeVsIntensity(ProtLeakage, type = "sample")


# Collapse and give statistics options:
#1 = median relative peptide, 2 = maxLFQ... soon to be directLFQwhich is much faster

ProtLeakage <- CollapseToProtein(ProtLeakage, 1)



# plot protein numbers
PlotProtAndPep(ProtLeakage)



PlotDataComplete(ProtLeakage)



# Compute correlations between peptides mapping to prot
ProtLeakage <- SharedPeptideCor(ProtLeakage)

PlotPepCor(ProtLeakage)


# Impute, still old imputation very slow, soon will update with C++ implementation
ProtLeakage <- KNN_impute(ProtLeakage)


## currently does label and LC/MS run but options comming soon
ProtLeakage <- BatchCorrect(ProtLeakage,run = T,labels = T)


ProtLeakage <- ComputePCA(ProtLeakage,imputed = F)

## plot PCA options are "Run order" "Total protein" "Condition" "Label"
PlotPCA(ProtLeakage, by = "Run order")

## also assigns louvain clusters
ProtLeakage <- ComputeUMAP(ProtLeakage)

## plots by cluster
PlotUMAP(ProtLeakage)
PlotUMAP(ProtLeakage, by = 'Condition')
PlotUMAP(ProtLeakage, by = 'Total protein')
PlotUMAP(ProtLeakage, by = 'Label')
PlotUMAP(ProtLeakage, by = 'Run order')


FeatureUMAP(ProtLeakage, prot = 'P24472', imputed = F)
FeatureUMAP(ProtLeakage, prot = 'P00329', imputed = F)
FeatureUMAP(ProtLeakage, prot = 'P56395', imputed = F)

FeatureUMAP(ProtLeakage, prot = 'P20152', imputed = F) #Fib
FeatureUMAP(ProtLeakage, prot = 'Q61233', imputed = F) #immune

FeatureUMAP(ProtLeakage, prot = 'P33267', imputed = F)



test2 <- ProtLeakage@reductions$UMAP
test2$ID <- rownames(test2)

testt <- ProtLeakage@meta.data %>% dplyr::select(ID,Intensity,Stain_Diameter)

test <- test %>% left_join(testt, by = c('ID'))
test$dead <- 'no'
test$dead[test$Stain_Diameter > 10] <- 'yes'



ggplot(test, aes(x = umap_1, y = umap_2, color = dead))+ geom_point() + dot_plot + theme_classic()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())

ggplot(test, aes(x = umap_1, y = umap_2, color = cluster))+ geom_point() + dot_plot + theme_classic()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())


write.csv(ProtLeakage@matricies@protein,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/protein.csv')
write.csv(ProtLeakage@matricies@protein.imputed,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/protein.imputed.csv')
write.csv(ProtLeakage@meta.data,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/meta.csv')


#------------------------------Frozen------------------------------------------

#Generate nPOP object from raw data
ProtLeakage <- MQ_to_QQC(data_path_frozen,linker_frozen, plex = 29,PIF_in = .7, PEP_in = 1)


# Normalize single cell runs to reference channel,
# filter out data points over twice reference
# Generate turn long table format to Peptide X Cell matrix

ProtLeakage <- TMT_Reference_channel_norm(ProtLeakage)


## Mapping cellenONE meta data to raw data
ProtLeakage <- link_cellenONE_Raw(ProtLeakage,Frozen_cells)

#plot exp design on glass slide
PlotSlideLayout_celltype(ProtLeakage)

PlotSlideLayout_label(ProtLeakage)




## Calculate statistics to look for decrease in LC/MS performance
# over time
ProtLeakage <- Calculate_run_order_statistics(ProtLeakage)

# Plot MS Intensity Statistics
PlotIntensityDrift(ProtLeakage)

#Plot Retention time Statistics
PlotRTDrift(ProtLeakage)


# Test negative controls, i.e. samples with no cell
ProtLeakage <- EvaluateNegativeControls(ProtLeakage)



# Make the classic neg ctrl plots
PlotNegCtrl(ProtLeakage,CV_thresh = .45)


# filter bad cells based off above, put in log10 intensity
ProtLeakage <- FilterBadCells(ProtLeakage, CV_thresh = .45)



# Compute the size relative to the carrier of each cell
PlotSCtoCarrierRatio(ProtLeakage)


# Plot cell size vs intensity in MS, options to color code by" "Run order" or "sample"
PlotCellSizeVsIntensity(ProtLeakage, type = "sample")


# Collapse and give statistics options:
#1 = median relative peptide, 2 = maxLFQ... soon to be directLFQwhich is much faster

ProtLeakage <- CollapseToProtein(ProtLeakage, 1)



# plot protein numbers
PlotProtAndPep(ProtLeakage)



PlotDataComplete(ProtLeakage)



# Compute correlations between peptides mapping to prot
ProtLeakage <- SharedPeptideCor(ProtLeakage)

PlotPepCor(ProtLeakage)


# Impute, still old imputation very slow, soon will update with C++ implementation
ProtLeakage <- KNN_impute(ProtLeakage)


## currently does label and LC/MS run but options comming soon
ProtLeakage <- BatchCorrect(ProtLeakage,run = T,labels = T)


ProtLeakage <- ComputePCA(ProtLeakage,imputed = F)

## plot PCA options are "Run order" "Total protein" "Condition" "Label"
PlotPCA(ProtLeakage, by = "Run order")

## also assigns louvain clusters
ProtLeakage <- ComputeUMAP(ProtLeakage)

## plots by cluster
PlotUMAP(ProtLeakage)
PlotUMAP(ProtLeakage, by = 'Condition')
PlotUMAP(ProtLeakage, by = 'Total protein')
PlotUMAP(ProtLeakage, by = 'Label')
PlotUMAP(ProtLeakage, by = 'Run order')


FeatureUMAP(ProtLeakage, prot = 'P24472', imputed = F)
FeatureUMAP(ProtLeakage, prot = 'P00329', imputed = F)
FeatureUMAP(ProtLeakage, prot = 'P56395', imputed = F)

FeatureUMAP(ProtLeakage, prot = 'P20152', imputed = F) #Fib
FeatureUMAP(ProtLeakage, prot = 'Q61233', imputed = F) #immune

FeatureUMAP(ProtLeakage, prot = 'P33267', imputed = F)



test2 <- ProtLeakage@reductions$UMAP
test2$ID <- rownames(test2)

testt <- ProtLeakage@meta.data %>% dplyr::select(ID,Intensity,Stain_Diameter)

test <- test %>% left_join(testt, by = c('ID'))
test$dead <- 'no'
test$dead[test$Stain_Diameter > 10] <- 'yes'



ggplot(test, aes(x = umap_1, y = umap_2, color = dead))+ geom_point() + dot_plot + theme_classic()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())

ggplot(test, aes(x = umap_1, y = umap_2, color = cluster))+ geom_point() + dot_plot + theme_classic()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())


write.csv(ProtLeakage@matricies@protein,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/protein.csv')
write.csv(ProtLeakage@matricies@protein.imputed,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/protein.imputed.csv')
write.csv(ProtLeakage@meta.data,'/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/meta.csv')


