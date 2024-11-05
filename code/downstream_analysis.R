## For getting file path
# install.packages("rstudioapi")
library(rstudioapi)

# Source all the functions and packages needed
source(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),'/func.R'))




# Read data
Frozen_anno <- read.csv('https://drive.google.com/uc?export=download&id=18YD1Cc9fz-UKlQwZWKhcdpX5ASrcUG5L')
Frozen_data <- read.csv('https://drive.google.com/uc?export=download&id=1uuSxzHfgPcNNYjpaAZIVPQ48i5R-6Gmq',row.names = 1,header = T)

Fresh_anno <- read.csv('https://drive.google.com/uc?export=download&id=1XEn5JY8Vmd_mSQmHXMUt0pe6_3wtGadm')
Fresh_data <- read.csv('https://drive.google.com/uc?export=download&id=1-f-gtu8cdC4OuWphLIxpJzLQ8AWDILU2',row.names = 1,header = T)







# Figure 2 - Cell type live/dead plot

Frozen_anno_lim <- Frozen_anno %>% select(dead,Cell_type)
Frozen_anno_lim$Cond <- 'Frozen'

Fresh_anno$dead <- 'no'
Fresh_anno$dead[Fresh_anno$Stain_Diameter > 10] <- 'yes'

Fresh_anno_lim <- Fresh_anno %>% select(dead,Cell_type)
Fresh_anno_lim$Cond <- 'Fresh'


plot_ct <- rbind(Fresh_anno_lim,Frozen_anno_lim)

plot_ct <- plot_ct %>% filter(Cell_type != 'NA')
plot_ct <- plot_ct %>% filter(Cell_type != 'Chondrocyte')
plot_ct <- plot_ct %>% filter(Cell_type != 'Ciliated')

ggplot(plot_ct, aes(x = Cell_type,fill = dead)) + geom_bar(position = 'dodge') +
  facet_wrap(~Cond,nrow = 2) + dot_plot + ylab('# single cells') + xlab('')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))

ggplot(plot_ct, aes(x = Cell_type)) + geom_bar(position = 'dodge') +
  dot_plot + ylab('# single cells') + xlab('')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))

Fresh_anno_lim <- Fresh_anno_lim %>% group_by(Cell_type) %>% summarise(fract = sum(dead =='yes',na.rm = T)/sum(is.na(dead) ==F))
Fresh_anno_lim <- Fresh_anno_lim %>% filter(is.na(Cell_type) == F)
Frozen_anno_lim <- Frozen_anno_lim %>% group_by(Cell_type) %>% summarise(fract = sum(dead =='yes',na.rm = T)/sum(is.na(dead) ==F))
Frozen_anno_lim <- Frozen_anno_lim %>% filter(is.na(Cell_type) == F)
Fresh_anno_lim$Cond <- 'Fresh'
Frozen_anno_lim$Cond <- 'Frozen'

plot_dead <- rbind(Fresh_anno_lim,Frozen_anno_lim)
plot_dead <- plot_dead %>% filter(Cell_type != 'Ciliated')

ggplot(plot_dead, aes(x = Cell_type,y = fract,color = Cond,group = Cond )) + geom_line() + geom_point(size = 4)+
  dot_plot + ylab(' fraction single cells') + xlab('')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1)) +
  scale_color_manual(values = c('Orange','Red'))

plot_ct$Cell_type[plot_ct$Cell_type == 'Chondrocyte'] <- ' Chondrocyte'
plot_ct$Cell_type[plot_ct$Cell_type == 'Immune'] <- ' Immune'
plot_ct$Cell_type[plot_ct$Cell_type == 'Chondrocyte'] <- ' Chondrocyte'
plot_ct$Cell_type[plot_ct$Cell_type == 'Club'] <- 'zClub'


# Figure 3 - Leaking vs Intact protein abundance
Frozen_data <- Frozen_data[rowSums(is.na(Frozen_data)==F) > 500,]

club <- Frozen_anno %>% filter(Cell_type == 'Club')
alive_club <- club %>% filter(dead == 'no')
dead_club <- club %>% filter(dead == 'yes')

bas <- Frozen_anno %>% filter(Cell_type == 'Basal')
alive_bas <- bas %>% filter(dead == 'no')
dead_bas <- bas %>% filter(dead == 'yes')

fib <- Frozen_anno %>% filter(Cell_type == 'Fibroblast')
alive_fib <- fib %>% filter(dead == 'no')
dead_fib <- fib %>% filter(dead == 'yes')

imm <- Frozen_anno %>% filter(Cell_type == 'Immune')
alive_imm <- imm %>% filter(dead == 'no')
dead_imm <- imm %>% filter(dead == 'yes')

fib_prot <- rowMeans(Frozen_data[,dead_fib$ID],na.rm = T) - rowMeans(Frozen_data[,alive_fib$ID],na.rm = T)

club_prot <- rowMeans(Frozen_data[,dead_club$ID],na.rm = T) - rowMeans(Frozen_data[,alive_club$ID],na.rm = T)

imm_prot <- rowMeans(Frozen_data[,dead_imm$ID],na.rm = T) - rowMeans(Frozen_data[,alive_imm$ID],na.rm = T)

bas_prot <- rowMeans(Frozen_data[,dead_bas$ID],na.rm = T) - rowMeans(Frozen_data[,alive_bas$ID],na.rm = T)


plot(club_prot,fib_prot,xlim = c(-1.5,.5),ylim = c(-1,.5))
cor(club_prot,fib_prot)

plot_diff <- as.data.frame(club_prot)
plot_diff$fib_prot <- fib_prot

plot_diff$club_prot <- plot_diff$club_prot - median(club_prot$club_prot[club_prot$compartment == 'Mitochondria'])
plot_diff$fib_prot <- plot_diff$fib_prot - median(club_prot$club_prot[club_prot$compartment == 'Mitochondria'])


ggplot(plot_diff, aes(x = fib_prot, y = club_prot)) + geom_point(alpha = .4)+ dot_plot + ggtitle('log2(Permeable/Intact), Cor = 0.66')+
  xlab('Fibroblast') + ylab('Club')

ggplot(plot_diff, aes(x = fib_prot, y = club_prot)) + ggpointdensity::geom_pointdensity() + dot_plot + ggtitle('log2(Permeable/Intact), Cor = 0.66')+
  xlab('Fibroblast') + ylab('Club') + scale_color_viridis() +theme_classic()


mat <- cbind(fib_prot,club_prot,imm_prot,bas_prot)

my_palette <- colorRampPalette(c("white", "red"))(100)

mat_plot <- cor(mat,use = 'pairwise.complete.obs')
diag(mat_plot) <- NA
pheatmap::pheatmap(cor(mat,use = 'pairwise.complete.obs'),color = my_palette)
pheatmap::pheatmap(mat_plot,color = my_palette)


mat$fib_prot <- mat$fib_prot - median(club_prot$club_prot[club_prot$compartment == 'Mitochondria'],na.rm = T)

write.csv(mat,'/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/SuppData/2024_Leduc_CellDeath/SupplementalFile3/leakage_by_celltype.csv')


# Figure 4 - Leaking protein compartment enrichment


Skel <- read.delim('/Users/andrewleduc/Desktop/Projects/AliveDead/SubCellLoc/Cytoskel.tsv')
Mito <- read.delim('/Users/andrewleduc/Desktop/Projects/AliveDead/SubCellLoc/Mito.tsv')
Cytosol <- read.delim('/Users/andrewleduc/Desktop/Projects/AliveDead/SubCellLoc/Cytosol.tsv')
Membrane <- read.delim('/Users/andrewleduc/Desktop/Projects/AliveDead/SubCellLoc/Membrane.tsv')
Nuc <- read.delim('/Users/andrewleduc/Desktop/Projects/AliveDead/SubCellLoc/Nuc.tsv')


Cytosol_unique <- Cytosol %>% filter(!Entry %in% c(Mito$Entry,Skel$Entry,Membrane$Entry,Nuc$Entry))
Mito_unique <- Mito %>% filter(!Entry %in% c(Cytosol$Entry,Skel$Entry,Membrane$Entry,Nuc$Entry))
#Skel_unique <- Skel %>% filter(!Entry %in% c(Mito$Entry,Cytosol$Entry,Membrane$Entry,Nuc$Entry))
Membrane_unique <- Membrane %>% filter(!Entry %in% c(Mito$Entry,Skel$Entry,Cytosol$Entry,Nuc$Entry))
Nuc_unique <- Nuc %>% filter(!Entry %in% c(Mito$Entry,Skel$Entry,Membrane$Entry,Cytosol$Entry))



club_prot <- as.data.frame(club_prot)
club_prot$prot <- rownames(club_prot)

club_prot$compartment <- NA

club_prot$compartment[club_prot$prot %in% Cytosol_unique$Entry] <- '  Cytosol'
club_prot$compartment[club_prot$prot %in% Membrane_unique$Entry] <- 'Membrane'
club_prot$compartment[club_prot$prot %in% Nuc_unique$Entry] <- ' Nuclear'
club_prot$compartment[club_prot$prot %in% Mito_unique$Entry] <- 'Mitochondria'

club_prot <- club_prot %>% filter(is.na(compartment)==F)

ggplot(club_prot,aes(x = compartment,y = club_prot)) + geom_boxplot() + dot_plot+ xlab('')+ylab('log2(Permeable - Intact)')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1)) + coord_cartesian(ylim = c(-1.5,.25))

club_prot$club_prot <- club_prot$club_prot - median(club_prot$club_prot[club_prot$compartment == 'Mitochondria'],na.rm = T)



# Sup Fig 1b

Fresh_data <- read.csv('/Users/andrewleduc/Desktop/Projects/AliveDead/Fresh/protein.csv',row.names = 1,header = T)
Frozen_data <- read.csv('/Users/andrewleduc/Desktop/Projects/AliveDead/Frozen/protein.csv',row.names = 1,header = T)

Fresh_sec <- Fresh_anno %>% filter(Cell_type %in% c('Club','Basal'))
Frozen_sec <- Frozen_anno %>% filter(Cell_type %in% c('Club','Basal'))

Frozen_data <- Frozen_data[,Frozen_sec$ID]
Fresh_data<- Fresh_data[,Fresh_sec$ID]

Frozen_data <- Frozen_data[rowSums(is.na(Frozen_data)==F) > 50,]

sect <- intersect(rownames(Fresh_data),rownames(Frozen_data))

Fresh_data <- Fresh_data[sect,]
Frozen_data <- Frozen_data[sect,]


Find_dif <- GetCorVect(Fresh_data,Frozen_data)

hist(Find_dif$cors,30)

Find_dif_bad <- Find_dif %>% filter(cors > .5)
Find_dif_bad <- Find_dif_bad %>% filter(cors <.5)



cyt <- intersect(rownames(Frozen_data),Cytosol_unique$Entry)

a <- Heatmap(cor(t(Frozen_data[cyt,]),use = 'pairwise.complete.obs'))
Heatmap(cor(t(Fresh_data[cyt,]),use = 'pairwise.complete.obs'))

order_plot <- rownames(Frozen_data[cyt,])[row_order(a)]

Heatmap(cor(t(Fresh_data[order_plot,]),use = 'pairwise.complete.obs'),cluster_rows = F,cluster_columns = F)
Heatmap(cor(t(Frozen_data[order_plot,]),use = 'pairwise.complete.obs'),cluster_rows = F,cluster_columns = F)


b <- Heatmap(cor(t(Frozen_data[Find_dif_bad$genes,]),use = 'pairwise.complete.obs'))
order_plot2 <- rownames(Frozen_data[Find_dif_bad$genes,])[row_order(b)]
Heatmap(cor(t(Fresh_data[order_plot2,]),use = 'pairwise.complete.obs'),cluster_rows = F,cluster_columns = F)
Heatmap(cor(t(Frozen_data[order_plot2,]),use = 'pairwise.complete.obs'),cluster_rows = F,cluster_columns = F)




### Sup figure 2 plots

# continuation of code from https://github.com/Andrew-Leduc/QuantQC/tree/main/AnalysisFromPaper/Analysis.R for pSCoPE data


dead <- AppNote@reductions$UMAP %>% filter(cluster %in% c(5,6,7))

alive <- AppNote@reductions$UMAP %>% filter(!cluster %in% c(5,6,7))

dead <- rowMeans(AppNote@matricies@protein[,rownames(alive)],na.rm = T) - rowMeans(AppNote@matricies@protein[,rownames(dead)],na.rm = T)


dead <- as.data.frame(dead)

FeatureUMAP(AppNote, prot = 'P20290') + theme_classic()

dead$prot <- rownames(dead)


Hum <- Proc_fasta('/Users/andrewleduc/Desktop/Github/QuantQC/inst/extdata/Human.fasta')
Mouse <- Proc_fasta('/Users/andrewleduc/Desktop/Github/QuantQC/inst/extdata/Mouse.fasta')
Mouse$split_gene <- toupper(Mouse$split_gene)


Mouse <- Mouse %>% filter(split_gene %in% Hum$split_gene)
Hum <- Hum %>% filter(split_gene %in% Mouse$split_gene)

Mouse <- Mouse[order(Mouse$split_gene),]
Hum <- Hum[order(Hum$split_gene),]

convert <- Mouse %>% left_join(Hum, by = c('split_gene'))
convert <- convert %>% filter(split_prot.y %in% dead$prot)
convert <- convert %>% filter(split_prot.x %in% mat_comp$prot)
convert <- convert %>% distinct(split_prot.y,.keep_all = T)

mat_comp <- mat %>% filter(prot %in% convert$split_prot.x)
dead <- dead %>% filter(prot %in% convert$split_prot.y)

dead <- dead[convert$split_prot.y,]
mat_comp <- mat_comp[convert$split_prot.x,]

dead$dead <- -dead$dead - median(dead$dead)

cor(mat_comp$bas_prot,dead$dead -1)

both_spec <- cbind(mat_comp$bas_prot,dead$dead -1)

colnames(both_spec) <- c('Mouse_Primary_Cells','Human_CellLines')
both_spec <- as.data.frame(both_spec)

ggplot(both_spec, aes(x = Mouse_Primary_Cells,y = Human_CellLines)) + geom_point() + dot_plot +
  ggtitle('Cor = 0.50')







####
# test <- test %>% left_join(Frozen_anno, by = ('ID'))
# test$Cell_type[test$Cell_type == 'Ciliated'] <- 'Club'
# test$Cell_type[test$cluster == '1'] <- 'Fibroblast'
# test$Cell_type[test$cluster == '6'] <- 'Chondrocyte'
# test <- test %>% filter(is.na(Cell_type)== F)
# test <- test %>% filter(is.na(dead)== F)
# 
# ggplot(test,aes(x = umap_1,y = umap_2,color = Cell_type))+geom_point() + theme_classic()
# ggplot(test,aes(x = umap_1,y = umap_2,color = dead))+geom_point() + theme_classic() +
#   scale_color_manual(values = c('gray60', 'black'))
# 
# test2 <- test2 %>% left_join(Fresh_anno, by = ('ID'))
# test2$Cell_type[test2$Cell_type == 'Ciliated'] <- 'Club'
# test2$Cell_type[test2$cluster == '2'] <- 'Club'
# test2$Cell_type[test2$cluster == '6'] <- 'Chondrocyte'
# test2 <- test2 %>% filter(is.na(Cell_type)== F)
# test2 <- test2 %>% filter(is.na(dead)== F)
# 
# ggplot(test2,aes(x = umap_1,y = umap_2,color = Cell_type))+geom_point() + theme_classic()
# 
# test2$dead <- 'no'
# test2$dead[test2$Stain_Diameter >10] <- 'yes'
# 
# test2$Cell_type[test2$cluster == '3'] <- 'Club'
# test2$Cell_type[test2$cluster == '6'] <- 'Chondrocyte'
# 
# ggplot(test2,aes(x = umap_1,y = umap_2,color = dead))+geom_point() + theme_classic()+
#   scale_color_manual(values = c('gray60', 'black'))
# 
# 





