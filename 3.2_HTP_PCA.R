################################################################################
### Author: Nguyen Trung Duc
### Email: ntduc11@gmail.com
################################################################################
#===============================================================================
###                            DATA VISUALIZATION
#===============================================================================
# Remove previous work
rm(list=ls())
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
#1. Load library
library("easypackages") #Load multiple libraries
libraries("gplots", "ggplot2", "dplyr","lavaan", "plyr", "cowplot", "rmarkdown", "rmarkdown",
          "readr", "caTools", "bitops","Rmisc", "PupillometryR", "FactoMineR",
          "factoextra","corrplot","ade4","ComplexHeatmap", "ggrepel","dabestr")
libraries( "rstatix", "ggpubr","stringr", "ggplot2", "tidyverse",
           "data.table", "readr","plotly", "DT", "pheatmap", "VennDiagram","patchwork",
           "heatmaply", "ggcorrplot", "RColorBrewer", "hrbrthemes", "tm", "proustr")
libraries("plyr", "dplyr", "reshape2", "readxl","stringr", "lubridate")                       #Data transformation
libraries("FactoMineR", "factoextra","corrplot","ade4","ComplexHeatmap", "ggrepel","dabestr") #data-visualization 
libraries("xlsx", "rJava", "fmsb","rlang", "car")      #Read data file
libraries("VIM", "outliers", "mice")                   #Outliers detection and management
libraries("asremlPlus", "nlme","lme4", "FW", "sommer") #Linear mixed models
libraries("metan", "SpATS","statgenGxE", "statgenSTA") #GxExM analysis - Multi-traits Selection

#2. Set working directory
#setwd("F:/DUC_Code/HTP_NUE")

#===============================================================================
###                        Correlation heatmap
#===============================================================================
data1 <-read.csv("Phenotypes/68DAS.csv")
data2 <-read.csv("Phenotypes/75DAS.csv")
data3 <-read.csv("Phenotypes/83DAS.csv")
data4 <-read.csv("Phenotypes/manual.csv")

#Convert variables into appropriate data types
data1$Group     <-as.factor(data1$Group)      # Group (ind, aus, japonica, basmati, admix) as factor
data1$Genotype  <-as.factor(data1$Genotype)  # Genotype as factor
data1$Treatment <-as.factor(data1$Treatment)  # Treatment as factor
data1$PlantID   <-as.factor(data1$PlantID)    # PlantID as factor
data1$Rep       <-as.factor(data1$Rep)        # Replication as factor

#Convert variables into appropriate data types
data2$Group     <-as.factor(data2$Group)      # Group (ind, aus, japonica, basmati, admix) as factor
data2$Genotype  <-as.factor(data2$Genotype)  # Genotype as factor
data2$Treatment <-as.factor(data2$Treatment)  # Treatment as factor
data2$PlantID   <-as.factor(data2$PlantID)    # PlantID as factor
data2$Rep       <-as.factor(data2$Rep)        # Replication as factor

#Convert variables into appropriate data types
data3$Group     <-as.factor(data3$Group)      # Group (ind, aus, japonica, basmati, admix) as factor
data3$Genotype  <-as.factor(data3$Genotype)  # Genotype as factor
data3$Treatment <-as.factor(data3$Treatment)  # Treatment as factor
data3$PlantID   <-as.factor(data3$PlantID)    # PlantID as factor
data3$Rep       <-as.factor(data3$Rep)        # Replication as factor

#Convert variables into appropriate data types
data4$Group     <-as.factor(data4$Group)      # Group (ind, aus, japonica, basmati, admix) as factor
data4$Genotype  <-as.factor(data4$Genotype)  # Genotype as factor
data4$Treatment <-as.factor(data4$Treatment)  # Treatment as factor
data4$PlantID   <-as.factor(data4$PlantID)    # PlantID as factor
data4$Rep       <-as.factor(data4$Rep)        # Replication as factor
#===============================================================================
#           Heatmap correlation between i-traits and manual traits
#===============================================================================
#Define subset data
Control   <- subset(data1, Treatment == "Control")
Nstress   <- subset(data1, Treatment == "Nstress")
traits <- colnames(data1)[-c(1:5)]
C <- Control[c(traits)]
Ns <-Nstress[c(traits)]
### For Control treatment (NStress can do the same)
CorC <-cor(C, method = "pearson", use = "pairwise.complete.obs")
CorNs <-cor(Ns, method = "pearson", use = "pairwise.complete.obs")
### Heatmap---------------------------------------------------------------------
libraries("cluster","circlize", "ComplexHeatmap", "dendextend",
          "factoextra",  "NbClust")

#Computing the number of clusters-----------------------------------------------
# Elbow method
fviz_nbclust(CorC, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(CorC, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(CorC, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

#Define group
row_dend = as.dendrogram(hclust(dist(CorC)))
row_dend = color_branches(row_dend, k = 4)
col_dend = as.dendrogram(hclust(dist(CorC)))
col_dend = color_branches(col_dend, k = 4)
#Set annotation-----------------------------------------------------------------
set.seed(123)
Traits <-read.csv("Phenotypes/Traits.csv")
ann <- data.frame(Traits$Group)
colnames(ann) <- c('Group')
colors <- list('Group' = c('Manual' = '#2cef3b',
                  'Biomass-related' = '#ef3b2c',
                    'Architectural' = '#3b2cef', 
                    'Physiological' = '#ef9d2c'))
colAnn <-HeatmapAnnotation(df = ann,
                        which = 'col',
                          col = colors)

png("04.2_Heatmap/CorC_manual_83DAS.png", width=14, height=10.5, units="in", res=800)
hmap <-Heatmap(CorC,
               name = "Corr",
               show_row_names = TRUE,
               row_names_gp = gpar(fontsize = 8),
               show_column_names = TRUE,
               column_names_gp = gpar(fontsize = 9),
               cluster_rows = row_dend,
               cluster_columns = col_dend,
               show_column_dend = TRUE,
               show_row_dend = TRUE,
               row_dend_reorder = TRUE,
               column_dend_reorder = TRUE,
               top_annotation = colAnn)
draw(hmap,heatmap_legend_side = "left", annotation_legend_side= "right")
dev.off()

png("04.2_Heatmap/CorNs_manual_83DAS.png", width=14, height=10.5, units="in", res=800)
hmap <-Heatmap(CorNs,
               name = "Corr",
               show_row_names = TRUE,
               row_names_gp = gpar(fontsize = 8),
               show_column_names = TRUE,
               column_names_gp = gpar(fontsize = 9),
               cluster_rows = row_dend,
               cluster_columns = col_dend,
               show_column_dend = TRUE,
               show_row_dend = TRUE,
               row_dend_reorder = TRUE,
               column_dend_reorder = TRUE,
               top_annotation = colAnn)
draw(hmap,heatmap_legend_side = "left", annotation_legend_side= "right")
dev.off()


### Correlation manual traits---------------------------------------------------
#Define subset data
w = 13
h = 10

Control_manual   <- subset(data4, Treatment == "Control")
Nstress_manual   <- subset(data4, Treatment == "Nstress")
traits <- colnames(data4)[-c(1:5)]

C_manual <- Control_manual[c(traits)]
Ns_manual <-Nstress_manual[c(traits)]
### For Control treatment (NStress can do the same)
CorC_manual <-cor(C_manual, method = "pearson", use = "pairwise.complete.obs")
CorNs_manual <-cor(Ns_manual, method = "pearson", use = "pairwise.complete.obs")

png("04.2_Heatmap/Cor_Control_manual.png", width=10.5, height=10, units="in", res=800)
c1 <-corr_coef(Control_manual)
plot(c1,
     size.text.cor = 3,
     size.text.signif = 3,
     size.text.lab = 14,
     digits.cor = 2,
     digits.pval = 2,
     caption = TRUE,
     reorder = FALSE)
dev.off()

png("04.2_Heatmap/Cor_Nstress_manual.png", width=10.5, height=10, units="in", res=800)
c2 <-corr_coef(Nstress_manual)
plot(c2,
     size.text.cor = 3,
     size.text.signif = 3,
     size.text.lab = 14,
     digits.cor = 2,
     digits.pval = 2,
     caption = TRUE,
     reorder = FALSE)
dev.off()

#===============================================================================
#Principal component analysis
#===============================================================================
#==============================  68DAS   =======================================
#Define subset data
Control   <- subset(data1, Treatment == "Control")
Nstress   <- subset(data1, Treatment == "Nstress")
traits <- colnames(data1)[-c(1:5)]

C <- Control[c(traits)]
Ns <-Nstress[c(traits)]
### For Control treatment (NStress can do the same)
CorC <-cor(C, method = "pearson", use = "pairwise.complete.obs")
CorNs <-cor(Ns, method = "pearson", use = "pairwise.complete.obs")
#PCA analysis
###Control
res.pca <- PCA(CorC,  graph = FALSE)
#02. Extract and visualize eigenvalues/variances:
### Extract eigenvalues/variances
get_eig(res.pca)
### Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
#03.Extract and visualize results for variables:
# Extract the results for variables
var <- get_pca_var(res.pca)
# Coordinates of variables
head(var$coord)
# Contribution of variables
head(var$contrib)
# Control variable colors using their contributions
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)             # Avoid text overlapping
# Contributions of variables to PC1
w = 12
h = 4.5
PC1 <-fviz_contrib(res.pca, choice = "var", axes = 1, top = 120) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('04.1_Figures/PC1-68DAS_Control.png', width = w, height = h)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 120)

### Nstress
res.pca <- PCA(Ns,  graph = FALSE)
#02. Extract and visualize eigenvalues/variances:
### Extract eigenvalues/variances
get_eig(res.pca)
### Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
#03.Extract and visualize results for variables:
# Extract the results for variables
var <- get_pca_var(res.pca)
# Coordinates of variables
head(var$coord)
# Contribution of variables
head(var$contrib)
# Control variable colors using their contributions
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)             # Avoid text overlapping
# Contributions of variables to PC1
w = 12
h = 4.5
PC2 <-fviz_contrib(res.pca, choice = "var", axes = 1, top = 120) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('04.1_Figures/PC1-68DAS_Nstress.png', width = w, height = h)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 120)

#===============================================================================
#==============================  75DAS   =======================================
#Define subset data
Control   <- subset(data2, Treatment == "Control")
Nstress   <- subset(data2, Treatment == "Nstress")
traits <- colnames(data2)[-c(1:5)]

C <- Control[c(traits)]
Ns <-Nstress[c(traits)]

### For Control treatment (NStress can do the same)
CorC <-cor(C, method = "pearson", use = "pairwise.complete.obs")
CorNs <-cor(Ns, method = "pearson", use = "pairwise.complete.obs")

#===============================================================================
#PCA analysis
###Control
res.pca <- PCA(CorC,  graph = FALSE)
#02. Extract and visualize eigenvalues/variances:
### Extract eigenvalues/variances
get_eig(res.pca)
### Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
#03.Extract and visualize results for variables:
# Extract the results for variables
var <- get_pca_var(res.pca)
# Coordinates of variables
head(var$coord)
# Contribution of variables
head(var$contrib)
# Control variable colors using their contributions
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)             # Avoid text overlapping
# Contributions of variables to PC1
w = 12
h = 4.5
PC1 <-fviz_contrib(res.pca, choice = "var", axes = 1, top = 120) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('04.1_Figures/PC1-75DAS_Control.png', width = w, height = h)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 120)

### Nstress
res.pca <- PCA(Ns,  graph = FALSE)
#02. Extract and visualize eigenvalues/variances:
### Extract eigenvalues/variances
get_eig(res.pca)
### Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
#03.Extract and visualize results for variables:
# Extract the results for variables
var <- get_pca_var(res.pca)
# Coordinates of variables
head(var$coord)
# Contribution of variables
head(var$contrib)
# Control variable colors using their contributions
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)             # Avoid text overlapping
# Contributions of variables to PC1
w = 12
h = 4.5
PC2 <-fviz_contrib(res.pca, choice = "var", axes = 1, top = 120) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('04.1_Figures/PC1-75DAS_Nstress.png', width = w, height = h)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 120)

#===============================================================================
#==============================  83DAS   =======================================
#Define subset data
Control   <- subset(data3, Treatment == "Control")
Nstress   <- subset(data3, Treatment == "Nstress")
traits <- colnames(data3)[-c(1:5)]

C <- Control[c(traits)]
Ns <-Nstress[c(traits)]
### For Control treatment (NStress can do the same)
CorC <-cor(C, method = "pearson", use = "pairwise.complete.obs")
CorNs <-cor(Ns, method = "pearson", use = "pairwise.complete.obs")

#PCA analysis
###Control
res.pca <- PCA(CorC,  graph = FALSE)
#02. Extract and visualize eigenvalues/variances:
### Extract eigenvalues/variances
get_eig(res.pca)
### Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
#03.Extract and visualize results for variables:
# Extract the results for variables
var <- get_pca_var(res.pca)
# Coordinates of variables
head(var$coord)
# Contribution of variables
head(var$contrib)
# Control variable colors using their contributions
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)             # Avoid text overlapping
# Contributions of variables to PC1
w = 12
h = 4.5
PC1 <-fviz_contrib(res.pca, choice = "var", axes = 1, top = 120) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('04.1_Figures/PC1-83DAS_Control.png', width = w, height = h)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 120)

### Nstress
res.pca <- PCA(Ns,  graph = FALSE)
#02. Extract and visualize eigenvalues/variances:
### Extract eigenvalues/variances
get_eig(res.pca)
### Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
#03.Extract and visualize results for variables:
# Extract the results for variables
var <- get_pca_var(res.pca)
# Coordinates of variables
head(var$coord)
# Contribution of variables
head(var$contrib)
# Control variable colors using their contributions
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)             # Avoid text overlapping
# Contributions of variables to PC1
w = 12
h = 4.5
PC2 <-fviz_contrib(res.pca, choice = "var", axes = 1, top = 120) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('04.1_Figures/PC1-83DAS_Nstress.png', width = w, height = h)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 120)

#===============================================================================
#==============================  Manual   ======================================
#Define subset data
Control   <- subset(data4, Treatment == "Control")
Nstress   <- subset(data4, Treatment == "Nstress")
traits <- colnames(data4)[-c(1:5)]

C <- Control[c(traits)]
Ns <-Nstress[c(traits)]
### For Control treatment (NStress can do the same)
CorC <-cor(C, method = "pearson", use = "pairwise.complete.obs")
CorNs <-cor(Ns, method = "pearson", use = "pairwise.complete.obs")

#PCA analysis
###Control
res.pca <- PCA(CorC,  graph = FALSE)
#02. Extract and visualize eigenvalues/variances:
### Extract eigenvalues/variances
get_eig(res.pca)
### Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
#03.Extract and visualize results for variables:
# Extract the results for variables
var <- get_pca_var(res.pca)
# Coordinates of variables
head(var$coord)
# Contribution of variables
head(var$contrib)
# Control variable colors using their contributions
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) +             # Avoid text overlapping
             theme_bw()
w = 8
h = 8
ggsave('Figures/PCA_manual_Control.png', width = w, height = h)

# Contributions of variables to PC1
w = 12
h = 4.5
PC1 <-fviz_contrib(res.pca, choice = "var", axes = 1, top = 120) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Figures/PC1-Manual_Control.png', width = w, height = h)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 120)

### Nstress
res.pca <- PCA(Ns,  graph = FALSE)
#02. Extract and visualize eigenvalues/variances:
### Extract eigenvalues/variances
get_eig(res.pca)
### Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
#03.Extract and visualize results for variables:
# Extract the results for variables
var <- get_pca_var(res.pca)
# Coordinates of variables
head(var$coord)
# Contribution of variables
head(var$contrib)
# Control variable colors using their contributions
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) +             # Avoid text overlapping
             theme_bw() 
w = 8
h = 8
ggsave('Figures/PCA_manual_Nstress.png', width = w, height = h)

# Contributions of variables to PC1
w = 12
h = 4.5
PC2 <-fviz_contrib(res.pca, choice = "var", axes = 1, top = 120) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Figures/PC1-Manual_Nstress.png', width = w, height = h)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 120)

