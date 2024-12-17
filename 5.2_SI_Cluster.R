################################################################################
### Author: Nguyen Trung Duc
### Institute:ICAR- Indian Agricultural Research Institute, New Delhi, India
###           Nanaji Deshmukh Plant Phenomics Centre, ICAR-IARI,New Delhi, India         
###           Vietnam National University of Agriculture, Hanoi, Vietnam
### Email: ntduc11@gmail.com
################################################################################
#===============================================================================
#                       Paper: GWAS - HTP - NUE - RICE
#===============================================================================
###                           NUE  Donors selection
#===============================================================================

#Load Library-------------------------------------------------------------------
library(easypackages)
libraries("FactoMineR", "factoextra","corrplot", "ade4","RColorBrewer", "pheatmap",
          "ComplexHeatmap", "heatmaply", "metan")
#####Hierarchical clustering (for diversity analysis)
library("magrittr")
library("ape")
# Helper packages
library("dplyr")       # for data manipulation
library("tidyverse") # data manipulation
library("ggplot2")     # for data visualization
library("dendextend") # for comparing two dendrograms
# Modeling packages
library("cluster")     # for general clustering algorithms
library("factoextra")  # for visualizing cluster results
#Set working directory
setwd("F:/DUC_Code/HTP_NUE/04_SI-MGIDI/Data")

# Remove previous work
rm(list=ls())

# Load data
dat.control <- read.csv("SI_Control_mean.csv", header = TRUE, row.names = 1)
dat.Nstress <- read.csv("SI_Nstress_mean.csv", header = TRUE, row.names = 1)

#Compute optimum k cluster------------------------------------------------------
#Control treatment ---------------------
# Elbow method
fviz_nbclust(dat.control, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(dat.control, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(dat.control, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

#Nstress treatment ----------------------
# Elbow method
fviz_nbclust(dat.Nstress, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(dat.Nstress, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(dat.Nstress, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")
       
       
# Visualize K-means clustering------------------------------------------------------------
km.res.control <- eclust(dat.control, "kmeans", k = 4, nstart = 25, graph = FALSE)
km.res.Nstress <- eclust(dat.Nstress, "kmeans", k = 4, nstart = 25, graph = FALSE)

# Visualize k-means clusters Method 1 ------------------------------------------
fviz_cluster(km.res.control,  frame.type = "norm", frame.level = 0.68)+ theme_bw() # the default
fviz_cluster(km.res.Nstress,  frame.type = "norm", frame.level = 0.68)+ theme_bw() # the default

# Visualize k-means clusters Method 2 ------------------------------------------
fviz_cluster(km.res.control, ellipse.type = "euclid", star.plot = TRUE, repel = TRUE, ggtheme = theme_bw() )
fviz_cluster(km.res.Nstress, ellipse.type = "euclid", star.plot = TRUE, repel = TRUE, ggtheme = theme_bw() )

