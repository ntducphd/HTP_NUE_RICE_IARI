################################################################################
### Author: Nguyen Trung Duc
### Institute:ICAR - Indian Agricultural Research Institute, New Delhi, India
###           Nanaji Deshmukh Plant Phenomics Centre, ICAR-IARI,New Delhi, India         
###           Vietnam National University of Agriculture, Hanoi, Vietnam
### Email: ntduc11@gmail.com
################################################################################
#===============================================================================
#                       Paper: GWAS - HTP - NUE - RICE
#===============================================================================
###                  Genetic and phenomics diversity analysis
#===============================================================================
#Clean the global environment
rm(list = ls())

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
setwd("E:/HTP_GWAS_NUE/Codes/01_Genetic_Diversity/Phenotypes")

# Remove previous work
rm(list=ls())

data.cor68 <- read.csv("68DAS_control.csv", header = TRUE, row.names = 1)
data.cor75 <- read.csv("75DAS_control.csv", header = TRUE, row.names = 1)
data.cor83 <- read.csv("83DAS_control.csv", header = TRUE, row.names = 1)
# K-means clustering
km.res68 <- eclust(data.cor68, "kmeans", k = 4, nstart = 25, graph = FALSE)
km.res75 <- eclust(data.cor75, "kmeans", k = 4, nstart = 25, graph = FALSE)
km.res83 <- eclust(data.cor83, "kmeans", k = 4, nstart = 25, graph = FALSE)

# Visualize k-means clusters method 1 ------------------------------------------
fviz_cluster(km.res68,  frame.type = "norm", frame.level = 0.68)+ theme_bw() # the default
fviz_cluster(km.res75,  frame.type = "norm", frame.level = 0.68)+ theme_bw() # the default
fviz_cluster(km.res83,  frame.type = "norm", frame.level = 0.68)+ theme_bw() # the default

# Visualize k-means clusters method 2 ------------------------------------------
fviz_cluster(km.res68, ellipse.type = "euclid", star.plot = TRUE, repel = TRUE, ggtheme = theme_bw() )
fviz_cluster(km.res75, ellipse.type = "euclid", star.plot = TRUE, repel = TRUE, ggtheme = theme_bw() )
fviz_cluster(km.res83, ellipse.type = "euclid", star.plot = TRUE, repel = TRUE, ggtheme = theme_bw() )

###Extract cluster result-------------------------------------------------------
km.res68$size
km.res75$size
km.res83$size

###Compare two clusters
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
# Compute distance matrix
res.dist.68 <- dist(data.cor68, method = "euclidean")
res.dist.75 <- dist(data.cor75, method = "euclidean")
# Compute hierarchical clustering
hc1 <- hclust(res.dist.68, "ward.D")
hc2 <- hclust(res.dist.75, "ward.D")
# Create two dendrograms
dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)
tanglegram(hc1, hc2)

#===============================================================================
#Compute optimum k cluster------------------------------------------------------
### 68DAS-----------------------------------------------------------------------
# Elbow method
fviz_nbclust(data.cor68, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(data.cor68, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(data.cor68, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

### 75DAS-----------------------------------------------------------------------
# Elbow method
fviz_nbclust(data.cor75, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(data.cor75, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(data.cor75, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

### 83DAS-----------------------------------------------------------------------
# Elbow method
fviz_nbclust(data.cor83, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(data.cor83, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot = 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(data.cor83, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

