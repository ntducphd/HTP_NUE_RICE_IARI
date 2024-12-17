################################################################################
### Author: Nguyen Trung Duc
### Institute: ICAR- Indian Agricultural Research Institute, New Delhi, India
###            Nanaji Deshmukh Plant Phenomics Centre, ICAR-IARI, New Delhi, India 
###            Vietnam National University of Agriculture, Hanoi, Vietnam
### Email: ntduc11@gmail.com
################################################################################
#===============================================================================
###                            DATA VISUALIZATION
#===============================================================================
# Remove previous work
rm(list=ls())

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
libraries("metan", "SpATS","statgenGxE", "statgenSTA") #GxExM analysis - Multi-traits Selectio

#2. Set working directory
setwd("F:/DUC.Code/HTP_NUE")

#3. Load data file
dat <-read.csv("06.1_Figures/phenoRatio.csv")
attach(dat)

#4. Define Variables
factors <-c('PlantID', 'Treatment')
traits <- colnames(dat)[-c(1:5)]

#5. Convert variables into appropriate data types
dat$Group     <-as.factor(dat$Group)     # Group (ind, aus, japonica, admix) as factor
dat$Genotype  <-as.factor(dat$Genotype) # Treatment as factor
dat$Treatment <-as.factor(dat$Treatment) # Treatment as factor
dat$PlantID   <-as.factor(dat$PlantID)   # PlantID as factor
dat$Rep       <-as.factor(dat$Rep)       # Replication as factor
#dat$DAS       <-as.factor(dat$DAS)       # DAS as factor

### Histogram
library("dplyr")
library("ggplot2")
library("hrbrthemes")

w = 4
h = 5
###Full code for create histogram using ggpubr package
gghistogram(
  dat,
  x = "NUE.SI",
  y = "..count..",
  combine = FALSE,
  merge = FALSE,
  weight = NULL,
  color = "black",
  fill = "blue",
  palette = c("Dark2"),
  size = 0.5,
  linetype = "solid",
  alpha = 0.5,
  bins = 30,
  binwidth = NULL,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  facet.by = NULL,
  panel.labs = NULL,
  short.panel.labs = TRUE,
  add = "mean",               #("none", "mean", "median")
  add.params = list(linetype = "dashed"),
  rug = TRUE,
  add_density = TRUE,
  label = NULL,
  font.label = list(size = 16, color = "black"),
  label.select = NULL,
  repel = FALSE,
  label.rectangle = FALSE,
  position = position_identity(),
  ggtheme = theme_pubr()
) + xlab("NUE Stress Index")
ggsave('06.1_Figures/NUE.SI.png', width = w, height = h)

