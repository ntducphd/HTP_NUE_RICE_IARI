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
libraries("metan", "SpATS","statgenGxE", "statgenSTA") #GxExM analysis - Multi-traits Selection

#2. Set working directory
setwd("F:/DUC.Code/HTP_NUE")

#3. Load data file
dat <-read.csv("Phenotypes/manual.68DAS.csv")
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


## Load functions
source("Functions/R_rainclouds.R")
source("Functions/R_summarySE.R")
source("Functions/pca.plot.R")
source("Functions/analysis.R")
source("Functions/modeling.models.R")

#Gardner-Altman estimation plot-------------------------------------------------
#2. Gardner-Altman estimation plot and Box plot "dabest" package
#Nature Methods 2019, 1548-7105. 10.1038/s41592-019-0470-3
w = 5
h = 7
two.group.unpaired <- dat %>%dabest(Treatment, NUEg, idx = c("Control", "Nstress"), paired = FALSE)
two.group.unpaired.meandiff <- mean_diff(two.group.unpaired)
plot(two.group.unpaired.meandiff, color.column = Treatment,
                                       palette = "Dark2",
                                 tick.fontsize = 16.5,
                           axes.title.fontsize = 16.5,
                                  rawplot.type = "swarmplot",
                                   show.legend = FALSE,
                                         theme = ggplot2::theme_classic())
ggsave('04_Figures/NUEg-.png', width = w, height = h)

