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
###PART 01. HTP DATA ANALYSIS (REML, BLUP, MODELING, DONORS SELECTION)
#===============================================================================

# Remove previous work----------------------------------------------------------
rm(list=ls())

#1. Load library----------------------------------------------------------------
library("easypackages") #Load multiple libraries
libraries( "rstatix", "ggpubr","stringr", "ggplot2", "tidyverse",
          "data.table", "readr","plotly", "DT", "pheatmap", "VennDiagram","patchwork",
          "heatmaply", "ggcorrplot", "RColorBrewer", "hrbrthemes", "tm", "proustr")
libraries("plyr", "dplyr", "reshape2", "readxl","stringr", "lubridate") #Data transformation
libraries("FactoMineR", "factoextra","corrplot","ade4","ComplexHeatmap", "ggrepel","dabestr")#data-visualization 
libraries("xlsx", "rJava", "fmsb","rlang", "car")      #Read data file
libraries("VIM", "outliers", "mice")                   #Outliers detection and management
libraries("asremlPlus", "nlme","lme4", "FW", "sommer") #Linear mixed models
libraries("metan", "SpATS","statgenGxE", "statgenSTA") #GxExM analysis - Multi-traits Selection

#2. Set working directory
setwd("F:/DUC_Code/HTP_NUE")

#3. Load data file
dat <-read.csv("Phenotypes/manual.csv")
attach(dat)

#4. Define Variables------------------------------------------------------------
#df <- to_factor(dat, 1:5) # 5 first columns as factor
#str(df)
#factors <-c('PlantID', 'Treatment')
#traits <- colnames(dat)[-c(1:5)]
#dat.new <- dat[c(factors, traits)]
#attach(dat.new)
#Convert variables into appropriate data types
dat$Group    <-as.factor(dat$Group)      # Group (ind, aus, japonica, basmati, admix) as factor
dat$Genotype <-as.factor(dat$Genotype)  # Genotype as factor
dat$Treatment <-as.factor(dat$Treatment)  # Treatment as factor
dat$PlantID   <-as.factor(dat$PlantID)    # PlantID as factor
dat$Rep       <-as.factor(dat$Rep)        # Replication as factor

#===============================================================================
###                       Descriptive Statistics
#===============================================================================
# Whole data frame - rstatix package
##Grouped data by treatment
treatment.sum <-dat %>%
  group_by(Treatment) %>% 
  get_summary_stats(type = "common")
write.xlsx(treatment.sum, file = "03_REML/Treatment_Sum.xlsx",
           sheetName = "TreatmentSum", append = FALSE)

#===============================================================================
###                    ANOVA - MIX MODEL with REML/BLUP
#===============================================================================
### ANOVA - Joint analysis of variance------------------------------------------
## Performs a joint analysis of variance to check for the presence of genotype-vs-environment interactions using both randomized complete block and alpha-lattice designs.
# ANOVA for all variables in data
anova.dat <- anova_joint(dat,    #Dataset containing the columns related to Treatment, Genotypes, replication/block and response variable(s)
           env = Treatment,      #Levels of the environments
           gen = PlantID,        #Levels of the genotypes.
           rep = Rep,            #Levels of the replications/blocks.
          resp = everything(),   #The response variable(s)
        #block = NULL            #All effects, except the error, are assumed to be fixed
          )
plot(anova.dat)
# Report ANOVA result
ANOVA_report<-gmd(anova.dat, "details")
### Write the result to excel file
write.xlsx(ANOVA_report, file = "03_REML/ANOVA_Sum.xlsx",
           sheetName = "ANOVA_Sum", append = FALSE)

## REML/BLUP - Mixed-effect models Variance components, genetic parameters, random effects and BLUP prediction in multi-environment trials
## Mixed-effect or random-effect models
REML.dat <-  gamem_met(dat,        # dataset:Environments, Genotypes, replication/block and response variable(s)
                    env = Treatment,     # levels of the environments/treatment
                    gen = PlantID,       # levels of the genotypes
                    rep = Rep,           # levels of the replications/blocks
                    resp = everything(), # response variable(s) resp = c(var1, var2, var3) or resp = everything()
                    # block = NULL,      # defaults to NULL. In this case, a randomized complete block design is considered. If block is informed, then an alpha-lattice design is employed considering block as random to make use of inter-block information, whereas the complete replicate effect is always taken as fixed, as no inter-replicate information was to be recovered (Mohring et al., 2015)
                    # by = NULL,         # One variable (factor) to compute the function by. It is a shortcut to group_by().This is especially useful, for example, when the researcher want to analyze environments within mega-environments. In this case, an object of class waasb_grouped is returned.
                    random = "all",      # random = "gen" (Model 1), "env" (Model 2), "all" (Model 3), 
                    # prob = 0.05,       # The probability for estimating confidence interval for BLUP's prediction
                    verbose = TRUE)      # Logical argument. If verbose = FALSE the code will run silently

#Variance contribution plot
w = 12.8
h = 4
p <-plot(REML.dat, type = "vcomp")+ geom_hline(yintercept = 0.5, linetype = 2)
p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('03_REML/Manual-REML.png', width = w, height = h)

#Genetic parameters
genpar <-gmd(REML.dat, "genpar", type = "ENV")
genpar
write.xlsx(genpar, file = "03_REML/genpar_manual.xlsx",
           sheetName = "genpar", append = FALSE)
#Variance components
vcomp <- gmd(REML.dat, what = "vcomp")
vcomp
write.xlsx(vcomp, file = "03_REML/vcomp_manual.xlsx",
           sheetName = "vcomp", append = FALSE)
