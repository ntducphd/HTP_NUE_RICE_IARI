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
dat <-read.csv("Phenotypes/68-83DAS.csv")
attach(dat)

#4. Define Variables
factors <-c('Group', 'Treatment')
traits <- colnames(dat)[-c(1:6)]

#5. Convert variables into appropriate data types
dat$Group     <-as.factor(dat$Group)     # Group (ind, aus, japonica, admix) as factor
dat$Genotype  <-as.factor(dat$Genotype) # Treatment as factor
dat$Treatment <-as.factor(dat$Treatment) # Treatment as factor
dat$PlantID   <-as.factor(dat$PlantID)   # PlantID as factor
dat$Rep       <-as.factor(dat$Rep)       # Replication as factor
dat$DAS       <-as.factor(dat$DAS)       # DAS as factor

#Rainclouds Plot----------------------------------------------------------------
## Load functions
source("Functions/R_rainclouds.R")
source("Functions/R_summarySE.R")
source("Functions/pca.plot.R")
source("Functions/analysis.R")
source("Functions/modeling.models.R")

# width and height variables for saved plots
library(ggplot2)
w = 12
h = 4.5

# Biomass-related traits--------------------------------------------------------
### PSA
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "PSA.P", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p1 <- ggplot(dat, aes(x = DAS, y = PSA.P, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = PSA.P, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = PSA.P, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = PSA.P_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = PSA.P_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = PSA.P_mean, group = Treatment, colour = Treatment, ymin = PSA.P_mean-se, ymax = PSA.P_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("PSA.P (cm^2)")
p1
ggsave('Figures/PSA.P.png', width = w, height = h)

#Bio.volume
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "E.Biovolume", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p2 <- ggplot(dat, aes(x = DAS, y = E.Biovolume, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = E.Biovolume, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = E.Biovolume, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = E.Biovolume_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = E.Biovolume_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = E.Biovolume_mean, group = Treatment, colour = Treatment, ymin = E.Biovolume_mean-se, ymax = E.Biovolume_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("E.Biovolume (pixels)")
p2
ggsave('Figures/E.Biovolume.png', width = w, height = h)

#AreaSV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "AreaSV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p3 <- ggplot(dat, aes(x = DAS, y = AreaSV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = AreaSV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = AreaSV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = AreaSV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = AreaSV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = AreaSV_mean, group = Treatment, colour = Treatment, ymin = AreaSV_mean-se, ymax = AreaSV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("AreaSV (pixels)")
p3
ggsave('Figures/AreaSV.png', width = w, height = h)

#AreaTV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "AreaTV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p3 <- ggplot(dat, aes(x = DAS, y = AreaTV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = AreaTV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = AreaTV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = AreaTV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = AreaTV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = AreaTV_mean, group = Treatment, colour = Treatment, ymin = AreaTV_mean-se, ymax = AreaTV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("AreaTV (pixels)")
p3
ggsave('Figures/AreaTV.png', width = w, height = h)

#Architectural traits-----------------------------------------------------------
#iPH
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "iPH", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p4 <- ggplot(dat, aes(x = DAS, y = iPH, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = iPH, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = iPH, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = iPH_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = iPH_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = iPH_mean, group = Treatment, colour = Treatment, ymin = iPH_mean-se, ymax = iPH_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("iPH (cm)")
p4
ggsave('Figures/iPH.png', width = w, height = h)

#CHASV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "CHASV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p5 <- ggplot(dat, aes(x = DAS, y = CHASV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = CHASV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = CHASV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = CHASV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = CHASV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = CHASV_mean, group = Treatment, colour = Treatment, ymin = CHASV_mean-se, ymax = CHASV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("CHASV (pixels)")
p5
ggsave('Figures/CHASV.png', width = w, height = h)

#CHATV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "CHATV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p6 <- ggplot(dat, aes(x = DAS, y = CHATV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = CHATV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = CHATV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = CHATV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = CHATV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = CHATV_mean, group = Treatment, colour = Treatment, ymin = CHATV_mean-se, ymax = CHATV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("CHATV (pixels)")
p6
ggsave('Figures/CHATV.png', width = w, height = h)

#Physiological traits-----------------------------------------------------------
#Color.GPA.SV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "Color.GPA.SV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p7 <- ggplot(dat, aes(x = DAS, y = Color.GPA.SV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = Color.GPA.SV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = Color.GPA.SV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.GPA.SV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.GPA.SV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.GPA.SV_mean, group = Treatment, colour = Treatment, ymin = Color.GPA.SV_mean-se, ymax = Color.GPA.SV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("Color.GPA.SV (pixels)")
p7
ggsave('Figures/Color.GPA.SV.png', width = w, height = h)

#Color.GPA.TV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "Color.GPA.TV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p8 <- ggplot(dat, aes(x = DAS, y = Color.GPA.TV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = Color.GPA.TV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = Color.GPA.TV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.GPA.TV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.GPA.TV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.GPA.TV_mean, group = Treatment, colour = Treatment, ymin = Color.GPA.TV_mean-se, ymax = Color.GPA.TV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("Color.GPA.TV (pixels)")
p8
ggsave('Figures/Color.GPA.TV.png', width = w, height = h)

#Color.YPA.SV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "Color.YPA.SV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p9 <- ggplot(dat, aes(x = DAS, y = Color.YPA.SV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = Color.YPA.SV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = Color.YPA.SV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.YPA.SV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.YPA.SV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.YPA.SV_mean, group = Treatment, colour = Treatment, ymin = Color.YPA.SV_mean-se, ymax = Color.YPA.SV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("Color.YPA.SV (pixels)")
p9
ggsave('Figures/Color.YPA.SV.png', width = w, height = h)

#Color.YPA.TV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "Color.YPA.TV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p10 <- ggplot(dat, aes(x = DAS, y = Color.YPA.TV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = Color.YPA.TV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = Color.YPA.TV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.YPA.TV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.YPA.TV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = Color.YPA.TV_mean, group = Treatment, colour = Treatment, ymin = Color.YPA.TV_mean-se, ymax = Color.YPA.TV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("Color.YPA.TV (pixels)")
p10
ggsave('Figures/Color.YPA.TV.png', width = w, height = h)

#IR.SV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "IR.SV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p11 <- ggplot(dat, aes(x = DAS, y = IR.SV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = IR.SV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = IR.SV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = IR.SV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = IR.SV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = IR.SV_mean, group = Treatment, colour = Treatment, ymin = IR.SV_mean-se, ymax = IR.SV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("IR.SV")
p11
ggsave('Figures/IR.SV.png', width = w, height = h)

#IR.TV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "IR.TV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p12 <- ggplot(dat, aes(x = DAS, y = IR.TV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = IR.TV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = IR.TV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = IR.TV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = IR.TV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = IR.TV_mean, group = Treatment, colour = Treatment, ymin = IR.TV_mean-se, ymax = IR.TV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("IR.TV")
p12
ggsave('Figures/IR.TV.png', width = w, height = h)

#NIR.SV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "NIR.SV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p13 <- ggplot(dat, aes(x = DAS, y = NIR.SV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = NIR.SV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = NIR.SV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = NIR.SV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = NIR.SV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = NIR.SV_mean, group = Treatment, colour = Treatment, ymin = NIR.SV_mean-se, ymax = NIR.SV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("NIR.SV")
p13
ggsave('Figures/NIR.SV.png', width = w, height = h)

#NIR.TV
#load the repeated measures factorial data
sumrepdat <- summarySE(dat, measurevar = "NIR.TV", groupvars=c("DAS", "Treatment"))
head(sumrepdat)
p14 <- ggplot(dat, aes(x = DAS, y = NIR.TV, fill = Treatment)) +
  geom_flat_violin(aes(fill = Treatment),position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(DAS)-.15, y = NIR.TV, colour = Treatment),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(aes(x = DAS, y = NIR.TV, fill = Treatment),outlier.shape = NA, alpha = .5, width = .2, colour = "black")+
  geom_line(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = NIR.TV_mean, group = Treatment, colour = Treatment), linetype = 5)+
  geom_point(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = NIR.TV_mean, group = Treatment, colour = Treatment),size = 2, shape = 20) +
  geom_errorbar(data = sumrepdat, aes(x = as.numeric(DAS)+.1, y = NIR.TV_mean, group = Treatment, colour = Treatment, ymin = NIR.TV_mean-se, ymax = NIR.TV_mean+se), width = .05)+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size=14),
        axis.text = element_text(size=14, colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle=0, size=14, hjust=0.5),
        axis.text.y = element_text(angle=0, size=14, hjust=0.5)) +
  xlab("Day After Sowing") + ylab("NIR.TV")
p14
ggsave('Figures/NIR.TV.png', width = w, height = h)
