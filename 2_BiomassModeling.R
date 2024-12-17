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
###                         Biomass modeling
#===============================================================================

# Remove previous work----------------------------------------------------------
rm(list=ls())
setwd("E:/MSc.Thesis/Codes/HTP_NUE")

#Method 1. Single variable model
library("ggstatsplot")
library("ggplot2")
library("metan")
#Load data
dat <-read.csv("02_Biomass_model/Data/Rice19/RIL35DAT.csv")

#Method 1. Single variable model------------------------------------------------
#Klukas et al. (2014)-----------------------------------------------------------
#Fresh weight
p1 <-ggstatsplot::ggscatterstats(
  data = dat, 
  x = FW, 
  y = E.Biovolume,
  title = "Fresh weight prediction model",
  messages = FALSE,
  #label.var = "PlantID"
)
#Dry weight
p2 <-ggstatsplot::ggscatterstats(
  data = dat, 
  x = DW, 
  y = E.Biovolume,
  title = "Dry weight prediction model",
  messages = FALSE
)
#Berger et al. (2012)-----------------------------------------------------------
#Fresh weight
p3 <-ggstatsplot::ggscatterstats(
  data = dat, 
  x = FW, 
  y = PSA.B,
  title = "Fresh weight prediction model",
  messages = FALSE,
  #label.var = "PlantID"
)
#Dry weight
p4 <-ggstatsplot::ggscatterstats(
  data = dat, 
  x = DW, 
  y = PSA.B,
  title = "Dry weight prediction model",
  messages = FALSE
)
#Parent et al. (2015)-----------------------------------------------------------
#Fresh weight
p5 <-ggstatsplot::ggscatterstats(
  data = dat, 
  x = FW, 
  y = PSA.P,
  title = "Fresh weight prediction model",
  messages = FALSE
)
#Dry weight
p6 <-ggstatsplot::ggscatterstats(
  data = dat, 
  x = DW, 
  y = PSA.P,
  title = "Dry weight prediction model",
  messages = FALSE
)
arrange_ggplot(p1, p2)
arrange_ggplot(p3, p4)
arrange_ggplot(p5, p6)
#arrange_ggplot(p1, p3, p5, p2, p4, p6)

#Chen et al. (2018)-------------------------------------------------------------
