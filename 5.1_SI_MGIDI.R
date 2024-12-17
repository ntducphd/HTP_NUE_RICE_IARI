################################################################################
### Author: Nguyen Trung Duc
### Institute:ICAR- Indian Agricultural Research Institute, New Delhi, India
###           Nanaji Deshmukh Plant Phenomics Centre, ICAR-IARI,New Delhi, India         
###           Vietnam National University of Agriculture, Hanoi, Vietnam
### Email: ntduc11@gmail.com
################################################################################
#===============================================================================
###Selection index based on MGIDI - multi-trait genotype-ideotype distance index
#===============================================================================

#Clean the global environment
rm(list = ls())

#1. Load library----------------------------------------------------------------
library("easypackages") #Load multiple libraries
libraries("metan", "tidyverse","future.apply","FactoMineR","factoextra","corrplot", "ade4","RColorBrewer",
         "pheatmap", "ComplexHeatmap","heatmaply","ggplot2","ggrepel", "rstatix","ggpubr",
         "FWDselect","ggforce","GGally","lmerTest")
libraries("xlsx", "rJava", "fmsb") #Export the results
libraries("grid","gridExtra","gtable","ggplot2","lattice") #Arrange plot
library("gghighlight") #Easy Way to Highlight a GGPlot in R
#2.1. Set working directory-----------------------------------------------------
setwd("F:/DUC_Code/HTP_NUE/05_SI-MGIDI")

#2.2. Load functions to_factor--------------------------------------------------
to_factor <- function(.data, ...){
  return(mutate(.data, across(c(...), as.factor)))
}

#3. Load dataset----------------------------------------------------------------
dat1 <-read.csv("Data/SI_Control.csv")
dat2 <-read.csv("Data/SI_Nstress.csv")
#Define factors (metan package)
df1 <- to_factor(dat1, 1:2) # 2 first columns as factor
df2 <- to_factor(dat2, 1:2) # 2 first columns as factor
str(df1)
str(df2)

#Correlation of dataset---------------------------------------------------------
p1 <-corr_coef(df1) %>% plot()
p2 <-corr_coef(df2) %>% plot()
arrange_ggplot(p1, p2)
library(GGally)

#Correlogram
dat3 <-read.csv("Data/SI.csv")
dat3$Treatment <- factor(dat3$Treatment, exclude = c("", NA))
library(GGally)
png("cor_all.png", width=10, height=6, units="in", res=600)
dat3 %>%
  drop_na(GrainWt, Biomass, NUpE, NUtE, NHI, Treatment) %>%
  select(GrainWt, Biomass, NUpE, NUtE, NHI, Treatment) %>%
  ggpairs(progress = FALSE,
          ggplot2::aes(alpha = 0.65, color = Treatment))+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")+
  theme_bw()
dev.off()

#MGIDI--------------------------------------------------------------------------
#Control
mod1 <- gamem(df1,
             gen = Genotype,
             rep = Rep,
             resp = everything())
mgidi_index1 <- mgidi(mod1, SI = 10) # Selection intensity
mgidi_index1
mgidi_index1$sel_gen
#Genetic parameters
gen_par <- gmd(mod1) # what = "genpar" is the default value
gen_par

###Stress
mod2 <- gamem(df2,
              gen = Genotype,
              rep = Rep,
              resp = everything())
mgidi_index2 <- mgidi(mod2, SI = 10) # Selection intensity
mgidi_index2
mgidi_index2$sel_gen
#Genetic parameters
gen_par <- gmd(mod2) # what = "genpar" is the default value
gen_par
#Exploiting the strengths and weaknesses of selected genotypes------------------
#Control
plot(mgidi_index1,
     radar = TRUE,
     type = "contribution",
     genotypes = "selected", #c ("selected", "all")
     width.bar = 1,
     arrange.label = TRUE,
     check.overlap = TRUE,
     size.text = 11,
     rotate = TRUE,
     size.line = 1)

#Nstress
plot(mgidi_index2,
     radar = TRUE,
     type = "contribution",
     width.bar = 1,
     arrange.label = TRUE,
     check.overlap = FALSE,
     size.text = 11,
     rotate = TRUE,
     size.line = 1)

#Venn plot to show the relationships between the indexes--------------
Control <- gmd(mgidi_index1, "sel_gen")
NStress <- gmd(mgidi_index2, "sel_gen")
# Create the plot
venn_plot(Control, NStress, show_elements = FALSE,
          name_size = 8,
          text_size = 8)

#Plot MGIDI selection index-----------------------------------------------------
# Load function
# source("MGIDI_functions.R")
# Plot MGIDI index
plot(mgidi_index1,
           arrange.label = TRUE,
           size.point = 2.5,
           size.line = 0.7,
           size.text = 8,
           col.sel = "red",
           col.nonsel = "black")
plot(mgidi_index2,
           arrange.label = TRUE,
           size.point = 2.5,
           size.line = 0.7,
           size.text = 8,
           col.sel = "red",
           col.nonsel = "black")

# Plot Proportion of phenotypic variance----------------------------------------
#Control 
p <-plot(mod1, type = "vcomp") + geom_hline(yintercept = 0.5, linetype = 2)
# Vertical rotation of x axis text
p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#Nstress
p <-plot(mod2, type = "vcomp") + geom_hline(yintercept = 0.5, linetype = 2)
# Vertical rotation of x axis text
p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Phenotypic and genotypic variance-covariance correlation matrices-----------
### Control
gcov <- gmd(mod1, what = "gcov")
round_cols(gcov, digits = 0) # to fit in the page width
pcov <- gmd(mod1, what = "pcov")
round_cols(pcov, digits = 0) # to fit in the page width
gcor <- gmd(mod1, what = "gcor") # genotypic correlation
pcor <- gmd(mod1, what = "pcor") # phenotypic correlation
pcor_gcor1 <- make_lower_upper(pcor, gcor, diag = 1)
p1 <-corrplot(pcor_gcor1,
              method = 'pie',
              order = "AOE",
              addCoef.col = "black",
              insig = "p-value",
              diag = FALSE)

### Nstress
gcov <- gmd(mod2, what = "gcov")
round_cols(gcov, digits = 0) # to fit in the page width
pcov <- gmd(mod2, what = "pcov")
round_cols(pcov, digits = 0) # to fit in the page width
gcor <- gmd(mod2, what = "gcor") # genotypic correlation
pcor <- gmd(mod2, what = "pcor") # phenotypic correlation
pcor_gcor2 <- make_lower_upper(pcor, gcor, diag = 1)
p2 <-corrplot(pcor_gcor2,
              method = 'pie',
              order = "AOE",
              addCoef.col = "black",
              insig = "p-value",
              diag = FALSE)

#===============================================================================
### Compare with FAI-BLUP and SH index
### Control condition-----------------------------------------------------------
#MGIDI index
mgidi_index1 <- mgidi(mod1, SI = 10) #Selection intensity 10%

#FAI-BLUP
fai1 <-fai_blup(mod1, SI = 10) #Selection intensity
fai1$sel_gen

##Smith-Hazel with all variables
smith1 <-Smith_Hazel(mod1, SI = 10)#Selection intensity
smith1$sel_gen
#Comparing the indexes
coincidence1 <- coincidence_index(mgidi_index1, fai1, smith1, total = 300)
#Venn plot to show the relationships between the MGIDI, FAI, SH indexes
MGIDI <- gmd(mgidi_index1, "sel_gen")
FAI_BLUP <- gmd(fai1, "sel_gen")
Smith_Hazel <- gmd(smith1, "sel_gen")
venn_plot(MGIDI, FAI_BLUP, Smith_Hazel, show_elements = FALSE,
          name_size = 8,
          text_size = 8)

##Join all selection indexes under Control condition --------------------------
sg_mgidi1 <-
  mgidi_index1$sel_dif %>%
  select_cols(Factor, VAR, Xo, sense, h2, SGperc) %>%
  rename(MGIDI = SGperc)
sg_fai1 <-
  fai1$sel_dif_trait$ID1 %>%
  select_cols(VAR, SGperc) %>%
  rename(FAI_BLUP = SGperc)
sg_hz_all1 <-
  smith1$sel_dif %>%
  select_cols(VAR, SGperc) %>%
  rename(HZ1 = SGperc)
# BLUPs
blups <-
  gmd(mod1, "blupg") %>%
  desc_stat(stats = c("mean, se")) %>%
  round_cols() %>%
  add_cols(mse = paste(mean, "+-", se, sep = "")) %>%
  select(VAR = variable, mse = mse)
# Selection gain
selection_gain <-
  left_join(sg_mgidi1, blups, by = "VAR") %>%
  left_join(sg_fai1, by = "VAR") %>%
  left_join(sg_hz_all1, by = "VAR") %>%
  select_cols(Factor, VAR, mse, MGIDI, FAI_BLUP, HZ1) %>%
  tidy_strings(Factor, sep = "")
selection_gain

### Nstress condition-----------------------------------------------------------
#MGIDI index
mgidi_index2 <- mgidi(mod2, SI = 10) #Selection intensity
#FAI-BLUP
fai2 <-fai_blup(mod2, SI = 10) #Selection intensity
fai2$sel_gen
##Smith-Hazel with all variables
smith2 <-Smith_Hazel(mod2, SI = 10)#Selection intensity
smith2$sel_gen
#Comparing the indexes
coincidence2 <- coincidence_index(mgidi_index2, fai2, smith2, total = 300)
coincidence2
#Venn plot to show the relationships between the MGIDI, FAI, SH indexes
MGIDI <- gmd(mgidi_index2, "sel_gen")
FAI_BLUP <- gmd(fai2, "sel_gen")
Smith_Hazel <- gmd(smith2, "sel_gen")
venn_plot(MGIDI, FAI_BLUP, Smith_Hazel, show_elements = FALSE,
          name_size = 8,
          text_size = 8)
##Join all selection indexes under Nstress condition ---------------------------
sg_mgidi2 <-
  mgidi_index2$sel_dif %>%
  select_cols(Factor, VAR, Xo, sense, h2, SGperc) %>%
  rename(MGIDI = SGperc)
sg_fai2 <-
  fai2$sel_dif_trait$ID1 %>%
  select_cols(VAR, SGperc) %>%
  rename(FAI_BLUP = SGperc)
sg_hz_all2 <-
  smith2$sel_dif %>%
  select_cols(VAR, SGperc) %>%
  rename(HZ1 = SGperc)
# BLUPs
blups <-
  gmd(mod2, "blupg") %>%
  desc_stat(stats = c("mean, se")) %>%
  round_cols() %>%
  add_cols(mse = paste(mean, "+-", se, sep = "")) %>%
  select(VAR = variable, mse = mse)
# Selection gain
selection_gain <-
  left_join(sg_mgidi2, blups, by = "VAR") %>%
  left_join(sg_fai2, by = "VAR") %>%
  left_join(sg_hz_all2, by = "VAR") %>%
  select_cols(Factor, VAR, mse, MGIDI, FAI_BLUP, HZ1) %>%
  tidy_strings(Factor, sep = "")
selection_gain

#===============================================================================
#Visual phenotype---------------------------------------------------------------
#Generate plot function
plot_selected <- function(data,
                          gen,
                          var,
                          xlab = "Genotype",
                          ylab = expression(paste("Grain yield"))){
  ovmean <- desc_stat(data, {{var}}, stats = "mean") %>% pull()
  plot_bars(data, {{gen}}, {{var}},
            xlab = xlab,
            ylab = ylab,
            size.text = 10, #6 for pdf file
            size.line = 0.2,
            width.bar = 1,
            n.dodge = 2,
            order = "desc") +
    geom_hline(yintercept = ovmean, linetype = 2, color = "red", size = 0.2)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

#Control condition
plot_selected(dat1, Genotype, GrainWt, ylab = "Grain weight (g)")
plot_selected(dat1, Genotype, Biomass, ylab = "Biomass (g)")
plot_selected(dat1, Genotype, NUpE, ylab = "Nitrogen Uptake Efficiency")
plot_selected(dat1, Genotype, NUtE, ylab = "Nitrogen Utilization Efficiency")
plot_selected(dat1, Genotype, NHI, ylab = "Nitrogen Harvest Index")

#Nstress condition
plot_selected(dat2, Genotype, GrainWt, ylab = "Grain weight (g)")
plot_selected(dat2, Genotype, Biomass, ylab = "Biomass (g)")
plot_selected(dat2, Genotype, NUpE, ylab = "Nitrogen Uptake Efficiency")
plot_selected(dat2, Genotype, NUtE, ylab = "Nitrogen Utilization Efficiency")
plot_selected(dat2, Genotype, NHI, ylab = "Nitrogen Harvest Index")

#Plot FAI-BLUP and SH 
#Control
pdf('FAI-BLUP-Control.pdf', width=8, height=8, onefile=F)
plot(fai1,
     arrange.label = TRUE,
     size.point = 2.5,
     size.line = 0.7,
     size.text = 6,
     col.sel = "red",
     col.nonsel = "black")
dev.off()

plot(smith1,
     arrange.label = TRUE,
     size.point = 2.5,
     size.line = 0.7,
     size.text = 8,
     col.sel = "red",
     col.nonsel = "black")

#Nstress

plot(fai2,
     arrange.label = TRUE,
     size.point = 2.5,
     size.line = 0.7,
     size.text = 8,
     col.sel = "red",
     col.nonsel = "black")

plot(smith2,
     arrange.label = TRUE,
     size.point = 2.5,
     size.line = 0.7,
     size.text = 8,
     col.sel = "red",
     col.nonsel = "black")

