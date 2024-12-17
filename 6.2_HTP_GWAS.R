################################################################################
### Author: Nguyen Trung Duc
### Institute: ICAR- Indian Agricultural Research Institute, New Delhi, India
###            Nanaji Deshmukh Plant Phenomics Centre, ICAR-IARI, New Delhi, India 
###            Vietnam National University of Agriculture, Hanoi, Vietnam
### Email: ntduc11@gmail.com
################################################################################
#===============================================================================
#                       Paper: GWAS - HTP - NUE - RICE
#===============================================================================
###PART 02. GENOME-WIDE ASSOCIATION STUDY AND QTL MAPPING
#===============================================================================

# Remove previous work
rm(list=ls())

#1.Preparation
library("multtest")
library("gplots")
library("LDheatmap")
library("genetics")
library("MASS")
library("compiler")
library("RColorBrewer")
library("scatterplot3d") 

#2. Install packages
install.packages("rMVP")
install.packages("GAPIT3")
install.packages("CMplot")
#install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
#devtools::install_github("xiaolei-lab/rMVP")

## Load library
library("rMVP")
library("GAPIT3")
library("CMplot")

# Loading packages for GAPIT and GAPIT functions
source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
#source("GAPIT3/GAPIT.library.R")
#source("GAPIT3/gapit_functions.txt")

#3. Set working directory
setwd("D:/Projects/2024_HTP_NUE/3_Codes/06.2_GWAS")
#Step 1: Set data directory and import files
myY <- read.csv("K50_Phenotypes.csv", head = TRUE)
myG <- read.csv("147SNPgenoypes.csv", head = FALSE)

#4. Convert HapMap format to numerical
myGAPIT <- GAPIT(G=myG, output.numerical=TRUE)
myGD= myGAPIT$GD
myGM= myGAPIT$GM

#5. GWAS analysis
## Model selection
myGAPIT <- GAPIT(
  Y=myY[,c(1,2)],
  G=myG,
  SNP.MAF=0.05,
  PCA.total=3,
  NJtree.group=4,
  Model.selection = TRUE)

## GWAS with five methods
myGAPIT_MLM <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  SNP.MAF=0.05,
  PCA.total=3,
  model=c("MLM", "CMLM", "FarmCPU", "Blink", "SUPER"),
  Multiple_analysis=TRUE)

#Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,c(1,2)],
  G=myG,
  SNP.MAF=0.05,
  PCA.total=3,
  NJtree.group=4,
  kinship.cluster=c("average", "complete", "ward"),
  kinship.group=c("Mean", "Max"),
  group.from=200,
  group.to=1000000,
  group.by=10)

#Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model=c("GLM","MLM","SUPER","MLMM","FarmCPU","Blink"),# choose model
  PCA.total=5,                                          # set total PCAs
  NJtree.group=8,                                       # set the number of clusting group in Njtree plot
  Inter.Plot=TRUE,                                      # perform interactive plot
  Multiple_analysis=TRUE,                               # perform multiple analysis
  PCA.3d=TRUE,                                          # plot 3d interactive PCA
  file.output=T)
  #,QTN.position=mysimulation$QTN.position)

#GLM - The GAPIT use Least Squares to solve the module
myGAPIT_GLM <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="GLM",
  PCA.total=5,
  SNP.MAF=0.05,
  file.output=T)

#MLM - EMMA method is used in GAPIT,
myGAPIT_MLM <- GAPIT(
  Y=myY,#[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="MLM",
  PCA.total=5,
  NJtree.group=8,
  SNP.MAF=0.05,
  file.output=T
)

#CMLM - Compress Mixed Linear Model is published by Zhang in 2010.
myGAPIT_CMLM <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="CMLM",
  PCA.total=5,
  SNP.MAF=0.05,
  file.output=T)

#MLMM - Multiple Loci Mixied linear Model is published by Segura in 2012
myGAPIT_MLMM <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="MLMM",
  PCA.total=5,
  SNP.MAF=0.05,
  file.output=T)

#SUPER - Settlement of MLM Under Progressively Exclusive Relation- ship is published by Qishan in 2014
myGAPIT_SUPER <- GAPIT(
  Y=myY,#[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="SUPER",
  PCA.total=5,
  SNP.MAF=0.05,
  file.output=T)

# Farm-CPU - Fixed and random model Circulating Probability Unification (FarmCPU) is published by Xiaolei in 2016
myGAPIT_FarmCPU <- GAPIT(
  Y=myY,#[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="FarmCPU",
  PCA.total=5,
  SNP.MAF=0.05,
  file.output=T)

# Blink - Bayesian-information and Linkage-disequilibrium 34 Iteratively Nested Keyway
myGAPIT_FarmCPU <- GAPIT(
  Y=myY,#[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="FarmCPU",
  PCA.total=5,
  SNP.MAF=0.05,
  Inter.Plot=TRUE,           # perform interactive plot
  Multiple_analysis=TRUE,    # perform multiple analysis
  PCA.3d=TRUE,               # plot 3d interactive PCA
  file.output=T)

#6. Reports
library("CMplot")
#latest version of CMPlot
source("https://raw.githubusercontent.com/YinLiLin/CMplot/master/R/CMplot.r")

## SNP-density plot
CMplot(pig60K,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
       file.output=TRUE,verbose=TRUE,width=9,height=6)

## PCA plot
pca <- attach.big.matrix("mvp.pc.desc")[, 1:3]
#pca <- prcomp(t(as.matrix(genotype)))$x[, 1:3]
MVP.PCAplot(PCA=pca, Ncluster=3, class=NULL, col=c("red", "green", "yellow"), file.type="jpg")

## CMpot - Circular and Rectangular Manhattan Plot
# Circular-Manhattan plot
# Genome-wide association study(GWAS)
CMplot(pig60K,type="p",plot.type="c",chr.labels=paste("Chr",c(1:18,"X","Y"),sep=""),r=0.4,cir.legend=TRUE,
         outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
         memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)
#Attach chromosome density on the bottom of Manhattan plot
CMplot(pig60K, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(1e-6,1e-4),threshold.lty=c(1,2),
         threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
         chr.den.col=c("darkgreen", "yellow", "red"),signal.col=c("red","green"),signal.cex=c(1.5,1.5),
         signal.pch=c(19,19),file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
         width=14,height=6)
#Add genes or SNP names around the highlighted SNPs
SNPs <- pig60K[pig60K[,5] < (0.05 / nrow(pig60K)), 1]
genes <- paste("GENE", 1:length(SNPs), sep="_")
set.seed(666666)
CMplot(pig60K[,c(1:3,5)], plot.type="m",LOG10=TRUE,col=c("grey30","grey60"),highlight=SNPs,
         highlight.col=c("red","blue","green"),highlight.cex=1,highlight.pch=c(15:17), highlight.text=genes,      
         highlight.text.col=c("red","blue","green"),threshold=0.05/nrow(pig60K),threshold.lty=2,   
         amplify=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)

#Multi_tracks Rectangular-Manhattan plot
SNPs <- list(
  pig60K$SNP[pig60K$trait1<1e-6],
  pig60K$SNP[pig60K$trait2<1e-6],
  pig60K$SNP[pig60K$trait3<1e-6])
CMplot(pig60K, plot.type="m",multracks=TRUE,threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green","blue"),
       signal.cex=1, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1.4)
#Note: if you are not supposed to change the color of signal, 
#      please set signal.col=NULL and highlight.col=NULL.

## Q-Q plot for single trait or method
MVP.Report(pig60K,plot.type="q",conf.int.col=NULL,box=TRUE,file.type="jpg",memo="",dpi=300)
## Q-Q plot for multiple traits or methods
MVP.Report(imMVP,plot.type="q",col=c("dodgerblue1", "olivedrab3", "darkgoldenrod1"),threshold=1e6,
             signal.pch=19,signal.cex=1.5,signal.col="red",conf.int.col="grey",box=FALSE,multracks=
               TRUE,file.type="jpg",memo="",dpi=300)

#7. QTL mapping from GWAS


#8. Genomics Prediction/Selection 
set.seed(198521)
nfold=5
sets=sample(cut(1:nrow(myY ),nfold,labels=FALSE),nrow(myY ))
j=1
training=myY[,c(1,2)]
training[sets==j,2]=NA
training_index=is.na(training[,2])
testing=myY[training_index,c(1,2)]
myGAPIT<- GAPIT(
  Y=training,
  GD=myGD,
  GM=myGM,
  PCA.total=3,
  model=c("gBLUP")# could be "gBLUP", "cBLUP" or "sBLUP"
)
inf.pred=merge(myGAPIT$Pred[myGAPIT$Pred[,3]==2,c(1,3,5,8)],myY,by.x="Taxa",by.y="Taxa")
cor(inf.pred[,4],inf.pred[,5])


#gBLUP - gBLUP used marker kinship to replace the pedgree relationship matrix.
myGAPIT_gBLUP <- GAPIT(
  Y=myY,#[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="gBLUP",
  PCA.total=5,
  file.output=T
)

#cBLUP - cBLUP used group kinship to replace the individual matrix.
myGAPIT_cBLUP <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="cBLUP",
  PCA.total=5,
  file.output=T
)
#sBLUP - sBLUP used SUPER method to build psedue QTN kinship matrix.
myGAPIT_sBLUP <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="sBLUP",
  PCA.total=5,
  file.output=T
)

#===============================================================================
###                 rMVP - Genome-Wide Association Study
#===============================================================================
#1. Hapmap file
#1.1 Load data file
MVP.Data(fileHMP="147geno_50K_hapmap.txt",
         filePhe="50K.Manual.Stress.txt",
         sep.hmp="\t",
         sep.phe="\t",
         SNP.effect="Add",
         fileKin=TRUE,
         filePC=TRUE,
         #priority="memory",
         #maxLine=10000,
         out="mvp.hmp"
)

#1.2. Population structure
pca <- attach.big.matrix("mvp.hmp.pc.desc")[, 1:3]
#pca <- prcomp(t(as.matrix(genotype)))$x[, 1:3]
MVP.PCAplot(PCA=pca, Ncluster=3, class=NULL, col=c("red", "green", "yellow"), file.type="jpg")

#1.3. Start GWAS
#1.3.1. Load needed file 
genotype <- attach.big.matrix("mvp.hmp.geno.desc")
phenotype <- read.table("50K.StressIndex.txt",head=TRUE)
map <- read.table("mvp.hmp.geno.map" , head = TRUE)

#If you have more than one phenotype
for(i in 2:ncol(phenotype)){
  imMVP <- MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
    #K=Kinship,
    #CV.GLM=Covariates,
    #CV.MLM=Covariates,
    #CV.FarmCPU=Covariates,
    nPC.GLM=5,
    nPC.MLM=3,
    nPC.FarmCPU=3,
    priority="speed",
    #ncpus=10,
    vc.method="BRENT",
    maxLoop=10,
    method.bin="static",
    #permutation.threshold=TRUE,
    #permutation.rep=100,
    threshold=0.05,
    method=c("GLM", "MLM", "FarmCPU")
  )
  gc()
}
#1. Prepare data set (plink, vcf, halfpmap)
MVP.Data(fileBed="rMVP/VTrice_thinned_maf",
         filePhe=NULL,
         fileKin=TRUE,
         filePC=TRUE,       
         #priority="speed",
         #maxLine=10000,
         out="mvp.plink")


#Population structure
pca <- attach.big.matrix("mvp.plink.pc.desc")[, 1:3]
#pca <- prcomp(t(as.matrix(genotype)))$x[, 1:3]
MVP.PCAplot(PCA=pca, Ncluster=3, class=NULL, col=c("red", "green", "yellow"), file.type="jpg")

#2. Start GWAS
#2.1.Load needed file 
genotype <- attach.big.matrix("mvp.plink.geno.desc")
map <- read.table("mvp.plink.geno.map" , head = TRUE)
phenotype <- read.csv("3K.68DAS.Control.csv",head=TRUE)


#2.2. Run GWAS
imMVP <- MVP(
  phe=phenotype,
  geno=genotype,
  map=map,
  #K=Kinship,
  #CV.GLM=Covariates,       ## if you have additional covariates, please keep there open.
  #CV.MLM=Covariates,
  #CV.FarmCPU=Covariates,
  nPC.GLM=5,                ## if you have added PC into covariates, please keep there closed.
  nPC.MLM=3,
  nPC.FarmCPU=3,
  priority="speed",         ## for Kinship construction
  #ncpus=10,
  vc.method="BRENT",        ## only works for MLM
  maxLoop=10,
  method.bin="static",      ## "FaST-LMM", "static" (#only works for FarmCPU)
  #permutation.threshold=TRUE,
  #permutation.rep=100,
  threshold=0.05,
  method=c("GLM", "MLM", "FarmCPU"))

#If you have more than one phenotype
for(i in 2:ncol(phenotype)){
  imMVP <- MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
    #K=Kinship,
    #CV.GLM=Covariates,
    #CV.MLM=Covariates,
    #CV.FarmCPU=Covariates,
    nPC.GLM=5,
    nPC.MLM=3,
    nPC.FarmCPU=3,
    priority="speed",
    #ncpus=10,
    vc.method="BRENT",
    maxLoop=10,
    method.bin="static",
    #permutation.threshold=TRUE,
    #permutation.rep=100,
    threshold=0.05,
    method=c("GLM", "MLM", "FarmCPU")
  )
  gc()
}

#===============================================================================
###                   statgenGWAS - Single trait GWAS 
#===============================================================================
#1.Install package
install.packages("statgenGWAS")
#remotes::install_github("Biometris/statgenGWAS", ref = "develop", dependencies = TRUE)

#2. Load library
library("statgenGWAS")

#3.Example
library(statgenGWAS)
## Read data.
data("dropsMarkers")
data("dropsMap")
data("dropsPheno")

## Create gData object
## Add genotypes as row names of dropsMarkers and drop Ind column.
rownames(dropsMarkers) <- dropsMarkers[["Ind"]]
dropsMarkers <- dropsMarkers[colnames(dropsMarkers) != "Ind"]

## Add genotypes as row names of dropsMap.
rownames(dropsMap) <- dropsMap[["SNP.names"]]
## Rename Chomosome and Position columns.
colnames(dropsMap)[match(c("Chromosome", "Position"), colnames(dropsMap))] <- c("chr", "pos")

## Rename Variety_ID in phenotypic data to genotype.
colnames(dropsPheno)[colnames(dropsPheno) == "Variety_ID"] <- "genotype"
## Select relevant columns and convert data to a list.
dropsPhenoList <- split(x = dropsPheno[c("genotype", "grain.yield",
                                         "grain.number", "seed.size",
                                         "anthesis", "silking", "plant.height",
                                         "tassel.height", "ear.height")], 
                        f = dropsPheno[["Experiment"]])

## Create a gData object all data.
gDataDrops <- createGData(geno = dropsMarkers, map = dropsMap, pheno = dropsPhenoList)

#=======================================================================
#=================== Recoding and cleaning of markers ==================
#=======================================================================
## Remove duplicate SNPs from gDataDrops.
gDataDropsDedup <- codeMarkers(gDataDrops, impute = FALSE, verbose = TRUE) 

## Copy gData object.
gDataDropsMiss <- gDataDrops
## Add random missing values to 1% of the values in the marker matrix.
set.seed(1)
nVal <- nrow(gDataDropsMiss$markers) * ncol(gDataDropsMiss$markers)
gDataDropsMiss$markers[sample(x = 1:nVal, size = nVal / 100)] <- NA

## Impute missing values with random value.
## Remove SNPs and genotypes with proportion of NA larger than 0.01.
gDataDropsImputed <- codeMarkers(gData = gDataDropsMiss, 
                                 nMissGeno = 0.01, 
                                 nMiss = 0.01, 
                                 impute = TRUE, 
                                 imputeType = "random", 
                                 verbose = TRUE)

gDataDropsImputedBeagle <- codeMarkers(gData = gDataDropsMiss, 
                                       impute = TRUE,
                                       imputeType = "beagle",
                                       verbose = TRUE)

#=======================================================================
#========================== Single trait GWAS ==========================
#=======================================================================
## Run single trait GWAS for traits 'grain.yield' and 'anthesis' for trial Mur13W.
GWASDrops <- runSingleTraitGwas(gData = gDataDropsDedup,
                                trials = "Mur13W",
                                traits = c("grain.yield", "anthesis"))
print(head(GWASDrops$GWAResult$Mur13W), row.names = FALSE)
print(GWASDrops$signSnp$Mur13W, row.names = FALSE)

## GWAS Summary
## Create summary of GWASDrops.
summary(GWASDrops)

## QQ plot of GWAS Drops.
plot(GWASDrops, plotType = "qq", trait = "grain.yield")

## Manhattan plot of GWAS Drops.
plot(GWASDrops, plotType = "manhattan", trait = "grain.yield")
## Plot a manhattan plot of GWAS Drops.
## Set significance threshold to 4 and only plot chromosomes 6 to 8.
plot(GWASDrops, plotType = "manhattan", trait = "grain.yield", yThr = 4, chr = 6:8)

## Plot a manhattan plot of GWAS Drops.
## Plot only 5% of SNPs with a LOD below 3.
set.seed(1)
plot(GWASDrops, plotType = "manhattan", trait = "grain.yield", lod = 3)

#=======================================================================
#============================= QTL plots ===============================
#=======================================================================
## Qtl plot of GWAS Drops.
## Set significance threshold to 4 and normalize effect estimates.
plot(GWASDrops, plotType = "qtl", yThr = 4, normalize = TRUE)

## Plot a qtl plot of GWAS Drops for Mur13W.
plot(GWASDrops, plotType = "qtl")

## Plot a qtl plot of GWAS Drops for Mur13W.
## Set significance threshold to 4.
plot(GWASDrops, plotType = "qtl", yThr = 4)

## Plot a qtl plot of GWAS Drops for Mur13W.
## Set significance threshold to 4 and normalize effect estimates.
plot(GWASDrops, plotType = "qtl", yThr = 4, normalize = TRUE)

#=======================================================================
#========================== Kinship matrices ===========================
#=======================================================================
## Run single trait GWAS for trial 'Mur13W' and trait 'grain.yield'
## Use chromosome specific kinship matrices computed using method of van Raden.
GWASDropsChrSpec <- runSingleTraitGwas(gData = gDataDropsDedup, 
                                       traits = "grain.yield",
                                       trials = "Mur13W",
                                       GLSMethod = "multi",
                                       kinshipMethod = "vanRaden")

#=======================================================================
#========================== Further options ============================
#=======================================================================
## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Use a fixed significance threshold of 4.
GWASDropsFixThr <- runSingleTraitGwas(gData = gDataDropsDedup,
                                      trials = "Mur13W",
                                      traits = "grain.yield",
                                      thrType = "fixed",
                                      LODThr = 4)

## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Use the Newton Raphson algorithm for computing the variance components.
GWASDropsNR <- runSingleTraitGwas(gData = gDataDropsDedup,
                                  trials = "Mur13W",
                                  traits = "grain.yield",
                                  remlAlgo = "NR")

## Genomic control correction
GWASDrops$GWASInfo$inflationFactor$Mur13W
## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Perform genomic correction on the p-Values.
GWASDropsGenControl <- runSingleTraitGwas(gData = gDataDropsDedup,
                                          trials = "Mur13W",
                                          traits = "grain.yield",
                                          genomicControl = TRUE)

## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Use PZE-106021410, the most significant SNP, a SNP covariate.
GWASDropsSnpCov <- runSingleTraitGwas(gData = gDataDropsDedup,
                                      trials = "Mur13W",
                                      traits = "grain.yield",
                                      snpCov = "PZE-106021410")

## Minor Allele Frequency
## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Only include SNPs that have a MAC of at least 20
GWASDropsMAC <- runSingleTraitGwas(gData = gDataDropsDedup,
                                   trials = "Mur13W",
                                   traits = "grain.yield",
                                   useMAF = FALSE,
                                   MAC = 20)
## SNPs close to significant SNPs
## Run single trait GWAS for trait 'grain.yield' for Mur13W.
## Include SNPs within 200000 centimorgan of significant SNPs with a minimum LD of 0.1.
GWASDropsInclClose <- runSingleTraitGwas(gData = gDataDropsDedup,
                                         trials = "Mur13W",
                                         traits = "grain.yield",
                                         sizeInclRegion = 200000,
                                         minR2 = 0.1)
## Check signSnp in output.
print(head(GWASDropsInclClose$signSnp$Mur13W), row.names = FALSE)
