library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(Vennerable)
source("helper_exploreDEGs.R")
source("auxFunctions.R")
source('~/Desktop/BioTools/RBioFunctions/metaFunctions_forNetworkAnalysis.R')

#### Check the expression of those genes
normWide <- get(load("../Results/data/normData_contrast_meanWide.RData"))
normLong<- get(load("../Results/data/normData_contrast_meanLong.RData"))
###############################
bw2Genome_modelNO3 <- get(load("../Results/bwt2_genome/DEanalysis_anova_NO3.RData"))
signifNO3 <- bw2Genome_modelNO3[bw2Genome_modelNO3$adj.P.Val < 0.05,]
dim(signifNO3)


### Enrichment of Network
YNM <- read.table("../meta/YNM_S03.txt",stringsAsFactors = F,header = T,sep="\t")
NitroGenes <- unique(c(YNM$TF_AGI,YNM$Promoter_AGI))
length(NitroGenes)

#################################################################
#NxG <- read.csv("../meta/Gaudinier2018/SingleMutants_logFC - shoots_singlemutants.csv",row.names = 1,header = T)
NxG <- read.delim("../meta/Gaudinier2018/SupplementalTable15_NxG.txt",row.names = 1,header = T,sep = "\t")
dim(NxG)

head(NxG)

NxG_signif <- NxG[,grep('wt1mM_vs_wt10mM.adj.P.Val',colnames(NxG)),drop=F]
NxG_signif <- NxG_signif[apply(NxG_signif,1,function(x){any(x < 0.01)}),,drop=F]
dim(NxG_signif)

NxG_signif <- NxG[rownames(NxG_signif),c('wt1mM_vs_wt10mM.logFC','InNetwork','SymbolsColumn')]
Gaudinier2018 <- rownames(NxG_signif)
#################################################################


####################################
############
library(Vennerable)

############ YNM vs krouk vs significant
vennList <- list("YNM"=NitroGenes,
                 "signifNO3"=rownames(signifNO3),
                 "Gaudinier2018_1vs10mM"=Gaudinier2018)

outFile <- makePath('Gaudinier2018_1vs10mM_vs_YNM','pdf')
pdf(outFile,paper = "a4r")
plot(Venn(vennList), doWeights = T,type="circles")
plot(Venn(vennList), doWeights = F,type="circles")
dev.off()


#####
#################### Enrichment
Universe <- unique(rownames(NxG))

notDE <- Universe[!Universe %in% Gaudinier2018]
DE <- Gaudinier2018

#tst <- intersect(rownames(signifNO3),NitroGenes)
tst <- intersect(intersect(NitroGenes,Gaudinier2018),rownames(signifNO3))
DE_inNet <- DE[DE %in% tst]
DE_notNet <- DE[!DE %in% tst]
notDE_inNet <- notDE[notDE %in% tst]
notDE_notNet <- notDE[!notDE %in% tst]

sum(
  length(DE_inNet),
  length(DE_notNet),
  length(notDE_inNet),
  length(notDE_notNet)
) == length(Universe)

contTable <- matrix( c(length(DE_inNet),
                       length(DE_notNet),
                       length(notDE_inNet),
                       length(notDE_notNet)),nrow=2,ncol=2)

fisher.test(contTable,
            alternative="greater")


genesInContingencyTable <- list(
  "DE_inNet"=DE_inNet,
  "DE_notNet"=DE_notNet,
  "notDE_inNet"=notDE_inNet,
  "notDE_notNet"=notDE_notNet)



Gaud_NxT <- condenseGeneList_toMatrix(list('Gaudinier_NxT'=DE_inNet,
'Gaudinier_YNM'=intersect(NitroGenes,Gaudinier2018),
'Gaudinier_YNM_NxT'=intersect(intersect(NitroGenes,Gaudinier2018),rownames(signifNO3))))

View(Gaud_NxT)



