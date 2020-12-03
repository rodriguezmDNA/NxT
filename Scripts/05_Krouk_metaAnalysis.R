library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(Vennerable)
source("helper_exploreDEGs.R")
source("auxFunctions.R")


#####
dir.create("../Results/04_Krouk",showWarnings = F)

outPattern <- paste0('','')
makePath <- function (title,ext){paste0("../Results/04_Krouk/",title,'.',ext)}
makePath('','')


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


####################################
logFC_krouk2010 <- read.table("https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2010-11-12-r123/MediaObjects/13059_2010_2499_MOESM8_ESM.txt",sep = "\t",header = T,row.names = 1,
                              stringsAsFactors = F, skip=1,
                              col.names = c("Gene","min03","min06","min09","min12","min15","min20"))
rownames(logFC_krouk2010) <- toupper(rownames(logFC_krouk2010))
krouk2010 <- rownames(logFC_krouk2010)
length(krouk2010)

####################################
############
library(Vennerable)

############ YNM vs krouk vs significant
vennList <- list("YNM"=NitroGenes,
                 "signifNO3"=rownames(signifNO3),
                 "Krouk2010"=krouk2010)

outFile <- makePath('Krouk2010_vs_YNM','pdf')
pdf(outFile,paper = "a4r")
plot(Venn(vennList), doWeights = T,type="circles")
plot(Venn(vennList), doWeights = F,type="circles")
dev.off()


#### Enrichment 

#################### Enrichment
Universe <- unique(rownames(bw2Genome_modelNO3))

notDE <- Universe[!Universe %in% krouk2010]
DE <- krouk2010

DE_inNet <- DE[DE %in% NitroGenes]
DE_notNet <- DE[!DE %in% NitroGenes]
notDE_inNet <- notDE[notDE %in% NitroGenes]
notDE_notNet <- notDE[!notDE %in% NitroGenes]

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

View(tibble('Gene'=DE_inNet) %>% left_join(t2gtair10,by='Gene'))
NxT_inNet <- genesInContingencyTable$DE_inNet













################# Expression of common genes between Krouk and YNM
KroukYNM <- intersect(NitroGenes,krouk2010)
KroukYNM <- setdiff(KroukYNM,rownames(signifNO3))

outFile <- makePath('Krouk2010_overlap_YNM','pdf')
pdf(outFile,width = 20,height = 15)

hmData <- normWide[KroukYNM,order]

hmData <- hmData[6:10]-hmData[1:5] 

hmPalette <- colorRampPalette(c('cyan','black','yellow'))(200)
paletteLength <- length(hmPalette)
boundPalette <- max(abs(range(hmData)))
myBreaks <- c(seq(-boundPalette, -0.0001, length.out=ceiling(paletteLength/2)), 
              seq(boundPalette/paletteLength, boundPalette, 
                  length.out=floor(paletteLength/2)))

hmData <- hmData %>% mutate(Gene=rownames(hmData)) %>% 
  left_join(t2gtair10,by='Gene') %>% select(- (ncol(hmData) + 1))


p01 <- ifelse(bw2Genome_modelNO3[KroukYNM,9] < .1,"*","")
hmData$Name <- paste(hmData$Name,p01)

wh <- 8
pheatmap(hmData[,- (ncol(hmData))], 
         border_color = NA, 
         #annotation_row = signifDPFmembership[,'cluster',drop=F],
         cluster_cols = F,  cluster_rows = T,labels_row = hmData$Name,
         fontsize_row = 6, fontsize_col = 8,
         cellheight = wh,cellwidth = wh,
         color = hmPalette,breaks = myBreaks)

pVal <- KroukYNM[bw2Genome_modelNO3[KroukYNM,9] < .1]
makeLinePlots(normLong,pVal,title = 'YNM overlap Krouk p<0.1')

pVal <- KroukYNM[!bw2Genome_modelNO3[KroukYNM,9] < .1]
makeLinePlots(normLong,pVal,title = 'YNM overlap Krouk - p>0.1')
dev.off()



bw2Genome_modelNO3 %>%
  ggplot(aes(x=F,y=adj.P.Val)) + 
  geom_point(size=0.5) + 
  ###
  geom_hline(yintercept = .1,color='skyblue') +
  geom_hline(yintercept = .2,color='orange') +
  geom_vline(xintercept = min(bw2Genome_modelNO3[bw2Genome_modelNO3$adj.P.Val < 0.1,'F']),color='skyblue') +
  geom_vline(xintercept = min(bw2Genome_modelNO3[bw2Genome_modelNO3$adj.P.Val < 0.2,'F']),color='orange') +
  theme_light()


signifNO3_p1 <- bw2Genome_modelNO3[bw2Genome_modelNO3$adj.P.Val < 0.1,]
dim(signifNO3_p1)
vennList <- list("YNM"=NitroGenes,
                 "signifNO3"=rownames(signifNO3_p1),
                 "Krouk2010"=krouk2010)
plot(Venn(vennList), doWeights = T,type="circles")
plot(Venn(vennList), doWeights = F,type="circles")


p1_YNM <- intersect(rownames(signifNO3_p1),NitroGenes)
makeLinePlots(normLong,p1_YNM,title = 'p < 0.1',nrow = )


# #### Other 
# ################################################

NitroGenes_expression <- bw2Genome_modelNO3[intersect(rownames(bw2Genome_modelNO3),NitroGenes),]
dim(NitroGenes_expression)
head(NitroGenes_expression)

NitroGenes_expression %>%
  ggplot(aes(x=F,y=adj.P.Val)) +
  geom_point(size=0.5) +
  ###
  geom_hline(yintercept = .1,color='skyblue') +
  geom_hline(yintercept = .2,color='orange') +
  geom_vline(xintercept = min(NitroGenes_expression[NitroGenes_expression$adj.P.Val < 0.1,'F']),color='skyblue') +
  geom_vline(xintercept = min(NitroGenes_expression[NitroGenes_expression$adj.P.Val < 0.2,'F']),color='orange') +
  theme_light()

NitroGenes_expression_p2 <- NitroGenes_expression[NitroGenes_expression$adj.P.Val < 0.2,]
dim(NitroGenes_expression_p2)

vennList <- list("YNM"=NitroGenes,
                 "krouk2010"=krouk2010,
                 "DEGs"=rownames(NitroGenes_expression_p2))

# outFile <- makePath('Krouk2010_vs_YNM','pdf')
# pdf(outFile,paper = "a4r")
plot(Venn(vennList), doWeights = T,type="circles")
plot(Venn(vennList), doWeights = F,type="circles")
#dev.off()





 
# #### Other 
# ################################################
# kroukOnly <- setdiff(krouk2010,union(NitroGenes,rownames(signifNO3)))
# kroukOnly_expression <- bw2Genome_modelNO3[intersect(rownames(bw2Genome_modelNO3),kroukOnly),]
# dim(kroukOnly_expression)
# head(kroukOnly_expression)
# 
# kroukOnly_expression %>%
#   ggplot(aes(x=F,y=adj.P.Val)) + 
#   geom_point(size=0.5) + 
#   ###
#   geom_hline(yintercept = .1,color='skyblue') +
#   geom_hline(yintercept = .2,color='orange') +
#   geom_vline(xintercept = min(kroukOnly_expression[kroukOnly_expression$adj.P.Val < 0.1,'F']),color='skyblue') +
#   geom_vline(xintercept = min(kroukOnly_expression[kroukOnly_expression$adj.P.Val < 0.2,'F']),color='orange') +
#   theme_light()
# 
# kroukOnly_expression_p2 <- kroukOnly_expression[kroukOnly_expression$adj.P.Val < 0.2,]
# 
# 
# 
# vennList <- list("YNM"=NitroGenes,
#                  "signifNO3"=rownames(signifNO3),
#                  "Krouk2010"=rownames(kroukOnly_expression_p2))
# 
# # outFile <- makePath('Krouk2010_vs_YNM','pdf')
# # pdf(outFile,paper = "a4r")
# plot(Venn(vennList), doWeights = T,type="circles")
# plot(Venn(vennList), doWeights = F,type="circles")
# #dev.off()
