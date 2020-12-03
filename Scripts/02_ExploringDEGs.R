source("auxFunctions.R")
library(tidyverse)
library(pheatmap)

dir.create("../Results/data")
dir.create("../Results/01_DEgenes_vs_network")

#######################
t2gtair10 <- read.table("~/Desktop/Spring2020/scRNASeq/t2gtair10.txt",
                        header = F,col.names = c("","Gene","Name"),
                        sep = "\t",quote = "",row.names = 1,stringsAsFactors = F)
t2gtair10$Name <- ifelse(t2gtair10$Name=="",t2gtair10$Gene,t2gtair10$Name)
head(t2gtair10)

t2gtair10 <- t2gtair10[!duplicated(t2gtair10$Gene),]


###### Read expression data
normData_contrast <- get(load("../Results/bwt2_genome_contrast/combatData.RData"))
normData_contrast_meanLong <- MeanDataLong(normData_contrast)
normData_contrast_meanWide <- normData_contrast_meanLong %>% pivot_wider(id_cols = Gene,values_from = meanExpr,names_from = Group)
normData_contrast_meanWide <- data.frame(normData_contrast_meanWide[,-1],
                                         row.names = normData_contrast_meanWide$Gene)


save(normData_contrast_meanWide,file = "../Results/data/normData_contrast_meanWide.RData")
save(normData_contrast_meanLong,file = "../Results/data/normData_contrast_meanLong.RData")

###### Read DE data
bw2Genome_modelNO3 <- get(load("../Results/bwt2_genome/DEanalysis_anova_NO3.RData"))
bw2Genome_modelFull <- get(load("../Results/bwt2_genome/DEanalysis_anova.RData"))
bw2Genome_contrast <- get(load("../Results/bwt2_genome_contrast/DEanalysis_contrasts.RData"))

signifNO3 <- bw2Genome_modelNO3[bw2Genome_modelNO3$adj.P.Val < 0.05,]
signifModel <- bw2Genome_modelFull[bw2Genome_modelFull$adj.P.Val < 0.05,]
signifContrast <- bw2Genome_contrast[bw2Genome_contrast$adj.P.Val < 0.05,]

library(Vennerable)
vennList <- list("signifModel"=rownames(signifModel),
                 "signifNO3"=rownames(signifNO3),
                 "signifContrast"=rownames(signifContrast))
plot(Venn(vennList), doWeights = F,type="circles")


#### Network
YNM <- read.table("../meta/YNM_S03.txt",stringsAsFactors = F,header = T,sep="\t")
NitroGenes <- unique(c(YNM$TF_AGI,YNM$Promoter_AGI))
length(NitroGenes)

pdf("../Results/01_DEgenes_vs_network/DEvsNetwork.pdf",paper = 'a4r')
vennList <- list("signifModel"=rownames(signifModel),
                 "signifNO3"=rownames(signifNO3),
                 "signifContrast"=rownames(signifContrast),
                 "Gaudinier2018"=NitroGenes)
plot(Venn(vennList), doWeights = F,type="ellipses")


vennList <- list(
                 "signifNO3"=rownames(signifNO3),
                 "Gaudinier2018"=NitroGenes)
plot(Venn(vennList), doWeights = T,type="circles")
dev.off()

### Aggregating NxT genes

NO3Genes <- rownames(signifNO3)
modelGenes <- rownames(signifModel)
contrastGenes <- rownames(signifContrast)
allGenes <- unique(c(modelGenes,contrastGenes,NO3Genes))


vennList <- list("NxT"=allGenes,
                 "Gaudinier2018"=NitroGenes)
plot(Venn(vennList), doWeights = T,type="circles")


#################### Enrichment
Universe <- unique(rownames(bw2Genome_contrast))

notDE <- Universe[!Universe %in% allGenes]
DE <- allGenes

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

NxT_inNet <- genesInContingencyTable$DE_inNet


#############
library(RColorBrewer)
#hmPalette <- colorRampPalette(brewer.pal(name = "Blues",9))(500)
####

order<- c(
  #"t0min_1mM",
  "t15min_10mM","t45min_10mM","t90min_10mM","t180min_10mM","t1200min_10mM",
  "t15min_1mM","t45min_1mM","t90min_1mM","t180min_1mM","t1200min_1mM")
normData_contrast_meanWide <- normData_contrast_meanWide[,order]

hmData <- normData_contrast_meanWide[NxT_inNet,6:10]-normData_contrast_meanWide[NxT_inNet,1:5] 

hmPalette <- colorRampPalette(c('cyan','black','yellow'))(200)
paletteLength <- length(hmPalette)
boundPalette <- max(abs(range(hmData)))
myBreaks <- c(seq(-boundPalette, -0.0001, length.out=ceiling(paletteLength/2)), 
              seq(boundPalette/paletteLength, boundPalette, 
                  length.out=floor(paletteLength/2)))

##
hmData <- hmData %>% mutate(Gene=rownames(hmData)) %>% 
  left_join(t2gtair10,by='Gene') %>% select(-6)

###
wh <- 13

pdf("../Results/01_DEgenes_vs_network/DE_networkGenes.pdf",paper = 'a4r')
pheatmap(hmData[,-6], 
         border_color = NA, 
         #annotation_row = signifDPFmembership[,'cluster',drop=F],
         cluster_cols = F,  cluster_rows = T,labels_row = hmData$Name,
         fontsize_row = 6, fontsize_col = 8,
         cellheight = wh,cellwidth = wh,
         color = hmPalette,breaks = myBreaks)

makeLinePlots(normData_contrast_meanLong,NxT_inNet,nrow = 3)

####################################
### Breaking down the network & DEGs by TF or promoter
ynmTF <- unique(YNM$TF_AGI)
ynmPromoter <- unique(YNM$Promoter_AGI)

vennList <- list("NxT"=NxT_inNet,
                 "TFs"=ynmTF,
                 "Promoters"=ynmPromoter)
plot.new()
plot(Venn(vennList), doWeights = F,type="circles")
dev.off()

pdf("../Results/01_DEgenes_vs_network/DE_networkGenes_TF_Promoters.pdf",paper = 'a4r')
NxT_tfs <- intersect(NxT_inNet,unique(YNM$TF_AGI))
length(NxT_tfs)
makeLinePlots(normData_contrast_meanLong,NxT_tfs,nrow = 2,title = "TFs")

NxT_prom <- intersect(NxT_inNet,unique(YNM$Promoter_AGI))
length(NxT_prom)
makeLinePlots(normData_contrast_meanLong,NxT_prom,nrow = 2,title = "Promoters")
dev.off()
####


#### Save this
nxt_networkType <- rbind( data.frame('Gene'=NxT_inNet[NxT_inNet %in% ynmTF],'Type'='TF'),
       data.frame('Gene'=NxT_inNet[NxT_inNet %in% ynmPromoter],'Type'='Target'))

nxt_networkType <- nxt_networkType %>% left_join(t2gtair10,by='Gene')
write.table(file='../Results/01_DEgenes_vs_network/DEvsNetwork.txt',
            x=nxt_networkType,quote = F,sep = '\t',row.names = F)


#####################
nxt_network <- YNM[YNM$TF_AGI %in% NxT_inNet | YNM$Promoter_AGI %in% NxT_inNet,]
###
 # Save this, send to cytoscape
write.table(file='../Results/01_DEgenes_vs_network/YNM_subset.txt',
            x=nxt_network,quote = F,sep = '\t',row.names = F)

# 
# outdegree
# for (tmpTF in intersect(ynmTF,NxT_inNet)){
#   print(tmpTF)
#   tmpNet <- YNM[YNM$TF_AGI %in% tmpTF,]
#   outdegree <- c(outdegree,nrow(tmpNet))
#   gens <- c(tmpTF,tmpNet$Promoter_AGI)
#   targets <- c(tmpNet$Promoter_AGI)
#   
#   makeLinePlots(normData_contrast_meanLong,gens,title = tmpTF)
#   if (nrow(tmpNet) > 1) makeHeatMap(normData_contrast_meanWide,gens,order,hmPalette,show_rownames = T)
# }







hmData <- normData_contrast_meanWide[rownames(signifNO3),order]
hmData <- hmData[,6:10]-hmData[,1:5] 

hmPalette <- colorRampPalette(c('cyan','black','yellow'))(200)
paletteLength <- length(hmPalette)
boundPalette <- max(abs(range(hmData)))
myBreaks <- c(seq(-boundPalette, -0.0001, length.out=ceiling(paletteLength/2)), 
              seq(boundPalette/paletteLength, boundPalette, 
                  length.out=floor(paletteLength/2)))

##
hmData <- hmData %>% mutate(Gene=rownames(hmData)) %>% 
  left_join(t2gtair10,by='Gene') %>% select(-6)

###
wh <- 2
pdf("../Results/01_DEgenes_vs_network/DE_significant.pdf",height = 20,width = 10)
pheatmap(hmData[,-6], 
         border_color = NA, 
         #annotation_row = signifDPFmembership[,'cluster',drop=F],
         cluster_cols = F,  cluster_rows = T,labels_row = hmData$Name,
         fontsize_row = 2, fontsize_col = 8, #cutree_rows = 5,
         cellheight = wh,cellwidth = 8, treeheight_col = 200,
         color = hmPalette,breaks = myBreaks)

pheatmap(hmData[,-6], scale = 'row',
         border_color = NA, 
         #annotation_row = signifDPFmembership[,'cluster',drop=F],
         cluster_cols = F,  cluster_rows = T,labels_row = hmData$Name,
         fontsize_row = 2, fontsize_col = 8, #cutree_rows = 5,
         cellheight = wh,cellwidth = 8, treeheight_col = 200,
         color = hmPalette,breaks = myBreaks)

hmPalette <- colorRampPalette(c('cyan','white','yellow'))(200)
pheatmap(hmData[,-6], 
         border_color = NA, 
         #annotation_row = signifDPFmembership[,'cluster',drop=F],
         cluster_cols = F,  cluster_rows = T,labels_row = hmData$Name,
         fontsize_row = 2, fontsize_col = 8, #cutree_rows = 5,
         cellheight = wh,cellwidth = 8, treeheight_col = 200,
         color = hmPalette,breaks = myBreaks)


hmPalette <- colorRampPalette(c('cyan','gray80','yellow'))(200)
pheatmap(hmData[,-6], 
         border_color = NA, 
         #annotation_row = signifDPFmembership[,'cluster',drop=F],
         cluster_cols = F,  cluster_rows = T,labels_row = hmData$Name,
         fontsize_row = 2, fontsize_col = 8, #cutree_rows = 5,
         cellheight = wh,cellwidth = 8, treeheight_col = 200,
         color = hmPalette,breaks = myBreaks)

hmPalette <- colorRampPalette(RColorBrewer::brewer.pal(n=7,"PuOr"))(200)
pheatmap(hmData[,-6], scale = 'row',
         border_color = NA, 
         #annotation_row = signifDPFmembership[,'cluster',drop=F],
         cluster_cols = F,  cluster_rows = T,labels_row = hmData$Name,
         fontsize_row = 2, fontsize_col = 8, #cutree_rows = 5,
         cellheight = wh,cellwidth = 8, treeheight_col = 200,
         color = hmPalette,breaks = myBreaks)
dev.off()
