library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

source("~/Desktop/Nitr_timeSeries/TimeSeries/Scripts/auxFunctions.R")
source("~/Desktop/BioTools/RBioFunctions/metaFunctions_forNetworkAnalysis.R")


####################################################
### Enrichment of Network
YNM <- read.table("../meta/YNM_S03.txt",stringsAsFactors = F,header = T,sep="\t")
NitroGenes <- unique(c(YNM$TF_AGI,YNM$Promoter_AGI))
length(NitroGenes)

### Expression data
normLong <- get(load("../Results/data/normData_contrast_meanLong.RData"))
normWide <- get(load("../Results/data/normData_contrast_meanWide.RData"))

### Arabidopsis annotation
t2gtair10 <- read.table("~/Desktop/Spring2020/scRNASeq/t2gtair10.txt",
                        header = F,col.names = c("","Gene","Name"),
                        sep = "\t",quote = "",row.names = 1,stringsAsFactors = F)
t2gtair10$Name <- ifelse(t2gtair10$Name=="",t2gtair10$Gene,t2gtair10$Name)
head(t2gtair10)

t2gtair10 <- t2gtair10[!duplicated(t2gtair10$Gene),]


### Universe is the genes used in the contrast (full)set
# all(rownames(bw2Genome_modelFull) %in% rownames(bw2Genome_contrast)) The full model (sans time 0) is contained within this set
Universe <- unique(rownames(normWide))



order<- c(
  #"t0min_1mM",
  "t15min_10mM","t45min_10mM","t90min_10mM","t180min_10mM","t1200min_10mM",
  "t15min_1mM","t45min_1mM","t90min_1mM","t180min_1mM","t1200min_1mM")

dirsK10 <- list.files("../DPF/a_combat/",pattern = "NitrTimeSeries_k10",full.names = T)

fileName <- gsub("NitrTimeSeries_",'',basename(dirsK10))
kMembFile <- list.files(dirsK10,pattern = "DPFgenesToPatterns.txt",full.names = T)  
Clustresults <- read.delim(kMembFile,sep = "\t",as.is = T)

# outPDF <- paste0("Summary_",fileName,'.pdf')
# pdf(outPDF,paper = "a4")

####################################################
Clustresults_lst <- lapply(Clustresults, function(x){x[x!='']})

## Get average cluster expression
GenesXClusterDF <- data.frame(do.call("rbind",lapply(names(Clustresults_lst), 
                                                     function(x){cbind("Gene"=Clustresults_lst[[x]],"cluster"=x)})),stringsAsFactors = F)

exprData <- normLong %>% left_join(GenesXClusterDF,by = "Gene") %>% filter(!is.na(cluster))

exprDataXDPF <- exprData %>% group_by(Group,cluster,Time,NO3) %>% 
  summarise("meanXcluster"=mean(meanExpr)) %>% ungroup() 

exprDataXDPF <- exprDataXDPF %>% 
  mutate(cluster=factor(exprDataXDPF$cluster,
                        levels = gtools::mixedsort(unique(exprDataXDPF$cluster))))

GenesXCluster <- sapply(Clustresults_lst,length)

#### Genes per cluster
GxP <- data.frame("Pattern"=gsub("Pattern_","",names(GenesXCluster)),
                  "n"=GenesXCluster) %>%
  ggplot(aes(x=factor(gtools::mixedorder(Pattern,decreasing = T)),y=n)) +
  geom_bar(stat="identity",width = 0.75,fill="skyblue") +
  geom_text(aes(x=factor(gtools::mixedorder(Pattern,decreasing = T)),y=median(n),label=n),size=3,color="black") +
  xlab("Cluster") +
  ylab("Genes assigned") +
  coord_flip() +
  theme_light()

print(GxP)

#### Avg cluster expression
####################################################
linPlot <- exprDataXDPF %>% mutate(Time=as.numeric(gsub("t|min","",exprDataXDPF$Time))) %>%
  #levels = gtools::mixedsort( unique(gsub("t|min","",exprDataXDPF$Time)) 
  #))) %>%
  ggplot(aes(x=Time,y=meanXcluster)) +
  geom_point(aes(color=NO3,group=NO3,shape=NO3) ) +
  geom_line(aes(color=NO3,group=NO3) ) +
  theme(aspect.ratio=16/9) +
  facet_wrap(as.factor(cluster)~.,scales = "free_x",nrow = 3) +
  theme_classic()
print(linPlot)

meanXcluster <- exprDataXDPF %>% mutate(Time=as.numeric(gsub("t|min","",exprDataXDPF$Time)))  %>% 
  pivot_wider(names_from = Group,values_from = meanXcluster,id_cols = cluster) %>% data.frame(.,row.names = .$cluster)
meanXcluster <- meanXcluster[,-1]

# Reorder & get average change per timepoint
meanXcluster <- meanXcluster[,order]
hmData <- meanXcluster[6:10]-meanXcluster[1:5] 

#hmPalette <- colorRampPalette(brewer.pal(n=7,'PuOr'))(20)
hmPalette <- colorRampPalette(c('cyan','black','yellow'))(200)
paletteLength <- length(hmPalette)
boundPalette <- max(abs(range(hmData)))
myBreaks <- c(seq(-boundPalette, -0.0001, length.out=ceiling(paletteLength/2)), 
              seq(boundPalette/paletteLength, boundPalette, 
                  length.out=floor(paletteLength/2)))


wh <- 20
pheatmap(hmData, 
         border_color = NA, 
         #annotation_row = signifDPFmembership[,'cluster',drop=F],
         cluster_cols = F,  cluster_rows = T,
         fontsize_row = 6, fontsize_col = 8,
         cellheight = wh,cellwidth = wh,
         color = hmPalette,breaks = myBreaks)

#### FET x Cluster
####################################################
FET_cluster <- lapply(Clustresults_lst,function(x){
  inCluster <- x
  notCluster <- Universe[!Universe %in% inCluster]
  
  clust_inNet <- inCluster[inCluster %in% NitroGenes]
  clust_notNet <- inCluster[!inCluster %in% NitroGenes]
  notClust_inNet <- notCluster[notCluster %in% NitroGenes]
  notClust_notNet <- notCluster[!notCluster %in% NitroGenes]
  
  sum(
    length(clust_inNet),
    length(clust_notNet),
    length(notClust_inNet),
    length(notClust_notNet)
  ) == length(Universe)
  
  
  
  tmpContingency <- matrix( c(length(clust_inNet),
                              length(clust_notNet),
                              length(notClust_inNet),
                              length(notClust_notNet)),nrow=2,ncol=2)
  
  
  
  tmpMembers <- condenseGeneList_toMatrix(list(
    "clust_inNet"=clust_inNet,
    "clust_notNet"=clust_notNet,
    "notClust_inNet"=notClust_inNet,
    "notClust_notNet"=notClust_notNet))
  
  tmpFET <- fisher.test(tmpContingency,
                        alternative="greater")    
  returnThis <- list("contingency"=tmpContingency,"FET"=tmpFET,"contMatMembers"=tmpMembers)
  return(returnThis)
})

## Extract 
contingencyMatrices <- lapply(FET_cluster,"[[",1)
contMatMembers <- lapply(FET_cluster,"[[",3)

# Correct p value
FETresults <- do.call("rbind",lapply(FET_cluster,function(x){
  data.frame("OR"=x[["FET"]]$estimate,
             "pVal"=x[["FET"]]$p.value)}))
FETresults$FDR <- p.adjust(FETresults$pVal,method = "fdr")


stars <- ifelse(FETresults$FDR < 0.01,"*","")
hmData <- FETresults[,1,drop=F]
wh <- 14
hmPalette <- colorRampPalette(brewer.pal(n=7,"Blues"))(20)
pheatmap(hmData,border_color = NA, display_numbers = cbind(stars), number_color = "white",
         cluster_cols = F,  cluster_rows = F,
         #fontsize_row = 6, fontsize_col = 8,
         cellheight = wh,cellwidth = wh, main='enrichment of clusters in the network',
         color = hmPalette)


############################################
#### Network genes in these clusters
signDPFMembers <- lapply(contMatMembers[FETresults$FDR < 0.01],
                         function(x){tmp=x[,1];
                         tmp[tmp != ''] })

genesNetworkSignificantClusters <- Reduce(union,signDPFMembers)

signifDPFmembership <- do.call('rbind',lapply(names(signDPFMembers),function(x){
  data.frame(cbind('cluster'=x,'gene'=signDPFMembers[[x]]),
             row.names = signDPFMembers[[x]],stringsAsFactors = F)
}))

commonGenes <- intersect(rownames(normWide),rownames(signifDPFmembership))
signifDPFmembership <- signifDPFmembership[commonGenes,] #Subset by genes in table and network

# Order by cluster
signifDPFmembership <- signifDPFmembership[gtools::mixedorder(signifDPFmembership$cluster),] 
hmData <- normWide[rownames(signifDPFmembership),]

# Get fold change
hmData <- hmData[,order]
hmData_FC <- hmData[6:10]-hmData[1:5]
head(hmData_FC)

# Get gene aliases
geneNames <- t2gtair10 %>% filter(Gene %in% rownames(hmData_FC)) %>% data.frame(.,row.names = .$Gene)
geneNames <- geneNames[rownames(hmData_FC),] #Reorder

clustMemb <- table(signifDPFmembership$cluster)
clustMemb <- clustMemb[gtools::mixedsort(names(clustMemb))]

hmPalette <- colorRampPalette(c('cyan','black','yellow'))(200)
paletteLength <- length(hmPalette)
boundPalette <- max(abs(range(hmData_FC)))
myBreaks <- c(seq(-boundPalette, -0.0001, length.out=ceiling(paletteLength/2)), 
              seq(boundPalette/paletteLength, boundPalette, 
                  length.out=floor(paletteLength/2)))

title <- paste0("Total of network genes: ", nrow(hmData_FC))
wh <- 5
pheatmap(hmData_FC, labels_row = geneNames$Name,
         border_color = NA, gaps_row = cumsum(clustMemb),
         annotation_row = signifDPFmembership[,'cluster',drop=F],
         cluster_cols = F,  cluster_rows = F,
         fontsize_row = 6, fontsize_col = 8,
         cellheight = wh,cellwidth = wh, main=title,
         color = hmPalette,breaks = myBreaks)
dev.off()

#### Check genes in clusters vs DE genes


  
# exprData %>% filter(Gene %in% genesNetworkSignificantClusters) %>%
#   mutate(Time=factor(gsub("t|min","",Time),levels = gtools::mixedsort( unique(gsub("t|min","",Time))))) %>%
#   #mutate(Time=as.numeric(gsub("t|min","",Time))) %>%
#   mutate("tmp"=paste0(Gene," (",Name,")")) %>%
#   ggplot(aes(x=Time,y=meanExpr)) +
#   geom_point(aes(color=NO3,group=Gene,shape=NO3) ) +
#   geom_line(aes(color=NO3,group=NO3) ) +
#   theme(aspect.ratio=16/9) +
#   #ggtitle(title) +
#   facet_wrap(as.factor(tmp)~cluster,scales = "free_x",nrow = 4) +
#   theme_classic()
