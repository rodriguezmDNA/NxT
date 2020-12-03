library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

source("helper_exploreDEGs.R")
source("auxFunctions.R")
#######################
chooseK <- 1

### Enrichment of Network
YNM <- read.table("../meta/YNM_S03.txt",stringsAsFactors = F,header = T,sep="\t")
NitroGenes <- unique(c(YNM$TF_AGI,YNM$Promoter_AGI))
length(NitroGenes)


#####
dir.create("../Results/03b_motifEnrichment")

outPattern <- paste0('Cluster_',chooseK)
makePath <- function (title,ext){paste0("../Results/03b_motifEnrichment/",title,"_",outPattern,'.',ext)}

kDir <- list.files("../Results/02b_DEgenes_ClusterAnalysis/k_4/AMEoutput/",pattern = paste0("_",chooseK),full.names = T)
#######################
ameResults <- read.delim(list.files(kDir,pattern = 'ame.tsv',full.names = T),comment.char = "#",stringsAsFactors = F)
ameSeqs <- read.delim(list.files(kDir,pattern = 'sequences.tsv',full.names = T),comment.char = "#",stringsAsFactors = F)

head(ameResults)
dim(ameResults)
all(ameResults$adj_p.value < 0.05) # All are significant after multiple testing correction


################ Get the predicted TFs and targets
motif2seq <- split(ameSeqs,ameSeqs$motif_ID)
motif2seq <- lapply(motif2seq,function(x){
  x %>% filter(class == "tp") %>% select(motif_ID,seq_ID)  
})
motif2seq <- do.call('rbind',motif2seq)

## Dictionary to add motif names to IDs
motifNames <- toupper(ameResults$motif_alt_ID)
names(motifNames) <- ameResults$motif_ID
##########################################

# Adjacency table for network building
motif2seq$motif_alt_ID <- motifNames[motif2seq$motif_ID]
rownames(motif2seq) <- NULL
head(motif2seq)
motif2seq <- motif2seq %>% select(motif_alt_ID,seq_ID) %>% rename('target'='seq_ID')


######

motif2seq <- motif2seq %>% 
              ### Add names to the genes
              left_join(t2gtair10,by = c("target"="Gene")) %>% rename('target_name'='Name') %>%
              # and to the motifs, if possible 
              left_join(t2gtair10,by = c("motif_alt_ID"="Gene")) %>% 
              rename('motif_name'='Name') %>%
              mutate(motif_name = ifelse(is.na(motif_name),motif_alt_ID,motif_name))

######

load("../meta/AGI_Dictionary/AGI_list_uniq.RData")
load("../meta/AGI_Dictionary/AGI_table_uniq.RData")
source('../meta/AGI_Dictionary/AthGenNoms.R')

getAGI <- data.frame('name'=unique(motif2seq$motif_name[grep("^[^AT.G.*$]",motif2seq$motif_name)]),
                     'AGI'='',
                     stringsAsFactors = F)

for (i in seq(nrow(getAGI))){
  print(i)
  tmp <- gene_to_AGI(getAGI[i,"name"])
  if (!is.null(tmp)){ getAGI$AGI[i] <- rownames(tmp) }

}

fillAGI <- getAGI[getAGI$AGI == '',]
fillAGI <- fillAGI %>% left_join(t2gtair10, by=c('name'='Name'))

getAGI <- getAGI %>% left_join(fillAGI) %>% mutate('AGI'=ifelse(is.na(Gene),AGI,Gene)) %>% select(-c(Gene))
head(getAGI)

getAGI[6,2] <- 'AT1G43160' #gene_to_AGI('RAP2')
getAGI[17,2] <- 'AT4G11080' #gene_to_AGI('BOX')


motif2seq <- motif2seq %>% left_join(getAGI,by=c('motif_name'='name')) %>% 
  mutate('AGI'=ifelse(is.na(AGI),motif_name,AGI)) %>%
  rename('motif_AGI'='AGI')

unique(motif2seq$motif_AGI)
motif2seq[motif2seq$motif_AGI == "AGL13",'motif_AGI'] <- 'AT3G61120'  #gene_to_AGI("AGL13")
motif2seq[motif2seq$motif_AGI == "TF3A",'motif_AGI'] <- 'AT1G72050' #gene_to_AGI("TF3")


outFile <- makePath('motif2seq_adjacencytable','txt')
write.table(file=outFile,motif2seq,quote = F,row.names = F,sep = "\t")


nodeAttribute <- data.frame(rbind(cbind('node'=motif2seq$target_name,'AGI'=motif2seq$target,'type'='target'),
                                  cbind('node'=motif2seq$motif_name,'AGI'=motif2seq$motif_AGI,'type'='TF')),stringsAsFactors = F)
nodeAttribute$YNM <- ifelse(nodeAttribute$AGI %in% NitroGenes,"YNM","")

outFile <- makePath('motif2seq_attributeNode','txt')
write.table(file=outFile,nodeAttribute,quote = F,row.names = F,sep = "\t")

############################################################
##### degree in/out

TF_outdegree <- sapply(split(motif2seq,motif2seq$motif_AGI),nrow)
TF_outdegree <- sort(TF_outdegree,decreasing = T)

##
targets_indegree <- sapply(split(motif2seq,motif2seq$target),nrow)
targets_indegree <- sort(targets_indegree,decreasing = T)



############################################################
targets <- unique(motif2seq$target)
predictedTFs <- unique(motif2seq$motif_AGI)



#to check the fraction of the genes that had a significant motif upstream
#clust3 <- read.delim("../Results/DEgenes_ClusterAnalysis/k_4/cluster_3.txt",header = F)
library(Vennerable)
vennList <- list("targets"=targets,
                 "predictedTFs"=predictedTFs,
                 #"clust3"=clust3$V1, #Fraction of the genes that had a significant motif upstream
                 "YNM"=NitroGenes)
outFile <- makePath('TF_target_YNM','pdf')
pdf(outFile,paper = "a4r")
plot(Venn(vennList), doWeights = F,type="circles")



outedgreeTable <- tibble('TF'=names(TF_outdegree),'outdegree'=TF_outdegree) %>% left_join(t2gtair10,by = c("TF"="Gene")) %>% 
  mutate(Name = fct_reorder(Name, (outdegree))) %>%
  mutate('YNM'=ifelse(TF %in% NitroGenes,"*",""))

inedgreeTable <- tibble('target'=names(targets_indegree),'indegree'=targets_indegree) %>% left_join(t2gtair10,by = c("target"="Gene")) %>% 
  mutate(Name = fct_reorder(Name, (indegree)))


####
print(
outedgreeTable %>%
  ggplot(aes(x=Name,y=outdegree)) + 
  geom_bar(stat="identity",width = 0.75,fill="skyblue") +
  geom_text(aes(x=Name,y=median(outdegree),label=outdegree),size=3,color="black") +
    geom_text(aes(x=Name,y=-0.65,label=YNM),size=6,color="orange",nudge_x = -0.25) +
  xlab("TF") +
  ylab("outdegree") +
  coord_flip() +
  theme_light()
)

print(
inedgreeTable %>%
  ggplot(aes(x=Name,y=indegree)) +
  geom_bar(stat="identity",width = 0.75,fill="skyblue") +
  geom_text(aes(x=Name,y=median(indegree),label=indegree),size=3,color="black") +
  xlab("target") +
  ylab("indegree") +
  coord_flip() +
  theme_light()
)
dev.off()

###############################
bw2Genome_modelNO3 <- get(load("../Results/bwt2_genome/DEanalysis_anova_NO3.RData"))
signifNO3 <- bw2Genome_modelNO3[bw2Genome_modelNO3$adj.P.Val < 0.05,]
dim(signifNO3)

### 
predictedTF_vs_DEGs <- bw2Genome_modelNO3[intersect(rownames(bw2Genome_modelNO3), 
                                                    predictedTFs),]

## 23 out of 28 were detected, but none significantly DE after multiple test correction
nrow(predictedTF_vs_DEGs) 
length(predictedTFs)



#### Check the expression of those genes
normWide <- get(load("../Results/data/normData_contrast_meanWide.RData"))
normLong<- get(load("../Results/data/normData_contrast_meanLong.RData"))

########################################################################################
############################################ Check expression of predicted TFs
common <- intersect(rownames(normWide), predictedTFs)


normWide_predictedTFs <- normWide[,order]
normWide_predictedTFs <- normWide_predictedTFs[names(TF_outdegree[names(TF_outdegree) %in% common]),] #Arrange by outdegree

normWide_predictedTFs <- normWide_predictedTFs[rowMeans(normWide_predictedTFs[,-11]) > 2,]
predictedTFs_expessed <- rownames(normWide_predictedTFs)
####
normWide_predictedTFs <- normWide_predictedTFs %>% mutate(Gene=rownames(normWide_predictedTFs)) %>% 
  left_join(t2gtair10,by='Gene') %>% select(-11)

hmData <- normWide_predictedTFs[6:10]-normWide_predictedTFs[1:5] 

hmPalette <- colorRampPalette(c('cyan','black','yellow'))(200)
paletteLength <- length(hmPalette)
boundPalette <- max(abs(range(hmData)))
myBreaks <- c(seq(-boundPalette, -0.0001, length.out=ceiling(paletteLength/2)), 
              seq(boundPalette/paletteLength, boundPalette, 
                  length.out=floor(paletteLength/2)))
wh <- 13

outFile <- makePath('predictedTF_expressed_greaterthan2','pdf')
pdf(outFile,paper = "a4r")
pheatmap(hmData, 
               border_color = NA, 
               #annotation_row = signifDPFmembership[,'cluster',drop=F],
               cluster_cols = F,  cluster_rows = F,labels_row = normWide_predictedTFs$Name,
               fontsize_row = 6, fontsize_col = 8,
               cellheight = wh,cellwidth = wh, main='predicted TFs - ordered by outdegree',
               color = hmPalette,breaks = myBreaks)


#### Line plots
makeLinePlots(normLong,predictedTFs,title = 'predicted TFs')
dev.off()
########################################################################################
############################################ Check expression of targets

common <- intersect(rownames(normWide), targets)

normWide_predictedTargets <- normWide[,order]
normWide_predictedTargets <- normWide_predictedTargets[names(targets_indegree[names(targets_indegree) %in% common]),] #Arrange by outdegree

normWide_predictedTargets <- normWide_predictedTargets %>% mutate(Gene=rownames(normWide_predictedTargets)) %>% 
  left_join(t2gtair10,by='Gene') %>% select(-11)

hmData <- normWide_predictedTargets[6:10]-normWide_predictedTargets[1:5] 

hmPalette <- colorRampPalette(c('cyan','black','yellow'))(200)
paletteLength <- length(hmPalette)
boundPalette <- max(abs(range(hmData)))
myBreaks <- c(seq(-boundPalette, -0.0001, length.out=ceiling(paletteLength/2)), 
              seq(boundPalette/paletteLength, boundPalette, 
                  length.out=floor(paletteLength/2)))
wh <- 13

outFile <- makePath('target_expression','pdf')
pdf(outFile,paper = "a4r")
pheatmap(hmData, 
         border_color = NA, 
         #annotation_row = signifDPFmembership[,'cluster',drop=F],
         cluster_cols = F,  cluster_rows = F,labels_row = normWide_predictedTargets$Name,
         fontsize_row = 6, fontsize_col = 8,
         cellheight = wh,cellwidth = wh, main='targets - ordered by indegree',
         color = hmPalette,breaks = myBreaks)


#### Line plots
makeLinePlots(normLong,targets,title='targets')
dev.off()

##################################################################

NxTwork <- c(predictedTFs_expessed,targets)
length(NxTwork)

NxTwork <- intersect(rownames(bw2Genome_modelNO3),NxTwork)
NxTwork_DE <- bw2Genome_modelNO3[NxTwork,]
table(NxTwork_DE$adj.P.Val < 0.05)

#table(bw2Genome_modelNO3[rownames(bw2Genome_modelNO3) %in% targets,]$adj.P.Val < 0.05)
#table(bw2Genome_modelNO3[rownames(bw2Genome_modelNO3) %in% predictedTFs,]$adj.P.Val < 0.05)

outFile <- makePath('NxT_Cluster1_FDRvsF','pdf')
pdf(outFile,paper = "a4r")
NxTwork_DE %>%
  ggplot(aes(x=F,y=adj.P.Val)) +
  geom_point(size=0.5) +
  ###
  geom_hline(yintercept = .1,color='skyblue') +
  geom_hline(yintercept = .2,color='orange') +
  geom_hline(yintercept = .05,color='gray') +
  geom_vline(xintercept = min(NxTwork_DE[NxTwork_DE$adj.P.Val < 0.1,'F']),color='skyblue') +
  geom_vline(xintercept = min(NxTwork_DE[NxTwork_DE$adj.P.Val < 0.2,'F']),color='orange') +
  geom_vline(xintercept = min(NxTwork_DE[NxTwork_DE$adj.P.Val < 0.02,'F']),color='gray') +
  theme_classic()
dev.off()

#NxTwork_DE_p2 <- NxTwork_DE[NxTwork_DE$adj.P.Val < 0.2,]
#dim(NxTwork_DE_p2)


#### HM
normWide_nxt <- normWide[rownames(NxTwork_DE),order]

hmData <- normWide_nxt[6:10]-normWide_nxt[1:5] 

hmPalette <- colorRampPalette(c('cyan','black','yellow'))(200)
paletteLength <- length(hmPalette)
boundPalette <- max(abs(range(hmData)))
myBreaks <- c(seq(-boundPalette, -0.0001, length.out=ceiling(paletteLength/2)), 
              seq(boundPalette/paletteLength, boundPalette, 
                  length.out=floor(paletteLength/2)))

hmData <- hmData %>% mutate(Gene=rownames(hmData)) %>% 
  left_join(t2gtair10,by='Gene') %>% select(-6)

rownames(hmData) <- rownames(normWide_nxt)


outFile <- makePath('NXT_Cluster1_patterns','pdf')
pdf(outFile,paper = "a4r")
pheatmap(hmData[,-6], 
       border_color = NA, cutree_rows = 6,
       #annotation_row = signifDPFmembership[,'cluster',drop=F],
       cluster_cols = F,  cluster_rows = T,labels_row = hmData$Name,
       fontsize_row = 6, fontsize_col = 8,
       cellheight = wh,cellwidth = wh, main='predicted TFs - ordered by outdegree',
       color = hmPalette,breaks = myBreaks)


clustH <- hclust(dist(hmData[,-6]))
hMembers <- cutree(clustH,k = 6)

for( k in unique(hMembers)){
  tmp <- names(hMembers[hMembers==k])
  print(makeLinePlots(normLong,tmp,title = k))
}
dev.off()

