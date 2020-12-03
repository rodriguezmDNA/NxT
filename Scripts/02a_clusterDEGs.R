library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

source("~/Desktop/Nitr_timeSeries/TimeSeries/Scripts/auxFunctions.R")
source("~/Desktop/BioTools/RBioFunctions/metaFunctions_forNetworkAnalysis.R")
source("helper_exploreDEGs.R")

####################
set.seed(42)
dir.create("../Results/02_DEG_Clustering")
####################################################
### Enrichment of Network
YNM <- read.table("../meta/YNM_S03.txt",stringsAsFactors = F,header = T,sep="\t")
NitroGenes <- unique(c(YNM$TF_AGI,YNM$Promoter_AGI))
length(NitroGenes)

### Expression data
normLong <- get(load("../Results/data/normData_contrast_meanLong.RData"))
normWide <- get(load("../Results/data/normData_contrast_meanWide.RData"))

################
bw2Genome_modelNO3 <- get(load("../Results/bwt2_genome/DEanalysis_anova_NO3.RData"))
signifNO3 <- bw2Genome_modelNO3[bw2Genome_modelNO3$adj.P.Val < 0.05,]
Universe <- rownames(bw2Genome_modelNO3)

################
clusDat <- normWide[rownames(signifNO3),]
dim(clusDat)

wssplot <- function(data, nc=15, seed=42){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

## Use maximum number of clusters possible
pdf("../Results/02_DEG_Clustering/elbow.pdf",paper = 'a4r')
wssplot(clusDat, nc=10) 
dev.off()


library(cluster)
kChoice <- 4
k.means.fit <- kmeans(clusDat, kChoice,iter.max = 1000) ##Maximum is nrow-1
table(k.means.fit$cluster)
### Analysis of cluster membership
k <- clustAnalysis()
print(k$clustsize)



pdf("../Results/02_DEG_Clustering/k4_clusters_StdDev_PerNxT.pdf",paper = "a4r")
#### Coeff of variation
exprData %>% group_by(Group,cluster,Time,NO3) %>% 
  summarise("sd"=sd(meanExpr)) %>% ungroup()  %>%
  ggplot(aes(x=factor(cluster),y=sd)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(aes(color=Group),width = 0.15) +
  facet_wrap(as.factor(cluster)~.,scales = "free_x",nrow = 1) +
  theme_classic()
dev.off()

pdf("../Results/02_DEG_Clustering/k4_clusters_avg.pdf",paper = "a4r")
genemembership <- k$genemembership
print(k$hm)
print(k$lineplot)
print(k$clustsize)
dev.off()

pdf("../Results/02_DEG_Clustering/k4_clusters_heatmap.pdf",width = 20,height = 10)
hmByCluster(genemembership)
dev.off()


pdf("../Results/02_DEG_Clustering/avg_expression_Cluster_k4.pdf",paper = "a4r")
exprDataXDPF %>% #filter(cluster == 1) %>%
  #mutate(Time=as.numeric(gsub("t|min","",exprDataXDPF$Time))) %>%
  mutate(Time=factor(gsub("t|min","",Time),levels = gtools::mixedsort( unique(gsub("t|min","",Time))))) %>%  
  ggplot(aes(x=Time,y=meanXcluster)) +
  geom_point(aes(color=NO3,group=NO3,shape=NO3) ) +
  geom_line(aes(color=NO3,group=NO3) ) +
  theme(aspect.ratio=6/24) +
  facet_wrap(.~as.factor(cluster),scales = "free_y",nrow = 4) +
  theme_classic()
dev.off()


#### Fisher exact test
fishRes <- FET(k$genemembership)
outFET <- analyzeFET(fishRes)
#outFET$hm
outFET$contingencyMatrices

pdf("../Results/02_DEG_Clustering/k4_fisherET.pdf",paper = 'a4r')
ClusterXNetwork <- lapply(outFET$contMatMembers,function(x){tmp <- x[,1]; tmp[tmp!='']})
  print(
  tibble('cluster'=names(ClusterXNetwork),'networkGenesInCluster'=sapply(ClusterXNetwork, length)) %>%
    ggplot(aes(x=cluster,y=networkGenesInCluster)) +
    geom_bar(stat="identity",width = 0.75,fill="skyblue") +
    geom_text(aes(x=cluster,y=median(networkGenesInCluster),label=networkGenesInCluster),size=3,color="black") +
    xlab("Cluster") +
    ylab("YNM genes in cluster") +
    coord_flip() +
    theme_light()
  )
### Heatmap of DE genes overlapping with YNM genes
NitroGeneMembers <- lapply(k$genemembership,function(x){ x[x%in% NitroGenes] })
NitroGeneMembers <- NitroGeneMembers[sapply(NitroGeneMembers,length) > 0]

plot.new()
outFET$hm
hmByCluster(NitroGeneMembers)
dev.off()

########### Further explore genes in cluster

genemembership_matrix <- condenseGeneList_toMatrix(genemembership)

pathOut <- paste0("../Results/02b_DEgenes_ClusterAnalysis/","k_",kChoice)
dir.create(pathOut,recursive = T,showWarnings = F)

lapply(names(genemembership),function(x){
  outFile <- paste(pathOut,'/cluster_',x,'.txt',sep='')
  outThis <- data.frame(genemembership[[x]])
  write.table(outThis,file = outFile,quote = F,sep = "\t",row.names = F,col.names = F)
  return()
})



####  Code to call python and extract promoter sequence for genes in each cluster
# for each in `ls ../k_4/*.txt`;
# do
# python filter_fastx.py -l $each -f ../Arabidopsis_thaliana.TAIR10.47.500bp_Upstream.fa
# done

####  Code to perform motif analysis on each set of promoters
####################################
# 
# for each in `ls GeneLists/*`;
# do
# fName=$(basename $each)
# fName=${fName%.fa}
# outMEME="MEMEoutput/"$fName
# outAME="AMEoutput/"$fName
# echo $each $fName $outName
# echo -e "Processing" $fName "to" $outName;
# #meme $each -dna -oc . -nostatus -mod anr -nmotifs 8 -minw 6 -maxw 50 -objfun classic -markov_order 0 -revcomp -o $outMEME;
# ame --verbose 1 --oc . --scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 --control --shuffle-- --kmer 2 --o $outAME $each ArabidopsisDAPv1.meme
# done

################