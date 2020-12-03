order<- c(
  #"t0min_1mM",
  "t15min_10mM","t45min_10mM","t90min_10mM","t180min_10mM","t1200min_10mM",
  "t15min_1mM","t45min_1mM","t90min_1mM","t180min_1mM","t1200min_1mM")

############################################################
t2gtair10 <- read.table("~/Desktop/Spring2020/scRNASeq/t2gtair10.txt",
                        header = F,col.names = c("","Gene","Name"),
                        sep = "\t",quote = "",row.names = 1,stringsAsFactors = F)
t2gtair10$Name <- ifelse(t2gtair10$Name=="",t2gtair10$Gene,t2gtair10$Name)
head(t2gtair10)

t2gtair10 <- t2gtair10[!duplicated(t2gtair10$Gene),]
############################################################

clustAnalysis <- function(){
 ############### Analysis
  GenesXClusterDF <- data.frame('cluster'=k.means.fit$cluster,'Gene'=names(k.means.fit$cluster),stringsAsFactors = F)
  exprData <- normLong %>% left_join(GenesXClusterDF,by = "Gene") %>% filter(!is.na(cluster))
  
  exprDataXDPF <- exprData %>% group_by(Group,cluster,Time,NO3) %>% 
    summarise("meanXcluster"=mean(meanExpr)) %>% ungroup() 
  
  exprDataXDPF <- exprDataXDPF %>% 
    mutate(cluster=factor(exprDataXDPF$cluster,
                          levels = gtools::mixedsort(unique(exprDataXDPF$cluster))))
  
  GenesXCluster <- table(k.means.fit$cluster)
  
  GenesXCluster <- as.data.frame(GenesXCluster)
  colnames(GenesXCluster) <- c("Pattern",'n')
  
  #### Genes per cluster
  GxP <- GenesXCluster %>% 
    ggplot(aes(x=factor(gtools::mixedorder(Pattern,decreasing = T)),y=n)) +
    geom_bar(stat="identity",width = 0.75,fill="skyblue") +
    geom_text(aes(x=factor(gtools::mixedorder(Pattern,decreasing = T)),y=median(n),label=n),size=3,color="black") +
    xlab("Cluster") +
    ylab("Genes assigned") +
    coord_flip() +
    theme_light()
  print(GxP)
  
  ########################################################################
  #### Explore avg expression by cluster
  linPlot <- exprDataXDPF %>% 
    #mutate(Time=as.numeric(gsub("t|min","",exprDataXDPF$Time))) %>%
    mutate(Time=factor(gsub("t|min","",Time),levels = gtools::mixedsort( unique(gsub("t|min","",Time))))) %>%  
    ggplot(aes(x=Time,y=meanXcluster)) +
    geom_point(aes(color=NO3,group=NO3,shape=NO3) ) +
    geom_line(aes(color=NO3,group=NO3) ) +
    theme(aspect.ratio=16/9) +
    facet_wrap(as.factor(cluster)~.,scales = "free_x",nrow = 1) +
    theme_classic()
  print(linPlot)
  
  #### Heatmap
  meanXcluster <- exprDataXDPF %>% mutate(Time=as.numeric(gsub("t|min","",exprDataXDPF$Time)))  %>% 
    pivot_wider(names_from = Group,values_from = meanXcluster,id_cols = cluster) %>% data.frame(.,row.names = .$cluster) %>% .[,-1]
  
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
  hm <- pheatmap(hmData, 
           border_color = NA, 
           #annotation_row = signifDPFmembership[,'cluster',drop=F],
           cluster_cols = F,  cluster_rows = T,
           fontsize_row = 6, fontsize_col = 8,
           cellheight = wh,cellwidth = wh,
           color = hmPalette,breaks = myBreaks)
  dev.off()
  
  Clustresults_lst <- lapply(split(GenesXClusterDF,GenesXClusterDF$cluster),function(x){x[,2]})
  
  ######## 
  outList <- list()
  outList[['clustsize']] = GxP
  outList[['lineplot']] = linPlot
  outList[['hm']] = hm
  outList[['genemembership']]=Clustresults_lst
  return(outList)
  
}


#### FET x Cluster
####################################################
FET <- function(kMembers){
  out <- lapply(kMembers,function(x){
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
  return(out)
}


####

analyzeFET <- function(fishRes){
  
  contingencyMatrices <- lapply(fishRes,"[[",1)
  contMatMembers <- lapply(fishRes,"[[",3)
  
  # Correct p value
  ####################
  FETresults <- do.call("rbind",lapply(fishRes,function(x){
    data.frame("OR"=x[["FET"]]$estimate,
               "pVal"=x[["FET"]]$p.value)}))
  FETresults$FDR <- p.adjust(FETresults$pVal,method = "fdr")
  
  ####################
  stars <- ifelse(FETresults$FDR < 0.01,"*","")
  hmData <- FETresults[,1,drop=F]
  wh <- 14
  hmPalette <- colorRampPalette(brewer.pal(n=7,"Blues"))(20)
  hm <- pheatmap(hmData,border_color = NA, display_numbers = cbind(stars), number_color = "white",
                 cluster_cols = F,  cluster_rows = F,
                 #fontsize_row = 6, fontsize_col = 8,
                 cellheight = wh,cellwidth = wh, main='enrichment of clusters in the network',
                 color = hmPalette)
  dev.off()
  out <- list()
  out[['hm']] <- hm
  out[['FETresults']] <- FETresults
  out[['contingencyMatrices']] <- contingencyMatrices
  out[['contMatMembers']] <- contMatMembers
  ####################
  return(out)
}


hmByCluster <- function(genemembership){
  signifDPFmembership <- do.call('rbind',lapply(names(genemembership),function(x){
    data.frame(cbind('cluster'=x,'gene'=genemembership[[x]]),
               row.names = genemembership[[x]],stringsAsFactors = F) }))
  
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
  
  wh <- 8
  pheatmap(hmData_FC, labels_row = geneNames$Name,
           border_color = NA, gaps_row = cumsum(clustMemb),
           annotation_row = signifDPFmembership[,'cluster',drop=F],
           cluster_cols = F,  cluster_rows = F,
           fontsize_row = 3, fontsize_col = 8,
           #cellheight = wh,
           cellwidth = wh,
           color = hmPalette,breaks = myBreaks)
}