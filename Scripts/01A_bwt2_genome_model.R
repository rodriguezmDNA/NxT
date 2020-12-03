#BiocManager::install("tximport")
#BiocManager::install("Glimma")
#BiocManager::install("sva")
library(tidyverse)
library(edgeR)
library(Glimma)
library(sva)
library(reshape)
#http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#downstream_dge_in_bioconductor

dir.create("../Results/bwt2_genome",recursive = T,showWarnings =F)

t2gtair10 <- read.table("~/Desktop/Spring2020/scRNASeq/t2gtair10.txt",
                        header = F,col.names = c("","Gene","Name"),
                        sep = "\t",quote = "",row.names = 1,stringsAsFactors = F)
t2gtair10$Name <- ifelse(t2gtair10$Name=="",t2gtair10$Gene,t2gtair10$Name)
head(t2gtair10)

t2gtair10 <- t2gtair10[!duplicated(t2gtair10$Gene),]

### Read experiment metadata
meta <- read.delim("../meta/TimeSeriesExperiment_fixStrips4and5.txt",
                   header = T,sep = "\t",stringsAsFactors = F)
meta$Strip <- gsub(" ","",meta$Strip)
meta$Time <- paste0("t",meta$Time)
meta$sample <- paste0(meta$Time,"_",meta$NO3,"_",meta$Strip,"_set",meta$Set)
meta$group <- paste0(meta$Time,"_",meta$NO3)

### Read counts

### Prep edgeR object
counts <- read.delim("../RawCounts/RawCounts_bwt2_genome.csv",sep = ",",
                     row.names = 1,header = T,stringsAsFactors = F)
colnames(counts) <- gsub("\\.","-",colnames(counts))
head(counts)

#### Filter
## This
removeMeta <- meta$sample[meta$Time == "t0min"]
removeCounts <- meta$File[meta$Time == "t0min"]
meta <- meta[!meta$sample%in%removeMeta,,drop=F]
counts <- counts[,!colnames(counts)%in%removeCounts,drop=F]

meta$Time <- factor(meta$Time)
meta$Time <- relevel(meta$Time,ref = "t15min")
## Or
#meta$Time <- factor(meta$Time)
#meta$Time <- relevel(meta$Time,ref = "t0min")
meta$NO3 <- factor(meta$NO3)
########

nrow(meta) == ncol(counts) & all(meta$File %in% colnames(counts))

##

design <- model.matrix(~Time*NO3, data =meta)#+Lane
colnames(design) <- gsub("meta|Time|NO3|Strip|:|-|/|group","",colnames(design))
head(design)

dge <- DGEList(counts=counts,samples = meta,remove.zeros = T)
# Filter: Genes with total counts more than 

sampleMin <- min(table(meta$group))
minCPM <- 1
isexpr <- rowSums(cpm(dge) > minCPM) >= sampleMin #Make sure to use the minimum 
dge <- dge[isexpr,,keep.lib.size = FALSE]
dim(dge)

dge <- calcNormFactors(dge)
v <- voomWithQualityWeights(dge, design=design, normalize.method = "quantile", plot=TRUE)
dim(v)


boxplot(v$E, range=0,
        ylab="log2[counts]", xlab="sample", main="voom normalized counts",
        cex.axis=0.5,las=2)


glMDSPlot(v, labels=meta$sample, groups=meta,
          folder="../glimma/bwt2_genome_sanst0_NxT",
          launch=T)

modcombat = model.matrix(~Time*NO3, data =meta)#+Lane
batch = meta$Strip
combat_edata <- ComBat(dat=v$E, batch=batch,mod = modcombat)

glMDSPlot(combat_edata, labels=meta$sample, groups=meta,
          folder="../glimma/bwt2_genome_combat_sanst0_NxT",
          launch=T)

fit <- lmFit(combat_edata, design)

fit <- eBayes(fit)
cbind(colnames(fit$coefficients))
DEG <- topTable(fit,coef=7:10,number=Inf,sort.by="none")
sum(DEG$adj.P.Val<0.05)

save(DEG,file="../Results/bwt2_genome/DEanalysis_anova.RData")
signifDEG <- DEG[DEG$adj.P.Val<0.05,]


voomData <- v$E
combatData <- combat_edata
colnames(voomData) <- meta$sample
colnames(combatData) <- meta$sample

save(voomData,file="../Results/bwt2_genome/voomData.RData")
save(combatData,file="../Results/bwt2_genome/combatData.RData")
# colnames(voomData) <- paste0("v_",colnames(voomData))
# colnames(combatData) <- paste0("c_",colnames(combatData))
# 
# pdf("../06_Results/bwt2_genome/Corr_Voom-Combat.pdf",width = 20,height = 15)
# pheatmap::pheatmap(cor(cbind(voomData,combatData)),border_color = NA)
# dev.off()

source("auxFunctions.R")

mean_voomData_long <- MeanDataLong(voomData)
mean_combatData_long <- MeanDataLong(combatData)
save(mean_voomData_long,file="../Results/bwt2_genome/mean_voomData_long.RData")
save(mean_combatData_long,file="../Results/bwt2_genome/mean_combatData_long.RData")


DEG_NO3 <- topTable(fit,coef=6:10,number=Inf,sort.by="none")
sum(DEG_NO3$adj.P.Val<0.05)
save(DEG_NO3,file="../Results/bwt2_genome/DEanalysis_anova_NO3.RData")


#### Making sense of the results
####################################################################
#### Making sense of the results
# pltMe = "AT5G67300"
# print(DEG[pltMe,c(1,4,5,3,2,6,7,8,9)])
# 
# wide_combatExpression <- mean_combatData_long %>% filter(Gene %in% pltMe) %>%
#   select(Gene,Group,meanExpr) %>% spread(key=Group,value=meanExpr)
# 
# wide_combatExpression <- data.frame(wide_combatExpression[,-1],row.names = wide_combatExpression$Gene)
# 
# wide_combatExpression <- wide_combatExpression[,c(3,7,9,5,1,
#                                                   4,8,10,6,2)]
#                                                   
# ctrst <- tibble("gen"=pltMe,"logFC"=c( ### Over time vs reference at 15min:
#   (wide_combatExpression[,6]-wide_combatExpression[,1]),
#   (wide_combatExpression[,7]-wide_combatExpression[,6])-(wide_combatExpression[,2]-wide_combatExpression[,1]),
#   (wide_combatExpression[,8]-wide_combatExpression[,6])-(wide_combatExpression[,3]-wide_combatExpression[,1]),
#   (wide_combatExpression[,9]-wide_combatExpression[,6])-(wide_combatExpression[,4]-wide_combatExpression[,1]),
#   (wide_combatExpression[,10]-wide_combatExpression[,6])-(wide_combatExpression[,5]-wide_combatExpression[,1])),
#   "time"=factor(c("15","45","90","180","1200"),levels = c("15","45","90","180","1200")))
# print(ctrst)
# ctrst %>%
#   ggplot(aes(x=time,y=logFC)) +
#   geom_point() + ggtitle(pltMe) + geom_line(aes(group=NA)) +
#   geom_text(aes(x=time,y=max(logFC)*1.3,label=signif(logFC,3))) +
#   theme_bw()
# 
# mean_combatData_long %>% filter(Gene %in% pltMe) %>%
#   mutate("t"=factor(gsub("t|min","",Time),levels=c("0","15","45","90","180","1200"))) %>%
#   ggplot(aes(x=t,y=meanExpr)) +
#   ggtitle(t2gtair10[t2gtair10$Gene==pltMe,"Name"]) +
#   geom_line(aes(color=NO3,group=NO3)) +
#   geom_point(aes(color=NO3)) +
#   facet_grid(Name~.,scales = "free_y") +
#  theme_minimal()
###################################################################