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

design <- model.matrix(~Time+NO3, data =meta)#+Lane
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
          folder="../glimma/bwt2_genome_sanst0_NplusT",
          launch=T)

modcombat = model.matrix(~Time+NO3, data =meta)#+Lane
batch = meta$Strip
combat_edata <- ComBat(dat=v$E, batch=batch,mod = modcombat)

glMDSPlot(combat_edata, labels=meta$sample, groups=meta,
          folder="../glimma/bwt2_genome_combat_sanst0_NplusT",
          launch=T)

fit <- lmFit(combat_edata, design)

fit <- eBayes(fit)
cbind(colnames(fit$coefficients))
DEG <- topTable(fit,coef=6,number=Inf,sort.by="none") #N only
sum(DEG$adj.P.Val<0.05)

save(DEG,file="../Results/bwt2_genome/onlyNO3model.RData")
signifDEG <- DEG[DEG$adj.P.Val<0.05,]
