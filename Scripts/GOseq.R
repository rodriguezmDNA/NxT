
source("~/Desktop/Nitr_timeSeries/TimeSeries/Scripts/auxFunctions.R")
source("~/Desktop/BioTools/RBioFunctions/metaFunctions_forNetworkAnalysis.R")
source("~/Desktop/BioTools/RBioFunctions/generic_functions.R")

differentGenomicBuilds <- T
## https://pathwaycommons.github.io/guide/primers/statistics/fishers_exact_test/

## Load libraries and Functions
{
  ########
  #biocLite(c("GenomicRanges"))
  #biocLite(c("rtracklayer"))
  library("goseq")
  library(rtracklayer)
  library(GenomicRanges)
  library(Rsamtools)
  library(purrrlyr) #install.packages("purrrlyr")
}

shortname <- "TFmutants_wtNtreat"
DEList <- get(load("../Results/data/normData_contrast_meanWide.RData"))

nrow(DEList)


## Read files
##GO
GFFfile = "../meta/TAIR10_withTransposons.gff"
GOFile = "../meta/ATH_GO_GOSLIM.txt"

source("~/Desktop/BioTools/RBioFunctions/metaFunctions_forNetworkAnalysis.R")

kMembFile <- list.files("../DPF/combat/NitrTimeSeries_k15/",full.names = T,pattern = "Patterns")
fileName <- gsub("../DPF/combat/NitrTimeSeries_|//DPFgenesToPatterns.txt",'',kMembFile)
Clustresults <- read.delim(kMembFile,sep = "\t",as.is = T)


##########################################################################################
##########################################################################################
# ##### Run this chunk if testing genes from a cluster analysis (output from DPF)
{
Clustresults_lst <- lapply(Clustresults,function(x){x[x!='']})
GenesPerCategory <- sapply(Clustresults_lst,length)
}

barplot(GenesPerCategory,ylim = c(0,max(GenesPerCategory)*1.25),las=2)


## All genes
assayed.genes <- rownames(DEList)

###### Prep GO
{ # If reading the 3.2 version
  go <- read.delim("../meta/ATH_GO_GOSLIM.txt",header = F,stringsAsFactors = F)
  
  
  
  colnames(go) <- c(
    "locusName", "TAIRAccession", "AGI", "GOrelationship", "GOterm",
    "GOID", "TAIRKeywordID", "Aspect", "GOslimTerm", "EvidenceCode",
    "EvidenceDesc", "EvidenceWith", "Reference", "Annotator", "Date"
  )
  
  ## -- Select only AGI and GO id
  go.goseq <- go[,c("AGI", "GOID")]
  head(go.goseq)
  
  colnames(go) <- c("ITAG",
                    "GOID")
  
  go.goseq <- go[,c("ITAG", "GOID")]
  head(go.goseq)
}

###### Prep GFF
######### Prepare GFF data
####################################
GFFfile = "../meta/TAIR10_withTransposons.gff"
GFF <- import.gff(GFFfile,version="3",feature.type="gene")
grl <- IRanges::reduce(split(GFF, mcols(GFF)$Name))
reducedGTF <- unlist(grl, use.names=T)
mcols(reducedGTF)$Name <- rep(names(grl), elementNROWS(grl))
reducedGTF
AllLengths <- width(reducedGTF)
names(AllLengths) <- mcols(reducedGTF)$Name
#
head(reducedGTF)
head(AllLengths)


GOList <- list()
maxN <- ncol(GeneList)

for( each in seq(1,maxN)){
  print (each) 
  #tmp <- GeneList[GenePatterns[,each]!="",each]
  
  #assayed.genes <- rownames(tmp)
  ### Since the GeneList matrix is uneven, the function fills in blanks with empty quotes "". 
  de.genes <- GeneList[GeneList[,each]!="",each] 
  
  # goseq needs a binary vector with 0 corresponding to genes not DE/interest and 1 being a DE gene (or gene of interest)
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  table(gene.vector)
  
  ## When using annotations and GO lists from different genome builds, we can remove the info belonging to the build and just leaving the gene stable ID
  
  if (differentGenomicBuilds) names(AllLengths) <- removeDotGene(names(AllLengths))
  
  ## Get the gene lengths only of the genes assayed 
  # (the AllLengths contains information for ALL the genes in the genome but not all of them are *normally* used)
  GeneLengths <- AllLengths[removeDotGene(names(gene.vector))]
  ##
  #cbind(GeneLengths,gene.vector)
  
  all.genes <- names(gene.vector)
  
  ########## GOseq functions
  ##################################################
  ## nullp creates a weight probability of a gene based DE based on it's length
  pwf=nullp(gene.vector, bias.data=GeneLengths)
  head(pwf)
  
  if (differentGenomicBuilds) rownames(pwf) <- removeDotGene( rownames(pwf) )
  
  go.all <- goseq(pwf, "tair10", gene2cat=go.goseq)
  ################################################################################
  
  ### Filter significant categories
  #table(go.all$over_represented_pvalue < 0.05)
  #go.sign <- go.all[go.all$over_represented_pvalue < 0.05,]
  
  # Save to a list
  #View(go.sign)
  GOList[[as.character(colnames(GeneList)[each])]] <- go.all
}


######### Filter by Significant Categories
####################################################################################################
sapply(GOList, dim)
signGOList <- lapply(GOList, function(x){  x[x$over_represented_pvalue < 0.05,] })
sapply(signGOList, nrow)

# ##### Get DE genes on each category
# #########################
# library(tidyverse)
# getGenesPerCategory <- lapply(seq_along(signGOList),function(DP){
#   cat("Genes per category in pattern", DP,"\n")
#   
#   ## Get DE genes
#   de.genes <- GeneList[GeneList[,DP]!="",DP]
#   
#   ### Get categories
#   categories <- signGOList[[DP]][,"category"]
#   
#   ## From the full GO seq list, get the genes that overlap with the DE
#   DEgenesInCategory <- go.goseq %>%
#     filter(.$GOID %in% categories) %>% group_by(GOID) %>% filter(ITAG %in% removeDotGene(de.genes)) %>% nest()
#   ### Convert the table with significant GO categories to a tibble and join the list-columns table
#   test <- as.tibble(signGOList[[DP]])
#   test <- test %>% 
#     left_join(DEgenesInCategory,by=c("category"="GOID")) %>% mutate("Set"=DP) #%>% rename("DEGenesinCategory"="data") 
#   return(test)
# })
# 
# 
# names(getGenesPerCategory) <- names(signGOList)
# names(getGenesPerCategory) 
# 
# ### Make a single table
# #####################
# significantGOtibble <- bind_rows(getGenesPerCategory) %>% 
#   mutate(uniqueID = row_number())
# 
# nrow(significantGOtibble) == sum(sapply(getGenesPerCategory, nrow)) # Should be true
# 
# ## Find the orthologues, not for Ath
# ################################
# # significantGOtibble <- significantGOtibble %>%  # To this table
# #   left_join( # Join the following transformation:
# #     significantGOtibble %>% select(uniqueID,data) %>%
# #       # Group by category (maybe this isn't necessary)
# #       mutate(                  # Create a new column that is
# #         "Sly2AGI_orthologues" = map(data, .f= findOrtho_Sly2AGI)) %>% 
# #       select(uniqueID,3) 
# #     , by = c("uniqueID" = "uniqueID")
# #   ) %>% rename("DE_ITAGsinCats"="data")
# # nrow(significantGOtibble) == sum(sapply(getGenesPerCategory, nrow)) # Should be true
# ################
# 
# ## Fix terms with no name associated
# significantGOtibble <- significantGOtibble %>%  
#   mutate("term"=ifelse( is.na(term),category,term  ) )
# 
# ## Save
# save(significantGOtibble,file = "significantGOtibble.RData")
# ################################
# 
# ## I made pretty graphs
# significantGOtibble %>% filter(!is.na(term)) %>%
#   group_by(ontology,term) %>% summarise("Count"=n()) %>% arrange(ontology,desc(Count)) %>% filter(Count > 1) %>%
#   ggplot(.,aes(x=reorder(term,Count), # Order X axis labels by Count
#                y=Count, #Counts go in the Y axis
#                fill=ontology)) + #Color by ontology
#   geom_bar(stat="identity") + coord_flip() + # Do a bar plot 
#   xlab("GO category") + 
#   facet_grid(ontology~.,scales="free") #Separate by ontology
# 
# significantGOtibble %>% select(Set) %>% group_by(Set) %>% count() %>%
#   ggplot(.,aes(x=as.factor(Set),y=n)) +
#   geom_bar(stat="identity",width = 10/length(unique(significantGOtibble$Set)) ) + # Do a bar plot 
#   xlab("Set") + ylab("# of significant Categories") + theme_light()
# 
# 
# ######### Write genes from GO categories to files
# ####################################################################################################
# ##### Create directories to save outputs
# {
#   resultsPath <- "GenesPerCategory/"
#   itagPath <- paste0(resultsPath,"ITAG_byCategory/")
#   orthoPath <- paste0(resultsPath,"AthOrthologues_byCategory/")
#   dir.create(itagPath,recursive = T) 
#   dir.create(orthoPath,recursive = T)
# }
# 
# unite_(significantGOtibble,"tmpName",c("term","ontology","Set"),sep="_",remove=F) %>%  ## Create new names for the
#   mutate("itagFile"=paste0(itagPath,"/",gsub("[[:space:]]|/|>","-",tmpName),".txt")) %>%
#   by_row(~write.table(.$DE_ITAGsinCats, file = .$itagFile,sep = "\t",row.names = F,quote = F)) %>%
#   mutate("orthoFile"=paste0(orthoPath,"/",gsub("[[:space:]]|/|>","-",tmpName),".txt")) %>%
#   by_row(~write.table(.$Sly2AGI_orthologues, file = .$orthoFile,sep = "\t",row.names = F,quote = F))
# 
# 
# ################################################################################################################################
# ################################################################################################################################
# 
# 
# ########### Subset the p values of overrepresented categories and split by ontology
# ############################################
# {    
#   ## Get only two columns
#   go2 <- lapply(signGOList, function(go){
#     #go <- signGOList[[1]]
#     go[is.na(go$term),"term"] <- go[is.na(go$term),"category"]
#     rownames(go) <- go$term
#     go <- go[,c("ontology","over_represented_pvalue")]
#   })
#   
#   ## Split by category
#   splitted <- lapply(seq_along(go2), function(x){
#     tmpName <- names(go2)[[x]]
#     go3 <- go2[[x]]
#     go3
#     if (nrow(go3) != 0){
#       go3 <- split(go3,go3$ontology)
#       names(go3) <- paste0(tmpName,".",names(go3))
#       return(go3)
#     }
#   })
#   
#   splitted <- unlist(splitted,recursive = F)
#   names(splitted)
#   
#   BioProcess <- splitted[grep("BP",names(splitted))]
#   MolFunction <- splitted[grep("MF",names(splitted))]
#   CellComponent <- splitted[grep("CC",names(splitted))]
#   
#   ##
#   BioProcess <- lapply(seq_along(BioProcess), function(x){
#     DF <- BioProcess[[x]][,2,drop=F]
#     colnames(DF) <- gsub("\\...","",names(BioProcess)[[x]])
#     return(DF)
#   })
#   names(BioProcess) <- lapply(BioProcess,colnames)
#   ##
#   MolFunction <- lapply(seq_along(MolFunction), function(x){
#     DF <- MolFunction[[x]][,2,drop=F]
#     colnames(DF) <- gsub("\\...","",names(MolFunction)[[x]])
#     return(DF)
#   })
#   names(MolFunction) <- lapply(MolFunction,colnames)
#   ##
#   CellComponent <- lapply(seq_along(CellComponent), function(x){
#     DF <- CellComponent[[x]][,2,drop=F]
#     colnames(DF) <- gsub("\\...","",names(CellComponent)[[x]])
#     return(DF)
#   })
#   names(CellComponent) <- lapply(CellComponent,colnames)
#   
#   ##### Convert to tables
#   ###############
#   BioProcessDF <- condenseListTables(BioProcess)
#   MolFunctionDF <- condenseListTables(MolFunction)
#   CellComponentDF <- condenseListTables(CellComponent)
# }
# dir.create("Results",showWarnings = T)
# write.table(BioProcessDF,file = "Results/GOenrichment_BioProcess.tsv",quote = F,row.names = T, sep="\t")
# write.table(MolFunctionDF,file = "Results/GOenrichment_MolFunction.tsv",quote = F,row.names = T, sep="\t")
# write.table(CellComponentDF,file = "Results/GOenrichment_CellComponent.tsv",quote = F,row.names = T, sep="\t")
# 
# 
# 
# 
# 
# 
# # 
# # 
# # 
# # 
# # ######### Prepare GFF data
# # ####################################
# # GFFfile = "meta/ITAG3.2_gene_models.gff"
# # GFF <- import.gff(GFFfile,version="3",feature.type="gene")
# # grl <- IRanges::reduce(split(GFF, mcols(GFF)$Name))
# # reducedGTF <- unlist(grl, use.names=T)
# # mcols(reducedGTF)$Name <- rep(names(grl), elementNROWS(grl))
# # reducedGTF
# # AllLengths <- width(reducedGTF)
# # names(AllLengths) <- mcols(reducedGTF)$Name
# # #
# # head(reducedGTF)
# # head(AllLengths)
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ## Filter logFC of significant genes                             
# # logfcExprs <- lapply(DEList,function(x){x[x[,grep("adj.P.Val",colnames(x))] < 0.05,grep("logFC",colnames(x)),drop=F]})
# # deGenesList <- lapply(logfcExprs, rownames)
# # 
# # 
# # ### GO terms of DE genes in data
# # ## Get GO terms enriched within each group
# # background <- Reduce(union,lapply(DEList,rownames))
# # par(mfrow=c(3,3))
# # GOenrich_DE <- findGOTerms(deGenesList,background,GFFfile,GOFile)
# # par(mfrow=c(1,1))
# # lapply(GOenrich_DE,dim)
# # ## --
# # 
# # # Filter by pValue and only Biological Process
# # sapply(GOenrich_DE, nrow)
# # signGOList <- lapply(GOenrich_DE, function(x){  x[x$over_represented_pvalue < 0.05 & x$ontology == "BP",] })
# # ## Use only overrepresentated pValue
# # for (each in names(signGOList)){
# #   rownames(signGOList[[each]]) <- signGOList[[each]]$term
# #   signGOList[[each]] <- signGOList[[each]][,"over_represented_pvalue",drop=F]
# #   colnames(signGOList[[each]]) <- each
# # }
# # sapply(signGOList, nrow)
# # 
# # ## Heatmap of GO categories per cluster
# # library(RColorBrewer)
# # library(gplots)
# # ## --
# # titulo <- "GOEnrichment_from_DE_genes"
# # colors <- colorRampPalette(c("skyblue","steelblue2","steelblue4"))
# # path <- paste("",titulo,"_",shortname,".pdf",sep="")
# # pdf(path,paper = "a4")
# # #
# # hmData <- as.matrix(condenseListTables(signGOList))
# # head(hmData)
# # hmData[hmData==0] <- NaN
# # head(hmData)
# # hmData <- -(log10(hmData))
# # #######
# # heatmap.2(hmData,col=colors(120),
# #           keysize = 1.5,
# #           symkey = F,
# #           #key.par=list(mar=c(3.5,0,3,0)),
# #           na.color = "white",
# #           margins = c(2,6),
# #           
# #           # lmat -- added 2 lattice sections (5 and 6) for padding
# #           lmat=rbind(c(3, 8, 2), 
# #                      c(7, 4, 9),
# #                      c(6, 5, 1)), 
# #           #lhei=c(1.5, 0,5), 
# #           #lwid=c(2.5, 4, 1),
# #           #lmat = lmat, 
# #           lhei=c(0.15,0.22,0.9),
# #           lwid=c(0.8,0.4,0.3),
# #           ##
# #           scale="none",
# #           density.info = "none", 
# #           key.xlab = "-log10(pVal)", 
# #           trace = "none",
# #           dendrogram = "none",
# #           cexRow = 0.064,
# #           cexCol = 0.15,
# #           Rowv = F,Colv = F,
# #           main="GO Enrichment",cex.main=0.5)
# # dev.off()
# # 
# # 
# # titulo <- paste("",titulo,"_",shortname,".csv",sep="")
# # write.csv(titulo,x = hmData)