  setwd("~/Rwork/1803-SCourbier/GoSeq/")

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

source("~/Rwork/JoelsRcode/GOSeq_jrmScript/omicFunctions/metaFunctions_forNetworkAnalysis.R")
source("~/Rwork/JoelsRcode/GOSeq_jrmScript/omicFunctions/generic_functions.R")
source("~/Rwork/JoelsRcode/GOSeq_jrmScript/omicFunctions/orthologueFindingFunctions.R")
source("~/Rwork/JoelsRcode/GOSeq_jrmScript/omicFunctions/GOseq_functions.R")

differentGenomicBuilds <- T
######


## load DE genes
load ("DEList_2018.05.30.Rdata")
## Load genes to test for enrichment
###  Either from DPF
#setwd("C:/Users/courb001/surfdrive/surfdrive/Utrecht/GOSeq_Joel2/") 
  ##########################################################################################
  ##########################################################################################
  # ##### Run this chunk if testing genes from a cluster analysis (output from DPF)
  # {  
  # GeneList <- read.table("DPFgenesToPatterns.txt",sep="\t",header = T,as.is = T)
  # GenesPerCategory <- sapply(apply(GeneList,2,function(x){table(x!="")}),"[",1)
  # names(GenesPerCategory) <- gsub(".FALSE","",names(GenesPerCategory))
  # }
  # #########################

  
  ## Or from a list of differential expression tables (output from limma/voom)
  {
    GeneList <- lapply(DEList,function(x){
    ## Get genes that are significant
    significant <- x[x[,grep("adj.P",colnames(x))] <= 0.05,]
    return ( rownames(significant))
  })
    GenesPerCategory <- sapply(GeneList, length)
    GeneList <- condenseGeneList_toMatrix(GeneList)
  }
  #############################################
  
  
barplot(GenesPerCategory,ylim = c(0,max(GenesPerCategory)*1.25),las=2)
  
  ## All genes
  assayed.genes <- Reduce(union,lapply(DEList,rownames))

  ######### Prepare GFF data
  ####################################
  GFFfile = "meta/ITAG3.2_gene_models.gff"
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

#########


######### Prepare list of GO terms
####################################

{ # If reading the 3.2 version
  go <- read.delim("meta/ITAG3.2_protein_go_Split_jrm.tsv",header = F,stringsAsFactors = F)
  colnames(go) <- c("ITAG",
                    "GOID")
  
  go.goseq <- go[,c("ITAG", "GOID")]
  head(go.goseq)
}
#########

##################
## The function readGO takes a GO dataframe of two columns, where each row is a single gene (first column)
# second column containes the categories associated with the gene. 
# readGO reshapes the table to have one category per row, which is then passed to GOseq
##################
# {
# # Read the plant GOslim annotation
# #######
#   GOitag <- "meta/0728fas.blast.map.annot.interpro.annex.plant_slim.GOstat.txt"
#   go.goseq <- readGO(GOitag);head(go.goseq)
#   go.goseq
#   nrow(go.goseq)
# }
# ## Or the merge
# {
#   GOitag <- "meta/0728fas.blast.map.annot.interpro.annex.GOstatMerge.txt"
#   go.goseq <- readGO(GOitag);head(go.goseq)
#   nrow(go.goseq)
# }
# ####

######### GO terms for each cluster
####################################################################################################

### This script assumes you've created a matrix of genes of interest (GeneList) using the condenseGeneList_toMatrix() function.
# This function essentially creates a table with each column contains genes that were significant or of interest for a particular experiment (ie, a contrast or the result of clustering)
# The function takes a list of genes and creates an uneven (not symmetrical) matrix. 

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
  if (differentGenomicBuilds) go.goseq$ITAG <- removeDotGene ( go.goseq$ITAG)  
  
  
  go.all <- goseq(pwf, "ITAG3.1", gene2cat=go.goseq)
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

##### Get DE genes on each category
#########################
library(tidyverse)
getGenesPerCategory <- lapply(seq_along(signGOList),function(DP){
  cat("Genes per category in pattern", DP,"\n")
  
  ## Get DE genes
  de.genes <- GeneList[GeneList[,DP]!="",DP]
  
  ### Get categories
  categories <- signGOList[[DP]][,"category"]
  
  ## From the full GO seq list, get the genes that overlap with the DE
  DEgenesInCategory <- go.goseq %>%
    filter(.$GOID %in% categories) %>% group_by(GOID) %>% filter(ITAG %in% removeDotGene(de.genes)) %>% nest()
  ### Convert the table with significant GO categories to a tibble and join the list-columns table
  test <- as.tibble(signGOList[[DP]])
  test <- test %>% 
    left_join(DEgenesInCategory,by=c("category"="GOID")) %>% mutate("Set"=DP) #%>% rename("DEGenesinCategory"="data") 
  return(test)
})


names(getGenesPerCategory) <- names(signGOList)
names(getGenesPerCategory) 

### Make a single table
#####################
significantGOtibble <- bind_rows(getGenesPerCategory) %>% 
  mutate(uniqueID = row_number())

nrow(significantGOtibble) == sum(sapply(getGenesPerCategory, nrow)) # Should be true

## Find the orthologues
################################
significantGOtibble <- significantGOtibble %>%  # To this table
  left_join( # Join the following transformation:
    significantGOtibble %>% select(uniqueID,data) %>%
      # Group by category (maybe this isn't necessary)
      mutate(                  # Create a new column that is
        "Sly2AGI_orthologues" = map(data, .f= findOrtho_Sly2AGI)) %>% 
      select(uniqueID,3) 
    , by = c("uniqueID" = "uniqueID")
  ) %>% rename("DE_ITAGsinCats"="data")
nrow(significantGOtibble) == sum(sapply(getGenesPerCategory, nrow)) # Should be true
################

## Fix terms with no name associated
significantGOtibble <- significantGOtibble %>%  
  mutate("term"=ifelse( is.na(term),category,term  ) )

## Save
save(significantGOtibble,file = "significantGOtibble.RData")
################################

## I made pretty graphs
significantGOtibble %>% filter(!is.na(term)) %>%
  group_by(ontology,term) %>% summarise("Count"=n()) %>% arrange(ontology,desc(Count)) %>% filter(Count > 1) %>%
  ggplot(.,aes(x=reorder(term,Count), # Order X axis labels by Count
               y=Count, #Counts go in the Y axis
               fill=ontology)) + #Color by ontology
  geom_bar(stat="identity") + coord_flip() + # Do a bar plot 
  xlab("GO category") + 
  facet_grid(ontology~.,scales="free") #Separate by ontology

significantGOtibble %>% select(Set) %>% group_by(Set) %>% count() %>%
  ggplot(.,aes(x=as.factor(Set),y=n)) +
  geom_bar(stat="identity",width = 10/length(unique(significantGOtibble$Set)) ) + # Do a bar plot 
  xlab("Set") + ylab("# of significant Categories") + theme_light()


######### Write genes from GO categories to files
####################################################################################################
##### Create directories to save outputs
{
  resultsPath <- "GenesPerCategory/"
  itagPath <- paste0(resultsPath,"ITAG_byCategory/")
  orthoPath <- paste0(resultsPath,"AthOrthologues_byCategory/")
  dir.create(itagPath,recursive = T) 
  dir.create(orthoPath,recursive = T)
}

unite_(significantGOtibble,"tmpName",c("term","ontology","Set"),sep="_",remove=F) %>%  ## Create new names for the
  mutate("itagFile"=paste0(itagPath,"/",gsub("[[:space:]]|/|>","-",tmpName),".txt")) %>%
  by_row(~write.table(.$DE_ITAGsinCats, file = .$itagFile,sep = "\t",row.names = F,quote = F)) %>%
  mutate("orthoFile"=paste0(orthoPath,"/",gsub("[[:space:]]|/|>","-",tmpName),".txt")) %>%
  by_row(~write.table(.$Sly2AGI_orthologues, file = .$orthoFile,sep = "\t",row.names = F,quote = F))


################################################################################################################################
################################################################################################################################


########### Subset the p values of overrepresented categories and split by ontology
############################################
{    
  ## Get only two columns
  go2 <- lapply(signGOList, function(go){
    #go <- signGOList[[1]]
    go[is.na(go$term),"term"] <- go[is.na(go$term),"category"]
    rownames(go) <- go$term
    go <- go[,c("ontology","over_represented_pvalue")]
  })
  
  ## Split by category
  splitted <- lapply(seq_along(go2), function(x){
    tmpName <- names(go2)[[x]]
    go3 <- go2[[x]]
    go3
    if (nrow(go3) != 0){
      go3 <- split(go3,go3$ontology)
      names(go3) <- paste0(tmpName,".",names(go3))
      return(go3)
    }
  })
  
  splitted <- unlist(splitted,recursive = F)
  names(splitted)
  
  BioProcess <- splitted[grep("BP",names(splitted))]
  MolFunction <- splitted[grep("MF",names(splitted))]
  CellComponent <- splitted[grep("CC",names(splitted))]
  
  ##
  BioProcess <- lapply(seq_along(BioProcess), function(x){
    DF <- BioProcess[[x]][,2,drop=F]
    colnames(DF) <- gsub("\\...","",names(BioProcess)[[x]])
    return(DF)
  })
  names(BioProcess) <- lapply(BioProcess,colnames)
  ##
  MolFunction <- lapply(seq_along(MolFunction), function(x){
    DF <- MolFunction[[x]][,2,drop=F]
    colnames(DF) <- gsub("\\...","",names(MolFunction)[[x]])
    return(DF)
  })
  names(MolFunction) <- lapply(MolFunction,colnames)
  ##
  CellComponent <- lapply(seq_along(CellComponent), function(x){
    DF <- CellComponent[[x]][,2,drop=F]
    colnames(DF) <- gsub("\\...","",names(CellComponent)[[x]])
    return(DF)
  })
  names(CellComponent) <- lapply(CellComponent,colnames)
  
  ##### Convert to tables
  ###############
  BioProcessDF <- condenseListTables(BioProcess)
  MolFunctionDF <- condenseListTables(MolFunction)
  CellComponentDF <- condenseListTables(CellComponent)
}
dir.create("Results",showWarnings = T)
write.table(BioProcessDF,file = "Results/GOenrichment_BioProcess.tsv",quote = F,row.names = T, sep="\t")
write.table(MolFunctionDF,file = "Results/GOenrichment_MolFunction.tsv",quote = F,row.names = T, sep="\t")
write.table(CellComponentDF,file = "Results/GOenrichment_CellComponent.tsv",quote = F,row.names = T, sep="\t")

# BioProcessDF[grep("protein",rownames(BioProcessDF)),,drop=F]
# BioProcessDF[grep("trans",rownames(BioProcessDF)),,drop=F]
# BioProcessDF[grep("translation",rownames(BioProcessDF)),,drop=F]
# BioProcessDF[grep("ribo",rownames(BioProcessDF)),,drop=F]
# 
# BioProcessDF[grep("translation|ribo",rownames(BioProcessDF)),,drop=F]
# 
# BioProcessDF[grep("photosyn",rownames(BioProcessDF)),,drop=F]
# BioProcessDF[grep("nitr",rownames(BioProcessDF)),,drop=F]


test <- -log10(BioProcessDF)+100
test[is.na(test)] <- 0

range(test)
names(test) <- names(DEList)
#test <- test[grep("protein|photo|transl|nitr|ribo",rownames(test)),]


#color <- colorRampPalette(c("white","skyblue","purple")) #No
#color <- colorRampPalette(c("white",RColorBrewer::brewer.pal(2,"YlGn"))) #No
#color <- colorRampPalette(c("white",RColorBrewer::brewer.pal(2,"YlGnBu"))) #No
color <- colorRampPalette(c("white",RColorBrewer::brewer.pal(9,"Blues"))) #Yes
#color <- colorRampPalette(c("white",RColorBrewer::brewer.pal(9,"Oranges"))) #Yes
#color <- colorRampPalette(c("white",RColorBrewer::brewer.pal(9,"PuBuGn"))) #Yes
#color <- colorRampPalette(c("white",cm.colors(20)))
#dev.off()
pdf("GOSeq_filter.pdf",width = 8,height = 10)
gplots::heatmap.2(as.matrix(test), dendrogram = "both", margins = c(12,26), 
                  trace = "none", symm = F, keysize = 0.8, key.title = "-log10(pVal)+100",
                  col = color(100))
dev.off()

