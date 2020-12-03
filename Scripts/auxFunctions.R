
require(tidyverse)
MeanDataLong <- function(expressionData){

  expressionData_long <- data.frame(expressionData) %>% mutate("Gene"=gsub("\\..$*","",rownames(expressionData))) %>% 
    pivot_longer(cols = colnames(expressionData))  %>% mutate("Group"=gsub("_Strip.*$","",name))
  
  expressionData_long <- expressionData_long %>% group_by(Gene,Group) %>%
    summarise("meanExpr"=mean(value)) %>% ungroup() %>%
    mutate("Time"=gsub("_1.*mM.*$","",Group),
           "NO3"=gsub("t.*min_|_Strip.*$","",Group)) %>% left_join(t2gtair10,by = "Gene")
  return(expressionData_long)
}


order<- c(
          #"t0min_1mM",
          "t15min_10mM","t45min_10mM","t90min_10mM","t180min_10mM","t1200min_10mM",
          "t15min_1mM","t45min_1mM","t90min_1mM","t180min_1mM","t1200min_1mM")


makeHeatMap <- function(dataWide,sub=NULL,newOrder=NULL,hmPalette,show_rownames = F){
  if (is.null(sub)) {sub <- rownames(dataWide)}
  if (is.null(newOrder)) {newOrder <- colnames(dataWide)}
  hmDat <- dataWide[sub,newOrder]
  #hmPalette <- colorRampPalette(brewer.pal(name = "PuOr",9))(500)
  pheatmap::pheatmap(hmDat, cluster_cols = F, show_rownames = show_rownames,
                     cellwidth = 8, fontsize_row = 6, gaps_col = c(1,6),
                     border_color = NA,
                     color = hmPalette)
}


makeLinePlots <- function(dataLong,sub=NULL,title="",nrow = 4){
  if (is.null(sub)) {sub <- rownames(dataWide)}
  
  dataLong %>% filter(Gene %in% sub) %>%
    mutate(Time=factor(gsub("t|min","",Time),levels = gtools::mixedsort( unique(gsub("t|min","",Time))))) %>%
    #mutate(Time=as.numeric(gsub("t|min","",Time))) %>%
    mutate("tmp"=paste0(Gene," (",Name,")")) %>%
    ggplot(aes(x=Time,y=meanExpr)) +
    geom_point(aes(color=NO3,group=NO3,shape=NO3) ) +
    geom_line(aes(color=NO3,group=NO3) ) +
    theme(aspect.ratio=16/9) +
    ggtitle(title) +
    facet_wrap(as.factor(tmp)~.,scales = "free_x",nrow = nrow) +
    theme_classic()
}


# 
# pltMe <- "AT1G65483"
# expressionData_long %>% filter(Gene %in% pltMe) %>%
#   mutate("t"=factor(gsub("t|min","",Time),levels=c("15","45","90","180","1200"))) %>%
#   ggplot(aes(x=t,y=meanExpr)) +
#   ggtitle(t2gtair10[t2gtair10$Gene==pltMe,"Name"]) +
#   geom_line(aes(color=NO3,group=NO3)) +
#   geom_point(aes(color=NO3)) +
#   facet_grid(Name~.,scales = "free_y") +
#   theme_minimal()
# 
