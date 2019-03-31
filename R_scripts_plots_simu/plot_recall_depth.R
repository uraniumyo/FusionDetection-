suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
})

##### recall vs sequencing depth at coverages : 5,10,20,30 ####

# set colors
colorlist <- brewer.pal(7, "Set1")[c(1,2,3)]
cols <- c("last"=colorlist[2],"minimap2"=colorlist[1],"ngmlr"=colorlist[3])

# read data

dir_pb <- '/Users/yaoyao/Desktop/FusionDetection/R_scripts_plots_simu/recall_vs_depth_pb.txt'
dir_ont <- '/Users/yaoyao/Desktop/FusionDetection/R_scripts_plots_simu/recall_vs_depth_ont.txt'

reformat_recall<-  function(dir){
  recall <- read.table(dir,sep = '',header=FALSE)
  rownames(recall)<-recall[,1] 
  recall <-recall[,-1] 
  recall <- as.data.frame(t(recall))
  return(recall)
}
 
recall_pb <-  reformat_recall(dir_pb)
recall_ont <-  reformat_recall(dir_ont)

# recall vs sequencing depth pacbio 
ggplot(recall_pb)+theme_bw()+
  geom_point(aes(x=coverage,y=ngmlr,color='ngmlr'))+
  geom_line(aes(x=coverage,y=ngmlr,color='ngmlr'),size=1.1)+
  geom_point(aes(x=coverage,y=minimap2,color='minimap2'))+
  geom_line(aes(x=coverage,y=minimap2,color='minimap2'),size=1.1)+
  geom_point(aes(x=coverage,y=last,color='last'))+
  geom_line(aes(x=coverage,y=last,color='last'),size=1.1)+
  ylab('Recall %')+
  #scale_x_continuous(limits=c(5, 30))+ 
  scale_x_discrete(limits=c(5,10,20,30))+
  scale_y_continuous(limits = c(0,100))+
  theme(legend.background = element_rect(fill="white",size=.2, color='black',linetype="solid"))+
  scale_color_manual(name  ="aligners",values = cols,labels=c("last","minimap2","ngmlr"))


# recall vs sequencing depth nanopore
ggplot(recall_ont)+theme_bw()+
  geom_point(aes(x=coverage,y=ngmlr,color='ngmlr'))+
  geom_line(aes(x=coverage,y=ngmlr,color='ngmlr'),size=1.1)+
  geom_point(aes(x=coverage,y=minimap2,color='minimap2'))+
  geom_line(aes(x=coverage,y=minimap2,color='minimap2'),size=1.1)+
  geom_point(aes(x=coverage,y=last,color='last'))+
  geom_line(aes(x=coverage,y=last,color='last'),size=1.1)+
  ylab('Recall %')+
  #scale_x_continuous(limits=c(5, 30))+ 
  scale_x_discrete(limits=c(5,10,20,30))+
  scale_y_continuous(limits = c(0,100))+
  theme(legend.background = element_rect(fill="white",size=.2, color='black',linetype="solid"))+
  scale_color_manual(name  ="aligners",values = cols,labels=c("last","minimap2","ngmlr"))




