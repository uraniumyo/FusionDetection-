suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(Biostrings)
  library(VennDiagram)
  library(reshape2)
  library(ggthemes)
  library(grid)
  library(scales) 
  library(RColorBrewer)
  library(ROCR)
  library(caTools)
})

######################### TRUE fusion result #######################

# read the actual fusion info
fusion_truth <- read.csv("/Users/yaoyao/Desktop/FusionDetection/convert_name/fusion_summary.csv", header = TRUE,sep='')
# add breakpoint to end distance (0 end, 1 middle) ##min doesnot work with mutate, use pmin
fusion_truth <- fusion_truth %>% mutate(fustxlen=(gene_A_length+gene_B_length),
                                        dist = 2*pmin(gene_A_length,gene_B_length)/fustxlen)

# plot: length distribution of simulated fusion transcripts
dftxlen <- data.frame(name='simulated_transcript',length=fusion_truth$fustxlen)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
max(dftxlen$length) #44747
min(dftxlen$length) #221
median(dftxlen$length)#1811.5
mean(dftxlen$length) #2452.27

# only 2 reads over 20000bp 
'
ggplot(dftxlen,aes(x=name,length))+geom_violin()+scale_y_continuous(limits = c(0,20000))+
  stat_summary(fun.data=data_summary)+
  theme(axis.title.x = element_blank())+ylab('length(bp)')
'
# plot: fusion transcript length distribution in histgram
ggplot(dftxlen,aes(x=length))+geom_histogram(binwidth=200,alpha=0.85)+
  #scale_x_continuous(limits = c(0,40000))+
  xlab('length(bp)')+
  geom_vline(aes(xintercept=mean(length)),
              color="blue", linetype="dashed", size=0.5,alpha=0.7)+
  geom_rug(aes(x=length),alpha=0.15)+theme_bw() 

########################### set color for aligners #################

colorlist <- brewer.pal(7, "Set1")[c(1,2,3)]
cols <- c("last"=colorlist[2],"minimap2"=colorlist[1],"ngmlr"=colorlist[3])


######################## read fusion files ##########################

# read TPFP results at covergae of 30
ngmlrtpfp1 <- read.table('/Users/yaoyao/Desktop/testforpython/ngmlr_pb30.txt')
mnp2tpfp1 <- read.table('/Users/yaoyao/Desktop/testforpython/mnp2_pb30.txt')
lasttpfp1 <- read.table('/Users/yaoyao/Desktop/testforpython/last_pb30.txt')

# number of reads aligned to at least 2 different genes
#nrow(ngmlrtpfp1) #1241
#nrow(mnp2tpfp1) #4198
#nrow(lasttpfp1) #10759

# remove reads containing NULLs (genes without annotation)
ngmlrtpfp <- ngmlrtpfp1 %>% filter(V2 != "NULL" & V3 != "NULL"& !grepl('NULL',V10))
mnp2tpfp <- mnp2tpfp1 %>% filter(V2 != "NULL" & V3 != "NULL" & !grepl('NULL',V10))
lasttpfp <- lasttpfp1 %>% filter(V2 != "NULL" & V3 != "NULL" & !grepl('NULL',V10))

# number of reads - all alignments within annotation ranges 
# number and % of alignments with NULL(out of gene annotation range/map to contigs without annotation)
#nrow(ngmlrtpfp) 1026 #215 17.32%
#nrow(mnp2tpfp) 3687 #511 12.17%
#nrow(lasttpfp) 9543 #1216 11.30%

# porprotion of reads with annotation for all alignments
a <- vector()
a[1] <- nrow(ngmlrtpfp)/nrow(ngmlrtpfp1) 
a[2] <- nrow(mnp2tpfp)/nrow(mnp2tpfp1) 
a[3] <- nrow(lasttpfp)/nrow(lasttpfp1)
dfwithanno <-  data.frame(aligners=c('ngmlr','minimap2','last'), anno=a)

# plot: proportion of reads with all alignemnts strictly fall in annotations
ggplot(dfwithanno,aes(aligners,anno,fill=aligners))+theme_bw()+
  geom_bar(width=0.65,stat="identity",position="dodge",color='black',alpha=0.8)+
  scale_fill_manual(values=cols)+
  theme(legend.background = element_rect(fill="white",size=.4, color='black',linetype="solid"))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits = c(0,1))+
  theme(axis.title.y = element_text(size=14))+
  ylab('porprotion of reads')+
  #ggtitle("porprotion of reads with annotation for all alignments")+
  scale_color_manual(name  ="aligners",values = cols,labels=c("last","minimap2","ngmlr"))+
  geom_text(aes(y=anno, label = paste(round(100*anno, 2), "%", sep="")),vjust =2, colour = "white", position = position_dodge(.9),size=5)



######################### TP/FP results #############################

# number of fusion genes output 
length(unique(ngmlrtpfp$V10)) #201
length(unique(mnp2tpfp$V10)) #473
length(unique(lasttpfp$V10)) #1189

# number of TRUE fusion genes output
length(unique((ngmlrtpfp[ngmlrtpfp$V7=='TP',])$V10)) #104 51.74%
length(unique((mnp2tpfp[mnp2tpfp$V7=='TP',])$V10)) #318 67.23%
length(unique((lasttpfp[lasttpfp$V7=='TP',])$V10)) #353 29.69%

# number of FALSE fusion genes output
length(unique((ngmlrtpfp[ngmlrtpfp$V7 =='FP',])$V10)) #97 48.26%
length(unique((mnp2tpfp[mnp2tpfp$V7 =='FP',])$V10)) #155 32.77%
length(unique((lasttpfp[lasttpfp$V7 =='FP',])$V10)) #836 70.31%


### number of supporting reads ###
# i - threshold for number of supporting reads
# for detected fusions with freq lower than i, remove from TP, considered as FP
# precision here = TP/(TP+FP) = number of true detected fusions/ all detected fusions above threshold
cal_precision <- function(u,i){
  true <- unique(u[u$V7=='TP',]$V10)
  Freq <- plyr::count(u$V10)
  freq_rm <- Freq %>% filter(Freq$freq>=i)
  precision <- sum(freq_rm$x %in% true)/nrow(freq_rm)
  precision[is.na(precision)]<-0
  precision
}

a=vector()
b=vector()
c=vector()
index=vector()

for (i in 1:31){
  index[i] <- i
  a[i] <- cal_precision(ngmlrtpfp,i)
  b[i] <- cal_precision(mnp2tpfp,i)
  c[i] <- cal_precision(lasttpfp,i)
}
dfprecision <- data.frame(coverage=index,ngmlr=a,minimap2=b,last=c)

# plot: number of supporting reads vs precision of detected fusions
ggplot(dfprecision)+theme_bw()+
  geom_point(aes(x=coverage,y=ngmlr,color='ngmlr'))+
  geom_line(aes(x=coverage,y=ngmlr,color='ngmlr'))+
  geom_point(aes(x=coverage,y=minimap2,color='minimap2'))+
  geom_line(aes(x=coverage,y=minimap2,color='minimap2'))+
  geom_point(aes(x=coverage,y=last,color='last'))+
  geom_line(aes(x=coverage,y=last,color='last'))+
  xlab('number of supporting reads')+ylab('precision')+
  scale_x_continuous(breaks=seq(0, 30, 2))+scale_y_continuous(limits = c(0,1),breaks=seq(0,1,0.1))+
  theme(legend.background = element_rect(fill="white",size=.2, color='black',linetype="solid"))+
  scale_color_manual(name  ="aligners",values = cols,labels=c("last","minimap2","ngmlr"))
  #ggtitle("Precision with increasing number of supporting reads as threshold")

# i - threshold for number of supporting reads
# for detected fusions with freq lower than i, remove from TP, considered as FP
# recall here = TP/(TP+FN)= number of detected fusions/ 500 simulated fusions
cal_recall <- function(u,i){
  true <- unique(u[u$V7=='TP',]$V10)
  Freq <- plyr::count(u$V10)
  freq_rm <- Freq %>% filter(Freq$freq>=i)
  recall <- sum(freq_rm$x %in% true)/500
  recall
}

d=vector()
e=vector()
f=vector()
index2=vector()
for (i in 1:31){
  index2[i] <- i
  d[i] <- cal_recall(ngmlrtpfp,i)
  e[i] <- cal_recall(mnp2tpfp,i)
  f[i] <- cal_recall(lasttpfp,i)
}
dfrecall <- data.frame(coverage=index2,ngmlr=d,minimap2=e,last=f)

# plot: number of supporting reads vs recall of detected fusions
ggplot(dfrecall)+theme_bw()+
  geom_point(aes(x=coverage,y=ngmlr,color='ngmlr'))+
  geom_line(aes(x=coverage,y=ngmlr,color='ngmlr'))+
  geom_point(aes(x=coverage,y=minimap2,color='minimap2'))+
  geom_line(aes(x=coverage,y=minimap2,color='minimap2'))+
  geom_point(aes(x=coverage,y=last,color='last'))+
  geom_line(aes(x=coverage,y=last,color='last'))+
  xlab('number of supporting reads')+ylab('recall')+
  scale_x_continuous(breaks=seq(0, 30, 2))+scale_y_continuous(limits = c(0,0.8),breaks=seq(0,0.8,0.1))+
  theme(legend.background = element_rect(fill="white",size=.2, color='black',linetype="solid"))+
  scale_color_manual(name  ="aligners",values = cols,labels=c("last","minimap2","ngmlr"))
  #ggtitle("Recall with increasing number of supporting reads as threshold")


# get true fusion name list 
get_fusion <- function(u){
  tmp <- list()
  for (i in 1:nrow(u)){ tmp[i] <- paste(u[i,],collapse = "_")} 
  tmp}
true <- get_fusion(as.matrix(fusion_truth[17:18]))
true <- data.frame(x=unlist(true))

# calculate performance for prcurve and auc value
cal_PR_curve <- function(u,true){
  result <- plyr::count(u$V10)
  tmp <-  full_join(result,true,by='x')
  tmp[is.na(tmp$freq),]$freq <- -1
  freq <- tmp$freq
  label <- tmp$x %in% true$x 
  pred <- prediction(freq,label)
  perf <- performance(pred,"prec", "rec")
  # remove the last result when cut off below 0 and recall=1 
  perf@x.values[[1]][length(perf@x.values[[1]])] <- NaN
  perf@y.values[[1]][length(perf@y.values[[1]])] <- NaN
  perf@alpha.values[[1]][length(perf@alpha.values[[1]])] <- NaN
  #save auc value(roc..)
  auc <-  performance(pred,"auc")@y.values
  return(list(perf=perf,AUC=auc))
}

ngmlrpr <- cal_PR_curve(ngmlrtpfp,true)
mnppr <- cal_PR_curve(mnp2tpfp,true)
lastpr <- cal_PR_curve(lasttpfp,true)


# manually calculate PRcurve results for debug
dfPRcurve <- data.frame(coverage=dfrecall$coverage,
                        ngmlr.p=dfprecision$ngmlr,ngmlr.r=dfrecall$ngmlr,
                        minimap2.p=dfprecision$minimap2,minimap2.r=dfrecall$minimap2,
                        last.p=dfprecision$last,last.r=dfrecall$last)
# Plot: Precision-Recall curve 
ggplot(dfPRcurve)+theme_bw()+
  geom_point(aes(x=ngmlr.r,y=ngmlr.p,color='ngmlr'), size=0.5)+geom_line(aes(x=ngmlr.r,y=ngmlr.p,color='ngmlr'), size=0.7)+
  geom_point(aes(x=minimap2.r,y=minimap2.p,color='minimap2'), size=0.5)+geom_line(aes(x=minimap2.r,y=minimap2.p,color='minimap2'), size=0.7)+
  geom_point(aes(x=last.r,y=last.p,color='last'), size=0.5)+geom_line(aes(x=last.r,y=last.p,color='last'), size=0.7)+
  ylab('Precision')+xlab('Recall')+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,1))+
  theme(legend.background = element_rect(fill="white",size=.2, color='black',linetype="solid"))+
  scale_color_manual(name  ="aligners",values = cols,labels=c("last","minimap2","ngmlr"))


# calculate AUC values
ngmlr_auc <- trapz(rev(dfrecall$ngmlr),rev(dfprecision$ngmlr))
mnp2_auc <- trapz(rev(dfrecall$minimap2),rev(dfprecision$minimap2))
last_auc <- trapz(rev(dfrecall$last),rev(dfprecision$last))


# Plot: AUC values 
dfauc <- data.frame(aligners=c('ngmlr','minimap2','last'),auc=c(ngmlr_auc,mnp2_auc,last_auc))

ggplot(dfauc,aes(aligners,auc,fill=aligners))+theme_bw()+
  geom_bar(stat="identity",position="dodge",color='black',alpha=0.8)+
  scale_fill_manual(values=cols)+
  theme(legend.position=c(0.85,0.75))+
  theme(legend.background = element_rect(fill="white",size=.4, color='black',linetype="solid"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_blank())+
  scale_y_continuous(limits = c(0,0.8))+
  theme(axis.title.y = element_text(size=14))+
  #ggtitle("AUC values")+
  ylab('AUC values')+
  geom_text(aes(y=auc, label = paste(round(auc, 4))),vjust =2, colour = "white", position = position_dodge(.9),size=5)

# F1 accuracy measure
# F1 score = 2*precision*recall/(precision+recall)
dfF1<- dfPRcurve %>%mutate(ngmlr.f1=2*ngmlr.p*ngmlr.r/(ngmlr.p+ngmlr.r),
                           minimap2.f1=2*minimap2.p*minimap2.r/(minimap2.p+minimap2.r),
                           last.f1=2*last.p*last.r/(last.p+last.r))
###Plot: F1 score versus number of supporting reads
ggplot(dfF1)+theme_bw()+
  geom_point(aes(x=coverage,y=ngmlr.f1,color='ngmlr'))+
  geom_line(aes(x=coverage,y=ngmlr.f1,color='ngmlr'))+
  geom_point(aes(x=coverage,y=minimap2.f1,color='minimap2'))+
  geom_line(aes(x=coverage,y=minimap2.f1,color='minimap2'))+
  geom_point(aes(x=coverage,y=last.f1,color='last'))+
  geom_line(aes(x=coverage,y=last.f1,color='last'))+
  xlab('number of supporting reads')+ylab('F1 score')+
  scale_x_continuous(breaks=seq(0, 30, 2))+scale_y_continuous(limits = c(0,0.8),breaks=seq(0,0.8,0.1))+
  theme(legend.background = element_rect(fill="white",size=.2, color='black',linetype="solid"))+
  scale_color_manual(name  ="aligners",values = cols,labels=c("last","minimap2","ngmlr"))
  #ggtitle("F1 score with increasing number of supporting reads as threshold")



################## effects of allowing duplicated gene families ##################

# number of TP reads
nrow(ngmlrtpfp[ngmlrtpfp$V7=='TP',]) #614
nrow(mnp2tpfp[mnp2tpfp$V7=='TP',]) #2786
nrow(lasttpfp[lasttpfp$V7=='TP',]) #5612

# number of FP reads
nrow(ngmlrtpfp[ngmlrtpfp$V7=='FP',]) #412
nrow(mnp2tpfp[mnp2tpfp$V7=='FP',]) #901
nrow(lasttpfp[lasttpfp$V7=='FP',]) #3931

# number of TP reads [allowing duplicated gene families]
nrow(ngmlrtpfp[ngmlrtpfp$V9=='TP',]) #616
nrow(mnp2tpfp[mnp2tpfp$V9=='TP',]) #2809
nrow(lasttpfp[lasttpfp$V9=='TP',]) #5970

# proportion of reads that are correctly aligned (TP/TP+FP) (precision in reads not fusions)
nrow(ngmlrtpfp[ngmlrtpfp$V7=='TP',])/nrow(ngmlrtpfp) #0.5984405
nrow(mnp2tpfp[mnp2tpfp$V7=='TP',])/nrow(mnp2tpfp) #0.7556279
nrow(lasttpfp[lasttpfp$V7=='TP',])/nrow(lasttpfp) #0.588075

# proportion of reads that are correctly aligned [allowing duplicated gene families]
nrow(ngmlrtpfp[ngmlrtpfp$V9=='TP',])/nrow(ngmlrtpfp) #0.6003899
nrow(mnp2tpfp[mnp2tpfp$V9=='TP',])/nrow(mnp2tpfp) #0.761866
nrow(lasttpfp[lasttpfp$V9=='TP',])/nrow(lasttpfp) #0.6255894


# TP FP differences(number of reads) when considering duplicated gene families
tp <- vector()
fp <- vector()
tp[1]<- nrow(ngmlrtpfp[ngmlrtpfp$V7=='TP',])
tp[2]<- nrow(ngmlrtpfp[ngmlrtpfp$V9=='TP',]) 
tp[3]<- nrow(mnp2tpfp[mnp2tpfp$V7=='TP',])
tp[4]<- nrow(mnp2tpfp[mnp2tpfp$V9=='TP',])
tp[5]<- nrow(lasttpfp[lasttpfp$V7=='TP',])
tp[6]<- nrow(lasttpfp[lasttpfp$V9=='TP',])

fp[1]<- nrow(ngmlrtpfp[ngmlrtpfp$V7=='FP',])
fp[2]<- nrow(ngmlrtpfp[ngmlrtpfp$V9=='FP',]) 
fp[3]<- nrow(mnp2tpfp[mnp2tpfp$V7=='FP',])
fp[4]<- nrow(mnp2tpfp[mnp2tpfp$V9=='FP',])
fp[5]<- nrow(lasttpfp[lasttpfp$V7=='FP',])
fp[6]<- nrow(lasttpfp[lasttpfp$V9=='FP',])
grp <- rep(c('original','gene_family'),3)

dfduptpfp <-data.frame(aligners=c('ngmlr','ngmlr','minimap2','minimap2','last','last'),
                       group=grp,
                       TP=tp,FP=fp) 
# plot: arrows, TP FP differences(number of reads) when considering duplicated gene families
ggplot(dfduptpfp)+theme_bw()+
  geom_point(aes(FP,TP,color=aligners))+
  #scale_y_continuous(trans=log2_trans())+ scale_x_continuous(trans=log2_trans())+
  #scale_y_log10()+scale_x_log10()+
  geom_path(aes(x=FP, y=TP, group = aligners,color=aligners), size=1.3,
            arrow = arrow(length = unit(0.35, "cm"),ends='last',type='open'))+
  scale_color_manual(name  ="aligners",values = cols,labels=c("last","minimap2","ngmlr"))
  #ggtitle("Changes in number of TP/FP reads when considering duplicated gene families")


# plot: boxplot, changes in precision when considering duplicated gene families
dfduptpfp <- dfduptpfp %>% mutate(precision=TP/(TP+FP))
dfdup <-data.frame(aligners=dfduptpfp$aligners,group=dfduptpfp$group,precision=dfduptpfp$precision)

ggplot(dfdup,aes(aligners,precision,fill=group))+theme_bw()+
  geom_bar(stat="identity",position="dodge",color='black')+
  scale_fill_manual(values=c("#FB8072", "#80B1D3"))+
  theme(legend.position=c(0.85,0.885))+
  theme(legend.background = element_rect(fill="white",size=.4, color='black',linetype="solid"))+
  theme(axis.text.x = element_text(size = 13))+theme(axis.title.x = element_blank())+
  scale_y_continuous(limits = c(0,0.8))+
  theme(axis.title.y = element_text(size=14))+
  #ggtitle("Changes in precision when considering duplicated gene families")+
  geom_text(aes(y=precision, label = paste(round(100*precision, 2), "%", sep="")),vjust =2, colour = "white", position = position_dodge(.9),size=5)

 

##################### venn diagrams of detected fusions #######################


# number of TP fusions 
length(unique((ngmlrtpfp[ngmlrtpfp$V7 =='TP',])$V10)) #104
length(unique((mnp2tpfp[mnp2tpfp$V7=='TP',])$V10)) #318
length(unique((lasttpfp[lasttpfp$V7 =='TP',])$V10)) #353

# number of FP fusions 
length(unique((ngmlrtpfp[ngmlrtpfp$V7 =='FP',])$V10)) #97
length(unique((mnp2tpfp[mnp2tpfp$V7=='FP',])$V10)) #155
length(unique((lasttpfp[lasttpfp$V7 =='FP',])$V10)) #836


# venn diagram of all
vennall <-list(minimap2=unique(mnp2tpfp$V10),ngmlr=unique(ngmlrtpfp$V10),last=unique(lasttpfp$V10))
venn.diagram(vennall,filename='/Users/yaoyao/Desktop/vennall.tiff',
             #lty = rep("solid", 3),lwd = 2,
             main="Intersection of all detected gene fusions",
             main.cex = 1.7, main.fontface = "plain",
             col = "transparent",fill = c(colorlist[1], colorlist[3], colorlist[2]),
             alpha = 0.60,cex = 1.5,fontfamily = "serif",fontface = "bold",
             lty = 2, cat.cex = 1.2,
             cat.col = c("red", "darkgreen", "darkblue"),cat.fontface="bold",margin = 0.04)

# venn diagram of TP
vennTP <- list(minimap2 = unique((mnp2tpfp[mnp2tpfp$V7=='TP',])$V10),
               ngmlr = unique((ngmlrtpfp[ngmlrtpfp$V7 =='TP',])$V10),
               last = unique((lasttpfp[lasttpfp$V7 =='TP',])$V10))
venn.diagram(vennTP,filename='/Users/yaoyao/Desktop/vennTP.tiff',
             #lty = rep("solid", 3),lwd = 2,
             main="Intersection of true discoveries",
             main.cex = 1.7, main.fontface = "plain",
             col = "transparent",fill = c(colorlist[1], colorlist[3], colorlist[2]),
             alpha = 0.60,cex = 1.5,fontfamily = "serif",fontface = "bold",
             lty = 2, cat.cex = 1.2,
             cat.col = c("red", "darkgreen", "darkblue"),cat.fontface="bold",margin = 0.04)


# venn diagram of FP
vennFP <- list(minimap2 = unique((mnp2tpfp[mnp2tpfp$V7=='FP',])$V10),
               ngmlr = unique((ngmlrtpfp[ngmlrtpfp$V7 =='FP',])$V10),
               last = unique((lasttpfp[lasttpfp$V7 =='FP',])$V10))
venn.diagram(vennFP,filename='/Users/yaoyao/Desktop/vennFP.tiff',
             #lty = rep("solid", 3),lwd = 2,
             main="Intersection of false discoveries",
             main.cex = 1.7, main.fontface = "plain",
             col = "transparent",fill = c(colorlist[1], colorlist[3], colorlist[2]),
             alpha = 0.60,cex = 1.5,fontfamily = "serif",fontface = "bold",
             lty = 2, cat.cex = 1.2,
             cat.col = c("red", "darkgreen", "darkblue"),cat.fontface="bold",margin = 0.04)



################ breakpoint analysis #############################

cal_true_break <- function(u){
  utp <- u %>% filter(u$V7=='TP')
  len <- nrow(utp)
  a <- vector()
  for (i in 1:21){
   a[i] <- sum(abs(u[u$V7=='TP',]$V5) < i & abs(u[u$V7=='TP',]$V4) < i)/len
   }
  return(a)
}

a<- cal_true_break(ngmlrtpfp)
b<- cal_true_break(mnp2tpfp)
c<- cal_true_break(lasttpfp)
dfbreakpoint <- data.frame(distance=c(0:20),ngmlr=a,minimap2=b,last=c)

# plot: increase distance range vs accuracy of breakpoint detection 
ggplot(dfbreakpoint)+theme_bw()+
  geom_point(aes(x=distance,y=ngmlr,color='ngmlr'))+
  geom_line(aes(x=distance,y=ngmlr,color='ngmlr'))+
  geom_point(aes(x=distance,y=minimap2,color='minimap2'))+
  geom_line(aes(x=distance,y=minimap2,color='minimap2'))+
  geom_point(aes(x=distance,y=last,color='last'))+
  geom_line(aes(x=distance,y=last,color='last'))+
  xlab('allowed distance to true breakpoint(bp)')+ylab('accuracy')+
  scale_x_continuous(breaks=seq(0, 30, 2))+scale_y_continuous(limits = c(0,1),breaks=seq(0,1,0.1))+
  theme(legend.background = element_rect(fill="white",size=.2, color='black',linetype="solid"))+
  scale_color_manual(name  ="aligners",values = cols,labels=c("last","minimap2","ngmlr"))
  #ggtitle("Breakpoint accuracy with increasing distance ranges")


### breakpoint accuracy when set threshold to 20bp
bp20 <- vector() 
bp20[1] <- nrow(ngmlrtpfp[ngmlrtpfp$V8=='TP',])/nrow(ngmlrtpfp[ngmlrtpfp$V7=='TP',]) #41.21%
bp20[2]<- nrow(mnp2tpfp[mnp2tpfp$V8=='TP',])/nrow(mnp2tpfp[mnp2tpfp$V7=='TP',]) #58.00%
bp20[3] <- nrow(lasttpfp[lasttpfp$V8=='TP',])/nrow(lasttpfp[lasttpfp$V7=='TP',]) #76.62%

dfbreak20bp <- data.frame(aligners=c('ngmlr','minimap2','last'),accuracy=bp20)
# plot: accuracy of breakpoint detection when dist within 20bp
ggplot(dfbreak20bp,aes(aligners,accuracy,fill=aligners))+theme_bw()+
  geom_bar(stat="identity",position="dodge",color='black',alpha=0.8)+
  scale_fill_manual(values=cols)+
  theme(legend.position=c(0.8,0.8))+
  theme(legend.background = element_rect(fill="white",size=.4, color='black',linetype="solid"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(size=14))+
  scale_y_continuous(limits = c(0,1),breaks=seq(0,1,0.1))+
  ylab('accuracy')+
  #ggtitle("Breakpoint accuracy when allow distance range at 20bp")+
  geom_text(aes(y=accuracy, label = paste(round(100*accuracy, 2), "%", sep="")),vjust =2, colour = "white", position = position_dodge(.9),size=5)


### breakexons(BE) vs intactexons(IE) ###
BE <- as.matrix(fusion_truth[fusion_truth$type=='break',17:18])
IE <- as.matrix(fusion_truth[fusion_truth$type=='intact',17:18])

get_fusion <- function(u){
tmp <- list()
for (i in 1:nrow(u)){ tmp[i] <- paste(u[i,],collapse = "_")} 
tmp}

BEfusion <- get_fusion(BE)
IEfusion <- get_fusion(IE)

'''
nrow(ngmlrtpfp[ngmlrtpfp$V8=='TP'&ngmlrtpfp$V10 %in% IEfusion,])/nrow(ngmlrtpfp[ngmlrtpfp$V10 %in% IEfusion,])
nrow(ngmlrtpfp[ngmlrtpfp$V8=='TP'&ngmlrtpfp$V10 %in% BEfusion,])/nrow(ngmlrtpfp[ngmlrtpfp$V10 %in% BEfusion,]) #24/95
nrow(mnp2tpfp[mnp2tpfp$V8=='TP'&mnp2tpfp$V10 %in% IEfusion,])/nrow(mnp2tpfp[mnp2tpfp$V10 %in% IEfusion,])
nrow(mnp2tpfp[mnp2tpfp$V8=='TP'&mnp2tpfp$V10 %in% BEfusion,])/nrow(mnp2tpfp[mnp2tpfp$V10 %in% BEfusion,])
nrow(lasttpfp[lasttpfp$V8=='TP'&lasttpfp$V10 %in% IEfusion,])/nrow(lasttpfp[lasttpfp$V10 %in% IEfusion,])
nrow(lasttpfp[lasttpfp$V8=='TP'&lasttpfp$V10 %in% BEfusion,])/nrow(lasttpfp[lasttpfp$V10 %in% BEfusion,])
'''
# number of correct detections of breakpoint/ number of correct mapped reads in IE/BE type 
breakpoint_IE_percent <- function(u,v,w){
  return(list(ngmlr=nrow(u[u$V8=='TP'&u$V10 %in% IEfusion,])/nrow(u[u$V10 %in% IEfusion,]),
              minimap2= nrow(v[v$V8=='TP'&v$V10 %in% IEfusion,])/nrow(v[v$V10 %in% IEfusion,]),
              last = nrow(w[w$V8=='TP'&w$V10 %in% IEfusion,])/nrow(w[w$V10 %in% IEfusion,])))}
breakpoint_BE_percent <- function(u,v,w){
  return(list(ngmlr=nrow(u[u$V8=='TP'&u$V10 %in% BEfusion,])/nrow(u[u$V10 %in% BEfusion,]),
              minimap2= nrow(v[v$V8=='TP'&v$V10 %in% BEfusion,])/nrow(v[v$V10 %in% BEfusion,]),
              last = nrow(w[w$V8=='TP'&w$V10 %in% BEfusion,])/nrow(w[w$V10 %in% BEfusion,])))}
ie <- breakpoint_IE_percent(ngmlrtpfp,mnp2tpfp,lasttpfp)
be <- breakpoint_BE_percent(ngmlrtpfp,mnp2tpfp,lasttpfp)
dfIEBE <- data.frame(aligners=c('ngmlr','minimap2','last'),IE=unlist(ie),BE=unlist(be))

# plot: accuracy of breakpoint detection in IE/BE type fusions
dfIEBE <- melt(dfIEBE,id.vars="aligners",variable.name="type",value.name="accuracy")

ggplot(dfIEBE,aes(aligners,accuracy,fill=type))+theme_bw()+
  geom_bar(stat="identity",position="dodge",color='black')+
  scale_fill_manual(values=c("#FB8072", "#80B1D3"))+
  theme(legend.position=c(0.8,0.8))+
  theme(legend.background = element_rect(fill="white",size=.4, color='black',linetype="solid"))+
  theme(axis.text.x = element_text(size = 13))+theme(axis.title.x = element_blank())+
  #scale_y_continuous(limits = c(0,1))+
  expand_limits(y=c(0,1))+ 
  theme(axis.title.y = element_text(size=14))+
  #ggtitle("Accuracy of breakpoint detection in IE/BE type fusions")
  geom_text(aes(y=accuracy, label = paste(round(100*accuracy, 2), "%", sep="")),vjust =2, colour = "white", position = position_dodge(.9),size=5)


'''
# maybe violin plot for distances ?
# number of reads different& distrange too large[0,40k]...diffcult to check on violin plot
get_dist <-function(u,name){
  C1 <- data.frame(dist = abs(u[u$V7=='TP',][,4]),aligner = name)
  C2 <- data.frame(dist = abs(u[u$V7=='TP',][,5]),aligner = name)
  dat <- rbind(C1, C2)
  return(dat)
}
a <- get_dist(ngmlrtpfp,'ngmlr')
b <- get_dist(mnp2tpfp,'minimap2')
c <-  get_dist(mnp2tpfp,'last')
dat <- rbind(a,b,c)
ggplot(dat) + geom_violin(aes(x=aligner,y=dist))+scale_y_continuous(limits = c(10000,20000))
'''


### proportion of true discoveries but possibly missing one/more exons 
# median length of human intron 1334bp
# number of reads missing exon/ number of correctly mapped reads
intronlen <- 1000
introndist <- vector()
#136/614
introndist[1] <- sum(ngmlrtpfp[ngmlrtpfp$V7=='TP',]$V4 > intronlen|ngmlrtpfp[ngmlrtpfp$V7=='TP',]$V5 > intronlen)/nrow(ngmlrtpfp[ngmlrtpfp$V7=='TP',])
#483/2786
introndist[2] <- sum(mnp2tpfp[mnp2tpfp$V7=='TP',]$V4 > intronlen|mnp2tpfp[mnp2tpfp$V7=='TP',]$V5 > intronlen)/nrow(mnp2tpfp[mnp2tpfp$V7=='TP',])
#416/5612
introndist[3] <- sum(lasttpfp[lasttpfp$V7=='TP',]$V4 > intronlen|lasttpfp[lasttpfp$V7=='TP',]$V5 > intronlen)/nrow(lasttpfp[lasttpfp$V7=='TP',])
dfmissexon <- data.frame(aligners=c('ngmlr','minimap2','last'),missexon = introndist)

ggplot(dfmissexon,aes(aligners,missexon,fill=aligners))+theme_bw()+
  geom_bar(stat="identity",position="dodge",color='black',alpha=0.8)+
  scale_fill_manual(values=cols)+
  theme(legend.position=c(0.8,0.8))+
  theme(legend.background = element_rect(fill="white",size=.4, color='black',linetype="solid"))+
  theme(axis.text.x =  element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(size=14))+
  expand_limits(y=c(0,0.5))+
  ylab('proportion of reads')+
  #ggtitle("Proportion of reads missing one or more exons in alignments")+
  geom_text(aes(y=missexon, label = paste(round(100*missexon, 2), "%", sep="")),vjust =2, colour = "white", position = position_dodge(.9),size=5)


# relationship between breakpoint location vs number of detected reads (0 close to ends; 1 in the middle)
dfcountngmlr <- plyr::count(ngmlrtpfp$V6)[-nrow(plyr::count(ngmlrtpfp$V6)),]
dfcountmnp2 <- plyr::count(mnp2tpfp$V6)[-nrow(plyr::count(mnp2tpfp$V6)),]
dfcountlast <- plyr::count(lasttpfp$V6)[-nrow(plyr::count(lasttpfp$V6)),]

g <- ggplot() + theme_bw()+
  geom_point(aes(x=dfcountngmlr$x,y=dfcountngmlr$freq,color="ngmlr"),alpha=0.35)+
  stat_smooth(method = 'loess',aes(x=dfcountngmlr$x,y=dfcountngmlr$freq,color="ngmlr"),alpha=0.1,se = FALSE,size=0.4)+
  geom_point(aes(x=dfcountmnp2$x,y=dfcountmnp2$freq,color='minimap2'),alpha=0.35)+
  stat_smooth(method = 'loess',aes(x=dfcountmnp2$x,y=dfcountmnp2$freq,color='minimap2'),alpha=0.1,se = FALSE,size=0.4)+
  geom_point(aes(x=dfcountlast$x,y=dfcountlast$freq,color='last'),alpha=0.35)+
  stat_smooth(method = 'loess',aes(x=dfcountlast$x,y=dfcountlast$freq,color='last'),alpha=0.1,se = FALSE,size=0.4)+
  #labs(title = "relation of breakpoint location and number of detected reads")+
  xlab('distance to read ends')+ylab('number of supporting reads')+
  stat_smooth(method = 'loess',aes(x=dfcountngmlr$x,y=dfcountngmlr$freq,color="ngmlr"),alpha=0.1,se = FALSE,size=0.4)+
  scale_color_manual(name  ="aligners",
                     values = cols,
                     labels=c("last","minimap2","ngmlr"))+theme(legend.position=c(0.2,0.77))+
  theme(legend.background = element_rect(fill="white",size=.2, color='black',linetype="solid"))+
  # the actual simulated fusion breakpoints
  geom_rug(aes(x=fusion_truth$dist),alpha=0.15)


