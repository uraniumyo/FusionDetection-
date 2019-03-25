suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
})


# read runtime data
time <- read.csv('/Users/yaoyao/Desktop/runtime.csv',sep = '',header=FALSE)
rownames(time)<-time[,1] 
time<-time[,-1] 
t <- as.data.frame(t(time))

# plot
colorlist <- brewer.pal(7, "Set1")[c(1,2,3)]
cols <- c("last"=colorlist[2],"minimap2"=colorlist[1],"ngmlr"=colorlist[3])

#CPU run time for aligners             
ggplot(t)+ theme_bw()+
  geom_bar(aes(x = factor(coverage), y = size),stat= 'identity',width=0.65,alpha=0.6)+
  geom_point(aes(x=factor(coverage),y = last* 1 / 200,color="last"))+
  geom_line(aes(x=factor(coverage),y = last* 1 / 200,group = 1,color="last"),size=0.9)+
  geom_point(aes(x=factor(coverage),y = minimap2* 1 / 200,color="minimap2"))+ 
  geom_line(aes(x=factor(coverage),y = minimap2* 1 / 200,group = 1,color="minimap2"),color=colorlist[1],size=0.9)+
  geom_point(aes(x=factor(coverage),y = ngmlr* 1 / 200,color="ngmlr"))+ 
  geom_line(aes(x=factor(coverage),y = ngmlr* 1 / 200,group = 1,color="ngmlr"),size=0.9)+
  scale_y_continuous(name = 'size(Mb)', sec.axis = sec_axis(~ . * 200, name = 'runtime(sec)'), limits = c(0, 110))+
  scale_color_manual(name  ="aligners",
                     values = cols,
                       labels=c("last","minimap2","ngmlr"))+
  theme(legend.position=c(0.2,0.83))+ 
  theme(legend.background = element_rect(fill="white",size=.2, color='black',linetype="solid"))+
  xlab('coverages')+ ylab('runtime(sec)')

