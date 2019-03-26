suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(Biostrings)
  library(reshape2)
  library(ggthemes)
  library(grid)
  library(scales) 
  library(RColorBrewer)
})

########################### set color for aligners #################

colorlist <- brewer.pal(7, "Set1")[c(1,2,3)]
cols <- c("last"=colorlist[2],"minimap2"=colorlist[1],"ngmlr"=colorlist[3])


######################## read fusion files ##########################

# read filtered fusion reads
read_file<- function(u){
  no_col <- max(count.fields(u, sep = ""))
  reads <- read.table(u,sep="",fill=TRUE,header = F,col.names=c(1:no_col),na.strings = '',stringsAsFactors = FALSE)
  # swap fusiontype to last column 
  for (i in 1:dim(reads)[1]){
    count <- plyr::count(as.character(is.na(reads[i,])))$freq[1]
    if (is.na(reads[i,no_col])){
      reads[i,no_col] <- reads[i,count]
      reads[i,count] <- NA
    }
  }
  names(reads)<- c('qread',paste(c('chr','pos_start','pos_end','strand','gene','name','feature'),rep(1:((no_col-2)/7),each=7),sep=''),'type')
  return(reads)
}

# read bam files (skbr3)
ngmlrskbr3_1 <- read_file('/Users/yaoyao/Desktop/PB_bamfiles/ngmlr_skbr3.txt')
mnp2skbr3_1 <- read_file('/Users/yaoyao/Desktop/PB_bamfiles/mnp_skbr3.txt')
lastskbr3_1 <- read_file('/Users/yaoyao/Desktop/PB_bamfiles/last_skbr3.txt')


######################## extract fusion names ##########################


get_fusion_names <- function(u){
  namecol <- grep('name',colnames(u),value = TRUE)
  u2 <- u[c('qread',namecol)]
  u3 <- u2[c(namecol)]
  # replace NULL with NA
  # u3[apply(u3, 2, function(x) x=="NULL")] = NA
  mat <- as.matrix(u3)
  
  # get gene names 
  list <- list()
  for (i in 1:dim(mat)[1]){
    list[[i]] <- as.vector(unique(na.omit(mat[i,])))
  }
  # keep reads with more than 2 alignments
  a <- list[lengths(list) > 1]
  
  # combine gene names as fusions
  list2 <- list()
  for (i in 1:length(a)){
    list2[[i]] <- paste(a[[i]],collapse = "_")
  }
  count_freq <- plyr::count(unlist(list2))
  # remove reads with alignments to NULL
  # only reads with all alignments annotated are kept 
  count_freq <- count_freq %>% filter(!grepl('NULL',x))
  return(count_freq)
} 

ngmlr_skbr3_freq <- get_fusion_names(ngmlrskbr3_1)
mnp2_skbr3_freq <- get_fusion_names(mnp2skbr3_1)
last_skbr3_freq <- get_fusion_names(lastskbr3_1)

# number of suspicious fusions 
nrow(ngmlr_skbr3_freq) #10408
nrow(mnp2_skbr3_freq) #16157
nrow(last_skbr3_freq) #44061

# when number of supporting reads is 2
dim(ngmlr_skbr3_freq[ngmlr_skbr3_freq $freq>=2,]) #693
dim(mnp2_skbr3_freq[mnp2_skbr3_freq $freq>=2,]) #245
dim(last_skbr3_freq[last_skbr3_freq $freq>=2,]) #2180

# more than 2
dim(ngmlr_skbr3_freq[ngmlr_skbr3_freq $freq>2,]) #343
dim(mnp2_skbr3_freq[mnp2_skbr3_freq $freq>2,]) #60
dim(last_skbr3_freq[last_skbr3_freq $freq>2,]) #812


######################## filter fusions with thresholds ##########################


# filter results with number of supporting reads as threshold
filter_with_freq <- function(u,v,w){
  a <- vector()
  b <- vector()
  c <- vector()
  inatleast2 <- vector()
  inall3 <- vector()
  # get number of fusions with number i supporting reads
  for (i in 1:5){
   a[i]<-dim(u[u$freq>=i,])[1]
   b[i]<-dim(v[v$freq>=i,])[1]
   c[i]<-dim(w[w$freq>=i,])[1]
   #intersect 2 aligners
   nandm <- intersect(u[u$freq>=i,]$x,v[v$freq>=i,]$x)
   nandl <- intersect(u[u$freq>=i,]$x,w[w$freq>=i,]$x)
   mandl <- intersect(v[v$freq>=i,]$x,w[w$freq>=i,]$x) 
   #number of fusions detected by atleast 2 aligners
   inatleast2[i] <- length(union(mandl,union(nandm,nandl)))
   #number of fusions detected by all 3 aligners
   inall3[i] <- length(intersect(intersect(nandm,nandl),mandl))
  }
  return(list(ngmlr=a,minimap2=b,last=c,inatleast2=inatleast2,inall3=inall3))
}

nfusion_with_threshold <- filter_with_freq(ngmlr_skbr3_freq,mnp2_skbr3_freq,last_skbr3_freq)
nfusion_with_threshold <- as.data.frame(nfusion_with_threshold)
coverage <- data.frame(coverage=c(1:5))
nfusion_with_threshold <- cbind(coverage,nfusion_with_threshold)

# plot: number of fusions 
ggplot(nfusion_with_threshold)+
geom_point(aes(x=coverage,y=ngmlr,color="#4DAF4A" ))+geom_line(aes(x=coverage,y=ngmlr,color="#4DAF4A" ), size=0.5)+
geom_point(aes(x=coverage,y=minimap2,color="#E41A1C"))+geom_line(aes(x=coverage,y=minimap2,color="#E41A1C"), size=0.5)+
geom_point(aes(x=coverage,y=last,color="#377EB8" ))+geom_line(aes(x=coverage,y=last,color="#377EB8" ), size=0.5)+
geom_point(aes(x=coverage,y=inatleast2,color='orange'))+geom_line(aes(x=coverage,y=inatleast2,color='orange'), size=0.7)+
geom_point(aes(x=coverage,y=inall3,color='yellow'))+geom_line(aes(x=coverage,y=inall3,color='yellow'), size=0.7)+
  scale_color_manual(name  ="",
                     values = c("#4DAF4A","#E41A1C","#377EB8" ,'orange','yellow'),
                     labels=c("last","minimap2","ngmlr",'in2aligners','in3aligners'))+
  #ggtitle('number of fusion detected at each threshold')+
ylab('number of fusions') + xlab('number of supporting reads')


# check results when threshold=2 (keep fusions>=2 reads)
oneinall3 <-intersect(ngmlr_skbr3_freq[ngmlr_skbr3_freq $freq==1,]$x,
                      intersect(mnp2_skbr3_freq[mnp2_skbr3_freq $freq==1,]$x,
                                last_skbr3_freq[last_skbr3_freq $freq==1,]$x)) #2629

nandm <- intersect(ngmlr_skbr3_freq[ngmlr_skbr3_freq $freq>=2,]$x,
                   mnp2_skbr3_freq[mnp2_skbr3_freq $freq>=2,]$x) #86 
nandl <- intersect(ngmlr_skbr3_freq[ngmlr_skbr3_freq $freq>=2,]$x,
                   last_skbr3_freq[last_skbr3_freq $freq>=2,]$x) #290 
mandl <- intersect(mnp2_skbr3_freq[mnp2_skbr3_freq $freq>=2,]$x,
                   last_skbr3_freq[last_skbr3_freq $freq>=2,]$x) #84 

inatleast2 <- union( mandl,union(nandm,nandl)) #360 

inall3 <- intersect(intersect(nandm,nandl),mandl) #50 

output <- data.frame(fusion=inall3)

write.table(output, file = "/Users/yaoyao/Desktop/PB_bamfiles/skbr3_fusion_result.txt",row.names=FALSE,sep="\t", quote = FALSE)
