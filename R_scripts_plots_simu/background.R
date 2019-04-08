# reads with noise (background reads randomly sampled from transcriptome, no fusion)

# read TPFP results of noise files
ngmlrnoise_pb <- read.table('/Users/yaoyao/Desktop/PB_bamfiles/background/tpfp/ngmlr_pb_noise_3.txt')
ngmlrnoise_ont <- read.table('/Users/yaoyao/Desktop/PB_bamfiles/background/tpfp/ngmlr_ont_noise_3.txt')
mnp2noise_pb <- read.table('/Users/yaoyao/Desktop/PB_bamfiles/background/tpfp/mnp2_pb_noise_3.txt')
mnp2noise_ont <- read.table('/Users/yaoyao/Desktop/PB_bamfiles/background/tpfp/mnp2_ont_noise_3.txt')
lastnoise_pb <- read.table('/Users/yaoyao/Desktop/PB_bamfiles/background/tpfp/last_pb_noise.txt')
lastnoise_ont <- read.table('/Users/yaoyao/Desktop/PB_bamfiles/background/tpfp/last_ont_noise.txt')


# remove reads containing NULLs (genes without annotation)
ngmlrnoise_pb <- ngmlrnoise_pb %>% filter(V2 != "NULL" & V3 != "NULL"& !grepl('NULL',V10))
ngmlrnoise_ont <- ngmlrnoise_ont %>% filter(V2 != "NULL" & V3 != "NULL"& !grepl('NULL',V10))
mnp2noise_pb <- mnp2noise_pb %>% filter(V2 != "NULL" & V3 != "NULL" & !grepl('NULL',V10))
mnp2noise_ont <- mnp2noise_ont %>% filter(V2 != "NULL" & V3 != "NULL" & !grepl('NULL',V10))
lastnoise_pb <- lastnoise_pb %>% filter(V2 != "NULL" & V3 != "NULL" & !grepl('NULL',V10))
lastnoise_ont <- lastnoise_ont %>% filter(V2 != "NULL" & V3 != "NULL" & !grepl('NULL',V10))


# number of fusion genes output 
length(unique(ngmlrnoise_pb$V10)) #361
length(unique(ngmlrnoise_ont$V10)) #319
length(unique(mnp2noise_pb$V10)) #485
length(unique(mnp2noise_ont$V10)) #580
length(unique(lastnoise_pb$V10)) #2171
length(unique(lastnoise_ont$V10)) #1438

# number of TRUE fusion genes output
length(unique((ngmlrnoise_pb[ngmlrnoise_pb$V7=='TP',])$V10)) #104
length(unique((ngmlrnoise_ont[ngmlrnoise_ont$V7=='TP',])$V10)) #99
length(unique((mnp2noise_pb[mnp2noise_pb$V7=='TP',])$V10)) #318
length(unique((mnp2noise_ont[mnp2noise_ont$V7=='TP',])$V10)) #375
length(unique((lastnoise_pb[lastnoise_pb$V7=='TP',])$V10)) #353
length(unique((lastnoise_ont[lastnoise_ont$V7=='TP',])$V10)) #358

# number of FALSE fusion genes output
length(unique((ngmlrnoise_pb[ngmlrnoise_pb$V7=='FP',])$V10)) #257
length(unique((ngmlrnoise_ont[ngmlrnoise_ont$V7=='FP',])$V10)) #220
length(unique((mnp2noise_pb[mnp2noise_pb$V7=='FP',])$V10)) #167
length(unique((mnp2noise_ont[mnp2noise_ont$V7=='FP',])$V10)) #205
length(unique((lastnoise_pb[lastnoise_pb$V7=='FP',])$V10)) #1818
length(unique((lastnoise_ont[lastnoise_ont$V7=='FP',])$V10)) #1080



a<-vector()
b<-vector()
c<-vector()
# FP reads/ all reads [after adding background noise]
a[2]<-sum(ngmlrnoise_pb$V7=='FP')/nrow(ngmlrnoise_pb) # 0.54444
a[4]<-sum(ngmlrnoise_ont$V7=='FP')/nrow(ngmlrnoise_ont) #0.5053717
b[2]<-sum(mnp2noise_pb$V7=='FP')/nrow(mnp2noise_pb) # 0.2494612
b[4]<-sum(mnp2noise_ont$V7=='FP')/nrow(mnp2noise_ont)  # 0.2539998
c[2]<-sum(lastnoise_pb$V7=='FP')/nrow(lastnoise_pb)  # 0.552184
c[4]<-sum(lastnoise_ont$V7=='FP')/nrow(lastnoise_ont) # 0.476945

# FP reads/ all reads [before, without background noise]
a[1]<- sum(ngmlrtpfp$V7=='FP')/nrow(ngmlrtpfp) #0.4015595
a[3]<- sum(ngmlront$V7=='FP')/nrow(ngmlront) #0.372752
b[1]<-sum(mnp2tpfp$V7=='FP')/nrow(mnp2tpfp) #0.2443721
b[3]<-sum(mnp2ont$V7=='FP')/nrow(mnp2ont) #0.252987
c[1]<-sum(lasttpfp$V7=='FP')/nrow(lasttpfp) #0.411925
c[3]<-sum(lastont$V7=='FP')/nrow(lastont) # 0.3597183
dfFDRread <- data.frame(dataset=c('pacbio','pacbio_noise','ont','ont_noise'),ngmlr=a,minimap2=b,last=c)

e<-vector()
f<-vector()
g<-vector()
# FP genes /all fusion candidates [with background]
e[2]<- length(unique((ngmlrnoise_pb[ngmlrnoise_pb$V7=='FP',])$V10))/length(unique(ngmlrnoise_pb$V10)) #0.7119114
e[4]<- length(unique((ngmlrnoise_ont[ngmlrnoise_ont$V7=='FP',])$V10))/length(unique(ngmlrnoise_ont$V10)) #0.6896552
f[2]<- length(unique((mnp2noise_pb[mnp2noise_pb$V7=='FP',])$V10))/length(unique(mnp2noise_pb$V10)) #0.3443299
f[4]<- length(unique((mnp2noise_ont[mnp2noise_ont$V7=='FP',])$V10))/length(unique(mnp2noise_ont$V10))  #0.3534483
g[2]<- length(unique((lastnoise_pb[lastnoise_pb$V7=='FP',])$V10))/length(unique(lastnoise_pb$V10)) #0.8374021
g[4]<- length(unique((lastnoise_ont[lastnoise_ont$V7=='FP',])$V10))/length(unique(lastnoise_ont$V10))#0.7510431


# FP genes /all fusion candidates [before, without background]
e[1]<- length(unique((ngmlrtpfp[ngmlrtpfp$V7 =='FP',])$V10))/length(unique(ngmlrtpfp$V10)) #0.4825871
e[3]<- length(unique((ngmlront[ngmlront$V7 =='FP',])$V10))/length(unique(ngmlront$V10)) # 0.484375
f[1]<- length(unique((mnp2tpfp[mnp2tpfp$V7=='FP',])$V10))/length(unique(mnp2tpfp$V10)) #0.3276956
f[3]<- length(unique((mnp2ont[mnp2ont$V7 =='FP',])$V10))/length(unique(mnp2ont$V10)) #0.3444056
g[1]<- length(unique((lasttpfp[lasttpfp$V7 =='FP',])$V10))/length(unique(lasttpfp$V10)) #0.7031119
g[3]<- length(unique((lastont[lastont$V7 =='FP',])$V10))/length(unique(lastont$V10)) #0.6057269

dfFDRgene<- data.frame(dataset=c('pacbio','pacbio_noise','ont','ont_noise'),ngmlr=e,minimap2=f,last=g)

