suppressPackageStartupMessages({
library(dplyr)
library(rtracklayer)
  })

# read fusion summary table
sumtable <- read.csv("summary.csv", header = TRUE, stringsAsFactors=FALSE)

# read gtf file ,select gene id/symbol
gtf1 <- rtracklayer::import('gencode.v29.annotation.gtf')
gtf_df <- as.data.frame(gtf1)
geneid_df <- dplyr::select(gtf_df,c(gene_id,gene_name))
colnames(geneid_df) <- c('geneid','genename')
# merge replicates
geneid_df <- geneid_df[!duplicated(geneid_df$geneid),]

# add corresponding gene symbol
sumtable <- left_join(sumtable, geneid_df,by = c('gene_A'='geneid'))
colnames(sumtable)[ncol(sumtable)] <- c('genename_A')

sumtable <- left_join(sumtable, geneid_df,by = c('gene_B'='geneid'))
colnames(sumtable)[ncol(sumtable)] <- c('genename_B')

# write to csv
write.table(x = sumtable ,file = "fusion_summary.csv",row.names = TRUE,
          col.names = TRUE)

