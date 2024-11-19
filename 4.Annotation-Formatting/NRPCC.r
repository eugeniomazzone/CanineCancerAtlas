cov <- read.table('Depth.txt', sep='\t', header=T)

pur <- read.table('GoFresults.txt', sep='\t', header=T)
pur <- pur[c(2,6,9, 11)]
colnames(pur) <- c('Ploidy','Purity','Path','Sample')

df <- merge(pur, cov, by='Sample', all.x=T, all.y=T)
df$EstimatedBy <- 'ASCAT'
df$EstimatedBy[is.na(df$Ploidy)] = 'Median'
df$Ploidy[is.na(df$Ploidy)] =2
df$Purity[is.na(df$Purity)] =median(df$Purity[!is.na(df$Purity)])
df$NRPCC <- (df$Tumor.Depth * df$Purity) / df$Ploidy
df$minCov <- (15 * df$Ploidy) / df$Purity

df$Excluded <- 


write.table(file='NRPCC.tsv', df, row.names=F, quote=F, sep='\t')

