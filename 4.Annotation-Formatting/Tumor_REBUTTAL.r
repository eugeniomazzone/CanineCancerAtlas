library(xlsx)

cov <- read.table('Coverage_Tumor.tsv', sep=',')
colnames(cov) <- c('Sample','Coverage','BioProject')

pur <- read.xlsx('TumorPurity.xlsx',1)
pur <- pur[c(9, 2, 6, 10)]
colnames(pur) <- c('Sample','Ploidy','Purity','BioProject')

pur$Sample <- do.call(rbind,strsplit(pur$Sample, '[/]'))[,4]
pur$Sample <- do.call(rbind,strsplit(pur$Sample, '[_]'))[,2]

df <- merge(pur, cov, by='Sample', all.x=T, all.y=T)
df <- df[c(1:3,5:6)]
colnames(df)[5] = 'BioProject'
df$EstimatedBy <- 'ASCAT'
df$EstimatedBy[is.na(df$Ploidy)] = 'Median'
df$Ploidy[is.na(df$Ploidy)] =2
df$Purity[is.na(df$Purity)] =median(df$Purity[!is.na(df$Purity)])
df$NRPCC <- (df$Coverage * df$Purity) / df$Ploidy
df$minCov <- (15 * df$Ploidy) / df$Purity

wb <- createWorkbook()  

sheet <- createSheet(wb, sheetName = 'NRPCC')
addDataFrame(df, sheet, row.names = T)

saveWorkbook(wb, "NRPCC_REBUTTAL.xlsx") 

write.table(file='NRPCC_REBUTTAL.tsv', df, row.names=F, quote=F, sep='\t')

