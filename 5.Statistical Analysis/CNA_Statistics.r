options(java.parameters = "-Xmx20000m")

library(sigminer)
library(dplyr)

biop <- c('/Annotation/TCN/')
tt <- c('B-cell Lymphoma')

df <- data.frame()
for (n in 1:length(biop)){
	subFiles <- list.files(biop[n])
	for (samp in subFiles) {
		sample <- paste0(biop[n],samp)
		openStat <- read.table(sample, header=T)
		colnames(openStat)[1]='Sample'
		openStat$Istotype <- tt[n]
		df <- rbind(df,openStat)
	}
}
colnames(df)[9]='minor_cn'
colnames(df)[1]='sample'

results <- data.frame()

dfToWork <- df
isto <- 'ALL'
dfToWork <-dfToWork[dfToWork$TCN!=2,]

delSize <- median(dfToWork$Width[dfToWork$TCN<2])
gainSize <- median(dfToWork$Width[dfToWork$TCN>2])
wholeSize <- median(dfToWork$Width[dfToWork$TCN!=2])

delN <- nrow(dfToWork[dfToWork$TCN<2,])
gainN <- nrow(dfToWork[dfToWork$TCN>2,])
wholeN <- nrow(dfToWork[dfToWork$TCN!=2,])

countDf <- dfToWork %>% count(Chr)
cyto <- read.table('/refFiles/cytoBandIdeo.txt', header=T)[1:38,c('chrom','chromEnd')]
cyto$chrom <- as.numeric(sapply(cyto$chrom, function(x) unlist(strsplit(x, 'chr'))[2]))
cyto <- cyto[order(as.numeric(cyto$chrom)),]
colnames(cyto)[2]<- 'ChrLen'

countDf <- merge(countDf, cyto, by.y='chrom', by.x='Chr')

statTest <- cor.test(countDf$n,countDf$ChrLen)
pval_all <- statTest$p.value
estim_all <- statTest$estimate

countDf <- dfToWork[dfToWork$TCN<2,] %>% count(Chr)

for (i in 1:38){
	if (!(i %in% countDf$Chr)){ countDf <- rbind(countDf, c(i,0)) }
}
countDf <- countDf[order(countDf$Chr),]
cyto <- read.table('/refFiles/cytoBandIdeo.txt', header=T)[1:38,c('chrom','chromEnd')]
cyto$chrom <- as.numeric(sapply(cyto$chrom, function(x) unlist(strsplit(x, 'chr'))[2]))
cyto <- cyto[order(as.numeric(cyto$chrom)),]
colnames(cyto)[2]<- 'ChrLen'

countDf <- merge(countDf, cyto, by.y='chrom', by.x='Chr')
statTest <- cor.test(countDf$n,countDf$ChrLen)
pval_del <- statTest$p.value
estim_del <- statTest$estimate

countDf <- dfToWork[dfToWork$TCN>2,] %>% count(Chr)
for (i in 1:38){
	if (!(i %in% countDf$Chr)){ countDf <- rbind(countDf, c(i,0)) }
}
countDf <- countDf[order(countDf$Chr),]
cyto <- read.table('/refFiles/cytoBandIdeo.txt', header=T)[1:38,c('chrom','chromEnd')]
cyto$chrom <- as.numeric(sapply(cyto$chrom, function(x) unlist(strsplit(x, 'chr'))[2]))
cyto <- cyto[order(as.numeric(cyto$chrom)),]
colnames(cyto)[2]<- 'ChrLen'
countDf <- merge(countDf, cyto, by.y='chrom', by.x='Chr')
statTest <- cor.test(countDf$n,countDf$ChrLen)
pval_gain <- statTest$p.value
estim_gain <- statTest$estimate

df2 <- data.frame(wholeSize, delSize, gainSize, wholeN, delN, gainN, pval_all, corr_all=estim_all, pval_gain, corr_gain=estim_gain, pval_del, corr_del=estim_del)
rownames(df2)=isto
results <- rbind(results,df2)

for (isto in unique(df$Istotype)){
	dfToWork <- subset(df, Istotype==isto)
	dfToWork <-dfToWork[dfToWork$TCN!=2,]

	delSize <- median(dfToWork$Width[dfToWork$TCN<2])
	gainSize <- median(dfToWork$Width[dfToWork$TCN>2])
	wholeSize <- median(dfToWork$Width[dfToWork$TCN!=2])

	delN <- nrow(dfToWork[dfToWork$TCN<2,])
	gainN <- nrow(dfToWork[dfToWork$TCN>2,])
	wholeN <- nrow(dfToWork[dfToWork$TCN!=2,])

	countDf <- dfToWork %>% count(Chr)
cyto <- read.table('/refFiles/cytoBandIdeo.txt', header=T)[1:38,c('chrom','chromEnd')]
cyto$chrom <- as.numeric(sapply(cyto$chrom, function(x) unlist(strsplit(x, 'chr'))[2]))
cyto <- cyto[order(as.numeric(cyto$chrom)),]
colnames(cyto)[2]<- 'ChrLen'
countDf <- merge(countDf, cyto, by.y='chrom', by.x='Chr')
	statTest <- cor.test(countDf$n,countDf$ChrLen)
	pval_all <- statTest$p.value
	estim_all <- statTest$estimate

	countDf <- dfToWork[dfToWork$TCN<2,] %>% count(Chr)

	for (i in 1:38){
		if (!(i %in% countDf$Chr)){ countDf <- rbind(countDf, c(i,0)) }
	}
	countDf <- countDf[order(countDf$Chr),]
cyto <- read.table('/refFiles/cytoBandIdeo.txt', header=T)[1:38,c('chrom','chromEnd')]
cyto$chrom <- as.numeric(sapply(cyto$chrom, function(x) unlist(strsplit(x, 'chr'))[2]))
cyto <- cyto[order(as.numeric(cyto$chrom)),]
colnames(cyto)[2]<- 'ChrLen'
countDf <- merge(countDf, cyto, by.y='chrom', by.x='Chr')
	statTest <- cor.test(countDf$n,countDf$ChrLen)
	pval_del <- statTest$p.value
	estim_del <- statTest$estimate

	countDf <- dfToWork[dfToWork$TCN>2,] %>% count(Chr)
	for (i in 1:38){
		if (!(i %in% countDf$Chr)){ countDf <- rbind(countDf, c(i,0)) }
	}
	countDf <- countDf[order(countDf$Chr),]
cyto <- read.table('/refFiles/cytoBandIdeo.txt', header=T)[1:38,c('chrom','chromEnd')]
cyto$chrom <- as.numeric(sapply(cyto$chrom, function(x) unlist(strsplit(x, 'chr'))[2]))
cyto <- cyto[order(as.numeric(cyto$chrom)),]
colnames(cyto)[2]<- 'ChrLen'
countDf <- merge(countDf, cyto, by.y='chrom', by.x='Chr')
	statTest <- cor.test(countDf$n,countDf$ChrLen)
	pval_gain <- statTest$p.value
	estim_gain <- statTest$estimate

	df2 <- data.frame(wholeSize, delSize, gainSize, wholeN, delN, gainN, pval_all, corr_all=estim_all, pval_gain, corr_gain=estim_gain, pval_del, corr_del=estim_del)
	rownames(df2)=isto
	results <- rbind(results,df2)
}

write.table(file='CorrelationTCN.tsv', results, sep='\t',quote=F)



