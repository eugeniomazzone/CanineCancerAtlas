options(java.parameters = "-Xmx20000m")

library(xlsx)
library(sigminer)
library(dplyr)

biop <- c('PRJNA752630')
tt <- c('B-cell Lymphoma')

df <- data.frame()
for (n in 1:length(biop)){
	subFiles <- list.files(paste0(biop[n],'/TCN/'))
	for (samp in subFiles) {
		sample <- paste0(biop[n],'/TCN/',samp)
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
cnaBurden <- 

countDf <- dfToWork %>% count(Chr)
countDf$ChrLen <- read.table('canine_cytoband_for_gistic_order.txt')[1:38,'V5']
statTest <- cor.test(countDf$n,countDf$ChrLen)
pval_all <- statTest$p.value
estim_all <- statTest$estimate

countDf <- dfToWork[dfToWork$TCN<2,] %>% count(Chr)

for (i in 1:38){
	if (!(i %in% countDf$Chr)){ countDf <- rbind(countDf, c(i,0)) }
}
countDf <- countDf[order(countDf$Chr),]
countDf$ChrLen <- read.table('canine_cytoband_for_gistic_order.txt')[1:38,'V5']
statTest <- cor.test(countDf$n,countDf$ChrLen)
pval_del <- statTest$p.value
estim_del <- statTest$estimate

countDf <- dfToWork[dfToWork$TCN>2,] %>% count(Chr)
for (i in 1:38){
	if (!(i %in% countDf$Chr)){ countDf <- rbind(countDf, c(i,0)) }
}
countDf <- countDf[order(countDf$Chr),]
countDf$ChrLen <- read.table('canine_cytoband_for_gistic_order.txt')[1:38,'V5']
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
	countDf$ChrLen <- read.table('canine_cytoband_for_gistic_order.txt')[1:38,'V5']
	statTest <- cor.test(countDf$n,countDf$ChrLen)
	pval_all <- statTest$p.value
	estim_all <- statTest$estimate

	countDf <- dfToWork[dfToWork$TCN<2,] %>% count(Chr)

	for (i in 1:38){
		if (!(i %in% countDf$Chr)){ countDf <- rbind(countDf, c(i,0)) }
	}
	countDf <- countDf[order(countDf$Chr),]
	countDf$ChrLen <- read.table('canine_cytoband_for_gistic_order.txt')[1:38,'V5']
	statTest <- cor.test(countDf$n,countDf$ChrLen)
	pval_del <- statTest$p.value
	estim_del <- statTest$estimate

	countDf <- dfToWork[dfToWork$TCN>2,] %>% count(Chr)
	for (i in 1:38){
		if (!(i %in% countDf$Chr)){ countDf <- rbind(countDf, c(i,0)) }
	}
	countDf <- countDf[order(countDf$Chr),]
	countDf$ChrLen <- read.table('canine_cytoband_for_gistic_order.txt')[1:38,'V5']
	statTest <- cor.test(countDf$n,countDf$ChrLen)
	pval_gain <- statTest$p.value
	estim_gain <- statTest$estimate

	df2 <- data.frame(wholeSize, delSize, gainSize, wholeN, delN, gainN, pval_all, corr_all=estim_all, pval_gain, corr_gain=estim_gain, pval_del, corr_del=estim_del)
	rownames(df2)=isto
	results <- rbind(results,df2)
}

wb <- createWorkbook()
sheet <- createSheet(wb, sheetName = paste0('cn_','List'))
addDataFrame(results, sheet, row.names = T)
saveWorkbook(wb, "CNA_Correlation.xlsx") 



