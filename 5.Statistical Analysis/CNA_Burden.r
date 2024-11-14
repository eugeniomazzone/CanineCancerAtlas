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

isto <- 'ALL'
df <- df[df$TCN!=2,] ### Loss-Gain will not be considered
alteredSample <- df %>% group_by(sample, Istotype) %>% summarize( ALT_L = sum(Width))
genomeLeng <- sum(read.table('canine_cytoband_for_gistic_order.txt')[1:38,'V5'])
alteredSample$cnaBurden <- alteredSample$ALT_L/genomeLeng

df2 <- quantile(alteredSample$cnaBurden)
df2 <- t(data.frame(df2))
rownames(df2)=isto
results <- rbind(results,df2)

for (isto in unique(alteredSample$Istotype)){
	dfToWork <- subset(alteredSample, Istotype==isto)
	df2 <- quantile(dfToWork$cnaBurden)
	df2 <- t(data.frame(df2))
	rownames(df2)=isto
	results <- rbind(results,df2)
}


wb <- createWorkbook()
sheet <- createSheet(wb, sheetName = 'CNA Burden')
addDataFrame(results, sheet, row.names = T)
saveWorkbook(wb, "CNA_Burden.xlsx") 



