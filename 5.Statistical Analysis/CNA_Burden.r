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

isto <- 'ALL'
df <- df[df$TCN!=2,] ### Loss-Gain will not be considered
alteredSample <- df %>% group_by(sample, Istotype) %>% summarize( ALT_L = sum(Width))
cyto <- read.table('/refFiles/cytoBandIdeo.txt', header=T)[1:38,c('chrom','chromEnd')]
cyto$chrom <- as.numeric(sapply(cyto$chrom, function(x) unlist(strsplit(x, 'chr'))[2]))
cyto <- cyto[order(as.numeric(cyto$chrom)),]
colnames(cyto)[2]<- 'ChrLen'
genomeLeng <- sum(cyto$ChrLen)
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

write.table(file='CNA_Burden.tsv', results, sep='\t',quote=F)



