tt <- c('B-Cell Lymphoma')
path <- '/WorkDir/CNA/'
subFiles <- list.files(path)
subFiles <- subFiles[grep(subFiles, pattern='rda', invert=T)]

df <- data.frame()
for (samp in subFiles) {
	sample <- paste0(path,samp,'/ASCAT/ASCN/')
	sampleStat <- list.files(sample, pattern='txt')
	openStat <- read.table(paste0(sample,'/',sampleStat), header=T)
	toAna <- openStat[which.max(openStat$GoF),]
	toAna$sample <- paste0(sample,'gamma',format(round(toAna$gamma,2),nsmall=2),'/',list.files(paste0(sample,'/gamma',format(round(toAna$gamma,2),nsmall=2)), pattern='.cn'))
	toAna$Histotype <- tt[1]
	toAna$Sample <- samp
	df <- rbind(df,toAna)
}


df <- df[grepl(df$sample, pattern='cn'),]

write.table(file='GoFresults.txt', row.names=FALSE,sep="\t", quote = FALSE, df)

