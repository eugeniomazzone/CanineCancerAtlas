library(xlsx)
bp <- c('PRJNA752630')
tt <- c('B-Cell Lymphoma')
path <- 'CNA/'
subFiles <- list.files(path)

df <- data.frame()
for (samp in subFiles) {
	sample <- paste0(path,samp)
	sampleStat <- list.files(sample, pattern='txt')
	openStat <- read.table(paste0(sample,'/',sampleStat), header=T)
	toAna <- openStat[which.max(openStat$GoF),]
	toAna$sample <- paste0(sample,'/gamma',toAna$gamma,'/',list.files(paste0(sample,'/gamma',toAna$gamma), pattern='.cn'))
	toAna$BioProject <- bp[1]
	toAna$Histotype <- tt[1]
	df <- rbind(df,toAna)
}


df <- df[grepl(df$sample, pattern='cn'),]

wb <- createWorkbook()  
message("Creating sheet", 'Purity_Details')
sheet <- createSheet(wb, sheetName = 'Purity Details')
message("Adding data frame", 'Purity Details')
addDataFrame(df, sheet, row.names = F)

mydiff <- data.frame()
for(tt in unique(df$Histotype)){
newdf=subset(df, Histotype==tt)
results <- t(data.frame(quantile(newdf$aberrant.cell.fraction)))
rownames(results) <- tt
mydiff <- rbind(mydiff, results)
}

write.table(file='GoFresults.txt', row.names=FALSE,sep="\t", quote = FALSE, mydiff)

message("Creating sheet", 'Purity_Details')
sheet <- createSheet(wb, sheetName = 'Purity by Histo')
message("Adding data frame", 'Purity Details')
addDataFrame(mydiff, sheet, row.names = T)

saveWorkbook(wb, "TumorPurity.xlsx") 
