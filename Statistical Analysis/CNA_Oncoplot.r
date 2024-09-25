library(maftools)
library(dplyr)
library(writexl)
library(xlsx)
library(ComplexHeatmap)

cn <- read.table('CNA_Annotated.csv', header=T)
metadataSample <- unique(cn %>% subset(select=c(sample, isto)))
cn$isto <- NULL

pv <- tidyr::pivot_wider(cn, names_from='sample', values_from='Status')

M <- as.data.frame(pv)
rownames(M) <- M$genes
M$genes <- NULL
M[M=='Del'] ='Del;'
M[M=='Amp'] ='Amp;'
for (i in 1:nrow(M)){
	M[i,sapply(M[i,],function(x) identical(x,list(c('Amp', 'Del'))))]='MultiHit'
}

countAmp <- rowSums(M == "Amp;")
countDel <- rowSums(M == "Del;")
countAmp <- sort(countAmp, decreasing=T)
countDel <- sort(countDel, decreasing=T)

top=30
AmpToShow <- names(countAmp[1:top]) 
countDel <- countDel[!(names(countDel)=='IGLV2-33')]
DelToShow <- names(countDel[1:top]) 

dfAmp <- M[AmpToShow,]
dfDel <- M[DelToShow,]
dfAmp[!(dfAmp == 'Amp;')]=''
dfDel[!(dfDel == 'Del;')]=''

col = c("Del" = "blue", "Amp" = "red", 'MultiHit'='purple')

fabcolors = c('#F19A94','#3E6DCA','#27DF7A','#A646A4','#E88ABD','#B0D022','#AF541F','#FFDC40','#80C5ED','#FF3953')
col.assign <- setNames(fabcolors, unique(metadataSample$isto))
rownames(metadataSample) <- metadataSample$sample
match.id <- match(metadataSample$sample, colnames(dfAmp))
df <- data.frame(metadataSample[match.id, 'isto'])
rownames(df) <- metadataSample[match.id, 'sample']
mix.ha <- HeatmapAnnotation(df = df, col = list(isto = col.assign))

sample.order <- NULL
for (tp in unique(metadataSample$isto))
{
	isto_sp <- subset(metadataSample, isto==tp)$sample
	dfW <- dfDel[isto_sp]
	sampleC <- colSums(dfW=='Del;')[order(colSums(dfW=='Del;'),decreasing=T)]
	sampleC <- sampleC[sampleC!=0]
	sample.order <- c(sample.order, names(sampleC))
}



alter_fun = function(x, y, w, h, v) {
  # background
  grid.rect(x, y, w-unit(2, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  # graphics for alterations
  if(v["Del"]) grid.rect(x, y, w-unit(2, "mm"), h-unit(1, "mm"), gp = gpar(fill = col["Del"], col = NA))
}


ht1 = oncoPrint(dfDel,
          alter_fun = alter_fun, col = col, 
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          pct_side = "right", row_names_side = "left",
          use_raster=TRUE, bottom_annotation = mix.ha, , column_order=sample.order)
  
jpeg(filename='Oncoplot_DEL.jpeg', width=6000, height=3000, res=300)
draw(ht1)
dev.off()

alter_fun = function(x, y, w, h, v) {
  # background
  grid.rect(x, y, w-unit(2, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  # graphics for alterations
  if(v["Amp"])    grid.rect(x, y, w-unit(2, "mm"), h-unit(1, "mm"), gp = gpar(fill = col["Amp"], col = NA))
}

sample.order <- NULL
for (tp in unique(metadataSample$isto))
{
	isto_sp <- subset(metadataSample, isto==tp)$sample
	dfW <- dfAmp[isto_sp]
	sampleC <- colSums(dfW=='Amp;')[order(colSums(dfW=='Amp;'),decreasing=T)]
	sampleC <- sampleC[sampleC!=0]
	sample.order <- c(sample.order, names(sampleC))
}

ht1 = ComplexHeatmap::oncoPrint(dfAmp,
          alter_fun = alter_fun, col = col, 
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          pct_side = "right", row_names_side = "left",
          use_raster=TRUE, bottom_annotation = mix.ha, column_order=sample.order)
  
jpeg(filename='Oncoplot_AMP.jpeg', width=6000, height=3000, res=300)
draw(ht1)
dev.off()


