options(java.parameters = "-Xmx25000m")

library(maftools)
library(viridis)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(ComplexHeatmap)

biop <- c('/Annotation/Delly_sansa_tsv/')
tt <- c('B-cell Lymphoma')

bt.genes <- rtracklayer::import('/refFiles/Canis_lupus_familiaris.CanFam3.1.104.gtf.gz')
bt.genes <- subset(bt.genes, select=c(gene_id,gene_name))
bt.genes <- data.frame(bt.genes)
bt.genes <- unique(bt.genes[c('gene_id','gene_name','start','end','seqnames')])
bt.genes <-bt.genes[!is.na(bt.genes$gene_name),]
bt.genes <-bt.genes[(bt.genes$seqnames %in% 1:38),]
rownames(bt.genes)<-NULL
bt.genes <- makeGRangesFromDataFrame(bt.genes, keep.extra.columns=TRUE)
bt.genes <- unique(data.frame(bt.genes)[,c('gene_id', 'gene_name')])
rownames(bt.genes) <- bt.genes$gene_id
bt.genes$gene_id <- NULL

#### READING AND LOADING ALL MAF FILES

sv <- data.frame()
for (n in 1:length(biop)) {
	print(biop[n])
	path2<-biop[n]
	files2 <- list.files(path2)
	for (f in files2) { 
		sv2 <- read.table(paste0(path2,f), sep='\t',header=T, row.names=NULL)
		if (nrow(sv2)>0){
			sv2$Histotype <- tt[n]
			sv2$sample <- f
			sv <- rbind(sv,sv2)
		}
	}
}


sv <- sv[!is.na(sv$query.startfeature) | !is.na(sv$query.endfeature) | !is.na(sv$query.containedfeature),]

mySvDf <- data.frame()
for(tp in unique(tt)){ 
	sv2 <- subset(sv, Histotype==tp)
	for (samp in unique(sv2$sample)){
		print(samp)
		sv3 <- subset(sv2 , sample == samp)
		
		for (qtp in unique(sv3$query.svtype)){
		sv4 <- subset(sv3 , query.svtype == qtp)
		sv_start=unique(sort(unlist(strsplit(sv4$query.startfeature, '\\(|\\)|\\.'))))
		sv_start <- sv_start[nchar(sv_start)>7]
		sv_end=unique(sort(unlist(strsplit(sv4$query.endfeature, '\\(|\\)|\\.'))))
		sv_end <- sv_end[nchar(sv_end)>7]
		sv_ctn <- unique(sort(unlist(strsplit(sv4$query.containedfeature, '\\(|\\)|\\,|\\.'))))
		sv_ctn <- sv_ctn[nchar(sv_ctn)>7]
		svSample <- data.frame(genes = unique(sort(c(sv_start, sv_end, sv_ctn))))
		svSample$SVtype <- qtp
		svSample$Histotype <- tp
		svSample$sample <- unlist(strsplit(samp,'\\.'))[1]
		mySvDf <- rbind(mySvDf,svSample)
		}
	}
}


mySvDf$genes <- bt.genes[mySvDf$genes,]
mySvDf <- mySvDf[!is.na(mySvDf$genes),]

sv2 <- mySvDf
sv2$Histotype <- NULL ### Remove Histotypes
pv <- tidyr::pivot_wider(sv2, names_from='sample', values_from='SVtype')
M <- as.data.frame(pv)
rownames(M) <- M$genes
M$genes <- NULL 
for (i in 1:nrow(M)){
	M[i,sapply(M[i,],function(x) length(unlist(x))>1)]='MultiHit'
}

countAmp <- rowSums(M=='DUP' | M=='MultiHit' | M=='BND' | M=='INS' | M=='INV')
countAmp <- sort(countAmp, decreasing=T)

top=30
AmpToShow <- names(countAmp[1:top]) 

dfSV <- M[AmpToShow,]
dfSV[(dfSV!='DUP' & dfSV!='MultiHit' & dfSV!='BND' & dfSV!='INS' & dfSV!='INV')]=''
dfSV[is.na(dfSV)]=''

sample.order <- NULL
for (tp in unique(metadataSample$Histotype))
{
	isto_sp <- subset(df, Histotype==tp)$sample
	dfW <- dfSV[isto_sp]
	sampleC <- colSums(dfW!='')[order(colSums(dfW!=''),decreasing=T)]
	sampleC <- sampleC[sampleC!=0]
	sample.order <- c(sample.order, names(sampleC))
}


col = c("INV" = "green", "DUP" = "red", "BND" = "blue", 'MultiHit'='purple')
alter_fun = function(x, y, w, h, v) {
  # background
  grid.rect(x, y, w-unit(1, "mm"), h-unit(2, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  # graphics for alterations
  if(v["INV"])    grid.rect(x, y, w-unit(2, "mm"), h-unit(1, "mm"), gp = gpar(fill = col["INV"], col = NA))
  if(v["BND"])    grid.rect(x, y, w-unit(2, "mm"), h-unit(1, "mm"), gp = gpar(fill = col["BND"], col = NA))
  if(v["DUP"])    grid.rect(x, y, w-unit(2, "mm"), h-unit(1, "mm"), gp = gpar(fill = col["DUP"], col = NA))
  if(v["MultiHit"])    grid.rect(x, y, w-unit(2, "mm"), h-unit(1, "mm"), gp = gpar(fill = col["MultiHit"], col = NA))
}


ht1 = oncoPrint(dfSV,
          alter_fun = alter_fun, col = col, 
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          pct_side = "right", row_names_side = "left",
          use_raster=TRUE, column_order=sample.order)
  
tiff(filename='oncoplotSV.tiff', width=6000, height=3000, res=300)
draw(ht1)
dev.off()






