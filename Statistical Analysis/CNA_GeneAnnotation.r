options(java.parameters = "-Xmx7000m")

library(xlsx)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)


biop <- c('PRJNA752630')
tt <- c('B-cell Lymphoma')

#### READING DATA
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
rownames(df)<-NULL
colnames(df)[9]='minor_cn'
colnames(df)[1]='sample'

#### READING GTF
bt.genes <- rtracklayer::import('Canis_lupus_familiaris.CanFam3.1.104.gtf.gz')
bt.genes <- subset(bt.genes, select=c(gene_id,gene_name))
bt.genes <- data.frame(bt.genes)
bt.genes <- unique(bt.genes[c('gene_id','gene_name','start','end','seqnames')])
bt.genes <-bt.genes[!is.na(bt.genes$gene_name),]
bt.genes <-bt.genes[(bt.genes$seqnames %in% 1:38),]
rownames(bt.genes)<-NULL
bt.genes <- makeGRangesFromDataFrame(bt.genes, keep.extra.columns=TRUE)


#### collecting GAINS
Gdf <- df[df$TCN>2,]
rownames(Gdf)=NULL
Gdf <- makeGRangesFromDataFrame(Gdf,seqnames.field='Chr', keep.extra.columns=TRUE)

overlappp <- findOverlaps(Gdf, bt.genes,
         maxgap=-1, minoverlap=100,
         type="any",    
         select="all",
         ignore.strand=FALSE)

indexes<-1:length(Gdf)
CNA <- data.frame()
for (seq in unique(Gdf$sample)) {
	modifG <- unique(bt.genes[subjectHits(overlappp[queryHits(overlappp) %in% indexes[Gdf$sample==seq]])]$gene_name)
	print(seq)
	if (length(modifG)>0){
		supDf <- data.frame(genes=modifG)
		supDf$sample <- seq
		supDf$isto <- unique(Gdf[Gdf$sample==seq]$Istotype)
		supDf$Status <- 'Amp'
		CNA <- rbind(CNA, supDf)
		
	}
}

Gdf <- df[df$TCN<2,]
rownames(Gdf)=NULL
Gdf <- makeGRangesFromDataFrame(Gdf,seqnames.field='Chr', keep.extra.columns=TRUE)

overlappp <- findOverlaps(Gdf, bt.genes,
         maxgap=-1, minoverlap=100,
         type="any",    
         select="all",
         ignore.strand=FALSE)

indexes<-1:length(Gdf)
for (seq in unique(Gdf$sample)) {
	modifG <- unique(bt.genes[subjectHits(overlappp[queryHits(overlappp) %in% indexes[Gdf$sample==seq]])]$gene_name)
	print(seq)
	if (length(modifG)>0){
		supDf <- data.frame(genes=modifG)
		supDf$sample <- seq
		supDf$isto <- unique(Gdf[Gdf$sample==seq]$Istotype)
		supDf$Status <- 'Del'
		CNA <- rbind(CNA, supDf)
		
	}
}


write.table(file='CNA_Annotated.csv', CNA, row.names=F)


