options(java.parameters = "-Xmx7000m")

library(dplyr)
library(GenomicRanges)
library(rtracklayer)

tt <- c('B-Cell Lymphoma')

#### READING DATA
df <- data.frame()

subFiles <- list.files(paste0('TCN/'))
for (samp in subFiles) {
	sample <- paste0('TCN/',samp)
	openStat <- read.table(sample, header=T)
	colnames(openStat)[1]='Sample'
	openStat$Istotype <- tt
	df <- rbind(df,openStat)
}

rownames(df)<-NULL
colnames(df)[9]='minor_cn'
colnames(df)[1]='sample'

#### READING GTF
bt.genes <- rtracklayer::import('/refFiles/Canis_lupus_familiaris.CanFam3.1.104.gtf.gz')
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

