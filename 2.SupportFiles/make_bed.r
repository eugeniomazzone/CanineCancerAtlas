library(GenomicRanges)
library(CNVRanger)
#library(AnnotationHub)
library(writexl) 
library(rtracklayer)
library(regioneR)
library(BSgenome.Cfamiliaris.UCSC.canFam3.masked)


gtf_my <- rtracklayer::import('canFam3.ensGene.gtf.gz')
gtf_my <- subset(gtf_my, type=='exon')
gtf_my <- subset(data.frame(gtf_my), select=c(seqnames, start, end))
gtf_my <- gtf_my[!grepl(gtf_my$seqnames, pattern='chrUn|chrM|chrX'),]
gtf_my <- gtf_my[-gtf_my$start + gtf_my$end !=0,]
rownames(gtf_my) <- NULL

gtf_my <- GenomicRanges::reduce(GenomicRanges::GRanges(gtf_my), min.gapwidth=100L)
gtf_my <- subset(data.frame(gtf_my), select=c(seqnames, start, end))

c=0
for (i in 1:nrow(gtf_my)){
c=c+gtf_my[i,'end'] - gtf_my[i,'start']
}
print(c/1000000)



write.table(gtf_my, file='refined_bed.bed', row.names=FALSE, sep='\t', quote = FALSE)
