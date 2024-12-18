library(maftools)
library(viridis)

biop <- c('/Annotation/SNP/')
tt <- c('B-cell Lymphoma')

mafs <- data.frame()
clinical <- data.frame()

for (n in 1:length(biop)) {

path2<-biop[n]
files2 <- list.files(path2, pattern='canFam3')
files2=files2[!grepl(files2, pattern='outQC')]
files2=files2[!grepl(files2, pattern='meta')]
files2=files2[!grepl(files2, pattern='dups')]
mafs2 <- annovarToMaf(paste(path2,files2,sep=''), table='ensGene',ens2hugo = FALSE, refBuild="canFam3")
mafs2$Hugo_Symbol <- mafs2$Gene.refGene

cli.data2 <- table(rep(biop[n], length(files2)) , unlist(strsplit(files2, split='.vcf.canFam3_multianno.txt', fixed=TRUE)), rep(tt[n], length(files2)))
cli.data2 <- as.data.frame(cli.data2)[1:3]
colnames(cli.data2) = c('Project','Tumor_Sample_Barcode', 'Histotypes')

mafs <- rbind(mafs,mafs2)
clinical <- rbind(clinical, cli.data2)
}

mafs <- mafs[!grepl(mafs$Chromosome, pattern='chrUn|chrM|chrX'),]

##--------------------------------------------------------------------------------------

mafs <-mafs[!grepl('ENSCAF',mafs$Hugo_Symbol),]
objMaf <- read.maf(maf=mafs, clinicalData=clinical)

tiff(filename='OncoplotSNP.tiff', width=6000, height=3000, res=300)
oncoplot(maf=objMaf,gene_mar = 10, top=30, fontSize = 1.1, clinicalFeature =c('Histotypes'), sortByAnnotation = TRUE, genesToIgnore='TTN', showTitle=FALSE,  legendFontSize = 2, annotationFontSize = 2, sepwd_samples=0, drawColBar = FALSE, draw_titv = TRUE)
dev.off()


