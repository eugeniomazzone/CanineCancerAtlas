library(maftools)
library(viridis)
library(ggplot2)
library(dplyr)

#### CREATING DIR FOR LATER

dir.create("CoocSmallMut")
dir.create("CoocMAF")

#### DEFINING THE BIOPROJECTS

biop <- c('/Annotation/SNP/')
tt <- c('B-cell Lymphoma')
cds <- c(57)

#### READING AND LOADING ALL MAF FILES

mafs <- data.frame()
clinical <- data.frame()
for (n in 1:length(biop)) {

path2<-biop[n]
files2 <- list.files(path2)
files2=files2[!grepl(files2, pattern='outQC|meta|dups')]
mafs2 <- annovarToMaf(paste(path2,files2,sep=''), table='ensGene',ens2hugo = FALSE, refBuild="canFam3")
mafs2$Hugo_Symbol <- mafs2$Gene.refGene
mafs2 <- mafs2[mafs2$Func.refGene!='intronic']
mafs2 <- mafs2[mafs2$Func.refGene!='intergenic']
mafs2 <- mafs2[mafs2$Func.refGene!='UTR3']
mafs2 <- mafs2[mafs2$Func.refGene!='UTR5']
mafs2 <- mafs2[mafs2$Func.refGene!='downstream']
mafs2 <- mafs2[mafs2$Func.refGene!='upstream']
mafs2 <- mafs2[mafs2$Func.refGene!='nRNA_intronic']
mafs2 <- mafs2[mafs2$Func.refGene!='nRNA_splicing']
mafs2 <- mafs2[mafs2$Func.refGene!='nRNA_exonic']
mafs2 <- mafs2[mafs2$Func.refGene!='ncRNA_exonic']
mafs2 <- mafs2[mafs2$Func.refGene!='ncRNA_intronic']
mafs2 <- mafs2[mafs2$ExonicFunc.refGene!='synonymous SNV']
mafs2 <- mafs2[!grepl(mafs2$Chromosome, pattern='chrUn|chrM|chrX'),]

cli.data2 <- table(rep(biop[n], length(files2)) , unlist(strsplit(files2, split='.vcf.canFam3_multianno.txt', fixed=TRUE)), rep(tt[n], length(files2)))
cli.data2 <- as.data.frame(cli.data2)[1:3]
colnames(cli.data2) = c('Project','Tumor_Sample_Barcode', 'TumorType')

mafs <- rbind(mafs,mafs2)
clinical <- rbind(clinical, cli.data2)
}

mafs <-mafs[!grepl('ENSCAF',mafs$Hugo_Symbol),]
mafs <-mafs[!grepl('TTN',mafs$Hugo_Symbol),]

full_df <- merge(mafs, clinical, by='Tumor_Sample_Barcode')
#### Making coocurrence for each istotype

for(tt in unique(full_df$TumorType)) {
mafsToWork <- subset(full_df, TumorType==tt)
objMaf <- read.maf(mafsToWork)
jpeg(filename=paste0('CoocSmallMut/CoocSmall',tt,'.jpeg'), width=3000, height=3000, res=300)
soma <- somaticInteractions(maf = objMaf, pvalue = c(0.05, 0.1), fontSize=0.6)
dev.off()
write.csv(file=paste0('CoocMAF/cooc_',tt,'.csv'), soma, row.names=FALSE)

}






