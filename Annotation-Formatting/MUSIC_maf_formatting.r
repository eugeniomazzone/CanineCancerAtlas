library(maftools)

#### READING AND LOADING ALL MAF FILES

mafs <- data.frame()
clinical <- data.frame()

path2<-paste('SNP/',sep='')
files2 <- list.files(path2)
files2=files2[!grepl(files2, pattern='outQC')]
files2=files2[!grepl(files2, pattern='meta')]
files2=files2[!grepl(files2, pattern='dups')]
mafs2 <- annovarToMaf(paste(path2,files2,sep=''), table='ensGene',ens2hugo = FALSE, refBuild="canFam3")
mafs2$Hugo_Symbol <- mafs2$Gene.refGene
mafs2 <- mafs2[!grepl(mafs2$Chromosome, pattern='chrUn|chrM|chrX'),]

cli.data2 <- table(rep("PRJNA752630",length(files2)),rep("B-Cell Lymphomas", length(files2)))
cli.data2 <- as.data.frame(cli.data2)[1:3]
colnames(cli.data2) = c('Project','Tumor_Sample_Barcode', 'TumorType')

mafs <- rbind(mafs,mafs2)
clinical <- rbind(clinical, cli.data2)

mafs <-mafs[!grepl('ENSCAF',mafs$Hugo_Symbol),]
mafs <-mafs[!grepl('TTN',mafs$Hugo_Symbol),]

col16 <- do.call(rbind,strsplit(mafs$Tumor_Sample_Barcode, '-'))[,2]
col1 <- mafs[,7]
col9 <- mafs[,8]

mafs[,1] <- col1
mafs[,9] <- col9
mafs[,16] <- col16

write.table(mafs, file='music_maf.txt', row.names=F, quote=F, sep="\t", col.names=F)

