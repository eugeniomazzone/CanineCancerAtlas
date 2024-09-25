library(maftools)
library(viridis)
library(ggplot2)
library(dplyr)
library(writexl)

#### DEFINING THE BIOPROJECTS

biop <- c('PRJNA752630')
tt <- c('B-cell Lymphoma')
cds <- c(57)

#### READING AND LOADING ALL MAF FILES

mafs <- data.frame()
clinical <- data.frame()
for (n in 1:length(biop)) {

path2<-paste(biop[n],'/SNP/',sep='')
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

#### EXTRACTING HISTOTYPE MUTATIONS

mafs <-mafs[!grepl('ENSCAF',mafs$Hugo_Symbol),]
mafs <-mafs[!grepl('TTN',mafs$Hugo_Symbol),]
pp <- dplyr::count(mafs, Hugo_Symbol, group_by=Tumor_Sample_Barcode, sort = TRUE)
pp$Histotype <- sapply(pp$group_by, function(samp) clinical$TumorType[clinical$Tumor_Sample_Barcode==samp])
pp1 <- dplyr::count(pp, Hugo_Symbol, group_by=Histotype, sort = TRUE)
pp_lost <- pp1[pp1$group_by %in% c('Pulmonary Adenocarcinoma', 'Urinary Carcinoma')]
for_coo <- rbind(pp1[pp1$n>=5], pp_lost)
colnames(for_coo)[2]='Histotype'

gene_order <- for_coo %>% dplyr::count(by=Hugo_Symbol, sort=T)
write.table(file='gene_by_histo.csv',gene_order, row.names=FALSE)

#### COUNT MUT*HISTO

gene_by_histo = dplyr::count(for_coo, group_by = Histotype)
colnames(gene_by_histo)[1]='Histotype'

#### MAKE TMB

meansss <- vector()

for (l in levels(clinical$TumorType)) {
	new_samples <- as.vector(clinical$Tumor_Sample_Barcode[clinical$TumorType==l])
	new_bp <- as.vector(clinical$Project[clinical$TumorType==l])
	new_count <- apply(as.data.frame(new_samples), 1, function(r) sum(mafs$Tumor_Sample_Barcode==r))

	new_frame1 <- data.frame(new_count,new_samples, new_bp)

	for (bp in levels(as.factor(new_bp))) {
		new_frame1$new_count[new_frame1$new_bp==bp] <- sapply(new_frame1$new_count[new_frame1$new_bp==bp], function(r) (r/cds[biop==bp]))
	}
	meansss <- c(meansss, median(new_frame1$new_count))
}

newf <- data.frame(TumorType=levels(clinical$TumorType), means=meansss)
colnames(newf) = c('Histotype', 'TMB')

#### CORRELATION 

cor.test(newf$TMB, gene_by_histo$n)
scat_pt <- data.frame('Histotypes'=newf$Histotype, 'TMB'=newf$TMB, 'mutated genes'=gene_by_histo$n)

gg <-ggplot(scat_pt, aes(y=TMB, x=mutated.genes)) + geom_point() + geom_text_repel(
    aes(label=Histotypes), size=4, box.padding = unit(0.5, "lines")
  ) + labs(x="Number of mutated genes", y="Tumor Mutational Burden")+ theme(axis.text.x=element_text(size=15,face = "bold"),
        axis.text.y =element_text(size=15,face = "bold"),
        axis.title.y =element_text(size=15,face = "bold"), axis.title.x = element_text(size=15,face = "bold"))

tiff(filename='corr.tiff', height=2000, width=2000, res=300)
gg
dev.off()

#### HEATMAP

last_df <- for_coo %>% pivot_wider(names_from=Histotype, values_from=n )
last_df <- replace(last_df, is.na(last_df), 0)
last_df <- data.frame(last_df)
rownames(last_df) <- last_df$Hugo_Symbol
last_df$Hugo_Symbol <- NULL
support <- dplyr::count(clinical, group_by=TumorType)
colnames(support)[1]='Histotype'
rownames(support) <- support$Histotype

gene_order <- dplyr::count(for_coo, group_by=Hugo_Symbol)
gene_order <- arrange(gene_order, -n)

to_plot <- last_df[gene_order$group_by[1:20],]

to_plot$genes <- rownames(to_plot)
melted_plot <- melt(to_plot)
melted_plot$variable <- as.character(melted_plot$variable)

melted_plot[grepl(melted_plot$variable, pattern='B.cell.Lymphoma'),]$variable = 'B-cell Lymphoma'

melted_plot$variable <- factor(melted_plot$variable)

for (i in 1:nrow(support)){
	melted_plot[grepl(melted_plot$variable, pattern=support$Histotype[i]),]$value <- melted_plot[grepl(melted_plot$variable, pattern=support$Histotype[i]),]$value/support[i,]$n
}


gg <- ggplot(data = melted_plot, aes(x=variable, y=genes, fill=value)) + geom_tile() + labs(x="Histotypes", y="Genes") + theme(axis.text.x=element_text(size=15,face = "bold",angle=45, hjust=1),
        axis.text.y =element_text(size=15,face = "bold"),
        axis.title.y =element_text(size=15,face = "bold"), axis.title.x = element_text(size=15,face = "bold"))

tiff(filename='heat.tiff', height=2000, width=2000, res=300)
gg
dev.off()



