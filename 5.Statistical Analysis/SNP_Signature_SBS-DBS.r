options(java.parameters = "-Xmx20000m")

library(maftools)
library(viridis)
library(ggplot2)
library(dplyr)
library(sigminer)
library(BSgenome.Cfamiliaris.UCSC.canFam3)


#### CREATING DIR FOR LATER

dir.create("SignatureDBS")
dir.create("SignatureSBS")

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

#### EXTRACTING HISTOTYPE MUTATIONS

mafs <-mafs[!grepl('ENSCAF',mafs$Hugo_Symbol),]
mafs <-mafs[!grepl('TTN',mafs$Hugo_Symbol),]
full_df <- merge(mafs, clinical, by='Tumor_Sample_Barcode')

#### We'll work only with BCL data

sigDf <- data.frame()
for (tp in unique(full_df$TumorType)){
	df <- subset(full_df, TumorType==tp)
	df$TumorType <- NULL
	df$Project <- NULL
	mafDf <- read_maf(df)

	mt_tally <- sig_tally(
	  mafDf,
	  ref_genome = "BSgenome.Cfamiliaris.UCSC.canFam3",
	  useSyn = TRUE, genome_build='hg19'
	)

	mt_sig <- sig_extract(mt_tally$nmf_matrix,
	  n_sig = 1,
	  nrun = 30,
	  cores = 2,
	)

	act_refit = sig_fit(t(mt_tally$all_matrices$SBS_96), sig_index = "ALL")
	act_refit2 = act_refit[apply(act_refit, 1, function(x) sum(x) > 0.1),]

	sim <- get_sig_similarity(mt_sig, sig_db='latest_SBS_GRCh37')
	sim <- sim$similarity
	rownames(sim) <- tp
	sigDf <- rbind(sigDf, sim)

	jpeg(filename=paste0('SignatureSBS/SignatureSBS_',tp,'.jpeg'), width=6000, height=3000, res=300)
	p<- show_sig_profile(mt_sig, mode = "SBS", paint_axis_text = FALSE, x_label_angle = 90)
	show(p)
	dev.off()
}

write.table(file='Sig_SBS_Similarity.csv',sigDf, row.names=FALSE)

###### DBS
sigDf <- data.frame()
for (tp in unique(full_df$TumorType)[c(1:6,8,9,10)]){
	df <- subset(full_df, TumorType==tp)
	df$TumorType <- NULL
	df$Project <- NULL
	mafDf <- read_maf(df)

	mt_tally <- sig_tally(
	  mafDf,
	  mode = "DBS",
	  ref_genome = "BSgenome.Cfamiliaris.UCSC.canFam3",
	  useSyn = TRUE, genome_build='hg19'
	)

	mt_sig <- sig_extract(mt_tally$nmf_matrix,
	  n_sig = 1,
	  nrun = 30,
	  cores = 2,
	)

	sim <- get_sig_similarity(mt_sig, sig_db='latest_DBS_GRCh38')
	sim <- sim$similarity
	rownames(sim) <- tp
	sigDf <- rbind(sigDf, sim)

	jpeg(filename=paste0('SignatureDBS/SignatureDBS_',tp,'.jpeg'), width=6000, height=3000, res=300)
	p<- show_sig_profile(mt_sig, mode = "DBS", paint_axis_text = FALSE, x_label_angle = 90)
	show(p)
	dev.off()
}


write.table(file='Sig_DBS_Similarity.csv',sigDf, row.names=FALSE)

