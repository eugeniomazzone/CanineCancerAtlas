options(java.parameters = "-Xmx28000m")

library(maftools)
library(viridis)
library(ggplot2)
library(dplyr)
library(writexl)
library(xlsx)
library(dndscv)

load('RefCDS_dog_CanFam3.1.rda')

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

full_df2 <- full_df[full_df$ExonicFunc.refGene!='synonymous SNV']
pp <- dplyr::count(full_df2, Hugo_Symbol, group_by=Tumor_Sample_Barcode, sort = TRUE)
pp$Histotype <- sapply(pp$group_by, function(samp) clinical$TumorType[clinical$Tumor_Sample_Barcode==samp])
pp1 <- dplyr::count(pp, Hugo_Symbol, group_by=Histotype, sort = TRUE)
pp_lost <- pp1[pp1$group_by %in% c('Pulmonary Adenocarcinoma', 'Urinary Carcinoma')]
pp1 <- rbind(pp1[pp1$n>=5], pp_lost)

wb <- createWorkbook()  
for(tp in unique(pp1$group_by))
{
	print(tp)
	new_df <- subset(pp1, group_by==tp)
	message("Creating sheet", tp)
	sheet <- createSheet(wb, sheetName = tp)
	message("Adding data frame", tp)
	addDataFrame(new_df, sheet, row.names = F)
}
saveWorkbook(wb, "Recuccrent_Small_Mut_by_Histo.xlsx")  

full_df <- subset(full_df, select=c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'Reference_Allele','Tumor_Seq_Allele2', 'TumorType'))
colnames(full_df) = c("sampleID", "chr", "pos", "ref", "mut", "TumorType")
full_df$chr <- do.call(rbind,strsplit(full_df$chr, 'chr'))[,2]


wb <- createWorkbook()  
for(tp in unique(tt[!grepl(tt, pattern='Pulmonary Adenocarcinoma|Urinary Carcinoma')]))
{
	print(tp)
	new_df <- subset(full_df, TumorType==tp)[,1:5]
	dndsout = dndscv(new_df, refdb='RefCDS_dog_CanFam3.1.rda', cv=NULL, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample=50000)
	message("Creating sheet", tp)
	sheet <- createSheet(wb, sheetName = tp)
	message("Adding data frame", tp)
	addDataFrame(subset(dndsout$sel_cv, pglobal_cv<=0.01), sheet, row.names = F)
}
saveWorkbook(wb, "Small_Driver_Mut_Histo.xlsx")  


wb <- createWorkbook()  
new_df <- full_df[,1:5]
dndsout = dndscv(new_df, refdb='RefCDS_dog_CanFam3.1.rda', cv=NULL, , max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample=50000)
message("Creating sheet", 'Candidate Driver')
sheet <- createSheet(wb, sheetName = 'Candidate Driver')
message("Adding data frame", 'Candidate Driver')
addDataFrame(subset(dndsout$sel_cv, pglobal_cv<=0.05), sheet, row.names = F)

dndsout <- subset(dndsout$sel_cv, qglobal_cv<=0.05)
dndsout$gene_name <- do.call(rbind,strsplit(dndsout$gene_name, '\\.'))[,2]

dndsout$gene_name[dndsout$gene_name=='K-RAS']='KRAS'
dndsout$gene_name[dndsout$gene_name=='N-RAS']='NRAS'
dndsout$gene_name[dndsout$gene_name=='P53']='TP53'
dndsout$gene_name[dndsout$gene_name=='ENSCAFG00000013622']='IGLV2-33'
results <- data.frame(dndsout$gene_name)

for (i in 1:nrow(results)){
for (j in unique(tt)){
print(i)
print(j)
print(sum(grepl(full_df2[full_df2$TumorType==j,]$Hugo_Symbol, pattern=results[i,'dndsout.gene_name'])))
results[i,j] = sum(grepl(full_df2[full_df2$TumorType==j,]$Hugo_Symbol, pattern=results[i,]))
}
}

message("Creating sheet", 'Driver Details')
sheet <- createSheet(wb, sheetName = 'Driver Details')
message("Adding data frame", 'Driver Details')
addDataFrame(results, sheet, row.names = F)
saveWorkbook(wb, "Small_Driver_Mut.xlsx") 

