options(java.parameters = "-Xmx20000m")

library(xlsx)
library(sigminer)

#### CREATING DIR FOR LATER
dir.create("CNA_out")

#### DEFINING THE BIOPROJECTS
biop <- c('PRJNA752630')
tt <- c('B-cell Lymphoma')

sigDf <- data.frame()
wb <- createWorkbook()
for (n in 1:length(unique(tt))){
	df <- data.frame()
	for (Project in biop[tt==tt[n]]) {
		subFiles <- list.files(paste0(Project,'/TCN/'))
		for (samp in subFiles) {
			sample <- paste0(biop[n],'/TCN/',samp)
			openStat <- read.table(sample, header=T)
			colnames(openStat)[1]='Sample'
			openStat$Istotype <- tt[n]
			df <- rbind(df,openStat)
		}
	}
	colnames(df)[9]='minor_cn'
	colnames(df)[1]='sample'

	cn <- read_copynumber(df,
	  seg_cols = c("Chr", "Start", "End", "TCN"),
	  genome_measure = "wg", complement = TRUE, add_loh = TRUE
	)

	tally_s <- sig_tally(cn, method = "S")
	sig_denovo = sig_extract(tally_s$all_matrices$CN_48, n_sig=1)

	sim <- get_sig_similarity(sig_denovo, sig_db='latest_CN_GRCh37')
	sim <- sim$rss
	rownames(sim) <- tt[n]
	sigDf <- rbind(sigDf, sim)

	act_refit = sig_fit(t(tally_s$all_matrices$CN_48), sig_index = "ALL", sig_db = "CNS_TCGA")

	act_refit2 = act_refit[apply(act_refit, 1, function(x) sum(x) > 0.1),]

	jpeg(file=paste0('CNA_out/DeNovo_',tt[n],'.jpeg'), res=150, width=2400, height=2000)
	p<-show_sig_profile(sig_denovo, mode = "copynumber", method = "S", style = "cosmic")
	show(p)
	dev.off()

	jpeg(file=paste0('CNA_out/Exposure_',tt[n],'.jpeg'), res=150, width=2400, height=2000)
	p <- show_sig_exposure(act_refit2)
	show(p)
	dev.off()
	p
}

sheet <- createSheet(wb, sheetName = paste0('cn',tt[n],'List'))
addDataFrame(sigDf, sheet, row.names = T)
saveWorkbook(wb, "CNA_Sig_Similarity.xlsx") 


