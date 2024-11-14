library(maftools)
library(viridis)
library(ggplot2)
library(ggwordcloud)

#### DEFINING THE DATA FOLDERS

biop <- c('PRJNA752630')
tt <- c('B-cell Lymphoma')
cds <- c(57)

#### READING AND LOADING ALL MAF FILES

mafs <- data.frame()
clinical <- data.frame()
for (n in 1:length(biop)) {
	path2<-paste(biop[n],'/SNP/',sep='')
	files2 <- list.files(path2)
	files2=files2[!grepl(files2, pattern='outQC')]
	files2=files2[!grepl(files2, pattern='meta')]
	files2=files2[!grepl(files2, pattern='dups')]
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
	mafs2 <- mafs2[mafs2$ExonicFunc.refGene!='synonymous SNV']
	mafs2 <- mafs2[!grepl(mafs2$Chromosome, pattern='chrUn|chrM|chrX'),]

	cli.data2 <- table(rep(biop[n], length(files2)) , unlist(strsplit(files2, split='.vcf.canFam3_multianno.txt', fixed=TRUE)), rep(tt[n], length(files2)))
	cli.data2 <- as.data.frame(cli.data2)[1:3]
	colnames(cli.data2) = c('Project','Tumor_Sample_Barcode', 'TumorType')

	mafs <- rbind(mafs,mafs2)
	clinical <- rbind(clinical, cli.data2)
}

#### EVALUATING TMB

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


#### PREPARING THE PLOT

new_TT <- vector()
new_frame <- data.frame()
new_pos <- vector()
start_pt <- vector()
end_pt <- vector()

c=0
meansss <- vector()
for (l in newf[sort(newf$means, index.return=TRUE)$ix,]$TumorType) {

	new_samples <- as.vector(clinical$Tumor_Sample_Barcode[clinical$TumorType==l])
	new_bp <- as.vector(clinical$Project[clinical$TumorType==l])
	new_count <- apply(as.data.frame(new_samples), 1, function(r) sum(mafs$Tumor_Sample_Barcode==r))

	new_frame1 <- data.frame(new_count,new_samples, new_bp)

	for (bp in levels(as.factor(new_bp))) {
		new_frame1$new_count[new_frame1$new_bp==bp] <- sapply(new_frame1$new_count[new_frame1$new_bp==bp], function(r) (r/cds[biop==bp]))
	}
	new_frame <- rbind(new_frame,  new_frame1[order(new_frame1$new_count),])
	new_pos <- c(new_pos, seq(from = c+0.25, to = c+0.75, length.out = length(new_frame1$new_count)))
	new_TT <- c( new_TT , rep(l, length(new_count)))
	meansss <- c(meansss, median(new_frame1$new_count))
	start_pt <- c(start_pt, c + 0.25)
	end_pt <- c(end_pt, c+0.75)
	c=c+1
}

df_tmb <- data.frame(x=new_pos ,y=new_frame$new_count, tum=new_TT)
n <- 1:length(tt)

rects <- data.frame(xstart=seq(0,9), xend=seq(1,10), col=factor(rep(c(0,1),5)))

### CONSTRUCTIONG THE PLOT

gg <- ggplot() 
gg <- gg + geom_rect(data=rects, aes(ymin=0, ymax=Inf, xmin=xstart, xmax=xend, fill=col), alpha=0.4) + scale_fill_manual(values=c('white', 'grey50')) + theme_bw()
gg <- gg + geom_point(data=df_tmb, aes(x, y)) + scale_y_log10() + scale_x_continuous( breaks= 1:length(unique(newf$TumorType)) + 0.5, labels=newf[sort(newf$means, index.return=TRUE)$ix,]$TumorType)  + labs(x=element_blank(), y='Tumor Mutational Burden') 
gg <- gg + geom_segment(aes(x = start_pt[n], xend = end_pt[n] , y = meansss[n] , yend = meansss[n]))  + theme(axis.text.x =element_text(size=15,angle = 45, hjust=1,face = "bold"),
	axis.text.y =element_text(size=15,face = "bold"),
	axis.title.y =element_text(size=20, face='bold')) + coord_cartesian(xlim=c(0.5,9.5)) + theme(legend.position='none')

### WRITING THE PLOT

tiff(filename='TMB.tiff', height=1200, width=3000, res=180)
gg
dev.off()

