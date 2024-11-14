options(java.parameters = "-Xmx20000m")

library(xlsx)
library(dplyr)
library(ggplot2)

biop <- c('PRJNA752630')
tt <- c('B-cell Lymphoma')

df <- data.frame()
for (n in 1:length(biop)){
	subFiles <- list.files(paste0(biop[n],'/TCN/'))
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

df <- df[df$TCN!=2,]
alteredSample <- df %>% group_by(sample, Istotype) %>% summarize( ALT_L = sum(Width))
genomeLeng <- sum(read.table('canine_cytoband_for_gistic_order.txt')[1:38,'V5'])
alteredSample$cnaBurden <- alteredSample$ALT_L/genomeLeng

results <- data.frame()
for (isto in unique(alteredSample$Istotype)){
	dfToWork <- subset(alteredSample, Istotype==isto)
	df2 <- quantile(dfToWork$cnaBurden)
	df2 <- t(data.frame(df2))
	rownames(df2)=isto
	results <- rbind(results,df2)
}

meansss <- results['50%']
meansss$isto <- rownames(meansss)
colnames(meansss)[1] <- 'mean'
meansss <- meansss[order(meansss$mean),]
alteredSample <- as.data.frame(alteredSample)

#### PREPARING THE PLOT

new_TT <- vector()
new_frame <- data.frame()
new_pos <- vector()
start_pt <- vector()
end_pt <- vector()

c=0
for (l in meansss$isto) {
	new_frame1 <- subset(alteredSample, Istotype==l)
	new_frame1 <- new_frame1[order(new_frame1$cnaBurden),]
	new_frame <- rbind(new_frame, new_frame1)
	new_pos <- c(new_pos, seq(from = c+0.25, to = c+0.75, length.out = length(new_frame1$Istotype)))
	new_TT <- c( new_TT , rep(l, length(new_frame1$Istotype)))
	start_pt <- c(start_pt, c + 0.25)
	end_pt <- c(end_pt, c+0.75)
	c=c+1
}

df_tmb <- data.frame(x=new_pos ,y=new_frame$cnaBurden, tum=new_TT)
n <- 1:length(meansss$isto)

rects <- data.frame(xstart=seq(0,9), xend=seq(1,10), col=factor(rep(c(0,1),5)))

gg <- ggplot() 
gg <- gg + geom_rect(data=rects, aes(ymin=0, ymax=Inf, xmin=xstart, xmax=xend, fill=col), alpha=0.4) + scale_fill_manual(values=c('white', 'grey50')) + theme_bw()
gg <- gg + geom_point(data=df_tmb, aes(x, y)) + scale_y_log10() + scale_x_continuous( breaks= 1:length(unique(df_tmb$tum)) + 0.5, labels=meansss$isto)  + labs(x=element_blank(), y='CNA Burden') 
gg <- gg + geom_segment(aes(x = start_pt[n], xend = end_pt[n] , y = meansss$mean[n] , yend = meansss$mean[n]))  + theme(axis.text.x =element_text(size=15,angle = 45, hjust=1,face = "bold"),
	axis.text.y =element_text(size=15,face = "bold"),
	axis.title.y =element_text(size=20, face='bold')) + coord_cartesian(xlim=c(0.5,9.5)) + theme(legend.position='none')

jpeg(filename='cnaBurden2.jpeg', height=1200, width=3000, res=180)
gg
dev.off()



