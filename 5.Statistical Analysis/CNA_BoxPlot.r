options(java.parameters = "-Xmx20000m")

library(xlsx)
library(dplyr)
library(ggplot2)
library(viridis)

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

results <- data.frame()

dfToWork <- df
isto <- 'ALL'
dfToWork <-dfToWork[dfToWork$TCN!=2,]

delss <- dfToWork[dfToWork$TCN<2,] %>% group_by(sample) %>% summarise(mean = mean(Width))
gainss <- dfToWork[dfToWork$TCN>2,] %>% group_by(sample) %>% summarise(mean = mean(Width))

var.test(delss$mean, gainss$mean, alternative="less")

gainss$Status <- 'Gain'
delss$Status <- 'Del'
my_data <- rbind(gainss, delss)
my_data <- as.data.frame(my_data)

jpeg(filename='WidthCNA.jpeg', height=2000, width=1400, res=180)
my_data %>%
  ggplot( aes(x=Status, y=mean, fill=Status)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("CNA Width Distribution") +
    xlab("") + theme(axis.text.x =element_text(size=15,face = "bold"),
	axis.text.y =element_text(size=15,face = "bold"),
	axis.title.y =element_text(size=20, face='bold'))
dev.off()


