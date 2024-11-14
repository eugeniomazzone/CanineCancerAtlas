library(BSgenome)
library(EaCoN)

dir.create("CNA")
setwd("CNA")
args = commandArgs(trailingOnly=TRUE)
biop='../BQRS'

bed <- "bed_refined.bed"
processed_bed <- paste0(unlist(strsplit(bed, '.bed'))[1],'_canFam3_b150.GC.rda')
if (file.exists(processed_bed)){
	print('BED already processed')
} else
	print('BED will be processed')
	BINpack.Maker(bed.file = paste0('../../',bed), bin.size = 150, genome.pkg = "BSgenome.Cfamiliaris.UCSC.canFam3") ### RUN ONE TIME FOR NEW BEDFILE
}

print(paste0('Start processing: ', args[2]))
try(WES.Bin(testBAM = paste(biop,'/', args[2], sep=''), refBAM = paste(biop,'/', args[1], sep=''), BINpack = paste0('../',processed_bed), samplename = args[2]))

try(WES.Normalize.ff(BIN.RDS.file = paste(args[2],'/',args[2],'_canFam3_b150_binned.RDS', sep=''), BINpack = paste0('../',processed_bed)))

try(Segment.ff(RDS.file = paste(args[2],'/',args[2],'_canFam3_b150_processed.RDS', sep=''), segmenter = "ASCAT"))

try(ASCN.ff(RDS.file = paste(args[2],'/ASCAT/L2R/',args[2],'.SEG.ASCAT.RDS', sep='')))
