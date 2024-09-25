library(BSgenome)
library(EaCoN)

setwd("CNA")
args = commandArgs(trailingOnly=TRUE)
biop='../BQRS'

#BINpack.Maker(bed.file = "../../bed_refined.bed", bin.size = 150, genome.pkg = "BSgenome.Cfamiliaris.UCSC.canFam3") ### RUN ONE TIME FOR NEW BEDFILE

try(WES.Bin(testBAM = paste(biop,'/', args[2], sep=''), refBAM = paste(biop,'/', args[1], sep=''), BINpack = "../refined_bed_canFam3_b150.GC.rda", samplename = args[2]))

try(WES.Normalize.ff(BIN.RDS.file = paste(args[2],'/',args[2],'_canFam3_b150_binned.RDS', sep=''), BINpack = "../refined_bed_canFam3_b150.GC.rda"))

try(Segment.ff(RDS.file = paste(args[2],'/',args[2],'_canFam3_b150_processed.RDS', sep=''), segmenter = "ASCAT"))

try(ASCN.ff(RDS.file = paste(args[2],'/ASCAT/L2R/',args[2],'.SEG.ASCAT.RDS', sep='')))
