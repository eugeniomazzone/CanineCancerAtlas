install.packages('devtools')
library(devtools)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("affxparser", "Biostrings", "aroma.light", "BSgenome", "copynumber", "GenomicRanges", "limma", "rhdf5","BSgenome.Cfamiliaris.UCSC.canFam3"))

install.packages(c('pbapply', 'squash', 'iotools', 'readr', 'seqminer'))

devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")
devtools::install_github("mskcc/facets")
install_github("aroneklund/copynumber")

install.packages('sequenza_3.0.0.tar.gz', repos = NULL, type="source")

devtools::install_github("gustaveroussy/EaCoN")


