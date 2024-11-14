NewsomaticInteractions <- function(maf, top = 25, genes = NULL, pvalue = c(0.05, 0.01), 
    returnAll = TRUE, geneOrder = NULL, fontSize = 0.8, leftMar = 4, 
    topMar = 4, showSigSymbols = TRUE, showCounts = FALSE, countStats = "all", 
    countType = "all", countsFontSize = 0.8, countsFontColor = "black", 
    colPal = "BrBG", revPal = FALSE, showSum = TRUE, plotPadj = FALSE, 
    colNC = 9, nShiftSymbols = 5, sigSymbolsSize = 2, sigSymbolsFontSize = 0.9, 
    pvSymbols = c(46, 42), limitColorBreaks = TRUE, isAmp=TRUE, isDel=TRUE) 
{
    if (is.null(genes)) {
        genes = getGeneSummary(x = maf)[1:top, Hugo_Symbol]
    }
    if (length(genes) < 2) {
        stop("Minimum two genes required!")
    }
    om = createOncoMatrix(m = maf, g = genes, isAmp=isAmp, isDel=isDel)
    #print(om)
    all.tsbs = levels(getSampleSummary(x = maf)[, Tumor_Sample_Barcode])
    mutMat = t(om$numericMatrix)
    missing.tsbs = all.tsbs[!all.tsbs %in% rownames(mutMat)]
    if (nrow(mutMat) < 2) {
        stop("Minimum two genes required!")
    }
    mutMat[mutMat > 0] = 1
    if (length(missing.tsbs) > 0) {
        missing.tsbs = as.data.frame(matrix(data = 0, nrow = length(missing.tsbs), 
            ncol = ncol(mutMat)), row.names = missing.tsbs)
        colnames(missing.tsbs) = colnames(mutMat)
        mutMat = rbind(mutMat, missing.tsbs)
    }
    interactions = sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat), 
        function(j) {
            f = try(fisher.test(mutMat[, i], mutMat[, j]), silent = TRUE)
            if (class(f) == "try-error") {
                if (all(mutMat[, i] == mutMat[, j])) {
                  if (colnames(mutMat)[i] != colnames(mutMat)[j]) {
                    warning("All the samples are in the same direction for the genes ", 
                      colnames(mutMat)[i], " and ", colnames(mutMat)[j], 
                      "! Could not perform Fisher test.")
                  }
                  NA
                }
                else {
                  if (colnames(mutMat)[i] != colnames(mutMat)[j]) {
                    warning("Contigency table could not created for the genes ", 
                      colnames(mutMat)[i], " and ", colnames(mutMat)[j], 
                      "! Could not perform Fisher test.")
                  }
                  NA
                }
            }
            else {
                ifelse(f$estimate > 1, -log10(f$p.val), log10(f$p.val))
            }
        }))
    oddsRatio <- oddsGenes <- sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat), 
        function(j) {
            f = try(fisher.test(mutMat[, i], mutMat[, j]), silent = TRUE)
            if (class(f) == "try-error") 
                if (all(mutMat[, i] == mutMat[, j])) {
                  NA
                }
                else {
                  NA
                }
            else f$estimate
        }))
    rownames(interactions) = colnames(interactions) = rownames(oddsRatio) = colnames(oddsRatio) = colnames(mutMat)
    sigPairs = which(x = 10^-abs(interactions) < 1, arr.ind = TRUE)
    sigPairs2 = which(x = 10^-abs(interactions) >= 1, arr.ind = TRUE)
    if (nrow(sigPairs) < 1) {
        stop("No meaningful interactions found.")
    }
    sigPairs = rbind(sigPairs, sigPairs2)
    sigPairsTbl = data.table::rbindlist(lapply(X = seq_along(1:nrow(sigPairs)), 
        function(i) {
            x = sigPairs[i, ]
            g1 = rownames(interactions[x[1], x[2], drop = FALSE])
            g2 = colnames(interactions[x[1], x[2], drop = FALSE])
            tbl = as.data.frame(table(factor(apply(X = mutMat[, 
                c(g1, g2), drop = FALSE], 1, paste, collapse = ""), 
                levels = c("00", "01", "11", "10"))))
            combn = data.frame(t(tbl$Freq))
            colnames(combn) = tbl$Var1
            pval = 10^-abs(interactions[x[1], x[2]])
            fest = oddsRatio[x[1], x[2]]
            d = data.table::data.table(gene1 = g1, gene2 = g2, 
                pValue = pval, oddsRatio = fest)
            d = cbind(d, combn)
            d
        }), fill = TRUE)
    sigPairsTbl[, `:=`(pAdj, p.adjust(pValue, method = "fdr"))]
    sigPairsTbl[is.na(sigPairsTbl)] = 0
    sigPairsTbl$Event = ifelse(test = sigPairsTbl$oddsRatio > 
        1, yes = "Co_Occurence", no = "Mutually_Exclusive")
    sigPairsTbl$pair = apply(X = sigPairsTbl[, .(gene1, gene2)], 
        MARGIN = 1, FUN = function(x) paste(sort(unique(x)), 
            collapse = ", "))
    sigPairsTbl[, `:=`(event_ratio, `01` + `10`)]
    sigPairsTbl[, `:=`(event_ratio, paste0(`11`, "/", event_ratio))]
    sigPairsTblSig = sigPairsTbl[order(as.numeric(pValue))][!duplicated(pair)]
        sigPairsTblSig = sigPairsTblSig[!gene1 == gene2]
    if (!returnAll) {
        sigPairsTblSig = sigPairsTblSig[pValue < min(pvalue)]
    }
    return(sigPairsTblSig)
}


createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE, cbio = FALSE, isAmp=TRUE, isDel=TRUE){

  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }

  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)

  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

      vc = c("")
      names(vc) = 0

      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }

  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
  }


  cnv_events = c(c("Amp", "Del"), as.character(subMaf[Variant_Type == "CNV"][, .N, Variant_Classification][, Variant_Classification]))
  cnv_events = unique(cnv_events)

  if(cbio){
    vc = c("Nonstop_Mutation", "Frame_Shift_Del", "Missense_Mutation",
           "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins")
    vc.cbio = c("Truncating", "Truncating", "Missense", "Truncating", "Truncating", "Truncating",
                "In-frame", "In-frame")
    names(vc.cbio) = vc
    subMaf[,Variant_Classification_temp := vc.cbio[as.character(subMaf$Variant_Classification)]]
    subMaf$Variant_Classification_temp = ifelse(test = is.na(subMaf$Variant_Classification_temp), yes = as.character(subMaf$Variant_Classification), no = subMaf$Variant_Classification_temp)
    subMaf[,Variant_Classification := as.factor(as.character(Variant_Classification_temp))]
    subMaf[,Variant_Classification_temp := NULL]
  }

if (isAmp && isDel){
  oncomat = data.table::dcast(data = subMaf[subMaf$Variant_Classification %in% c('Amp','Del'),.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x, cnv = cnv_events){
                                #x = unique(as.character(x)) #>=2 distinct variant classification = Multi_Hit
                                x = as.character(x) # >= 2 same/distinct variant classification = Multi_Hit See #347
                                xad = x[x %in% cnv]
                                xvc = x[!x %in% cnv]

                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                }

                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)
                              }
if (isAmp && !isDel){
  oncomat = data.table::dcast(data = subMaf[subMaf$Variant_Classification %in% c('Amp'),.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x, cnv = cnv_events){
                                #x = unique(as.character(x)) #>=2 distinct variant classification = Multi_Hit
                                x = as.character(x) # >= 2 same/distinct variant classification = Multi_Hit See #347
                                xad = x[x %in% cnv]
                                xvc = x[!x %in% cnv]

                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                }

                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)
                              }

if (!isAmp && isDel){
  oncomat = data.table::dcast(data = subMaf[subMaf$Variant_Classification %in% c('Del'),.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x, cnv = cnv_events){
                                #x = unique(as.character(x)) #>=2 distinct variant classification = Multi_Hit
                                x = as.character(x) # >= 2 same/distinct variant classification = Multi_Hit See #347
                                xad = x[x %in% cnv]
                                xvc = x[!x %in% cnv]

                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                }

                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)
                              }
  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])

  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)

  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)

  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }

  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }

  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)


  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId

    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]

    mdf = mdf[, -ncol(mdf)]

    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort

    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy

    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes, cnvc = cnv_events))
  }
}
