#' Compute similarity between two mutation data
#'
#' @keywords internal
#' @param mut1 Input sample mutation data.
#' @param mut2 Refrence mutation data.
#'
#' @return Score of similarity

computeSim = function(mut1, mut2){
  mut1 = as.numeric(mut1)
  mut2 = as.numeric(mut2)
  similarity = sum(mut1*mut2)/(sqrt(sum(mut1**2)) * sqrt(sum(mut2**2)))
}

#'  Transform the input data to simSig input data
#'
#' @param tmMut Location of the mutation file that is to be converted or name of data frame in environment.
#' Use data(sample.mut.ref) to see the example data.
#' @param bsg Location of the mutation file that is to be converted or name of data frame in environment.
#' The default value is BSgenome.Hsapiens.UCSC.hg19, and Only set if another genome build is required. Must be a BSgenome object.
#'
#' @return A data frame that contains sample IDs for the rows and trinucleotide contexts for the columns.
#' Each entry is the count of how many times a mutation with that trinucleotide context is seen in the sample.
#' @export
#'
#' @examples testTrans = transMut(tmMut = sample.mut.ref)
transMut = function(tmMut, bsg = NULL){
  tmMut = mut.to.sigs.input(mut.ref = tmMut,
                            bsg = bsg)
}
#' Get the similarity score of a tumor sample with mSignaturedb
#'
#' @param trTmMut Either a data frame or location of input text file, where rows are samples, columns are trinucleotide contexts.
#' Use data(sample.mut.trans) to see the example data.
#' @param tmMutRef Location of the mutation file that is to be converted or name of data frame in environment.
#' "tmMutRef" is the standard signature value for one site and one country.
#' Default value is "mSigdb" that convert from "mSignatureDB" database.
#' USe data("mSigdb") to see the example dasta.
#' @param clinicalData Location of the mutation file that is to be converted or name of data frame in environment.
#' Include country and site information of every sample in trTmMut.
#' Use data("sample.clinical") to see the example data.
#' @param contexts.needed FALSE if tumor.file is a context file, TRUE if it is only mutation counts.
#' @param tri.counts.method Set to either:
#' \itemize{
#'  \item 'default' -- no further normalization \item 'exome' -- normalized by
#'   number of times each trinucleotide context is observed in the exome \item
#'   'genome' -- normalized by number of times each trinucleotide context is
#'   observed in the genome \item 'exome2genome' -- multiplied by a ratio of that
#'   trinucleotide's occurence in the genome to the trinucleotide's occurence in
#'   the exome \item 'genome2exome' -- multiplied by a ratio of that
#'   trinucleotide's occurence in the exome to the trinucleotide's occurence in
#'   the genome \item data frame containing user defined scaling factor -- count
#'   data for each trinucleotide context is multiplied by the corresponding value
#'   given in the data frame }
#'
#' @return The similarity score of a tumor sample with mSignaturedb
#' @export
#'
#' @examples testSim = sigSim(trTmMut = sample.mut.trans, clinicalData = sample.clinical, contexts.needed = TRUE)
sigSim = function(trTmMut, clinicalData, tmMutRef = mSigdb, contexts.needed = FALSE, tri.counts.method = 'default'){
  sampleList = rownames(trTmMut)
  bsample = nrow(trTmMut)
  clinicalData = clinicalData[(clinicalData$sampleId%in%sampleList), ]
  asample = nrow(clinicalData)
  if (bsample == asample) {
  }else{
    warning(paste(bsample - asample,"samples were not calculated because they don't have clinical data."), call. = FALSE)
  }
  sampleList = names(table(clinicalData$sampleID))
  siteList = names(table(mSigdb$tmSite))
  brow = nrow(clinicalData)
  clinicalData = clinicalData[(clinicalData$tmSite%in%siteList), ]
  arow = nrow(clinicalData)
  if (brow == arow) {
  }else{
    warning(paste(brow - arow,"samples were not calculated because their sites are not part of the tmMutRef."), call. = FALSE)
  }
  deconSigs = list()
  sampleList = as.character(clinicalData$sampleId)
  for (sample_id in sampleList) {
    deconSigs[[sample_id]] = deconstructSigs::whichSignatures(tumor.ref = trTmMut,
                                             sample.id = sample_id,
                                             contexts.needed = contexts.needed,
                                             tri.counts.method = tri.counts.method,
                                             signatures.ref = signatures.cosmic)
  }
  simScores = list()
  misCountry = 0
  for (sample_id in names(deconSigs)) {
    sample_mut = deconSigs[[sample_id]]
    sample_mut = as.numeric(unname(sample_mut$weights))
    sample_cli = unname(clinicalData[as.character(clinicalData$sampleId) == as.character(sample_id), ])
    sample_country = as.character(sample_cli[[3]])
    sample_Site = as.character(sample_cli[[2]])
    tmMutRef = tmMutRef[tmMutRef$tmSite == sample_Site, ]
    refCountry = as.character(tmMutRef$Country)
    if(sample_country%in%refCountry){
    }else{
      sample_country = refCountry[1]
      misCountry =  misCountry + 1
    }
    ref_mut = unname(unlist(tmMutRef[(tmMutRef$Country == sample_country & tmMutRef$tmSite == sample_Site), ]))[c(3:32)]  
    simScores[[sample_id]] = computeSim(sample_mut, ref_mut)
  }
  if (misCountry == 0) {
  }else{
    warning(paste(misCountry,"samples use the default country param because their countries are not part of the tmMutRef."), call. = FALSE)
  }
  simScores
}

#' Title Get the similarity scores of every sample with every tumor site
#'
#' @param contexts.needed FALSE if tumor.file is a context file, TRUE if it is only mutation counts.
#' @param tri.counts.method Set to either:
#' \itemize{
#'  \item 'default' -- no further normalization \item 'exome' -- normalized by
#'   number of times each trinucleotide context is observed in the exome \item
#'   'genome' -- normalized by number of times each trinucleotide context is
#'   observed in the genome \item 'exome2genome' -- multiplied by a ratio of that
#'   trinucleotide's occurence in the genome to the trinucleotide's occurence in
#'   the exome \item 'genome2exome' -- multiplied by a ratio of that
#'   trinucleotide's occurence in the exome to the trinucleotide's occurence in
#'   the genome \item data frame containing user defined scaling factor -- count
#'   data for each trinucleotide context is multiplied by the corresponding value
#'   given in the data frame }
#' @param trTmMut Either a data frame or location of input text file, where rows are samples, columns are trinucleotide contexts.
#' Use data(sample.mut.trans) to see the example data.
#' @param predictRef Location of the mutation file that is to be converted or name of data frame in environment.
#' "tmMutRef" is the standard signature value for one site.
#' Default value is "mSigsitedb" that convert from "mSignatureDB" database.
#' USe data("mSigsitedb") to see the example dasta.
#'
#' @return A dataframe contain similarity scores of every sample with every tumor site
#' @export
#'
#' @examples testPredict = predictTm(trTmMut = sample.mut.trans, contexts.needed = TRUE)
predictTm = function(trTmMut, predictRef = mSigsitedb, contexts.needed = FALSE, tri.counts.method = 'default'){
    sampleList = rownames(trTmMut)
    deconSigs = list()
    
    for (sample_id in sampleList) {
      deconSigs[[sample_id]] = deconstructSigs::whichSignatures(tumor.ref = trTmMut,
                                                                sample.id = sample_id,
                                                                contexts.needed = contexts.needed,
                                                                tri.counts.method = tri.counts.method,
                                                                signatures.ref = signatures.cosmic)
    }
    
    siteList = names(table(predictRef$tmSite))
    predicScores = as.data.frame(setNames(replicate(nrow(predictRef),numeric(0), simplify = F), siteList))
    
    for (sample_id in sampleList) {
      samplePriScore = sapply(siteList, FUN = function(site){
        sampleSig = deconSigs[[sample_id]]
        sampleSig = as.numeric(unname(sampleSig$weights))
        siteSig = unname(unlist(predictRef[(predictRef$tmSite == site), ]))[c(2:31)]  
        sampleScore = computeSim(sampleSig, siteSig) 
      })
      samplePriScore = as.list(samplePriScore)
      predicScores = rbind(predicScores, samplePriScore, stringsAsFactors = FALSE)
    }
    
    row.names(predicScores) = sampleList
    predicScores
}


#' Title Plot the similarity
#'
#' @param contexts.needed FALSE if tumor.file is a context file, TRUE if it is only mutation counts.
#' @param tri.counts.method Set to either:
#' \itemize{
#'  \item 'default' -- no further normalization \item 'exome' -- normalized by
#'   number of times each trinucleotide context is observed in the exome \item
#'   'genome' -- normalized by number of times each trinucleotide context is
#'   observed in the genome \item 'exome2genome' -- multiplied by a ratio of that
#'   trinucleotide's occurence in the genome to the trinucleotide's occurence in
#'   the exome \item 'genome2exome' -- multiplied by a ratio of that
#'   trinucleotide's occurence in the exome to the trinucleotide's occurence in
#'   the genome \item data frame containing user defined scaling factor -- count
#'   data for each trinucleotide context is multiplied by the corresponding value
#'   given in the data frame }
#' @param trTmMut Either a data frame or location of input text file, where rows are samples, columns are trinucleotide contexts.
#' Use data(sample.mut.trans) to see the example data.
#' @param predictRef Location of the mutation file that is to be converted or name of data frame in environment.
#' "tmMutRef" is the standard signature value for one site.
#' Default value is "mSigsitedb" that convert from "mSignatureDB" database.
#' USe data("mSigsitedb") to see the example dasta.
#' @param sampleId Name of sample.
#'
#' @return Plot a histogram of similarity scores for one sample with every tumor.
#' @export
#'
#' @examples testPriPlot = plotPredict(trTmMut = sample.mut.trans, contexts.needed = TRUE, sampleId = 1)
plotPredict = function(trTmMut, predictRef = mSigsitedb, sampleId, contexts.needed = FALSE, tri.counts.method = 'default'){
  mSigsitedb$tmSite = gsub(" ","", mSigsitedb$tmSite)
  sampleSig = deconstructSigs::whichSignatures(tumor.ref = trTmMut,
                                               sample.id = sampleId,
                                               contexts.needed = contexts.needed,
                                               tri.counts.method = tri.counts.method,
                                               signatures.ref = signatures.cosmic)
  siteList = names(table(predictRef$tmSite))
  samplePriScore = sapply(siteList, FUN = function(site){
    sampleSig = as.numeric(unname(sampleSig$weights))
    siteSig = unname(unlist(predictRef[(predictRef$tmSite == site), ]))[c(2:31)]  
    sampleScore = computeSim(sampleSig, siteSig) 
  })
  sortSite = siteList[order(samplePriScore)]
  layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2), nrow = 3))
  par(mar=c(3,8,3,0),las = '2') 
  predictbar = barplot(samplePriScore, col=RColorBrewer::brewer.pal(10, "Set3"), names.arg = siteList, horiz = TRUE, xlim = c(0,1))
  text(cex=1, y = predictbar, x=unname(unlist(samplePriScore))+0.1, labels = round(unname(unlist(samplePriScore)),4))
  par(mar=c(1,0,5,0),las = '2')
  plot(x=1:3,y=1:3,type="n",bty="n",ann=F,axes=F)
  legend("topright", legend = sortSite[1:3], col = RColorBrewer::brewer.pal(3, "Set3"), pch = 15,cex = 1.2,pt.cex = 4 ,xpd=TRUE,bty = "n", title = "TOP 3 SITES")
}

#' Plot the signature ratio
#'
#' @param tmMutRef Location of the mutation file that is to be converted or name of data frame in environment.
#' "tmMutRef" is the standard signature value for one site and one country.
#' Default value is "mSigdb" that convert from "mSignatureDB" database.
#' USe data("mSigdb") to see the example dasta.
#' @param sampleId Name of sample.
#' @param tmSite Site of sample.
#' @param country Country of sample.
#' @param contexts.needed FALSE if tumor.file is a context file, TRUE if it is only mutation counts.
#' @param trTmMut Either a data frame or location of input text file, where rows are samples, columns are trinucleotide contexts.
#' Use data(randomly.generated.tumors) to see the example data.
#' @param tri.counts.method Set to either:
#' \itemize{
#'  \item 'default' -- no further normalization \item 'exome' -- normalized by
#'   number of times each trinucleotide context is observed in the exome \item
#'   'genome' -- normalized by number of times each trinucleotide context is
#'   observed in the genome \item 'exome2genome' -- multiplied by a ratio of that
#'   trinucleotide's occurence in the genome to the trinucleotide's occurence in
#'   the exome \item 'genome2exome' -- multiplied by a ratio of that
#'   trinucleotide's occurence in the exome to the trinucleotide's occurence in
#'   the genome \item data frame containing user defined scaling factor -- count
#'   data for each trinucleotide context is multiplied by the corresponding value
#'   given in the data frame }
#'
#' @return Plot a histogram to compare sample and ref.
#' @export
#'
#' @examples testPlot = plotSig(trTmMut = sample.mut.trans, tmSite = "bone", country = "UK", sampleId = "1", contexts.needed = TRUE)
plotSig = function(trTmMut, tmSite, country, tmMutRef = mSigdb, sampleId, contexts.needed = FALSE, tri.counts.method = 'default'){
  mSigsitedb$tmSite = gsub(" ","", mSigsitedb$tmSite)
  siteList = names(table(tmMutRef$tmSite))
  
  if (!(tmSite%in%siteList)) {
    stop("The tmSite is not in tmMutref")
  }
  
  tmRef = tmMutRef[tmMutRef$tmSite==tmSite, ]
  countryList = as.character(tmRef$Country)
  
  if (!(country%in%countryList)) {
    stop(paste("The country should be",countryList))
  }
  
  sampleSig = deconstructSigs::whichSignatures(tumor.ref = trTmMut,
                              sample.id = sampleId,
                              contexts.needed = contexts.needed,
                              tri.counts.method = tri.counts.method,
                              signatures.ref = signatures.cosmic)
  
  xlab = paste("signature",c(1:30),sep = "")
  campareDf =  rbind.data.frame(sampleSig$weights, unname(unlist(tmMutRef[(tmMutRef$tmSite==tmSite & tmMutRef$Country == country),]))[c(3:32)], make.row.names =  FALSE)
  boxColor = RColorBrewer::brewer.pal(3, "Set3")
  layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2), nrow = 3))
  par(mar=c(10,3,3,0),las = '2')
  barplot(as.matrix(campareDf), col = boxColor[c(1,2)], ylim = c(0,1), cex.names = 1.5)
  par(mar=c(1,0,5,0),las = '2')
  plot(x=1:3,y=1:3,type="n",bty="n",ann=F,axes=F)
  legend("topright" , legend = c("Sample","Reference"), col = boxColor[c(1,2)], pch = 15, cex = 1.2,pt.cex = 4 ,xpd=TRUE,bty = "n")
}

