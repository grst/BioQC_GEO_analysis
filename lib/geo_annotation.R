#' Functions related to extracting information from GEO studies.
#' 
#' Information retrieved from the GEO needs to be preprocessed
#' and santized before we send it into our database. 

geo_annotation.tissueMap = read.csv("lib/res/map_tissue_annotation.csv", strip.white=TRUE, stringsAsFactors=FALSE)


#' Retrieve Gene Symbols for BioQC
#' 
#' @param eset
#' @return eset with addtional BioqcGeneSymbol column. 
attachGeneSymbols = function(eset) {
  annotation = annotateProbesets(fData(eset)$ID, orthologue=TRUE)
  gene.symbol.col = annotation[,"GeneSymbol",drop=FALSE]
  colnames(gene.symbol.col) = c("BioqcGeneSymbol")
  fData(eset) = cbind(fData(eset), gene.symbol.col)
  return(eset)
}

#' Extract the tissue name from sample annotation 
#'
#' In GEO Series, tissue information can be found
#' in one of the characteristics_ch1 column.  
#' 
#' @param characteristics_ch1 column content
extractTissue = function(characteristics) {
  characteristics = unlist(strsplit(characteristics, split=";"))
  for(char in characteristics) {
    char = trim(char)
    if(substring(char, 0, 7) == "tissue:"){
      return(trim(substring(char, 8)))
    }
  }
}

#' Map the often very specific tissue annotation from GEO
#' to one tissue type, e.g. liver. 
#'
#' @param tissue.orig tissue annotation as found in the GEO
#'    Series and retrieved by \code{\link{extractTissue}}
sanitizeTissue = function(tissue.orig) {
  if(tissue.orig %in% geo_annotation.tissueMap$annotation) {
    geo_annotation.tissueMap[geo_annotation.tissueMap[,'annotation'] == tissue.orig,'tissue']    
  } else {
    return("other")
  }
}

geoIdFromPath = function(path) {
  pat = "((GSE|GDS)\\d+)(_\\d+)?.Rdata"  
  id = sub(pat, "\\1", basename(path))
  return(id)
}


