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
  library(ribiosAnnotation)
  ### We focus on GeneExpression Data atm
  #   # Unfortunately, the anntation of Gene symbols is not consistent within 
  #   # GEO Series. These names should cover most of it, though
  #   gene.symbols = c("Gene symbol", 
  #                    "Gene Symbol",
  #                    "GENE_SYMBOL",
  #                    "Gene_symbol", 
  #                    "GeneSymbol",
  #                    "geneSymbol",
  #                    "G_Symbol",
  #                    "Symbol",
  #                    "SYMBOL",
  #                    "Gene.Symbol",
  #                    "GENE")
  annotation = annotateProbesets(fData(eset)$ID, orthologue=TRUE)
  gene.symbol.col = annotation[,"GeneSymbol",drop=FALSE]
  colnames(gene.symbol.col) = c("BioqcGeneSymbol")
  fData(eset) = cbind(fData(eset), gene.symbol.col)
  return(eset)
}

#' Extract the tissue name from sample annotation 
#'
#' In GEO Series, tissue information can be found
#' in one of the characteristics_ch?? columns starting
#' with "tissue: ..." 
#' We simply search all columns for this pattern. 
#' 
#' @param pdata.row row of pData(eset)
extractTissue = function(pdata.row) {
  if("tissue" %in% colnames(pdata.row)) {
    return(pdata.row$tissue)
  } else {
    pat = "tissue: (.*)"
    return(tolower(sub(pat, "\\1", pdata.row[grepl(pat, pdata.row)])[1]))
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


