#' Functions related to extracting information from GEO studies.
#' 
#' Information retrieved from the GEO needs to be preprocessed
#' and santized before we send it into our database. 


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
#'
extractTissue = function(pdata.row) {
  if("tissue" %in% colnames(pdata.row)) {
    return(pdata.row$tissue)
  } else {
    pat = "tissue: (.*)"
    return(tolower(sub(pat, "\\1", pdata.row[grepl(pat, pdata.row)])[1]))
  } 
}

#' 
#'
#'
sanitizeTissue = function(tissue.orig) {
  
}


