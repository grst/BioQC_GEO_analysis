#' Functions related to extracting information from GEO studies.
#' 
#' Information retrieved from the GEO needs to be preprocessed
#' and santized before we send it into our database. 

library(assertthat)
library(ribiosAnnotation)
library(ribiosUtils)
geo_annotation.tissueMap = read.csv("lib/res/map_tissue_annotation.csv", strip.white=TRUE, stringsAsFactors=FALSE)


#' Retrieve Gene Symbols for BioQC
#' 
#' @param eset
#' @param platform.id
#' @return eset with addtional BioqcGeneSymbol column. 
attachGeneSymbols = function(eset, platform.id=NULL) {
  annotation.package = dbGetQuery(mydb, "select bioc_package from gpl where gpl = ?", platform.id)
  if(nrow(annotation.package) > 0 && !is.na(annotation.package[1,1])) {
    assert_that(nrow(annotation.package) == 1)
    package.name = sprintf("%s.db", annotation.package[1,1])
    require(package.name, character.only=TRUE)
    fdata = fData(eset)
    gene.ids = select(get(package.name), keys=as.character(fdata$ID), columns=c("ENTREZID"), keytype="PROBEID")
    # returns a 1:many mapping. Use matchColumn to resolve that
    gene.ids.matched = matchColumn(fdata$ID, gene.ids, "PROBEID", multi=FALSE)
    ortholog.res = annotateHumanOrthologs(gene.ids.matched$ENTREZID)
    gene.symbols.orth = matchColumn(gene.ids.matched$ENTREZID, ortholog.res, "OrigGeneID", multi=FALSE)
    # save back to fData
    fdata = cbind(fdata, data.frame(BioqcGeneSymbol=gene.symbols.orth$GeneSymbol))
    fData(eset) = fdata
    return(eset)
  } else {
    # return all NA's if cannot annotate. 
    return cbind(fData(eset), rep(NA, ncol(fData(eset))))
  } 
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
  return("no_tissue_annotated")
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


