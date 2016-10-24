###################
# Library for BioQC GEO analysis. 
# 
# Collection of functions that are used in 
# multiple analyses. 
##################


#' Normalize a vector between 0 and 1
#'
#' @param x A vector.
norm01 = function(x){
  (x-min(x))/(max(x)-min(x))
}


#' Collapse multiple signatures that belong to the same tissue into one.
#' 
#' BioQC comes with multiple signatures per tissue e.g. Liver, Liver_fetal
#' Liver_NGS, ...
#' The annotation in the GEO is more coarse than that, only stating 'liver'. 
#' This function therefore collapses multiple signatures into one by taking the
#' maximum score of the selected signatures. 
#' 
#' @param table
#' @param sig.list list of the signature names corresponding to the rownames of the table
#' @param new.name new row name for the collapsed rows. 
#' @param method function to be applied to each column. Defaults to max.  
collapseSignatures = function(table, sig.list, new.name, method=max) {
  row.collapsed = apply(table[sig.list, ], 2, method)
  row.collapsed.df = data.frame(t(unlist(row.collapsed)))
  rownames(row.collapsed.df) = c(new.name)
  table.collapsed = rbind(table, row.collapsed.df)
  return(table.collapsed[-which(rownames(table.collapsed) %in% sig.list), ,drop=FALSE])
}
