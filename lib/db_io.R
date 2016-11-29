###############################
# Library to read and write data from and to the database. 
###############################

library(stringr)

#' Escape a string for SQL
#' 
#' @note this is not to be understood as being safe in any way. Use only when prepared 
#' statements are note an option
#' @param text the string to escape.
db_escape = function(text) {
  return(gsub("'", "''", text))
}

#' Read a GMT file and write it to the signatures database
#' 
#' @param gmt_file path to gmt file
gmt2db = function(gmt_file) {
  filename = basename(gmt_file)
  gmt <- readGmt(gmt_file)
  gmt_list = lapply(gmt, function(line) {
    return(list(id=NA, name=line$name, source=filename, desc=line$desc, genes=paste(line$genes, collapse=',')))
  })
  gmt_table = do.call(rbind.data.frame, gmt_list)
  dbAppendDf("BIOQC_SIGNATURES", gmt_table)
}

#' Download all signatures from the database and combine them in one gmt. 
#' 
#' Uses the database id as signature identifier. 
#' @param output_file
db2gmt = function(output_file) {
  signatuers = dbGetQuery(mydb, "select * from bioqc_signatures
                           order by source, name")
  if(file.exists(output_file)) {
    # need to clear the file, as we are appending later. 
    file.remove(output_file)
  }
  for (i in 1:nrow(signatures)) {
    row = signatures[i, ]
    name = as.character(row$ID)
    desc = str_c(row$SOURCE, row$NAME, sep=":")
    genes = str_c(unlist(str_split(row$GENE_SYMBOLS, ",")), collapse="\t") 
    cat(paste(str_c(name, desc, genes, sep="\t"), "\n"), file=output_file, append=TRUE)
  }
}

#' Read a BioQC result file (pvalue-matrix) and write it to 
#' the BIOQC_RES table. 
#' 
#' Signature names need to be the numeric ids from the signatures
#' table. Ideally, create your gmt file with db2gmt
bioqc2db = function(bioqc_res_file) {
  
}