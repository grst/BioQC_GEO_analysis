#!Rscript

########
# USAGE:
#   Rscript run_bioqc.R <esetName.Rdata> <outputPath>
# where esetName.Rdata contains a biobase ExpressionSet with variable
# name eset. outputPath will be appended with <esetName_bioqc_res.tab>
#
# The script applies BioQC to the ExpressionSet and stores the
# raw p-values as a data table.
########

library(tools)

args = commandArgs(trailingOnly = TRUE)
esetFile = args[1]
outputPath = args[2]
geo.id = file_path_sans_ext(basename(esetFile))
print(geo.id)

outFile = file.path(outPath, paste(geo.id, "_bioqc_res.tab", sep=""))

# get the first colname of the list that is in fData
col.name = gene.symbols[which(gene.symbols %in% colnames(fData(eset)))[1]]

if(is.na(col.name)) {
  print(sprintf("%s: Study does not provide Gene Symbols.", geo.id))
} else{
  library(BioQC)
  
  # do BioQC analysis and save pValues to table 
  gmtFile = system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
  gmt <- readGmt(gmtFile)
  load(esetFile)
  
  #run BioQC
  bioqcRes = wmwTest(eset, gmt, valType="p.greater", col="BioqcGeneSymbol")
  
  #save result to table
  write.table(bioqcRes, file=outFile)
}


