#!/bin/env Rscript 

#############
# USAGE:
#   run_bioqc.R <outputDir> <chunkFile> 
# where chunkFile is a file containing paths to Rdata objects storing
# ExpressionSets, one file per line. The ExpressionSet in the Rdata object
# must be named 'eset'. 
#
# The script runs BioQC on each expression set and stores the raw 
# p-values as data tables. 
#############


library(tools)
library(Biobase)
library(BioQC)
library(ribiosAnnotation)
source("lib/lib.R")
source("lib/geo_annotation.R")

# options(error = quote({
#   dump.frames("ribios.dump", to.file = TRUE)
#   quit(save = "no", status = 1L)
# }))

args = commandArgs(trailingOnly=TRUE)

chunkFile = args[2]
outDir = args[1]
outPath = file.path(outDir, "%s_bioqc_res.tab")
esetFiles = readLines(chunkFile)

gmtFile = system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
gmt <- readGmt(gmtFile)

runFile = function(esetFile) {
  load(esetFile)
  eset = attachGeneSymbols(eset)
  #run BioQC
  bioqcRes = wmwTest(eset, gmt, valType="p.greater", col="BioqcGeneSymbol")
  return(bioqcRes)
}

for (esetFile in esetFiles) {
  print(sprintf("Working on %s", esetFile))
  tryCatch ({
    res = runFile(esetFile)
    outFile = sprintf(outPath, tools::file_path_sans_ext(basename(esetFile)))
    print(sprintf("Writing to %s", outFile))
    write.table(res, file=outFile)
  }, 
  error=function(cond) {
    print(sprintf("%s failed: ", esetFile))
    message(cond)
  })
}


