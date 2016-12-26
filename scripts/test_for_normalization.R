#!/bin/env Rscript 

#############
# USAGE:
#   test_for_normalization.R <outputDir> <chunkFile> 
# where chunkFile is a file containing paths to Rdata objects storing
# ExpressionSets, one file per line. The ExpressionSet in the Rdata object
# must be named 'eset_res'. 
#
# If gmtFile is omitted, the default BioQC gmt file will be used. 
#
# The script runs BioQC on each expression set and stores the raw 
# p-values as data tables. 
#############


stopifnot(suppressPackageStartupMessages(require(tools)))
stopifnot(suppressPackageStartupMessages(require(readr)))
stopifnot(suppressPackageStartupMessages(require(Biobase)))
source("lib/lib.R")
source("lib/db_io.R")


args = commandArgs(trailingOnly=TRUE)

chunkFile = args[2]
outDir = args[1]
outPath = file.path(outDir, "%s_stats.txt")
esetFiles = readLines(chunkFile)

runFile = function(esetFile) {
  load(esetFile)
  s = summary(apply(exprs(eset_res), 1, mean))
  out_file = sprintf(outPath, tools::file_path_sans_ext(basename(esetFile)))
  capture.output(s, file=out_file)
  print(sprintf("%s written to %s.", esetFile, out_file))
}

for (esetFile in esetFiles) {
  print(sprintf("%s started.", esetFile))
  tryCatch ({
    runFile(esetFile)
  }, 
  error=function(cond) {
    print(sprintf("%s failed: ", esetFile))
    print(cond)
  })
}


