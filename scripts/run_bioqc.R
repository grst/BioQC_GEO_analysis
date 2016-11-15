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


stopifnot(suppressPackageStartupMessages(require(tools)))
stopifnot(suppressPackageStartupMessages(require(Biobase)))
stopifnot(suppressPackageStartupMessages(require(BioQC)))
stopifnot(suppressPackageStartupMessages(require(ribiosAnnotation)))
stopifnot(suppressPackageStartupMessages(require(assertthat)))
source("lib/lib.R")
source("lib/geo_annotation.R")
# source("lib/db.R")

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
  # assert_that(length(levels(as.factor(pData(eset)$platform_id))) == 1)
  # platform.id = as.character(pData(eset)[1, 'platform_id'])
  
  eset = attachOrthologousSymbols(eset)
  
  # run BioQC
  if(!all(is.na(fData(eset)[["BioqcGeneSymbol"]]))) {
      bioqcRes = wmwTest(eset, gmt, valType="p.greater", col="BioqcGeneSymbol")
      return(bioqcRes)
  } else {
      stop("Gene Symbols could not be annotated.")
  }
}

for (esetFile in esetFiles) {
  print(sprintf("%s started.", esetFile))
  tryCatch ({
    res = runFile(esetFile)
    outFile = sprintf(outPath, tools::file_path_sans_ext(basename(esetFile)))
    write.table(res, file=outFile)
    print(sprintf("%s written to: %s", esetFile, outFile))
  }, 
  error=function(cond) {
    print(sprintf("%s failed: ", esetFile))
    print(cond)
  })
}


