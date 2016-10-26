#!/bin/env Rscript 

#############
# USAGE:
#   mk_sample_db.R <outputDir> <chunkFile> 
# where chunkFile is a file containing paths to Rdata objects storing
# ExpressionSets, one file per line. The ExpressionSet in the Rdata object
# must be named 'eset'. 
#
# The script extracts meta-information of Series and Samples and writes them
# to the Database
#############


library(tools)
library(Biobase)
source("lib/lib.R")
source("lib/geo_annotation.R")

args = commandArgs(trailingOnly=TRUE)

chunkFile = args[2]
outDir = args[1]
outPath = file.path(outDir, "%s_samples.tab")
esetFiles = readLines(chunkFile)

runFile = function(esetFile) {
  load(esetFile)
  series.id = geoIdFromPath(esetFile)
  tissue.orig = apply(pData(eset), 1, extractTissue)
  tissue = sapply(tissue.orig, sanitizeTissue)
  samples = data.frame(id=pData(eset)$geo_accession,
                        series=rep(series.id, nrow(pData(eset))),   
                        platform=pData(eset)$platform_id,
                        tissue=tissue,
                        organism=rep(NA, nrow(pData(eset))),
                        tissue_orig=tissue.orig)
  return(samples)
}

for (esetFile in esetFiles) {
  samples = runFile(esetFile)
  outFile = sprintf(outPath, tools::file_path_sans_ext(basename(esetFile)))
  print(sprintf("Writing to %s", outFile))
  write.table(samples, file=outFile)
}

