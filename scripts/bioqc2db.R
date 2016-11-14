#!/bin/env Rscript

###################
# Script reads BioQC result files and writes them into the 
# postgreSQL database. 
# 
# USAGE:
#   bioqc2db.R <chunkfile>
# where <chunkfile> contains the path of one bioqc_res table file per line. 
# 
# the table files are the result of write.table(wmwTest(...))
##################

stopifnot(suppressPackageStartupMessages(require(data.table)))
source("lib/db.R")
source("lib/geo_annotation.R")

args = commandArgs(trailingOnly=TRUE)
chunkFile = args[1]
bqcFiles = readLines(chunkFile)

for(bqcFile in bqcFiles) {
    gse = geoIdFromPath(bqcFile)
    bioqcRes = data.table(read.table(bqcFile), keep.rownames=TRUE)
    print(sprintf("Processing %s with %d rows and %d cols", bqcFile, nrow(bioqcRes), ncol(bioqcRes)))
    res.molten = melt(bioqcRes, id.vars="rn")
    setcolorder(res.molten, c(2,1,3))
    table = cbind(rep(gse, nrow(res.molten)), res.molten)
    dbAppendDf("bioqc_res", table)
}
