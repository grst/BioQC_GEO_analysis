#!Rscript

########
# USAGE:
#   Rscript run_bioqc.R <eset.Rdata> <output_file.tab>
# where eset.Rdata contains a biobase ExpressionSet with variable
# name eset. 
#
# The script applies BioQC to the ExpressionSet and stores the
# raw p-values as a data table.
########

args = commandArgs(trailingOnly = TRUE)
esetFile = args[1]
outputFile = args[2]

#esetFile = "/homebasel/biocomp/sturmg/projects/GEO_BioQC/GDS_GPL570/GDS4074.Rdata"
#outputFile = "/homebasel/biocomp/sturmg/projects/GEO_BioQC/BioQC_GEO_analysis/plots/GDS4074"

library(BioQC)

gmtFile = system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
gmt <- readGmt(gmtFile)
load(esetFile)

#run BioQC
bioqcRes = wmwTest(eset, gmt, valType="p.greater", col="Gene symbol")

#save result to table
write.table(bioqcRes, file=outputFile)
