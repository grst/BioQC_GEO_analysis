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
geo_id = file_path_sans_ext(basename(esetFile))
print(geo_id)

outFile = file.path(outPath, paste(geoID, "_bioqc_res.tab", sep=""))

#esetFile = "/homebasel/biocomp/sturmg/projects/GEO_BioQC/GDS_GPL570/GDS4074.Rdata"
#outFile = "/homebasel/biocomp/sturmg/projects/GEO_BioQC/BioQC_GEO_analysis/plots/GDS4074"

library(BioQC)

gmtFile = system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
gmt <- readGmt(gmtFile)
load(esetFile)

#run BioQC
bioqcRes = wmwTest(eset, gmt, valType="p.greater", col="Gene symbol")

#save result to table
write.table(bioqcRes, file=outFile)
