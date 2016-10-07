#!Rscript

# GDS4882
#        GDS2526
#        GDS3315
#        GDS4551
#        GDS3481
#        GDS3290
#        GDS4074
#        GDS2453
#        GDS3836
#        GDS4252

args = commandArgs(trailingOnly = TRUE)
esetFile = args[1]
outputBasename = args[2]

#esetFile = "/homebasel/biocomp/sturmg/projects/GEO_BioQC/GDS_GPL570/GDS4074.Rdata"
#qoutputBasename = "/homebasel/biocomp/sturmg/projects/GEO_BioQC/BioQC_GEO_analysis/plots/GDS4074"

  
library(RColorBrewer)
library(BioQC)
library(ggplot2)
library(reshape2)
gmtFile = system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
gmt <- readGmt(gmtFile)
load(esetFile)

testEset = function(eset) {
  bioqcRes = wmwTest(eset, gmt, valType="p.greater", col="Gene symbol")
  # experimentData(eset)
  bioqcResFil <- filterPmat(bioqcRes, 1E-6)
  bioqcAbsLogRes <- absLog10p(bioqcResFil)
  write.table(bioqcRes, file=paste(outputBasename, "_bioqc_res.tab", sep=""))
  return(bioqcAbsLogRes)
}

bioqcAbsLogRes = testEset(eset)

# do heatmap
hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')  
mat.melted = melt(bioqcAbsLogRes)
print(sprintf("writing to %s", outputBasename))
pdf(paste(outputBasename, "_heatmap.pdf", sep=""))
ggplot(data=mat.melted, aes(x=Var2, y=Var1, fill=value)) + geom_tile() +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(basename(esetFile))
dev.off()


# do PCA
pca = prcomp(t(exprs(eset)))
expVar <- function(pcaRes, n) {vars <- pcaRes$sdev^2; (vars/sum(vars))[n]}
pdf(paste(outputBasename, "_pca.pdf", sep=""))
biplot(pca, col=c("#335555dd", "transparent"), cex=1.15,
       xlab=sprintf("Principal component 1 (%1.2f%%)", expVar(pca,1)*100),
       ylab=sprintf("Principal component 1 (%1.2f%%)", expVar(pca,2)*100),
       main=basename(esetFile))
dev.off()