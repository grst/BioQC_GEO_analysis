library(RColorBrewer)
library(ggplot2)
library(reshape2)

###################
# Library for BioQC GEO analysis. 
# 
# Collection of functions that are used in 
# multiple analyses. 
##################


#' Normalize a vector between 0 and 1
#'
#' @param x A vector.
norm01 = function(x){
       (x-min(x))/(max(x)-min(x))
}


#' Perform and Plot a PCA of an ExpressionSet
#' 
#' @param eset An ExpressionSet
esetPca = function(eset, title) {
    pca = prcomp(t(exprs(eset)))
    expVar <- function(pcaRes, n) {vars <- pcaRes$sdev^2; (vars/sum(vars))[n]}
    biplot(pca, col=c("#335555dd", "transparent"), cex=1.15,
            xlab=sprintf("Principal component 1 (%1.2f%%)", expVar(pca,1)*100),
            ylab=sprintf("Principal component 1 (%1.2f%%)", expVar(pca,2)*100),
            main=title)
}


#' Create a heatmap from a Matrix
#'
#' The matrix contains samples as columns and
#' tissue signatures as rows. 
#'
#' @param bioqc_res A Matrix containing the (transformed) p-values. 
bioqcHeatmap = function(bioqc_res, title) {
    hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')  
    mat.melted = melt(bioqc_res)
    ggplot(data=mat.melted, aes(x=Var2, y=Var1, fill=value)) + geom_tile() +
            coord_equal() +
            scale_fill_gradientn(colours = hm.palette(100)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            ggtitle(title)
}


#' Collapse multiple signatures that belong to the same tissue into one. 
collapseSignatures = function(table, sig_list, method=max) {
  
}
