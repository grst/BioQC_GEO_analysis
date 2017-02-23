#!/bin/env Rscript

##########################################################
# Generate heatmaps for all samples, grouped by tissue. 
# 
# Saves output to results/heatmaps_db. 
##########################################################

source("lib/knitr.R")
source("lib/plots.R")
source("lib/db.R")
library(data.table)
library(BioQC)
library(scales)

prepend_control = Vectorize(function(str) {
  return(str_c("0", str, sep = "_"))
})

tissues = dbGetQuery(mydb, "
  select distinct tgroup
  from bioqc_tissue_set 
  where tissue_set = 'gtex_all'")

getTissueSamples = function(tissue) {
  # get all samples belonging to one tissue and filter for siginifant 
  # signatures in one sql query! 
  query = "
  select /*+ parallel(16) */ bsst.gsm
         , bsst.tgroup
         , found_sig_name as SIGNATURE
         , found_sig_pvalue as PVALUE
    from bioqc_selected_samples_tset bsst
    join bioqc_res_tset brt 
      on brt.gsm = bsst.gsm
      and brt.tissue_set = bsst.tissue_set
    where bsst.tgroup = ?
    and bsst.tissue_set = 'gtex_all'
    order by gsm
  "
  data = data.table(dbGetQuery(mydb, query, tissue))
  data[,pvalue.log:=absLog10p(as.numeric(PVALUE))]
  data[,GSM:=as.character(GSM)]
  return(data)
}

getReferenceSamples = function(tissue) { 
  query = "
  select /*+ parallel(16) */ bsst.gsm
                           , bsst.tgroup
                           , bs.name as SIGNATURE
                           , br.pvalue as PVALUE
  from bioqc_selected_samples_tset bsst
  join bioqc_res br
    on br.gsm = bsst.gsm
  join bioqc_signatures bs
    on bs.id = br.signature
  where bsst.tgroup = ?
  and bsst.tissue_set = 'gtex_all'
  and bs.source = 'baseline_signatures.gmt'
  and bs.name in ('random_10_0', 'random_10_1', 'random_100_0', 'random_100_1', 'awesome_housekeepers', 'enzyme_goslim')
  "
  data = data.table(dbGetQuery(mydb, query, tissue))
  data[,pvalue.log:=absLog10p(as.numeric(PVALUE))]
  data[,GSM:=as.character(GSM)]
  data[,SIGNATURE:=prepend_control(SIGNATURE)]
  return(data)  
}

for(tissue in tissues$TGROUP) {
  data = getTissueSamples(tissue)
  data = rbind(data, getReferenceSamples(tissue))
  data[,GSM:=factor(GSM, levels=unique(GSM))]
  data[,SIGNATURE:=factor(SIGNATURE, levels=sort(unique(SIGNATURE), decreasing = TRUE))]
  print(tissue)
  hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')  
  sampids = levels(data$GSM)
 
  pdf(file=sprintf("results/heatmaps_db/%s.pdf", tissue),
      width=min(nrow(data)*.3 + 5, 30),
      height=length(levels(data$SIGNATURE))*.33+2)
  for (i in seq(1, length(sampids), 70)) {
    print(ggplot(data=data[GSM %in% sampids[i:(i+70)], ], aes(x=GSM, y=SIGNATURE, fill=pvalue.log)) + 
            geom_tile() +
            coord_equal() +
            scale_fill_gradientn(colours = hm.palette(100), limits=c(0, 30), oob=squish) +
            scale_y_discrete(drop=FALSE) + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            ggtitle(tissue))
  }
  dev.off()
}
