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

tissues = dbGetQuery(mydb, "
  select distinct bnt.tissue
  from bioqc_gsm bg
  join bioqc_normalize_tissues bnt on bnt.tissue_orig = lower(bg.tissue_orig)
  join bioqc_tissue_set bts on bts.tissue = bnt.tissue
  where bts.tissue_set = 'gtex_all'")

getTissueSamples = function(tissue) {
  # get all samples belonging to one tissue and filter for siginifant 
  # signatures in one sql query! 
  query = "
  with gmt_signatures as (
    select * from bioqc_signatures 
    where source = 'gtex_ngs_0.85_5.gmt'
  )
  select /*+ parallel(16) */ br.gsm
         , br.tissue
        -- , ur.signature
         , gs.name as SIGNATURE
         , br.pvalue
    from bioqc_res_fil br
    join gmt_signatures gs on gs.id = br.signature
    where tissue = ?
    order by gsm
  "
  data = dbGetQuery(mydb, query, tissue)
  data = data.frame(data, pvalue.log=absLog10p(as.numeric(data[,"PVALUE"])))
  data$GSM = as.character(data$GSM)
  return(data)
}

for(tissue in tissues$TISSUE) {
  data = data.table(getTissueSamples(tissue))
  data$GSM = factor(data$GSM, levels=unique(data$GSM))
  data$SIGNATURE = factor(data$SIGNATURE, levels=sort(unique(data$SIGNATURE), decreasing = TRUE))
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
