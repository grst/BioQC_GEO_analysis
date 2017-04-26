# Sample Selection and processing {#sample-selection}



Here, we document the sample selection process before running *BioQC*. 

## Required annotation
A sample is *usable* for this study, if 

* the gene symbos are annotated (requirement to run BioQC)
* the tissue of origin is annotated (requirment to draw conclusions about contamination)

We consider two approaches for annotating gene symbols: 

* Using the Bioconductor [AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html) package. The GEOmetadb provides a mapping of the GPL identifier to these packages. 
* Using the `annotGPL=TRUE` option of [GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html)'s `getGEO`. This requires an annotation file being available for the respective platform. We retrieved a list of the available annoataion files in an earlier step. 

We compare the two methods with respect to the amount of usable samples that we can get. 

First we need to run the sql-script [annotation_stats](https://github.com/grst/BioQC_GEO_analysis/blob/master/db/views/annotation_stats.sql) which will create views summarizing the amount of usable samples. 
From these views, we can calculate a Venn diagram: 

 
<img src="02_select_and_get_samples_files/figure-html/sample_filtering-1.png" width="768" style="display:block; margin: auto" style="display: block; margin: auto;" />

The `getGEO` method appears to be the more powerful method. Ideal would be a combination of the two, however, for the sake of simplicity, we stick to `getGEO`, loosing 499 studies (35602 samples). 

This leaves us with the following filtering result: 

comment                           GSM     GSE
---------------------------  --------  ------
total                         1945417   73719
tissue annotated               760798   24267
annotation file available      768346   31579
tissue and annotation file     275206    9632



### Export list of samples {#sample-list}
We store the respective gse identifiers in `results/gse_lists/gse_tissue_annotation.txt`: 

```r
sqlGse = "
select distinct gse
from (
  select *
  from bioqc_studies_has_tissue
  intersect
    select *
    from bioqc_studies_has_annot) u"
gse = dbGetQuery(mydb, sqlGse)
writeLines(gse$GSE, file(gse.file))
```
