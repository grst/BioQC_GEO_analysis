## Sample Processing with BioQC



The following processes are ressource intensive, therefore we execute them on a high performance cluster (HPC). We use [chunksub](https://github.com/grst/chunksub) to distribute the [list of sample ids](#sample-list) to the workers. This involves four major steps which are also documented in the project's [Makefile](https://github.com/grst/BioQC_GEO_analysis/blob/master/Makefile). 

1. We download the studies with [GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html) and store them as R [ExpressionSet](https://bioconductor.org/packages/devel/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) using the R script [geo_to_eset.R](https://github.com/grst/BioQC_GEO_analysis/blob/master/scripts/geo_to_eset.R). For some series, the download is not successful. 
2. We annotated human orthologous genes for all studies using [ribiosAnnotation](https://github.com/Accio/ribios) in [annotate_eset.R](https://github.com/grst/BioQC_GEO_analysis/blob/master/scripts/annotate_eset.R). This is necessary as the tissue signatures are built on a human dataset. The annotation failes for species which are not in the *ribios* database. 
3. We run *BioQC* on these studies use [run_bioqc.R](https://github.com/grst/BioQC_GEO_analysis/blob/master/scripts/run_bioqc.R). 
4. Finally, we prefilter BioQC's results for having a p-value < 0.05 and import them into the database. 

## Update database
we aggregate the bioqc results and import them manually in the DBS. 
we collate all samples on which we successfully ran BioQC. This is our Background. 


## Sample Post-selection

The failures during download and annotation reduce the number of samples available to our study.

<!-- since there are muliple gpl and gsm in a gse, there might've been some gsm 
inserted, that don't have a tissue annotated, albeit the gse was selected. -->


```r
sql_select = "
select /*+ USE_HASH(bs, bgg) parallel(16) */ count(distinct bg.gsm) as GSM
                                           , count(distinct bgg.gse) as GSE"
sql_from = "
from bioqc_bioqc_success bs
join bioqc_gsm bg
  on bg.gsm = bs.gsm
left outer join bioqc_gse_gsm bgg
  on bgg.gsm = bs.gsm 
"
sql_where = ""
res = dbGetQuery(mydb, str_c(sql_select, sql_from, sql_where, sep="\n"))
kable(res)
```

    GSM    GSE
-------  -----
 253714   8083


### Excluding multi-channel microarrays

Multi channel microarrays date back to the early age of gene expression studies. They don't provide absolute gene expression data and are not meaningful outside their experimental context. We therefore exclude these experiments:

```r
sql_select2 = sql_select
sql_from2 = sql_from
sql_where2 = str_c(sql_where, "where channel_count = 1", sep="\n")
res = dbGetQuery(mydb, str_c(sql_select2, sql_from2, sql_where2, sep="\n"))
kable(res)
```

    GSM    GSE
-------  -----
 235741   7561

### Exclude non-mapped tissues
We exclude samples that have a tissue annotated, but it is not mapped to a [normalized tissue](#tissue-normalization). 


```r
sql_select3 = sql_select2
sql_from3 = str_c(sql_from2, "
join bioqc_normalize_tissues bnt
  on bnt.tissue_orig = lower(bg.tissue_orig)", sep="\n")
sql_where3 = sql_where2
res = dbGetQuery(mydb, str_c(sql_select3, sql_from3, sql_where3, sep="\n"))
kable(res)
```

    GSM    GSE
-------  -----
 135877   3770

### Select organisms
We were interested in the organism distribution.

```r
sql_select4 = str_c(sql_select3, ", bg.organism_ch1", sep="\n")
sql_from4 = sql_from3
sql_where4 = sql_where3
res = dbGetQuery(mydb, str_c(sql_select4, sql_from4, sql_where4, "
group by organism_ch1
order by gsm desc", sep="\n"))
kable(res)
```



   GSM    GSE  ORGANISM_CH1                             
------  -----  -----------------------------------------
 65944   1201  Homo sapiens                             
 38128   2267  Mus musculus                             
 29909    278  Rattus norvegicus                        
  1082     24  Macaca mulatta                           
   259      7  Macaca fascicularis                      
   202      2  Mus musculus musculus x M. m. domesticus 
    80      2  Cercocebus atys                          
    57      1  Oryctolagus cuniculus                    
    36      1  Chlorocebus aethiops                     
    32      3  Mus musculus domesticus                  
    25      1  Pan troglodytes                          
    23      1  Papio cynocephalus                       
    19      1  Mus spretus                              
    18      1  Capra hircus                             
    16      1  Phodopus sungorus                        
    12      1  Mus musculus musculus x M. m. castaneus  
    12      1  Papio hamadryas                          
     8      1  Macaca nemestrina                        
     6      1  Mus musculus musculus                    
     6      1  Mus musculus castaneus                   
     3      1  Mus sp.                                  

Results suggest that it makes sense to limit the analysis to the three main organisms: *H. sapiens*, *M. musculus*, *R. norvegicus*. This makes also sense as these species are closesly related and therefore the signatures are more likely to translate within these species. 
We are left with the following amount of samples: 


```r
sql_select5 = sql_select3
sql_from5 = sql_from4
sql_where5 = str_c(sql_where4, "
  and organism_ch1 in ('Homo sapiens', 'Mus musculus', 'Rattus norvegicus')", sep="\n")
res = dbGetQuery(mydb, str_c(sql_select5, sql_from5, sql_where5, sep="\n"))
kable(res)
```

    GSM    GSE
-------  -----
 133981   3727

We create a materialized view `BIOQC_RES_FIL` with bioqc results that are filted according 
to these criteria in `db/views/bioqc_res_fil.sql`, except that it also contains the GTEx signatures. ``


### Select tissues

Of which tissues are enough samples available that we can make a meaningful statement about contamination? 

```r
sqlTissue = "
select /*+ parallel(16) */ tissue, count(distinct gsm) as samples from bioqc_res_fil
group by tissue
order by samples desc"
resTissue = dbGetQuery(mydb, sqlTissue)
kable(resTissue)
```



TISSUE               SAMPLES
------------------  --------
blood                  32548
liver                  26732
lung                    9553
kidney                  7507
bone marrow             6536
brain                   5719
heart                   4813
breast tumor            3990
spleen                  3506
adipose                 2857
skeletal muscle         2445
skin                    2271
hippocampus             2053
colon                   1914
hepatocyte              1870
tumor                   1627
white blood cells       1535
lymph node              1335
breast                  1205
cerebellum              1175
testis                  1128
pbmc                     936
frontal cortex           805
placenta                 766
pancreas                 746
retina                   711
pancreatic islets        645
thymus                   638
ovary                    549
prostate                 522
hypothalamus             468
mammary gland            464
jejunum                  446
prefrontal cortex        409
cortex                   328
embryo                   257
uterus                   239
cervix                   147
stomach                  139
salivary gland           128
bladder                  123
ventral midbrain          78
eye                       69
adrenal gland             67
monocyte                  59
neuroblastoma             27
plasma                     9
