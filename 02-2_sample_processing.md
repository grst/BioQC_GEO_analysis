## Processing Samples with BioQC



The following processes are ressource intensive, therefore we execute them on high performance cluster (HPC). We use [chunksub](https://github.com/grst/chunksub) to distribute the [list of sample ids](#sample-list) to the workers. This involves four major steps which are also documented in the project's [Makefile](https://github.com/grst/BioQC_GEO_analysis/blob/master/Makefile). 

1. We download the studies with [GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html) and store them as R [ExpressionSet](https://bioconductor.org/packages/devel/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) using [geo_to_eset.R](https://github.com/grst/BioQC_GEO_analysis/blob/master/scripts/geo_to_eset.R). For some packages, the download is not successful. 
2. We annotated human orthologous genes for all studies using [ribiosAnnotation](https://github.com/Accio/ribios) in [annotate_eset.R](https://github.com/grst/BioQC_GEO_analysis/blob/master/scripts/annotate_eset.R). This is necessary as the tissue signatures are built on a human dataset. The annotation failes for species which are not in the `ribios` database. 
3. We run *BioQC* on these studies use [run_bioqc.R](https://github.com/grst/BioQC_GEO_analysis/blob/master/scripts/run_bioqc.R). 
4. Finally, we prefilter BioQC's results for having a p-value of at most 0.05 and import them into the database. 


After failures, we are left with the following number of studies/samples:



Although we lost a considerable amount of studies, we only lost a few samples. Apparently, the respective studies were relatively small. 
<!-- i also found empty studies -->

## Excluding multi channel microarrays

Multi channel microarrays date back to the early age of gene expression studies. They don't provide absolute gene expression data and are not meaningful outside their experimental context. We therefore exclude these experiments


## Organism selection. 
We were interested in the organism distribution.


Results suggest that it makes sense to limit the samples to the three main organisms: *H. sapiens*, *M. musculus*, *R. norvegicus*. This makes also sense as these species are closesly related and therefore the signatures are more likely to translate within these species. 
We are left with the following amount of samples: 



We create a materialized view `BIOQC_RES_FIL` with bioqc results that are filted according 
to these criteria in `db/views/bioqc_res_fil.sql`. Make sure the results are identical: 



## Tissue selection
Number of samples per tissue:

