# Validating Tissue Signatures

The authors of *BioQC* have taken three independent approaches to show that their signatures are valid and biologically meaningful. However, they did not break down the predictive performance (*i.e.* is the signature able to identify its tissue) of each signature with quantitative performance measures. 

To address this, we independently derived signatures on the GTEx dataset using *gini-index* and performed both a 10-fold cross validation on the same dataset and a cross-species, cross-platform validation on a different dataset. To this end, we created the python package [pygenesig](https://github.com/grst/pygenesig), a framework to create and validate signatures. 


## Data
 * GTEx
 * GNF Mouse GeneAtlas V3 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10246
    - platform Affymetrix Mouse Genome 430 2.0 Array (GPL1261)
