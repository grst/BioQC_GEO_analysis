# BioQC GEO Analysis

Systematically testing the GEO for *tissue heterogeneity*. 

*Manuscript in preparation*

## Introduction
Recently we created a software tool, BioQC, that detects tissue heterogeneity in gene expression data and shared it with the community of genome researchers via Bioconductor. The concept of tissue heterogeneity stems from our observations that gene expression data is often compromised by cells originating from other tissues than the target tissue of profiling. Tissue heterogeneity can be caused by physiological or pathological processes, such as immune cell infiltration. Alternatively, they can be caused by technical imperfection to separate complex or nearby tissues or even human errors. Failures in detecting tissue heterogeneity may have profound implications on data interpretation and reproducibility. 

As bioinformaticians working on drug discovery in the pharma industry, we are convinced that gene expression data available in publicly available databases such as NCBI Gene Expression Omnibus (GEO) or EBI ArrayExpress has great potential to catalyse new therapeutic agents. Disease signatures derived from disease models or patient biopsies, for instance, can be used to assess cellular models used for discovery and to guide compound selection. Molecular phenotypes of compounds, in another example, can be used to validate both efficacy and pre-clinical safety of compounds. Apparently all such applications depend critically on the quality of gene expression data. Several groups have scrutinised publicly available datasets and have identified deleterious factors of data quality such as batch effects, human error, and even data manipulation and faking. However, tissue heterogeneity has not been explicitly addressed so far and there is neither data nor knowledge about its prevalence. To fill this gap, we undertake a systematic investigation of publicly available gene expression datasets.

We first explain how we curated samples from GEO and mapped the tissue annotation to signatures. Second, we show how we use *pygenesig* to independently create and validate signatures. Third, we provide technical information about our data storage stragegy using a DBS. 

## BioQC
BioQC implements a computationally efficient Wilcoxon-Mann-Whitney (WMW) test. 

## Signatures
Tissue signatures are essentially a list of genes which are enriched in a certain tissue. 

## The Experiment in Brief

We take all samples from GEO which have the required annotation. We run BioQC on each of the samples. 
We compare the tissue predicted by BioQC with the annotated tissue and call a sample 'contamined' 
if the prediction does not match the annotation. 
