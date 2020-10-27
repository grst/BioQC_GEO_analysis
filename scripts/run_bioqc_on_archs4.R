#!/usr/bin/env Rscript
#$ -pe smp 2
#$ -cwd
#$ -V
#$ -l h_vmem 100G

###############################################################
# Run BioQC on Archs 4 data
# ./run_bioqc_on_archs4.R SPECIES
# where SPECIES = human|mouse
###############################################################

args  = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("Specify either 'human' or 'mouse' as command line argument. ")
} else {
  species = args[1]
}

if(species == "human") {

  ARCHS4_META = "data/archs4/human_gsm_meta.rda"
  ARCHS4_COUNTS = "./data/archs4/human_matrix.rda"
  GENE_LENGTHS = "./data/ensembl_v90/gene_length_human.bed"
  ENSG_TO_SYMBOL = "./data/ensembl_v90/ensg_to_symbol_human.tsv"

  NORMALIZE_TISSUE_TABLE = "./data/bioqc_geo_oracle_dump/BIOQC_NORMALIZE_TISSUES_DATA_TABLE.csv"
  TISSUE_SET_TABLE = "./data/bioqc_geo_oracle_dump/BIOQC_TISSUE_SET_DATA_TABLE.csv"
  CHUNK_SIZE = 200
  OUT_DIR = "./data/archs4/processed/human_chunks"
  MOUSE = FALSE
  HOMOLOGENE = "./manual_annotation/homologene.data"

} else if(species == "mouse") {
  ARCHS4_META = "data/archs4/mouse_gsm_meta.rda"
  ARCHS4_COUNTS = "./data/archs4/mouse_matrix.rda"
  GENE_LENGTHS = "./data/ensembl_v90/gene_length_mouse.bed"
  ENSG_TO_SYMBOL = "./data/ensembl_v90/ensg_to_symbol_mouse.tsv"

  NORMALIZE_TISSUE_TABLE = "./data/bioqc_geo_oracle_dump/BIOQC_NORMALIZE_TISSUES_DATA_TABLE.csv"
  TISSUE_SET_TABLE = "./data/bioqc_geo_oracle_dump/BIOQC_TISSUE_SET_DATA_TABLE.csv"
  CHUNK_SIZE = 200
  OUT_DIR = "./data/archs4/processed/mouse_chunks"
  MOUSE = TRUE
  HOMOLOGENE = "./manual_annotation/homologene.data"

} else {
  stop("Invalid arguments. ")
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive=TRUE)

#' Convert counts to TPM, the simple way (not as accurate)
#'
#' From https://gist.github.com/slowkow/6e34ccb4d1311b8fe62e
#'
#' @param counts a numeric matrix of raw counts (gene x sample)
#' @param lengths a vector of gene lenghts (same order and dimention as `gene`)
#' @return tpm a numeric matrix normalized by library size and gene length
counts_to_tpm_simple = function(counts, lengths) {
  apply(counts, 2, function(x) {
    rate <- x / lengths
    rate / sum(rate) * 1e6
  })
}

library(BioQC)
library(dplyr)
library(tibble)
library(readr)
library(Biobase)
library(stringr)
library(tidyr)
library(parallel)


# Load input data
gene_lengths = read_tsv(GENE_LENGTHS)
ensg_to_symbol = read_tsv(ENSG_TO_SYMBOL)

normalize_tissue = read_csv(NORMALIZE_TISSUE_TABLE)
bioqc_tissue_set = read_csv(TISSUE_SET_TABLE) %>%
  select(TISSUE, TGROUP, TISSUE_SET) %>%
  distinct()

load(ARCHS4_META)
load(ARCHS4_COUNTS)

# get metadata
sample_df = lapply(gsmMeta, function(x) {
  tibble_row(
    GSM = x$Sample_geo_accession,
    GPL = x$Sample_platform_id,
    GSE = x$Sample_series_id[1],
    source_name_ch1 = x$Sample_source_name_ch1,
    library_strategy = x$Sample_library_strategy,
    library_source = x$Sample_library_source,
    instrument_model = x$Sample_instrument_model,
    contact_country = x$Sample_contact_country,
    submission_date = x$Sample_submission_date,
    last_update_date = x$Sample_last_update_date,
    molecule_ch1 = x$Sample_molecule_ch1,
    sample_title = x$Sample_title,
    sample_characteristics_ch1 = paste0(x$Sample_characteristics_ch1, collapse=";"),
    sample_description = paste0(x$Sample_description, collapse=";"),
    sample_data_processing = paste0(x$Sample_data_processing, collapse=";"),
    sample_extraction_protocol = paste0(x$Sample_extract_protocol_ch1, collapse=";"),
    sample_growth_protocol = paste0(x$Sample_growth_protocol_ch1, collapse=";")
  )
}) %>% bind_rows()

filtering_stats = list()
filtering_stats[[1]] = tibble_row(step="unfiltered", n_samples=nrow(sample_df), n_studies=length(unique(sample_df$GSE)))

sample_df = sample_df %>% filter(library_strategy == "RNA-Seq", library_source == "transcriptomic", molecule_ch1 %in% c("total RNA", "polyA RNA"))

filtering_stats[[2]] = tibble_row(step="library_type", n_samples=nrow(sample_df), n_studies=length(unique(sample_df$GSE)))

sample_df = sample_df %>% filter(!str_detect(str_to_lower(sample_df$sample_title), "single cell|single-cell|smart-seq|smartseq"),
                                 !str_detect(str_to_lower(sample_df$sample_characteristics_ch1), "single cell|single-cell|smart-seq|smartseq"),
                                 !str_detect(str_to_lower(sample_df$sample_description), "single cell|single-cell|smart-seq|smartseq"),
                                 !str_detect(str_to_lower(sample_df$sample_data_processing), "single cell|single-cell|smart-seq|smartseq"),
                                 !str_detect(str_to_lower(sample_df$sample_extraction_protocol), "single cell|single-cell|smart-seq|smartseq"),
                                 !str_detect(str_to_lower(sample_df$sample_growth_protocol), "single cell|single-cell|smart-seq|smartseq"))

filtering_stats[[3]] = tibble_row(step="no_single_cells", n_samples=nrow(sample_df), n_studies=length(unique(sample_df$GSE)))

read_counts = colSums(exp)
sample_df = sample_df %>% filter(GSM %in% names(read_counts[read_counts > 500000]))

filtering_stats[[4]] = tibble_row(step="200k_reads", n_samples=nrow(sample_df), n_studies=length(unique(sample_df$GSE)))


archs4_meta = sample_df %>%
  mutate(source_name_ch1 = tolower(source_name_ch1)) %>%
  inner_join(normalize_tissue, by=c("source_name_ch1"="TISSUE_ORIG")) %>%
  inner_join(bioqc_tissue_set)

filtering_stats[[5]] = tibble_row(step="tissue_in_cv", n_samples=length(unique(archs4_meta$GSM)), n_studies=length(unique(archs4_meta$GSE)))

archs4_meta %>%
  write_tsv(file.path(OUT_DIR, "archs4_meta.tsv"))

filtering_stats %>% bind_rows() %>% write_tsv(file.path(OUT_DIR, "filtering_stats.tsv"))

# Select samples with tissue from archs4 matrix.
exp = exp[, archs4_meta$GSM %>% unique()]

# normalize to TPM
gene_lengths = gene_lengths %>%
  inner_join(ensg_to_symbol, by=c("gene"="Gene stable ID"))

tmp_gene_lengths = gene_lengths %>%
  filter(`Gene name` %in% rownames(exp)) %>%
  group_by(`Gene name`) %>%
  summarise(gene_length = max(longest_isoform)) %>%
  arrange(`Gene name`) %>%
  as.data.frame() %>%
  column_to_rownames("Gene name")

exp_tpm = counts_to_tpm_simple(exp[rownames(tmp_gene_lengths),], tmp_gene_lengths$gene_length)

# load homologene
if (MOUSE) {
  homologene = read_tsv("./manual_annotation/homologene.data",
                        col_names=c("HID", "Taxonomy ID", "Gene ID", "Gene Symbol", "Protein gi", "Protein accession"))
  homologene_hsa_mmu = homologene %>% filter(`Taxonomy ID` %in% c(10090, 9606)) %>%
    mutate(species = if_else(`Taxonomy ID` == 10090, "mmu", "hsa")) %>%
    select(`HID`, species, `Gene Symbol`) %>%
    pivot_wider(names_from="species", values_from="Gene Symbol", values_fn=function(x) x[[1]])
}


## Export expression sets
lapply(seq(1, ncol(exp_tpm), CHUNK_SIZE), function(i) {
  start = i
  stop = min(ncol(exp_tpm), (i+CHUNK_SIZE-1))
  message(paste0("Processing range ", start, ":", stop))
  tmp_exp_tpm = exp_tpm[, start:stop]

  # Build eset
  fdata = data.frame(gene_symbol=rownames(tmp_exp_tpm))
  rownames(fdata) = fdata$gene_symbol

  if (MOUSE) {
    fdata = fdata %>% left_join(homologene_hsa_mmu, by=c("gene_symbol"="mmu"))

    tmp_exp_tpm = tmp_exp_tpm %>%
      as_tibble()

    tmp_exp_tpm$gene_symbol = fdata$hsa

    tmp_exp_tpm = tmp_exp_tpm %>%
      filter(!is.na(gene_symbol))

    tmp_exp_tpm = tmp_exp_tpm %>% as.data.frame() %>% arrange(gene_symbol) %>% column_to_rownames("gene_symbol") %>% as.matrix()
    fdata = data.frame(gene_symbol=rownames(tmp_exp_tpm)) %>% arrange(gene_symbol)
    rownames(fdata) = fdata$gene_symbol
  }

  tmp_eset = ExpressionSet(tmp_exp_tpm, featureData = AnnotatedDataFrame(fdata))

  save(tmp_eset, file = file.path(OUT_DIR, paste0(species, "_chunk_", i, ".rda")), compress = FALSE)
})

# Run this in a jobscript:
# wmwTest(tmp_eset, signatures_gmt, col="gene_symbol") %>% as_tibble(rownames="signature") %>% pivot_longer(cols=-signature, names_to="GSM", values_to="pvalue")


