###############################################################
# Run BioQC on Archs 4 data
###############################################################

ARHCS4_META = "data/archs4/human_gsm_meta.rda"
ARCHS4_COUNTS = "./data/archs4/human_matrix.rda"
GENE_LENGTHS = "./data/ensembl_v90/gene_length_human.bed"
ENSG_TO_SYMBOL = "./data/ensembl_v90/ensg_to_symbol_human.tsv"

NORMALIZE_TISSUE_TABLE = "./data/bioqc_geo_oracle_dump/BIOQC_NORMALIZE_TISSUES_DATA_TABLE.csv"
TISSUE_SET_TABLE = "./data/bioqc_geo_oracle_dump/BIOQC_TISSUE_SET_DATA_TABLE.csv"
CHUNK_SIZE = 200
OUT_DIR = "./data/archs4/processed/human_chunks"

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

# met metadata
sample_df = lapply(gsmMeta, function(x) {
  tibble_row(GSM=x$Sample_geo_accession, GPL=x$Sample_platform_id, source_name_ch1=x$Sample_source_name_ch1)
}) %>% bind_rows()

archs4_meta = sample_df %>%
  mutate(source_name_ch1 = tolower(source_name_ch1)) %>%
  inner_join(normalize_tissue, by=c("source_name_ch1"="TISSUE_ORIG")) %>%
  inner_join(bioqc_tissue_set) %>%
  select(-source_name_ch1)

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

# Build eset
fdata = data.frame(gene_symbol=rownames(exp_tpm))
rownames(fdata) = fdata$gene_symbol
eset = ExpressionSet(exp_tpm, featureData = AnnotatedDataFrame(fdata))


## run BioQC
lapply(seq(1, ncol(eset), CHUNK_SIZE), function(i) {
  tmp_eset = eset[, i:min(ncol(eset), (i+CHUNK_SIZE))]
  save(tmp_eset, file = file.path(OUT_DIR, paste0("chunk_", i, ".rda")), compress = FALSE)
})

# Run this in a jobscript:
# wmwTest(tmp_eset, signatures_gmt, col="gene_symbol") %>% as_tibble(rownames="signature") %>% pivot_longer(cols=-signature, names_to="GSM", values_to="pvalue")


