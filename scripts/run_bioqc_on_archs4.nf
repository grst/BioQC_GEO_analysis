#!/usr/bin/env nextflow

params.chunk_path = "./data/archs4/processed/**/chunk*.rda"
params.signature_table = "./data/bioqc_geo_oracle_dump/BIOQC_SIGNATURES_DATA_TABLE.csv"
params.out_dir = "./data/archs4/processed"

process run_bioqc {
    conda "bioconductor-bioqc r-readr r-dplyr r-tidyr"

    input:
        file eset_chunk from Channel.fromPath(params.chunk_path)
        file "signature_table.csv" from Channel.fromPath(params.signature_table).collect()

    output:
        file "${eset_chunk}.tsv" into bioqc_res_chunk

    script:
        """#!/usr/bin/env Rscript
        library(dplyr)
        library(BioQC)
        library(readr)
        library(tidyr)

        load("${eset_chunk}")

        signatures = read_csv("signature_table.csv")

        signatures_gmt = lapply(1:nrow(signatures), function(i) {
          id = signatures\$ID[i]
          source = signatures\$SOURCE[i]
          symbols = signatures\$GENE_SYMBOLS[i]
          genes = strsplit(symbols, ",")[[1]]
          genes = genes[genes != ""]
          tmp_l = list(name=id, desc=source, genes=genes)
          tmp_l
        }) %>% GmtList()

        wmwTest(tmp_eset, signatures_gmt, col="gene_symbol") %>%
            as_tibble(rownames="signature") %>%
            pivot_longer(cols=-signature, names_to="GSM", values_to="pvalue") %>%
            write_tsv("${eset_chunk}.tsv")
        """
}

process concat {
    publishDir params.out_dir, mode: "link"
    input:
        file chunks from bioqc_res_chunk.collect()

    output:
        file "bioqc_res_all.tsv"

    script:
        """
        head -n1 <(cat *.rda.tsv) > bioqc_res_all.tsv
        for f in *.rda.tsv ; do
          tail -n+2 \$f >> bioqc_res_all.tsv
        done
        """

}
