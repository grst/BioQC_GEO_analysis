import pandas as pd

sig_table = pd.read_csv("../data/bioqc_geo_oracle_dump/BIOQC_SIGNATURES_DATA_TABLE.csv")

sig_table


def make_gmt(df, out_path):
    gmt = df.loc[:, ["NAME", "DESCRIPTION", "GENE_SYMBOLS"]]
    gmt["GENE_SYMBOLS"] = [",".join(g for g in x.split(",") if g !="") for x in gmt["GENE_SYMBOLS"]]
    gmt.to_csv(out_path, header=False, index=False)


for filename, df_group in sig_table.groupby("SOURCE"):
    make_gmt(df_group, filename)    


