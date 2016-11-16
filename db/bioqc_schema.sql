/** Table Tissues: holds manually selected tissues **/
create table bioqc_tissues(id varchar2(80) not null primary key);

create table bioqc_signatures(id varchar2(80) not null primary key,
    description clob null,
    gene_symbols clob null);

create table bioqc_tissues_signatures(tissue varchar2(80),
    signature varchar2(80), 
    constraint pk_tissue_signature
      primary key(tissue, signature),
    constraint fk_tissue
      foreign key (tissue)
      references bioqc_tissues("ID"),
    constraint fk_signature
      foreign key (signature)
      references bioqc_signatures("ID"));

create table bioqc_res(gse varchar2(10) not null, 
    gsm varchar2(10) not null,
    signature varchar2(80),
    pvalue binary_double,
    primary key(gse, gsm, signature),
    constraint fk_gse
      foreign key (gse)
      references bioqc_gse(gse),
    constraint fk_gsm
      foreign key (gsm)
      references bioqc_gsm(gsm),
    constraint fk_bioqc_res_signature
      foreign key (signature)
      references bioqc_signatures("ID"));

